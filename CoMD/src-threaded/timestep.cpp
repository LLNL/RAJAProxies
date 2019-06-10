/// \file
/// Leapfrog time integrator

#include "timestep.h"

#include "CoMDTypes.h"
#include "linkCells.h"
#include "parallel.h"
#include "performanceTimers.h"

static void advanceVelocity(SimFlat* s, RAJA::TypedIndexSet<RAJA::RangeSegment> *extent, real_t dt);
static void advancePosition(SimFlat* s, RAJA::TypedIndexSet<RAJA::RangeSegment> *extent, real_t dt);

/* This structure is redefined here in order to inline some functions from
 * linkCells.cpp here.
 */
typedef struct AtomMsgSt
{
   int gid;
   int type;
   real_t rx, ry, rz;
   real_t px, py, pz;
}
AtomMsg;

/// Advance the simulation time to t+dt using a leap frog method
/// (equivalent to velocity verlet).
///
/// Forces must be computed before calling the integrator the first time.
///
///  - Advance velocities half time step using forces
///  - Advance positions full time step using velocities
///  - Update link cells and exchange remote particles
///  - Compute forces
///  - Update velocities half time step using forces
///
/// This leaves positions, velocities, and forces at t+dt, with the
/// forces ready to perform the half step velocity update at the top of
/// the next call.
///
/// After nSteps the kinetic energy is computed for diagnostic output.
double timestep(SimFlat* s, int nSteps, real_t dt)
{
   for (int ii=0; ii<nSteps; ++ii)
   {
      startTimer(velocityTimer);
      advanceVelocity(s, s->isLocal, 0.5*dt);
      stopTimer(velocityTimer);

      startTimer(positionTimer);
      advancePosition(s, s->isLocal, dt);
      stopTimer(positionTimer);

      startTimer(redistributeTimer);
      redistributeAtoms(s);
      stopTimer(redistributeTimer);

      startTimer(computeForceTimer);
      computeForce(s);
      stopTimer(computeForceTimer);

      startTimer(velocityTimer);
      advanceVelocity(s, s->isLocal, 0.5*dt);
      stopTimer(velocityTimer);
   }

   startTimer(kineticEnergyTimer);
   kineticEnergy(s);
   stopTimer(kineticEnergyTimer);

   return ePotential;
}

void computeForce(SimFlat* s)
{
   s->pot->force(s);
}


void advanceVelocity(SimFlat* s, RAJA::TypedIndexSet<RAJA::RangeSegment> *extent, real_t dt)
{
  RAJA::kernel<atomWorkKernel>(
  RAJA::make_tuple(
    RAJA::RangeSegment(0, s->boxes->nLocalBoxes),
    RAJA::RangeSegment(0, MAXATOMS) ),
    [=] COMD_DEVICE (int iBox, int iOffLocal) {
      const int nIBox = s->boxes->nAtoms[iBox];
      if(iOffLocal < nIBox) {
        const int iOff = iOffLocal + (iBox * MAXATOMS);
        real3_ptr p = s->atoms->p;
        real3_ptr f = s->atoms->f;

        p[iOff][0] += dt*f[iOff][0];
        p[iOff][1] += dt*f[iOff][1];
        p[iOff][2] += dt*f[iOff][2];
      }
    } );
}

void advancePosition(SimFlat* s, RAJA::TypedIndexSet<RAJA::RangeSegment> *extent, real_t dt)
{
  RAJA::kernel<atomWorkKernel>(
  RAJA::make_tuple(
  RAJA::RangeSegment(0, s->boxes->nLocalBoxes),
    RAJA::RangeSegment(0, MAXATOMS) ),
    [=] COMD_DEVICE (int iBox, int iOffLocal) {
      const int nIBox = s->boxes->nAtoms[iBox];
      if(iOffLocal < nIBox) {
        const int iOff = iOffLocal + (iBox * MAXATOMS);
        const int iSpecies = s->atoms->iSpecies[iOff];
        const real_t invMass = 1.0/s->species[iSpecies].mass;
        real3_ptr p = s->atoms->p;
        real3_ptr r = s->atoms->r;

        r[iOff][0] += dt*p[iOff][0]*invMass;
        r[iOff][1] += dt*p[iOff][1]*invMass;
        r[iOff][2] += dt*p[iOff][2]*invMass;
      }
    } ) ;
}

/// Calculates total kinetic and potential energy across all tasks.  The
/// local potential energy is a by-product of the force routine.
void kineticEnergy(SimFlat* s)
{
   real_t eLocal[2];
   eLocal[0] = ePotential;

   rajaReduceSumRealKernel kenergy(0.0);

   eLocal[1] = 0.0;

  RAJA::kernel<atomWorkKernel>(
  RAJA::make_tuple(
    RAJA::RangeSegment(0, s->boxes->nLocalBoxes),
    RAJA::RangeSegment(0, MAXATOMS) ),
    [=] COMD_DEVICE (int iBox, int iOffLocal) {
      const int nIBox = s->boxes->nAtoms[iBox];
      if(iOffLocal < nIBox) {
        const int iOff = iOffLocal + (iBox * MAXATOMS);
        const int iSpecies = s->atoms->iSpecies[iOff];
        const real_t invMass = 0.5/s->species[iSpecies].mass;
        real3_ptr p = s->atoms->p;

        kenergy += ( p[iOff][0] * p[iOff][0] +
                     p[iOff][1] * p[iOff][1] +
                     p[iOff][2] * p[iOff][2] )*
                     invMass ;
      }
    } ) ;

   eLocal[1] = kenergy;
// printf("kenergy: %f\n", eLocal[1]);

   real_t eSum[2];
   startTimer(commReduceTimer);
   addRealParallel(eLocal, eSum, 2);
   stopTimer(commReduceTimer);

   ePotential = eSum[0];
   eKinetic = eSum[1];
}

/// \details
/// The force exchange assumes that the atoms are in the same order in
/// both a given local link cell and the corresponding remote cell(s).
/// However, the atom exchange does not guarantee this property,
/// especially when atoms cross a domain decomposition boundary and move
/// from one task to another.  Trying to maintain the atom order during
/// the atom exchange would immensely complicate that code.  Instead, we
/// just sort the atoms after the atom exchange.
COMD_HOST_DEVICE void sortAtomsInCell(Atoms* atoms, LinkCell* boxes, int iBox)
{
   int nAtoms = boxes->nAtoms[iBox];

   //AtomMsg tmp[nAtoms];
   AtomMsg tmp[MAXATOMS];

   int begin = iBox*MAXATOMS;
   int end = begin + nAtoms;

   int sorted = 1;
   for (int ii=begin, iTmp=0; ii<end; ++ii, ++iTmp)
   {
      tmp[iTmp].gid  = atoms->gid[ii];
      tmp[iTmp].type = atoms->iSpecies[ii];
      tmp[iTmp].rx =   atoms->r[ii][0];
      tmp[iTmp].ry =   atoms->r[ii][1];
      tmp[iTmp].rz =   atoms->r[ii][2];
      tmp[iTmp].px =   atoms->p[ii][0];
      tmp[iTmp].py =   atoms->p[ii][1];
      tmp[iTmp].pz =   atoms->p[ii][2];
      if(iTmp > 0) {
        if(tmp[iTmp].gid < tmp[iTmp-1].gid)
          sorted = 0;
      }
   }
   if(!sorted) {
     return;
   }
   //qsort(&tmp, nAtoms, sizeof(AtomMsg), sortAtomsById);
   /* TODO: Make this into a function like qsort to clean things up.
    * Note: Switching from qsort to insertion sort roughly doubles sorting performance
    */
     /* Begin Insertion Sort */
     /* This has been changed to an insertion sort instead of a quick sort because the
      * elements of this array are already mostly sorted and insertion sort works well
      * on mostly sorted data.
      */
     int i, j;
     AtomMsg key;
     for (i = 1; i < nAtoms; i++) {
       key.gid  = tmp[i].gid;
       key.type = tmp[i].type;
       key.rx   = tmp[i].rx;
       key.ry   = tmp[i].ry;
       key.rz   = tmp[i].rz;
       key.px   = tmp[i].px;
       key.py   = tmp[i].py;
       key.pz   = tmp[i].pz;
       
       j = i-1;
 
       /* Move elements of tmp[0..i-1], that are
          greater than key, to one position ahead
          of their current position */
       while (j >= 0 && tmp[j].gid > key.gid)
       {
         tmp[j+1].gid  = tmp[j].gid;
         tmp[j+1].type = tmp[j].type;
         tmp[j+1].rx   = tmp[j].rx;
         tmp[j+1].ry   = tmp[j].ry;
         tmp[j+1].rz   = tmp[j].rz;
         tmp[j+1].px   = tmp[j].px;
         tmp[j+1].py   = tmp[j].py;
         tmp[j+1].pz   = tmp[j].pz;
         j = j-1;
       }
       tmp[j+1].gid  = key.gid;
       tmp[j+1].type = key.type;
       tmp[j+1].rx   = key.rx;
       tmp[j+1].ry   = key.ry;
       tmp[j+1].rz   = key.rz;
       tmp[j+1].px   = key.px;
       tmp[j+1].py   = key.py;
       tmp[j+1].pz   = key.pz;
     }
     /* End Insertion Sort */

   for (int ii=begin, iTmp=0; ii<end; ++ii, ++iTmp)
   {
      atoms->gid[ii]   = tmp[iTmp].gid;
      atoms->iSpecies[ii] = tmp[iTmp].type;
      atoms->r[ii][0]  = tmp[iTmp].rx;
      atoms->r[ii][1]  = tmp[iTmp].ry;
      atoms->r[ii][2]  = tmp[iTmp].rz;
      atoms->p[ii][0]  = tmp[iTmp].px;
      atoms->p[ii][1]  = tmp[iTmp].py;
      atoms->p[ii][2]  = tmp[iTmp].pz;
   }

}

/// \details
/// This function provides one-stop shopping for the sequence of events
/// that must occur for a proper exchange of halo atoms after the atom
/// positions have been updated by the integrator.
///
/// - updateLinkCells: Since atoms have moved, some may be in the wrong
///   link cells.
/// - haloExchange (atom version): Sends atom data to remote tasks.
/// - sort: Sort the atoms.
///
/// \see updateLinkCells
/// \see initAtomHaloExchange
/// \see sortAtomsInCell
void redistributeAtoms(SimFlat* sim)
{
   startTimer(updateLinkCellsTimer);
   updateLinkCells(sim->boxes, sim->atoms);
   stopTimer(updateLinkCellsTimer);

   startTimer(atomHaloTimer);
   haloExchange(sim->atomExchange, sim, 0);
   stopTimer(atomHaloTimer);

   startTimer(atomSortTimer);
   RAJA::kernel<redistributeKernel>(
   RAJA::make_tuple(
   RAJA::RangeSegment(0, sim->boxes->nTotalBoxes)),
   [=] COMD_DEVICE (int iBox) {
     sortAtomsInCell(sim->atoms, sim->boxes, iBox);
   } );
   stopTimer(atomSortTimer);
}

