/// \file
/// Leapfrog time integrator

#include "timestep.h"

#include <omp.h>

#include "CoMDTypes.h"
#include "linkCells.h"
#include "parallel.h"
#include "performanceTimers.h"

static void advanceVelocity(SimFlat* s, RAJA::TypedIndexSet<RAJA::RangeSegment> *extent, real_t dt);
static void advancePosition(SimFlat* s, RAJA::TypedIndexSet<RAJA::RangeSegment> *extent, real_t dt);

/* This structure is redefined here in order to inline some functions from
 * linkCells.cpp here.  TODO: Mark the functions needed as both host and
 * device functions with a macro if DO_CUDA is defined.
 */
typedef struct AtomMsgSt
{
   int gid;
   int type;
   real_t rx, ry, rz;
   real_t px, py, pz;
}
AtomMsg;

/* This variable acts as a flag as to whether the timestep function has been
 * called yet or not.  This changes the functionality of the kineticEnergy
 * function depending on whether initalization is finished or not.
 */
int inTimestep = 0;

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
  if(!inTimestep)
    inTimestep = 1;
   for (int ii=0; ii<nSteps; ++ii)
   {
      startTimer(velocityTimer);
#ifdef DO_CUDA
      advanceVelocity(s, globalSim->isLocal, 0.5*dt);
#else
      advanceVelocity(s, s->isLocal, 0.5*dt);
#endif
      stopTimer(velocityTimer);

      startTimer(positionTimer);
#ifdef DO_CUDA
      advancePosition(s, globalSim->isLocal, dt);
#else
      advancePosition(s, s->isLocal, dt);
#endif
      stopTimer(positionTimer);

      startTimer(redistributeTimer);
      redistributeAtoms(s);
      stopTimer(redistributeTimer);

      startTimer(computeForceTimer);
      computeForce(s);
      stopTimer(computeForceTimer);

      startTimer(velocityTimer);
#ifdef DO_CUDA
      advanceVelocity(s, globalSim->isLocal, 0.5*dt);
#else
      advanceVelocity(s, s->isLocal, 0.5*dt);
#endif
      stopTimer(velocityTimer);
   }

   kineticEnergy(s);

#ifdef DO_CUDA
   return globalSim->ePotential;
#else
   return s->ePotential;
#endif
}

void computeForce(SimFlat* s)
{
#ifdef DO_CUDA
   globalSim->pot->force(s);
#else
   s->pot->force(s);
#endif
}


void advanceVelocity(SimFlat* s, RAJA::TypedIndexSet<RAJA::RangeSegment> *extent, real_t dt)
{
#ifdef DO_CUDA
  RAJA::kernel<atomWorkGPU>(
  RAJA::make_tuple(
    RAJA::RangeSegment(0, globalSim->boxes->nLocalBoxes),
    RAJA::RangeSegment(0, MAXATOMS) ),
    [=] RAJA_DEVICE (int iBox, int iOffLocal) {
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
#else
   RAJA::forall<atomWork>(*extent, [=] (int iOff) {
      s->atoms->p[iOff][0] += dt*s->atoms->f[iOff][0];
      s->atoms->p[iOff][1] += dt*s->atoms->f[iOff][1];
      s->atoms->p[iOff][2] += dt*s->atoms->f[iOff][2];
   } ) ;
#endif
}

void advancePosition(SimFlat* s, RAJA::TypedIndexSet<RAJA::RangeSegment> *extent, real_t dt)
{
#ifdef DO_CUDA
  RAJA::kernel<atomWorkGPU>(
  RAJA::make_tuple(
  RAJA::RangeSegment(0, globalSim->boxes->nLocalBoxes),
    RAJA::RangeSegment(0, MAXATOMS) ),
    [=] RAJA_DEVICE (int iBox, int iOffLocal) {
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
#else
   RAJA::forall<atomWork>(*extent, [=] (int iOff) {
      int iSpecies = s->atoms->iSpecies[iOff];
      real_t invMass = 1.0/s->species[iSpecies].mass;
      s->atoms->r[iOff][0] += dt*s->atoms->p[iOff][0]*invMass;
      s->atoms->r[iOff][1] += dt*s->atoms->p[iOff][1]*invMass;
      s->atoms->r[iOff][2] += dt*s->atoms->p[iOff][2]*invMass;
   } ) ;
#endif
}

/// Calculates total kinetic and potential energy across all tasks.  The
/// local potential energy is a by-product of the force routine.
void kineticEnergy(SimFlat* s)
{
   real_t eLocal[2];
#ifdef DO_CUDA
   rajaReduceSumRealCUDA kenergy(0.0);
   eLocal[0] = globalSim->ePotential;
#else
   rajaReduceSumReal kenergy(0.0);
   eLocal[0] = s->ePotential;
#endif

   eLocal[1] = 0;
#ifdef DO_CUDA
   int localBoxes;
   if(inTimestep)
     localBoxes = globalSim->boxes->nLocalBoxes;
   else
     localBoxes = s->boxes->nLocalBoxes;
  RAJA::kernel<atomWorkGPU>(
  RAJA::make_tuple(
    //RAJA::RangeSegment(0, s->boxes->nLocalBoxes),
    RAJA::RangeSegment(0, localBoxes),
    RAJA::RangeSegment(0, MAXATOMS) ),
    [=] RAJA_DEVICE (int iBox, int iOffLocal) {
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
#else
   RAJA::forall<atomWork>(*s->isLocal, [=] (int iOff) {
      const int iSpecies = s->atoms->iSpecies[iOff];
      const real_t invMass = 0.5/s->species[iSpecies].mass;
      kenergy += ( s->atoms->p[iOff][0] * s->atoms->p[iOff][0] +
                   s->atoms->p[iOff][1] * s->atoms->p[iOff][1] +
                   s->atoms->p[iOff][2] * s->atoms->p[iOff][2] )*
                  invMass ;
   } ) ;
#endif
   eLocal[1] = kenergy;

   real_t eSum[2];
   startTimer(commReduceTimer);
   addRealParallel(eLocal, eSum, 2);
   stopTimer(commReduceTimer);

#ifdef DO_CUDA
   globalSim->ePotential = eSum[0];
   globalSim->eKinetic = eSum[1];
#else
   s->ePotential = eSum[0];
   s->eKinetic = eSum[1];
#endif
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
#ifdef DO_CUDA
   /* The updateLinkCells function has been altered to only accept sim
    * in order to avoid dereferencing any of the members of sim in CPU
    * code which could trigger an unnecessary page fault.
    */
   updateLinkCells(sim);
#else
   updateLinkCells(sim->boxes, sim->atoms);
#endif
   stopTimer(updateLinkCellsTimer);

   startTimer(atomHaloTimer);
#ifdef DO_CUDA
   haloExchange(globalSim->atomExchange, sim);
#else
   haloExchange(sim->atomExchange, sim);
#endif
   stopTimer(atomHaloTimer);

#ifdef DO_CUDA
   RAJA::kernel<redistributeGPU>(
   RAJA::make_tuple(
   RAJA::RangeSegment(0, globalSim->boxes->nTotalBoxes)),
   [=] RAJA_DEVICE (int iBox) {
     Atoms *atoms = sim->atoms;
     LinkCell *boxes = sim->boxes;
     int nAtoms = boxes->nAtoms[iBox];

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
       if(iTmp > 0){
         if(tmp[iTmp].gid < tmp[iTmp-1].gid)
           sorted = 0;
       }
     }
     if(!sorted){
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
   } );
#else
   RAJA::forall<linkCellTraversal>(RAJA::RangeSegment(0,sim->boxes->nTotalBoxes), [=] (int ii) {
     sortAtomsInCell(sim->atoms, sim->boxes, ii);
   });
#endif

}

void updateIndexSets(SimFlat *s)
{
//TO DO rewrite this
/*
#pragma omp parallel for
   for (int i=0; i<s->boxes->nTotalBoxes; ++i) {
      RAJA::RangeSegment *seg =
         static_cast<RAJA::RangeSegment *>(s->isTotal->getSegment(i)) ;

      seg->setEnd(i*MAXATOMS + s->boxes->nAtoms[i]) ;
   }
*/
}

