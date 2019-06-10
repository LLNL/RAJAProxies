/// \file
/// Computes forces for the 12-6 Lennard Jones (LJ) potential.
///
/// The Lennard-Jones model is not a good representation for the
/// bonding in copper, its use has been limited to constant volume
/// simulations where the embedding energy contribution to the cohesive
/// energy is not included in the two-body potential
///
/// The parameters here are taken from Wolf and Phillpot and fit to the
/// room temperature lattice constant and the bulk melt temperature
/// Ref: D. Wolf and S.Yip eds. Materials Interfaces (Chapman & Hall
///      1992) Page 230.
///
/// Notes on LJ:
///
/// http://en.wikipedia.org/wiki/Lennard_Jones_potential
///
/// The total inter-atomic potential energy in the LJ model is:
///
/// \f[
///   E_{tot} = \sum_{ij} U_{LJ}(r_{ij})
/// \f]
/// \f[
///   U_{LJ}(r_{ij}) = 4 \epsilon
///           \left\{ \left(\frac{\sigma}{r_{ij}}\right)^{12}
///           - \left(\frac{\sigma}{r_{ij}}\right)^6 \right\}
/// \f]
///
/// where \f$\epsilon\f$ and \f$\sigma\f$ are the material parameters in the potential.
///    - \f$\epsilon\f$ = well depth
///    - \f$\sigma\f$   = hard sphere diameter
///
///  To limit the interation range, the LJ potential is typically
///  truncated to zero at some cutoff distance. A common choice for the
///  cutoff distance is 2.5 * \f$\sigma\f$.
///  This implementation can optionally shift the potential slightly
///  upward so the value of the potential is zero at the cuotff
///  distance.  This shift has no effect on the particle dynamics.
///
///
/// The force on atom i is given by
///
/// \f[
///   F_i = -\nabla_i \sum_{jk} U_{LJ}(r_{jk})
/// \f]
///
/// where the subsrcipt i on the gradient operator indicates that the
/// derivatives are taken with respect to the coordinates of atom i.
/// Liberal use of the chain rule leads to the expression
///
/// \f{eqnarray*}{
///   F_i &=& - \sum_j U'_{LJ}(r_{ij})\hat{r}_{ij}\\
///       &=& \sum_j 24 \frac{\epsilon}{r_{ij}} \left\{ 2 \left(\frac{\sigma}{r_{ij}}\right)^{12}
///               - \left(\frac{\sigma}{r_{ij}}\right)^6 \right\} \hat{r}_{ij}
/// \f}
///
/// where \f$\hat{r}_{ij}\f$ is a unit vector in the direction from atom
/// i to atom j.
///
///

#include "ljForce.h"

#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "constants.h"
#include "mytype.h"
#include "parallel.h"
#include "linkCells.h"
#include "memUtils.h"
#include "CoMDTypes.h"

#define POT_SHIFT 1.0

/// Derived struct for a Lennard Jones potential.
/// Polymorphic with BasePotential.
/// \see BasePotential
typedef struct LjPotentialSt
{
   real_t cutoff;          //!< potential cutoff distance in Angstroms
   real_t mass;            //!< mass of atoms in intenal units
   real_t lat;             //!< lattice spacing (angs) of unit cell
   char latticeType[8];    //!< lattice type, e.g. FCC, BCC, etc.
   char  name[3];	   //!< element name
   int	 atomicNo;	   //!< atomic number
   int  (*force)(SimFlat* s); //!< function pointer to force routine
   void (*print)(FILE* file, BasePotential* pot);
   void (*destroy)(BasePotential** pot); //!< destruction of the potential
   real_t sigma;
   real_t epsilon;
} LjPotential;

static int ljForce(SimFlat* s);
static void ljPrint(FILE* file, BasePotential* pot);

void ljDestroy(BasePotential** inppot)
{
   if ( ! inppot ) return;
   LjPotential* pot = (LjPotential*)(*inppot);
   if ( ! pot ) return;
   comdFree(pot);
   *inppot = NULL;

   return;
}

/// Initialize an Lennard Jones potential for Copper.
BasePotential* initLjPot(void)
{
   LjPotential *pot = (LjPotential*)comdMalloc(sizeof(LjPotential));
   pot->force = ljForce;
   pot->print = ljPrint;
   pot->destroy = ljDestroy;
   pot->sigma = 2.315;	                  // Angstrom
   pot->epsilon = 0.167;                  // eV
   pot->mass = 63.55 * amuToInternalMass; // Atomic Mass Units (amu)

   pot->lat = 3.615;                      // Equilibrium lattice const in Angs
   strcpy(pot->latticeType, "FCC");       // lattice type, i.e. FCC, BCC, etc.
   pot->cutoff = 2.5*pot->sigma;          // Potential cutoff in Angs

   strcpy(pot->name, "Cu");
   pot->atomicNo = 29;

   return (BasePotential*) pot;
}

void ljPrint(FILE* file, BasePotential* pot)
{
   LjPotential* ljPot = (LjPotential*) pot;
   fprintf(file, "  Potential type   : Lennard-Jones\n");
   fprintf(file, "  Species name     : %s\n", ljPot->name);
   fprintf(file, "  Atomic number    : %d\n", ljPot->atomicNo);
   fprintf(file, "  Mass             : " FMT1 " amu\n", ljPot->mass / amuToInternalMass); // print in amu
   fprintf(file, "  Lattice Type     : %s\n", ljPot->latticeType);
   fprintf(file, "  Lattice spacing  : " FMT1 " Angstroms\n", ljPot->lat);
   fprintf(file, "  Cutoff           : " FMT1 " Angstroms\n", ljPot->cutoff);
   fprintf(file, "  Epsilon          : " FMT1 " eV\n", ljPot->epsilon);
   fprintf(file, "  Sigma            : " FMT1 " Angstroms\n", ljPot->sigma);
}

int ljForce(SimFlat* s)
{
   LjPotential* pot = (LjPotential *) s->pot;
   const real_t sigma = pot->sigma;
   const real_t epsilon = pot->epsilon;
   const real_t rCut = pot->cutoff;
   const real_t rCut2 = rCut*rCut;
   const real_t s6 = sigma*sigma*sigma*sigma*sigma*sigma;
   const real_t rCut6 = s6 / (rCut2*rCut2*rCut2);
   const real_t eShift = POT_SHIFT * rCut6 * (rCut6 - 1.0);

   // zero forces and energy
   rajaReduceSumRealKernel ePot(0.0);

   ePotential = 0.0;

   profileStart(forceZeroingTimer);
  RAJA::kernel<atomWorkKernel>(
  RAJA::make_tuple(
    RAJA::RangeSegment(0, s->boxes->nLocalBoxes),
    RAJA::RangeSegment(0, MAXATOMS) ),
    [=] COMD_DEVICE (int iBox, int iOffLocal) {
      const int nIBox = s->boxes->nAtoms[iBox];
      if(iOffLocal < nIBox) {
        const int iOff = iOffLocal + (iBox * MAXATOMS);
        real3_ptr f = s->atoms->f;
        real_ptr U = s->atoms->U;
        f[iOff][0] = 0.0;
        f[iOff][1] = 0.0;
        f[iOff][2] = 0.0;
        U[iOff] = 0.0;
      }
    } ) ;

   profileStop(forceZeroingTimer);
   {
   profileStart(forceFunctionTimer);
     RAJA::kernel<forcePolicyKernel>(
       RAJA::make_tuple(
         *s->isLocalSegment,                // local boxes
         RAJA::RangeSegment(0,27),          // 27 neighbor boxes
         RAJA::RangeSegment(0, MAXATOMS),   // atoms i in local box
         RAJA::RangeSegment(0, MAXATOMS) ), // atoms j in neighbor box
       [=] COMD_DEVICE (int iBoxID, int nghb, int iOff, int jOff) {
         const int nLocalBoxes = s->boxes->nLocalBoxes;
         const int nIBox = s->boxes->nAtoms[iBoxID];
         const int jBoxID = s->boxes->nbrBoxes[iBoxID][nghb];
         const int nJBox = s->boxes->nAtoms[jBoxID];
         const int iOffLocal = iOff;
         const int jOffLocal = jOff;

         iOff += iBoxID*MAXATOMS;
         jOff += jBoxID*MAXATOMS;
         const int iGid = s->atoms->gid[iOff];
         const int jGid = s->atoms->gid[jOff];

         if( (iOffLocal < nIBox && jOffLocal < nJBox) && !(jBoxID < nLocalBoxes && jGid <= iGid)) {
           real3 dr;
           real_t r2 = 0.0;
           real3_ptr r =  s->atoms->r;
           for (int m=0; m<3; m++)
           {
             dr[m] = r[iOff][m] - r[jOff][m];
             r2 += dr[m]*dr[m];
           }
           if ( r2 <= rCut2 && r2 > 0.0)
           {
             // Important note:
             // from this point on r actually refers to 1.0/r
             real_ptr U = s->atoms->U ;
             real3_ptr f = s->atoms->f ;
             r2 = 1.0/r2;
             const real_t r6 = s6 * (r2*r2*r2);
             const real_t eLocal = r6 * (r6 - 1.0) - eShift;
             U[iOff] += 0.5*eLocal; // Shouldn't this be atomic too?

             if (jBoxID < nLocalBoxes)
               ePot += eLocal;
             else
               ePot += 0.5*eLocal;

             // different formulation to avoid sqrt computation
             const real_t fr = - 4.0*epsilon*r6*r2*(12.0*r6 - 6.0);

             for (int m=0; m<3; m++)
             {
               dr[m] *= fr;
#ifdef DO_CUDA
               atomicAdd(&f[iOff][m], -dr[m]);
               atomicAdd(&f[jOff][m], dr[m]);
#else
               f[iOff][m] -= dr[m];
               f[jOff][m] += dr[m];
#endif
             }
           }  //end if within cutoff
         }//end if atoms exist
       });
   profileStop(forceFunctionTimer);
   }

   ePotential = ePot*4.0*epsilon;

   return 0;
}
