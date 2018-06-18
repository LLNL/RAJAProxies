/// \file
/// CoMD data structures.

#ifndef __COMDTYPES_H_
#define __COMDTYPES_H_

#include <stdio.h>
#include "mytype.h"
#include "haloExchange.h"
#include "linkCells.h"
#include "decomposition.h"
#include "initAtoms.h"
#include "RAJA/RAJA.hpp"

struct SimFlatSt;

/// The base struct from which all potentials derive.  Think of this as an
/// abstract base class.
///
/// CoMD uses the following units:
///  - distance is in Angstroms
///  - energy is in eV
///  - time is in fs
///  - force in in eV/Angstrom
///
///  The choice of distance, energy, and time units means that the unit
///  of mass is eV*fs^2/Angstrom^2.  Hence, we must convert masses that
///  are input in AMU (atomic mass units) into internal mass units.
typedef struct BasePotentialSt
{
   real_t cutoff;          //!< potential cutoff distance in Angstroms
   real_t mass;            //!< mass of atoms in intenal units
   real_t lat;             //!< lattice spacing (angs) of unit cell
   char latticeType[8];    //!< lattice type, e.g. FCC, BCC, etc.
   char  name[3];	   //!< element name
   int	 atomicNo;	   //!< atomic number
   int  (*force)(struct SimFlatSt* s); //!< function pointer to force routine
   void (*print)(FILE* file, struct BasePotentialSt* pot);
   void (*destroy)(struct BasePotentialSt** pot); //!< destruction of the potential
} BasePotential;


/// species data: chosen to match the data found in the setfl/funcfl files
typedef struct SpeciesDataSt
{
   char  name[3];   //!< element name
   int	 atomicNo;  //!< atomic number
   real_t mass;     //!< mass in internal units
} SpeciesData;

/// Simple struct to store the initial energy and number of atoms.
/// Used to check energy conservation and atom conservation.
typedef struct ValidateSt
{
   double eTot0; //<! Initial total energy
   int nAtoms0;  //<! Initial global number of atoms
} Validate;

///
/// The fundamental simulation data structure with MAXATOMS in every
/// link cell.
///
typedef struct SimFlatSt
{
   int nSteps;            //<! number of time steps to run
   int printRate;         //<! number of steps between output
   double dt;             //<! time step

   Domain* domain;        //<! domain decomposition data

   LinkCell* boxes;       //<! link-cell data

   Atoms* atoms;          //<! atom data (positions, momenta, ...)

   SpeciesData* species;  //<! species data (per species, not per atom)

   real_t ePotential;     //!< the total potential energy of the system
   real_t eKinetic;       //!< the total kinetic energy of the system

   BasePotential *pot;	  //!< the potential

   HaloExchange* atomExchange;

   RAJA::TypedIndexSet<RAJA::RangeSegment> *isTotal ;
   RAJA::TypedIndexSet<RAJA::RangeSegment> *isLocal ;
  RAJA::RangeSegment *isLocalSegment ;

} SimFlat;


//
// RAJA execution policies
//
//typedef RAJA::omp_parallel_for_segit linkCellTraversal ;
//typedef RAJA::TypedIndexSet<RAJA::RangeSegment>::ExecPolicy<RAJA::seq_segit, RAJA::simd_exec> linkCellWork;
//typedef RAJA::TypedIndexSet<RAJA::RangeSegment>::ExecPolicy<RAJA::omp_parallel_for_segit, RAJA::simd_exec> atomWork;
//typedef RAJA::omp_parallel_segit task_graph_policy;
#ifdef ENABLE_OPENMP
typedef RAJA::omp_parallel_for_segit linkCellTraversal ;
typedef RAJA::ExecPolicy<RAJA::seq_segit, RAJA::simd_exec> linkCellWork;
typedef RAJA::ExecPolicy<RAJA::omp_parallel_for_segit, RAJA::simd_exec> atomWork;
typedef RAJA::omp_parallel_segit task_graph_policy;

typedef RAJA::KernelPolicy<
  RAJA::statement::For<0, RAJA::omp_parallel_for_segit,
  RAJA::statement::For<1, RAJA::seq_exec,
  RAJA::statement::For<2, RAJA::seq_exec,
  RAJA::statement::For<3, RAJA::simd_exec,
  RAJA::statement::Lambda<0> > > > > > ljForcePolicy;

typedef RAJA::ReduceSum<RAJA::omp_reduce, real_t> rajaReduceSumReal;
#else
typedef RAJA::seq_segit linkCellTraversal ;
typedef RAJA::ExecPolicy<RAJA::seq_segit, RAJA::simd_exec> linkCellWork;
typedef RAJA::ExecPolicy<RAJA::seq_segit, RAJA::simd_exec> atomWork;
typedef RAJA::seq_segit task_graph_policy;

typedef RAJA::KernelPolicy<
  RAJA::statement::For<0, RAJA::simd_exec,
  RAJA::statement::For<1, RAJA::simd_exec,
  RAJA::statement::For<2, RAJA::simd_exec,
  RAJA::statement::For<3, RAJA::simd_exec,
  RAJA::statement::Lambda<0> > > > > > ljForcePolicy;

typedef RAJA::ReduceSum<RAJA::seq_reduce, real_t> rajaReduceSumReal;
#endif
// "task-graph" requires segment dependency graph is set up...not yet...
//typedef RAJA::omp_taskgraph_segit task_graph_policy;
//typedef RAJA::seq_segit task_graph_policy;

#endif
