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
#include "performanceTimers.h"
#include "RAJA/RAJA.hpp"
//#include "RAJA/policy/cuda/kernel/CudaKernel.hpp"

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

#ifdef DO_CUDA
/* This will become essentially a copy of the SimFlat structure that is passed throughout many of
 * the CoMD functions.  Only the constant values and variables only needed on the CPU will be
 * copied into this structure, none of the atom data or large arrays.
 *
 * The purpose of this is to avoid accessing the simulation structure on the CPU side after it has
 * been transferred to GPU memory.  This creates unnecessary page faults and drastically reduces
 * the performance of this code.
 */
extern SimFlat *globalSim;
#endif

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
  RAJA::statement::Lambda<0> > > > > > forcePolicy;

typedef RAJA::ReduceSum<RAJA::omp_reduce, real_t> rajaReduceSumReal;
#endif

#ifdef DO_CUDA
#ifndef ENABLE_OPENMP
typedef RAJA::seq_segit linkCellTraversal;
typedef RAJA::ExecPolicy<RAJA::seq_segit, RAJA::simd_exec> linkCellWork;
typedef RAJA::ExecPolicy<RAJA::seq_segit, RAJA::simd_exec> atomWork;
typedef RAJA::seq_segit task_graph_policy;

typedef RAJA::KernelPolicy<
  RAJA::statement::For<0, RAJA::seq_exec,
  RAJA::statement::For<1, RAJA::seq_exec,
  RAJA::statement::For<2, RAJA::seq_exec,
  RAJA::statement::For<3, RAJA::simd_exec,
  RAJA::statement::Lambda<0> > > > > > forcePolicy;

typedef RAJA::ReduceSum<RAJA::seq_reduce, real_t> rajaReduceSumReal;
#endif

#define CUDA_ASYNC

typedef RAJA::KernelPolicy<
#ifdef CUDA_ASYNC
  RAJA::statement::CudaKernelAsync<
#else
  RAJA::statement::CudaKernel<
#endif
    RAJA::statement::For<0, RAJA::cuda_block_exec,
    RAJA::statement::For<1, RAJA::cuda_thread_exec,
    RAJA::statement::Lambda<0> > > > > atomWorkGPU;

typedef RAJA::KernelPolicy<
#ifdef CUDA_ASYNC
  RAJA::statement::CudaKernelAsync<
#else
  RAJA::statement::CudaKernel<
#endif
    RAJA::statement::For<0, RAJA::cuda_threadblock_exec<128>,
    RAJA::statement::Lambda<0> > > > redistributeGPU;

typedef RAJA::KernelPolicy<
#ifdef CUDA_ASYNC
  RAJA::statement::CudaKernelAsync<
#else
  RAJA::statement::CudaKernel<
#endif
    RAJA::statement::For<0, RAJA::cuda_block_exec,
    RAJA::statement::For<1, RAJA::cuda_thread_exec,
    RAJA::statement::Lambda<0> > > > > atomPackGPU;

typedef RAJA::KernelPolicy<
#ifdef CUDA_ASYNC
  RAJA::statement::CudaKernelAsync<
#else
  RAJA::statement::CudaKernel<
#endif
    RAJA::statement::For<0, RAJA::cuda_block_exec,
    RAJA::statement::For<1, RAJA::cuda_block_exec,
    RAJA::statement::For<2, RAJA::cuda_thread_exec,
    RAJA::statement::For<3, RAJA::cuda_thread_exec,
    RAJA::statement::Lambda<0> > > > > > > forcePolicyGPU;

typedef RAJA::ReduceSum<RAJA::cuda_reduce<27>, real_t> rajaReduceSumRealCUDA;
typedef RAJA::ReduceSum<RAJA::cuda_reduce<27>, int> rajaReduceSumIntCUDA;
#endif

#ifndef ENABLE_OPENMP
#ifndef DO_CUDA
typedef RAJA::seq_segit linkCellTraversal;
typedef RAJA::ExecPolicy<RAJA::seq_segit, RAJA::simd_exec> linkCellWork;
typedef RAJA::ExecPolicy<RAJA::seq_segit, RAJA::simd_exec> atomWork;
typedef RAJA::seq_segit task_graph_policy;

typedef RAJA::KernelPolicy<
  RAJA::statement::For<0, RAJA::seq_exec,
  RAJA::statement::For<1, RAJA::seq_exec,
  RAJA::statement::For<2, RAJA::seq_exec,
  RAJA::statement::For<3, RAJA::simd_exec,
  RAJA::statement::Lambda<0> > > > > > forcePolicy;

typedef RAJA::ReduceSum<RAJA::seq_reduce, real_t> rajaReduceSumReal;
#endif
#endif

#endif
