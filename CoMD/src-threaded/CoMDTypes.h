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

#ifdef ENABLE_OPENMP
#include <omp.h>
#endif

#ifdef ENABLE_CALIPER
#include <caliper/cali.h>
#endif

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

/* This is essentially a copy of the SimFlat structure that is passed throughout many of
 * the CoMD functions.  Only the constant values and variables only needed on the CPU will be
 * copied into this structure, none of the atom data or large arrays.
 *
 * The purpose of this is to avoid accessing the simulation structure on the CPU side after it has
 * been transferred to GPU memory.  This creates unnecessary page faults and drastically reduces
 * the performance of the code.
 */
//extern SimFlat *globalSim;
extern real_t ePotential;     //!< the total potential energy of the system
extern real_t eKinetic;       //!< the total kinetic energy of the system
#ifdef DO_CUDA
extern int nLocal;    //!< total number of atoms on this processor
extern int nGlobal;   //!< total number of atoms in simulation
// Allows for __device__ to be put before kernels only when CUDA is enabled
#define COMD_DEVICE __device__
#else
#define COMD_DEVICE
#endif

/*
########################################################################
########################## OpenMP Portion ##############################
########################################################################
*/

#ifdef ENABLE_OPENMP
typedef RAJA::omp_parallel_for_segit linkCellTraversal ;
typedef RAJA::ExecPolicy<RAJA::seq_segit, RAJA::simd_exec> linkCellWork;
typedef RAJA::ExecPolicy<RAJA::omp_parallel_for_segit, RAJA::simd_exec> atomWork;
typedef RAJA::omp_parallel_segit task_graph_policy;

typedef RAJA::KernelPolicy<
    RAJA::statement::For<0, RAJA::omp_parallel_for_segit,
    RAJA::statement::For<1, RAJA::simd_exec,
    RAJA::statement::Lambda<0> > > > atomWorkKernel;

typedef RAJA::KernelPolicy<
    RAJA::statement::For<0, RAJA::omp_parallel_for_segit,
    RAJA::statement::Lambda<0> > > redistributeKernel;

// Used for eam
typedef RAJA::KernelPolicy<
  RAJA::statement::For<0, RAJA::omp_parallel_for_segit,
  RAJA::statement::For<1, RAJA::seq_exec,
  RAJA::statement::For<2, RAJA::seq_exec,
  RAJA::statement::For<3, RAJA::simd_exec,
  RAJA::statement::Lambda<0> > > > > > eamforcePolicyKernel;

// Used for ljForce
typedef RAJA::KernelPolicy<
  RAJA::statement::For<0, RAJA::omp_parallel_for_segit,
  RAJA::statement::For<1, RAJA::seq_exec,
  RAJA::statement::For<2, RAJA::seq_exec,
  RAJA::statement::For<3, RAJA::simd_exec,
  RAJA::statement::Lambda<0> > > > > > forcePolicyKernel;

typedef RAJA::ReduceSum<RAJA::omp_reduce, real_t> rajaReduceSumReal;
typedef RAJA::ReduceSum<RAJA::seq_reduce, real_t> rajaReduceSumRealKernel;
typedef RAJA::ReduceSum<RAJA::seq_reduce, int> rajaReduceSumInt;
#endif

/*
########################################################################
########################### CUDA Portion ###############################
########################################################################
*/

#ifdef DO_CUDA
//#define CUDA_ASYNC

#ifndef ENABLE_OPENMP
typedef RAJA::seq_segit linkCellTraversal;
typedef RAJA::ExecPolicy<RAJA::seq_segit, RAJA::simd_exec> linkCellWork;
typedef RAJA::ExecPolicy<RAJA::seq_segit, RAJA::simd_exec> atomWork;
typedef RAJA::seq_segit task_graph_policy;
#endif

typedef RAJA::KernelPolicy<
#ifdef CUDA_ASYNC
  RAJA::statement::CudaKernelAsync<
#else
  RAJA::statement::CudaKernel<
#endif
    RAJA::statement::For<0, RAJA::cuda_block_x_loop,
    RAJA::statement::For<1, RAJA::cuda_thread_x_loop,
    RAJA::statement::Lambda<0> > > > > atomWorkKernel;

typedef RAJA::KernelPolicy<
#ifdef CUDA_ASYNC
  RAJA::statement::CudaKernelAsync<
#else
  RAJA::statement::CudaKernel<
#endif
    RAJA::statement::Tile<0, RAJA::statement::tile_fixed<128>, RAJA::cuda_block_x_loop,
    RAJA::statement::For<0, RAJA::cuda_thread_x_loop,
    RAJA::statement::Lambda<0> > > > > redistributeKernel;

// Used for ljForce
typedef RAJA::KernelPolicy<
#ifdef CUDA_ASYNC
  RAJA::statement::CudaKernelAsync<
#else
  RAJA::statement::CudaKernel<
#endif
    RAJA::statement::For<0, RAJA::cuda_block_x_loop,
    RAJA::statement::For<1, RAJA::cuda_block_y_loop,
    RAJA::statement::For<2, RAJA::cuda_thread_x_loop,
    RAJA::statement::For<3, RAJA::cuda_thread_y_loop,
    RAJA::statement::Lambda<0> > > > > > > forcePolicyKernel;

typedef RAJA::ReduceSum<RAJA::seq_reduce, real_t> rajaReduceSumReal;
typedef RAJA::ReduceSum<RAJA::cuda_reduce, real_t> rajaReduceSumRealKernel;
typedef RAJA::ReduceSum<RAJA::cuda_reduce, int> rajaReduceSumInt;
#endif

/*
########################################################################
########################## Serial Portion ##############################
########################################################################
*/

#ifndef ENABLE_OPENMP
#ifndef DO_CUDA
typedef RAJA::seq_segit linkCellTraversal;
typedef RAJA::ExecPolicy<RAJA::seq_segit, RAJA::simd_exec> linkCellWork;
typedef RAJA::ExecPolicy<RAJA::seq_segit, RAJA::simd_exec> atomWork;
typedef RAJA::seq_segit task_graph_policy;

typedef RAJA::KernelPolicy<
    RAJA::statement::For<0, RAJA::seq_segit,
    RAJA::statement::For<1, RAJA::simd_exec,
    RAJA::statement::Lambda<0> > > > atomWorkKernel;

typedef RAJA::KernelPolicy<
    RAJA::statement::For<0, RAJA::seq_exec,
    RAJA::statement::Lambda<0> > > redistributeKernel;

// Used for eam
typedef RAJA::KernelPolicy<
  RAJA::statement::For<0, RAJA::seq_exec,
  RAJA::statement::For<1, RAJA::seq_exec,
  RAJA::statement::For<2, RAJA::seq_exec,
  RAJA::statement::For<3, RAJA::simd_exec,
  RAJA::statement::Lambda<0> > > > > > eamforcePolicyKernel;

// Used for ljForce
typedef RAJA::KernelPolicy<
  RAJA::statement::For<0, RAJA::seq_exec,
  RAJA::statement::For<1, RAJA::seq_exec,
  RAJA::statement::For<2, RAJA::seq_exec,
  RAJA::statement::For<3, RAJA::simd_exec,
  RAJA::statement::Lambda<0> > > > > > forcePolicyKernel;

// TODO Remove these...
/* Temporary Sequential Versions */
typedef RAJA::KernelPolicy<
    RAJA::statement::For<0, RAJA::seq_segit,
    RAJA::statement::For<1, RAJA::simd_exec,
    RAJA::statement::Lambda<0> > > > atomWorkKernelSeq;

// Used for eam
typedef RAJA::KernelPolicy<
  RAJA::statement::For<0, RAJA::seq_exec,
  RAJA::statement::For<1, RAJA::seq_exec,
  RAJA::statement::For<2, RAJA::seq_exec,
  RAJA::statement::For<3, RAJA::simd_exec,
  RAJA::statement::Lambda<0> > > > > > eamforcePolicyKernelSeq;

typedef RAJA::ReduceSum<RAJA::seq_reduce, real_t> rajaReduceSumReal;
typedef RAJA::ReduceSum<RAJA::seq_reduce, real_t> rajaReduceSumRealKernel;
typedef RAJA::ReduceSum<RAJA::seq_reduce, int> rajaReduceSumInt;
#endif
#endif

#endif
