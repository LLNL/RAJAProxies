/// \file
/// Frequently needed typedefs.

#ifndef __MYTYPE_H_
#define __MYTYPE_H_

/// \def SINGLE determines whether single or double precision is built
// NOTE: Switching to SINGLE will require some code changes at this time.
//       Usage of floor() in CUDA kernels will need to become the float
//       version, and the MPI functionality will need to use MPI_FLOAT.
#ifdef SINGLE
typedef float real_t;  //!< define native type for CoMD as single precision
  #define FMT1 "%g"    //!< /def format argument for floats 
  #define EMT1 "%e"    //!< /def format argument for eng floats
#else
typedef double real_t; //!< define native type for CoMD as double precision
  #define FMT1 "%lg"   //!< \def format argument for doubles 
  #define EMT1 "%le"   //!< \def format argument for eng doubles 
#endif
typedef real_t * __restrict__ /* __attribute__((align_value(64))) */ real_ptr ;

typedef real_t real3[3]; //!< a convenience vector with three real_t 
typedef real3 * __restrict__ /* __attribute__((align_value(64))) */ real3_ptr ;

#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
#define COMD_HOST_DEVICE __host__ __device__
#else
#define COMD_HOST_DEVICE
#endif

static void zeroReal3(real3 a)
{
   a[0] = 0.0;
   a[1] = 0.0;
   a[2] = 0.0;
}

#define screenOut stdout

#endif
