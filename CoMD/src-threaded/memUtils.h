/// \file
/// Wrappers for memory allocation.

#ifndef _MEMUTILS_H_
#define _MEMUTILS_H_

#include <stdlib.h>
#include <string.h>

#if defined(ENABLE_HIP)
#include <hip/hip_runtime.h>
#endif

#define freeMe(s,element) {if(s->element) comdFree(s->element);  s->element = NULL;}

static void* comdMalloc(size_t num, size_t iSize)
{
   void *localPtr ;
#ifdef ENABLE_CUDA
   cudaMallocManaged(&localPtr, num*iSize);
#elif defined(ENABLE_HIP)
   hipMalloc(&localPtr, num*iSize);
#else
   posix_memalign(&localPtr, 64, num*iSize) ;
#endif
   return localPtr ;
}

static void* comdCalloc(size_t num, size_t iSize)
{
   void *localPtr ;
#ifdef ENABLE_CUDA
   cudaMallocManaged(&localPtr, num*iSize);
   cudaMemset(localPtr, 0, num*iSize);
#elif defined(ENABLE_HIP)
   hipMalloc(&localPtr, num*iSize);
   hipMemset(localPtr, 0, num*iSize);
#else
   posix_memalign(&localPtr, 64, num*iSize) ;
   memset(localPtr, 0, num*iSize) ;
#endif
   return localPtr ;
}

static void comdFree(void *ptr)
{
#ifdef ENABLE_CUDA
  cudaFree(ptr);
#elif defined(ENABLE_HIP)
  hipFree(ptr);
#else
  free(ptr);
#endif
}
#endif
