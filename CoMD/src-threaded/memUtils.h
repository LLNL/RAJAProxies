/// \file
/// Wrappers for memory allocation.

#ifndef _MEMUTILS_H_
#define _MEMUTILS_H_

#include <stdlib.h>
#include <string.h>

#define freeMe(s,element) {if(s->element) comdFree(s->element);  s->element = NULL;}

static void* comdMalloc(size_t iSize)
{
   void *localPtr ;
#ifdef DO_CUDA
   cudaMallocManaged(&localPtr, iSize);
#else
   posix_memalign(&localPtr, 64, iSize) ;
#endif
   return localPtr ;
}

static void* comdCalloc(size_t num, size_t iSize)
{
   void *localPtr ;
#ifdef DO_CUDA
   cudaMallocManaged(&localPtr, iSize);
   cudaMemset(localPtr, 0, num*iSize);
#else
   posix_memalign(&localPtr, 64, num*iSize) ;
   memset(localPtr, 0, num*iSize) ;
#endif
   return localPtr ;
}

static void comdFree(void *ptr)
{
#ifdef DO_CUDA
  cudaFree(ptr);
#else
  free(ptr);
#endif
}
#endif
