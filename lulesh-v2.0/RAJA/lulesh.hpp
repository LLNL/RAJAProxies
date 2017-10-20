
#include "RAJA/RAJA.hpp"
#include "luleshPolicy.hpp"
#include "luleshMemory.hpp"

#define MAX(a, b) ( ((a) > (b)) ? (a) : (b))

/* if luleshPolicy.hxx USE_CASE >= 9, must use lulesh_ptr.h */
#if USE_CASE >= LULESH_CUDA_CANONICAL
#if defined(LULESH_HEADER)
#undef LULESH_HEADER
#endif
#define LULESH_HEADER 1
#endif

#if !defined(LULESH_HEADER)
#include "lulesh_stl.hpp"
#elif (LULESH_HEADER == 1)
#include "lulesh_ptr.h"
#else
#include "lulesh_tuple.h"
#endif
