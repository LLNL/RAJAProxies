
#include "RAJA/RAJA.hpp"

//
//   RAJA IndexSet type used in loop traversals.
//
using LULESH_ISET = RAJA::TypedIndexSet<RAJA::RangeSegment,
                                        RAJA::ListSegment,
                                        RAJA::RangeStrideSegment>;

#include "luleshPolicy.hpp"
#ifdef RAJA_ENABLE_CHAI
#include "luleshMemory_chai.hpp"
#else
#include "luleshMemory.hpp"
#endif

#define MAX(a, b) ( ((a) > (b)) ? (a) : (b))

/* if luleshPolicy.hxx USE_CASE >= 9, must use lulesh_ptr.h */
#if USE_CASE >= LULESH_CUDA_CANONICAL
#if defined(LULESH_HEADER)
#undef LULESH_HEADER
#endif
#define LULESH_HEADER 1
#endif

#ifdef RAJA_ENABLE_CHAI
#include "lulesh_chai.hpp"
#else
//#if !defined(LULESH_HEADER)
//#include "lulesh_stl.hpp"
//#elif (LULESH_HEADER == 1)
#include "lulesh_ptr.hpp"
//#else
//#include "lulesh_tuple.hpp"
//#endif
#endif
