// This work was performed under the auspices of the U.S. Department of Energy by
// Lawrence Livermore National Laboratory under Contract DE-AC52-07NA27344.
//

//
// ALLOCATE/RELEASE FUNCTIONS 
//
/**********************************/
/* Memory Pool                    */
/**********************************/

namespace lulesh2 {

template <typename VARTYPE >
struct MemoryPool {
public:
   typedef chai::ManagedArray<VARTYPE>  POOLTYPE;
   typedef POOLTYPE&                    RETTYPE;

   MemoryPool()
   {
      for (int i=0; i<32; ++i) {
         lenType[i] = 0 ;
      }
   }

   RETTYPE allocate(int len) {
      int i = 0;
      for (i=0; i<32; ++i) {
         if (lenType[i] == len) {
            lenType[i] = -lenType[i] ;
            break ;
         }
         else if (lenType[i] == 0) {
            lenType[i] = -len ;
#ifdef RAJA_ENABLE_CUDA
            ptr[i].allocate(len, chai::GPU) ;
#else
            ptr[i].allocate(len, chai::CPU) ;
#endif
            break ;
         }
      }
      return ptr[i];
   }

   bool release(const RETTYPE oldPtr) {
      int i ;
      bool success = true ;
      for (i=0; i<32; ++i) {
         if (ptr[i] == oldPtr) {
            lenType[i] = -lenType[i] ;
            break ;
         }
      }
      if (i == 32) {
         success = false ; /* error -- not found */
      }
      return success ;
   }

   POOLTYPE ptr[32] ; 
   int lenType[32] ;
} ;

}

