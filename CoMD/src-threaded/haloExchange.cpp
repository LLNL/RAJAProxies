/// \file
/// Communicate halo data such as "ghost" atoms with neighboring tasks.
/// In addition to ghost atoms, the EAM potential also needs to exchange
/// some force information.  Hence this file implements both an atom
/// exchange and a force exchange, each with slightly different
/// properties due to their different roles.
///
/// The halo exchange in CoMD 1.1 takes advantage of the Cartesian domain
/// decomposition as well as the link cell structure to quickly
/// determine what data needs to be sent.
///
/// This halo exchange implementation is able to send data to all 26
/// neighboring tasks using only 6 messages.  This is accomplished by
/// sending data across the x-faces, then the y-faces, and finally
/// across the z-faces.  Some of the data that was received from the
/// x-faces is included in the y-face sends and so on.  This
/// accumulation of data allows data to reach edge neighbors and corner
/// neighbors by a two or three step process.
///
/// The advantage of this type of structured halo exchange is that it
/// minimizes the number of MPI messages to send, and maximizes the size
/// of those messages.
///
/// The disadvantage of this halo exchange is that it serializes message
/// traffic.  Only two messages can be in flight at once. The x-axis
/// messages must be received and processed before the y-axis messages
/// can begin.  Architectures with low message latency and many off node
/// network links would likely benefit from alternate halo exchange
/// strategies that send independent messages to each neighbor task.

#include "haloExchange.h"

#include <assert.h>

#include "CoMDTypes.h"
#include "decomposition.h"
#include "parallel.h"
#include "linkCells.h"
#include "eam.h"
#include "memUtils.h"
#include "performanceTimers.h"

#define MAX(A,B) ((A) > (B) ? (A) : (B))

/// Don't change the order of the faces in this enum.
enum HaloFaceOrder {HALO_X_MINUS, HALO_X_PLUS,
                    HALO_Y_MINUS, HALO_Y_PLUS,
                    HALO_Z_MINUS, HALO_Z_PLUS};

/// Don't change the order of the axes in this enum.
enum HaloAxisOrder {HALO_X_AXIS, HALO_Y_AXIS, HALO_Z_AXIS};

/// Extra data members that are needed for the exchange of atom data.
/// For an atom exchange, the HaloExchangeSt::parms will point to a
/// structure of this type.
typedef struct AtomExchangeParmsSt
{
   int nCells[6];        //!< Number of cells in cellList for each face.
   int* cellList[6];     //!< List of link cells from which to load data for each face.
   real_t* pbcFactor[6]; //!< Whether this face is a periodic boundary.
   char* sendBufM;
   char* sendBufP;
   char* recvBufM;
   char* recvBufP;
}
AtomExchangeParms;

/// Extra data members that are needed for the exchange of force data.
/// For an force exchange, the HaloExchangeSt::parms will point to a
/// structure of this type.
typedef struct ForceExchangeParmsSt
{
   int nCells[6];     //!< Number of cells to send/recv for each face.
   int* sendCells[6]; //!< List of link cells to send for each face.
   int* recvCells[6]; //!< List of link cells to recv for each face.
   char* sendBufM;
   char* sendBufP;
   char* recvBufM;
   char* recvBufP;
}
ForceExchangeParms;

/// Package data for the force exchange.
typedef struct ForceMsgSt
{
   real_t dfEmbed;
}
ForceMsg;

static HaloExchange* initHaloExchange(Domain* domain);
static void exchangeData(HaloExchange* haloExchange, void* data, int iAxis, int source);

static int* mkAtomCellList(LinkCell* boxes, enum HaloFaceOrder iFace, const int nCells);
static int loadAtomsBuffer(void* vparms, void* data, int face, char* charBuf);
// This has been moved to linkCells.cpp
extern void unloadAtomsBuffer(void* vparms, void* data, int face, int bufSize, char* charBuf);
static void destroyAtomsExchange(void* vparms);

static int* mkForceSendCellList(LinkCell* boxes, int face, int nCells);
static int* mkForceRecvCellList(LinkCell* boxes, int face, int nCells);
static int loadForceBuffer(void* vparms, void* data, int face, char* charBuf);
static void unloadForceBuffer(void* vparms, void* data, int face, int bufSize, char* charBuf);
static void destroyForceExchange(void* vparms);
extern "C" int sortAtomsById(const void* a, const void* b);

/// \details
/// When called in proper sequence by redistributeAtoms, the atom halo
/// exchange helps serve three purposes:
/// - Send ghost atom data to neighbor tasks.
/// - Shift atom coordinates by the global simulation size when they cross
///   periodic boundaries.  This shift is performed in loadAtomsBuffer.
/// - Transfer ownership of atoms between tasks as the atoms move across
///   spatial domain boundaries.  This transfer of ownership occurs in
///   two places.  The former owner gives up ownership when
///   updateLinkCells moves a formerly local atom into a halo link cell.
///   The new owner accepts ownership when unloadAtomsBuffer calls
///   putAtomInBox to place a received atom into a local link cell.
///
/// This constructor does the following:
///
/// - Sets the bufCapacity to hold the largest possible number of atoms
///   that can be sent across a face.
/// - Initialize function pointers to the atom-specific versions
/// - Sets the number of link cells to send across each face.
/// - Builds the list of link cells to send across each face.  As
///   explained in the comments for mkAtomCellList, this list must
///   include any link cell, local or halo, that could possibly contain
///   an atom that needs to be sent across the face.  Atoms that need to
///   be sent include "ghost atoms" that are located in local link
///   cells that correspond to halo link cells on receiving tasks as well as
///   formerly local atoms that have just moved into halo link cells and
///   need to be sent to the rank that owns the spatial domain the atom
///   has moved into.
/// - Sets a coordinate shift factor for each face to account for
///   periodic boundary conditions.  For most faces the factor is zero.
///   For faces on the +x, +y, or +z face of the simulation domain
///   the factor is -1.0 (to shift the coordinates by -1 times the
///   simulation domain size).  For -x, -y, and -z faces of the
///   simulation domain, the factor is +1.0.
///
/// \see redistributeAtoms
HaloExchange* initAtomHaloExchange(Domain* domain, LinkCell* boxes)
{
   HaloExchange* hh = initHaloExchange(domain);

   int size0 = (boxes->gridSize[1]+2)*(boxes->gridSize[2]+2);
   int size1 = (boxes->gridSize[0]+2)*(boxes->gridSize[2]+2);
   int size2 = (boxes->gridSize[0]+2)*(boxes->gridSize[1]+2);
   int maxSize = MAX(size0, size1);
   maxSize = MAX(size1, size2); // Shouldn't this be MAX(maxSize, size2)?

   hh->bufCapacity = maxSize*2*MAXATOMS*sizeof(AtomMsg);

   hh->loadBuffer = loadAtomsBuffer;
   hh->unloadBuffer = unloadAtomsBuffer;
   hh->destroy = destroyAtomsExchange;

   AtomExchangeParms* parms = (AtomExchangeParms*) comdMalloc(1, sizeof(AtomExchangeParms));

   parms->nCells[HALO_X_MINUS] = 2*(boxes->gridSize[1]+2)*(boxes->gridSize[2]+2);
   parms->nCells[HALO_Y_MINUS] = 2*(boxes->gridSize[0]+2)*(boxes->gridSize[2]+2);
   parms->nCells[HALO_Z_MINUS] = 2*(boxes->gridSize[0]+2)*(boxes->gridSize[1]+2);
   parms->nCells[HALO_X_PLUS]  = parms->nCells[HALO_X_MINUS];
   parms->nCells[HALO_Y_PLUS]  = parms->nCells[HALO_Y_MINUS];
   parms->nCells[HALO_Z_PLUS]  = parms->nCells[HALO_Z_MINUS];

   /* Make sure there is enough space in the buffer for the integer array of
    * offsets to allow for parallel packing/unpacking.
    */
   maxSize = parms->nCells[HALO_X_MINUS];
   if(parms->nCells[HALO_Y_MINUS] > maxSize) maxSize = parms->nCells[HALO_Y_MINUS];
   if(parms->nCells[HALO_Z_MINUS] > maxSize) maxSize = parms->nCells[HALO_Z_MINUS];
   if(parms->nCells[HALO_X_PLUS] > maxSize) maxSize = parms->nCells[HALO_X_PLUS];
   if(parms->nCells[HALO_Y_PLUS] > maxSize) maxSize = parms->nCells[HALO_Y_PLUS];
   if(parms->nCells[HALO_Z_PLUS] > maxSize) maxSize = parms->nCells[HALO_Z_PLUS];

   int intsSize = (maxSize+1) * sizeof(int);
   if( (intsSize % sizeof(AtomMsg)) != 0 )
     intsSize += sizeof(AtomMsg) - (intsSize % sizeof(AtomMsg));

   hh->bufCapacity += intsSize;

   for (int ii=0; ii<6; ++ii)
      parms->cellList[ii] = mkAtomCellList(boxes, HaloFaceOrder(ii), parms->nCells[ii]);

   for (int ii=0; ii<6; ++ii)
   {
      parms->pbcFactor[ii] = (real_t *)comdMalloc(3, sizeof(real_t));
      for (int jj=0; jj<3; ++jj)
         parms->pbcFactor[ii][jj] = 0.0;
   }
   int* procCoord = domain->procCoord; //alias
   int* procGrid  = domain->procGrid; //alias
   if (procCoord[HALO_X_AXIS] == 0)                       parms->pbcFactor[HALO_X_MINUS][HALO_X_AXIS] = +1.0;
   if (procCoord[HALO_X_AXIS] == procGrid[HALO_X_AXIS]-1) parms->pbcFactor[HALO_X_PLUS][HALO_X_AXIS]  = -1.0;
   if (procCoord[HALO_Y_AXIS] == 0)                       parms->pbcFactor[HALO_Y_MINUS][HALO_Y_AXIS] = +1.0;
   if (procCoord[HALO_Y_AXIS] == procGrid[HALO_Y_AXIS]-1) parms->pbcFactor[HALO_Y_PLUS][HALO_Y_AXIS]  = -1.0;
   if (procCoord[HALO_Z_AXIS] == 0)                       parms->pbcFactor[HALO_Z_MINUS][HALO_Z_AXIS] = +1.0;
   if (procCoord[HALO_Z_AXIS] == procGrid[HALO_Z_AXIS]-1) parms->pbcFactor[HALO_Z_PLUS][HALO_Z_AXIS]  = -1.0;

   parms->sendBufM = (char *) comdMalloc(hh->bufCapacity, sizeof(char));
   parms->sendBufP = (char *) comdMalloc(hh->bufCapacity, sizeof(char));
   parms->recvBufM = (char *) comdMalloc(hh->bufCapacity, sizeof(char));
   parms->recvBufP = (char *) comdMalloc(hh->bufCapacity, sizeof(char));

   hh->parms = parms;
   return hh;
}

/// The force exchange is considerably simpler than the atom exchange.
/// In the force case we only need to exchange data that is needed to
/// complete the force calculation.  Since the atoms have not moved we
/// only need to send data from local link cells and we are guaranteed
/// that the same atoms exist in the same order in corresponding halo
/// cells on remote tasks.  The only tricky part is the size of the
/// plane of local cells that needs to be sent grows in each direction.
/// This is because the y-axis send must send some of the data that was
/// received from the x-axis send, and the z-axis must send some data
/// from the y-axis send.  This accumulation of data to send is
/// responsible for data reaching neighbor cells that share only edges
/// or corners.
///
/// \see eam.c for an explanation of the requirement to exchange
/// force data.
HaloExchange* initForceHaloExchange(Domain* domain, LinkCell* boxes)
{
   HaloExchange* hh = initHaloExchange(domain);

   hh->loadBuffer = loadForceBuffer;
   hh->unloadBuffer = unloadForceBuffer;
   hh->destroy = destroyForceExchange;

   int size0 = (boxes->gridSize[1])*(boxes->gridSize[2]);
   int size1 = (boxes->gridSize[0]+2)*(boxes->gridSize[2]);
   int size2 = (boxes->gridSize[0]+2)*(boxes->gridSize[1]+2);
   int maxSize = MAX(size0, size1);
   maxSize = MAX(size1, size2);
   hh->bufCapacity = (maxSize)*MAXATOMS*sizeof(ForceMsg);

   ForceExchangeParms* parms = (ForceExchangeParms*) comdMalloc(1, sizeof(ForceExchangeParms));

   parms->nCells[HALO_X_MINUS] = (boxes->gridSize[1]  )*(boxes->gridSize[2]  );
   parms->nCells[HALO_Y_MINUS] = (boxes->gridSize[0]+2)*(boxes->gridSize[2]  );
   parms->nCells[HALO_Z_MINUS] = (boxes->gridSize[0]+2)*(boxes->gridSize[1]+2);
   parms->nCells[HALO_X_PLUS]  = parms->nCells[HALO_X_MINUS];
   parms->nCells[HALO_Y_PLUS]  = parms->nCells[HALO_Y_MINUS];
   parms->nCells[HALO_Z_PLUS]  = parms->nCells[HALO_Z_MINUS];

   for (int ii=0; ii<6; ++ii)
   {
      parms->sendCells[ii] = mkForceSendCellList(boxes, ii, parms->nCells[ii]);
      parms->recvCells[ii] = mkForceRecvCellList(boxes, ii, parms->nCells[ii]);
   }

   parms->sendBufM = (char *) comdMalloc(hh->bufCapacity, sizeof(char));
   parms->sendBufP = (char *) comdMalloc(hh->bufCapacity, sizeof(char));
   parms->recvBufM = (char *) comdMalloc(hh->bufCapacity, sizeof(char));
   parms->recvBufP = (char *) comdMalloc(hh->bufCapacity, sizeof(char));
   hh->parms = parms;

   return hh;
}

void destroyHaloExchange(HaloExchange** haloExchange)
{
   (*haloExchange)->destroy((*haloExchange)->parms);
   comdFree((*haloExchange)->parms);
   comdFree(*haloExchange);
   *haloExchange = NULL;
}

// haloExchange is used by both atom and force
// the extra argument source is to indicate which one 
// invokes this function so that we can cast the parms to 
// according types
void haloExchange(HaloExchange* haloExchange, void* data, int source)
{
   for (int iAxis=0; iAxis<3; ++iAxis) {
      exchangeData(haloExchange, data, iAxis, source);
   }
}

/// Base class constructor.
HaloExchange* initHaloExchange(Domain* domain)
{
   HaloExchange* hh = (HaloExchange*) comdMalloc(1, sizeof(HaloExchange));

   // Rank of neighbor task for each face.
   hh->nbrRank[HALO_X_MINUS] = processorNum(domain, -1,  0,  0);
   hh->nbrRank[HALO_X_PLUS]  = processorNum(domain, +1,  0,  0);
   hh->nbrRank[HALO_Y_MINUS] = processorNum(domain,  0, -1,  0);
   hh->nbrRank[HALO_Y_PLUS]  = processorNum(domain,  0, +1,  0);
   hh->nbrRank[HALO_Z_MINUS] = processorNum(domain,  0,  0, -1);
   hh->nbrRank[HALO_Z_PLUS]  = processorNum(domain,  0,  0, +1);
   hh->bufCapacity = 0; // will be set by sub-class.

   return hh;
}

/// This is the function that does the heavy lifting for the
/// communication of halo data.  It is called once for each axis and
/// sends and receives two message.  Loading and unloading of the
/// buffers is in the hands of the sub-class virtual functions.
///
/// \param [in] iAxis     Axis index.
/// \param [in, out] data Pointer to data that will be passed to the load and
///                       unload functions

// source - 0: exchange atom data
//        - 1: exchange force data
void exchangeData(HaloExchange* haloExchange, void* data, int iAxis, int source)
{
   enum HaloFaceOrder faceM = HaloFaceOrder(2*iAxis);
   enum HaloFaceOrder faceP = HaloFaceOrder(faceM+1);

   char* sendBufM;
   char* sendBufP;
   char* recvBufM;
   char* recvBufP;
   if (source == 0) {
      AtomExchangeParms* parms = (AtomExchangeParms*) haloExchange->parms;
      sendBufM = parms->sendBufM;
      sendBufP = parms->sendBufP;
      recvBufM = parms->recvBufM;
      recvBufP = parms->recvBufP; 
   } else {
      ForceExchangeParms* parms = (ForceExchangeParms*) haloExchange->parms;
      sendBufM = parms->sendBufM;
      sendBufP = parms->sendBufP;
      recvBufM = parms->recvBufM;
      recvBufP = parms->recvBufP; 
   }

   startTimer(atomPackTimer);
   int nSendM = haloExchange->loadBuffer(haloExchange->parms, data, faceM, sendBufM);
   int nSendP = haloExchange->loadBuffer(haloExchange->parms, data, faceP, sendBufP);
   stopTimer(atomPackTimer);

   int nbrRankM = haloExchange->nbrRank[faceM];
   int nbrRankP = haloExchange->nbrRank[faceP];

   int nRecvM, nRecvP;

   startTimer(atomCommTimer);
   nRecvP = sendReceiveParallel(sendBufM, nSendM, nbrRankM, recvBufP, haloExchange->bufCapacity, nbrRankP);
   nRecvM = sendReceiveParallel(sendBufP, nSendP, nbrRankP, recvBufM, haloExchange->bufCapacity, nbrRankM);
   stopTimer(atomCommTimer);

   startTimer(atomUnpackTimer);
   haloExchange->unloadBuffer(haloExchange->parms, data, faceM, nRecvM, recvBufM);
   haloExchange->unloadBuffer(haloExchange->parms, data, faceP, nRecvP, recvBufP);
   stopTimer(atomUnpackTimer);
}

/// Make a list of link cells that need to be sent across the specified
/// face.  For each face, the list must include all cells, local and
/// halo, in the first two planes of link cells.  Halo cells must be
/// included in the list of link cells to send since local atoms may
/// have moved from local cells into halo cells on this time step.
/// (Actual remote atoms should have been deleted, so the halo cells
/// should contain only these few atoms that have just crossed.)
/// Sending these atoms will allow them to be reassigned to the task
/// that covers the spatial domain they have moved into.
///
/// Note that link cell grid coordinates range from -1 to gridSize[iAxis].
/// \see initLinkCells for an explanation link cell grid coordinates.
///
/// \param [in] boxes  Link cell information.
/// \param [in] iFace  Index of the face data will be sent across.
/// \param [in] nCells Number of cells to send.  This is used for a
///                    consistency check.
/// \return The list of cells to send.  Caller is responsible to free
/// the list.
int* mkAtomCellList(LinkCell* boxes, enum HaloFaceOrder iFace, const int nCells)
{
   int* list = (int *) comdMalloc(nCells, sizeof(int));
   int xBegin = -1;
   int xEnd   = boxes->gridSize[0]+1;
   int yBegin = -1;
   int yEnd   = boxes->gridSize[1]+1;
   int zBegin = -1;
   int zEnd   = boxes->gridSize[2]+1;

   if (iFace == HALO_X_MINUS) xEnd = xBegin+2;
   if (iFace == HALO_X_PLUS)  xBegin = xEnd-2;
   if (iFace == HALO_Y_MINUS) yEnd = yBegin+2;
   if (iFace == HALO_Y_PLUS)  yBegin = yEnd-2;
   if (iFace == HALO_Z_MINUS) zEnd = zBegin+2;
   if (iFace == HALO_Z_PLUS)  zBegin = zEnd-2;

   int count = 0;
   for (int ix=xBegin; ix<xEnd; ++ix)
      for (int iy=yBegin; iy<yEnd; ++iy)
         for (int iz=zBegin; iz<zEnd; ++iz)
            list[count++] = getBoxFromTuple(boxes, ix, iy, iz);
   assert(count == nCells);
   return list;
}

/// The loadBuffer function for a halo exchange of atom data.  Iterates
/// link cells in the cellList and load any atoms into the send buffer.
/// This function also shifts coordinates of the atoms by an appropriate
/// factor if they are being sent across a periodic boundary.
///
/// \see HaloExchangeSt::loadBuffer for an explanation of the loadBuffer
/// parameters.
int loadAtomsBuffer(void* vparms, void* data, int face, char* charBuf)
{
   rajaReduceSumInt nBufReduce(0);

   AtomExchangeParms* parms = (AtomExchangeParms*) vparms;
   const int nCells = parms->nCells[face];
   int* cellList = parms->cellList[face];
   real_t* pbcFactor = parms->pbcFactor[face];

   int* bufferOffset = (int*)comdMalloc(nCells, sizeof(int));
   int size = (nCells+1) * sizeof(int);
   //Make sure the memory is alligned with AtomMsg size
#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
   if((size % sizeof(AtomMsg)) != 0)
     size += sizeof(AtomMsg) - (size % sizeof(AtomMsg));
#else
   // This is more efficient for the CPU-side versions since this work
   // only needs to be done once.  The CUDA/HIP versions need this to be on
   // the GPU to avoid page faults.
   SimFlat* s = (SimFlat*) data;
   int* offsetBuf = (int*)charBuf;
   real3 shift;
   shift[0] = pbcFactor[0] * s->domain->globalExtent[0];
   shift[1] = pbcFactor[1] * s->domain->globalExtent[1];
   shift[2] = pbcFactor[2] * s->domain->globalExtent[2];

   int sum = 0;
   offsetBuf[0] = nCells;
   for(int i = 0; i < nCells; i++) {
     bufferOffset[i] = sum;
     offsetBuf[i+1] = sum;
     sum += s->boxes->nAtoms[cellList[i]];
   }
#endif

   AtomMsg* buf = (AtomMsg*) (charBuf + size);
   // TODO: Try to optimize this.
   RAJA::kernel<redistributeKernel>(
     RAJA::make_tuple(
     RAJA::TypedRangeSegment<int>(0, nCells)),
   [=] RAJA_HOST_DEVICE (int iCell) {
       SimFlat* s = (SimFlat*) data;
       int* offsetBuf = (int*)charBuf;
       /* All accesses to the simulation data structure should be on the GPU to avoid page faults */
#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
       real3 shift;
       shift[0] = pbcFactor[0] * s->domain->globalExtent[0];
       shift[1] = pbcFactor[1] * s->domain->globalExtent[1];
       shift[2] = pbcFactor[2] * s->domain->globalExtent[2];

       int sum = 0;
       offsetBuf[0] = nCells;
       for(int i = 0; i < nCells; i++) {
         bufferOffset[i] = sum;
         offsetBuf[i+1] = sum;
         sum += s->boxes->nAtoms[cellList[i]];
       }
#endif

       int* gid = s->atoms->gid;
       int* iSpecies = s->atoms->iSpecies;
       real3_ptr r = s->atoms->r;
       real3_ptr p = s->atoms->p;

       int offset = bufferOffset[iCell];

       int iBox = cellList[iCell];
       int iOff = iBox*MAXATOMS;

       int ii = iOff;
       while(ii < iOff+s->boxes->nAtoms[iBox]) {
         buf[offset].gid  = gid[ii];
         buf[offset].type = iSpecies[ii];
         buf[offset].rx = r[ii][0] + shift[0];
         buf[offset].ry = r[ii][1] + shift[1];
         buf[offset].rz = r[ii][2] + shift[2];
         buf[offset].px = p[ii][0];
         buf[offset].py = p[ii][1];
         buf[offset].pz = p[ii][2];
         nBufReduce += 1;
         ii++;
         offset++;
       }
     });

   const int nBuf = (int)nBufReduce;

   comdFree(bufferOffset);

   return (nBuf*sizeof(AtomMsg)) + size;
}

/* ############################################################
   ################## unloadAtomsBuffer Moved #################
   ############################################################
   The unloadAtomsBuffer function has been moved to linkCells.cpp
   to avoid issues requiring separate compilation in the CUDA case.
 */

void destroyAtomsExchange(void* vparms)
{
   AtomExchangeParms* parms = (AtomExchangeParms*) vparms;

   for (int ii=0; ii<6; ++ii)
   {
      comdFree(parms->pbcFactor[ii]);
      comdFree(parms->cellList[ii]);
   }
   comdFree(parms->sendBufM);
   comdFree(parms->sendBufP);
   comdFree(parms->recvBufM);
   comdFree(parms->recvBufP);
}

/// Make a list of link cells that need to send data across the
/// specified face.  Note that this list must be compatible with the
/// corresponding recv list to ensure that the data goes to the correct
/// atoms.
///
/// \see initLinkCells for information about the conventions for grid
/// coordinates of link cells.
int* mkForceSendCellList(LinkCell* boxes, int face, int nCells)
{
   int* list = (int *)comdMalloc(nCells, sizeof(int));
   int xBegin, xEnd, yBegin, yEnd, zBegin, zEnd;

   int nx = boxes->gridSize[0];
   int ny = boxes->gridSize[1];
   int nz = boxes->gridSize[2];
   switch(face)
   {
     case HALO_X_MINUS:
      xBegin=0;    xEnd=1;    yBegin=0;    yEnd=ny;   zBegin=0;    zEnd=nz;
      break;
     case HALO_X_PLUS:
      xBegin=nx-1; xEnd=nx;   yBegin=0;    yEnd=ny;   zBegin=0;    zEnd=nz;
      break;
     case HALO_Y_MINUS:
      xBegin=-1;   xEnd=nx+1; yBegin=0;    yEnd=1;    zBegin=0;    zEnd=nz;
      break;
     case HALO_Y_PLUS:
      xBegin=-1;   xEnd=nx+1; yBegin=ny-1; yEnd=ny;   zBegin=0;    zEnd=nz;
      break;
     case HALO_Z_MINUS:
      xBegin=-1;   xEnd=nx+1; yBegin=-1;   yEnd=ny+1; zBegin=0;    zEnd=1;
      break;
     case HALO_Z_PLUS:
      xBegin=-1;   xEnd=nx+1; yBegin=-1;   yEnd=ny+1; zBegin=nz-1; zEnd=nz;
      break;
     default:
      assert(1==0);
   }

   int count = 0;
   for (int ix=xBegin; ix<xEnd; ++ix)
      for (int iy=yBegin; iy<yEnd; ++iy)
         for (int iz=zBegin; iz<zEnd; ++iz)
            list[count++] = getBoxFromTuple(boxes, ix, iy, iz);

   assert(count == nCells);
   return list;
}

/// Make a list of link cells that need to receive data across the
/// specified face.  Note that this list must be compatible with the
/// corresponding send list to ensure that the data goes to the correct
/// atoms.
///
/// \see initLinkCells for information about the conventions for grid
/// coordinates of link cells.
int* mkForceRecvCellList(LinkCell* boxes, int face, int nCells)
{
   int* list = (int *)comdMalloc(nCells, sizeof(int));
   int xBegin, xEnd, yBegin, yEnd, zBegin, zEnd;

   int nx = boxes->gridSize[0];
   int ny = boxes->gridSize[1];
   int nz = boxes->gridSize[2];
   switch(face)
   {
     case HALO_X_MINUS:
      xBegin=-1; xEnd=0;    yBegin=0;  yEnd=ny;   zBegin=0;  zEnd=nz;
      break;
     case HALO_X_PLUS:
      xBegin=nx; xEnd=nx+1; yBegin=0;  yEnd=ny;   zBegin=0;  zEnd=nz;
      break;
     case HALO_Y_MINUS:
      xBegin=-1; xEnd=nx+1; yBegin=-1; yEnd=0;    zBegin=0;  zEnd=nz;
      break;
     case HALO_Y_PLUS:
      xBegin=-1; xEnd=nx+1; yBegin=ny; yEnd=ny+1; zBegin=0;  zEnd=nz;
      break;
     case HALO_Z_MINUS:
      xBegin=-1; xEnd=nx+1; yBegin=-1; yEnd=ny+1; zBegin=-1; zEnd=0;
      break;
     case HALO_Z_PLUS:
      xBegin=-1; xEnd=nx+1; yBegin=-1; yEnd=ny+1; zBegin=nz; zEnd=nz+1;
      break;
     default:
      assert(1==0);
   }

   int count = 0;
   for (int ix=xBegin; ix<xEnd; ++ix)
      for (int iy=yBegin; iy<yEnd; ++iy)
         for (int iz=zBegin; iz<zEnd; ++iz)
            list[count++] = getBoxFromTuple(boxes, ix, iy, iz);

   assert(count == nCells);
   return list;
}

/// The loadBuffer function for a force exchange.
/// Iterate the send list and load the derivative of the embedding
/// energy with respect to the local density into the send buffer.
///
/// \see HaloExchangeSt::loadBuffer for an explanation of the loadBuffer
/// parameters.
int loadForceBuffer(void* vparms, void* vdata, int face, char* charBuf)
{
  // TODO: Implement this in parallel with RAJA
   ForceExchangeParms* parms = (ForceExchangeParms*) vparms;
   ForceExchangeData* data = (ForceExchangeData*) vdata;
   ForceMsg* buf = (ForceMsg*) charBuf;
   int nCells = parms->nCells[face];
   int* cellList = parms->sendCells[face];
   int nBuf = 0;
   for (int iCell=0; iCell<nCells; ++iCell)
   {
      int iBox = cellList[iCell];
      int iOff = iBox*MAXATOMS;
      for (int ii=iOff; ii<iOff+data->boxes->nAtoms[iBox]; ++ii)
      {
         buf[nBuf].dfEmbed = data->dfEmbed[ii];
         ++nBuf;
      }
   }
   return nBuf*sizeof(ForceMsg);
}

/// The unloadBuffer function for a force exchange.
/// Data is received in an order that naturally aligns with the atom
/// storage so it is simple to put the data where it belongs.
///
/// \see HaloExchangeSt::unloadBuffer for an explanation of the
/// unloadBuffer parameters.
void unloadForceBuffer(void* vparms, void* vdata, int face, int bufSize, char* charBuf)
{
  // TODO: Implement this in parallel with RAJA
  ForceExchangeParms* parms = (ForceExchangeParms*) vparms;
  ForceExchangeData* data = (ForceExchangeData*) vdata;
  ForceMsg* buf = (ForceMsg*) charBuf;
  assert(bufSize % sizeof(ForceMsg) == 0);

  int nCells = parms->nCells[face];
  int* cellList = parms->recvCells[face];
  int iBuf = 0;
  for (int iCell=0; iCell<nCells; ++iCell)
    {
      int iBox = cellList[iCell];
      int iOff = iBox*MAXATOMS;
      for (int ii=iOff; ii<iOff+data->boxes->nAtoms[iBox]; ++ii)
        {
          data->dfEmbed[ii] = buf[iBuf].dfEmbed;
          ++iBuf;
        }
    }
  assert(iBuf == bufSize/ sizeof(ForceMsg));
}

void destroyForceExchange(void* vparms)
{
   ForceExchangeParms* parms = (ForceExchangeParms*) vparms;

   for (int ii=0; ii<6; ++ii)
   {
      comdFree(parms->sendCells[ii]);
      comdFree(parms->recvCells[ii]);
   }
   comdFree(parms->sendBufM);
   comdFree(parms->sendBufP);
   comdFree(parms->recvBufM);
   comdFree(parms->recvBufP);
}

///  A function suitable for passing to qsort to sort atoms by gid.
///  Because every atom in the simulation is supposed to have a unique
///  id, this function checks that the atoms have different gids.  If
///  that assertion ever fails it is a sign that something has gone
///  wrong elsewhere in the code.
int sortAtomsById(const void* a, const void* b)
{
   int aId = ((AtomMsg*) a)->gid;
   int bId = ((AtomMsg*) b)->gid;
   assert(aId != bId);

   if (aId < bId)
      return -1;
   return 1;
}
