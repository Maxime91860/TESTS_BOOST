#if !defined(USE_MPI)
# error "You should specify USE_MPI=0 or USE_MPI=1 on the compile line"
#endif


// OpenMP will be compiled in if this flag is set to 1 AND the compiler beging
// used supports it (i.e. the _OPENMP symbol is defined)
#define USE_OMP 1

#if USE_MPI
#include <mpi.h>

/*
   define one of these three symbols:

   SEDOV_SYNC_POS_VEL_NONE
   SEDOV_SYNC_POS_VEL_EARLY
   SEDOV_SYNC_POS_VEL_LATE
*/

#define SEDOV_SYNC_POS_VEL_EARLY 1
#endif

#include <math.h>
#include <vector>

//********************
// Boost Serialization
//********************
#include <boost/serialization/vector.hpp>
#include <iostream>
#include <fstream>
#if _OPENMP
#include <omp.h>
#endif

//**************************************************
// Allow flexibility for arithmetic representations 
//**************************************************

#define MAX(a, b) ( ((a) > (b)) ? (a) : (b))


// Precision specification
typedef float        real4 ;
typedef double       real8 ;
typedef long double  real10 ;  // 10 bytes on x86

typedef int    Index_t ; // array subscript and loop index
typedef real8  Real_t ;  // floating point representation
typedef int    Int_t ;   // integer representation

enum { VolumeError = -1, QStopError = -2 } ;

inline real4  SQRT(real4  arg) { return sqrtf(arg) ; }
inline real8  SQRT(real8  arg) { return sqrt(arg) ; }
inline real10 SQRT(real10 arg) { return sqrtl(arg) ; }

inline real4  CBRT(real4  arg) { return cbrtf(arg) ; }
inline real8  CBRT(real8  arg) { return cbrt(arg) ; }
inline real10 CBRT(real10 arg) { return cbrtl(arg) ; }

inline real4  FABS(real4  arg) { return fabsf(arg) ; }
inline real8  FABS(real8  arg) { return fabs(arg) ; }
inline real10 FABS(real10 arg) { return fabsl(arg) ; }


// Stuff needed for boundary conditions
// 2 BCs on each of 6 hexahedral faces (12 bits)
#define XI_M        0x00007
#define XI_M_SYMM   0x00001
#define XI_M_FREE   0x00002
#define XI_M_COMM   0x00004

#define XI_P        0x00038
#define XI_P_SYMM   0x00008
#define XI_P_FREE   0x00010
#define XI_P_COMM   0x00020

#define ETA_M       0x001c0
#define ETA_M_SYMM  0x00040
#define ETA_M_FREE  0x00080
#define ETA_M_COMM  0x00100

#define ETA_P       0x00e00
#define ETA_P_SYMM  0x00200
#define ETA_P_FREE  0x00400
#define ETA_P_COMM  0x00800

#define ZETA_M      0x07000
#define ZETA_M_SYMM 0x01000
#define ZETA_M_FREE 0x02000
#define ZETA_M_COMM 0x04000

#define ZETA_P      0x38000
#define ZETA_P_SYMM 0x08000
#define ZETA_P_FREE 0x10000
#define ZETA_P_COMM 0x20000

// MPI Message Tags
#define MSG_COMM_SBN      1024
#define MSG_SYNC_POS_VEL  2048
#define MSG_MONOQ         3072

#define MAX_FIELDS_PER_MPI_COMM 6

// Assume 128 byte coherence
// Assume Real_t is an "integral power of 2" bytes wide
#define CACHE_COHERENCE_PAD_REAL (128 / sizeof(Real_t))

#define CACHE_ALIGN_REAL(n) \
   (((n) + (CACHE_COHERENCE_PAD_REAL - 1)) & ~(CACHE_COHERENCE_PAD_REAL-1))

//////////////////////////////////////////////////////
// Primary data structure
//////////////////////////////////////////////////////

/*
 * The implementation of the data abstraction used for lulesh
 * resides entirely in the Domain class below.  You can change
 * grouping and interleaving of fields here to maximize data layout
 * efficiency for your underlying architecture or compiler.
 *
 * For example, fields can be implemented as STL objects or
 * raw array pointers.  As another example, individual fields
 * m_x, m_y, m_z could be budled into
 *
 *    struct { Real_t x, y, z ; } *m_coord ;
 *
 * allowing accessor functions such as
 *
 *  "Real_t &x(Index_t idx) { return m_coord[idx].x ; }"
 *  "Real_t &y(Index_t idx) { return m_coord[idx].y ; }"
 *  "Real_t &z(Index_t idx) { return m_coord[idx].z ; }"
 */

class Domain {

   public:

   // Constructor
   Domain(Int_t numRanks, Index_t colLoc,
          Index_t rowLoc, Index_t planeLoc,
          Index_t nx, Int_t tp, Int_t nr, Int_t balance, Int_t cost);

   Domain ()   :
   m_e_cut(Real_t(1.0e-7)),
   m_p_cut(Real_t(1.0e-7)),
   m_q_cut(Real_t(1.0e-7)),
   m_v_cut(Real_t(1.0e-10)),
   m_u_cut(Real_t(1.0e-7)),
   m_hgcoef(Real_t(3.0)),
   m_ss4o3(Real_t(4.0)/Real_t(3.0)),
   m_qstop(Real_t(1.0e+12)),
   m_monoq_max_slope(Real_t(1.0)),
   m_monoq_limiter_mult(Real_t(2.0)),
   m_qlc_monoq(Real_t(0.5)),
   m_qqc_monoq(Real_t(2.0)/Real_t(3.0)),
   m_qqc(Real_t(2.0)),
   m_eosvmax(Real_t(1.0e+9)),
   m_eosvmin(Real_t(1.0e-9)),
   m_pmin(Real_t(0.)),
   m_emin(Real_t(-1.0e+15)),
   m_dvovmax(Real_t(0.1)),
   m_refdens(Real_t(1.0)) {};

   //
   // ALLOCATION
   //

   void AllocateNodePersistent(Int_t numNode) // Node-centered
   {
      m_x.resize(numNode);  // coordinates
      m_y.resize(numNode);
      m_z.resize(numNode);

      m_xd.resize(numNode); // velocities
      m_yd.resize(numNode);
      m_zd.resize(numNode);

      m_xdd.resize(numNode); // accelerations
      m_ydd.resize(numNode);
      m_zdd.resize(numNode);

      m_fx.resize(numNode);  // forces
      m_fy.resize(numNode);
      m_fz.resize(numNode);

      m_nodalMass.resize(numNode);  // mass
   }

   void AllocateElemPersistent(Int_t numElem) // Elem-centered
   {
      m_nodelist.resize(8*numElem);

      // elem connectivities through face
      m_lxim.resize(numElem);
      m_lxip.resize(numElem);
      m_letam.resize(numElem);
      m_letap.resize(numElem);
      m_lzetam.resize(numElem);
      m_lzetap.resize(numElem);

      m_elemBC.resize(numElem);

      m_e.resize(numElem);
      m_p.resize(numElem);

      m_q.resize(numElem);
      m_ql.resize(numElem);
      m_qq.resize(numElem);

      m_v.resize(numElem);

      m_volo.resize(numElem);
      m_delv.resize(numElem);
      m_vdov.resize(numElem);

      m_arealg.resize(numElem);

      m_ss.resize(numElem);

      m_elemMass.resize(numElem);
   }

   void AllocateGradients(Int_t numElem, Int_t allElem)
   {
      // Position gradients
      m_delx_xi.resize(numElem) ;
      m_delx_eta.resize(numElem) ;
      m_delx_zeta.resize(numElem) ;

      // Velocity gradients
      m_delv_xi.resize(allElem) ;
      m_delv_eta.resize(allElem);
      m_delv_zeta.resize(allElem) ;
   }

   void DeallocateGradients()
   {
      m_delx_zeta.clear() ;
      m_delx_eta.clear() ;
      m_delx_xi.clear() ;

      m_delv_zeta.clear() ;
      m_delv_eta.clear() ;
      m_delv_xi.clear() ;
   }

   void AllocateStrains(Int_t numElem)
   {
      m_dxx.resize(numElem) ;
      m_dyy.resize(numElem) ;
      m_dzz.resize(numElem) ;
   }

   void DeallocateStrains()
   {
      m_dzz.clear() ;
      m_dyy.clear() ;
      m_dxx.clear() ;
   }
   
   //
   // ACCESSORS
   //

   // Node-centered

   // Nodal coordinates
   Real_t& x(Index_t idx)    { return m_x[idx] ; }
   Real_t& y(Index_t idx)    { return m_y[idx] ; }
   Real_t& z(Index_t idx)    { return m_z[idx] ; }

   // Nodal velocities
   Real_t& xd(Index_t idx)   { return m_xd[idx] ; }
   Real_t& yd(Index_t idx)   { return m_yd[idx] ; }
   Real_t& zd(Index_t idx)   { return m_zd[idx] ; }

   // Nodal accelerations
   Real_t& xdd(Index_t idx)  { return m_xdd[idx] ; }
   Real_t& ydd(Index_t idx)  { return m_ydd[idx] ; }
   Real_t& zdd(Index_t idx)  { return m_zdd[idx] ; }

   // Nodal forces
   Real_t& fx(Index_t idx)   { return m_fx[idx] ; }
   Real_t& fy(Index_t idx)   { return m_fy[idx] ; }
   Real_t& fz(Index_t idx)   { return m_fz[idx] ; }

   // Nodal mass
   Real_t& nodalMass(Index_t idx) { return m_nodalMass[idx] ; }

   // Nodes on symmertry planes
   Index_t symmX(Index_t idx) { return m_symmX[idx] ; }
   Index_t symmY(Index_t idx) { return m_symmY[idx] ; }
   Index_t symmZ(Index_t idx) { return m_symmZ[idx] ; }
   bool symmXempty()          { return m_symmX.empty(); }
   bool symmYempty()          { return m_symmY.empty(); }
   bool symmZempty()          { return m_symmZ.empty(); }

   //
   // Element-centered
   //
   Index_t&  regElemSize(Index_t idx) { return m_regElemSize[idx] ; }
   Index_t&  regNumList(Index_t idx) { return m_regNumList[idx] ; }
   Index_t*  regNumList()            { return &m_regNumList[0] ; }
   Index_t*  regElemlist(Int_t r)    { return m_regElemlist[r] ; }
   Index_t&  regElemlist(Int_t r, Index_t idx) { return m_regElemlist[r][idx] ; }

   Index_t*  nodelist(Index_t idx)    { return &m_nodelist[Index_t(8)*idx] ; }

   // elem connectivities through face
   Index_t&  lxim(Index_t idx) { return m_lxim[idx] ; }
   Index_t&  lxip(Index_t idx) { return m_lxip[idx] ; }
   Index_t&  letam(Index_t idx) { return m_letam[idx] ; }
   Index_t&  letap(Index_t idx) { return m_letap[idx] ; }
   Index_t&  lzetam(Index_t idx) { return m_lzetam[idx] ; }
   Index_t&  lzetap(Index_t idx) { return m_lzetap[idx] ; }

   // elem face symm/free-surface flag
   Int_t&  elemBC(Index_t idx) { return m_elemBC[idx] ; }

   // Principal strains - temporary
   Real_t& dxx(Index_t idx)  { return m_dxx[idx] ; }
   Real_t& dyy(Index_t idx)  { return m_dyy[idx] ; }
   Real_t& dzz(Index_t idx)  { return m_dzz[idx] ; }

   // Velocity gradient - temporary
   Real_t& delv_xi(Index_t idx)    { return m_delv_xi[idx] ; }
   Real_t& delv_eta(Index_t idx)   { return m_delv_eta[idx] ; }
   Real_t& delv_zeta(Index_t idx)  { return m_delv_zeta[idx] ; }

   // Position gradient - temporary
   Real_t& delx_xi(Index_t idx)    { return m_delx_xi[idx] ; }
   Real_t& delx_eta(Index_t idx)   { return m_delx_eta[idx] ; }
   Real_t& delx_zeta(Index_t idx)  { return m_delx_zeta[idx] ; }

   // Energy
   Real_t& e(Index_t idx)          { return m_e[idx] ; }

   // Pressure
   Real_t& p(Index_t idx)          { return m_p[idx] ; }

   // Artificial viscosity
   Real_t& q(Index_t idx)          { return m_q[idx] ; }

   // Linear term for q
   Real_t& ql(Index_t idx)         { return m_ql[idx] ; }
   // Quadratic term for q
   Real_t& qq(Index_t idx)         { return m_qq[idx] ; }

   // Relative volume
   Real_t& v(Index_t idx)          { return m_v[idx] ; }
   Real_t& delv(Index_t idx)       { return m_delv[idx] ; }

   // Reference volume
   Real_t& volo(Index_t idx)       { return m_volo[idx] ; }

   // volume derivative over volume
   Real_t& vdov(Index_t idx)       { return m_vdov[idx] ; }

   // Element characteristic length
   Real_t& arealg(Index_t idx)     { return m_arealg[idx] ; }

   // Sound speed
   Real_t& ss(Index_t idx)         { return m_ss[idx] ; }

   // Element mass
   Real_t& elemMass(Index_t idx)  { return m_elemMass[idx] ; }

   Index_t nodeElemCount(Index_t idx)
   { return m_nodeElemStart[idx+1] - m_nodeElemStart[idx] ; }

   Index_t *nodeElemCornerList(Index_t idx)
   { return &m_nodeElemCornerList[m_nodeElemStart[idx]] ; }

   // Parameters 

   // Cutoffs
   Real_t u_cut() const               { return m_u_cut ; }
   Real_t e_cut() const               { return m_e_cut ; }
   Real_t p_cut() const               { return m_p_cut ; }
   Real_t q_cut() const               { return m_q_cut ; }
   Real_t v_cut() const               { return m_v_cut ; }

   // Other constants (usually are settable via input file in real codes)
   Real_t hgcoef() const              { return m_hgcoef ; }
   Real_t qstop() const               { return m_qstop ; }
   Real_t monoq_max_slope() const     { return m_monoq_max_slope ; }
   Real_t monoq_limiter_mult() const  { return m_monoq_limiter_mult ; }
   Real_t ss4o3() const               { return m_ss4o3 ; }
   Real_t qlc_monoq() const           { return m_qlc_monoq ; }
   Real_t qqc_monoq() const           { return m_qqc_monoq ; }
   Real_t qqc() const                 { return m_qqc ; }

   Real_t eosvmax() const             { return m_eosvmax ; }
   Real_t eosvmin() const             { return m_eosvmin ; }
   Real_t pmin() const                { return m_pmin ; }
   Real_t emin() const                { return m_emin ; }
   Real_t dvovmax() const             { return m_dvovmax ; }
   Real_t refdens() const             { return m_refdens ; }

   // Timestep controls, etc...
   Real_t& time()                 { return m_time ; }
   Real_t& deltatime()            { return m_deltatime ; }
   Real_t& deltatimemultlb()      { return m_deltatimemultlb ; }
   Real_t& deltatimemultub()      { return m_deltatimemultub ; }
   Real_t& stoptime()             { return m_stoptime ; }
   Real_t& dtcourant()            { return m_dtcourant ; }
   Real_t& dthydro()              { return m_dthydro ; }
   Real_t& dtmax()                { return m_dtmax ; }
   Real_t& dtfixed()              { return m_dtfixed ; }

   Int_t&  cycle()                { return m_cycle ; }
   Index_t&  numRanks()           { return m_numRanks ; }

   Index_t&  colLoc()             { return m_colLoc ; }
   Index_t&  rowLoc()             { return m_rowLoc ; }
   Index_t&  planeLoc()           { return m_planeLoc ; }
   Index_t&  tp()                 { return m_tp ; }

   Index_t&  sizeX()              { return m_sizeX ; }
   Index_t&  sizeY()              { return m_sizeY ; }
   Index_t&  sizeZ()              { return m_sizeZ ; }
   Index_t&  numReg()             { return m_numReg ; }
   Int_t&  cost()             { return m_cost ; }
   Index_t&  numElem()            { return m_numElem ; }
   Index_t&  numNode()            { return m_numNode ; }
   
   Index_t&  maxPlaneSize()       { return m_maxPlaneSize ; }
   Index_t&  maxEdgeSize()        { return m_maxEdgeSize ; }
   
   //
   // MPI-Related additional data
   //

#if USE_MPI   
   // Communication Work space 
   Real_t *commDataSend ;
   Real_t *commDataRecv ;
   
   // Maximum number of block neighbors 
   MPI_Request recvRequest[26] ; // 6 faces + 12 edges + 8 corners 
   MPI_Request sendRequest[26] ; // 6 faces + 12 edges + 8 corners 
#endif

  private:

   void BuildMesh(Int_t nx, Int_t edgeNodes, Int_t edgeElems);
   void SetupThreadSupportStructures();
   void CreateRegionIndexSets(Int_t nreg, Int_t balance);
   void SetupCommBuffers(Int_t edgeNodes);
   void SetupSymmetryPlanes(Int_t edgeNodes);
   void SetupElementConnectivities(Int_t edgeElems);
   void SetupBoundaryConditions(Int_t edgeElems);

   friend class boost::serialization::access;
   template <typename Archive>
   void serialize(Archive &ar, const unsigned int version){


      //Check de/serialization
      // if(Archive::is_loading::value){
      //    std::cout << "-------------------------\n";
      //    std::cout << "Start of deserialization.\n";
      //    std::cout << "-------------------------\n";
      // }
      // else {
      //    std::cout << "-------------------------\n";
      //    std::cout << "Start of serialization.\n";
      //    std::cout << "-------------------------\n";
      // }

      ar & m_x ;  /* coordinates */
      ar & m_y;
      ar & m_z;

      ar & m_xd ; /* velocities */
      ar & m_yd ;
      ar & m_zd ;

      ar & m_xdd ; /* accelerations */
      ar & m_ydd ;
      ar & m_zdd ;

      ar & m_fx ;  /* forces */
      ar & m_fy ;
      ar & m_fz ;

      ar & m_nodalMass ;  /* mass */

      ar & m_symmX ;  /* symmetry plane nodesets */
      ar & m_symmY ;
      ar & m_symmZ ;

      // Element-centered

      ar & m_numRanks ;
      ar & m_colLoc ;
      ar & m_rowLoc ;
      ar & m_planeLoc ;
      ar & m_tp ;

      ar & m_sizeX ;
      ar & m_sizeY ;
      ar & m_sizeZ ;
      ar & m_numElem ;
      ar & m_numNode ;

      ar & m_maxPlaneSize ;
      ar & m_maxEdgeSize ;

      // Region information
      ar &  m_numReg ;
      ar &  m_cost; //imbalance cost
      
      if(Archive::is_loading::value){
         m_regElemSize = new Index_t[m_numReg];
      }
      ar & boost::serialization::make_array <Index_t> (m_regElemSize, m_numReg); // Size of region sets

      if(Archive::is_loading::value){
         m_regNumList = new Index_t[m_numElem];
      }
      ar & boost::serialization::make_array <Index_t> (m_regNumList, m_numElem); // Region number per domain element

      if(Archive::is_loading::value){
         m_regElemlist = new Index_t*[m_numReg];
         for (int i = 0; i < m_numReg; i++){
            m_regElemlist[i] = new Index_t[m_regElemSize[i]];
         }
      }

      for (int i = 0; i < m_numReg; i++){
         ar & boost::serialization::make_array <Index_t> (m_regElemlist[i], m_regElemSize[i]);
      }

      ar &  m_nodelist ;     /* elemToNode connectivity */

      ar & m_lxim ;  /* element connectivity across each face */
      ar & m_lxip ;
      ar & m_letam ;
      ar & m_letap ;
      ar & m_lzetam ;
      ar & m_lzetap ;

      ar & m_elemBC ;  /* symmetry/free-surface flags for each elem face */

      ar & m_dxx ;  /* principal strains -- temporary */
      ar & m_dyy ;
      ar & m_dzz ;

      ar & m_delv_xi ;    /* velocity gradient -- temporary */
      ar & m_delv_eta ;
      ar & m_delv_zeta ;

      ar & m_delx_xi ;    /* coordinate gradient -- temporary */
      ar & m_delx_eta ;
      ar & m_delx_zeta ;

      ar & m_e ;   /* energy */

      ar & m_p ;   /* pressure */
      ar & m_q ;   /* q */
      ar & m_ql ;  /* linear term for q */
      ar & m_qq ;  /* quadratic term for q */

      ar & m_v ;     /* relative volume */
      ar & m_volo ;  /* reference volume */
      ar & m_vnew ;  /* new relative volume -- temporary */
      ar & m_delv ;  /* m_vnew - m_v */
      ar & m_vdov ;  /* volume derivative over volume */

      ar & m_arealg ;  /* characteristic length of an element */

      ar & m_ss ;      /* "sound speed" */

      ar & m_elemMass ;  /* mass */

      // Cutoffs (treat as constants)
      ar & const_cast<Real_t &>(m_e_cut);
      ar & const_cast<Real_t &>(m_p_cut);
      ar & const_cast<Real_t &>(m_q_cut);
      ar & const_cast<Real_t &>(m_v_cut);
      ar & const_cast<Real_t &>(m_u_cut);

      // Other constants (usually setable, but hardcoded in this proxy app)
      ar & const_cast<Real_t &>(m_hgcoef);
      ar & const_cast<Real_t &>(m_ss4o3);
      ar & const_cast<Real_t &>(m_qstop);
      ar & const_cast<Real_t &>(m_monoq_max_slope);
      ar & const_cast<Real_t &>(m_monoq_limiter_mult);  
      ar & const_cast<Real_t &>(m_qlc_monoq);
      ar & const_cast<Real_t &>(m_qqc_monoq);
      ar & const_cast<Real_t &>(m_qqc);
      ar & const_cast<Real_t &>(m_eosvmax);
      ar & const_cast<Real_t &>(m_eosvmin);
      ar & const_cast<Real_t &>(m_pmin);
      ar & const_cast<Real_t &>(m_emin);
      ar & const_cast<Real_t &>(m_dvovmax);
      ar & const_cast<Real_t &>(m_refdens);

      // Variables to keep track of timestep, simulation time, and cycle
      ar &  m_dtcourant ;         // courant constraint 
      ar &  m_dthydro ;           // volume change constraint 
      ar &   m_cycle ;             // iteration count for simulation 
      ar &  m_dtfixed ;           // fixed time increment 
      ar &  m_time ;              // current time 
      ar &  m_deltatime ;         // variable time increment 
      ar &  m_deltatimemultlb ;
      ar &  m_deltatimemultub ;
      ar &  m_dtmax ;             // maximum allowable time increment 
      ar &  m_stoptime ;          // end time for simulation 

      // OMP hack 
      #if _OPENMP
         Index_t numthreads = omp_get_max_threads();
      #else
         Index_t numthreads = 1;
      #endif

      if (numthreads > 1) {
         if(Archive::is_loading::value){
            m_nodeElemStart = new Index_t[m_numNode+1];
         }
         ar & boost::serialization::make_array <Index_t> (m_nodeElemStart, m_numNode+1);

         if(Archive::is_loading::value){
            m_nodeElemCornerList = new Index_t[m_nodeElemStart[m_numNode]];
         }
         ar & boost::serialization::make_array <Index_t> (m_nodeElemCornerList, m_nodeElemStart[m_numNode]);
      } else {
         m_nodeElemStart = NULL;
         m_nodeElemCornerList = NULL;
      }

      // Used in setup
      ar & m_rowMin;
      ar & m_rowMax;
      ar & m_colMin;
      ar & m_colMax;
      ar & m_planeMin;
      ar & m_planeMax; 

     #if USE_MPI   
     // account for face communication 
     Index_t comBufSize =
       (m_rowMin + m_rowMax + m_colMin + m_colMax + m_planeMin + m_planeMax) *
       m_maxPlaneSize * MAX_FIELDS_PER_MPI_COMM ;

     // account for edge communication 
     comBufSize +=
       ((m_rowMin & m_colMin) + (m_rowMin & m_planeMin) + (m_colMin & m_planeMin) +
        (m_rowMax & m_colMax) + (m_rowMax & m_planeMax) + (m_colMax & m_planeMax) +
        (m_rowMax & m_colMin) + (m_rowMin & m_planeMax) + (m_colMin & m_planeMax) +
        (m_rowMin & m_colMax) + (m_rowMax & m_planeMin) + (m_colMax & m_planeMin)) *
       m_maxEdgeSize * MAX_FIELDS_PER_MPI_COMM ;

     // account for corner communication 
     // factor of 16 is so each buffer has its own cache line 
     comBufSize += ((m_rowMin & m_colMin & m_planeMin) +
          (m_rowMin & m_colMin & m_planeMax) +
          (m_rowMin & m_colMax & m_planeMin) +
          (m_rowMin & m_colMax & m_planeMax) +
          (m_rowMax & m_colMin & m_planeMin) +
          (m_rowMax & m_colMin & m_planeMax) +
          (m_rowMax & m_colMax & m_planeMin) +
          (m_rowMax & m_colMax & m_planeMax)) * CACHE_COHERENCE_PAD_REAL ;


      // Communication Work space 
      if(Archive::is_loading::value){
         commDataSend = new Real_t[comBufSize];
         commDataRecv = new Real_t[comBufSize];
      }
      ar & boost::serialization::make_array <Real_t> (commDataRecv,comBufSize);
      ar & boost::serialization::make_array <Real_t> (commDataSend,comBufSize);

      #endif

      //Check de/serialization
      // if(Archive::is_loading::value){
      //    std::cout << "-------------------------\n";
      //    std::cout << "Deserialization finished.\n";
      //    std::cout << "-------------------------\n";
      // }
      // else {
      //    std::cout << "-------------------------\n";
      //    std::cout << "Serialization finished.\n";
      //    std::cout << "-------------------------\n";
      // }
   }

   //
   // IMPLEMENTATION
   //

   /* Node-centered */
   std::vector<Real_t> m_x ;  /* coordinates */
   std::vector<Real_t> m_y ;
   std::vector<Real_t> m_z ;

   std::vector<Real_t> m_xd ; /* velocities */
   std::vector<Real_t> m_yd ;
   std::vector<Real_t> m_zd ;

   std::vector<Real_t> m_xdd ; /* accelerations */
   std::vector<Real_t> m_ydd ;
   std::vector<Real_t> m_zdd ;

   std::vector<Real_t> m_fx ;  /* forces */
   std::vector<Real_t> m_fy ;
   std::vector<Real_t> m_fz ;

   std::vector<Real_t> m_nodalMass ;  /* mass */

   std::vector<Index_t> m_symmX ;  /* symmetry plane nodesets */
   std::vector<Index_t> m_symmY ;
   std::vector<Index_t> m_symmZ ;

   // Element-centered

   // Region information
   Int_t    m_numReg ;
   Int_t    m_cost; //imbalance cost
   Index_t *m_regElemSize ;   // Size of region sets
   Index_t *m_regNumList ;    // Region number per domain element
   Index_t **m_regElemlist ;  // region indexset 

   std::vector<Index_t>  m_nodelist ;     /* elemToNode connectivity */

   std::vector<Index_t>  m_lxim ;  /* element connectivity across each face */
   std::vector<Index_t>  m_lxip ;
   std::vector<Index_t>  m_letam ;
   std::vector<Index_t>  m_letap ;
   std::vector<Index_t>  m_lzetam ;
   std::vector<Index_t>  m_lzetap ;

   std::vector<Int_t>    m_elemBC ;  /* symmetry/free-surface flags for each elem face */

   std::vector<Real_t> m_dxx ;  /* principal strains -- temporary */
   std::vector<Real_t> m_dyy ;
   std::vector<Real_t> m_dzz ;

   std::vector<Real_t> m_delv_xi ;    /* velocity gradient -- temporary */
   std::vector<Real_t> m_delv_eta ;
   std::vector<Real_t> m_delv_zeta ;

   std::vector<Real_t> m_delx_xi ;    /* coordinate gradient -- temporary */
   std::vector<Real_t> m_delx_eta ;
   std::vector<Real_t> m_delx_zeta ;
   
   std::vector<Real_t> m_e ;   /* energy */

   std::vector<Real_t> m_p ;   /* pressure */
   std::vector<Real_t> m_q ;   /* q */
   std::vector<Real_t> m_ql ;  /* linear term for q */
   std::vector<Real_t> m_qq ;  /* quadratic term for q */

   std::vector<Real_t> m_v ;     /* relative volume */
   std::vector<Real_t> m_volo ;  /* reference volume */
   std::vector<Real_t> m_vnew ;  /* new relative volume -- temporary */
   std::vector<Real_t> m_delv ;  /* m_vnew - m_v */
   std::vector<Real_t> m_vdov ;  /* volume derivative over volume */

   std::vector<Real_t> m_arealg ;  /* characteristic length of an element */
   
   std::vector<Real_t> m_ss ;      /* "sound speed" */

   std::vector<Real_t> m_elemMass ;  /* mass */

   // Cutoffs (treat as constants)
   const Real_t  m_e_cut ;             // energy tolerance 
   const Real_t  m_p_cut ;             // pressure tolerance 
   const Real_t  m_q_cut ;             // q tolerance 
   const Real_t  m_v_cut ;             // relative volume tolerance 
   const Real_t  m_u_cut ;             // velocity tolerance 

   // Other constants (usually setable, but hardcoded in this proxy app)

   const Real_t  m_hgcoef ;            // hourglass control 
   const Real_t  m_ss4o3 ;
   const Real_t  m_qstop ;             // excessive q indicator 
   const Real_t  m_monoq_max_slope ;
   const Real_t  m_monoq_limiter_mult ;
   const Real_t  m_qlc_monoq ;         // linear term coef for q 
   const Real_t  m_qqc_monoq ;         // quadratic term coef for q 
   const Real_t  m_qqc ;
   const Real_t  m_eosvmax ;
   const Real_t  m_eosvmin ;
   const Real_t  m_pmin ;              // pressure floor 
   const Real_t  m_emin ;              // energy floor 
   const Real_t  m_dvovmax ;           // maximum allowable volume change 
   const Real_t  m_refdens ;           // reference density 

   // Variables to keep track of timestep, simulation time, and cycle
   Real_t  m_dtcourant ;         // courant constraint 
   Real_t  m_dthydro ;           // volume change constraint 
   Int_t   m_cycle ;             // iteration count for simulation 
   Real_t  m_dtfixed ;           // fixed time increment 
   Real_t  m_time ;              // current time 
   Real_t  m_deltatime ;         // variable time increment 
   Real_t  m_deltatimemultlb ;
   Real_t  m_deltatimemultub ;
   Real_t  m_dtmax ;             // maximum allowable time increment 
   Real_t  m_stoptime ;          // end time for simulation 


   Int_t   m_numRanks ;

   Index_t m_colLoc ;
   Index_t m_rowLoc ;
   Index_t m_planeLoc ;
   Index_t m_tp ;

   Index_t m_sizeX ;
   Index_t m_sizeY ;
   Index_t m_sizeZ ;
   Index_t m_numElem ;
   Index_t m_numNode ;

   Index_t m_maxPlaneSize ;
   Index_t m_maxEdgeSize ;

   // OMP hack 
   Index_t *m_nodeElemStart ;
   Index_t *m_nodeElemCornerList ;

   // Used in setup
   Index_t m_rowMin, m_rowMax;
   Index_t m_colMin, m_colMax;
   Index_t m_planeMin, m_planeMax ;

} ;

typedef Real_t &(Domain::* Domain_member )(Index_t) ;

struct cmdLineOpts {
   Int_t its; // -i 
   Int_t nx;  // -s 
   Int_t numReg; // -r 
   Int_t numFiles; // -f
   Int_t showProg; // -p
   Int_t quiet; // -q
   Int_t viz; // -v 
   Int_t cost; // -c
   Int_t balance; // -b
};



// Function Prototypes

// lulesh-par
Real_t CalcElemVolume( const Real_t x[8],
                       const Real_t y[8],
                       const Real_t z[8]);

// lulesh-util
void ParseCommandLineOptions(int argc, char *argv[],
                             Int_t myRank, struct cmdLineOpts *opts);
void VerifyAndWriteFinalOutput(Real_t elapsed_time,
                               Domain& locDom,
                               Int_t nx,
                               Int_t numRanks);

// lulesh-viz
void DumpToVisit(Domain& domain, int numFiles, int myRank, int numRanks);

// lulesh-comm
void CommRecv(Domain& domain, Int_t msgType, Index_t xferFields,
              Index_t dx, Index_t dy, Index_t dz,
              bool doRecv, bool planeOnly);
void CommSend(Domain& domain, Int_t msgType,
              Index_t xferFields, Domain_member *fieldData,
              Index_t dx, Index_t dy, Index_t dz,
              bool doSend, bool planeOnly);
void CommSBN(Domain& domain, Int_t xferFields, Domain_member *fieldData);
void CommSyncPosVel(Domain& domain);
void CommMonoQ(Domain& domain);

// lulesh-init
void InitMeshDecomp(Int_t numRanks, Int_t myRank,
                    Int_t *col, Int_t *row, Int_t *plane, Int_t *side);
