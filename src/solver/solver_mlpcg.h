/*!---------------------------------------------------------------------
\file
\brief data types for multilevel preconditioned cg

---------------------------------------------------------------------*/

/*! 
\addtogroup MLPCG 
*//*! @{ (documentation module open)*/

/*!------------------------------------------------------------------------
\brief presmoother enum

m.gee 9/02  

presmoother enum       

-------------------------------------------------------------------------*/
typedef enum _MLPCG_PRESMOOTH
{
               pre_none,
               pre_fwdGS,           /*!< Gauss-Seidel smoother */
               pre_Jacobi,          /*!< Jacobi smoother */
               pre_ilu              /*!< ilu(n) smoother with n=ilu_n */
} MLPCG_PRESMOOTH;
/*!------------------------------------------------------------------------
\brief presmoother enum

m.gee 9/02  

postmoother enum       

-------------------------------------------------------------------------*/
typedef enum _MLPCG_POSTMOOTH
{
               post_none,
               post_bckGS,           /*!< Gauss-Seidel smoother */
               post_Jacobi,          /*!< Jacobi smoother */
               post_ilu              /*!< ilu(n) smoother with n=ilu_n */
} MLPCG_POSTMOOTH;
/*!------------------------------------------------------------------------
\brief coarsesolver enum

m.gee 9/02  

coarsesolver enum       

-------------------------------------------------------------------------*/
typedef enum _MLPCG_COARSESOLVE
{
               co_none,
               co_ilu,           /*!< ilu solver */
               co_spooles,       /*!< lapack solver */
               co_lapack         /*!< spooles solver */
} MLPCG_COARSESOLVE;


/*!------------------------------------------------------------------------
\brief variables needed by the MLPCG solver

m.gee 9/02  

variables needed by the MLPCG solver       

-------------------------------------------------------------------------*/
typedef struct _MLPCGVARS
{
int                     numlev;      /*!< number of grids in the ml precond. */
double                  p_omega;     /*!< damping of the prolongator smoother */
int                     numdf;       /*!< number of coarse grid dofs per node */
int                     ilu_n;       /*!< ilu(n) */
double                  tol;         /*!< cg-tolerance */
int                     maxiter;     /*!< max number of cg iterations */
enum _MLPCG_COARSESOLVE coarsesolv;  /*!< coarsest level solver */
int                     co_ilu_n;    /*!< coarsest level ilu(n) in case of ilu solver */
enum _MLPCG_PRESMOOTH   presmoother; /*!< upgoing smoother on levels */
int                     presweep;    /*!< number of presmoothing sweeps */
enum _MLPCG_POSTMOOTH   postsmoother;/*!< upgoing smoother on levels */
int                     postsweep;   /*!< number of postsmoothing sweeps */
struct _DISCRET        *fielddis;    /*!< warning: this is a pointer to the field original, not a copy! */
struct _PARTDISCRET    *partdis;     /*!< warning: this is a pointer to the partition original, not a copy! */
} MLPCGVARS;
/*!------------------------------------------------------------------------
\brief matrix needed by the MLPCG solver on the finest grid only   

m.gee 6/01  

matrix needed by the MLPCG solver on the finest grid only      

-------------------------------------------------------------------------*/
typedef struct _DBCSR
{
int                     is_init;       /*!< was this matrix initialized ? */
int                     is_factored;   /*!< is this matrix already factored ? */
int                     ncall;         /*!< how often was this matrix solved */

int                     numeq_total;   /*!< total number of unknowns */
int                     numeq;         /*!< number of unknowns updated on this proc */ 
int                     nnz;           /*!< number of nonzeros on this proc */

int                     owner[MAXPROC][2]; /*!< contains for each proc the lowest and highest
                                                dof number 
                                                owner[proc][0] lowest dof on proc 
                                                owner[proc][1] highest dof on proc */
struct _ARRAY           update;        /*!< list of dofs updated on this proc */
struct _ARRAY           a;             /*!< the values of the sparse matrix */
struct _ARRAY           ja;            /*!< the column indizes of the sparse matrix */
struct _ARRAY           ia;            /*!< the row indizes of the sparse matrix */
struct _ARRAY           blocks;        /*!< nodal block information of csr matrix 
                                            block[i][0] = size of nodal block
                                            block[i][1..block[i][0]] = dofs to this block */
struct _ARRAY           gdofrecv;      /*!< ghost dof list to receive 
                                            gdofrecv.a.ia[0..fdim-1]    = dof numbers of external dofs needed 
                                            sorted in ascending order */
struct _ARRAY           recvbuff;
struct _ARRAY           computebuff;   /*!< after receiving all messages they are sorted to computebuff */
#ifdef PARALLEL
MPI_Status             *status;        /*!< receive status */
#endif
struct _ARRAY           gdofsend;      /*!< ghost dof list to send to other procs (3D integer)
                                            gdofsend.a.ia[proc][0] = number of values to send to proc proc 
                                            gdofsend.a.ia[proc][ 1..gdofsend.a.i3[proc][0] ] = 
                                            dof numbers to be send to proc proc in ascending order */
struct _ARRAY           sendbuff;
#ifdef PARALLEL
MPI_Request            *request;       /*!< send request */
#endif
int                     firstcoupledof;/*!< dof number of the first dof that has interproc coupling */

struct _DBCSR          *csc;           /*!< the treansposed matrix in compressed sparse column format */
struct _DBCSR          *ilu;           /*!< the ilu-decomposed matrix */
ARRAY                  *dense;         /*!< for dense solve on coarsest grid */
ARRAY                  *ipiv;          /*!< for dense solve on coarsest grid */

} DBCSR;
/*!------------------------------------------------------------------------
\brief one level of the multilevel preconditioner

m.gee 6/01  

one level of the multilevel preconditioner       

-------------------------------------------------------------------------*/
typedef struct _MLLEVEL
{
enum _MLPCG_COARSESOLVE coarsesolv;  /*!< coarsest level solver */
int                     co_ilu_n;    /*!< coarsest level ilu(n) in case of ilu solver */
enum _MLPCG_PRESMOOTH   presmoother; /*!< upgoing smoother on levels */
int                     presweep;    /*!< number of presmoothing sweeps */
enum _MLPCG_POSTMOOTH   postsmoother;/*!< upgoing smoother on levels */
int                     postsweep;   /*!< number of postsmoothing sweeps */

struct _DBCSR          *csr;         /*!< the sparse matrix of this level */

int                     nagg;        /*!< number of aggregates on this level */
struct _AGG            *agg;         /*!< vector of aggregates */

struct _DBCSR          *P;           /*!< the prolongator from this level+1 to this level */
} MLLEVEL;
/*!------------------------------------------------------------------------
\brief one level of the multilevel preconditioner

m.gee 6/01  

one level of the multilevel preconditioner       

-------------------------------------------------------------------------*/
typedef struct _AGG
{
int                     nblock;      /*!<  number of blocks in this aggregate */          
int                   **block;       /*!<  list of dof-blocks, which belong to this aggregate */
int                     numdf;       /*!<  number of dofs of this aggregate */
int                    *dof;         /*!<  dofs of the aggregate's supernode */
struct _ARRAY          *tentP;       /*!< this aggregates piece of tentative prolongator */
int                     tentP_nrow;  /*!< number of rows in this piece of tent. Prolongator */
int                    *tentP_rindex;/*!< the row indizes of this piece of tentative prolongator */
                                     /*   (the column indizes are the dofs) */
struct _ARRAY          *R;           /*!<  The R-part of the P=QR factorization */
enum
   {
    mlpcg_aggnone,
    mlpcg_agglocal,
    mlpcg_agginterproc
   }                    coupling;    /*!<  flag indicating whether this aggregate is near intprocessor bounday */
} AGG;
/*!------------------------------------------------------------------------
\brief the ml preconditioner

m.gee 6/01  

this structure holds the complete ml preconditioner       

-------------------------------------------------------------------------*/
typedef struct _MLPRECOND
{
int                     numlev;
struct _MLLEVEL        *level;           /*!< vector of levels */
struct _DISCRET        *fielddis;        /*!< warning: this is a pointer to the field original, not a copy! */
struct _PARTDISCRET    *partdis;         /*!< warning: this is a pointer to the partition original, not a copy! */
struct _ARRAY           director;        /*!< the nodal directors of the shell8 associated nodes */
NODE                  **node;            /*!< ptr to the nodes of the partition */
double                  omega;           /*!< damping factor of the prolongator smoother */
} MLPRECOND;
/*!------------------------------------------------------------------------
\brief the ml solver

m.gee 6/01  

this structure holds everything for the ml-cg algorithm       

-------------------------------------------------------------------------*/
typedef struct _MLSOLVER
{
double                     tol;         /*!< tolerance for the cg algorithm */
int                        maxiter;     /*!< max number of iterations for the cg algorithm */
struct _ARRAY              r;           /*!< iterate vector of the cg-algorithm */     
struct _ARRAY              z;           /*!< iterate vector of the cg-algorithm */
struct _ARRAY              p;           /*!< iterate vector of the cg-algorithm */
struct _ARRAY              q;           /*!< iterate vector of the cg-algorithm */
} MLSOLVER;




/*! @} (documentation module close)*/

