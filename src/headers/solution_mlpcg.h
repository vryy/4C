/*!---------------------------------------------------------------------
\file
\brief data types for multilevel preconditioned cg

---------------------------------------------------------------------*/

/*! 
\addtogroup MLPCG 
*//*! @{ (documentation module open)*/

/*!------------------------------------------------------------------------
\brief variables needed by the MLPCG solver

m.gee 6/01  

variables needed by the MLPCG solver       

-------------------------------------------------------------------------*/
typedef struct _MLPCGVARS
{
int                     numlev;     /*!< number of grids in the ml precond. */
double                  p_omega;    /*!< damping of the prolongator smoother */
int                     ilu_n;      /*!< ilu(n) */
double                  tol;        /*!< cg-tolerance */
int                     co_ilu_n;   /*!< coarsest level ilu(n) in case of ilu solver */
enum
   {
    co_ilu,                         /*!< ilu solver */
    co_spooles,                     /*!< lapack solver */
    co_lapack                       /*!< spooles solver */
   }                    coarsesolv; /*!< coarsest level solver */
enum
   {
    pre_fwdGS,                      /*!< Gauss-Seidel smoother */
    pre_Jacobi,                     /*!< Jacobi smoother */
    pre_ilu                         /*!< ilu(n) smoother with n=ilu_n */
   }                    presmoother;/*!< upgoing smoother on levels */
int                     presweep;   /*!< number of presmoothing sweeps */
enum
   {
    post_bckGS,                      /*!< Gauss-Seidel smoother */
    post_Jacobi,                     /*!< Jacobi smoother */
    post_ilu                         /*!< ilu(n) smoother with n=ilu_n */
   }                    postsmoother;/*!< upgoing smoother on levels */
int                     postsweep;   /*!< number of postsmoothing sweeps */
} MLPCGVARS;




/*!------------------------------------------------------------------------
\brief matrix needed by the MLPCG solver on the finest grid only   

m.gee 6/01  

matrix needed by the MLPCG solver on the finest grid only      

-------------------------------------------------------------------------*/
typedef struct _DBCSR_ROOT
{
int                     is_init;         /*!< was this matrix initialized ? */
int                     is_factored;     /*!< is this matrix already factored ? */
int                     ncall;           /*!< how often was this matrix solved */

int                     numeq_total;     /*!< total number of unknowns */
int                     numeq;           /*!< number of unknowns updated on this proc */ 
int                     nnz;             /*!< number of nonzeros on this proc */

struct _ARRAY           update;          /*!< list of dofs updated on this proc */
struct _ARRAY           a;               /*!< the values of the sparse matrix */
struct _ARRAY           ja;              /*!< the column indizes of the sparse matrix */
struct _ARRAY           ia;              /*!< the row indizes of the sparse matrix */
struct _ARRAY           blocks;          /*!< nodal block information of csr matrix 
                                              block[i][0] = size of nodal block
                                              block[i][1..numdf] = row-indizes in ia of row belonging
                                                                   to this block */
int                     firstcoupledof;  /*!< dof number of the first dof that has interproc coupling */

/* some arrays that are used for parallel assembly, mainly in the case of inter-proc-coupling conditions */
/* coupling is not supported by this format, but the data has to be here anyway to be compatibel with
   the assembly routines
*/
#ifdef PARALLEL 
int                     numcoupsend;     /*!< number of coupling information to be send by this proc */
int                     numcouprecv;     /*!< number of coupling information to be recv. by this proc */
struct _ARRAY          *couple_d_send;   /*!< send and receive buffers if necessary */
struct _ARRAY          *couple_i_send;
struct _ARRAY          *couple_d_recv;
struct _ARRAY          *couple_i_recv;
#endif
} DBCSR_ROOT;



/*!------------------------------------------------------------------------
\brief matrix needed by the MLPCG solver on each virtual grid

m.gee 6/01  

matrix needed by the MLPCG solver on each virtual grid       

-------------------------------------------------------------------------*/
typedef struct _DBCSR
{
int                     numeq_total;     /*!< total number of unknowns */
int                     numeq;           /*!< number of unknowns updated on this proc */ 
} DBCSR;




/*!------------------------------------------------------------------------
\brief one level of the multilevel preconditioner

m.gee 6/01  

one level of the multilevel preconditioner       

-------------------------------------------------------------------------*/
typedef struct _MLLEVEL
{
struct _DBCSR           csr;             /*!< the sparse matrix of this level */
} MLLEVEL;




/*!------------------------------------------------------------------------
\brief the ml preconditioner

m.gee 6/01  

matrix needed by the MLPCG solver on each virtual grid       

-------------------------------------------------------------------------*/
typedef struct _MLPRECOND
{
int                     numlev;
struct _MLLEVEL        *levels;          /*!< vector of levels */
} MLPRECOND;

























/*! @} (documentation module close)*/
