/*----------------------------------------------------------------------*
 | includes for solver package AZTEC 2.1                 m.gee 10/01    |
 *----------------------------------------------------------------------*/
#ifdef AZTEC_PACKAGE
#ifdef PARALLEL 
/*-------------------------------- with mpi parallel version of aztec2.1*/
#include </bau/stat33/users/statik/lib/AZTEC21_MPI/az_aztec.h>
#else
/*------------------------ without mpi , sequentiel version of aztec2.1 */
#include </bau/stat33/users/statik/lib/AZTEC21/az_aztec.h>
#endif
#endif /* end of ifdef AZTEC_PACKAGE */

/*----------------------------------------------------------------------*
 | includes for solver package HYPRE 1.6.0               m.gee 10/01    |
 *----------------------------------------------------------------------*/
#ifdef PARALLEL/*---------------- HYPRE exists only in parallel version */ 
#ifdef HYPRE_PACKAGE
#include "utilities.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_mv.h"
#include "HYPRE_IJ_mv.h"
#include "HYPRE_parcsr_ls.h"
#include "krylov.h"
#endif /* end of ifdef HYPRE_PACKAGE */
#endif/* end of ifdef PARALLEL */

/*----------------------------------------------------------------------*
 | includes for solver package SuperLU_DIST              m.gee 10/01    |
 *----------------------------------------------------------------------*/
#ifdef PARALLEL/*---------------- HYPRE exists only in parallel version */ 
#ifdef PARSUPERLU_PACKAGE
#include "superlu_ddefs.h"
#endif /* end of ifdef PARSUPERLU_PACKAGE */
#endif/* end of ifdef PARALLEL */

/*----------------------------------------------------------------------*
 | an enum which describes all types of sparse matrices   m.gee 5/01    |
 *----------------------------------------------------------------------*/
typedef enum _SPARSE_TYP
{
     sparse_none,            /* sparse matrix type not specified */
     msr,                    /* distributed modified sparse row format */
     parcsr,                 /* distributed compressed sparse row format */
     ucchb,                  /* unsymmetric column compressed Harwell-Boeing format */
     dense,                  /* dense matrix for Lapack */
     rc_ptr                  /* row column pointer format for mumps */
} SPARSE_TYP;

/*----------------------------------------------------------------------*
 | a union which holds all types of sparse matrices      m.gee 11/01    |
 *----------------------------------------------------------------------*/
typedef union _SPARSE_ARRAY
{

struct _AZ_ARRAY_MSR   *msr;    /* pointer to Aztec's DMSR matrix */
struct _H_PARCSR       *parcsr; /*         to HYPRE's ParCSR matrix */
struct _UCCHB          *ucchb;  /*         to Superlu's UCCHB matrix */
struct _DENSE          *dense;  /*         to dense matrix */
struct _RC_PTR         *rc_ptr; /*         to Mump's row/column ptr matrix */

} SPARSE_ARRAY;

/*----------------------------------------------------------------------*
 | distr. sparse matrices, vectors and general solver data    m.gee 5/01|
 | this is the main structure used by all types of solvers              |
 *----------------------------------------------------------------------*/
typedef struct _SOLVAR
{
enum   _FIELDTYP        fieldtyp;          /* type of field */
enum   _PART_TYP        parttyp;           /* typ of partition */
enum   _SOLVER_TYP      solvertyp;         /* typ of chosen solver */

struct _AZVAR          *azvar;             /* variables needed for solver aztec */
struct _HYPREVARS      *hyprevar;          /* variables needed for HYPRE EuclidCG */
struct _PSUPERLUVARS   *psuperluvars;      /* variables needed for Parallel SuperLU */
struct _LAPACKVARS     *lapackvars;        /* variables needed for Lapack */
struct _MUMPSVARS      *mumpsvars;         /* variables needed for MUMPS */

int                     nsysarray;         /* number of global sparse arrays for this field */   
enum  _SPARSE_TYP      *sysarray_typ;      /* vector of types for all sparse arrays */
union _SPARSE_ARRAY    *sysarray;          /* vector of sparse arrays */

int                     nrhs;              /* number of distributed rhs-vectors */
struct _DIST_VECTOR    *rhs;               /* vector of distributed rhs-vectors */
int                     nsol;              /* number of distributed solution vectors */
struct _DIST_VECTOR    *sol;               /* vector of dist. solution vectors */
} SOLVAR;

/*----------------------------------------------------------------------*
 | variables needed for solver MUMPS                     m.gee 01/02    |
 *----------------------------------------------------------------------*/
typedef struct _MUMPSVARS
{
int                     i;                   /* this is in progress.... */
} MUMPSVARS;

/*----------------------------------------------------------------------*
 | variables needed for solver aztec                     m.gee 10/01    |
 *----------------------------------------------------------------------*/
typedef struct _AZVAR
{
enum   _AZSOLVERTYP     azsolvertyp;        /* subtype of aztec solver, see enums.h */
enum   _AZPRECTYP       azprectyp;          /* type of aztec preconditioner, see enums.h */
int                     azreuse;            /* reuse of preconditioning feature, important, but not yet implemented */
int                     azgfill;            /* percentage fill in allowed */
int                     aziter;             /* maximum number of iterations allowed */
int                     azsub;              /* number of krylov subspaces for vertain solvers (e.g. gmres) */
int                     azgraph;            /* forgot about it.... */
int                     azpoly;             /* integer parameter with different meanings dependent on type of precond. */
double                  azdrop;             /* numerical drop tolerance for preconditioners using it, default: 0.0 */
double                  azfill;             /* allowed fill-in in percent of the memory used by the sparse matrix */             
double                  aztol;              /* tolerance */
double                  azomega;            /* relaxation parameter for some preconditioners */
} AZVAR;

/*----------------------------------------------------------------------*
 | variables needed for HYPRE solver package             m.gee 10/01    |
 *----------------------------------------------------------------------*/
typedef struct _HYPREVARS
{
enum _HYPREPRECTYP      hypre_prectyp;    /* type of hypre preconditioner */
int                     io;               /* flag to set solver quiet */
int                     maxiter;          /* max iterations allowed */
int                     numiter;
double                  resnorm;
int                     reuse;
double                  tol;
int                     kryldim;
double                  threshold;
int                     sweep[4];
int                     ifill;
double                  dfill;
int                     bj;

int                     parasymm;
int                     paralevel;
double                  parathresh;
double                  parafilter;
} HYPREVARS;

/*----------------------------------------------------------------------*
 | variables needed for solver ParSuperLU                m.gee 10/01    |
 *----------------------------------------------------------------------*/
typedef struct _PSUPERLUVARS
{
int                     reuse;
} PSUPERLUVARS;

/*----------------------------------------------------------------------*
 | variables needed for solver Lapack                    m.gee 11/01    |
 *----------------------------------------------------------------------*/
typedef struct _LAPACKVARS
{
int                     reuse;
} LAPACKVARS;

/*----------------------------------------------------------------------*
 | a sparse matrix in row/column pointer format          m.gee 12/02    |
 *----------------------------------------------------------------------*/
typedef struct _RC_PTR
{
int                     is_init;         /* was this matrix initialized ? */
int                     is_factored;     /* is this matrix already factored ? */
int                     ncall;           /* how often was this matrix solved */

int                     numeq_total;     /* total number of unknowns */
int                     numeq;           /* number of unknowns updated on this proc */ 
int                     nnz_total;       /* total number of nonzero entries */
int                     nnz;             /* number of nonzeros on this proc */

struct _ARRAY           update;          /* sorted list of dofs updated on this proc */
/*struct _ARRAY           irn;
struct _ARRAY           jcn;*/
struct _ARRAY           irn_loc;
struct _ARRAY           jcn_loc;
struct _ARRAY           A_loc;
struct _ARRAY           rowptr;
struct _ARRAY           bindx;

#ifdef PARALLEL 
int                     numcoupsend;     /* number of coupling information to be send by this proc */
int                     numcouprecv;     /* number of coupling information to be recv. by this proc */
struct _ARRAY          *couple_d_send;   /* send and receive buffers if necessary */
struct _ARRAY          *couple_i_send;
struct _ARRAY          *couple_d_recv;
struct _ARRAY          *couple_i_recv;
#endif
} RC_PTR;

/*----------------------------------------------------------------------*
 | a dense matrix                                        m.gee 11/01    |
 *----------------------------------------------------------------------*/
typedef struct _DENSE
{
int                     is_init;         /* was this matrix initialized ? */
int                     is_factored;     /* is this matrix already factored ? */
int                     ncall;           /* how often was this matrix solved */

int                     numeq_total;     /* total number of unknowns */
int                     numeq;           /* number of unknowns updated on this proc */ 
int                     nnz_total;       /* total number of nonzero entries */
int                     nnz;             /* number of nonzeros on this proc */

struct _ARRAY           update;          /* sorted list of dofs updated on this proc */
struct _ARRAY           A;               /* the dense matrix */
struct _ARRAY           ipiv;            /* pivoting information */
int                     lwork;           /* work space for sym. Lapack solver */
struct _ARRAY           work;            /* work space for sym. Lapack solver */


#ifdef PARALLEL 
int                     numcoupsend;     /* number of coupling information to be send by this proc */
int                     numcouprecv;     /* number of coupling information to be recv. by this proc */
struct _ARRAY          *couple_d_send;   /* send and receive buffers if necessary */
struct _ARRAY          *couple_i_send;
struct _ARRAY          *couple_d_recv;
struct _ARRAY          *couple_i_recv;
#endif
} DENSE;

/*----------------------------------------------------------------------*
 | a uccHB matrix                                        m.gee 11/01    |
 *----------------------------------------------------------------------*/
typedef struct _UCCHB
{
int                     is_init;          /* was this matrix initialized ? */
int                     is_factored;      /* is this matrix already factored ? */
int                     ncall;            /* how often was this matrix solved */

int                     numeq_total;      /* total number of unknowns */
int                     numeq;            /* number of unknowns updated on this proc */ 
int                     nnz_total;        /* total number of nonzero entries */
int                     nnz;              /* number of nonzeros on this proc */

struct _ARRAY           update;           /* list of dofs updated on this proc */
struct _ARRAY           a;                /* the ucchb matrix */
struct _ARRAY           asub;
struct _ARRAY           asub_backup;
struct _ARRAY           asub_perm_backup;
struct _ARRAY           xa;
struct _ARRAY           xa_backup;
struct _ARRAY           xa_perm_backup;

#ifdef PARALLEL 
#ifdef PARSUPERLU_PACKAGE
gridinfo_t              grid;             /* the 2D MPI-grid to solve on */
superlu_options_t       options;          /* options for superLU */
SuperLUStat_t           stat;             /* statistics */
ScalePermstruct_t       ScalePerm;        /* Scaling and permutation information */
LUstruct_t              LUstruct;         /* structure for the L & U decomposition */
SuperMatrix             A;                /* the matrix structure */
#endif

int                     numcoupsend;     /* number of coupling information to be send by this proc */
int                     numcouprecv;     /* number of coupling information to be recv. by this proc */
struct _ARRAY          *couple_d_send;   /* send and receive buffers if necessary */
struct _ARRAY          *couple_i_send;
struct _ARRAY          *couple_d_recv;
struct _ARRAY          *couple_i_recv;
#endif
} UCCHB;

/*----------------------------------------------------------------------*
 | a ParCSR matrix                                        m.gee 5/01    |
 *----------------------------------------------------------------------*/
typedef struct _H_PARCSR
{
int                     is_init;          /* was this matrix initialized ? */
int                     is_factored;      /* does precond. information exist ? */
int                     ncall;            /* how often was this matrix solved */
int                     numeq_total;      /* total number of unknowns */
int                     numeq;            /* number of unknowns updated on this proc */ 
int                     nnz;              /* number of nonzeros on this proc */

struct _ARRAY           perm;                /* permutation of update for each proc, type is IA*/
struct _ARRAY           perm_sizes;          /* size of perm on each proc, type is IV */
struct _ARRAY           update;              /* ascending list of dofs on all procs, type is IA*/
struct _ARRAY           bindx;               /* see documenation for DMSR format (Aztec2.1) */

#ifdef HYPRE_PACKAGE
HYPRE_IJMatrix          ij_matrix;           /* the matrix itself, hidden away */
HYPRE_Solver            solver;              /* hypre solver structure */
HYPRE_Solver            precond;             /* hypre preconditioner structure */
#endif

#ifdef PARALLEL 
int                     numcoupsend;     /* number of coupling information to be send by this proc */
int                     numcouprecv;     /* number of coupling information to be recv. by this proc */
struct _ARRAY          *couple_d_send;   /* send and receive buffers if necessary */
struct _ARRAY          *couple_i_send;
struct _ARRAY          *couple_d_recv;
struct _ARRAY          *couple_i_recv;
#endif
} H_PARCSR;

/*----------------------------------------------------------------------*
 | a DMSR Matrix (Distributed Modified Sparse Row format) m.gee 5/01    |
 *----------------------------------------------------------------------*/
typedef struct _AZ_ARRAY_MSR
{
int                     is_init;          /* was this matrix initialized ? */
int                     is_factored;      /* does precond. information exist ? */
int                     ncall;            /* how often was this matrix solved */
int                     numeq_total;      /* total number of unknowns */
int                     numeq;            /* number of unknowns updated on this proc */ 
int                     nnz;              /* number of nonzeros on this proc */

struct _ARRAY           update;           /* list of dofs updated on this proc */
int                     shift;            /* binary shift for searching in update */
int                    *bins;             /* binary mirror of update */
struct _ARRAY           bindx;            /* the sparse matrix */
struct _ARRAY           bindx_backup;
struct _ARRAY           val;
struct _ARRAY           val_backup;

#ifdef AZTEC_PACKAGE
double                  params[AZ_PARAMS_SIZE];    /* Aztec parameters */
double                  status[AZ_STATUS_SIZE];    /* Aztec return status */
int                     proc_config[AZ_PROC_SIZE]; /* MPI-configuration */
int                     options[AZ_OPTIONS_SIZE];  /* Aztec options */
int                    *data_org;                  /* Aztec internal data org. */
int                    *external;                  /* list of external dofs often needed by this proc */
int                    *update_index;              /* list of dofs updated on this proc */
int                    *extern_index;              /* list of external related dofs */
int                     N_external;                /* number of external dofs often needed by this proc */

AZ_MATRIX              *Amat;                      /* the matrix object */
AZ_PRECOND             *Aprec;                     /* the preconditioner object */
#endif

#ifdef PARALLEL 
int                     numcoupsend;     /* number of coupling information to be send by this proc */
int                     numcouprecv;     /* number of coupling information to be recv. by this proc */
struct _ARRAY          *couple_d_send;   /* send and receive buffers if necessary */
struct _ARRAY          *couple_i_send;
struct _ARRAY          *couple_d_recv;
struct _ARRAY          *couple_i_recv;
#endif
} AZ_ARRAY_MSR;

/*----------------------------------------------------------------------*
 | a distributed vector for solution                      m.gee 6/01    |
 *----------------------------------------------------------------------*/
typedef struct _DIST_VECTOR
{
int                     numeq_total;     /* total size of distr. vector */
int                     numeq;           /* local size of distr. vector */
struct _ARRAY           vec;             /* local piece of distr. vector */
} DIST_VECTOR;

/*----------------------------------------------------------------------*
 | put the prototypes of all                                            |
 | routines that have contact to solution.h in here                     |
 |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
#include "prototypes_sol.h"
/*----------------------------------------------------------------------*
 | define the global structure solv                                     |
 |                                                                      |
 |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
SOLVAR            *solv;
