/*----------------------------------------------------------------------*
 | includes for solver package AZTEC 2.1                 m.gee 10/01    |
 *----------------------------------------------------------------------*/
#ifdef AZTEC_PACKAGE
#ifdef PARALLEL 
/*-------------------------------- with mpi parallel version of aztec2.1*/
#ifndef SUN
#include </bau/stat33/users/statik/lib/AZTEC21_MPI/az_aztec.h>
#else
#include "../../../lib_sun/aztec21/lib/az_aztec.h"
#endif
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
 | includes for solver package Spooles                    m.gee 4/02    |
 *----------------------------------------------------------------------*/
#ifdef PARALLEL 
#ifdef SPOOLES_PACKAGE
#ifndef SUN
#include "/bau/stat16/users/statik/lib/spooles/MPI/spoolesMPI.h"
#else
#include "../../../lib_sun/spooles/MPI/spoolesMPI.h"
#endif
#endif
#endif

/*----------------------------------------------------------------------*
 | includes for solver package Spooles              s.offermanns 5/02    |
 *----------------------------------------------------------------------*/
#ifdef UMFPACK
#include "/bau/stat16/users/statik/lib/UMFPACK/include/umfpack.h"
#endif

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
     mds,                    /* mlib - direct - sparse (sym.&nonsym)*/
     msr,                    /* distributed modified sparse row format */
     vbr,                    /* distributed variable block row format */
     parcsr,                 /* distributed compressed sparse row format */
     ucchb,                  /* unsymmetric column compressed Harwell-Boeing format */
     dense,                  /* dense matrix for Lapack */
     rc_ptr,                 /* row column pointer format for mumps */
     skymatrix,              /* skyline format for solver colsol */  
     spoolmatrix,            /* matrix object for solver spooles */
     ccf,                    /* compressed column format for umfpack */
} SPARSE_TYP;

/*----------------------------------------------------------------------*
 | a union which holds all types of sparse matrices      m.gee 11/01    |
 *----------------------------------------------------------------------*/
typedef union _SPARSE_ARRAY
{

struct _ML_ARRAY_MDS   *mds;    /* mlib - symm. - unsymm. - sparse */
struct _AZ_ARRAY_MSR   *msr;    /* pointer to Aztec's DMSR matrix */
struct _AZ_ARRAY_VBR   *vbr;    /* pointer to Aztec's DVBR matrix */
struct _H_PARCSR       *parcsr; /*         to HYPRE's ParCSR matrix */
struct _UCCHB          *ucchb;  /*         to Superlu's UCCHB matrix */
struct _DENSE          *dense;  /*         to dense matrix */
struct _RC_PTR         *rc_ptr; /*         to Mump's row/column ptr matrix */
struct _CCF            *ccf;    /*         to Umfpack compressed column matrix */
struct _SKYMATRIX      *sky;    /*         to Colsol's skyline matrix */
struct _SPOOLMAT       *spo;    /*         to Spoole's matrix */

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



struct _MLVAR          *mlvar;             /* variables needed for hp's mlib solver */
struct _AZVAR          *azvar;             /* variables needed for solver aztec */
struct _HYPREVARS      *hyprevar;          /* variables needed for HYPRE EuclidCG */
struct _PSUPERLUVARS   *psuperluvars;      /* variables needed for Parallel SuperLU */
struct _LAPACKVARS     *lapackvars;        /* variables needed for Lapack */
struct _MUMPSVARS      *mumpsvars;         /* variables needed for MUMPS */
struct _COLSOLVARS     *colsolvars;        /* variables needed for colsol */

int                     nsysarray;         /* number of global sparse arrays for this field */   
enum  _SPARSE_TYP      *sysarray_typ;      /* vector of types for all sparse arrays */
union _SPARSE_ARRAY    *sysarray;          /* vector of sparse arrays */

int                     nrhs;              /* number of distributed rhs-vectors */
struct _DIST_VECTOR    *rhs;               /* vector of distributed rhs-vectors */
int                     nsol;              /* number of distributed solution vectors */
struct _DIST_VECTOR    *sol;               /* vector of dist. solution vectors */
} SOLVAR;

/*----------------------------------------------------------------------*
 | variables needed for solver mlib                         al 12/01    |
 *----------------------------------------------------------------------*/
typedef struct _MLVAR
{
int                     symm  ;/* 1 -> symmetric , 0 -> nonsymm.       */
int                     msglvl;/* 0..4 -> 4 complete debugging output  */
int                     maxzer;/* additional fill in ,= 0 no fill in   */
int                     order;

double                  pvttol;/* 0.0 reorder with minim. fill in      */
                               /* 1.0 best numerical stability         */
} MLVAR;

/*----------------------------------------------------------------------*
 | variables needed for solver colsol                    m.gee 01/02    |
 *----------------------------------------------------------------------*/
typedef struct _COLSOLVARS
{
int                     i;                   /* this is in progress.... */
} COLSOLVARS;


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
int                     blockdiag;
} AZVAR;

/*----------------------------------------------------------------------*
 | variables needed for HYPRE solver package             m.gee 10/01    |
 *----------------------------------------------------------------------*/
typedef struct _HYPREVARS
{
enum _HYPREPRECTYP      hypre_prectyp;    /* type of hypre preconditioner */
int                     io;               /* flag to set solver quiet */
int                     maxiter;          /* max iterations allowed */
int                     numiter;          /* number of iterations taken */
double                  resnorm;          /* residual norm achieved */
int                     reuse;            /* reuse feature (not yet impl.) */
double                  tol;              /* user-given tolerance */
int                     kryldim;          /* dimension of krylov subspace */
double                  threshold;        /* parameters for amg, see manual */
int                     sweep[4];         
int                     ifill;            /* fill in level for ilu */
double                  dfill;            /* fill in level in percent of the original matrix for ilu and parasails */
int                     bj;               /* ? */

int                     parasymm;         /* parasails preconditioner parameters */
int                     paralevel;
double                  parathresh;
double                  parafilter;
} HYPREVARS;

/*----------------------------------------------------------------------*
 | variables needed for solver ParSuperLU                m.gee 10/01    |
 *----------------------------------------------------------------------*/
typedef struct _PSUPERLUVARS
{
int                     reuse;            /* in progress.... */
} PSUPERLUVARS;

/*----------------------------------------------------------------------*
 | variables needed for solver Lapack                    m.gee 11/01    |
 *----------------------------------------------------------------------*/
typedef struct _LAPACKVARS
{
int                     reuse;            /* in progress.... */
} LAPACKVARS;

/*----------------------------------------------------------------------*
 | a skyline matrix                                      m.gee 01/02    |
 | this structure holds a skyline matrix to be solved with colsol       |
 *----------------------------------------------------------------------*/
typedef struct _SKYMATRIX
{
int                     is_init;         /* was this matrix initialized ? */
int                     is_factored;     /* is this matrix already factored ? */
int                     ncall;           /* how often was this matrix solved */

int                     numeq_total;     /* total number of unknowns */
int                     numeq;           /* number of unknowns updated on this proc */ 
int                     nnz_total;       /* total number of nonzero entries */
int                     nnz;             /* number of nonzeros on this proc */

struct _ARRAY           update;          /* sorted list of dofs updated on this proc */
struct _ARRAY           maxa;
struct _ARRAY           maxaf;
struct _ARRAY           A;
double                  det;
/* some arrays that are used for parallel assembly, mainly in the case of inter-proc-coupling conditions */
#ifdef PARALLEL 
int                     numcoupsend;     /* number of coupling information to be send by this proc */
int                     numcouprecv;     /* number of coupling information to be recv. by this proc */
struct _ARRAY          *couple_d_send;   /* send and receive buffers if necessary */
struct _ARRAY          *couple_i_send;
struct _ARRAY          *couple_d_recv;
struct _ARRAY          *couple_i_recv;
#endif
} SKYMATRIX;

/*----------------------------------------------------------------------*
 | a sparse matrix in row/column pointer format          m.gee 12/02    |
 | this structure holds a distributed sparse matrix for Mumps           |
 | it uses two integer vectors irn_loc,jcn_loc to hold indicees of an   |
 | entry in A_loc (see Mumps manual)
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

int                     icntl[20];
int                     comm;

struct _ARRAY           update;          /* sorted list of dofs updated on this proc */
struct _ARRAY           irn_loc;         /* proc-local row pointer vector */
struct _ARRAY           irn_locf;        /* fortran style pointer vector of irn_loc */
struct _ARRAY           jcn_loc;         /* proc-local column pointer vector */
struct _ARRAY           jcn_locf;        /* fortran style pointer vector of jcn_loc */
struct _ARRAY           A_loc;           /* values of the matrix */
struct _ARRAY           rowptr;          /* int vector holding the begin of each row in irn_loc */

struct _ARRAY           irn_glob;        /* on imyrank=0 the global row/column pointer arrays */
struct _ARRAY           jcn_glob;

/* some arrays that are used for parallel assembly, mainly in the case of inter-proc-coupling conditions */
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
 | a sparse matrix in compressed column format      s.offermanns 04/02  |
 | this structure holds a distributed sparse matrix for Umfpack 4beta   |
 | it uses two integer vectors Ap, Ap to hold indicees of an            |
 | entry in Ax (see Umpfack manual)
 *----------------------------------------------------------------------*/
typedef struct _CCF
{
int                     is_init;         /* was this matrix initialized ? */
int                     is_factored;     /* is this matrix already factored ? */
int                     ncall;           /* how often was this matrix solved */

int                     numeq_total;     /* total number of unknowns */
int                     numeq;           /* number of unknowns updated on this proc */ 
int                     nnz_total;       /* total number of nonzero entries */
int                     nnz;             /* number of nonzeros on this proc */

int                     comm;

struct _ARRAY           update;          /* sorted list of dofs updated on this proc */
struct _ARRAY           Ap;              /* column pointer vector */
struct _ARRAY           Ai;              /* row pointer vector */
struct _ARRAY           Ax;              /* values of the matrix */

/* some arrays that are used for parallel assembly, mainly in the case of inter-proc-coupling conditions */
#ifdef PARALLEL 
int                     numcoupsend;     /* number of coupling information to be send by this proc */
int                     numcouprecv;     /* number of coupling information to be recv. by this proc */
struct _ARRAY          *couple_d_send;   /* send and receive buffers if necessary */
struct _ARRAY          *couple_i_send;
struct _ARRAY          *couple_d_recv;
struct _ARRAY          *couple_i_recv;
#endif
} CCF;

/*----------------------------------------------------------------------*
 | a dense matrix                                        m.gee 11/01    |
 | this structure holds a dense matrix to be solved with lapack         |
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


/* some arrays that are used for parallel assembly, mainly in the case of inter-proc-coupling conditions */
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
 | unsym compressed harwell Boeing matrix, to be used with SuperLU
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
struct _ARRAY           asub;             /* pointer vector of the ucchb */
struct _ARRAY           asub_backup;      /* backup of pointer vector of the ucchb */
struct _ARRAY           asub_perm_backup; /* backup of permuted pointer vector of the ucchb */
struct _ARRAY           xa;               /* pointer vector of the ucchb */
struct _ARRAY           xa_backup;        /* backup of pointer vector of the ucchb */
struct _ARRAY           xa_perm_backup;   /* permuted backup of pointer vector of the ucchb */

#ifdef PARALLEL 
#ifdef PARSUPERLU_PACKAGE
gridinfo_t              grid;             /* the 2D MPI-grid to solve on */
superlu_options_t       options;          /* options for superLU */
SuperLUStat_t           stat;             /* statistics */
ScalePermstruct_t       ScalePerm;        /* Scaling and permutation information */
LUstruct_t              LUstruct;         /* structure for the L & U decomposition */
SuperMatrix             A;                /* the matrix structure */
#endif

/* some arrays that are used for parallel assembly, mainly in the case of inter-proc-coupling conditions */
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
 | matrix in distributed compressed sparse row format (PCSR)for HYPRE   |
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
HYPRE_IJMatrix          ij_matrix;           /* the matrix itself, hidden away by Hypre*/
HYPRE_Solver            solver;              /* hypre solver structure */
HYPRE_Solver            precond;             /* hypre preconditioner structure */
#endif

/* some arrays that are used for parallel assembly, mainly in the case of inter-proc-coupling conditions */
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
 | matrix in dist. modified sparse row format (DMSR) for Aztec2.1       |
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
struct _ARRAY           bindx_backup;     /* backup of bindx, as bindx is altered by solver */
struct _ARRAY           val;              /* values of matrix */
struct _ARRAY           val_backup;       /* backup of values of matrix as val is altered by solver */

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

/* some arrays that are used for parallel assembly, mainly in the case of inter-proc-coupling conditions */
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
 | a DMSR Matrix (Distributed Modified Sparse Row format) m.gee 5/01    |
 | matrix in dist. modified sparse row format (DMSR) for Aztec2.1       |
 *----------------------------------------------------------------------*/
typedef struct _AZ_ARRAY_VBR
{
int                     is_init;          /* was this matrix initialized ? */
int                     is_factored;      /* does precond. information exist ? */
int                     ncall;            /* how often was this matrix solved */
int                     numeq_total;      /* total number of unknowns */
int                     numeq;            /* number of unknowns updated on this proc */ 
int                     nnz;              /* number of nonzeros on this proc */

struct _ARRAY           update;           /* list of dofs updated on this proc */
struct _ARRAY           bupdate;          /* list of nodal blocks updated on this proc */
   /* vector of block connectivity */
   /* block_connect[actnodeId][0] = my own block size */
   /* block_connect[actnodeId][1] = number off-diagonal blocks in this row */
   /* block_connect[actnodeId][2..] = id_loc of off-diagonal blocks */
int                     block_consize;
int                   **block_connect;
struct _ARRAY           rpntr;
struct _ARRAY           cpntr;
struct _ARRAY           bpntr;
struct _ARRAY           bindx;
struct _ARRAY           indx;
struct _ARRAY           val;

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

/* some arrays that are used for parallel assembly, mainly in the case of inter-proc-coupling conditions */
#ifdef PARALLEL 
int                     numcoupsend;     /* number of coupling information to be send by this proc */
int                     numcouprecv;     /* number of coupling information to be recv. by this proc */
struct _ARRAY          *couple_d_send;   /* send and receive buffers if necessary */
struct _ARRAY          *couple_i_send;
struct _ARRAY          *couple_d_recv;
struct _ARRAY          *couple_i_recv;
#endif
} AZ_ARRAY_VBR;




/*----------------------------------------------------------------------*
 | a Spoole Matrix                                        m.gee 4/03    |
 *----------------------------------------------------------------------*/
typedef struct _SPOOLMAT
{
int                     is_init;          /* was this matrix initialized ? */
int                     is_factored;      /* does precond. information exist ? */
int                     ncall;            /* how often was this matrix solved */
int                     numeq_total;      /* total number of unknowns */
int                     numeq;            /* number of unknowns updated on this proc */ 
int                     nnz;              /* number of nonzeros on this proc */

struct _ARRAY           update;          /* list of dofs updated on this proc */
struct _ARRAY           irn_loc;         /* proc-local row pointer vector */
struct _ARRAY           jcn_loc;         /* proc-local column pointer vector */
struct _ARRAY           rowptr;          /* int vector holding the begin of each row in irn_loc */
struct _ARRAY           A_loc;           /* the values of the matrix */

#ifdef SPOOLES_PACKAGE
FrontMtx               *frontmtx;
InpMtx                 *mtxA;            /* the sparse matrix object */
DenseMtx               *mtxY;
DenseMtx               *mtxX;
Graph                  *graph;
IVL                    *adjIVL;
ETree                  *frontETree;
IVL                    *symbfacIVL;
SubMtxManager          *mtxmanager;
SubMtxManager          *solvemanager;
ChvManager             *chvmanager;
Chv                    *rootchv;
IV                     *oldToNewIV;
IV                     *newToOldIV;
IV                     *ownersIV;
IV                     *vtxmapIV;
DV                     *cumopsDV;
IV                     *ownedColumnsIV;
InpMtx                 *newA;    
DenseMtx               *newY;
SolveMap               *solvemap;
#endif

/* some arrays that are used for parallel assembly, mainly in the case of inter-proc-coupling conditions */
#ifdef PARALLEL 
int                     numcoupsend;     /* number of coupling information to be send by this proc */
int                     numcouprecv;     /* number of coupling information to be recv. by this proc */
struct _ARRAY          *couple_d_send;   /* send and receive buffers if necessary */
struct _ARRAY          *couple_i_send;
struct _ARRAY          *couple_d_recv;
struct _ARRAY          *couple_i_recv;
#endif
} SPOOLMAT;




/*----------------------------------------------------------------------*
 |  column pointer, row index sparse matrix representation  al  10/01   |
 |  for the lower triangle of the matrix                                |
 |                                                                      |
 |                    - HP's MLIB -                                     |
 *----------------------------------------------------------------------*/
typedef struct _ML_ARRAY_MDS
{
char               arrayname[50];
int                is_init;

/* input */

int                numeq;       /* number of equations                  */
int                nnz;         /* number of nonzeroes                  */
int                output;      /* fortran unit number =6 -> screen     */

double             cond;        /* condition number                     */

struct _ARRAY      colstr;      /* gives the index in rowind of the
                                   first nonzero in the lower triangular
                                   part of column j of the matrix       */
struct _ARRAY      rowind;      /* list of row indices for all nonzeros
                                   within each column                   */
 
/* output */

double             rcond;       /* estimate the reciprocal of the l-norm
                                   condition number                     */

int                inrtia[3];   /* number of positive, negative and 
                                   an indicator if there are zero
                                   eigenvalues                          */
                                   
double             global[150]; /* global communication array           */

int                ierr;        /* = 0; normal return                   */                                   
                                                   
} ML_ARRAY_MDS;




/*----------------------------------------------------------------------*
 | a distributed vector for solution                      m.gee 6/01    |
 | a vector distributed among processors.                               |
 | each processor holds a piece of size numeq of the vector, the total  |
 | size of vector is numeq_total                                        |
 | the layout of the vector in general suits the data format of one of  |
 | the sparse matrix formats above                                      |
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
 | the global variable solv and the head solution.h is not visible      |
 | from everywhere in the code, so the prototypes of all functions using|
 | structurs declared in solution.h are in here and not in prototypes.h |
 |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
#include "prototypes_sol.h"
/*----------------------------------------------------------------------*
 | define the global structure pointer *solv                            |
 |                                                                      |
 | global variable *solv, vector of lenght numfld of structures SOLVAR  |
 | is defined in solver_control.c
 |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
