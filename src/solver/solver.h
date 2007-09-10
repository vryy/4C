/*!---------------------------------------------------------------------
\file
\brief all data types for solver packages

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

---------------------------------------------------------------------*/

#ifndef SOLVER_H
#define SOLVER_H

/*----------------------------------------------------------------------*
 | includes for solver package AZTEC 2.1                 m.gee 10/01    |
 *----------------------------------------------------------------------*/
#ifdef AZTEC_PACKAGE
#ifdef TRILINOS_PACKAGE
#include <az_aztec.h>
#else
#include <aztec21/lib/az_aztec.h>
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

/* Avoid the stupid spooles timing header. It's never used anyway. */
#define _TIMINGS_

#ifdef HPUX10
#include <spooles/MPI/spoolesMPI.h>
#endif

#ifdef HPUX11
#include <spooles/MPI/spoolesMPI.h>
#endif

#ifdef HPUXITA
#include <spooles/MPI/spoolesMPI.h>
#endif

#if defined(LINUX) || defined(LINUX64)
#include <spooles/MPI/spoolesMPI.h>
#endif

#ifdef SX6
#include <spooles/MPI/spoolesMPI.h>
#endif

#ifdef SX8
#include <spooles/MPI/spoolesMPI.h>
#endif

#ifdef TX7
#include <spooles/MPI/spoolesMPI.h>
#endif

#ifdef AZUSA
#include "../../../../lib_ita1/spooles/MPI/spoolesMPI.h"
#endif

#ifdef SUN
#include "../../../lib_sun/spooles/MPI/spoolesMPI.h"
#endif

#ifdef LINUX_MUENCH
#include <spooles/MPI/spoolesMPI.h>
#endif

#ifdef HPUX_MUENCH
#include <spooles/MPI/spoolesMPI.h>
#endif

#endif
#endif

/*----------------------------------------------------------------------*
 | includes for solver package UMFPACK                      mn 05/02    |
 *----------------------------------------------------------------------*/
#ifdef UMFPACK

#if defined(LINUX) || defined(LINUX64)
#include <umfpack/umfpack.h>
#elif defined(LINUX_MUENCH)
#include <umfpack.h>
#elif defined(HPUX_MUENCH)
#include <umfpack.h>
#elif defined(WIN_MUENCH)
#include <umfpack.h>
#else
#include <umfpack/umfpack.h>
#endif

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
  parcsr,                 /* distributed compressed sparse row format */
  ucchb,                  /* unsymmetric column compressed Harwell-Boeing format */
  dense,                  /* dense matrix for Lapack */
  rc_ptr,                 /* row column pointer format for mumps */
  skymatrix,              /* skyline format for solver colsol */
  spoolmatrix,            /* matrix object for solver spooles */
  ccf,                    /* compressed column format for umfpack */
  bdcsr,                  /* block distributed compressed sparse row format */
  oll,                    /* orthogonal linked list matrix */
  trilinos                /* use Trilinos' Epetra_CrsMatrix matrix */
} SPARSE_TYP;


/*----------------------------------------------------------------------*
 | a union which holds all types of sparse matrices      m.gee 11/01    |
 *----------------------------------------------------------------------*/
typedef union _SPARSE_ARRAY
{
  struct _ML_ARRAY_MDS   *mds;      /* mlib - symm. - unsymm. - sparse */
  struct _AZ_ARRAY_MSR   *msr;      /* pointer to Aztec's DMSR matrix */
  struct _H_PARCSR       *parcsr;   /*         to HYPRE's ParCSR matrix */
  struct _UCCHB          *ucchb;    /*         to Superlu's UCCHB matrix */
  struct _DENSE          *dense;    /*         to dense matrix */
  struct _RC_PTR         *rc_ptr;   /*         to Mump's row/column ptr matrix */
  struct _CCF            *ccf;      /*         to Umfpack compressed column matrix */
  struct _SKYMATRIX      *sky;      /*         to Colsol's skyline matrix */
  struct _SPOOLMAT       *spo;      /*         to Spoole's matrix */
  struct _DBCSR          *bdcsr;    /* matrix needed by the MLPCG solver on the finest grid only */
  struct _OLL            *oll;      /* pointer to orthogonal linked list matrix */
  struct _TRILINOSMATRIX *trilinos; /* pointer to orthogonal linked list matrix */
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
  enum   _MATRIX_TYP      matrixtyp;         /* typ of chosen matrix */


  struct _MLVAR          *mlvar;             /* variables needed for hp's mlib solver */
  struct _AZVAR          *azvar;             /* variables needed for solver aztec */
  struct _HYPREVARS      *hyprevar;          /* variables needed for HYPRE EuclidCG */
  struct _PSUPERLUVARS   *psuperluvars;      /* variables needed for Parallel SuperLU */
  struct _LAPACKVARS     *lapackvars;        /* variables needed for Lapack */
  struct _MUMPSVARS      *mumpsvars;         /* variables needed for MUMPS */
  struct _COLSOLVARS     *colsolvars;        /* variables needed for colsol */

#ifdef MLPCG
  struct _MLPCGVARS      *mlpcgvars;         /* variables needed for MLPCG */
#endif

  INT                     nsysarray;         /* number of global sparse arrays for this field */
  enum  _SPARSE_TYP      *sysarray_typ;      /* vector of types for all sparse arrays */
  union _SPARSE_ARRAY    *sysarray;          /* vector of sparse arrays */

  INT                     nrhs;              /* number of distributed rhs-vectors */
  struct _DIST_VECTOR    *rhs;               /* vector of distributed rhs-vectors */
  INT                     nsol;              /* number of distributed solution vectors */
  struct _DIST_VECTOR    *sol;               /* vector of dist. solution vectors */
} SOLVAR;


/*----------------------------------------------------------------------*
 | variables needed for solver mlib                         al 12/01    |
 *----------------------------------------------------------------------*/
typedef struct _MLVAR
{
  INT                     symm  ;/* 1 -> symmetric , 0 -> nonsymm.       */
  INT                     msglvl;/* 0..4 -> 4 complete debugging output  */
  INT                     maxzer;/* additional fill in ,= 0 no fill in   */
  INT                     order;

  DOUBLE                  pvttol;/* 0.0 reorder with minim. fill in      */
  /* 1.0 best numerical stability         */
} MLVAR;


/*----------------------------------------------------------------------*
 | variables needed for solver colsol                    m.gee 01/02    |
 *----------------------------------------------------------------------*/
typedef struct _COLSOLVARS
{
  INT                     i;                   /* this is in progress.... */
} COLSOLVARS;


/*----------------------------------------------------------------------*
 | variables needed for solver MUMPS                     m.gee 01/02    |
 *----------------------------------------------------------------------*/
typedef struct _MUMPSVARS
{
  INT                     i;                   /* this is in progress.... */
} MUMPSVARS;


/*----------------------------------------------------------------------*
 | variables needed for solver aztec                     m.gee 10/01    |
 *----------------------------------------------------------------------*/
typedef struct _AZVAR
{
  enum   _AZSOLVERTYP     azsolvertyp;        /* subtype of aztec solver, see enums.h */
  enum   _AZPRECTYP       azprectyp;          /* type of aztec preconditioner, see enums.h */
  INT                     azreuse;            /* reuse of preconditioning feature, important,
                                                 but not yet implemented */
  INT                     azoutput;           /* output level for AztecOO 0=no output */
  INT                     azgfill;            /* percentage fill in allowed */
  INT                     aziter;             /* maximum number of iterations allowed */
  INT                     azsub;              /* number of krylov subspaces for certain solvers
                                                 (e.g. gmres) */
  INT                     azgraph;            /* forgot about it.... */
  INT                     azpoly;             /* integer parameter with different meanings
                                                 dependent on type of precond. */
  INT                     azoverlap;          /* amount of overlap for additive schwartz type preconditioners (e.g. ILU) */
  DOUBLE                  azdrop;             /* numerical drop tolerance for preconditioners
                                                 using it, default: 0.0 */
  DOUBLE                  azfill;             /* allowed fill-in in percent of the memory
                                                 used by the sparse matrix */
 INT                     azconv;             /*  Convergence check applied by
					      *  Aztec solver (for example
					      *  (||r|| unscaled,
					      *   ||r||/||r0||,
					      *   ||r||/||rhs||, ...)
					      *  The corresponding integer
					      *  numbers are defined in a
					      *  header and prescribed by the
					      *  AZTEC package */
  DOUBLE                  aztol;              /* tolerance */
  DOUBLE                  azomega;            /* relaxation parameter for some preconditioners */
  INT                     blockdiag;

#ifdef TRILINOS_PACKAGE
  INT                     azscal;             /* 0 = none, 1 = sym, 2 = infnorm */

  int                     mlprint;            /* ml print level 0 - 10 */
  int                     mlcsize;            /* size where to stop coarsening */
  int                     mlmaxlevel;         /* max no. of grids */
  int                     mlsmotimes[15];     /* no. smoothing steps on each level */
  int                     mlcoarsentype;      /* 0 UC 1 METIS 2 VBMETIS 3 MIS */
  int                     mlaggsize;
  int                     mlsmotype_fine;     /* 0 SGS 1 Jacobi 2 Chebychev 3 MLS 4 ILU 5 KLU */
  int                     mlsmotype_med;      /* 0 SGS 1 Jacobi 2 Chebychev 3 MLS 4 ILU 5 KLU */
  int                     mlsmotype_coarse;   /* 0 SGS 1 Jacobi 2 Chebychev 3 MLS 4 ILU 5 KLU 6 Superlu*/
  double                  mldamp_fine;        /* damping factor fine grid */
  double                  mldamp_med;         /* damping factor fine grid */
  double                  mldamp_coarse;      /* damping factor fine grid */
  double                  mldamp_prolong;     /* damping factor for prolongator smoother */
  double                  ml_threshold;       /* threshold for aggregation/prolongator smoother */
#endif
} AZVAR;


/*----------------------------------------------------------------------*
 | variables needed for HYPRE solver package             m.gee 10/01    |
 *----------------------------------------------------------------------*/
typedef struct _HYPREVARS
{
  enum _HYPREPRECTYP      hypre_prectyp;    /* type of hypre preconditioner */
  INT                     io;               /* flag to set solver quiet */
  INT                     maxiter;          /* max iterations allowed */
  INT                     numiter;          /* number of iterations taken */
  DOUBLE                  resnorm;          /* residual norm achieved */
  INT                     reuse;            /* reuse feature (not yet impl.) */
  DOUBLE                  tol;              /* user-given tolerance */
  INT                     kryldim;          /* dimension of krylov subspace */
  DOUBLE                  threshold;        /* parameters for amg, see manual */
  INT                     sweep[4];
  INT                     ifill;            /* fill in level for ilu */
  DOUBLE                  dfill;            /* fill in level in percent of the original
                                               matrix for ilu and parasails */
  INT                     bj;               /* ? */

  INT                     parasymm;         /* parasails preconditioner parameters */
  INT                     paralevel;
  DOUBLE                  parathresh;
  DOUBLE                  parafilter;
} HYPREVARS;


/*----------------------------------------------------------------------*
 | variables needed for solver ParSuperLU                m.gee 10/01    |
 *----------------------------------------------------------------------*/
typedef struct _PSUPERLUVARS
{
  INT                     reuse;            /* in progress.... */
} PSUPERLUVARS;


/*----------------------------------------------------------------------*
 | variables needed for solver Lapack                    m.gee 11/01    |
 *----------------------------------------------------------------------*/
typedef struct _LAPACKVARS
{
  INT                     reuse;            /* in progress.... */
} LAPACKVARS;


/*----------------------------------------------------------------------*
 | a skyline matrix                                      m.gee 01/02    |
 | this structure holds a skyline matrix to be solved with colsol       |
 *----------------------------------------------------------------------*/
typedef struct _SKYMATRIX
{
  INT                     is_init;         /* was this matrix initialized ? */
  INT                     is_factored;     /* is this matrix already factored ? */
  INT                     ncall;           /* how often was this matrix solved */

  INT                     numeq_total;     /* total number of unknowns */
  INT                     numeq;           /* number of unknowns updated on this proc */
  INT                     nnz_total;       /* total number of nonzero entries */
  INT                     nnz;             /* number of nonzeros on this proc */

  struct _ARRAY           update;          /* sorted list of dofs updated on this proc */
  struct _ARRAY           maxa;
  struct _ARRAY           maxaf;
  struct _ARRAY           A;
  DOUBLE                  det;
  /* some arrays that are used for parallel assembly, mainly in the case of
   * inter-proc-coupling conditions */
#ifdef PARALLEL
  INT                     numcoupsend;     /* number of coupling information
                                              to be send by this proc */
  INT                     numcouprecv;     /* number of coupling information
                                              to be recv. by this proc */
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
  INT                     is_init;         /* was this matrix initialized ? */
  INT                     is_factored;     /* is this matrix already factored ? */
  INT                     ncall;           /* how often was this matrix solved */

  INT                     numeq_total;     /* total number of unknowns */
  INT                     numeq;           /* number of unknowns updated on this proc */
  INT                     nnz_total;       /* total number of nonzero entries */
  INT                     nnz;             /* number of nonzeros on this proc */

  INT                     icntl[20];
  INT                     comm;

  struct _ARRAY           update;          /* sorted list of dofs updated on this proc */
  struct _ARRAY           irn_loc;         /* proc-local row pointer vector */
  struct _ARRAY           irn_locf;        /* fortran style pointer vector of irn_loc */
  struct _ARRAY           jcn_loc;         /* proc-local column pointer vector */
  struct _ARRAY           jcn_locf;        /* fortran style pointer vector of jcn_loc */
  struct _ARRAY           A_loc;           /* values of the matrix */
  struct _ARRAY           rowptr;          /* INT vector holding the begin of each row in irn_loc */

  struct _ARRAY           irn_glob;        /* on imyrank=0 the global row/column pointer arrays */
  struct _ARRAY           jcn_glob;

  /* some arrays that are used for parallel assembly, mainly in the
   * case of inter-proc-coupling conditions */
#ifdef PARALLEL
  INT                     numcoupsend;     /* number of coupling information
                                              to be send by this proc */
  INT                     numcouprecv;     /* number of coupling information
                                              to be recv. by this proc */
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
  INT                     is_init;         /* was this matrix initialized ? */
  INT                     is_factored;     /* is this matrix already factored ? */
  INT                     reuse;           /* is last factorization to be used again? */
  INT                     ncall;           /* how often was this matrix solved */

  INT                     numeq_total;     /* total number of unknowns */
  INT                     numeq;           /* number of unknowns updated on this proc */
  INT                     nnz_total;       /* total number of nonzero entries */
  INT                     nnz;             /* number of nonzeros on this proc */

  INT                     comm;

  struct _ARRAY           update;          /* sorted list of dofs updated on this proc */
  struct _ARRAY           Ap;              /* column pointer vector */
  struct _ARRAY           Ai;              /* row pointer vector */
  struct _ARRAY           Ax;              /* values of the matrix */

  /* some arrays that are used for parallel assembly, mainly in the
   * case of inter-proc-coupling conditions */
#ifdef PARALLEL
  INT                     numcoupsend;     /* number of coupling information
                                              to be send by this proc */
  INT                     numcouprecv;     /* number of coupling information
                                              to be recv. by this proc */
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
  INT                     is_init;         /* was this matrix initialized ? */
  INT                     is_factored;     /* is this matrix already factored ? */
  INT                     ncall;           /* how often was this matrix solved */

  INT                     numeq_total;     /* total number of unknowns */
  INT                     numeq;           /* number of unknowns updated on this proc */
  INT                     nnz_total;       /* total number of nonzero entries */
  INT                     nnz;             /* number of nonzeros on this proc */

  struct _ARRAY           update;          /* sorted list of dofs updated on this proc */
  struct _ARRAY           A;               /* the dense matrix */
  struct _ARRAY           ipiv;            /* pivoting information */
  INT                     lwork;           /* work space for sym. Lapack solver */
  struct _ARRAY           work;            /* work space for sym. Lapack solver */


  /* some arrays that are used for parallel assembly, mainly in the
   * case of inter-proc-coupling conditions */
#ifdef PARALLEL
  INT                     numcoupsend;     /* number of coupling information
                                              to be send by this proc */
  INT                     numcouprecv;     /* number of coupling information
                                              to be recv. by this proc */
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
  INT                     is_init;          /* was this matrix initialized ? */
  INT                     is_factored;      /* is this matrix already factored ? */
  INT                     ncall;            /* how often was this matrix solved */

  INT                     numeq_total;      /* total number of unknowns */
  INT                     numeq;            /* number of unknowns updated on this proc */
  INT                     nnz_total;        /* total number of nonzero entries */
  INT                     nnz;              /* number of nonzeros on this proc */

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

  /* some arrays that are used for parallel assembly, mainly in the
   * case of inter-proc-coupling conditions */
  INT                     numcoupsend;     /* number of coupling information
                                              to be send by this proc */
  INT                     numcouprecv;     /* number of coupling information
                                              to be recv. by this proc */
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
  INT                     is_init;          /* was this matrix initialized ? */
  INT                     is_factored;      /* does precond. information exist ? */
  INT                     ncall;            /* how often was this matrix solved */
  INT                     numeq_total;      /* total number of unknowns */
  INT                     numeq;            /* number of unknowns updated on this proc */
  INT                     nnz;              /* number of nonzeros on this proc */

  struct _ARRAY           perm;                /* permutation of update for each proc, type is IA*/
  struct _ARRAY           perm_sizes;          /* size of perm on each proc, type is IV */
  struct _ARRAY           update;              /* ascending list of dofs on all procs, type is IA*/
  struct _ARRAY           bindx;               /* see documenation for DMSR format (Aztec2.1) */

#ifdef HYPRE_PACKAGE
  HYPRE_IJMatrix          ij_matrix;           /* the matrix itself, hidden away by Hypre*/
  HYPRE_Solver            solver;              /* hypre solver structure */
  HYPRE_Solver            precond;             /* hypre preconditioner structure */
#endif

  /* some arrays that are used for parallel assembly, mainly in the
   * case of inter-proc-coupling conditions */
#ifdef PARALLEL
  INT                     numcoupsend;     /* number of coupling information
                                              to be send by this proc */
  INT                     numcouprecv;     /* number of coupling information
                                              to be recv. by this proc */
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
  INT                     is_init;          /* was this matrix initialized ? */
  INT                     is_factored;      /* does precond. information exist ? */
  INT                     is_transformed;   /* has the transformation been done? */
  INT                     ncall;            /* how often was this matrix solved */
  INT                     numeq_total;      /* total number of unknowns */
  INT                     numeq;            /* number of unknowns updated on this proc */
  INT                     nnz;              /* number of nonzeros on this proc */

  struct _ARRAY           update;           /* list of dofs updated on this proc */
  INT                     shift;            /* binary shift for searching in update */
  INT                    *bins;             /* binary mirror of update */
  struct _ARRAY           bindx;            /* the sparse matrix */
  struct _ARRAY           bindx_backup;     /* backup of bindx, as bindx is altered by solver */
  struct _ARRAY           val;              /* values of matrix */
  struct _ARRAY           val_backup;       /* backup of values of matrix as val is altered by solver */

  INT                    *invupdate;
  INT                    *invbindx;

#ifdef AZTEC_PACKAGE
  DOUBLE                  params[AZ_PARAMS_SIZE];    /* Aztec parameters */
  DOUBLE                  status[AZ_STATUS_SIZE];    /* Aztec return status */
  INT                     proc_config[AZ_PROC_SIZE]; /* MPI-configuration */
  INT                     options[AZ_OPTIONS_SIZE];  /* Aztec options */
  INT                    *data_org;                  /* Aztec internal data org. */
  INT                    *external;                  /* list of external dofs often needed
                                                        by this proc */
  INT                    *update_index;              /* list of dofs updated on this proc */
  INT                    *extern_index;              /* list of external related dofs */
  INT                     N_external;                /* number of external dofs often needed
                                                        by this proc */

  AZ_MATRIX              *Amat;                      /* the matrix object */
  AZ_PRECOND             *Aprec;                     /* the preconditioner object */
#endif

  /* some arrays that are used for parallel assembly, mainly in the
   * case of inter-proc-coupling conditions */
#ifdef PARALLEL
  INT                     numcoupsend;     /* number of coupling information
                                              to be send by this proc */
  INT                     numcouprecv;     /* number of coupling information
                                              to be recv. by this proc */
  struct _ARRAY          *couple_d_send;   /* send and receive buffers if necessary */
  struct _ARRAY          *couple_i_send;
  struct _ARRAY          *couple_d_recv;
  struct _ARRAY          *couple_i_recv;
#endif
} AZ_ARRAY_MSR;



/*----------------------------------------------------------------------*
 | a Spooles Matrix                                       m.gee 4/03    |
 *----------------------------------------------------------------------*/
typedef struct _SPOOLMAT
{
  INT                     is_init;         /* was this matrix initialized ? */
  INT                     is_factored;     /* does precond. information exist ? */
  INT                     ncall;           /* how often was this matrix solved */
  INT                     numeq_total;     /* total number of unknowns */
  INT                     numeq;           /* number of unknowns updated on this proc */
  INT                     nnz;             /* number of nonzeros on this proc */

  struct _ARRAY           update;          /* list of dofs updated on this proc */
  struct _ARRAY           irn_loc;         /* proc-local row pointer vector */
  struct _ARRAY           jcn_loc;         /* proc-local column pointer vector */
  struct _ARRAY           rowptr;          /* INT vector holding the begin of each
                                              row in irn_loc */
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

  /* some arrays that are used for parallel assembly, mainly in the
   * case of inter-proc-coupling conditions */
#ifdef PARALLEL
  INT                     numcoupsend;     /* number of coupling information
                                              to be send by this proc */
  INT                     numcouprecv;     /* number of coupling information
                                              to be recv. by this proc */
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
  INT                is_init;

  /* input */
  INT                numeq;       /* number of equations */
  INT                nnz;         /* number of nonzeroes */
  INT                output;      /* fortran unit number =6 -> screen */
  DOUBLE             cond;        /* condition number */
  struct _ARRAY      colstr;      /* gives the index in rowind of the
                                     first nonzero in the lower triangular
                                     part of column j of the matrix */
  struct _ARRAY      rowind;      /* list of row indices for all nonzeros */
  /* output */
  DOUBLE             rcond;       /* estimate the reciprocal of the l-norm
                                     condition number */
  INT                inrtia[3];   /* number of positive, negative and
                                     an indicator if there are zero
                                     eigenvalues */
  DOUBLE             global[150]; /* global communication array */
  INT                ierr;        /* = 0; normal return */
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
  INT                     numeq_total;     /* total size of distr. vector */
  INT                     numeq;           /* local size of distr. vector */
  struct _ARRAY           vec;             /* local piece of distr. vector */
} DIST_VECTOR;


/*----------------------------------------------------------------------*
  | define the global structure pointer *solv                            |
  |                                                                      |
  | global variable *solv, vector of lenght numfld of structures SOLVAR  |
  | is defined in solver_control.c
  |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/


/*!---------------------------------------------------------------------
  \brief a OLL Matrix

  <pre>                                                              mn 02/03
  This structure contains a matrix in OLL Format (orthopgonal linked list)
  </pre>
  -------------------------------------------------------------------------*/
typedef struct _OLL
{
  INT                     is_init;       /*!< was this matrix initialized ? */
  INT                     is_masked;     /*!< was this matrix masked ? */
  INT                     is_copied;     /*!< was this matrix copied ? */

  INT                     sparsepat;     /*!< sparsepattern was printed to gnu */

  INT                     numeq_total;   /*!< total number of unknowns */
  INT                     numeq;         /*!< number of unknowns updated on this proc */
  INT                     nnz;           /*!< number of nonzeros on this proc */

  struct _ARRAY           update;        /*!< list of dofs updated on this proc */

  INT                     rdim;          /*!< the local row dimension */
  INT                     cdim;          /*!< the local col dimension */
  struct _MATENTRY      **row;           /*!< pointer vector to the rows */
  struct _MATENTRY      **col;           /*!< pointer vector to the rows */

  INT                     total;         /*!< total size of spare */
  INT                     used;          /*!< used entries of spare */
  struct _MATENTRY       *spare;         /*!< vectore of spare MATENTRIES */

  union _SPARSE_ARRAY    *sysarray;      /*!< sparse array in solver format */
  enum  _SPARSE_TYP      *sysarray_typ;  /*!< typ for sparse array */

  struct _ARRAY           gdofrecv;      /*!< ghost dof list to receive
                                           gdofrecv.a.ia[0..fdim-1]    = dof
                                           numbers of external dofs needed
                                           sorted in ascending order */
  struct _ARRAY           recvbuff;
  struct _ARRAY           computebuff;   /*!< after receiving all messages
                                           they are sorted to computebuff */
#ifdef PARALLEL
  MPI_Status             *status;        /*!< receive status */
#endif
  struct _ARRAY           gdofsend;      /*!< ghost dof list to send to other
                                           procs (3D integer) gdofsend.a.ia[proc][0] =
                                           number of values to send to proc proc
                                           gdofsend.a.ia[proc][ 1..gdofsend.a.i3[proc][0] ] =
                                           dof numbers to be send to proc proc in ascending order */
  struct _ARRAY           sendbuff;
#ifdef PARALLEL
  MPI_Request            *request;       /*!< send request */
#endif

  /* some arrays that are used for parallel assembly, mainly in the
   * case of inter-proc-coupling conditions */
#ifdef PARALLEL
  INT                     numcoupsend;   /*!< number of coupling information
                                           to be send by this proc */
  INT                     numcouprecv;   /*!< number of coupling information
                                           to be recv. by this proc */
  struct _ARRAY          *couple_d_send; /*!< send and receive buffers if necessary */
  struct _ARRAY          *couple_i_send;
  struct _ARRAY          *couple_d_recv;
  struct _ARRAY          *couple_i_recv;
#endif
} OLL;


/*!----------------------------------------------------------------------
  \brief one entry to the oll-matrix

  <pre>                                                              mn 02/03
  This structure contains one non-zero entry of an OLL-matrix
  </pre>
  -------------------------------------------------------------------------*/
typedef struct _MATENTRY
{
  INT          r;             /*!< the true row index */
  INT          c;             /*!< the true col index */
  DOUBLE       val;           /*!< value at place r/c */
  struct _MATENTRY *rnext;    /*!< pointer to next value in row */
  struct _MATENTRY *cnext;    /*!< pointer to next value in column */
} MATENTRY;


/*!---------------------------------------------------------------------
  \brief a wrapper for Trilinos' Epetra_CrsMatrix

  <pre>                                                              gee 09/06
  This structure contains a Trilinos' Epetra_CrsMatrix
  </pre>
  -------------------------------------------------------------------------*/
typedef struct _TRILINOSMATRIX
{
  INT                     is_init;       /*!< was this matrix initialized ? */
  INT                     is_factored;   /*!< was this matrix factored ? */
  INT                     ncall;         /*!< no. of calls to solver */
  INT                     numeq_total;   /*!< total number of unknowns */
  INT                     numeq;         /*!< number of unknowns updated on this proc */

  struct _ARRAY           update;        /*!< list of dofs updated on this proc */

  void*                   epetracomm;    /*!< either Epetra_SerialComm or Epetra_MpiComm */
  void*                   rowmap;        /*!< rowmap of the matrix actually ptr to Epetra_Map */
  void*                   matrix;        /*!< Epetra_CrsMatrix */

  INT                     NumEntriesPerRow; /*!< trilinos Graph build info */
  INT                     StaticProfile; /*!< trilinos Graph build info */

  void*                   linearproblem; /*!< ptr to the Epetra_LinearProblem */
  void*                   solver;        /*!< ptr to any Trilinos solver */
  void*                   params;        /*!< ptr to Teuchos parameter list */

  void*                   prec;          /*!< ptr to any Trilinos preconditioner (type Epetra_Operator) */
  void*                   precmatrix;    /*!< copy of matrix, needed to reuse preconditioner */
  double*                 nullspace;     /*!< ptr to nullspace of this operator (used by ML) */

  union _SPARSE_ARRAY    *sysarray;      /*!< sparse array in solver format (used ofr spooles as subsolver only)*/
  enum  _SPARSE_TYP      *sysarray_typ;  /*!< typ for sparse array */

  /* some arrays that are used for parallel assembly, mainly in the
     case of inter-proc-coupling conditions */
#ifdef PARALLEL
  INT                     numcoupsend;   /*!< number of coupling information
                                           to be send by this proc */
  INT                     numcouprecv;   /*!< number of coupling information
                                           to be recv. by this proc */
  struct _ARRAY          *couple_d_send; /*!< send and receive buffers if necessary */
  struct _ARRAY          *couple_i_send;
  struct _ARRAY          *couple_d_recv;
  struct _ARRAY          *couple_i_recv;
#endif
} TRILINOSMATRIX;

/*!------------------------------------------------------------------------
  \brief matrix needed by the MLPCG solver on the finest grid only

  m.gee 6/01

  matrix needed by the MLPCG solver on the finest grid only
  -------------------------------------------------------------------------*/
typedef struct _DBCSR
{
  INT                     is_init;       /*!< was this matrix initialized ? */
  INT                     is_factored;   /*!< is this matrix already factored ? */
  INT                     ncall;         /*!< how often was this matrix solved */

  INT                     numeq_total;   /*!< total number of unknowns */
  INT                     numeq;         /*!< number of unknowns updated on this proc */
  INT                     nnz;           /*!< number of nonzeros on this proc */

  INT                     owner[MAXPROC][2]; /*!< contains for each proc the lowest and highest
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
                                           gdofrecv.a.ia[0..fdim-1]    = dof numbers
                                           of external dofs needed
                                           sorted in ascending order */
  struct _ARRAY           recvbuff;
  struct _ARRAY           computebuff;   /*!< after receiving all messages they
                                           are sorted to computebuff */
#ifdef PARALLEL
  MPI_Status             *status;        /*!< receive status */
#endif
  struct _ARRAY           gdofsend;      /*!< ghost dof list to send to other
                                           procs (3D integer) gdofsend.a.ia[proc][0] =
                                           number of values to send to proc proc
                                           gdofsend.a.ia[proc][ 1..gdofsend.a.i3[proc][0] ] =
                                           dof numbers to be send to proc proc in ascending order */
  struct _ARRAY           sendbuff;
#ifdef PARALLEL
  MPI_Request            *request;       /*!< send request */
#endif
  INT                     firstcoupledof;/*!< dof number of the first dof that
                                           has interproc coupling */

  struct _DBCSR          *csc;           /*!< the treansposed matrix in compressed
                                           sparse column format */
#ifdef MLPCG
  struct _DBCSR          *ilu;           /*!< the ilu-decomposed asm - matrix */
  struct _DBCSR          *stupid_asm;           /*!< the asm splitted matrix */
#endif
  ARRAY                  *dense;         /*!< for dense solve on coarsest grid */
  ARRAY                  *ipiv;          /*!< for dense solve on coarsest grid */

} DBCSR;



/*----------------------------------------------------------------------*
  | mlpcg                                                     m.gee9/02  |
 *----------------------------------------------------------------------*/
#ifdef MLPCG
#include "solver_mlpcg.h"
#endif


/*----------------------------------------------------------------------*
  | put the prototypes of all                                            |
  | routines that have contact to solver.h in here                       |
  | the global variable solv and the head solver.h is not visible        |
  | from everywhere in the code, so the prototypes of all functions using|
  | structurs declared in solver.h are in here and not in prototypes.h   |
  |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
#include "solver_prototypes.h"



#endif /* #ifndef SOLVER_H */
