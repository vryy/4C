/*----------------------------------------------------------------------*
 | NOTE:                                                                |
 | - if changes or additions are made to this file, a complete recompile|
 |   of the whole code is recommended                                   |
 | - if segmentation violation errors occure in runtime, check the      |
 |   values below, as some of them are the default sizes of arrays      | 
 | - please do not define thousands of all kinds of variables, because  |
 |   they are global, do only define globally important ones            |
 | - do NOT use common words (e.g. JACOBI, NODE, ELEMENT ...)           |
 | - always use strict upper case letters                               |
 |                                                                      |
 *----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
 | some definitions that are important for dynamic memory management    |
 *----------------------------------------------------------------------*/
#ifdef SIXTYFOUR /*--------------- a 64 bit pointer is of size long int */
typedef long int PTRSIZE;
#else/*-------------------------------- a 32 bit pointer is of size int */
typedef int PTRSIZE;
#endif

/*----------------------------------------------------------------------*
 | basic data types, do not use int, double or char !                   |
 *----------------------------------------------------------------------*/
#ifdef INT
#undef INT
#endif
typedef int       INT;
#ifdef DOUBLE
#undef DOUBLE
#endif
typedef double    DOUBLE;
#ifdef CHAR
#undef CHAR
#endif
typedef char      CHAR;

/*----------------------------------------------------------------------*
 | special definitions for special compilers.....                       |
 *----------------------------------------------------------------------*/
/* append underslashs, if necessary */
#undef CCA_APPEND_U
/* append underslash for gnu's linux compiler gcc and g77 */
#ifdef SUSE73 
#define CCA_APPEND_U (1)
#endif
/* append underslash for CUSS Sunfire */
#ifdef SUN 
#define CCA_APPEND_U (1)
#endif
/* append underslash for SIXTYFOUR flag */
#ifdef SIXTYFOUR 
#define CCA_APPEND_U (1)
#endif


#ifdef CCA_APPEND_U
#define dsytrf              dsytrf_
#define dsytri              dsytri_
#define dsytrs              dsytrs_
#define dgetrf              dgetrf_
#define dgetri              dgetri_
#define dgetrs              dgetrs_
#define dsygv               dsygv_
#define dsyevd              dsyevd_
#define dsyev               dsyev_
#define mydsyevx            mydsyevx_
#define colsol              colsol_
#define c1inv6              c1inv6_
#define c1ab                c1ab_
#define mxmatb              mxmatb_
#define mxmabt              mxmabt_
#define fortranpow          fortranpow_
#define mxmab               mxmab_
#define c1invf              c1invf_
#define solveq              solveq_
#define c1jacb              c1jacb_
#define mumps_interface     mumps_interface_
#define iluk                iluk_
#define lusol               lusol_
#define mlpcgveczero        mlpcgveczero_
#define mlpcgvecvec         mlpcgvecvec_
#define mlpcgupdupdvec      mlpcgupdupdvec_
#define mlpcgupdvec         mlpcgupdvec_
#define dveczero            dveczero_
#define iveczero            iveczero_
#define fsdoc               fsdoc_
#endif
#ifndef AZTEC_PACKAGE
void dsytrf(char *uplo, int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);
void dsytri(char *uplo, int *n, double *a, int *lda, int *ipiv, double *work, int *info);
void dsytrs(char *uplo, int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);
void dgetrf(int *m, int *n, double *a, int *lda, int *ipiv, int *info);
void dgetri(int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);
void dgetrs(char *trans, int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);
void dsygv(int *itype, char *jobz, char *uplo, int *n, double *a, int *lda, double *b, int *ldb, double *w, double *work, int *lwork, int *info);
void dsyevd(char *jobz, char *uplo, int *n, double *a, int *lda, double *w, double *work, int *lwork, int *iwork, int *liwork, int *info);
void dsyev(char *jobz, char *uplo, int *n, double *a, int *lda, double *w, double *work, int *lwork, int *info);
#endif
void colsol(double *a, double *v, int *maxa, int *nn, int *nrr, int *nrc, int *nwa, int *nqm, int *nr1, int *nr2, int *kkk, double *det, int *isc, int *nsch, int *ipr, int *info);
void iluk(int *n, double *a, int *ja, int *ia, int *lfil, double *alu, int *jlu, int *ju, int *levs, int *iwk, double *w, int *jw, int *ierr);
void lusol(int *n, double *y, double *x, double *alu, int *jlu, int *ju);
void mlpcgveczero(double *x, int *n);
void mlpcgvecvec(double *x, double *y, double *sum, int *n);
void mlpcgupdupdvec(double *a, double *y, double *facy, double *x, double *facx, int *init, int *n);
void mlpcgupdvec(double *y, double *x, double *fac, int *init, int *n);
void dveczero(double *x, int *n);
void iveczero(int *x, int *n);
void mydsyevx(char *jobz,char *range,char *uplo,int *n,double *a,int *lda,double *vl,double *vu,int *il,int *iu,double *abstol,int *m,double *w,double *z,int *ldz,double *work,int *lwork,int *iwork,int *ifail,int *info);
/*----------------------------------------------------------------------*
 | absolut value of an integer                                          |
 *----------------------------------------------------------------------*/
#define ABS(x)    ((x) <  0  ? (-x) : (x))

/*----------------------------------------------------------------------*
 | absolut value of a double                                            |
 *----------------------------------------------------------------------*/
#define FABS(x)   ((x) < 0.0 ? (-(x)) : (x))

/*----------------------------------------------------------------------*
 | square of a double                                                   |
 *----------------------------------------------------------------------*/
static DOUBLE dsqrarg;
#define DSQR(a) (a*a)

/*----------------------------------------------------------------------*
 | the larger of two doubles (fast version)                             |
 *----------------------------------------------------------------------*/
static DOUBLE dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b), (dmaxarg1) > (dmaxarg2) ? (dmaxarg1) : (dmaxarg2))

/*----------------------------------------------------------------------*
 | the smaller of two doubles (fast version)                            |
 *----------------------------------------------------------------------*/
static DOUBLE dminarg1,dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b), (dminarg1) < (dminarg2) ? (dminarg1) : (dminarg2))

/*----------------------------------------------------------------------*
 | the larger of two integer (fast version)                             |
 *----------------------------------------------------------------------*/
static INT imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b), (imaxarg1) > (imaxarg2) ? (imaxarg1) : (imaxarg2))

/*----------------------------------------------------------------------*
 | the smaller of two integer (fast version)                            |
 *----------------------------------------------------------------------*/
static INT iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b), (iminarg1) < (iminarg2) ? (iminarg1) : (iminarg2))

/*----------------------------------------------------------------------*
 | max number of processors                                             |
 *----------------------------------------------------------------------*/
#define MAXPROC          (16)

/*----------------------------------------------------------------------*
 | size of buffer to attach to intra-communicator in byte               |
 *----------------------------------------------------------------------*/
#define MPIBUFFSIZE      (52428800) /* this is 50 MB */

/*----------------------------------------------------------------------*
 | exact one RAD                                                        |
 *----------------------------------------------------------------------*/
#define RAD              (atan(1.0)/45.0)

/*----------------------------------------------------------------------*
 | exact PI                                                             |
 *----------------------------------------------------------------------*/
#define PI               (asin(1.0)*2.0)

/*----------------------------------------------------------------------*
 | maximum number columns in input file                                 |
 *----------------------------------------------------------------------*/
#define MAXNUMCOL        (500)

/*----------------------------------------------------------------------*
 | maximum number of nodes to an element                                |
 *----------------------------------------------------------------------*/
#define MAXNOD           (27)

/*----------------------------------------------------------------------*
 | maximum number of elements to a node                                 |
 *----------------------------------------------------------------------*/
#define MAXELE           (25)

/*----------------------------------------------------------------------*
 | maximum number of conditions (dirich/neumann) to a node              |
 *----------------------------------------------------------------------*/
#define MAXCONDPERNODE    (6)

/*----------------------------------------------------------------------*
 | maximum number of dofs to a node -> can be more than 6 for shell9    |
 *----------------------------------------------------------------------*/
#define MAXDOFPERNODE    (33)

/*----------------------------------------------------------------------*
 | maximum number of gaussian points in an element                      |
 *----------------------------------------------------------------------*/
#define MAXGAUSS         (360)  /*20 Layers 3x3x2 integrated*/

/*----------------------------------------------------------------------*
 | maximum number of dofs to an element                                 |
 *----------------------------------------------------------------------*/
#define MAXDOFPERELE     (MAXNOD*MAXDOFPERNODE)

/*----------------------------------------------------------------------*
 | maximum number nonzero entries in a row of a sparse system matrix    |
 | is number of nodes to an element *                                   |
 | number of dofs to a node *                                           |
 | number of elements to a node (8) *                                   |
 | 2 (unsymmetric case)                                                 |
 *----------------------------------------------------------------------*/
#define MAX_NNZPERROW     (MAXNOD*MAXDOFPERNODE*8*2)

/*----------------------------------------------------------------------*
 | maximum numberof fields                                              |
 *----------------------------------------------------------------------*/
#define MAXFIELD         (3)

/*----------------------------------------------------------------------*
 | maximum numberof fields                                              |
 *----------------------------------------------------------------------*/
#define MAXTIMECURVE         (5)
/*----------------------------------------------------------------------*
 | numbers                                                              |
 *----------------------------------------------------------------------*/
#define ZERO              (0.0)
#define ONE               (1.0)
#define TWO               (2.0)
#define THREE             (3.0)
#define FOUR              (4.0)
#define FIVE              (5.0)
#define SIX               (6.0)
#define SEVEN             (7.0)
#define EIGHT             (8.0)
#define NINE              (9.0)
#define TEN              (10.0)
#define ELEVEN           (11.0)
#define TWELVE           (12.0)
/*----------------------------------------------------------------------*
 | a set of different tolerances                                        |
 *----------------------------------------------------------------------*/
#define EPS5             (1.0E-05)
#define EPS6             (1.0E-06)
#define EPS7             (1.0E-07)
#define EPS8             (1.0E-08)
#define EPS9             (1.0E-09)
#define EPS10            (1.0E-10)
#define EPS11            (1.0E-11)
#define EPS12            (1.0E-12)
#define EPS13            (1.0E-13)
#define EPS14            (1.0E-14)
#define EPS15            (1.0E-15)
/*----------------------------------------------------------------------*
 | a set of numbers                                                     |
 *----------------------------------------------------------------------*/
#define VERYLARGEINT     (1000000000)

/*----------------------------------------------------------------------*
 | shell8                                                               |
 *----------------------------------------------------------------------*/
#define MAXNOD_SHELL8    (9)
#define NUMDOF_SHELL8    (6)
#define MAXHYB_SHELL8    (45)

/*----------------------------------------------------------------------*
 | shell9                                                               |
 |                                                                      |
 |    !!! set MAXDOFPERNODE according to number of kinematic Layers !!! |
 |                                                                      |
 *----------------------------------------------------------------------*/
#define MAXNOD_SHELL9    (9)
#define MAXKLAY_SHELL9   (10)                   /* max. nr of kin lay */
#define NUMDOF_SHELL9    (3+3*MAXKLAY_SHELL9)   /* numdf = 3 + 3*klay */
#define MAXHYB_SHELL9    (45)
#define A3FAC_SHELL9     (1.0) /* makes it possible to change the norm of a3L */
                               /* A3FAC_SHELL9 = 0.5 : |a3L| = 0.5 * hL -> as in shell8 */
                               /* A3FAC_SHELL9 = 1.0 : |a3L| = 1.0 * hL -> like in Dis. Braun */
/*----------------------------------------------------------------------*
 | wall1                                                                |
 *----------------------------------------------------------------------*/
#define MAXNOD_WALL1     (9)
/*----------------------------------------------------------------------*
 | brick8                                                               |
 *----------------------------------------------------------------------*/
#define MAXNOD_BRICK1    (20)
#define NUMDOF_BRICK1    (3)
/*----------------------------------------------------------------------*
 | fluid                                                              |
 *----------------------------------------------------------------------*/
#define MAXQINTC (6)   /* max. number of implemented integration cases for QUADS */
#define MAXQINTP (6)   /* max. number of integration parameters  for QUADS */ 
#define MAXTINTC (11)  /* max. number of implemented integration cases for TRIS */
#define MAXTINTP (13)  /* max. number of integration parameters  for TRIS */
/*----------------------------------------------------------------------*
 | fluid2                                                               |
 *----------------------------------------------------------------------*/
#define NUM_F2_VELDOF (2) /* number of velocity dofs per node */
#define MAXNOD_F2 (9)    /* max. number of nodes per element */
/*----------------------------------------------------------------------*
 | fluid3                                                               |
 *----------------------------------------------------------------------*/
#define NUM_F3_VELDOF (3) /* number of velocity dofs per node */
#define MAXNOD_F3 (27)    /* max. number of nodes per element */
