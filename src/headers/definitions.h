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
static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

/*----------------------------------------------------------------------*
 | the larger of two doubles (fast version)                             |
 *----------------------------------------------------------------------*/
static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b), (dmaxarg1) > (dmaxarg2) ? (dmaxarg1) : (dmaxarg2))

/*----------------------------------------------------------------------*
 | the smaller of two doubles (fast version)                            |
 *----------------------------------------------------------------------*/
static double dminarg1,dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b), (dminarg1) < (dminarg2) ? (dminarg1) : (dminarg2))

/*----------------------------------------------------------------------*
 | the larger of two integer (fast version)                             |
 *----------------------------------------------------------------------*/
static int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b), (imaxarg1) > (imaxarg2) ? (imaxarg1) : (imaxarg2))

/*----------------------------------------------------------------------*
 | the smaller of two integer (fast version)                            |
 *----------------------------------------------------------------------*/
static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b), (iminarg1) < (iminarg2) ? (iminarg1) : (iminarg2))

/*----------------------------------------------------------------------*
 | exact one RAD                                                        |
 *----------------------------------------------------------------------*/
#define RAD              (atan(1.0)/45.0)

/*----------------------------------------------------------------------*
 | exact PI                                                             |
 *----------------------------------------------------------------------*/
#define PI               (asin(1.0)*2.0)

/*----------------------------------------------------------------------*
 | maximum number of nodes to an element                                |
 *----------------------------------------------------------------------*/
#define MAXNOD           (27)

/*----------------------------------------------------------------------*
 | maximum number of elements to a node                                 |
 *----------------------------------------------------------------------*/
#define MAXELE           (25)

/*----------------------------------------------------------------------*
 | maximum number of dofs to a node                                     |
 *----------------------------------------------------------------------*/
#define MAXDOFPERNODE    (6)

/*----------------------------------------------------------------------*
 | maximum number of gaussian points in an element                      |
 *----------------------------------------------------------------------*/
#define MAXGAUSS         (27)

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

/*----------------------------------------------------------------------*
 | shell8                                                               |
 *----------------------------------------------------------------------*/
#define MAXNOD_SHELL8    (9)
#define NUMDOF_SHELL8    (6)
#define MAXHYB_SHELL8    (45)
