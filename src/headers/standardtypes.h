#include "math.h"
#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "ctype.h"
#ifdef PARALLEL 
   #include "mpi.h"
#endif
#include "definitions.h"
#include "am.h"
#include "enums.h"
#include "partition.h"
#include "geometrytypes.h"
#include "design.h"
#include "load.h"
#include "conditions.h"
#include "materials.h"
#include "bunker.h"
/*----------------------------------------------------------------------*
 | type definitions of basic structures                   m.gee 8/00    |
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | file I/O variables & fr-system                         m.gee 8/00    |
 *----------------------------------------------------------------------*/
typedef struct _FILES
{
/*------------------------------------------------------------- file I/O */
char             *inputfile_name;         /* input file name             */

char             *outputfile_kenner;      /* output file kenner          */
char              outputfile_name[100];   /* output file name            */
size_t            outlenght;              /* lenght of output file kenner*/
int               num_outputfiles;        /* number of output files      */
int               pss_counter;

FILE             *in_input;               /* file-pointer input file     */
FILE             *out_out;                /* file-pointer .out  file     */
FILE             *out_err;                /* file-pointer .err  file     */
FILE             *out_pss;                /* file-pointer .pss  file     */

FILE             *gidmsh;                 /* file pointer .flavia.msh    */
FILE             *gidres;                 /* file pointer .flavia.res    */
/*---------------------- variables needed by the free-field-input system */
char              title[5][500];          /* problem title                */
char              line[500];              /* linebuffer for reading       */
char            **input_file;             /* copy of the input file       */
char             *input_file_hook;
int               numcol;
int               numrows;                /* number of rows in inputfile  */
int               actrow;                 /* rowpointer used by fr        */
char             *actplace;               /* pointer to actual place in input-file */
} FILES;

/*----------------------------------------------------------------------*
 | tracing of time & array bugs                           m.gee 8/00    |
 *----------------------------------------------------------------------*/
typedef struct _TRACEROUT
{
struct _TRACEROUT   *prev;
struct _TRACEROUT   *next;
char                 name[50];
enum 
   {
    dsnone,
    dsin,
    dsout
   }                 dsroutcontrol;
   
} TRACEROUT;

typedef struct _TRACE
{
int                 trace_on;             /* switches trace on/off */
int                 num_arrays;           /* number of current int-arrays */
int                 size_arrays;          /* size of vector arrays    */
struct _ARRAY     **arrays;               /* pointer to the arrays    */
int                 deepness;
struct _TRACEROUT   routine[100];         /* chained list ring to trace routine */
struct _TRACEROUT   *actroutine;          /* ptr to actual routine */
} TRACE;

/*----------------------------------------------------------------------*
 | FIELD                                                   m.gee 8/00    |
 *----------------------------------------------------------------------*/
typedef struct _FIELD
{
enum _FIELDTYP    fieldtyp;       /* typ of field */

int               numnp;          /* number of nodes in this field */                          
int               numele;         /* number of elements in this field */
int               numeq;          /* number of unknowns */
int               numdf;          /* number of degrees of freedom */

int               db_bunker_pss_handle; /* pss-file handle for this bunker */
int               db_bunker_Id; /* Id of bunker, created for this field */

struct _ELEMENT   *element;       /* vector of elements */
struct _NODE      *node;          /* vector of nodes */

} FIELD;

/*----------------------------------------------------------------------*
 | general problem-variables                              m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _GENPROB
{
int               nele;          /* total number of elements */
int               nnode;         /* total number of nodes */
int               ndim;          /* dimension of problem (2 or 3) */
int               nmat;          /* total number of materials */
int               numfld;        /* number of fields */
int               numdf;         /* maximum number of dofs to one node (not used)*/
int               design;        /* need design information or not */
int               restart;       /* is restart or not (not used) */

enum _PROBLEM_TYP probtyp;       /* type of problem, see enum.h */
enum _TIME_TYP    timetyp;       /* type of time, see enum.h */

} GENPROB;

/*----------------------------------------------------------------------*
 | general IO-flags                                      m.gee 12/01    |
 *----------------------------------------------------------------------*/
typedef struct _IO_FLAGS
{
int               struct_disp_file;
int               struct_stress_file;
int               struct_disp_gid;
int               struct_stress_gid;
} IO_FLAGS;


/*----------------------------------------------------------------------*
 | general dynamic-variables                              m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _DYNAMIC                  /* not used */
{
char dyntyp[50];

int                nstep;
int                damp;
int                iter;
int                maxiter;

int                numcurve;
struct _CURVE     *curve;

double             dt;
double             maxtime;
double             beta;
double             gamma;
double             alpha_m;
double             alpha_f;
double             m_damp;
double             k_damp;
} DYNAMIC;

/*----------------------------------------------------------------------*
 | enum NR_CONTROLTYP                                    m.gee 11/01    |
 *----------------------------------------------------------------------*/
typedef enum _NR_CONTROLTYP         /* type of nonlinear static control */
{
                       control_none,
                       control_disp,
                       control_load,
                       control_arc
} NR_CONTROLTYP;                         
/*----------------------------------------------------------------------*
 | general static-variables                               m.gee 6/01    |
 *----------------------------------------------------------------------*/
typedef struct _STATIC_VAR               
{
int                 geolinear;          /* is linear calculation */
int                 geononlinear;       /* is nonlinear calculation (redundant!) */
enum _NR_CONTROLTYP nr_controltyp;      /* type of control */
int                 nstep;              /* number of steps */
int                 maxiter;            /* max number of iterations */
double              tolresid;           /* tolerance of residual forces */
double              toldisp;            /* tolerance of residual displacements */
double              stepsize;           /* ... */

struct _NODE       *controlnode;        /* ptr to control node */
int                 control_node_global;/* global control node Id (redundant) */
int                 control_dof;        /* dof to be controlled */
} STATIC_VAR;

/*----------------------------------------------------------------------*
 | general static-control-variables                       m.gee 6/01    |
 *----------------------------------------------------------------------*/
typedef struct _STANLN   /* under construction ....*/
{
double              sp1;
double              csp;                  /* current stiffness parameter */
double              rlold;
double              rlnew;
double              rlpre;

double              renorm;
double              rinorm;
double              rrnorm;

double              renergy;

struct _ARRAY       arcfac;               /* load factors of increments */
} STANLN;

/*----------------------------------------------------------------------*
 | Prototypes                                            m.gee 06/01    |
 *----------------------------------------------------------------------*/
#include "prototypes.h" /* prototypes to all routines */
/*----------------------------------------------------------------------*
 | define the global structures & structure pointers                    |
 |                                                                      |
 |                                                                      |
 |                                                        m.gee 8/00    |
 *----------------------------------------------------------------------*/
FILES             allfiles;
PAR               par;
PARTITION        *partition;
GENPROB           genprob;
#ifdef DEBUG
TRACE             trace;
#endif
DESIGN            *design;
FIELD             *field;
DYNAMIC           *dyn;
STATIC_VAR        *statvar;
MATERIAL          *mat;
IO_FLAGS           ioflags;
