/*----------------------------------------------------------------------*
 | includes of Ansi C standard headers                   m.gee 8/00    |
 *----------------------------------------------------------------------*/
#include "math.h"
#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "ctype.h"
/*----------------------------------------------------------------------*
 | includes of mpi of C                                   m.gee 8/00    |
 *----------------------------------------------------------------------*/
#ifdef PARALLEL 
   #include "mpi.h"
#endif
/*----------------------------------------------------------------------*
 | definitions file                                       m.gee 8/00    |
 *----------------------------------------------------------------------*/
#include "definitions.h"
/*----------------------------------------------------------------------*
 | structure types used by the array management           m.gee 8/00    |
 *----------------------------------------------------------------------*/
#include "am.h"
/*----------------------------------------------------------------------*
 | various types of enums                                 m.gee 8/00    |
 *----------------------------------------------------------------------*/
#include "enums.h"
/*----------------------------------------------------------------------*
 | structures concerning domain decomposition                m.gee 8/00 |
 | and intra-communicators                                              |
 *----------------------------------------------------------------------*/
#include "partition.h"
/*----------------------------------------------------------------------*
 | Nodes and Elements on finite element level             m.gee 8/00    |
 *----------------------------------------------------------------------*/
#include "geometrytypes.h"
/*----------------------------------------------------------------------*
 | Points, Lines, Surfaces and Volumes on design level    m.gee 8/00    |
 *----------------------------------------------------------------------*/
#include "design.h"
/*----------------------------------------------------------------------*
 | structures used by load cases                          m.gee 8/00    |
 *----------------------------------------------------------------------*/
#include "load.h"
/*----------------------------------------------------------------------*
 | standard conditions for nodes and elements on finite element level   |
 |                                                        m.gee 8/00    |
 *----------------------------------------------------------------------*/
#include "conditions.h"
/*----------------------------------------------------------------------*
 | material laws                                          m.gee 8/00    |
 *----------------------------------------------------------------------*/
#include "materials.h"
/*----------------------------------------------------------------------*/

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
int               pss_counter;            /* number of records on pss-file */

FILE             *in_input;               /* file-pointer input file     */
FILE             *out_out;                /* file-pointer .out  file     */
FILE             *out_err;                /* file-pointer .err  file     */
FILE             *out_pss;                /* file-pointer .pss  file     */

FILE             *gidmsh;                 /* file pointer .flavia.msh    */
FILE             *gidres;                 /* file pointer .flavia.res    */

/*---------------------- variables needed by the free-field-input system */
char              title[5][500];          /* problem title                */
char              line[500];              /* linebuffer for reading       */
char            **input_file;             /* copy of the input file in rows and columns */
char             *input_file_hook;        /* ptr the copy of the input file is allocated to */
int               numcol;                 /* number of cols in inputfile  */
int               numrows;                /* number of rows in inputfile  */
int               actrow;                 /* rowpointer used by fr        */
char             *actplace;               /* pointer to actual place in input-file */
} FILES;

/*----------------------------------------------------------------------*
 | tracing of time & array bugs                           m.gee 8/00    |
 | This structures is used by the chained list that keeps trrack of     |
 | the function calls                                                   |
 |                                                                      |
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
/*----------------------------------------------------------------------*
 | tracing of time & array bugs                           m.gee 8/00    |
 |                                                                      |
 *----------------------------------------------------------------------*/
typedef struct _TRACE
{
int                 trace_on;             /* switches trace on/off */
int                 num_arrays;           /* number of current int-arrays */
int                 size_arrays;          /* size of vector arrays    */
struct _ARRAY     **arrays;               /* pointer to the arrays    */
int                 deepness;             /* the actual deepness of the calling tree */
struct _TRACEROUT   routine[100];         /* chained list ring to trace routine */
struct _TRACEROUT   *actroutine;          /* ptr to actual routine */
} TRACE;

/*----------------------------------------------------------------------*
 | FIELD                                                   m.gee 8/00   |
 |                                                                      |
 | This is the Home-structure of each physical field on finite element  |
 | level. It holds the number of nodes and elements, and the nodes and  |
 | elements themselves. This structure is redundant on all procs, which |
 | means that each proc holds all nodes and elements, but not all       |
 | contents of all nodes and elements                                   |
 *----------------------------------------------------------------------*/
typedef struct _FIELD
{
enum _FIELDTYP    fieldtyp;       /* typ of field */

int               numnp;          /* number of nodes in this field */                          
int               numele;         /* number of elements in this field */
int               numeq;          /* number of unknowns */
int               numdf;          /* number of degrees of freedom */

struct _ELEMENT   *element;       /* vector of elements */
struct _NODE      *node;          /* vector of nodes */

} FIELD;

/*----------------------------------------------------------------------*
 | general problem-variables                              m.gee 4/01    |
 | General information is held here                                     |           
 *----------------------------------------------------------------------*/
typedef struct _GENPROB
{
int               nele;          /* total number of elements over all fields */
int               nnode;         /* total number of nodes over all fields*/
int               ndim;          /* dimension of problem (2 or 3) */
int               nmat;          /* total number of material laws */
int               numfld;        /* number of fields */
int               numdf;         /* maximum number of dofs to one node (not used, in progress)*/
int               design;        /* need design information or not */
int               restart;       /* is restart or not (not used yet) */

enum _PROBLEM_TYP probtyp;       /* type of problem, see enum.h */
enum _TIME_TYP    timetyp;       /* type of time, see enum.h */

} GENPROB;

/*----------------------------------------------------------------------*
 | general IO-flags                                      m.gee 12/01    |
 | flags to switch the output of certain data on or off, read from gid  |
 *----------------------------------------------------------------------*/
typedef struct _IO_FLAGS
{
int               struct_disp_file;    /* write displacements to .out */
int               struct_stress_file;  /* write structural stress to .out */
int               struct_disp_gid;     /* write structural displacements to .flavia.res */
int               struct_stress_gid;   /* write structural stresses to .flavia.res */
} IO_FLAGS;


/*----------------------------------------------------------------------*
 | general dynamic-variables                              m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _DYNAMIC                  /* not used */
{
char dyntyp[50];

int                nstep;    /* this all is in progress... */
int                damp;     /* some of these values are read from gid */
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
 | type of control algorithm for Newton-Raphson in nonlinear structural |
 | analysis                                                             |
 *----------------------------------------------------------------------*/
typedef enum _NR_CONTROLTYP         /* type of nonlinear static control */
{
                       control_none,
                       control_disp,     /* displacement control */
                       control_load,     /* not impl. yet */
                       control_arc       /* not implem. yet */
} NR_CONTROLTYP;                         

/*----------------------------------------------------------------------*
 | general static-variables                               m.gee 6/01    |
 | variables used by linear or nonlinear structural static analysis     |
 *----------------------------------------------------------------------*/
typedef struct _STATIC_VAR               
{
int                 geolinear;          /* is linear calculation */
int                 geononlinear;       /* is nonlinear calculation */
enum _NR_CONTROLTYP nr_controltyp;      /* type of control */
int                 nstep;              /* number of steps */
int                 maxiter;            /* max number of iterations in NR */
double              tolresid;           /* tolerance of residual forces */
double              toldisp;            /* tolerance of residual displacements */
double              stepsize;           /* steplenght */

struct _NODE       *controlnode;        /* ptr to control node */
int                 control_node_global;/* global control node Id (redundant) */
int                 control_dof;        /* dof of control node to be controlled */
} STATIC_VAR;

/*----------------------------------------------------------------------*
 | general static-control-variables                       m.gee 6/01    |
 | variables to perform Newton Raphson                                  |
 *----------------------------------------------------------------------*/
typedef struct _STANLN  
{
double              sp1;
double              csp;                  /* current stiffness parameter */
double              rlold;                /* load factor of last step */
double              rlnew;                /* load factor of actual step */
double              rlpre;                /* load factor from predictor */

double              renorm;               /* some norms */
double              rinorm;
double              rrnorm;

double              renergy;

struct _ARRAY       arcfac;               /* vector of load factors of increments */
} STANLN;

/*----------------------------------------------------------------------*
 | Prototypes                                            m.gee 06/01    |
 | prototypes of all files using standardtypes.h                        |
 *----------------------------------------------------------------------*/
#include "prototypes.h"
/*----------------------------------------------------------------------*
 | global variables                                                     |
 |                                                                      |
 |                                                                      |
 |                                                        m.gee 8/00    |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | ranks and communicators                                              |
 | This structure struct _PAR par; is defined in main_ccarat.c
 *----------------------------------------------------------------------*/
extern struct _PAR   par;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of partitions, size numfld                                    |
 | struct _PARTITION    *partition; defined in gobal_control.c          |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | struct _GENPROB       genprob; defined in global_control.c           |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure allfiles, which holds all file pointers                    |
 | is defined in input_control_global.c                                 |
 | struct _FILES  allfiles;                                             |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | global variable *solv, vector of lenght numfld of structures SOLVAR  |
 | defined in solver_control.c                                          |
 |                                                                      |
 |                                                       m.gee 11/00    |
 | struct _SOLVAR  *solv;                                               |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | enum _CALC_ACTION                                      m.gee 1/02    |
 | command passed from control routine to the element level             |
 | to tell element routines what to do                                  |
 | defined globally in global_calelm.c                                  |
 | enum _CALC_ACTION calc_action[MAXFIELD];                             |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | tracing variables                                                    |
 | defined in pss_ds.c                                                  |
 | #ifdef DEBUG                                                         |
 | extern struct _TRACE         trace;                                  |
 | #endif                                                               |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate design if needed                                 |
 | defined in global_control.c                                          |
 | struct _DESIGN       *design;                                        |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 | struct _FIELD         *field;                                        |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | struct _DYNAMIC      *dyn;                                           |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate static variables if needed                       |
 | defined in global_control.c                                          |
 | struct _STATIC_VAR   *statvar;                                       |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c                                          |
 | struct _MATERIAL     *mat;                                           |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 | struct _IO_FLAGS        ioflags;                                     |
 *----------------------------------------------------------------------*/
