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
 | structures concerning domain decomposition             m.gee 8/00    |
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
 | structures used by load curves                         m.gee 8/00    |
 *----------------------------------------------------------------------*/
#include "curve.h"
/*----------------------------------------------------------------------*
 | standard conditions for nodes and elements on finite element level   |
 |                                                        m.gee 8/00    |
 *----------------------------------------------------------------------*/
#include "conditions.h"
/*----------------------------------------------------------------------*
 | material laws                                          m.gee 8/00    |
 *----------------------------------------------------------------------*/
#include "materials.h"
/*----------------------------------------------------------------------*
 | structures used for static analysis                    m.gee 2/02    |
 *----------------------------------------------------------------------*/
#include "static_analysis.h"
/*----------------------------------------------------------------------*
 | structures used for dynamic analysis                   m.gee 2/02    |
 *----------------------------------------------------------------------*/
#include "dynamic.h"
/*----------------------------------------------------------------------*
 | structures used for bug and time tracing               m.gee 2/02    |
 *----------------------------------------------------------------------*/
#include "tracing.h"
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
enum _FIELDTYP     fieldtyp;      /* typ of field */
int                ndis;          /* number of discretizations in this field */
struct _DISCRET   *dis;           /* structure holding a number of discretisations */
} FIELD;
/*----------------------------------------------------------------------*
 | this structure hold one discretization of a field      m.gee 2/02    |
 *----------------------------------------------------------------------*/
typedef struct _DISCRET
{
int                numnp;          /* number of nodes in this field */                          
int                numele;         /* number of elements in this field */
int                numeq;          /* number of unknowns */
int                numdf;          /* number of degrees of freedom */

struct _ELEMENT   *element;       /* vector of elements */
struct _NODE      *node;          /* vector of nodes */

int                ngnode;        /* number of geometry points */        
int                ngline;        /* number of geometry lines */        
int                ngsurf;        /* number of geometry surfaces */        
int                ngvol;         /* number of geometry volumes */        
struct _GNODE     *gnode;         /* vector of geometry points */    
struct _GLINE     *gline;         /* vector of geometry lines */    
struct _GSURF     *gsurf;         /* vector of geometry surfaces */    
struct _GVOL      *gvol;          /* vector of geometry volumes */    
} DISCRET;




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
int               restart;       /* is restart or not (not used yet) */

enum _PROBLEM_TYP probtyp;       /* type of problem, see enum.h */
enum _TIME_TYP    timetyp;       /* type of time, see enum.h */

} GENPROB;






/*----------------------------------------------------------------------*
 | Prototypes                                            m.gee 06/01    |
 | prototypes of all files using standardtypes.h                        |
 *----------------------------------------------------------------------*/
#include "prototypes.h"
/*----------------------------------------------------------------------*
 | Prototypes                                            m.gee 06/01    |
 | Prototypes of routines which include the                             |
 | header solution.h are kept separately in prototypes_sol.h            |
 | This is done, because of 2 reasons:                                  |
 | 1. The structures defined in solution.h are very big                 |
 | 2. in solution.h the headers of all kinds of solver libraries are    |
 |    included. As we do not know exactly, what is defined in these     |
 |    headers, it is better to keep the knowledge of these headers as   |
 |    small as possible                                                 |
 *----------------------------------------------------------------------*/





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
 extern struct _PAR   par;                      
 *----------------------------------------------------------------------*/
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
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 |                                                       m.gee 02/02    |
 | number of load curves numcurve                                       |
 | vector of structures of curves                                       |
 | defined in input_curves.c                                            |
 | int                   numcurve;                                      |
 | struct _CURVE      *curve;                                           |
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
