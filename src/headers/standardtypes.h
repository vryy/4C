/*!---------------------------------------------------------------------
\file
\brief 

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

---------------------------------------------------------------------*/

#ifndef STANDARDTYPES_H
#define STANDARDTYPES_H

/*----------------------------------------------------------------------*
 | includes of Ansi C standard headers                   m.gee 8/00    |
 *----------------------------------------------------------------------*/
#include "math.h"
#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "ctype.h"
#include "time.h"
/*----------------------------------------------------------------------*
 | includes of mpi of C                                   m.gee 8/00    |
 *----------------------------------------------------------------------*/
#ifdef PARALLEL 
   #include "mpi.h"
#endif
/*----------------------------------------------------------------------*
 | check_defines file                                       mn 01/04    |
 *----------------------------------------------------------------------*/
#include "check_defines.h"
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

/*----------------------------------------------------------------------*/
#ifdef WALLCONTACT
#include "../wall1/wallcontact.h"
#endif
/*----------------------------------------------------------------------*/

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
 | structures used for eigensolution analysis              al 08/02     |
 *----------------------------------------------------------------------*/
#include "eigenanalysis.h"
/*----------------------------------------------------------------------*
 | structures used for dynamic analysis                   m.gee 2/02    |
 *----------------------------------------------------------------------*/
#include "dynamic.h"
/*----------------------------------------------------------------------*
 | structures used for fluid dynamic analysis             genk 05/02    |
 *----------------------------------------------------------------------*/
#include "fluid.h"
/*----------------------------------------------------------------------*
 | structures used for bug and time tracing               m.gee 2/02    |
 *----------------------------------------------------------------------*/
#include "tracing.h"
/*----------------------------------------------------------------------*
 | structures used for bug and time tracing               m.gee 10/02   |
 *----------------------------------------------------------------------*/
#include "container.h"
/*----------------------------------------------------------------------*/




/*! 
\addtogroup FRSYSTEM 
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h                                                  
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>

*----------------------------------------------------------------------*/
typedef struct _FILES
{
/*------------------------------------------------------------- file I/O */
char             *inputfile_name;         /* input file name             */

char             *outputfile_kenner;      /* output file kenner          */
char              outputfile_name[100];   /* output file name            */
char             *vispssfile_kenner;
char              vispssfile_name[100];
char             *visresfile_kenner;
char              visresfile_name[100];
size_t            outlenght;              /* lenght of output file kenner*/
size_t            vispsslength;
size_t            visreslength;
INT               num_outputfiles;        /* number of output files      */
INT               pss_counter;            /* number of records on pss-file */
INT               num_flaviaresfiles;     /*                               */

FILE             *in_input;               /* file-pointer input file     */
FILE             *out_out;                /* file-pointer .out  file     */
FILE             *out_err;                /* file-pointer .err  file     */
FILE             *out_pss;                /* file ptr to restart-pss file */
FILE             *in_pss;                 /* file-pointer .pss  file     */
FILE             *out_mon;                /* file-pointer .mon file     */
#ifdef GEMM
FILE             *out_gemm;               /* file-pointer .gemm file*/
#endif
FILE             *out_tur;                /* file-pointer .tur  file     */

FILE             *gidmsh;                 /* file pointer .flavia.msh    */
FILE             *gidres;                 /* file pointer .flavia.res    */
FILE             *gnu;                    /* file pointer .plt */

/*---------------------- variables needed by the free-field-input system */
char              title[5][500];          /* problem title                */
char              line[500];              /* linebuffer for reading       */
char            **input_file;             /* copy of the input file in rows and columns */
char             *input_file_hook;        /* ptr the copy of the input file is allocated to */
INT               numcol;                 /* number of cols in inputfile  */
INT               numrows;                /* number of rows in inputfile  */
INT               actrow;                 /* rowpointer used by fr        */
char             *actplace;               /* pointer to actual place in input-file */
} FILES;
/*! @} (documentation module close)*/





/*----------------------------------------------------------------------*
 | general IO-flags                                      m.gee 12/01    |
 | flags to switch the output of certain data on or off, read from gid  |
 *----------------------------------------------------------------------*/
typedef struct _IO_FLAGS
{
INT               struct_disp_file;       /* write displacements to .out */
INT               struct_stress_file;     /* write structural stress to .out */
INT               struct_disp_gid;        /* write structural displacements to .flavia.res */
INT               struct_stress_gid;      /* write structural stresses to .flavia.res */
INT               struct_stress_gid_smo;  /* write structural smoothed/unsmoothed stresses to .flavia.res */
INT               fluid_sol_file;         /* write vel/pre to .out */
INT               fluid_sol_gid;          /* write vel/pre to .flavia.res */
INT               fluid_stress_gid;       /* write stresses to .flavia.res */
INT               fluid_vis_file;         /* write solution to pss-file for VISUAL2 */
INT               ale_disp_file;          /* write ale displacements to .out */
INT               ale_disp_gid;           /* write ale displacement to .flavia.res */
INT               monitor;
INT               relative_displ;         /* write relative displacements to .err */
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
INT                ndis;          /* number of discretizations in this field */
struct _DISCRET   *dis;           /* structure holding a number of discretisations */
} FIELD;
/*----------------------------------------------------------------------*
 | this structure hold one discretization of a field      m.gee 2/02    |
 *----------------------------------------------------------------------*/
typedef struct _DISCRET
{
INT                numnp;          /* number of nodes in this field */ 
INT                numele;         /* number of elements in this field */
INT                numeq;          /* number of unknowns */
INT                numdf;          /* number of degrees of freedom */

struct _ELEMENT   *element;       /* vector of elements */
struct _NODE      *node;          /* vector of nodes */

INT                ngnode;        /* number of geometry points */        
INT                ngline;        /* number of geometry lines */        
INT                ngsurf;        /* number of geometry surfaces */        
INT                ngvol;         /* number of geometry volumes */        
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
INT               nele;          /* total number of elements over all fields */
INT               nnode;         /* total number of nodes over all fields*/
INT               maxnode;       /* largest node id in the problem */
INT               nodeshift;     /* shift of node ids for second dis */
INT               ndim;          /* dimension of problem (2 or 3) */
INT               nmat;          /* total number of material laws */
INT               numfld;        /* number of fields */
INT               numdf;         /* maximum number of dofs to one node (not used, in progress)*/
INT               restart;       /* is restart or not */
INT               visual;        /* flag for visualise mode or not */
INT               numsf;         /* actual number of struct-field */
INT               numff;         /* actual number of fluid field */
INT               numaf;         /* actual number of ale field */
INT               numls;         /* actual number of ls field */
INT               graderw;       /* flag is gradient enhanced material model */

#ifdef D_XFEM
INT               xfem_on_off;    /* flag to switch enriched formulation */
#endif
  
enum _PROBLEM_TYP probtyp;       /* type of problem, see enum.h */
enum _TIME_TYP    timetyp;       /* type of time, see enum.h */

#ifdef RESULTTEST
  INT             numresults;   /* number of known results waiting to
                                 * be compared to the calculated ones. */
#endif
  
NODE            **nodes;
INT               ngnode;
GNODE           **gnodes;
INT               ngline;
GLINE           **glines;
INT               ngsurf;
GSURF           **gsurfs;
INT               ngvol;
GVOL            **gvols;
} GENPROB;

/*---------------------------------------------------------------------*
 | monotoring informations                                  genk 01/03 |
 *---------------------------------------------------------------------*/
typedef struct _MONITOR
{
INT              numnp;
INT              numval;
INT              numstep;
ARRAY            onoff;
ARRAY            monnodes;
ARRAY            val;
ARRAY            time;
} MONITOR;


#ifdef RESULTTEST

/*
 * The results we expect. For testing.
 * Each record describes one value and where to find it.
 *
 * Be careful: Testing is only legal running one processor.
 */
typedef struct _RESULTDESCR
{
  FIELDTYP field;               /* the field the value steams from */
  INT dis;                      /* the discretisation */
  
  /* You can only have one. Either node or element. The other one must
   * equal -1. */
  INT node;                     /* the node the value belongs to */
  INT element;                  /* or the element */
  
  CHAR position[100];           /* The position of the checked value
                                 * inside its node or
                                 * element. Normally it looks like
                                 * "sol(row,col)" but for elements it
                                 * could be anything. */
  CHAR name[25];                /* the physical quantity, used for output */
  DOUBLE value;                 /* the expected value  */
  DOUBLE tolerance;             /* an exceptable tolerance */
} RESULTDESCR;

#endif

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


/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h                                                  
</pre>

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


/*----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 06/01
This structure struct _PAR par; is defined in input_control_global.c
and the type is in standardtypes.h                                                  
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
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
/*!----------------------------------------------------------------------
\brief the tracing variable

<pre>                                                         m.gee 8/00
defined in pss_ds.c, declared in tracing.h                                                  
</pre>
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
 | INT                   numcurve;                                      |
 | struct _CURVE      *curve;                                           |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 |                                                          mn 02/04    |
 | number of spatual load functions   numcurve                          |
 | vector of structures of functions                                    |
 | defined in input_funct.c                                             |
 | INT                 numfunct;                                        |
 | struct _FUNCT      *funct;                                           |
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
/*----------------------------------------------------------------------*
 |                                                          uk 06/04    |
 | An array of expected results. This is read from the input file.      |
 | defined in global_result_test.c                                      |
 | struct _RESULTDESCR      *resultdescr;                               |
 *----------------------------------------------------------------------*/

#endif
