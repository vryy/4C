/*!---------------------------------------------------------------------
\file
\brief contains bug and time tracing structures

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

---------------------------------------------------------------------*/
/*!
\addtogroup DSSYSTEM
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief tracing of routines

<pre>                                                         m.gee 8/00
This structures is used by the chained list that keeps track of
the function calls.
This chained list is organized as a ring of lenght 100 and is readily
preallocated. It can therefore trace routine calls up to a deepness
level of 100 routines before it starts overriding itself
</pre>

*----------------------------------------------------------------------*/
typedef struct _TRACEROUT
{
struct _TRACEROUT   *prev;          /*!< ptr to previous structure in chained list */
struct _TRACEROUT   *next;          /*!< ptr to next structure in chained list */
char                *name;          /*!< name of routine */
enum
   {
    dsnone,
    dsin,
    dsout
   }                 dsroutcontrol; /*!< status of routine, inside, outside, unknown */

} TRACEROUT;

/*!----------------------------------------------------------------------
\brief tracing of ARRAYs

<pre>                                                         m.gee 2/02
This structures is used by the chained list that keeps track of
the the ARRAYS which are allocated using the AM-System
The chained list is fully dynamic an creates one structure to point
to each ARRAY or ARRAY4D which is created. If an ARRAY is deleted
using amdel or am4del the structure is taken off the chain list and
deallocated.
A report about all ARRAYs and ARRAY4Ds can be written to .err to keep
to e.g. detect local or damaged ARRAYs which where not destroyed
</pre>

*----------------------------------------------------------------------*/
typedef struct _TRACEARRAY
{
   struct _TRACEARRAY   *prev;     /*!< ptr to previous structure in chained list */
   struct _TRACEARRAY   *next;     /*!< ptr to next structure in chained list */
   enum
   {
      array_none,
      array_2d,
      array_4d
   }                     arraytyp; /*!< type of array traced by this structure */
   union
   {
      struct _ARRAY      *a2;      /*!< ptr to the 2D array */
      struct _ARRAY4D    *a4;      /*!< ptr to the 3D array */
   }                     a;        /*!< name of union */
} TRACEARRAY;



/*!----------------------------------------------------------------------
\brief the tracing variable

<pre>                                                         m.gee 8/00
defined in pss_ds.c, declared in tracing.h
</pre>
*----------------------------------------------------------------------*/
typedef struct _TRACE
{
/* variables for watching the ARRAYS */
INT                 trace_on;             /*!< switches trace on/off */
INT                 num_arrays;           /*!< number of current INT-arrays */

struct _TRACEARRAY *arraychain;           /*!< start of the linear chained list */
struct _TRACEARRAY *endarraychain;        /*!< ptr to the actual end of the chain list */

/* variables for watching the routines */
INT                 deepness;             /*!< the actual deepness of the calling tree */
struct _TRACEROUT   routine[100];         /*!< chained list ring to trace routine */
struct _TRACEROUT   *actroutine;          /*!< ptr to actual routine */
} TRACE;


/*! @} (documentation module close)*/
