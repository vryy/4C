/*----------------------------------------------------------------------*
 | tracing of time & array bugs                           m.gee 8/00    |
 | This structures is used by the chained list that keeps track of      |
 | the function calls.                                                  |
 | This chained list is organized as a ring of lenght 100 and is readily|
 | preallocated. It can therefor trace routine calls up to a deepness   |
 | level of 100 routines before it starts overriding itself             |
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
 | tracing of time & array bugs                           m.gee 2/02    |
 | This structures is used by the chained list that keeps track of      |
 | the the ARRAYS which are allocated using the AM-System               |
 | The chained list is fully dynamic an creates one structure to point  |
 | to each ARRAY or ARRAY4D which is created. If an ARRAY is deleted    |
 | using amdel or am4del the structure is taken off the chain list and  |
 | deallocated.                                                         |
 | A report about all ARRAYs and ARRAY4Ds can be written to .err to keep|
 | to e.g. detect local or damaged ARRAYs which where not destroyed     |
 |                                                                      |
 *----------------------------------------------------------------------*/
typedef struct _TRACEARRAY
{
   struct _TRACEARRAY   *prev;
   struct _TRACEARRAY   *next;
   enum
   {
      array_none,
      array_2d,
      array_4d
   }                     arraytyp;
   union
   {
      struct _ARRAY      *a2;
      struct _ARRAY4D    *a4;
   }                     a;

} TRACEARRAY;



/*----------------------------------------------------------------------*
 | tracing of time & array bugs                           m.gee 8/00    |
 |                                                                      |
 *----------------------------------------------------------------------*/
typedef struct _TRACE
{
/* variables for watching the ARRAYS */
int                 trace_on;             /* switches trace on/off */
int                 num_arrays;           /* number of current int-arrays */

struct _TRACEARRAY *arraychain;           /* start of the linear chained list */
struct _TRACEARRAY *endarraychain;        /* ptr to the actual end of the chain list */

/* variables for watching the routines */
int                 deepness;             /* the actual deepness of the calling tree */
struct _TRACEROUT   routine[100];         /* chained list ring to trace routine */
struct _TRACEROUT   *actroutine;          /* ptr to actual routine */
} TRACE;
