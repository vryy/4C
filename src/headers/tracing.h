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
