/*----------------------------------------------------------------------*
 | general structural dynamic-variables                   m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef union _DYNAMIC                 
{
   struct _STRUCT_DYNAMIC    *sdyn;
   struct _FLUID_DYNAMIC     *fdyn;
} DYNAMIC;



/*----------------------------------------------------------------------*
 | general structural dynamic-variables                   m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _STRUCT_DYNAMIC                 
{
char               dyntyp[50];
int                nstep;    /* this all is in progress... */
int                damp;     
int                iter;
int                maxiter;
double             toldisp;
double             dt;
double             maxtime;
double             beta;
double             gamma;
double             alpha_m;
double             alpha_f;
double             m_damp;
double             k_damp;
} STRUCT_DYNAMIC;
/*----------------------------------------------------------------------*
 | general dynamic-variables for analysis                 m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _STRUCT_DYN_CALC                
{
double             rldfac;
double             rnorm;

double             constants[20];

double             epot;
double             eout;
double             etot;
double             ekin;
} STRUCT_DYN_CALC;



/*----------------------------------------------------------------------*
 | general fluid dynamic-variables                        m.gee 2/02    |
 *----------------------------------------------------------------------*/
typedef struct _FLUID_DYNAMIC                 
{
int                i;          /* not used at the moment */
} FLUID_DYNAMIC;
