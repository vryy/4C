/*----------------------------------------------------------------------*
 | general structural dynamic-variables                   m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef union _ALLDYNA                 
{
   struct _STRUCT_DYNAMIC    *sdyn; /* ptr for allocation of structural dynamic data */
   struct _FLUID_DYNAMIC     *fdyn; /* ptr for allocation of fluid dynamic data */
} ALLDYNA;



/*----------------------------------------------------------------------*
 | general structural dynamic-variables                   m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _STRUCT_DYNAMIC                 
{
char               dyntyp[50]; /* type of algorithm */
int                updevry_disp;/* write result very updevry step */
int                updevry_stress;/* write result very updevry step */
int                nstep;      /* number of steps */
int                step;       /* actual step */
int                damp;       /* flag to switch damping on/off */
int                iter;       /* number of active iteration */
int                maxiter;    /* maximum number of iterations */
double             toldisp;    /* displacement tolerance */
double             dt;         /* stepsize */
double             maxtime;    /* maximum total time */
double             time;       /* actual time */
double             beta;       /* time integration coefficients */
double             gamma;
double             alpha_m;
double             alpha_f;
double             m_damp;     /* factors for Raleigh damping */
double             k_damp;
} STRUCT_DYNAMIC;
/*----------------------------------------------------------------------*
 | general dynamic-variables for analysis                 m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _STRUCT_DYN_CALC                
{
double             rldfac;          /* ? */
double             rnorm;           /* ? */

double             constants[20];   /* derived constants for time integration */


double             epot;            /* potential energy */
double             eout;            /* external energy */
double             etot;            /* total energy */
double             ekin;            /* kinetic energy */

double             dinorm;  /* square of the L2-norm of the residual displacements */
double             dnorm;   /* square of the L2-norm of the displacements increment */

} STRUCT_DYN_CALC;



/*----------------------------------------------------------------------*
 | general fluid dynamic-variables                        m.gee 2/02    |
 *----------------------------------------------------------------------*/
typedef struct _FLUID_DYNAMIC                 
{
int                i;          /* not used at the moment */
} FLUID_DYNAMIC;
