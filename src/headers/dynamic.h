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
enum
   {
    gen_alfa,
    centr_diff
   }               Typ;         /* type of time integration algorithm */
int                updevry_disp;/* write result very updevry step */
int                updevry_stress;/* write result very updevry step */
int                res_write_evry;/* write restart every res_write_evry step */
int                nstep;      /* number of steps */
int                step;       /* actual step */
int                damp;       /* flag to switch damping on/off */
int                iter;       /* number of active iteration */
int                maxiter;    /* maximum number of iterations */
int                eigen;      /* flag for eigenvalue analysis */
int                contact;    /* flag to switch contact onoff */
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

int                timeadapt;  /* flag to switch adaptive time stepping on */
int                itwant;     /* requested number of newton iterations */
double             maxdt;      /* max allowed time step */
double             resultdt;   /* postprocessing time step */
int                writecounter; /* counter for output steps */
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
