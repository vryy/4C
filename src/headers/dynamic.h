#include "../fluid_full/fluid.h"

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
 | general fluid dynamic-variables for element evaluation               |
 |                                                        genk 03/02    |
 *----------------------------------------------------------------------*/
typedef struct _FLUID_DYN_CALC                
{
double dta;      /* actual time increment */
double thsl;     /* theta-s,l: const. for "stiffness" terms LHS */
double thsr;     /* theta-s,r: const. for "stiffness" terms RHS */
double thpl;     /* theta-p,l: const. for "pressure" terms LHS  */
double thpr;     /* theta-p,r: const. for "pressure" terms RHS  */
double velmax;   /* max. velocity, needed for stabilisaton parameter */
double tau[3];   /* array for stabilitity parameter */
double sigma;    /* const. for nonlinear iteration */   
int    itwost;   /* control variable for element evaluation */
int    isemim;   /* control variable for element evaluation */
int    iprerhs;  /* treatment of pressure in time discr. */
int    nik;
int    nic;
int    nir;
int    nie;
int    nil;
int    nif;
int    nii;
int    nis;     /* flags for nonlinear iteration */
int    ishape;  /* flag for new element shape */
union  _FLUID_DATA data;
} FLUID_DYN_CALC;


/*----------------------------------------------------------------------*
 | general fluid dynamic-variables from input                           |
 |                                                         genk 3/02    |
 *----------------------------------------------------------------------*/
typedef struct _FLUID_DYNAMIC    /* this is all in progress*/             
{
char               dyntyp[50];   /* dynamictype */
int                numdf;        /* number of dofs of the fluid elements */
int                iopfsi;       /* time integration method */
int                numcont;      /* number of continuation steps */
int                uppss;        /* update pss file every n steps */
int                idisp;        /* store results every n steps */      
int                nstep;        /* number of timesteps */
int                step;         
int                ite;          /* nonlinear iteration scheme */
int                itemax;       /* number of nonlin. iterations */
int                itchk;        /* convergence check during nonlin. iteration */
int                itnorm;       /* norm for conv. check d. nonlin. iteration */
int                stchk;        /* steady state check every n steps */
int                stnorm;       /* norm for steady state check */
int                iopfss;       /* starting algorithm */
int                numfss;       /* number of starting algorithm steps */
int                init;         /* initialisation of starting field */
int                iprerhs;      /* treatment of pressure in time discr. */
double             maxtime;      /* maximal simulation time */
double             time;         /* actual time */
double             dt;           /* time increment */
double             alpha;        /* time integration constant */
double             theta;        /* time integration constant */
double             gamma;        /* time integration constant */
double             ittol;        /* tolerance for iteration convergence check */
double             sttol;        /* tolerance for steady state check */
double             thetas;       /* constant for starting algorithm) */
struct _ARRAY      start;        /* starting field */
struct _FLUID_DYN_CALC dynvar;
} FLUID_DYNAMIC;
