/*----------------------------------------------------------------------*
 | general structural dynamic-variables                   m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef union _ALLDYNA                 
{
   struct _STRUCT_DYNAMIC    *sdyn;   /* ptr for allocation of structural dynamic data */
   struct _FLUID_DYNAMIC     *fdyn;   /* ptr for allocation of fluid dynamic data */
   struct _FSI_DYNAMIC       *fsidyn; /*ptr fo allocation of fsi dynamic data */
   struct _ALE_DYNAMIC       *adyn;   /* ptr for allocation of ale dynamic data */
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
INT                updevry_disp;/* write result very updevry step */
INT                updevry_stress;/* write result very updevry step */
INT                res_write_evry;/* write restart every res_write_evry step */
INT                nstep;      /* number of steps */
INT                step;       /* actual step */
INT                damp;       /* flag to switch damping on/off */
INT                iter;       /* number of active iteration */
INT                maxiter;    /* maximum number of iterations */
INT                eigen;      /* flag for eigenvalue analysis */
INT                contact;    /* flag to switch contact onoff */
DOUBLE             toldisp;    /* displacement tolerance */
DOUBLE             dt;         /* stepsize */
DOUBLE             maxtime;    /* maximum total time */
DOUBLE             time;       /* actual time */
DOUBLE             beta;       /* time integration coefficients */
DOUBLE             gamma;
DOUBLE             alpha_m;
DOUBLE             alpha_f;
DOUBLE             m_damp;     /* factors for Raleigh damping */
DOUBLE             k_damp;

INT                timeadapt;  /* flag to switch adaptive time stepping on */
INT                itwant;     /* requested number of newton iterations */
DOUBLE             maxdt;      /* max allowed time step */
DOUBLE             resultdt;   /* postprocessing time step */
INT                writecounter; /* counter for output steps */
} STRUCT_DYNAMIC;
/*----------------------------------------------------------------------*
 | general dynamic-variables for analysis                 m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _STRUCT_DYN_CALC                
{
DOUBLE             rldfac;          /* ? */
DOUBLE             rnorm;           /* ? */

DOUBLE             constants[20];   /* derived constants for time integration */


DOUBLE             epot;            /* potential energy */
DOUBLE             eout;            /* external energy */
DOUBLE             etot;            /* total energy */
DOUBLE             ekin;            /* kinetic energy */

DOUBLE             dinorm;  /* square of the L2-norm of the residual displacements */
DOUBLE             dnorm;   /* square of the L2-norm of the displacements increment */

} STRUCT_DYN_CALC;
/*----------------------------------------------------------------------*
 | general fsi variables                                  genk 09/02    |
 *----------------------------------------------------------------------*/
typedef struct _FSI_DYNAMIC                 
{
INT                ifsi;            /* coupling algorithm */
INT                ipre;            /* type of predictor */
INT                inrmfsi;         /* convergence criterion */
INT                ichecke;         /* energy check */
INT                inest;           /* nested iteration */
INT                ichopt;          /* optimal ordering for CHEBYCHEV parameter */
INT                iait;            /* Aitken iteration */
INT                itechapp;        /* No. of Iter. for approx. EW-Calculation */
INT                ichmax;          /* Max. No. of CHEBYCHEV iterations */
INT                isdmax;          /* Max. No. of steepest descent iterations */
INT                nstep;           /* number of steps */
INT                itemax;          /* max. number of iterations over fields */
INT                uppss;           
INT                step;
INT                iale;            
DOUBLE             time;
DOUBLE             dt;              /* time increment */
DOUBLE             maxtime;         /* total time */
DOUBLE             entol;           /* tolerance for energy check over fields */
DOUBLE             relax;           /* actual relaxation parameter */
DOUBLE             convtol;         /* tolerance for iteration over fields */ 
DOUBLE             deltaeint;       /* energy production at the interface */
ARRAY              sid;             /* structural interface dofs */
INT                numsid;          /* number of structural interface dofs */
INT                actpos;          /*  */
} FSI_DYNAMIC;

/*----------------------------------------------------------------------*
 | general ale dynamic variables                            ck 12/02    |
 *----------------------------------------------------------------------*/
typedef struct _ALE_DYNAMIC                
{
enum
   {
                  classic_lin
   } typ;                        /* switch dynamic algorithm */
INT                nstep;        /* number of steps */
INT                step;         /* actual step */
INT                updevry_disp; /* write result very updevry step */
DOUBLE             dt;           /* stepsize */
DOUBLE             maxtime;      /* maximum total time */
DOUBLE             time;         /* actual time */
} ALE_DYNAMIC;
