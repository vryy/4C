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
/*----------------------------------------------------------------------*
 | general structural dynamic-variables                   m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef union _ALLDYNA
{
   struct _STRUCT_DYNAMIC    *sdyn;   /* ptr for allocation of structural dynamic data */
   struct _FLUID_DYNAMIC     *fdyn;   /* ptr for allocation of fluid dynamic data */
   struct _FSI_DYNAMIC       *fsidyn; /*ptr for allocation of fsi dynamic data */
   struct _SSI_DYNAMIC       *ssidyn; /*ptr for allocation of ssi dynamic data */
   struct _ALE_DYNAMIC       *adyn;   /* ptr for allocation of ale dynamic data */
#ifdef D_LS
  struct _LS_DYNAMIC        *lsdyn;  /* ptr for allocation of ls dynamic data */
#endif
} ALLDYNA;



/*----------------------------------------------------------------------*
 | general structural dynamic-variables                   m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _STRUCT_DYNAMIC
{
enum
   {
    gen_alfa,
    centr_diff,
    Gen_EMM
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
#ifdef GEMM
DOUBLE             xsi;        /*  Parameter used by GEMM */
#endif
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
#ifdef D_WALL1
DOUBLE             total_linmom[2];
DOUBLE             total_angular_momentum;
DOUBLE             total_strain_energy;
DOUBLE             total_kinetic_energy;
DOUBLE             local_linmom[2];
DOUBLE             local_angular_momentum;
DOUBLE             local_strain_energy;
DOUBLE             local_kinetic_energy;
#endif
DOUBLE             dinorm;  /* square of the L2-norm of the residual displacements */
DOUBLE             dnorm;   /* square of the L2-norm of the displacements increment */

} STRUCT_DYN_CALC;
/*----------------------------------------------------------------------*
 | general fsi variables                                  genk 09/02    |
 *----------------------------------------------------------------------*/
typedef struct _FSI_DYNAMIC
{
INT                ifsi;            /*!< coupling algorithm */
INT                ipre;            /*!< type of predictor */
INT                inrmfsi;         /*!< convergence criterion */
INT                ichecke;         /*!< energy check */
INT                isdmax;          /*!< Max. No. of steepest descent iterations */
INT                nstep;           /*!< number of steps */
INT                itemax;          /*!< max. number of iterations over fields */
INT                uppss;           /*!<  */
INT                upres;           /*!< update .flavia.res every step */
INT                uprestart;       /*!< write restart every step */
INT                step;            /*!<  */
INT                iale;            /*!<  */
DOUBLE             time;            /*!<  */
DOUBLE             dt;              /*!< time increment */
DOUBLE             maxtime;         /*!< total time */
DOUBLE             entol;           /*!< tolerance for energy check over fields */
DOUBLE             relax;           /*!< actual relaxation parameter */
DOUBLE             convtol;         /*!< tolerance for iteration over fields */
DOUBLE             deltaeint;       /*!< energy production at the interface */
ARRAY              sid;             /*!< structural interface dofs */
INT                numsid;          /*!< number of structural interface dofs */
INT                actpos;          /*!<  */
INT                coupmethod;      /*!< flag, 0=mortar , 1=conforming */
} FSI_DYNAMIC;

/*----------------------------------------------------------------------*
 | general ssi variables                                  genk 10/03    |
 *----------------------------------------------------------------------*/
typedef struct _SSI_DYNAMIC                 
{
INT                ifsi;            /*!< coupling algorithm */
INT                ipre;            /*!< type of predictor */
INT                inrmfsi;         /*!< convergence criterion */
INT                ichecke;         /*!< energy check */
INT                inest;           /*!< nested iteration */
INT                ichopt;          /*!< optimal ordering for CHEBYCHEV parameter */
INT                iait;            /*!< Aitken iteration */
INT                itechapp;        /*!< No. of Iter. for approx. EW-Calculation */
INT                ichmax;          /*!< Max. No. of CHEBYCHEV iterations */
INT                isdmax;          /*!< Max. No. of steepest descent iterations */
INT                nstep;           /*!< number of steps */
INT                itemax;          /*!< max. number of iterations over fields */
INT                uppss;           /*!<  */
INT                upres;           /*!< update .flavia.res every step */
INT                res_write_evry;  /*!< write restart every step */
INT                step;            /*!<  */
INT                iale;            /*!<  */
DOUBLE             time;            /*!<  */
DOUBLE             dt;              /*!< time increment */
DOUBLE             maxtime;         /*!< total time */
DOUBLE             entol;           /*!< tolerance for energy check over fields */
DOUBLE             relax;           /*!< actual relaxation parameter */
DOUBLE             convtol;         /*!< tolerance for iteration over fields */ 
DOUBLE             deltaeint;       /*!< energy production at the interface */
ARRAY              sid;             /*!< structural interface dofs */
INT                numsid;          /*!< number of structural interface dofs */
INT                actpos;          /*!<  */
INT                conformmesh;     /*!< flag, 0=conf. discr., 1=nonconf. discr. */
INT                coupmethod;      /*!< flag, 0=mortar , 1=interpolation */
} SSI_DYNAMIC;


/*----------------------------------------------------------------------*
 | general ale dynamic variables                            ck 12/02    |
 *----------------------------------------------------------------------*/
typedef struct _ALE_DYNAMIC
{
enum
   {
                  classic_lin,   /*!< classic linear calculation */
		  min_Je_stiff,  /*!< incremental calculation
		                    stiffened with min J_element^2 */
                  two_step,      /*!< calculation in 2 steps per timestep */
		  springs,       /*!< springs rather than continous pseudo material */
		  laplace,       /*!< Laplace smoothing algorithm */
                  LAS            /*!< large amplitude sloshing */
   } typ;                        /*!< switch dynamic algorithm */

enum
   {
             no_quality,         /*!< no element quality monitoring */
             aspect_ratio,       /*!< aspect ratio element quality criterion */
             corner_angle,       /*!< corner angle element quality criterion */
	     min_detF            /*!< minimal elemental Jacobian determinant */
   } measure_quality;            /*!< switch, which quality measure to watch */

INT                nstep;        /*!< number of steps */
INT                step;         /*!< actual step */
INT                updevry_disp; /*!< write result very updevry step */
INT                num_initstep; /*!< number of initial steps with prestress */
DOUBLE             dt;           /*!< stepsize */
DOUBLE             maxtime;      /*!< maximum total time */
DOUBLE             time;         /*!< actual time */
INT                coupmethod;   /*!< flag, 0=mortar , 1=conforming */
} ALE_DYNAMIC;
