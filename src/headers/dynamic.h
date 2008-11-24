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
#ifndef CCADISCRET
   struct _STRUCT_DYNAMIC    *sdyn;   /* ptr for allocation of structural dynamic data */
#endif
   struct _FLUID_DYNAMIC     *fdyn;   /* ptr for allocation of fluid dynamic data */
   struct _FSI_DYNAMIC       *fsidyn; /* ptr for allocation of fsi dynamic data */
   struct _SSI_DYNAMIC       *ssidyn; /* ptr for allocation of ssi dynamic data */
   struct _ALE_DYNAMIC       *adyn;   /* ptr for allocation of ale dynamic data */
   struct _THERM_DYNAMIC     *tdyn;  /* ptr for allocation of THERMAL
                                          dynamic control data */
   struct _TSI_DYNAMIC       *tsidyn;  /* ptr for allocation of TSI
                                          dynamic control data */
} ALLDYNA;


/*----------------------------------------------------------------------*
 | time adaptivity only for read in                     bborn 10/07     |
 *----------------------------------------------------------------------*/
#ifndef CCADISCRET
typedef struct _TIMADA_DYNAMIC
{
  enum _timadakindenum {
    timada_kind_none,           /* no time adaptivity */
    timada_kind_zienxie,        /* Zienkiewicz-Xie indicator */
    timada_kind_ab2             /* Adams-Bahsforth2 indicator */
  } kind;                       /* type of adaptivity in time */
  DOUBLE dt_max;                /* maximally permitted step size */
  DOUBLE dt_min;                /* minimally permitted step size */
  DOUBLE dt_scl_min;            /* minimally permitted ratio of
                                 * new to last size */
  DOUBLE dt_scl_max;            /* maximally permitted ratio of
                                 * new to last size */
  DOUBLE dt_scl_saf;            /* safety scale of optimally predicted
                                 * new step size */
  enum {
    timada_err_norm_vague = 0,  /* undetermined norm */
    timada_err_norm_l1,         /* L1/linear norm */
    timada_err_norm_l2,         /* L2/Euclidean norm */
    timada_err_norm_rms,        /* root mean square (RMS) norm */
    timada_err_norm_inf         /* Maximum/infinity norm */
  } err_norm;                   /* error norm */
  DOUBLE err_tol;               /* error tolerance (target) */
  /* INT err_pow; */            /* order */
  INT adastpmax;                /* maximally permitted adaptive
                                 * step size iterations */
} TIMADA_DYNAMIC;
#endif


/*----------------------------------------------------------------------*
 | general structural dynamic-variables                   m.gee 4/01    |
 *----------------------------------------------------------------------*/
#ifndef CCADISCRET
typedef struct _STRUCT_DYNAMIC
{
enum
   {
    gen_alfa,                   /* generalised-alpha time integrator */
    gen_alfa_statics,           /* static analysis with gen_alfa */
    centr_diff,                 /* central differences (explicit) */
    Gen_EMM,                    /* generalised energy-momentum method */
    statics,                    /* static analysis */
    genalpha,                   /* generalised-alpha time integrator, 
                                 * new style (implicit) */
    onesteptheta,               /* one-step-theta time integrator (implicit) */
    gemm,                       /* generalised energy-momementum method */
    ab2                         /* Adams-Bashforth 2nd order (explicit) */
   }               Typ;         /* type of time integration algorithm */
enum _genavgtype
   {
    genavg_vague,               /* undefined mid-averaging type */
    genavg_imrlike,             /* alphaf-mid-averaging is done IMR-like, i.e.
                                 *    F_{int,m} = F_{int}(D_m)
                                 *              = F_{int}( (1-alpha_f)*D_{n+1} + alpha_f*D_n )
                                 * (IMR means implicit mid-point rule.) */
    genavg_trlike               /* alphaf-mid-averaging is done TR-like, i.e.
                                 *    F_{int,m} = (1-alpha_f)*F_{int,n+1} + alpha_f*F_{int,n}
                                 *              = (1-alpha_f)*F_{int}(D_{n+1}) + alpha_f*F_{int}(D_n)
                                 * (TR means trapezoidal rule.) */
   }               genavgtype;  /*  alphaf-mid-average type for generalised-alpha time integration */
enum _nlnSolvTyp
   {
    vague,                      /* undetermined sol.tech. type */
    fullnewton,                 /* full Newton-Raphson */
    lsnewton,                   /* line search Newton-Raphson */
    modnewton,                  /* modified Newton */
    matfreenewton,              /* matrix-free Newton iteration */
    nlncg,                      /* nonlinear CG iteration using Nox */
    ptc,                        /* pseudo transient continuation nonlinear iteration */
    newtonlinuzawa,             /* linear Uzawa solution for constraint system */
    augmentedlagrange,          /* non-linear Uzawa solution for constraint system */
    noxnewtonlinesearch         /* Line search Newton utilising NOX */
   }               nlnSolvTyp;  /* type of nonlinear solver to be used */
enum _PredType
   {
     pred_vague,                /* undetermined */
     pred_constdis,             /* constant displacements */
     pred_constdisvelacc        /* constant displacements, velocities and accelerations */
   }               predtype;    /* predictor type */
enum _convcheck
   {
     absres_or_absdis,          /* absolute norms of residual forces
                                 * OR iterative displacement increments */
     absres_and_absdis,         /* absolute norms of residual forces
                                 * AND iterative displacement increments */
     relres_or_absdis,          /* relative norm of residual forces
                                 * OR absolute norm if iterative
                                 * displacement increments */
     relres_and_absdis,         /* relative norm of residual forces
                                 * AND absolute norm if iterative
                                 * displacement increments */
     relres_or_reldis,          /* relative norms of residual forces
                                 * OR iterative displacement increments */
     relres_and_reldis,         /* relative norms of residual forces
                                 * AND iterative displacement increments */
     linuzawa,                  /* linear Uzawa for constraint system */
     nonlinuzawa                /* non-linear Uzawa for constraint system */
   }               convcheck;   /* convergence check of solution technique */
enum _dynkind
   {
     dynkind_deprecated,        /* see over left */
     dynkind_direct,            /* direct TIS */
     dynkind_directadaptive,    /* direct TIS with adaptivity */
     dynkind_invanalysis        /* inverse analysis */
   }               dynkind;     /* kind of dynamic analysis */
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
/* Generalised-alpha parameters */
DOUBLE             beta;       /* time integration coefficients */
DOUBLE             gamma;
DOUBLE             alpha_m;
DOUBLE             alpha_f;
/* Generalised energy-momentum parameter */
#ifdef GEMM
DOUBLE             xsi;        /*  Parameter used by GEMM */
#endif
/* damping */
enum _dampkind
   {
     damp_none = 0,            /* damping off */
     damp_rayleigh,            /* globally applied Rayleigh damping */
     damp_material             /* element-wise applied damping */
   }               dampkind;   /* type of damping */
DOUBLE             m_damp;     /* factor for Raleigh damping */
DOUBLE             k_damp;     /* factor for Raleigh damping */
/* time step size adaptivity --- depreciated??? */
INT                timeadapt;  /* flag to switch adaptive time stepping on */
INT                itwant;     /* requested number of newton iterations */
DOUBLE             maxdt;      /* max allowed time step */
DOUBLE             resultdt;   /* postprocessing time step */
/* time step size adaptivity --- new style */
#ifndef CCADISCRET
TIMADA_DYNAMIC     timada;     /* time adaptivity data */
#endif
/* output */
INT                writecounter; /* counter for output steps */
} STRUCT_DYNAMIC;
#endif
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
FSI_COUPLING       ifsi;            /*!< coupling algorithm */
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
enum {
      cf_none,       /*! No evaluation of coupling force                */
      cf_stress,     /*! Evaluation using stress values (derivatives)   */
      cf_nodeforce   /*! Evaluation via consistent nodal forces         */
     } coupforce;    /*!< how to calculate fsi coupling force */
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
                  incr_lin,
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
INT                coupmethod;   /*!< flag, 0=mortar , 1=conforming, 2=nonconforming */
} ALE_DYNAMIC;


/*----------------------------------------------------------------------*
 | general thermal dynamic control                        bborn 03/06   |
 *----------------------------------------------------------------------*/
typedef struct _THERM_DYNAMIC
{
  /* time thingies */
  DOUBLE             maxtime;    /* maximum total time */
  INT                nstep;      /* number of steps */
  DOUBLE             dt;         /* stepsize */
  DOUBLE             time;       /* actual time */
  INT                step;       /* actual step */
  /* integrator thingies */
  DOUBLE             gamma;      /* 'theta' in one-step-theta */
  /* NRI thingies */
  INT                iter;       /* number of active iteration */
  INT                maxiter;    /* maximum number of iterations */
  DOUBLE             toldisp;    /* displacement tolerance */
  /* output thingies */
  INT                out_res_ev; /* print result every */
  /* field thingies */
  DOUBLE             initmpr;    /* (flat) initial domain temperature */
} THERM_DYNAMIC;


/*----------------------------------------------------------------------*
 | general TSI dynamic control variables                  bborn 03/06   |
 *----------------------------------------------------------------------*/
typedef struct _TSI_DYNAMIC
{
  /* core algorithm choice */
  enum
  {
    tsi_none,                       /* default */
    tsi_full_fehlbg,                /* fully coupled analysis,
                                     * both fields vary in time
                                     * and are solved with Fehlberg4/5 */
    tsi_therm_ost_struct_genalp,    /* An instationary thermal field is solved
                                     * at each time step by prescribed
                                     * BCs by one-step-theta. Subsequently, the
                                     * structural dyamics is solved by
                                     * generalised-alpha introducing
                                     * the current temperature. */
    tsi_therm_presc_struct_genalp,  /* A steady thermal field is solved
                                     * at each time step by prescribed
                                     * heat load. Subsequently, the
                                     * structural dyamics is solved by
                                     * generalised-alpha introducing
                                     * the current temperature. */
    tsi_therm_stat_struct_genalp,   /* thermal field is initially
                                     * solved, then only the structure
                                     * is time-integrated with
                                     * generalised-alpha scheme */
    tsi_therm_stat_struct_cendif,   /* thermal field is initially
                                     * solved, then only the structure
                                     * is time-integrated with
                                     * central-differences scheme */
    tsi_therm_stat_struct_fehlbg,   /* thermal field is initially
                                     * solved, then only the structure
                                     * is time-integrated with
                                     * Fehlberg4 scheme */
    tsi_therm_pred_struct_dyn       /* thermal field is predefined
                                     * initially, ie it is not solved,
                                     * later on only the structural field
                                     * is integrated in time */
  }                kind;            /* kind of analysis */

  /* time constants */
  DOUBLE           maxtime;         /* total (final) time */
  INT              nstep;           /* number of steps */
  DOUBLE           dt;              /* time step size */

  /* time variables */
  DOUBLE           time;            /* time */
  INT              step;            /* current time step */

  /* iteration constants */
  INT              maxiter;         /* maximum number of iterations */
  DOUBLE           entol;           /* tolerance for energy check
                                       over fields */

  /* integrator constants */
  DOUBLE           th_gamma;        /* 'theta' of one-step-theta */
  DOUBLE           st_beta;         /* 0.0<'beta'<=0.5 in
                                     * generalised-alpha/Newmark's method */
  DOUBLE           st_gamma;        /* 0.0<'gamma'<=1.0 in
                                     * generalised-alpha/Newmark's method */
  DOUBLE           st_alpha_m;      /* 0.0<'alpha_m'<1.0 in
                                     * generalised-alpha method */
  DOUBLE           st_alpha_f;      /* 0.0<'alpha_f'<1.0 in
                                     * generalised-alpha method */

  /* integrator variables */
  INT              actstg;          /* RK current stage */

  /* output constants */
  INT              out_std_ev;      /* print to STDOUT every */
  INT              out_res_ev;      /* print result every */

  /* field settings */
  DOUBLE           th_initmpr;      /* (flat) initial domain temperature */
} TSI_DYNAMIC;
