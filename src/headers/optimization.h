/*!---------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Andreas Lipka
            lipka@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/lipka/
            0771 - 685-6575
</pre>

---------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | type definitions optimization                         a.lipka 5/01   |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 |  optimization data                                    a.lipka 5/01   |
 *----------------------------------------------------------------------*/
typedef struct _OPTI
{

     enum
     {                             /* type of structural optimization   */
        ot_none,
        ot_shape_optimization   ,
        ot_topology_optimization
     }                    opttype;

     INT                  numiter;  /* total number of iteration steps  */
     INT                  optstep;  /* current iteration step           */
     INT                  graph_out;/* graphical output every ? step AS */

     enum
     {                              /* type of optimization strategy    */
        os_none,
        os_fsd,
        os_nlp
     }                    strategy;

     enum
     {                              /* objective                        */
        oj_none,
        oj_structural_weight,
        oj_structural_volume,
        oj_strain_energy,
        oj_stress_function,
        oj_frequency ,
        oj_buckling,
        oj_ductility
     }                   objective;


     union
     {                              /* type of optimization strategy    */
       struct _OSFSD *fsd;
       struct _OSNLP *nlp;
     }                       strat;

     INT  numvar;                /* total number struct. opt. variables */
     INT  nvaind;                /* number of independent opt variables */
     INT  numeqc;                /* number equality     constraints     */
     INT  numiqc;                /* number inequality   constraints     */
     INT  numres;                /* number restriction data             */
     INT  numlin;                /* number linking     data             */
     struct _OBTV            *ovar;           /* optimization variables */
     struct _OEQC            *oeqc;           /* equality   constraints */
     struct _OIQC            *oiqc;           /* inequality constraints */
     struct _ORES            *ores;           /* variable restricitons  */
     struct _OLIN            *olin;           /* variable linking       */

     struct _OEIG      *oeig;/* struct contains eigenvalue analysis data */
     DOUBLE totvol;
     DOUBLE totwgt;
     DOUBLE totmas;
     DOUBLE nlndat;
     /*---- smoothing of objectives or densities or ... ---*/
     /* input file should contain:
     OPT_SMO       ON          : smooth gradients or other values
     SMO_TYPE      gradient
     SMO_ERAD      2.5         : radius
     SMO_EXPO      2.5         : exponent */
     /*----------------------------------------------------*/
     enum /* switch smoothing on/off */
     {
        sm_off,
        sm_on
     }                    optsmooth;
     enum /* what to smooth: */
     {
        sm_grad,
        sm_dens
     }                    smoothtype;
     /* modified sigmund smoothing function (diss: maute...)*/
     DOUBLE smoothrad; /* radius   */
     DOUBLE smoothexp; /* exponent */

} OPTI;
/*----------------------------------------------------------------------*
 | evaluation of function gradients                    a.lipka 5/01     |
 *----------------------------------------------------------------------*/
typedef enum _EVGRAD
{
   implicit,
   explicit,
   variational
}  EVGRAD;
/*----------------------------------------------------------------------*
 | linking rules of variable design coordinates        a.lipka 5/01     |
 *----------------------------------------------------------------------*/
typedef enum _LINKR
     {
        no_olin,         /* type of objects to be linked with           */
        dncl,            /* design nodal coordinates                    */
        dntl,            /* design nodal thickness                      */
        matl             /* linking of diffenrent materials             */

     } LINKR;
/*----------------------------------------------------------------------*
 | fully stressed design stragedy                      a.lipka 5/01     |
 *----------------------------------------------------------------------*/
typedef struct _OSFSD
{
     INT     numiter;   /* number of iteration steps with fsd           */
     INT      numvar;   /* number number of independent opt. variables  */
     INT  numvar_lin;   /* no. of elements with density as opt. var. when linking. AS */
     DOUBLE   grdinc;   /* forward difference step increment            */

     EVGRAD    fgrad;   /* evaluation of function gradients             */

     DOUBLE    acc  ;   /*  erforderliche genauigkeit                   */
     DOUBLE    alpha;   /*  genauigkeit fuer die gleichheitsrestriktion */
     DOUBLE    beta ;   /*  schreitweitenfaktor                         */
     DOUBLE    delta;   /*  schrittweitenbegrenzung                     */
     DOUBLE    gamma;   /*  gleichheitsrestriktionswert                 */

     DOUBLE    *grdobj;
     DOUBLE    *grdcon;
     DOUBLE    *var   ;
     DOUBLE    *resu  ;
     DOUBLE    *resl  ;
     DOUBLE    *grdobj_lin;     /* resized grdobj vector for linking AS*/
     DOUBLE    *grdcon_lin;     /* resized grdcon vector for linking AS*/
     DOUBLE    *var_lin   ;     /* resized var    vector for linking AS*/


} OSFSD;
/*----------------------------------------------------------------------*
 | nonlinear prog.       stragedy                      a.lipka 5/01     |
 *----------------------------------------------------------------------*/
typedef struct _OSNLP
{
     EVGRAD    fgrad;   /* evaluation of function gradients             */
     DOUBLE   grdinc;   /* forward difference step increment            */

     enum
     {                  /* type of scaling the variables                */
        no,
        initial,
        hessian
     }                         sclvar;

     DOUBLE    acc ;   /*  accuracy in program nlpql                    */
     DOUBLE sclfun ;    /* scaling factor for objective                 */
     DOUBLE sclcon ;    /* scaling factors for constraints              */
     DOUBLE scbou  ;    /* variable for autom. scaling of objective     */
     INT    numiter;    /* number of iteration steps with nlp           */
     INT    nlpldl ;    /* flag for control of nonlinear programming    */
     INT    nlluen ;    /* lueunberger selfscaling bfgs parameter       */
     INT    maxfun ;    /*  max. n. o. function calls in line search    */
     INT    lise   ;    /*  flag for control of nonlinear programming   */
     DOUBLE amue   ;    /*  factor for line search adjustment           */
     INT    merit  ;    /*  aug. lagrangian / l1-penalty merit function */



     INT numcon; /* max. number struct. opt. constraints      */
     INT numeqc; /* max. number equality     constraints      */
     INT iscvar; /* type of scaling the variables             */
     INT igrad ; /* define kind of gradient evaluation        */
     INT isqp  ; /* flag for control of nonlinear programming */
     INT mode  ; /* mode = 0 : normal execution.              */
     DOUBLE time; /* cpu time used for optimization with nlp  */
     INT nprint; /* write nlp process data to file            */
     DOUBLE *var;  /* initial guess for the optimal solution, on return: the last computed iterate */
     DOUBLE *resl; /* the upper bounds of the variables */
     DOUBLE *resu; /* the lower bounds of the variables */
     DOUBLE *vscale; /* scaling factors  for variables   */
     DOUBLE *rconst; /* equilibrium equality constraint array */
     DOUBLE *cscale; /* scaling factors  for constraints */
     DOUBLE *rumult; /* vector of lagrange multipliers  */
     INT    *active; /* active constraints (active(j) = .true.) */
     DOUBLE *df;     /* gradients of the objective function */
     DOUBLE *dg;     /* gradients of the active constraints */
     DOUBLE *rchess; /* approximation of the hessian matrix */
     DOUBLE *rdhess; /* hessian matrix ?*/
     DOUBLE *wa;     /* wa is a real working array of length lwa */
     INT    *kwa;    /* wa is a intg working array of length lwa */
     INT    *lwa;    /* wa is a intg working array of length lwa */
     DOUBLE *f; /* the objective function value of x */
     INT    nitstep; /* current iteration step */
     /* dimensions for wa's ... */
     INT n;
     INT nr;
     INT mnn2;
     INT mmax;   /* row dimension of dg >=1                    */
     INT nmax;   /* row dimension of c  >=1                    */
     INT nc1;    /* = nmax; first  dimension of c              */
     INT nc2;    /* = nmax; second dimension of c              */
     INT nlwa ;
     INT lkwa ;
     INT llwa ;
     DOUBLE *scg ;/* scaling factors for constraints */
     INT ixmax;
     INT ixmin;







} OSNLP;
/*----------------------------------------------------------------------*
 | optimization variables                              a.lipka 5/01     |
 *----------------------------------------------------------------------*/
typedef struct _OBTV
{
     INT objId;          /* object Id (design node, material number ... */
     INT posgl;          /* position in global variable vector          */
     /* type of object with variable attributes */
     enum
     {
        no_ovatt,
        dheight,                  /* variable design thickness          */
        dcoorx,                   /* variable design coordinates        */
        dcoory,
        dcoorz,
        eleofdesofmat,     /* variable material data of design elements */
        eleofmat,                 /* variable material data of elements */
        material                  /* variable material                  */
     }            ovatt;
     INT    numvar     ;/*  number of opt. variables                    */
     DOUBLE *var       ;/*  pointer to opt. variable                    */
     DOUBLE scva       ;/*  scaling factor                              */
     DOUBLE bupper     ;/*  upper limitation of opt. variable           */
     DOUBLE blower     ;/*  lower limitation of opt. variable           */

     INT    linked     ;/*  prescribed linking                          */
     DOUBLE  myweight   ;/*  weighting factor for linking               */
     DOUBLE  neweight   ;/* weighting factor for linking                */
     LINKR olin_type;          /* type of objects to be linked with              */
} OBTV;
/*----------------------------------------------------------------------*
 | equality constraints                                a.lipka 5/01     |
 *----------------------------------------------------------------------*/
typedef struct _OEQC
{
     INT objId;          /* object Id (design node, material number ... */

     enum
     {
        no_oeqco,                 /* object with   equality constraint  */
        eqdface,                  /* design face   equality constraint  */
        eqdvolu                   /* design volume equality constraint  */
     }                         oeqc_object;
     enum
     {
        no_oeqct,                  /* type of equality constraint       */
        volume,                    /* volume  equality constraint       */
        mass                       /* mass    equality constraint       */
     }                         oeqc_type;
     DOUBLE val        ;           /*  value  of equality constraint    */
     DOUBLE scl        ;           /*  scaling factor of constraint     */
} OEQC;
/*----------------------------------------------------------------------*
 | inequality constraints                              a.lipka 5/01     |
 *----------------------------------------------------------------------*/
typedef struct _OIQC
{
     INT objIds[2];   /* object Ids (design nodes, material numbers ... */

     enum
     {
        no_oiqco,                 /* object with inequality constraint  */
        iednode,                  /* design node inequality constraint  */
        iedface,                  /* design face inequality constraint  */
        iefeele                   /* finite ele. inequality constraint  */
     }                         oiqc_object;
     enum
     {
        no_oiqct,                  /* type of inequality constraint     */
        coorx,                     /* coord. constr. in z direction     */
        coory,                     /* coord. constr. in z direction     */
        coorz,                     /* coord. constr. in z direction     */
        sve                        /* v.mises element stress constraint */
     }                         oiqc_type;
     enum
     {
        no_opie,                   /* type of constraint operation      */
        lowt,                      /* >     */
        lart                       /* <     */
     }                         oiqc_operator;

     DOUBLE *val        ;         /* values  of inequality constraint   */
} OIQC;
/*----------------------------------------------------------------------*
 | variable restricitons                               a.lipka 5/01     |
 |                                                                      |
 | DNODE  1  ZL  -10  ZU 50.0   ZS 1.0 :                                |
 |                                                                      |
 | ores->ores_type = dnulr                                              |
 | ores->resuval   =  50.0                                              |
 | ores->reslval   = -10.0                                              |
 | ores->ressval   =   1.0                                              |
 *----------------------------------------------------------------------*/
typedef struct _ORES
{
     INT objId;          /* object Id (design node, material number ... */

     enum
     {
        no_orest,   /* type of restriction                              */
        dnthick,    /* design node thickness upper/lower restriction    */
        dncooz,     /* design coordinates: upper/lower restriction zdir */
        maulr       /*           material: upper/lower restriction      */

     } ores_type;

     DOUBLE resuval;                               /*  upper value */
     DOUBLE reslval;                               /*  lower value */
     DOUBLE ressval;                               /*  scale value */
} ORES;
/*----------------------------------------------------------------------*
 | variable linking                                    a.lipka 5/01     |
 |                                                                      |
 | DNODE  9  LIKE  3  :                                                 |
 |                                                                      |
 | olin->olin_type = dncl                                               |
 | olin->resuval   =  50.0                                              |
 | olin->reslval   = -10.0                                              |
 | olin->ressval   =   1.0                                              |
 *----------------------------------------------------------------------*/
typedef struct _OLIN
{
     INT objIds[2];   /* object Ids (design nodes, material numbers ... */
     DOUBLE  myweight   ;/* weighting factor for linking                   */
     DOUBLE  neweight   ;/* weighting factor for linking                   */

     LINKR olin_type;          /* type of objects to be linked with              */
} OLIN;
/*----------------------------------------------------------------------*
 | struct for eigenvalue analysis                      a.lipka 5/01     |
 |                                                                      |
 *----------------------------------------------------------------------*/
typedef struct _OEIG
{
     INT      numeigv; /* number of eigenvalues   */
     DOUBLE  eigv[10]; /* vector with eigenvalues */
     DOUBLE  eigf[10]; /* vector with frequencies */
     DOUBLE  eigs[10]; /* scaling vector          */
     DOUBLE  *scasens; /* vector with scaling factors for sens. ana.*/
     /*--- KREISSELMEIER-STEINHAUSER  ---*/
     DOUBLE  rhoks;
} OEIG;

typedef enum _OPT_GR_OUT
{
   gr_init,     /* initialize                       */
   gr_mesh,     /* write meshdata                   */
   gr_dens,     /* write density values of elements */
   gr_disp      /* write displacements */
} OPT_GR_OUT;
