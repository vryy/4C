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

     int                  numiter;  /* total number of iteration steps  */
                               
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
                                         
     int  numvar;                /* total number struct. opt. variables */
     int  nvaind;                /* number of independent opt variables */
     int  numeqc;                /* number equality     constraints     */
     int  numiqc;                /* number inequality   constraints     */
     int  numres;                /* number restriction data             */
     int  numlin;                /* number linking     data             */
     struct _OBTV            *ovar;           /* optimization variables */
     struct _OEQC            *oeqc;           /* equality   constraints */
     struct _OIQC            *oiqc;           /* inequality constraints */
     struct _ORES            *ores;           /* variable restricitons  */
     struct _OLIN            *olin;           /* variable linking       */
     
     double totvol;
     double totwgt;
     double totmas;
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
     double smoothrad; /* radius   */
     double smoothexp; /* exponent */

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
     int     numiter;   /* number of iteration steps with fsd           */
     int      numvar;   /* number number of independent opt. variables  */
     double   grdinc;   /* forward difference step increment            */
     
     EVGRAD    fgrad;   /* evaluation of function gradients             */

     double    acc  ;   /*  erforderliche genauigkeit                   */
     double    alpha;   /*  genauigkeit fuer die gleichheitsrestriktion */
     double    beta ;   /*  schreitweitenfaktor                         */
     double    delta;   /*  schrittweitenbegrenzung                     */
     double    gamma;   /*  gleichheitsrestriktionswert                 */

     double    *grdobj;
     double    *grdcon;
     double    *var   ;
     double    *resu  ;
     double    *resl  ;


} OSFSD;
/*----------------------------------------------------------------------*
 | nonlinear prog.       stragedy                      a.lipka 5/01     |
 *----------------------------------------------------------------------*/
typedef struct _OSNLP
{
     EVGRAD    fgrad;   /* evaluation of function gradients             */
     double   grdinc;   /* forward difference step increment            */
     
     enum 
     {                  /* type of scaling the variables                */
        no,
        initial,
        hessian
     }                         sclvar;
     
     double    acc ;   /*  accuracy in program nlpql                    */
     double sclfun ;    /* scaling factor for objective                 */
     double sclcon ;    /* scaling factors for constraints              */
     double scbou  ;    /* variable for autom. scaling of objective     */
     int    numiter;    /* number of iteration steps with nlp           */
     int    nlpldl ;    /* flag for control of nonlinear programming    */
     int    nlluen ;    /* lueunberger selfscaling bfgs parameter       */
     int    maxfun ;    /*  max. n. o. function calls in line search    */
     int    lise   ;    /*  flag for control of nonlinear programming   */
     double amue   ;    /*  factor for line search adjustment           */
     int    merit  ;    /*  aug. lagrangian / l1-penalty merit function */



     int numcon; /* max. number struct. opt. constraints      */
     int numeqc; /* max. number equality     constraints      */
     int iscvar; /* type of scaling the variables             */
     int igrad ; /* define kind of gradient evaluation        */
     int isqp  ; /* flag for control of nonlinear programming */
     int mode  ; /* mode = 0 : normal execution.              */
     double time; /* cpu time used for optimization with nlp  */   
     int nprint; /* write nlp process data to file            */
     double *var;  /* initial guess for the optimal solution, on return: the last computed iterate */
     double *resl; /* the upper bounds of the variables */
     double *resu; /* the lower bounds of the variables */
     double *vscale; /* scaling factors  for variables   */
     double *rconst; /* equilibrium equality constraint array */
     double *cscale; /* scaling factors  for constraints */
     double *rumult; /* vector of lagrange multipliers  */
     int    *active; /* active constraints (active(j) = .true.) */
     double *df;     /* gradients of the objective function */
     double *dg;     /* gradients of the active constraints */
     double *rchess; /* approximation of the hessian matrix */
     double *rdhess; /* hessian matrix ?*/
     double *wa;     /* wa is a real working array of length lwa */
     int    *kwa;    /* wa is a intg working array of length lwa */
     int    *lwa;    /* wa is a intg working array of length lwa */
     double *f; /* the objective function value of x */   
     int    nitstep; /* current iteration step */
     /* dimensions for wa's ... */
     int n;
     int nr;
     int mnn2;
     int mmax;   /* row dimension of dg >=1                    */
     int nmax;   /* row dimension of c  >=1                    */
     int nc1;    /* = nmax; first  dimension of c              */
     int nc2;    /* = nmax; second dimension of c              */
     int nlwa ;
     int lkwa ;
     int llwa ;
     double *scg ;/* scaling factors for constraints */ 
     int ixmax;
     int ixmin;
     






} OSNLP;
/*----------------------------------------------------------------------*
 | optimization variables                              a.lipka 5/01     |
 *----------------------------------------------------------------------*/
typedef struct _OBTV
{
     int objId;          /* object Id (design node, material number ... */
     int posgl;          /* position in global variable vector          */ 
     /* type of object with variable attributes */
     enum 
     {                  
        no_ovatt,
        dheight,                  /* variable design thickness          */
        dcoorx,                   /* variable design coordinates        */
        dcoory,
        dcoorz,
        eleofmat,                 /* variable material data of elements */
        material                  /* variable material                  */
     }            ovatt;
     int    numvar     ;/*  number of opt. variables                    */
     double *var       ;/*  pointer to opt. variable                    */
     double scva       ;/*  scaling factor                              */
     double bupper     ;/*  upper limitation of opt. variable           */
     double blower     ;/*  lower limitation of opt. variable           */
     
     int    linked     ;/*  prescribed linking                          */
     LINKR olin_type;          /* type of objects to be linked with              */
} OBTV;
/*----------------------------------------------------------------------*
 | equality constraints                                a.lipka 5/01     |
 *----------------------------------------------------------------------*/
typedef struct _OEQC
{
     int objId;          /* object Id (design node, material number ... */
    
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
     double val        ;           /*  value  of equality constraint    */
     double scl        ;           /*  scaling factor of constraint     */
} OEQC;
/*----------------------------------------------------------------------*
 | inequality constraints                              a.lipka 5/01     |
 *----------------------------------------------------------------------*/
typedef struct _OIQC
{
     int objIds[2];   /* object Ids (design nodes, material numbers ... */
     
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

     double *val        ;         /* values  of inequality constraint   */
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
     int objId;          /* object Id (design node, material number ... */
     
     enum 
     {                  
        no_orest,   /* type of restriction                              */
        dnthick,    /* design node thickness upper/lower restriction    */
        dncooz,     /* design coordinates: upper/lower restriction zdir */
        maulr       /*           material: upper/lower restriction      */
     
     } ores_type;

     double resuval;                               /*  upper value */
     double reslval;                               /*  lower value */
     double ressval;                               /*  scale value */
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
     int objIds[2];   /* object Ids (design nodes, material numbers ... */
     
     LINKR olin_type;          /* type of objects to be linked with              */
} OLIN;

typedef enum _OPT_GR_OUT
{
   gr_init,     /* initialize                       */
   gr_mesh,     /* write meshdata                   */
   gr_dens,     /* write density values of elements */
   gr_disp      /* write displacements */
} OPT_GR_OUT;
