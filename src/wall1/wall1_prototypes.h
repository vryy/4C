/*!----------------------------------------------------------------------
\file
\brief contains all prototypes for wall1 element

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0771 - 685-6122
</pre>

*----------------------------------------------------------------------*/
#ifdef D_WALL1

/*! 
\addtogroup WALL1
*//*! @{ (documentation module open)*/

/*-----------------------------------------------------------------------*
|  w1_bop.c                                                  al 9/01     |
|  calculate operator matrix at point r,s                                |
*-----------------------------------------------------------------------*/
void w1_bop(DOUBLE    **bop,   /* derivative operator                   */
            DOUBLE    **deriv, /* shape functions derivatives           */
            DOUBLE    **xjm,   /* jacobian matrix                       */
            DOUBLE      det,   /* determinant of jacobian               */
            INT         iel);  /* number of nodes at actual element     */
/*-----------------------------------------------------------------------*
|  w1_cal_deriv.c                                            al 9/01     |
|  compute displacement derivatives                                      |
*-----------------------------------------------------------------------*/
void w1_disd(ELEMENT   *ele,  /* actual element                         */
             DOUBLE   **bop,  /* derivative operator                    */
             DOUBLE    *gop,  /* additional derivative operator         */
             DOUBLE    *alpha,/*internal dof                            */
             WALL_TYPE wtype, /* plane stress/strain...                 */
             DOUBLE    *disd);/* displacement derivatives               */
/*-----------------------------------------------------------------------*
|  w1_cal_eps.c                                              al 9/01     |
|  evaluates linear/nonlinear strains from                               |
|  displacement derivatives                                              |
*-----------------------------------------------------------------------*/
void w1_eps( DOUBLE   *disd, /* displacement derivatives                */
             WALL_TYPE wtype,/* plane stress/strain...                  */
             DOUBLE   *eps); /* strains                                 */
/*-----------------------------------------------------------------------*
|  w1_cal_stress.c                                           al 9/01     |
|  evaluate element stresses                                             |
|  2-D isoparametric degenerated element                                 |
*-----------------------------------------------------------------------*/
void w1_cal_stress(ELEMENT   *ele, /* actual element                    */
                   W1_DATA   *data,/* wall1 data                        */
                   MATERIAL  *mat, /* actual material                   */
                   INT        kstep,  /* number of step      fh 06/02   */
                   INT        init);  /* initialize this function       */
/*-----------------------------------------------------------------------*
|  w1_cal_stress.c                                           al 9/01     |
|  principal stresses and directions                                     |
|  2D isoparametric degenerated element                                  |
| ---------------------------------------------------------------------- |
|       angle of principal stresses:                                     |
|                                                2*TAU XY                |
|                            TAN (2*ALFA) = -------------------          |
|                                           (SIGMA X - SIGMA Y)          |
*-----------------------------------------------------------------------*/
void w1_mami(DOUBLE *stress, /* vector of stresses                      */
             DOUBLE *fps,    /* first  principal stress                 */ 
             DOUBLE *sps,    /* second principal stress                 */ 
             DOUBLE *aps);   /* angle of principal direction            */ 
/*-----------------------------------------------------------------------*
|  w1_cal_stress.c                                           al 9/01     |
|       returns R/S coordinates of gauss integration points              |
|       2D isoparametric degenerated element                             |
| ---------------------------------------------------------------------- |
|       xr489 ---> R,S values for rectangles with 4,8,9 nodes            |
|       xr16  ---> R,S values for rectangles with 12,16 nodes            |
|       xt36  ---> R,S values for triangle   with 3,6   nodes            |
|       xt10  ---> R,S values for triangle   with 10    nodes            |
*-----------------------------------------------------------------------*/
DOUBLE w1rsn (INT node, /* number of actual integration point           */
              INT  irs, /* r/s identifier                               */
              INT  iel);/* number of nodes at actual element            */
/*-----------------------------------------------------------------------*
|  w1_cal_stress.c                                           al 9/01     |
|       extrapolation form gauss points for rectangles                   |
*-----------------------------------------------------------------------*/
void w1recs(DOUBLE *funval,/* function value                            */
            DOUBLE  r,     /* r/s coordinates                           */
            DOUBLE  s,     /* r/s coordinates                           */
            DOUBLE *fval,  /* function values at gauss points           */
            DOUBLE *fpar,  /* function parameters                       */
            INT     igauss,/* number of gauss points                    */
            INT     icode);/* ==1 initialize function parameters        */
                           /* > 1            function evaluation        */      
/*-----------------------------------------------------------------------*
|  w1_call_mat.c                                             al 9/01     |
|  select proper material law                                            |
*-----------------------------------------------------------------------*/
void w1_call_mat(ELEMENT   *ele, /* actual element                      */
                 MATERIAL  *mat, /* actual material                     */
                 WALL_TYPE wtype,/* plane stress/strain...              */  
                 DOUBLE **bop,   /* derivative operator                 */
                 DOUBLE  *gop,   /* additional derivative operator      */
                 DOUBLE  *alpha, /* internal dof                        */
                 INT ip,         /* integration point Id                */
                 DOUBLE *stress, /* vector of stresses                  */
                 DOUBLE **d,     /* constitutive matrix                 */
                 INT istore,     /* controls storing of new stresses    */
                 INT newval);    /* controls evaluation of new stresses */   
/*-----------------------------------------------------------------------*
|  w1_call_stiff.c                                           al 9/01     |
|  usual stiffness matrix total lagrangian formulation                   |
*-----------------------------------------------------------------------*/
void w1_keku(DOUBLE **s,    /* element stiffness matrix                 */
             DOUBLE **bs,   /* derivative operator                      */
             DOUBLE **d,    /* constitutive matrix                      */
             DOUBLE   fac,  /* integration factor                       */
             INT      nd,   /* total number degrees of freedom of ele.  */
             INT      neps);/* actual number of strain components       */
/*-----------------------------------------------------------------------*
|  w1_funcderiv.c                                            al 9/01     |
|  shape functions and derivatives                                       |
*-----------------------------------------------------------------------*/
void w1_funct_deriv(DOUBLE     *funct,  /* shape function values        */
                    DOUBLE    **deriv,  /* shape function derivatives   */
                    DOUBLE      r,      /* r/s coordinates              */
                    DOUBLE      s,      /* r/s coordinates              */
                    DIS_TYP     typ,    /* quad4, quad8, quad9 ...      */
                    INT         option);/* index for function evaluation*/ 
void w1_degfuncderiv(DOUBLE     *funct, 
                      DOUBLE    **deriv, 
                      DOUBLE      r, 
                      DIS_TYP     typ,
                      INT         option);
/*-----------------------------------------------------------------------*
|  w1_init.c                                                 al 9/01     |
|  initialize the element                                                |
*-----------------------------------------------------------------------*/
void w1init(PARTITION *actpart,
            MATERIAL    *mat );/* actual material                       */
/*-----------------------------------------------------------------------*
|  w1_intg.c                                                 al 9/01     |
|  integration points                                                    |
 ----------------------------------------------------------------------- |
|      ngp[0] --> number of integration points r-direction               |
|                 for triangular elements :                              |
|                 number of over all integration points                  |
|      ngp[1] --> number of integration points s-direction               |
|                 for triangular elements :                              |
|                 parameter for alternative integration                  |
|                 0 = standard    1 = gauss-radau                        |
|      ntyp   --> 1 = triangular, 0 = rectangular  element               |
|      xgrr   --> gauss sampling points        r-direction               |
|      xgss   --> gauss sampling points        s-direction               |
|      wgtr   --> weighting factors            r-direction               |
|      wgts   --> weighting factors            s-direction               |
*-----------------------------------------------------------------------*/
void w1intg(ELEMENT   *ele,    /* actual element                        */
            W1_DATA   *data,   /* wall1 data                            */
            INT        option);/* index for function evaluation         */
/*-----------------------------------------------------------------------*
|  w1_jaco.c                                                 al 9/01     |
|  calculate operator matrix at point r,s                                |
*-----------------------------------------------------------------------*/
void w1_jaco(DOUBLE    **deriv,/* shape function derivatives            */
             DOUBLE    **xjm,  /* jacobian matrix                       */
             DOUBLE     *det,  /* determinant of jacobian               */
             ELEMENT    *ele,  /* actual element                        */
             INT         iel); /* number of nodes at actual element     */
/*-----------------------------------------------------------------------*
|  w1_mat_linel.c                                            al 9/01     |
|  constitutive matrix - linear elastic - 2D                             |
|  plane stress, plane strain, rotational symmetry                       |
*-----------------------------------------------------------------------*/
void w1_mat_linel(DOUBLE ym,      /* young's modulus                    */     
                  DOUBLE pv,      /* poisson's ratio                    */     
                  WALL_TYPE wtype,/* plane stress/strain...             */   
                  DOUBLE **d);    /* constitutive matrix                */   
void w1_mat_stvpor(MATERIAL  *mat,
                   DOUBLE *matdata,
                   WALL_TYPE wtype,/* plane stress/strain...            */
                   DOUBLE **d);    /* constitutive matrix               */
/*-----------------------------------------------------------------------*
|  w1_mat_plast_dp.c                                         al 9/01     |
|  constitutive matrix - forces - Drucker Prager - 2D                    |
|  plane stress, plane strain, rotational symmetry                       |
*-----------------------------------------------------------------------*/
void w1_mat_plast_dp(DOUBLE ym,      /* young's modulus                 */     
                     DOUBLE pv,      /* poisson's ratio                 */     
                     DOUBLE sigy,    /* yield stresse                   */    
                     DOUBLE hard,    /* hardening modulus               */    
                     DOUBLE phi,     /* angle - drucker prager          */    
                     ELEMENT   *ele, /* actual element                  */ 
                     WALL_TYPE wtype,/* plane stress/strain...          */      
                     DOUBLE **bop,   /* derivative operator             */  
                     DOUBLE  *gop,   /* add. derivative operator        */  
                     DOUBLE  *alpha, /* internal dof                    */  
                     INT ip,         /* integration point Id            */         
                     DOUBLE *stress, /* vector of stresses              */
                     DOUBLE **d,     /* constitutive matrix             */    
                     INT istore,     /* controls storing of new stresses*/    
                     INT newval);    /* controls evaluation of stresses */ 
/*-----------------------------------------------------------------------*
|  w1_mat_plast_epc.c                                        al 9/01     |
|  const. matrix - forces - elastoplastic concrete - 2D                  |
|  plane stress, plane strain                                            |
*-----------------------------------------------------------------------*/
void w1_mat_plast_epc(DOUBLE emod      , /* young's modulus             */
                      DOUBLE vnu       , /* poisson's ratio             */
                      DOUBLE ftm       , /* tensile strenght            */                    
                      DOUBLE fcm       , /* compressive strenght        */        
                      DOUBLE gt        , /* tensile fracture energy     */        
                      DOUBLE gc        , /* compression fracture energy */        
                      DOUBLE gamma1    , /* fitting factor yield funct 1*/     
                      DOUBLE gamma2    , /* fitting factor yield funct 2*/     
                      INT    nstiff    , /* flag for tension stiffening */
                      INT    maxreb    , /* max. number of rebars       */             
                      DOUBLE *reb_area  ,/* rebar area                  */
                      DOUBLE *reb_ang   ,/* rebar angle                 */
                      DOUBLE *reb_so    ,/* minimum bond length         */
                      DOUBLE *reb_ds    ,/* diameter                    */
                      DOUBLE *reb_rgamma,/* 4=deformed bars 2=plane bar */
                      ELEMENT   *ele,    /* actual element              */
                      WALL_TYPE wtype,   /* plane stress/strain...      */          
                      DOUBLE **bop,      /* derivative operator         */
                      DOUBLE  *gop,      /* add derivative operator     */
                      DOUBLE  *alpha,    /* internal dof                */
                      INT ip,            /* integration point Id        */
                      DOUBLE *stressc,   /* vector of stresses          */
                      DOUBLE **d,        /* constitutive matrix         */
                      INT istore,        /* controls storing of stresses*/
                      INT newval);       /* controls eval. of stresses  */  
/*-----------------------------------------------------------------------|
| w1_mat_plast_epc_serv.c                                     al 9/01    |
|      initialize d-components d[0][3], d[1][3], d[2][3], d[3][3]        |
|      on elastic                                                        |
|-----------------------------------------------------------------------*/
void w1iwadi (DOUBLE ym,   /* young's modulus                           */                           
              DOUBLE pv,   /* poisson's ratio                           */ 
              DOUBLE *di); /* components of constitutive matrix         */
/*-----------------------------------------------------------------------|
| w1_mat_plast_epc_serv.c                                     al 9/01    |
|      topic : wall1 - yield criterion for plane strain                  |
|              yield - criterion:  drucker - prager                      |
|-----------------------------------------------------------------------*/
void w1yicsr (DOUBLE dev,      /* norm der deviatorischen spannungen    */
              DOUBLE hyd,      /* 1. invariante                         */
              DOUBLE sigym,    /* uniaxial yield stress                 */
              DOUBLE alpha,    /* neigungswinkel im invariantenraum     */
              DOUBLE *ft);     /* yield condition                       */ 
/*-----------------------------------------------------------------------|
| w1_mat_plast_epc_serv.c                                     al 9/01    |
|      topic : wall1 - yield criterion (cap - region)                    |
|              yield - criterion:  cap - model                           |
|-----------------------------------------------------------------------*/
void w1yiccap (DOUBLE dev,   /* norm of deviatoric predictor stresses   */
               DOUBLE hyd,   /* 1. invariant  of the predictor stresses */
               DOUBLE *sigym,/* yield stresses                          */
               DOUBLE alpha, /* factor for the first invariants         */
               DOUBLE *ft);  /* yield function                          */
/*-----------------------------------------------------------------------|
| w1_mat_plast_epc_serv.c                                     al 9/01    |
|    topic: calculate the gradients of the elastic predictor             |              
|           deviatoric/hydrostatic stresses                              |
|-----------------------------------------------------------------------*/
void w1preds(DOUBLE *sigma, /* elastic predictor projected              */
             DOUBLE *alpha, /* neigungswinkel im invariantenraum        */
             DOUBLE *devsig,/* deviatoric stresses                      */
             DOUBLE *dev,   /* norm der deviatorischen spannungen       */
             DOUBLE *dn,    /* gradient in deviatoric direction         */
             DOUBLE *dcom,  /* gradient in hdrostatic direction         */
             DOUBLE *grad); /* total gradient                           */
/*-----------------------------------------------------------------------|
| w1_mat_plast_epc_serv.c                                     al 9/01    |
|    topic: calculate the deviatoric/hydrostatic stress components       |
|           of the elastic predictor                                     |
|-----------------------------------------------------------------------*/
void w1pres (DOUBLE *sigma ,    /*  elastic predictor projected         */
             DOUBLE *devsig,    /*  deviatoric stresses                 */
             DOUBLE *sm    ,    /*  hydrrostatic stresses               */
             DOUBLE *dev   ,    /*  norm der deviatorischen spannungen  */
             DOUBLE *hyd);      /*  1. invariante                       */
/*-----------------------------------------------------------------------*
|  w1_mat_plast_mises.c                                      al 9/01     |
|  constitutive matrix - forces - linear elastic- von Mises - 2D a       |
|  plane stress, plane strain, rotational symmetry                       |
*-----------------------------------------------------------------------*/
void w1_mat_plast_mises(DOUBLE ym,      /* young's modulus              */
                        DOUBLE pv,      /* poisson's ratio              */
                        DOUBLE sigy,    /* yield stresse                */
                        DOUBLE hard,    /* hardening modulus            */       
                        DOUBLE gf,      /* fracture energy              */
                        DOUBLE betah,
                        ELEMENT   *ele, /* actual element               */
                        WALL_TYPE wtype,/* plane stress/strain...       */         
                        DOUBLE **bop,   /* derivative operator          */
                        DOUBLE  *gop,   /* add. derivative operator     */
                        DOUBLE  *alpha, /* internal dof                 */
                        INT ip,         /* integration point Id         */
                        DOUBLE *stress, /* vector of stresses           */
                        DOUBLE **d,     /* constitutive matrix          */
                        INT istore,     /* controls storing of stresses */
                        INT newval);    /* controls eval. of stresses   */ 
/*-----------------------------------------------------------------------*
|  w1_mat_plast_serv.c                                       al 9/01     |
|                                                                        |
|    topic : wall1 - yield criterion plane_stress                        |
|            -   -   -- -  --                                            |
|            yield - criterion:  phi > 0 drucker - prager                |
|                                phi = 0 mises                           |
|            for material model 3 - plasticity                           |
*-----------------------------------------------------------------------*/
void w1yilcr(DOUBLE E,    /* young's modulus                            */
             DOUBLE Eh,   /* hardening modulus                          */
             DOUBLE betah,
             DOUBLE sigy, /* yield stresse                              */
             DOUBLE epstn,/* equivalent uniaxial plastic strain         */    
             INT    isoft,/* softening                                  */
             DOUBLE dia,  /* equivalent element length                  */
             DOUBLE *tau, /* current stresses (local)                   */
             DOUBLE *ft); /* yield condition                            */     
/*-----------------------------------------------------------------------*
|  w1_mat_plast_serv.c                                       al 9/01     |
|  radial return for elements with mises material model                  |
*-----------------------------------------------------------------------*/
void w1radi(DOUBLE e,        /* young's modulus                         */
            DOUBLE eh,       /* hardening modulus                       */
            DOUBLE betah,
            DOUBLE sigy,     /* yield stresse                           */
            DOUBLE vnu,      /* poisson's ratio                         */
            DOUBLE dia,      /* equivalent element length               */
            DOUBLE *sigma,   /* elastic predicor projected onto yield   */
            DOUBLE *qn,      /* backstress vector                       */
            INT    isoft,    /* softening                               */
            DOUBLE *epstn,   /* equivalent uniaxial plastic strain      */       
            DOUBLE *dlam,    /* increment of plastic multiplier         */
            WALL_TYPE wtype);/* plane stress/strain...                  */
/*-----------------------------------------------------------------------*
|  w1_mat_plast_serv.c                                       al 9/01     |
|    forms the elasto-plastic consistent tangent material tensor         |
*-----------------------------------------------------------------------*/
void w1mapl(DOUBLE e,        /* young's modulus                         */
            DOUBLE eh,       /* hardening modulus                       */
            DOUBLE betah,
            DOUBLE sigy,     /* yield stresse                           */
            DOUBLE vnu,      /* poisson's ratio                         */
            DOUBLE dia,      /* equivalent element length               */
            DOUBLE *tau,     /* current stresses (local)                */
            INT    isoft,    /* softening                               */
            DOUBLE *epstn,   /* equivalent uniaxial plastic strain      */       
            DOUBLE *dlam,    /* increment of plastic multiplier         */
            DOUBLE **d,      /* constitutive matrix                     */
            WALL_TYPE wtype);/* plane stress/strain...                  */
/*-----------------------------------------------------------------------*
|  w1_mat_plast_serv.c                                       al 9/01     |
|    topic : wall1 - yield criterion plane_stress                        |
|            -   -   -- -  --                                            |
|            yield - criterion:  phi > 0 drucker - prager                |
|                                phi = 0 mises                           |
*-----------------------------------------------------------------------*/
void w1yilcr_dp(DOUBLE E,      /* young's modulus                       */
                DOUBLE Eh,     /* hardening modulus                     */
                DOUBLE phi,    /* angle - drucker prager                */
                DOUBLE sigy,   /* yield stresse                         */
                DOUBLE *sigym, /* uniaxial yield stress                 */
                DOUBLE epstn,  /* equivalent uniaxial plastic strain    */         
                DOUBLE *tau,   /* current stresses (local)              */
                DOUBLE *ft);   /* yield condition                       */          
/*-----------------------------------------------------------------------*
|  w1_mat_plast_serv.c                                       al 9/01     |
|    radial return for elements with mises material model                |
|                                                                        |
*-----------------------------------------------------------------------*/
void w1radi_dp(DOUBLE e,        /* young's modulus                      */
               DOUBLE eh,       /* hardening modulus                    */
               DOUBLE phi,      /* angle - drucker prager               */
               DOUBLE sigy,     /* yield stresse                        */
               DOUBLE vnu,      /* poisson's ratio                      */
               DOUBLE *sigma,   /* elastic predicor projected onto yield*/
               DOUBLE *qn,      /* backstress vector                    */
               DOUBLE *epstn,   /* equivalent uniaxial plastic strain   */          
               DOUBLE *sigym,   /* uniaxial yield stress                */ 
               DOUBLE *dlam,    /* increment of plastic multiplier      */
               WALL_TYPE wtype);/* plane stress/strain...               */ 
/*-----------------------------------------------------------------------*
|  w1_mat_plast_serv.c                                       al 9/01     |
|     calculate element diameter (equivalent length) for one element     |            
*-----------------------------------------------------------------------*/
void w1cdia(ELEMENT   *ele,    /* actual element                        */
            W1_DATA   *data,   /* wall1 data                            */
            DOUBLE    *funct_h,/* shape function values                 */
            DOUBLE   **deriv_h,/* shape function derivatives            */
            DOUBLE   **xjm_h); /* jacobian matrix                       */
/*-----------------------------------------------------------------------|
|  w1_mat_plast_serv.c                                       al 9/01     |
|     topic : wall1 - concrete material prevalues                        |
|             -   -   -                 -----                            |
|             yield - criterion:  drucker-prager / spherical cap         |
|-----------------------------------------------------------------------*/
void w1cpreva (DOUBLE *epst,     /* equivalent uniaxial plastic strain  */
               DOUBLE *sigym,    /* [0] uniaxial tension yield stress   */
                                 /* [1] yield stress "inverted cone"    */
                                 /* [2] uniaxial compr. yield stress    */
               DOUBLE *alpha,    /* factor for the first invariants     */
	       DOUBLE *hards,    /* hardening modulus                   */
               DOUBLE *e,        /* young modulus                       */
               DOUBLE *g,        /* shear modulus                       */
               DOUBLE *vnu,      /* poisson's ratio                     */
               DOUBLE *com,      /* bilk modulus                        */
               DOUBLE *fcm,      /* compressive strenght                */
               DOUBLE *gc,       /* compression fracture energy         */
               DOUBLE *ftm,      /* tensile strenght                    */
               DOUBLE *gt,       /* tensile fracture energy             */
               DOUBLE *gamma1,   /* fitting factor yield function 1     */
               DOUBLE *gamma2,   /* fitting factor yield function 2     */
               DOUBLE *dfac,     /* damage factor                       */
               DOUBLE *dia,      /* equivalent element length           */
	       DOUBLE *acrs,     /* average crack spacing               */
               DOUBLE *cappaet,  /* max. elastic tension strain         */
               DOUBLE *cappaut,  /* tensile fracture strain             */
               DOUBLE *cappae,   /* max. elastic compression strain     */
               DOUBLE *cappauc,  /* compressive fracture strain         */
               DOUBLE *sig3,     /* equivalent compressive stress       */
               DOUBLE *fbd,      /* tension stiffening stress           */
               DOUBLE **d,       /* elastic material matrix             */
               DOUBLE *sig,      /* stresses from last update           */
               DOUBLE *eps);     /* strains from last update            */
/*-----------------------------------------------------------------------|
|  w1_mat_plast_serv.c                                       al 9/01     |
|       topic : radial return for elements with elastoplastic material   |
|               (concrete)                                               |
|-----------------------------------------------------------------------*/
void w1cradi (DOUBLE *sigma,  /* elastic predictor projected onto yield */ 
              DOUBLE *epstn,  /* equivalent uniaxial plastic strain     */ 
              DOUBLE *dlam,   /* incremental plastic multiplier         */ 
              WALL_TYPE wtype,/* plane stress/strain...                 */
              INT    yip,     /* stress state   1 =elastic >=2 =plastic */ 
              DOUBLE *alpha,  /* factor for the first invariants        */ 
              DOUBLE *ft,     /* yield condition                        */ 
              DOUBLE *e,      /* young modulus                          */ 
              DOUBLE *g,      /* shear modulus                          */ 
	        DOUBLE *com,    /* bulk modulus                           */ 
              DOUBLE *sigym,  /* uniaxial predictor yield stress        */ 
              DOUBLE *hards,  /* plastic modulus                        */ 
              DOUBLE *sigy,   /* actual uniaxial yield stress           */ 
              DOUBLE *dn,     /* gradient components in deviatoric dir  */ 
              DOUBLE *dcom,   /* gradient components in hdrostatic dir  */ 
              DOUBLE *devsig, /* deviatoric predictor stresses          */ 
              DOUBLE *sm,     /* hydrostatic predictor stresses         */ 
              DOUBLE *fcm,    /* compressive strenght                   */ 
              DOUBLE *gc,     /* compression fracture energy            */ 
              DOUBLE *ftm,    /* tensile strenght                       */ 
	        DOUBLE *gt,     /* tensile fracture energy                */ 
              DOUBLE *gamma1, /* fitting factor yield function 1        */ 
              DOUBLE *gamma2, /* fitting factor yield function 2        */ 
              DOUBLE *dia,    /* equivalent element length              */ 
              DOUBLE *acrs);  /* average crack spacing                  */ 
/*-----------------------------------------------------------------------*
|  w1_mat_plast_serv.c                                       al 9/01     |
|                                                                        |
|    forms the elasto-plastic consistent tangent material tensor         |
*-----------------------------------------------------------------------*/
void w1mapl2(DOUBLE *tau,      /* current stresses (local)              */
             DOUBLE **d,       /* material matrix to be calculated      */
             DOUBLE *dlam,     /* increment of plastic multiplier       */
             WALL_TYPE wtype,  /* plane stress/strain...                */
             DOUBLE *alpha,    /* neigungswinkel der fliessflaechen     */
             DOUBLE *g,        /* schubmodul                            */
             DOUBLE *com,      /* kompressionsmodul                     */
             DOUBLE *hards,    /* plastic hardeningmodulus              */
             DOUBLE *dn,       /* gradient components in dev. direction */
             DOUBLE *grad);    /* total gradient                        */
/*-----------------------------------------------------------------------|
|  w1_mat_plast_serv.c                                       al 9/01     |
|       topic : radial return for multisurface problems                  |
|               (material model: concrete)                               |
|-----------------------------------------------------------------------*/
void w1cradms(DOUBLE *sigma,  /* elastic predictor projected onto yield */ 
              DOUBLE *epstn,  /* equivalent uniaxial plastic strain     */ 
              DOUBLE *dlam,   /* incremental plastic multiplier         */ 
              WALL_TYPE wtype,/* plane stress/strain...                 */
              INT     yip,    /* stress state   1 =elastic >=2 =plastic */ 
              DOUBLE *alpha,  /* factor for the first invariants        */ 
              DOUBLE *ft,     /* yield condition                        */ 
              DOUBLE  e,      /* young modulus                          */ 
              DOUBLE  g,      /* shear modulus                          */ 
	      DOUBLE  com,    /* bulk modulus                           */ 
              DOUBLE *sigym,  /* uniaxial predictor yield stress        */ 
              DOUBLE *hards,  /* plastic modulus                        */ 
              DOUBLE  dn[3][4],  /*gradient components in deviatoric dir*/ 
              DOUBLE  dcom[3][4],/*gradient components in hdrostatic dir*/ 
              DOUBLE *devsig, /* deviatoric predictor stresses          */ 
              DOUBLE *sm,     /* hydrostatic predictor stresses         */ 
              DOUBLE  fcm,    /* compressive strenght                   */ 
              DOUBLE  gc,     /* compression fracture energy            */ 
              DOUBLE  ftm,    /* tensile strenght                       */ 
	      DOUBLE  gt,     /* tensile fracture energy                */ 
              DOUBLE  gamma1, /* fitting factor yield function 1        */ 
              DOUBLE  gamma2, /* fitting factor yield function 2        */ 
              DOUBLE  dia,    /* equivalent element length              */ 
              DOUBLE  acrs);  /* average crack spacing                  */ 
/*-----------------------------------------------------------------------|
|  w1_mat_plast_serv.c                                       al 9/01     |
|     topic: convert some arrays                                         |
|-----------------------------------------------------------------------*/
      void w1conver(DOUBLE *dlam,       /*  incr. of plastic multiplier */             
                    DOUBLE *alpha,      /*  factor for the first inv.   */             
                    DOUBLE *hards,      /* hardening modulus            */             
                    DOUBLE dn[3][4],    /*  gradient comp. in dev. dir. */  
                    DOUBLE grad[3][4],  /*  total gradient              */             
                    DOUBLE *dlamc,      /*  -|                          */             
                    DOUBLE *alphac,     /*   |                          */             
                    DOUBLE *hardsc,     /*   |-  converted arrays       */             
                    DOUBLE dnc[2][4],   /*   |                          */             
                    DOUBLE gradc[2][4]);/*  -|                          */             
/*---------------------------------------------------------------------- |
|  w1_mat_plast_serv.c                                       al 9/01     |
|    forms the elasto-plastic consistent tangent material tensor         |
|    for wall element   (drucker - prager)                               |
|    (generalized for the multisurface problem: including the apex)      |
*-----------------------------------------------------------------------*/
void w1maplg(DOUBLE *tau,      /* current stresses (local)              */
             DOUBLE **d,       /* material matrix to be calculated      */
             DOUBLE  dlam,     /* increment of plastic multiplier       */
             WALL_TYPE wtype,  /* plane stress/strain...                */
             INT     yip,      /* stress state  1=elastic 2=plastic     */
             DOUBLE  emod,     /* elastizitaetsmodul                    */
             DOUBLE  g,        /* schubmodul                            */
             DOUBLE  com,      /* kompressionsmodul                     */
             DOUBLE *hards,    /* plastic hardeningmodulus              */
             DOUBLE  dn[2][4], /* gradient components in dev. direction */
             DOUBLE  grad[2][4]);/*norm of the dev. predictor stresses  */
/*-----------------------------------------------------------------------|
|  w1_mat_plast_serv.c                                        al 9/01    |
|       topic : radial return for elements with elastoplastic material   |
|               (concrete)      CAP-REGION                               |
|-----------------------------------------------------------------------*/
void w1radcap(DOUBLE *sigma,  /* elastic predictor projected onto yield */ 
              DOUBLE *epstn,  /* equivalent uniaxial plastic strain     */ 
              DOUBLE *dlam,   /* incremental plastic multiplier         */ 
              WALL_TYPE wtype,/* plane stress/strain...                 */
              DOUBLE *alpha,  /* factor for the first invariants        */ 
              DOUBLE *e,      /* young modulus                          */ 
              DOUBLE *g,      /* shear modulus                          */ 
	        DOUBLE *com,    /* bulk modulus                           */ 
              DOUBLE *sigym,  /* uniaxial predictor yield stress        */ 
              DOUBLE *hards,  /* plastic modulus                        */ 
              DOUBLE *grad,   /* total gradient                         */ 
              DOUBLE *devsig, /* deviatoric predictor stresses          */ 
              DOUBLE *dev,    /* norm of the dev. predictor stresses    */ 
              DOUBLE *hyd,    /* 1st invariant  of the predictor stress.*/ 
              DOUBLE *hydn,   /* 1st invariant  of the new stresses     */ 
              DOUBLE *fcm,    /* compressive strenght                   */ 
              DOUBLE *gc,     /* compression fracture energy            */ 
              DOUBLE *gamma2, /* fitting factor yield function 2        */ 
              DOUBLE *dia);   /* equivalent element length              */ 
/*-----------------------------------------------------------------------|
|  w1_mat_plast_serv.c                                        al 9/01    |
|       topic: forms the  e l a s t o - p l a s t i c                    |
|              c o n s i s t e n t  tangent material tensor              |
|              for wall element   (drucker - prager)                     |
|-----------------------------------------------------------------------*/
void w1mplcap(DOUBLE *tau,      /* current stresses (local)             */
              DOUBLE **d,       /* material matrix to be calculated     */
              DOUBLE *dlam,     /* increment of plastic multiplier      */
              WALL_TYPE wtype,  /* plane stress/strain...               */
              DOUBLE *alpha,    /* factor for the first invariants      */
              DOUBLE *emod,     /* young modulus                        */
              DOUBLE *vnu,      /* poisson's ratio                      */
              DOUBLE *hards,    /* plastic hardeningmodulus             */
              DOUBLE *sigym,    /* [0] uniaxial tension yield stress    */
              DOUBLE *grad,     /* total gradient                       */
              DOUBLE *cappae,   /* max. elastic compression strain      */
              DOUBLE *cappauc,  /* compressive fracture strain          */
              DOUBLE *epst,     /* equivalent uniaxial plastic strain   */
              DOUBLE *sig3);    /* equivalent compressive stress        */
/*-----------------------------------------------------------------------|
|  w1_mat_plast_serv.c                                        al 9/01    |
| compute average crack spacing for element wall1     (model: concrete)  |
|-----------------------------------------------------------------------*/
void w1acrs(DOUBLE  hte,       /* section properties                    */   
               INT  maxreb,    /* max. number of rebars                 */   
            DOUBLE *stress,    /* element stresses                      */   
            DOUBLE *angle,     /* angle of the principal concrete stress*/   
            DOUBLE *reb_area,  /* rebar area                            */
            DOUBLE *reb_ang,   /* rebar angle                           */
            DOUBLE *reb_so,    /* minimum bond length                   */
            DOUBLE *reb_ds,    /* diameter                              */
            DOUBLE *reb_rgamma,/* 4 := deformed bars = 2 := plane bar   */
            DOUBLE *thick,     /* thickness ot the structure            */   
            DOUBLE dia,        /* equivalent element length             */   
            DOUBLE *acrs);     /* average crack spacing                 */   
/*-----------------------------------------------------------------------|
|  w1_mat_plast_serv.c                                        al 9/01    |
|    topic: calculate the gradients of the elastic predictor             |              
|           deviatoric/hydrostatic stresses                              |
|-----------------------------------------------------------------------*/
void w1cpreds(DOUBLE *sigym, /* uniaxial predictor yield stress         */
              DOUBLE *alpha, /* factor for the first invariants         */
              DOUBLE *devsig,/* deviatoric predictor stresses           */
              DOUBLE  hyd,   /* 1. invariant of the predictor stresses  */
              DOUBLE *grad); /* total gradient                          */
/*-----------------------------------------------------------------------*
|  w1_mat_plast_serv.c                                       al 9/01     |
|    WALL1 : array = value                                               |
*-----------------------------------------------------------------------*/
void w1_init44(DOUBLE a[4][4],DOUBLE val);
/*-----------------------------------------------------------------------*
|  w1_mat_plast_serv.c                                       al 9/01     |
|    WALL1 : R(I,J) = A(I,K)*B(K,J) ---  R = A*B                         |
*-----------------------------------------------------------------------*/
void w1_AxB_444(DOUBLE a[4][4],DOUBLE b[4][4],DOUBLE r[4][4]);
/*-----------------------------------------------------------------------*
|  w1_mat_plast_serv.c                                       al 9/01     |
|    WALL1 : R(I,J) = A(I,K)*B(K,J) ---  R = A*B                         |
*-----------------------------------------------------------------------*/
void w1_AxBT_444(DOUBLE a[4][4],DOUBLE b[4][4],DOUBLE r[4][4]);
/*-----------------------------------------------------------------------*
|  w1_mat_plast_serv.c                                       al 9/01     |
|    WALL1 : R(I,J) = A(I,K)*B(K,J) ---  R = A*B                         |
*-----------------------------------------------------------------------*/
void w1_AxB_441(DOUBLE a[4][4],DOUBLE b[4],DOUBLE r[4]);
/*-----------------------------------------------------------------------*
|  w1_mat_plast_serv.c                                       al 9/01     |
|    WALL1 : R(I,J) = A(I,K)*B(K,J) ---  R = A*B                         |
*-----------------------------------------------------------------------*/
void w1_AxBT_414(DOUBLE a[4],DOUBLE b[4],DOUBLE r[4][4]);
/*-----------------------------------------------------------------------*
|  w1_mat_plast_serv.c                                       al 9/01     |
|    WALL1 : R(I,J) = A(I,K)*B(K,J) ---  R = A*B                         |
*-----------------------------------------------------------------------*/
void w1_AxFAC_44(DOUBLE a[4][4], DOUBLE r[4][4], DOUBLE fac);
/*-----------------------------------------------------------------------*
|  w1_mat_plast_serv.c                                       al 9/01     |
|    WALL1 : ADDITION   R = A+B   (FULL MATRICES)                        |
*-----------------------------------------------------------------------*/
void w1_AaddB_44(DOUBLE a[4][4], DOUBLE b[4][4], DOUBLE r[4][4]);
/*-----------------------------------------------------------------------*
|  w1_mat_plast_serv.c                                       al 9/01     |
|  constitutive matrix - forces - plastic - rebar                al 9/01 |
*-----------------------------------------------------------------------*/
void w1_mat_rebar(ELEMENT   *ele, /* actual element                     */
                  MATERIAL  *mat, /* actual material                    */
                  DOUBLE   **bop, /* derivative operator                */
                  DOUBLE    *gop, /* add. derivative operator           */
                  DOUBLE    *alpha,/* internal dof                      */
                  DOUBLE   **xjm, /* jacobian matrix                    */
                  DOUBLE *stress, /* vector of stresses                 */
                  DOUBLE     **d, /* constitutive matrix                */
                  INT         ip, /* integration point Id               */
                  INT       lanz, /* rebar id                           */
                  INT     istore);/* controls storing of new stress.2wa */
/*-----------------------------------------------------------------------|
|  w1_mat_serv.c                                              al 9/01    |
|      topic: condensed stress vector                                    |
|             plane strain --> plane stress                              |
|-----------------------------------------------------------------------*/
void w1consig (DOUBLE **d,/* current material matrix components d14-d44 */
               DOUBLE  *sigma,  /* current stresses (local)             */
               DOUBLE  *sigmac);/* condensed stress vector              */
/*-----------------------------------------------------------------------|
|  w1_mat_serv.c                                              al 9/01    |
|      topic : condense the constitutive tensor                          |
|              plane strain --> plane stress                             |
|-----------------------------------------------------------------------*/
void w1concep (DOUBLE **d);  /* material matrix to be calculated        */
/*-----------------------------------------------------------------------|
|  w1_mat_serv.c                                              al 9/01    |
|    topic: evaluate the incremental strain in thickness direction       |              
|           for wall element (plane stress)                              |
|-----------------------------------------------------------------------*/
void w1de33(DOUBLE *sigi,  /* stresses from last iteration step         */
            DOUBLE *epsi,  /* strains from last iteration step          */
            DOUBLE *di,    /* components d41,d42,d43,d44 of the         */
                           /* const. tensor from last iteration step    */
            DOUBLE *strain);/* current strains (local)                  */
/*-----------------------------------------------------------------------*
|  w1_static_ke.c                                            al 9/01     |
|  integration of linear stiffness ke for wall1 element                  |
*-----------------------------------------------------------------------*/
void w1static_ke(ELEMENT   *ele,    /* actual element                   */
                 W1_DATA   *data,   /* wall1 data                       */
                 MATERIAL  *mat,    /* actual material                  */
                 ARRAY     *estif_global, /* element stiffness matrix   */
                 ARRAY     *emass_global, /* element mass matrix        */
                 DOUBLE    *force,  /* global vector for internal forces*/  
                 INT        init);  /* initialize this function         */
/*-----------------------------------------------------------------------*
|  w1_static_keug.c                                        ah 06/02      |
|  integration of nonlinear stiffness keug for wall1 element             |
*-----------------------------------------------------------------------*/
void w1static_keug(ELEMENT  *ele,           /* actual element           */
                   W1_DATA   *data,         /* wall1 data               */
                   MATERIAL  *mat,          /* actual material          */
                   ARRAY     *estif_global, /* element stiffness matrix */
                   ARRAY     *emass_global, /* element mass matrix      */
                   DOUBLE    *force,/* global vector for internal forces*/
                   INT        init);        /* initialize this function */
/*-----------------------------------------------------------------------*
|  w1_static_ke.c                                            al 9/01     |
|  evaluates element forces                                              |
*-----------------------------------------------------------------------*/
void w1fi( DOUBLE  *F,    /* element stresses                           */
           DOUBLE   fac,  /* integration factor                         */
           DOUBLE **bop,  /* derivative operator                        */
           INT      nd,   /* total number degrees of freedom of element */
           DOUBLE  *fie); /* nodal forces                               */           
/*----------------------------------------------------------------------*/
/*  w1_call_mat.c                                         ah 06/02      */
/* get density out of material law                                      */
/*----------------------------------------------------------------------*/
void w1_getdensity(MATERIAL  *mat,        /* actual material            */
                   DOUBLE *density);      /* density of actual material */
/*----------------------------------------------------------------------*/
/*  w1_boplin.c                                           ah 06/02      */
/* evaluation of linear B-operator                                      */
/*----------------------------------------------------------------------*/
void w1_boplin(DOUBLE    **boplin,        /* Blin                       */
               DOUBLE    **deriv,     /* derivatives of ansatzfunctions */
               DOUBLE    **xjm,           /* jacobian matrix            */
               DOUBLE      det,           /* det of jacobian matrix     */
               INT         iel);          /* nodenumber of element      */

/*----------------------------------------------------------------------*/
/*  w1_b_barop.c                                                        */
/* evaluation of B_bar operator                                         */
/*----------------------------------------------------------------------*/
void w1_b_barop(ELEMENT *ele,                          /* actual element*/
                DOUBLE **b_bar,                        /* b_bar operator*/
		DOUBLE **int_b_bar, /*interpolated b-bar op.(E-M scheme)*/    
		DOUBLE **boplin,                       /*boplin operator*/
                DOUBLE *F,                        /*Deformation gradient*/
		INT numeps,               /* number of strain components*/
		INT nd,       /* total number degrees of freedom of ele.*/
		INT ip);                     /*Integration point counter*/
/*----------------------------------------------------------------------*/
/*  w1_defgrad.c                                          ah 06/02      */
/*  evaluation of deformation gradient F                                */
/*----------------------------------------------------------------------*/
void w1_defgrad(DOUBLE    *F,        /* Deformation gradient            */
                DOUBLE    *strain,   /* Green-Lagrange-strains          */
                DOUBLE    **xrefe,   /* coordinates in referenz-config. */
                DOUBLE    **xcure,   /* coordinates in current-config.  */
                DOUBLE    **boplin,  /* Blin                            */
                INT         iel);    /* nodenumber of element           */
/*----------------------------------------------------------------------*
|  w1_call_matgeononl.c                                     ah 06/02     |
|  select proper material law for large deformations                     |
*-----------------------------------------------------------------------*/
void w1_call_matgeononl(MATERIAL  *mat,   /* material of actual element */
                        WALL_TYPE wtype,/*Info about actual wall element*/
                        DOUBLE *strain,  /*strain at act. gaussian point*/
                        DOUBLE **stress, /*stress at act. gaussian point*/
                        DOUBLE **d,       /* material tangent           */
                        INT numeps);     /* number of strain components */
/*----------------------------------------------------------------------*
|  w1_mat_linelgeonon.c                                     ah 06/02     |
|  linear elastic material law for large deformations                    |
*-----------------------------------------------------------------------*/
void w1_mat_linelgeonon(DOUBLE ym,       /* Young's modulus             */
                        DOUBLE pv,       /* poisson's ratio             */
                        WALL_TYPE wtype,/*Info about actual wall element*/
                        DOUBLE *strain,  /*strain at act. gaussian point*/
                        DOUBLE **d,      /* material tangent            */
                        DOUBLE **stress, /*stress at act. gaussian point*/
                        INT numeps);     /* number of strain components */
/*-----------------------------------------------------------------------*
|  w1_cal_kg.c                                              ah 06/02     |
|  evaluation of geometric part of stiffness matrix in geononl. case     |
*-----------------------------------------------------------------------*/
void w1_kg(
           
	   ELEMENT  *ele,                /* active element pointer      */ 	   
           DOUBLE  **kg,                 /* geometric stiffness matrix  */
           DOUBLE  **boplin,             /* linear B-operator           */
           DOUBLE  **stress,             /*stress at act. gaussian point*/
           DOUBLE    fac,                /* integration factor          */
           INT       nd,                 /* dof's of element            */
           INT       neps,               /* number of strain components */ 
	   INT       ip);                /* active Gauss point          */
	                
/*-----------------------------------------------------------------------*
|  w1_cal_keu.c                                              ah 06/02    |
|  evaluation of elast+init. disp. part of stiffness in geononl. case    |
*-----------------------------------------------------------------------*/
void w1_keu(DOUBLE  **keu,    /* elastic + initial deformation stiffness*/
            DOUBLE  **b_bar,             /* b_bar operator              */            
            DOUBLE  **int_b_bar,         /* interpolated b_bar(for GEMM)*/
            DOUBLE  **D,                 /* material tangente           */
            DOUBLE    fac,               /* integration factor          */
            INT       nd,                /* dof's of element            */
            INT       neps);             /* number of strain components */
	                 
 /*----------------------------------------------------------------------*
 |  w1_cal_fint.c                                              ah 06/02   |
 | evaluate internal element forces for large def (total Lagr)           |
 *----------------------------------------------------------------------*/
void w1_fint(ELEMENT *ele,             /* active element pointer       */
             DOUBLE **stress,          /* 2.PK stresses                */ 
             DOUBLE **int_b_bar,       /* interpolated b_bar (for GEMM)*/ 
             DOUBLE  *fint,            /* internal forces              */ 
             DOUBLE   fac,             /* detJ*wr*ws*thickness         */ 
             INT      nd,              /* Element-DOF                  */
	       INT      ip);             /* Active Gauss point           */  
/*----------------------------------------------------------------------*
 | Transform stress and strain local-global                  fh 7/02    |
 | Local 3-direction is zero                                            |
 *----------------------------------------------------------------------*/ 
void w1_lss(DOUBLE    *a,    /* vector to be transformed                */
            DOUBLE    **G,   /* elements of transformation matrix       */
            DOUBLE    **GI,  /* inverse of G                            */
	    INT         it); /* flag for local/global strains/stresses  */                                                                                                                                              
/*----------------------------------------------------------------------*
 | Set Transformation Matrices G and GI                      fh 7/02    |
 *----------------------------------------------------------------------*/
void w1_sett(DOUBLE   **A,   /* matrix of direction cosines             */
            DOUBLE    **B,   /* G-Matrix     (L=1)                      */
            DOUBLE    **C);  /* inverse of G (L=2)                      */                                                         
/*----------------------------------------------------------------------*
 | Calculate Transformation Matrices G and G(Inv)            fh 7/02    |
 *----------------------------------------------------------------------*/
void w1_tram(DOUBLE   **xjm, /* Elements of Jacobian Matrix             */
            DOUBLE    **G,   /* G-Matrix                                */
            DOUBLE    **GI,  /* inverse of G                            */
	    DOUBLE    **dum);/* matrix of direction cosines             */                                                                     
/*----------------------------------------------------------------------*
 | w1_cal_fext.c                                             ah 07/02   |
 | element line and surface loads                                       |
 *----------------------------------------------------------------------*/
void w1_eleload(ELEMENT  *ele,     /* actual element                  */
                W1_DATA  *data,    /* wall1- Data                     */
                DOUBLE	 *loadvec, /* global element load vector fext */
                INT	  init,    /* flag if init or calculation     */
                INT       imyrank);        
void w1_fsiload(ELEMENT  *ele,     /* actual element                  */
                W1_DATA  *data,    /* wall1- Data                     */
                DOUBLE	 *loadvec, /* global element load vector fext */
                INT	  init,    /* flag if init or calculation     */
                INT       imyrank);
void w1_iedg(INT *iegnod, ELEMENT *ele, INT line, INT init);		
/*----------------------------------------------------------------------*
 | w1_cal_fext.c                                             ah 07/02   |
 | integration of element loads                                         |
 *----------------------------------------------------------------------*/
void w1_fextsurf(ELEMENT    *ele,        /* actuell element             */
                 DOUBLE    **eload,      /* element load vector         */
                 DOUBLE     *funct,      /* ansatz-functions            */
                 DOUBLE      fac,        /* integration factor          */
                 INT         iel);       /* element DOF                 */
/*----------------------------------------------------------------------*
 | w1_bop.c                                                  ah 9/02    |
 | additional operator matrix gop at r,s for incomp modes               |
 *----------------------------------------------------------------------*/
void w1_gop(DOUBLE    *gop,            /* operator matrix G             */
            DOUBLE    **xjm0,          /* jacobian matrix at r,s=0      */
            DOUBLE      det0,          /* det J at r,s=0                */
            DOUBLE      det,           /* det J at actual r,s           */
            DOUBLE      e1,            /* actual GP coordinate r        */
            DOUBLE      e2);           /* actual GP coordinate s        */
/*----------------------------------------------------------------------*
 | w1_incompmode.c                                              ah 9/02 |
 | stiffness matrix due to incompatible modes                           |
 *----------------------------------------------------------------------*/
void w1_knc(DOUBLE  **knc,             /* stiffness knc= BT C G         */
            DOUBLE  **bop,             /* operator matrix               */
            DOUBLE   *gop,             /* additional opperator matrix   */
            DOUBLE  **d,               /* constitutive matrix           */
            DOUBLE    fac);            /* integration factor            */
/*----------------------------------------------------------------------*
 | w1_incompmode.c                                              ah 9/02 |
 | stiffness matrix  knn      (4*4)                                     |
 *----------------------------------------------------------------------*/
void w1_knn(DOUBLE  **knn,             /* stiffness knn= GT C G         */
            DOUBLE   *gop,             /* additional opperator matrix   */
            DOUBLE  **d,               /* constitutive matrix           */
            DOUBLE    fac);            /* integration factor            */
/*----------------------------------------------------------------------*
 | w1_incompmode.c                                              ah 9/02 |
 | evaluate internal element forces due to incomp. modes                |
 *----------------------------------------------------------------------*/
void w1_fintn(DOUBLE  *F,              /* stress                        */
              DOUBLE   fac,            /* integration factor            */
              DOUBLE  *gop,            /* additional opperator matrix   */
              DOUBLE  *fintn);          /* INT forces due to inc modes  */
/*----------------------------------------------------------------------*
 | w1_math.c                                                    ah 9/02 |
 | calculate inverse of arbitrary NxN matrix                            |
 *----------------------------------------------------------------------*/
void w1_inverse_matrix
 (INT N,          /* I: size of matrix: N x N                           */
  DOUBLE **A,     /* I: Original Matrix                                 */
  DOUBLE **Y);    /* O: Inverse Matrix                                  */
/*----------------------------------------------------------------------*
 | w1_math.c                                                    ah 9/02 |
 | LU-Decomposition of  NxN - matrix ->"numerical recipes in C"         |
 *----------------------------------------------------------------------*/
void w1_ludcmp(DOUBLE **a,  /* I: Original Matrix (will be rearanged)   */
               INT n,       /* I: size of matrix: N x N                 */
               INT *indx,   /* O:row premutation pf pivoting            */
               DOUBLE *d);/* O:+1/-1,depending on even or odd row permut*/
/*----------------------------------------------------------------------*
 | w1_math.c                                                    ah 9/02 |
 | solving AX=B, A is LU decomposed  ->"numerical recipes in C"         |
 *----------------------------------------------------------------------*/
void w1_lubksb(DOUBLE **a,   /* I: LU-decomposed matrix                 */
               INT n,        /* I: size of matrix: N x N                */
               INT *indx,    /* I:row premutation pf pivoting           */
               DOUBLE *b);   /* In:B will be changed-> Out:Solution X   */
/*----------------------------------------------------------------------*
 | w1_math.c                                                    ah 9/02 |
 |   allocate a DOUBLE vector with subscript range v[nl..nh]            |
 *----------------------------------------------------------------------*/
DOUBLE *w1_vector(INT nl, INT nh);
/*----------------------------------------------------------------------*
 | w1_math.c                                                    ah 9/02 |
 |   free a DOUBLE vector allocated with w1_vector()                    |
 *----------------------------------------------------------------------*/
void w1_free_vector(DOUBLE *v, INT nl);
/*----------------------------------------------------------------------*
 | w1_incompmode.c                                              ah 9/02 |
 | static kondensation:     Kele =   K  - knc*inverse(knn)*kcn,         |
 |                      Fint ele = Fint - knc*inverse(knn)*fintn,       |
 *----------------------------------------------------------------------*/
void  w1_stat_cond(DOUBLE **knninv,/*I:stiffness inverse(knn) (4x4)     */
                   DOUBLE **knc,   /*I: knc = BT C G (8x4)              */
                   DOUBLE **deltak,/*O: knc*inverse(knn)*kcn            */ 
                   DOUBLE  *fintn, /*I: fintn=GT*singma                 */
                   DOUBLE  *deltaf,/*O: knc*inverse(knn)*fintn          */ 
                   ELEMENT *ele); /*I: actual element                  */
/*----------------------------------------------------------------------*
 | w1_incompmode.c                                              ah 9/02 |
 | update of internal dof alpha                                         |
 *----------------------------------------------------------------------*/
void  w1_updalpha(DOUBLE  *alpha,  /*O: internal dof for incomp modes   */
                  ELEMENT *ele,    /*I:actual element                   */
                  DOUBLE **knc,    /*I: mixed stiffness                 */
                  DOUBLE **knninv, /*I: inverse of incomp. stiffness    */
                  DOUBLE  *fintn,  /*I: INT. forces of inc. modes       */
                  INT      istore);/*I: update after loadstep->istore=1 */
/*----------------------------------------------------------------------*
 |  integration routine for WALL1 element                      al 6/01  |
 *----------------------------------------------------------------------*/
void w1_oint(
             ELEMENT   *ele, 
             W1_DATA   *data, 
             MATERIAL  *mat,
             DOUBLE    *retval,  /* return value */
             INT        init     /* ==2 calc.strain energy */
             );
/*-----------------------------------------------------------------------*
|  w1_mat_plast_mises3D.c                                      sh 8/02   |
|  constitutive matrix - forces - linear elastic- von Mises - 3D         |
|  uses a 3D-Material formulation; constitutive matrix and stresses      |
|  are condensed back to 2D conditions (plane stress, plane strain)      |
|  (rotational symmetry is not implemented)                              |
*-----------------------------------------------------------------------*/
void w1_mat_plast_mises_3D(DOUBLE     ym,
                           DOUBLE     pv,
                           DOUBLE     sigy,
                           DOUBLE     hard,
                           DOUBLE     gf,
                           DOUBLE     betah,
                           ELEMENT   *ele,
                           WALL_TYPE  wtype,
                           DOUBLE   **bop,
                           DOUBLE    *gop,
                           DOUBLE    *alpha,
                           INT        ip,
                           DOUBLE    *stressc,
                           DOUBLE   **d,
                           INT        istore,
                           INT        newval);
/*-----------------------------------------------------------------------|
| w1_mat_trans.c                                                         |
|      topic: blowing up plane stress/strain conditions           sh 7/02|
|             to 3D --> 3D-Material Law                                  |
|      topic: kondense 3D conditions                              sh 7/02|
|             to plane stress/strain conditions                          |
|-----------------------------------------------------------------------*/
void w1mat_trans_up (DOUBLE     ym,
                     DOUBLE     pv,
                     ELEMENT   *ele,                                        
                     WALL_TYPE  wtype,
                     DOUBLE   **bop,
                     DOUBLE    *gop,
                     DOUBLE    *alpha,
                     INT        ip,
                     DOUBLE     strain[4]); 
void w1mat_trans_down (DOUBLE   **d, 
                       ELEMENT   *ele,
                       WALL_TYPE  wtype,
                       INT        ip,
                       INT        yipc,
                       DOUBLE    *stressc, 
                       DOUBLE    *sig3D,
                       DOUBLE    *eps3D,
                       DOUBLE    *stress3D,
                       DOUBLE    *strain,  
                       DOUBLE    *qn3D);
/*-----------------------------------------------------------------------*
| w1_mat_trans.c                                                         |
|     changes 2 rows of a vector                               sh 8/02   |
|  vec[x,x,a,x,b,x,...] -> vec[x,x,b,x,a,x,...]                          |
*-----------------------------------------------------------------------*/
void w1_vec_switch(DOUBLE *vec,      /* vector do be modified           */
                   INT a,            /* row to be changed to b          */
                   INT b);           /* row to be changed to a          */
/*-----------------------------------------------------------------------*
| w1_mat_trans.c                                                         |
|     changes 2 rows & columns of a square matrix              sh 8/02   |
|     [-,a,-,b,-]        [-,b,-,a,-]                                     |
|     [a,a,a,c,a]        [b,b,b,c,b]                                     |
|  mat[-,a,-,b,-] ->  mat[-,b,-,a,-]                                     |
|     [b,c,b,b,b]        [a,c,a,a,a]                                     |
|     [-,a,-,b,-]        [-,b,-,a,-]                                     |
*-----------------------------------------------------------------------*/
void w1_matrix_switch(DOUBLE **mat,   /* matrix do be modified          */
                      INT a,          /* row & colum to be changed to b */     
                      INT b,          /* row & colum to be changed to a */
                      INT l);         /* length of row/column of matrix */
/*----------------------------------------------------------------------*/
void w1_mat_dam_mp(DOUBLE     youngs,
                   DOUBLE     nue,
                   DOUBLE     kappa_0,
                   DOUBLE     alph,
                   DOUBLE     beta,
                   ELEMENT   *ele,
                   WALL_TYPE  wtype,
                   DOUBLE   **bop,
                   DOUBLE    *gop,
                   DOUBLE    *alpha,
                   INT        ip,
                   DOUBLE    *stress,
                   DOUBLE   **D,
                   INT        istore,
                   INT        newval); 
/*----------------------------------------------------------------------*
 | constitutive matrix - forces - linear elastic- Damage - 2D   he 04/03|
 | plane stress, plane strain, rotational symmetry                      |
 *----------------------------------------------------------------------*/
void w1_mat_damage(DOUBLE ym,      /* young's modulus                   */
                   DOUBLE pv,      /* poisson's ratio                   */
                   INT    Equival, /* flag for equivalent strains       */
                   INT    Damtyp,  /* flag for Damage-Typ               */
                   DOUBLE Kappa_0, /* initial damage equivalent strain  */
                   DOUBLE Kappa_m, /* factor for damage-law             */
                   DOUBLE Alpha,   /* factor for expon. damage-function */       
                   DOUBLE Beta,    /* factor for expon. damage-function */
                   DOUBLE k_fac,   /* factor for de Vree                */       
                   ELEMENT   *ele, /* actual element                    */
                   WALL_TYPE wtype,/* plane stress/strain...            */         
                   DOUBLE **bop,   /* derivative operator               */
                   DOUBLE  *gop,
                   DOUBLE  *alpha,
                   INT ip,         /* integration point Id              */
                   DOUBLE *stress, /* vector of stresses                */
                   DOUBLE **d,     /* constitutive matrix               */
                   INT istore,     /* controls storing of stresses      */
                   INT newval);     /* controls eval. of stresses        */ 
/*----------------------------------------------------------------------*
 | create tensor from vector !!!!!!only for strains!!!!   he    04/03   |
 *----------------------------------------------------------------------*/
void w1_4to9(DOUBLE  *vector, 
             DOUBLE   tensor[3][3]);
/*----------------------------------------------------------------------*
 | compute Kroneker-Delta (2-stufiger Einheitstensor)      he    04/03   |
 *----------------------------------------------------------------------*/
void w1_kroneker(DOUBLE delta[3][3]);
/*----------------------------------------------------------------------*
 | compute ela. stresses with Hook                        he    04/03   |
 *----------------------------------------------------------------------*/
void w1_stress_ela(DOUBLE ym, 
                   DOUBLE pv, 
                   DOUBLE epsilon[3][3], 
                   DOUBLE sigma_el[3][3],
                   DOUBLE delta[3][3]);
/*----------------------------------------------------------------------*
 | compute actuel equivalent strains and                  he    04/03   |
 | their global derivatives of epsilon                                  |
 | Equival = 1 -  ENERGY RELEASE RATE CONCEPT - THERMODYNAMICS          |
 | Equival = 2 -  SIMO & JU (1987)                                      |
 | Equival = 3 -  JU (1989)                                             |
 | Equival = 4 -  de VREE (1995)                                        |
 *----------------------------------------------------------------------*/
void w1_equi_eps(DOUBLE   epsilon[3][3], 
                 DOUBLE   sigma_el[3][3],       
                 DOUBLE   delta[3][3], 
                 INT      Equival,
                 DOUBLE   *eta, 
                 DOUBLE   eta_der[3][3],
                 DOUBLE   pv,
                 DOUBLE   k);
/*----------------------------------------------------------------------*
 | BESTIMMUNG DER SCHAEDIGUNG UND IHRER PARTIELLEN        he    04/03   |
 | ABLEITUNG NACH KAPPA JE NACH DAMTYPE                                 |
 | DAMTYPE=1 - LINEAR-ENTFESTIGENDES DAMAGING                           |
 | DAMTYPE=2 - EXPONENTIELLES DAMAGING                                  |
 *----------------------------------------------------------------------*/
void w1_dam_typ(DOUBLE *damage, 
                DOUBLE *dam_deriv,
                DOUBLE kappa,
                DOUBLE Kappa_0,
                DOUBLE Kappa_m,
                DOUBLE Alpha,
                DOUBLE Beta,
                INT    Damtyp);
/*----------------------------------------------------------------------*
 | compute ela. stiffness tensor                           he    04/03   |
 *----------------------------------------------------------------------*/
void w1_mat_ela(DOUBLE   ym, 
                DOUBLE   pv,
                DOUBLE   delta[3][3], 
                DOUBLE   c_el[3][3][3][3]);
/*----------------------------------------------------------------------*
 | BESTIMMUNG DES SEKANTENTENSORS                         he    04/03   |
 *----------------------------------------------------------------------*/
void w1_sec(DOUBLE damage, 
            DOUBLE c_el[3][3][3][3],
            DOUBLE c_sec[3][3][3][3]);
/*----------------------------------------------------------------------*
 | BESTIMMUNG DES SPANNUNGSTENSORS                        he    04/03   |
 *----------------------------------------------------------------------*/
void w1_stress(DOUBLE c_sec[3][3][3][3], 
               DOUBLE epsilon[3][3],
               DOUBLE sigma[3][3]);
/*----------------------------------------------------------------------*
 | reduce tensor to vector                                 he    04/03   |
 *----------------------------------------------------------------------*/
void w1_9to4(DOUBLE   tensor[3][3], 
             DOUBLE  *vector);
/*----------------------------------------------------------------------*
 | reduce 4-stufigen tensor to 2-stufigen tensor           he    04/03   |
 *----------------------------------------------------------------------*/
void w1_81to16(DOUBLE   tensor4[3][3][3][3], 
               DOUBLE **tensor2);
/*----------------------------------------------------------------------*
 | reduce 4-stufigen tensor to 2-stufigen tensor           he    04/03   |
 *----------------------------------------------------------------------*/
void w1_81to16_1(DOUBLE   tensor4[3][3][3][3], 
                 DOUBLE   tensor2[4][4]);
/*----------------------------------------------------------------------*
 | BESTIMMUNG DER SPANNUNGSKOMPONENTEN sig                he    04/03   |
 *----------------------------------------------------------------------*/
void w1_cond(DOUBLE sig[4], 
             DOUBLE **d);
/*----------------------------------------------------------------------*/
void w1_strain_energy(ELEMENT *ele,                   /* actual element */
                      DOUBLE **stress,   /* 2PK str. at act. Gauss point*/
		      DOUBLE *strain,   /* Strains at actual Gauss point*/
		      DOUBLE  fac);        /* Int. factor at Gauss point*/		      		                                                                       
/*----------------------------------------------------------------------*/
/*  w1_kinetic_energy.c                                                 */
/*  Calculation of kinetic energy of an element                         */
/*----------------------------------------------------------------------*/                                                           
void w1_kinetic_energy(ELEMENT *ele,                  /* actual element */
                       DOUBLE **mass);                   /* mass matrix */                    		      		                                                                       
/*----------------------------------------------------------------------*/		      
/*  w1_update_history.c                                                 */
/*  Update of state variables for E-M Int. Scheme                       */
/*----------------------------------------------------------------------*/                                                           
void w1_update_history(ELEMENT *ele,                   /* actual element*/
                       W1_DATA *data,                      /* wall1 data*/
	               MATERIAL *mat);                 /*actual material*/
/*----------------------------------------------------------------------*/
/*  w1_history.c                                                        */
/*  Update of state variables for E-M Int. Scheme                       */
/* (called by w1_update_history)                                        */
/*----------------------------------------------------------------------*/                                                           
void w1_history(ELEMENT *ele,                           /*actual element*/
           DOUBLE **b_bar,                              /*b_bar operator*/
	   DOUBLE **boplin,                            /*boplin operator*/
	   DOUBLE **stress,              /* 2PK str. at act. Gauss point*/
	   DOUBLE  *F,                     /*Deformation gradient tensor*/
	   INT      numeps,               /* number of strain components*/
	   INT      nd,       /* total number degrees of freedom of ele.*/
	   INT      ip);                     /*Integration point counter*/	   
/*----------------------------------------------------------------------*
 |  w1_restart.c                                            sh 03/04    |
 *----------------------------------------------------------------------*/
void w1_write_restart(ELEMENT *actele, MATERIAL  *mat, INT nhandle, long int *handles, INT init);
void w1_read_restart(ELEMENT *actele, MATERIAL  *mat, long int *handles, INT init);
/*----------------------------------------------------------------------*/                                             
/*----------------------------------------------------------------------*/
/*  w1_tri_service.c                                                    */
/*----------------------------------------------------------------------*/                                                           
void w1_degrectri(DOUBLE     *funct, 
                  DOUBLE    **deriv, 
                  DOUBLE      r, 
                  DIS_TYP     typ,
                  INT         option);
void w1_edgejaco(ELEMENT    *ele, 
                 DOUBLE     *funct,    
                 DOUBLE    **deriv,   
                 DOUBLE    **xjm,     
                 DOUBLE     *det,          
                 INT         iel,
                 INT        *iedgnod);
/*----------------------------------------------------------------------*/                                             
#endif /*D_WALL1*/
/*! @} (documentation module close)*/
                                                                        
