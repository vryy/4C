/*!----------------------------------------------------------------------
\file
\brief contains all prototypes for wall1 element

*----------------------------------------------------------------------*/
#ifdef D_WALL1

/*! 
\addtogroup WALL1
*//*! @{ (documentation module open)*/

/*-----------------------------------------------------------------------*
|  w1_bop.c                                                  al 9/01     |
|  calculate operator matrix at point r,s                                |
*-----------------------------------------------------------------------*/
void w1_bop(double    **bop,   /* derivative operator                   */
            double    **deriv, /* shape functions derivatives           */
            double    **xjm,   /* jacobian matrix                       */
            double      det,   /* determinant of jacobian               */
            int         iel);  /* number of nodes at actual element     */
/*-----------------------------------------------------------------------*
|  w1_cal_deriv.c                                            al 9/01     |
|  compute displacement derivatives                                      |
*-----------------------------------------------------------------------*/
void w1_disd(ELEMENT   *ele,  /* actual element                         */
             double   **bop,  /* derivative operator                    */
             double    *gop,  /* additional derivative operator         */
             double    *alpha,/*internal dof                            */
             WALL_TYPE wtype, /* plane stress/strain...                 */
             double    *disd);/* displacement derivatives               */
/*-----------------------------------------------------------------------*
|  w1_cal_eps.c                                              al 9/01     |
|  evaluates linear/nonlinear strains from                               |
|  displacement derivatives                                              |
*-----------------------------------------------------------------------*/
void w1_eps( double   *disd, /* displacement derivatives                */
             WALL_TYPE wtype,/* plane stress/strain...                  */
             double   *eps); /* strains                                 */
/*-----------------------------------------------------------------------*
|  w1_cal_stress.c                                           al 9/01     |
|  evaluate element stresses                                             |
|  2-D isoparametric degenerated element                                 |
*-----------------------------------------------------------------------*/
void w1_cal_stress(ELEMENT   *ele, /* actual element                    */
                   W1_DATA   *data,/* wall1 data                        */
                   MATERIAL  *mat, /* actual material                   */
                   ARRAY     *estif_global, /* element stiffness matrix */ 
                   double    *force,  /* glob. vector - internal forces */
                   int        kstep,  /* number of step      fh 06/02   */
                   int        init);  /* initialize this function       */
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
void w1_mami(double *stress, /* vector of stresses                      */
             double *fps,    /* first  principal stress                 */ 
             double *sps,    /* second principal stress                 */ 
             double *aps);   /* angle of principal direction            */ 
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
double w1rsn (int node, /* number of actual integration point           */
              int  irs, /* r/s identifier                               */
              int  iel);/* number of nodes at actual element            */
/*-----------------------------------------------------------------------*
|  w1_cal_stress.c                                           al 9/01     |
|       extrapolation form gauss points for rectangles                   |
*-----------------------------------------------------------------------*/
void w1recs(double *funval,/* function value                            */
            double  r,     /* r/s coordinates                           */
            double  s,     /* r/s coordinates                           */
            double *fval,  /* function values at gauss points           */
            double *fpar,  /* function parameters                       */
            int     igauss,/* number of gauss points                    */
            int     icode);/* ==1 initialize function parameters        */
                           /* > 1            function evaluation        */      
/*-----------------------------------------------------------------------*
|  w1_call_mat.c                                             al 9/01     |
|  select proper material law                                            |
*-----------------------------------------------------------------------*/
void w1_call_mat(ELEMENT   *ele, /* actual element                      */
                 MATERIAL  *mat, /* actual material                     */
                 WALL_TYPE wtype,/* plane stress/strain...              */  
                 double **bop,   /* derivative operator                 */
                 double  *gop,   /* additional derivative operator      */
                 double  *alpha, /* internal dof                        */
                 double **xjm,   /* jacobian matrix                     */
                 int ip,         /* integration point Id                */
                 double *stress, /* vector of stresses                  */
                 double **d,     /* constitutive matrix                 */
                 int istore,     /* controls storing of new stresses    */
                 int newval);    /* controls evaluation of new stresses */   
/*-----------------------------------------------------------------------*
|  w1_call_stiff.c                                           al 9/01     |
|  usual stiffness matrix total lagrangian formulation                   |
*-----------------------------------------------------------------------*/
void w1_keku(double **s,    /* element stiffness matrix                 */
             double **bs,   /* derivative operator                      */
             double **d,    /* constitutive matrix                      */
             double   fac,  /* integration factor                       */
             int      nd,   /* total number degrees of freedom of ele.  */
             int      neps);/* actual number of strain components       */
/*-----------------------------------------------------------------------*
|  w1_funcderiv.c                                            al 9/01     |
|  shape functions and derivatives                                       |
*-----------------------------------------------------------------------*/
void w1_funct_deriv(double     *funct,  /* shape function values        */
                    double    **deriv,  /* shape function derivatives   */
                    double      r,      /* r/s coordinates              */
                    double      s,      /* r/s coordinates              */
                    DIS_TYP     typ,    /* quad4, quad8, quad9 ...      */
                    int         option);/* index for function evaluation*/ 
void w1_degfuncderiv(double     *funct, 
                      double    **deriv, 
                      double      r, 
                      DIS_TYP     typ,
                      int         option);
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
            int        option);/* index for function evaluation         */
/*-----------------------------------------------------------------------*
|  w1_jaco.c                                                 al 9/01     |
|  calculate operator matrix at point r,s                                |
*-----------------------------------------------------------------------*/
void w1_jaco(double     *funct,/* shape function values                 */
             double    **deriv,/* shape function derivatives            */
             double    **xjm,  /* jacobian matrix                       */
             double     *det,  /* determinant of jacobian               */
             ELEMENT    *ele,  /* actual element                        */
             int         iel); /* number of nodes at actual element     */
/*-----------------------------------------------------------------------*
|  w1_mat_linel.c                                            al 9/01     |
|  constitutive matrix - linear elastic - 2D                             |
|  plane stress, plane strain, rotational symmetry                       |
*-----------------------------------------------------------------------*/
void w1_mat_linel(double ym,      /* young's modulus                    */     
                  double pv,      /* poisson's ratio                    */     
                  WALL_TYPE wtype,/* plane stress/strain...             */   
                  double **d);    /* constitutive matrix                */   
/*-----------------------------------------------------------------------*
|  w1_mat_plast_dp.c                                         al 9/01     |
|  constitutive matrix - forces - Drucker Prager - 2D                    |
|  plane stress, plane strain, rotational symmetry                       |
*-----------------------------------------------------------------------*/
void w1_mat_plast_dp(double ym,      /* young's modulus                 */     
                     double pv,      /* poisson's ratio                 */     
                     double ALFAT,   /* temperature expansion factor    */  
                     double sigy,    /* yield stresse                   */    
                     double hard,    /* hardening modulus               */    
                     double phi,     /* angle - drucker prager          */    
                     ELEMENT   *ele, /* actual element                  */ 
                     WALL_TYPE wtype,/* plane stress/strain...          */      
                     double **bop,   /* derivative operator             */  
                     double  *gop,   /* add. derivative operator        */  
                     double  *alpha, /* internal dof                    */  
                     int ip,         /* integration point Id            */         
                     double *stress, /* vector of stresses              */
                     double **d,     /* constitutive matrix             */    
                     int istore,     /* controls storing of new stresses*/    
                     int newval);    /* controls evaluation of stresses */ 
/*-----------------------------------------------------------------------*
|  w1_mat_plast_epc.c                                        al 9/01     |
|  const. matrix - forces - elastoplastic concrete - 2D                  |
|  plane stress, plane strain                                            |
*-----------------------------------------------------------------------*/
void w1_mat_plast_epc(double dens      , /* density                     */
                      double emod      , /* young's modulus             */
                      double vnu       , /* poisson's ratio             */
                      double alfat     , /* temperature expansion factor*/ 
                      double xsi       ,                                
                      double sigyc     , /* yield stresse               */
                      double ftm       , /* tensile strenght            */                    
                      double fcm       , /* compressive strenght        */        
                      double gt        , /* tensile fracture energy     */        
                      double gc        , /* compression fracture energy */        
                      double gamma1    , /* fitting factor yield funct 1*/     
                      double gamma2    , /* fitting factor yield funct 2*/     
                      int    nstiff    , /* flag for tension stiffening */
                      int    maxreb    , /* max. number of rebars       */             
                      int    *rebar     ,/* rebar id                    */
                      double *reb_area  ,/* rebar area                  */
                      double *reb_ang   ,/* rebar angle                 */
                      double *reb_so    ,/* minimum bond length         */
                      double *reb_ds    ,/* diameter                    */
                      double *reb_rgamma,/* 4=deformed bars 2=plane bar */
                      double *reb_dens  ,/* density                     */
                      double *reb_alfat ,/* temperature expansion factor*/ 
                      double *reb_emod  ,/* young's modulus             */
                      double *reb_rebnue,/* poisson's ratio             */
                      double *reb_sigy  ,/* yield stresse               */
                      double *reb_hard  ,/* hardening modulus           */        
                      ELEMENT   *ele,    /* actual element              */
                      WALL_TYPE wtype,   /* plane stress/strain...      */          
                      double **bop,      /* derivative operator         */
                      double  *gop,      /* add derivative operator     */
                      double  *alpha,    /* internal dof                */
                      double **xjm,      /* jacobian matrix             */
                      int ip,            /* integration point Id        */
                      double *stressc,   /* vector of stresses          */
                      double **d,        /* constitutive matrix         */
                      int istore,        /* controls storing of stresses*/
                      int newval);       /* controls eval. of stresses  */  
/*-----------------------------------------------------------------------|
| w1_mat_plast_epc_serv.c                                     al 9/01    |
|      initialize d-components d[0][3], d[1][3], d[2][3], d[3][3]        |
|      on elastic                                                        |
|-----------------------------------------------------------------------*/
void w1iwadi (double ym,   /* young's modulus                           */                           
              double pv,   /* poisson's ratio                           */ 
              double *di); /* components of constitutive matrix         */
/*-----------------------------------------------------------------------|
| w1_mat_plast_epc_serv.c                                     al 9/01    |
|      topic : wall1 - yield criterion for plane strain                  |
|              yield - criterion:  drucker - prager                      |
|-----------------------------------------------------------------------*/
void w1yicsr (double dev,      /* norm der deviatorischen spannungen    */
              double hyd,      /* 1. invariante                         */
              double sigym,    /* uniaxial yield stress                 */
              double alpha,    /* neigungswinkel im invariantenraum     */
              double *ft);     /* yield condition                       */ 
/*-----------------------------------------------------------------------|
| w1_mat_plast_epc_serv.c                                     al 9/01    |
|      topic : wall1 - yield criterion (cap - region)                    |
|              yield - criterion:  cap - model                           |
|-----------------------------------------------------------------------*/
void w1yiccap (double dev,   /* norm of deviatoric predictor stresses   */
               double hyd,   /* 1. invariant  of the predictor stresses */
               double *sigym,/* yield stresses                          */
               double alpha, /* factor for the first invariants         */
               double *ft);  /* yield function                          */
/*-----------------------------------------------------------------------|
| w1_mat_plast_epc_serv.c                                     al 9/01    |
|    topic: calculate the gradients of the elastic predictor             |              
|           deviatoric/hydrostatic stresses                              |
|-----------------------------------------------------------------------*/
void w1preds(double *sigma, /* elastic predictor projected              */
             double *alpha, /* neigungswinkel im invariantenraum        */
             double *devsig,/* deviatoric stresses                      */
             double *dev,   /* norm der deviatorischen spannungen       */
             double *dn,    /* gradient in deviatoric direction         */
             double *dcom,  /* gradient in hdrostatic direction         */
             double *grad); /* total gradient                           */
/*-----------------------------------------------------------------------|
| w1_mat_plast_epc_serv.c                                     al 9/01    |
|    topic: calculate the deviatoric/hydrostatic stress components       |
|           of the elastic predictor                                     |
|-----------------------------------------------------------------------*/
void w1pres (double *sigma ,    /*  elastic predictor projected         */
             double *devsig,    /*  deviatoric stresses                 */
             double *sm    ,    /*  hydrrostatic stresses               */
             double *dev   ,    /*  norm der deviatorischen spannungen  */
             double *hyd);      /*  1. invariante                       */
/*-----------------------------------------------------------------------*
|  w1_mat_plast_mises.c                                      al 9/01     |
|  constitutive matrix - forces - linear elastic- von Mises - 2D a       |
|  plane stress, plane strain, rotational symmetry                       |
*-----------------------------------------------------------------------*/
void w1_mat_plast_mises(double ym,      /* young's modulus              */
                        double pv,      /* poisson's ratio              */
                        double ALFAT,   /* temperature expansion factor */
                        double sigy,    /* yield stresse                */
                        double hard,    /* hardening modulus            */       
                        double gf,      /* fracture energy              */
                        double betah,
                        ELEMENT   *ele, /* actual element               */
                        WALL_TYPE wtype,/* plane stress/strain...       */         
                        double **bop,   /* derivative operator          */
                        double  *gop,   /* add. derivative operator     */
                        double  *alpha, /* internal dof                 */
                        int ip,         /* integration point Id         */
                        double *stress, /* vector of stresses           */
                        double **d,     /* constitutive matrix          */
                        int istore,     /* controls storing of stresses */
                        int newval);    /* controls eval. of stresses   */ 
/*-----------------------------------------------------------------------*
|  w1_mat_plast_serv.c                                       al 9/01     |
|                                                                        |
|    topic : wall1 - yield criterion plane_stress                        |
|            -   -   -- -  --                                            |
|            yield - criterion:  phi > 0 drucker - prager                |
|                                phi = 0 mises                           |
|            for material model 3 - plasticity                           |
*-----------------------------------------------------------------------*/
void w1yilcr(double E,    /* young's modulus                            */
             double Eh,   /* hardening modulus                          */
             double betah,
             double sigy, /* yield stresse                              */
             double epstn,/* equivalent uniaxial plastic strain         */    
             int    isoft,/* softening                                  */
             double dia,  /* equivalent element length                  */
             double *tau, /* current stresses (local)                   */
             double *ft); /* yield condition                            */     
/*-----------------------------------------------------------------------*
|  w1_mat_plast_serv.c                                       al 9/01     |
|  radial return for elements with mises material model                  |
*-----------------------------------------------------------------------*/
void w1radi(double e,        /* young's modulus                         */
            double eh,       /* hardening modulus                       */
            double betah,
            double sigy,     /* yield stresse                           */
            double vnu,      /* poisson's ratio                         */
            double dia,      /* equivalent element length               */
            double *sigma,   /* elastic predicor projected onto yield   */
            double *qn,      /* backstress vector                       */
            int    isoft,    /* softening                               */
            double *epstn,   /* equivalent uniaxial plastic strain      */       
            double *dlam,    /* increment of plastic multiplier         */
            WALL_TYPE wtype);/* plane stress/strain...                  */
/*-----------------------------------------------------------------------*
|  w1_mat_plast_serv.c                                       al 9/01     |
|    forms the elasto-plastic consistent tangent material tensor         |
*-----------------------------------------------------------------------*/
void w1mapl(double e,        /* young's modulus                         */
            double eh,       /* hardening modulus                       */
            double betah,
            double sigy,     /* yield stresse                           */
            double vnu,      /* poisson's ratio                         */
            double dia,      /* equivalent element length               */
            double *tau,     /* current stresses (local)                */
            int    isoft,    /* softening                               */
            double *epstn,   /* equivalent uniaxial plastic strain      */       
            double *dlam,    /* increment of plastic multiplier         */
            double **d,      /* constitutive matrix                     */
            WALL_TYPE wtype);/* plane stress/strain...                  */
/*-----------------------------------------------------------------------*
|  w1_mat_plast_serv.c                                       al 9/01     |
|    topic : wall1 - yield criterion plane_stress                        |
|            -   -   -- -  --                                            |
|            yield - criterion:  phi > 0 drucker - prager                |
|                                phi = 0 mises                           |
*-----------------------------------------------------------------------*/
void w1yilcr_dp(double E,      /* young's modulus                       */
                double Eh,     /* hardening modulus                     */
                double phi,    /* angle - drucker prager                */
                double sigy,   /* yield stresse                         */
                double *sigym, /* uniaxial yield stress                 */
                double epstn,  /* equivalent uniaxial plastic strain    */         
                double *tau,   /* current stresses (local)              */
                double *ft);   /* yield condition                       */          
/*-----------------------------------------------------------------------*
|  w1_mat_plast_serv.c                                       al 9/01     |
|    radial return for elements with mises material model                |
|                                                                        |
*-----------------------------------------------------------------------*/
void w1radi_dp(double e,        /* young's modulus                      */
               double eh,       /* hardening modulus                    */
               double phi,      /* angle - drucker prager               */
               double sigy,     /* yield stresse                        */
               double vnu,      /* poisson's ratio                      */
               double *sigma,   /* elastic predicor projected onto yield*/
               double *qn,      /* backstress vector                    */
               double *epstn,   /* equivalent uniaxial plastic strain   */          
               double *sigym,   /* uniaxial yield stress                */ 
               double *dlam,    /* increment of plastic multiplier      */
               WALL_TYPE wtype);/* plane stress/strain...               */ 
/*-----------------------------------------------------------------------*
|  w1_mat_plast_serv.c                                       al 9/01     |
|     calculate element diameter (equivalent length) for one element     |            
*-----------------------------------------------------------------------*/
void w1cdia(ELEMENT   *ele,    /* actual element                        */
            W1_DATA   *data,   /* wall1 data                            */
            double    *funct_h,/* shape function values                 */
            double   **deriv_h,/* shape function derivatives            */
            double   **xjm_h); /* jacobian matrix                       */
/*-----------------------------------------------------------------------|
|  w1_mat_plast_serv.c                                       al 9/01     |
|     topic : wall1 - concrete material prevalues                        |
|             -   -   -                 -----                            |
|             yield - criterion:  drucker-prager / spherical cap         |
|-----------------------------------------------------------------------*/
void w1cpreva (double *epst,     /* equivalent uniaxial plastic strain  */
               double *sigym,    /* [0] uniaxial tension yield stress   */
                                 /* [1] yield stress "inverted cone"    */
                                 /* [2] uniaxial compr. yield stress    */
               double *alpha,    /* factor for the first invariants     */
	       double *hards,    /* hardening modulus                   */
               double *e,        /* young modulus                       */
               double *g,        /* shear modulus                       */
               double *vnu,      /* poisson's ratio                     */
               double *com,      /* bilk modulus                        */
               double *fcm,      /* compressive strenght                */
               double *gc,       /* compression fracture energy         */
               double *ftm,      /* tensile strenght                    */
               double *gt,       /* tensile fracture energy             */
               double *gamma1,   /* fitting factor yield function 1     */
               double *gamma2,   /* fitting factor yield function 2     */
               double *dfac,     /* damage factor                       */
               double *dia,      /* equivalent element length           */
	       double *acrs,     /* average crack spacing               */
               double *cappaet,  /* max. elastic tension strain         */
               double *cappaut,  /* tensile fracture strain             */
               double *cappae,   /* max. elastic compression strain     */
               double *cappauc,  /* compressive fracture strain         */
               double *sig3,     /* equivalent compressive stress       */
               double *fbd,      /* tension stiffening stress           */
               double **d,       /* elastic material matrix             */
               double *sig,      /* stresses from last update           */
               double *eps);     /* strains from last update            */
/*-----------------------------------------------------------------------|
|  w1_mat_plast_serv.c                                       al 9/01     |
|       topic : radial return for elements with elastoplastic material   |
|               (concrete)                                               |
|-----------------------------------------------------------------------*/
void w1cradi (double *sigma,  /* elastic predictor projected onto yield */ 
              double *epstn,  /* equivalent uniaxial plastic strain     */ 
              double *dlam,   /* incremental plastic multiplier         */ 
              WALL_TYPE wtype,/* plane stress/strain...                 */
              int    yip,     /* stress state   1 =elastic >=2 =plastic */ 
              double *alpha,  /* factor for the first invariants        */ 
              double *ft,     /* yield condition                        */ 
              double *e,      /* young modulus                          */ 
              double *g,      /* shear modulus                          */ 
	      double *com,    /* bulk modulus                           */ 
              double *sigym,  /* uniaxial predictor yield stress        */ 
              double *hards,  /* plastic modulus                        */ 
              double *sigy,   /* actual uniaxial yield stress           */ 
              double *dn,     /* gradient components in deviatoric dir  */ 
              double *dcom,   /* gradient components in hdrostatic dir  */ 
              double *grad,   /* total gradient                         */ 
              double *devsig, /* deviatoric predictor stresses          */ 
              double *sm,     /* hydrostatic predictor stresses         */ 
              double *fcm,    /* compressive strenght                   */ 
              double *gc,     /* compression fracture energy            */ 
              double *ftm,    /* tensile strenght                       */ 
	      double *gt,     /* tensile fracture energy                */ 
              double *gamma1, /* fitting factor yield function 1        */ 
              double *gamma2, /* fitting factor yield function 2        */ 
              double *dia,    /* equivalent element length              */ 
              double *acrs);  /* average crack spacing                  */ 
/*-----------------------------------------------------------------------*
|  w1_mat_plast_serv.c                                       al 9/01     |
|                                                                        |
|    forms the elasto-plastic consistent tangent material tensor         |
*-----------------------------------------------------------------------*/
void w1mapl2(double *tau,      /* current stresses (local)              */
             double **d,       /* material matrix to be calculated      */
             double *dlam,     /* increment of plastic multiplier       */
             WALL_TYPE wtype,  /* plane stress/strain...                */
             double *alpha,    /* neigungswinkel der fliessflaechen     */
             double *emod,     /* elastizitaetsmodul                    */
             double *g,        /* schubmodul                            */
             double *com,      /* kompressionsmodul                     */
             double *betah,    /* factor for isotrop/kinemat. hardening */
             double *hards,    /* plastic hardeningmodulus              */
             double *dn,       /* gradient components in dev. direction */
             double *grad,     /* total gradient                        */
             double *dev);     /* norm of the dev. predictor stresses   */
/*-----------------------------------------------------------------------|
|  w1_mat_plast_serv.c                                       al 9/01     |
|       topic : radial return for multisurface problems                  |
|               (material model: concrete)                               |
|-----------------------------------------------------------------------*/
void w1cradms(double *sigma,  /* elastic predictor projected onto yield */ 
              double *epstn,  /* equivalent uniaxial plastic strain     */ 
              double *dlam,   /* incremental plastic multiplier         */ 
              WALL_TYPE wtype,/* plane stress/strain...                 */
              int     yip,    /* stress state   1 =elastic >=2 =plastic */ 
              double *alpha,  /* factor for the first invariants        */ 
              double *ft,     /* yield condition                        */ 
              double  e,      /* young modulus                          */ 
              double  g,      /* shear modulus                          */ 
	      double  com,    /* bulk modulus                           */ 
              double *sigym,  /* uniaxial predictor yield stress        */ 
              double *hards,  /* plastic modulus                        */ 
              double  dn[3][4],  /*gradient components in deviatoric dir*/ 
              double  dcom[3][4],/*gradient components in hdrostatic dir*/ 
              double *devsig, /* deviatoric predictor stresses          */ 
              double *sm,     /* hydrostatic predictor stresses         */ 
              double  fcm,    /* compressive strenght                   */ 
              double  gc,     /* compression fracture energy            */ 
              double  ftm,    /* tensile strenght                       */ 
	      double  gt,     /* tensile fracture energy                */ 
              double  gamma1, /* fitting factor yield function 1        */ 
              double  gamma2, /* fitting factor yield function 2        */ 
              double  dia,    /* equivalent element length              */ 
              double  acrs);  /* average crack spacing                  */ 
/*-----------------------------------------------------------------------|
|  w1_mat_plast_serv.c                                       al 9/01     |
|     topic: convert some arrays                                         |
|-----------------------------------------------------------------------*/
      void w1conver(double *dlam,       /*  incr. of plastic multiplier */             
                    double *alpha,      /*  factor for the first inv.   */             
                    double *hards,      /* hardening modulus            */             
                    double dn[3][4],    /*  gradient comp. in dev. dir. */  
                    double grad[3][4],  /*  total gradient              */             
                    double *dlamc,      /*  -|                          */             
                    double *alphac,     /*   |                          */             
                    double *hardsc,     /*   |-  converted arrays       */             
                    double dnc[2][4],   /*   |                          */             
                    double gradc[2][4]);/*  -|                          */             
/*---------------------------------------------------------------------- |
|  w1_mat_plast_serv.c                                       al 9/01     |
|    forms the elasto-plastic consistent tangent material tensor         |
|    for wall element   (drucker - prager)                               |
|    (generalized for the multisurface problem: including the apex)      |
*-----------------------------------------------------------------------*/
void w1maplg(double *tau,      /* current stresses (local)              */
             double **d,       /* material matrix to be calculated      */
             double  dlam,     /* increment of plastic multiplier       */
             WALL_TYPE wtype,  /* plane stress/strain...                */
             int     yip,      /* stress state  1=elastic 2=plastic     */
             double  emod,     /* elastizitaetsmodul                    */
             double  g,        /* schubmodul                            */
             double  com,      /* kompressionsmodul                     */
             double *hards,    /* plastic hardeningmodulus              */
             double  dn[2][4], /* gradient components in dev. direction */
             double  grad[2][4]);/*norm of the dev. predictor stresses  */
/*-----------------------------------------------------------------------|
|  w1_mat_plast_serv.c                                        al 9/01    |
|       topic : radial return for elements with elastoplastic material   |
|               (concrete)      CAP-REGION                               |
|-----------------------------------------------------------------------*/
void w1radcap(double *sigma,  /* elastic predictor projected onto yield */ 
              double *epstn,  /* equivalent uniaxial plastic strain     */ 
              double *dlam,   /* incremental plastic multiplier         */ 
              WALL_TYPE wtype,/* plane stress/strain...                 */
              double *alpha,  /* factor for the first invariants        */ 
              double *e,      /* young modulus                          */ 
              double *g,      /* shear modulus                          */ 
	      double *com,    /* bulk modulus                           */ 
              double *sigym,  /* uniaxial predictor yield stress        */ 
              double *hards,  /* plastic modulus                        */ 
              double *grad,   /* total gradient                         */ 
              double *devsig, /* deviatoric predictor stresses          */ 
              double *sm,     /* hydrostatic predictor stresses         */ 
              double *dev,    /* norm of the dev. predictor stresses    */ 
              double *hyd,    /* 1st invariant  of the predictor stress.*/ 
              double *hydn,   /* 1st invariant  of the new stresses     */ 
              double *fcm,    /* compressive strenght                   */ 
              double *gc,     /* compression fracture energy            */ 
              double *ftm,    /* tensile strenght                       */ 
	      double *gt,     /* tensile fracture energy                */ 
              double *gamma2, /* fitting factor yield function 2        */ 
              double *dia);   /* equivalent element length              */ 
/*-----------------------------------------------------------------------|
|  w1_mat_plast_serv.c                                        al 9/01    |
|       topic: forms the  e l a s t o - p l a s t i c                    |
|              c o n s i s t e n t  tangent material tensor              |
|              for wall element   (drucker - prager)                     |
|-----------------------------------------------------------------------*/
void w1mplcap(double *tau,      /* current stresses (local)             */
              double **d,       /* material matrix to be calculated     */
              double *dlam,     /* increment of plastic multiplier      */
              WALL_TYPE wtype,  /* plane stress/strain...               */
              double *alpha,    /* factor for the first invariants      */
              double *emod,     /* young modulus                        */
              double *vnu,      /* poisson's ratio                      */
              double *hards,    /* plastic hardeningmodulus             */
              double *sigym,    /* [0] uniaxial tension yield stress    */
              double *grad,     /* total gradient                       */
              double *cappae,   /* max. elastic compression strain      */
              double *cappauc,  /* compressive fracture strain          */
              double *epst,     /* equivalent uniaxial plastic strain   */
              double *sig3);    /* equivalent compressive stress        */
/*-----------------------------------------------------------------------|
|  w1_mat_plast_serv.c                                        al 9/01    |
| compute average crack spacing for element wall1     (model: concrete)  |
|-----------------------------------------------------------------------*/
void w1acrs(double  hte,       /* section properties                    */   
               int  maxreb,    /* max. number of rebars                 */   
            double *stress,    /* element stresses                      */   
            double *angle,     /* angle of the principal concrete stress*/   
            double *reb_area,  /* rebar area                            */
            double *reb_ang,   /* rebar angle                           */
            double *reb_so,    /* minimum bond length                   */
            double *reb_ds,    /* diameter                              */
            double *reb_rgamma,/* 4 := deformed bars = 2 := plane bar   */
            double *thick,     /* thickness ot the structure            */   
            double dia,        /* equivalent element length             */   
            double **xjm,      /* jacobian matrix                       */
            double *acrs);     /* average crack spacing                 */   
/*-----------------------------------------------------------------------|
|  w1_mat_plast_serv.c                                        al 9/01    |
|    topic: calculate the gradients of the elastic predictor             |              
|           deviatoric/hydrostatic stresses                              |
|-----------------------------------------------------------------------*/
void w1cpreds(double *sigym, /* uniaxial predictor yield stress         */
              double *alpha, /* factor for the first invariants         */
              double *devsig,/* deviatoric predictor stresses           */
              double  hyd,   /* 1. invariant of the predictor stresses  */
              double *grad); /* total gradient                          */
/*-----------------------------------------------------------------------*
|  w1_mat_plast_serv.c                                       al 9/01     |
|    WALL1 : array = value                                               |
*-----------------------------------------------------------------------*/
void w1_init44(double a[4][4],double val);
/*-----------------------------------------------------------------------*
|  w1_mat_plast_serv.c                                       al 9/01     |
|    WALL1 : R(I,J) = A(I,K)*B(K,J) ---  R = A*B                         |
*-----------------------------------------------------------------------*/
void w1_AxB_444(double a[4][4],double b[4][4],double r[4][4]);
/*-----------------------------------------------------------------------*
|  w1_mat_plast_serv.c                                       al 9/01     |
|    WALL1 : R(I,J) = A(I,K)*B(K,J) ---  R = A*B                         |
*-----------------------------------------------------------------------*/
void w1_AxBT_444(double a[4][4],double b[4][4],double r[4][4]);
/*-----------------------------------------------------------------------*
|  w1_mat_plast_serv.c                                       al 9/01     |
|    WALL1 : R(I,J) = A(I,K)*B(K,J) ---  R = A*B                         |
*-----------------------------------------------------------------------*/
void w1_AxB_441(double a[4][4],double b[4],double r[4]);
/*-----------------------------------------------------------------------*
|  w1_mat_plast_serv.c                                       al 9/01     |
|    WALL1 : R(I,J) = A(I,K)*B(K,J) ---  R = A*B                         |
*-----------------------------------------------------------------------*/
void w1_AxBT_414(double a[4],double b[4],double r[4][4]);
/*-----------------------------------------------------------------------*
|  w1_mat_plast_serv.c                                       al 9/01     |
|    WALL1 : R(I,J) = A(I,K)*B(K,J) ---  R = A*B                         |
*-----------------------------------------------------------------------*/
void w1_AxFAC_44(double a[4][4], double r[4][4], double fac);
/*-----------------------------------------------------------------------*
|  w1_mat_plast_serv.c                                       al 9/01     |
|    WALL1 : ADDITION   R = A+B   (FULL MATRICES)                        |
*-----------------------------------------------------------------------*/
void w1_AaddB_44(double a[4][4], double b[4][4], double r[4][4]);
/*-----------------------------------------------------------------------*
|  w1_mat_plast_serv.c                                       al 9/01     |
|  constitutive matrix - forces - plastic - rebar                al 9/01 |
*-----------------------------------------------------------------------*/
void w1_mat_rebar(ELEMENT   *ele, /* actual element                     */
                  MATERIAL  *mat, /* actual material                    */
                  double   **bop, /* derivative operator                */
                  double    *gop, /* add. derivative operator           */
                  double    *alpha,/* internal dof                      */
                  double   **xjm, /* jacobian matrix                    */
                  double *stress, /* vector of stresses                 */
                  double     **d, /* constitutive matrix                */
                  int         ip, /* integration point Id               */
                  int       lanz, /* rebar id                           */
                  int     istore);/* controls storing of new stress.2wa */
/*-----------------------------------------------------------------------|
|  w1_mat_serv.c                                              al 9/01    |
|      topic: condensed stress vector                                    |
|             plane strain --> plane stress                              |
|-----------------------------------------------------------------------*/
void w1consig (double **d,/* current material matrix components d14-d44 */
               double  *sigma,  /* current stresses (local)             */
               double  *sigmac);/* condensed stress vector              */
/*-----------------------------------------------------------------------|
|  w1_mat_serv.c                                              al 9/01    |
|      topic : condense the constitutive tensor                          |
|              plane strain --> plane stress                             |
|-----------------------------------------------------------------------*/
void w1concep (double **d);  /* material matrix to be calculated        */
/*-----------------------------------------------------------------------|
|  w1_mat_serv.c                                              al 9/01    |
|    topic: evaluate the incremental strain in thickness direction       |              
|           for wall element (plane stress)                              |
|-----------------------------------------------------------------------*/
void w1de33(double *sigi,  /* stresses from last iteration step         */
            double *epsi,  /* strains from last iteration step          */
            double *di,    /* components d41,d42,d43,d44 of the         */
                           /* const. tensor from last iteration step    */
            double *strain);/* current strains (local)                  */
/*-----------------------------------------------------------------------*
|  w1_static_ke.c                                            al 9/01     |
|  integration of linear stiffness ke for wall1 element                  |
*-----------------------------------------------------------------------*/
void w1static_ke(ELEMENT   *ele,    /* actual element                   */
                 W1_DATA   *data,   /* wall1 data                       */
                 MATERIAL  *mat,    /* actual material                  */
                 ARRAY     *estif_global, /* element stiffness matrix   */
                 ARRAY     *emass_global, /* element mass matrix        */
                 double    *force,  /* global vector for internal forces*/  
                 int        init);  /* initialize this function         */
/*-----------------------------------------------------------------------*
|  w1_static_keug.c                                        ah 06/02      |
|  integration of nonlinear stiffness keug for wall1 element             |
*-----------------------------------------------------------------------*/
void w1static_keug(ELEMENT  *ele,           /* actual element           */
                   W1_DATA   *data,         /* wall1 data               */
                   MATERIAL  *mat,          /* actual material          */
                   ARRAY     *estif_global, /* element stiffness matrix */
                   ARRAY     *emass_global, /* element mass matrix      */
                   double    *force,/* global vector for internal forces*/
                   int        init);        /* initialize this function */
/*-----------------------------------------------------------------------*
|  w1_static_ke.c                                            al 9/01     |
|  evaluates element forces                                              |
*-----------------------------------------------------------------------*/
void w1fi( double  *F,    /* element stresses                           */
           double   fac,  /* integration factor                         */
           double **bop,  /* derivative operator                        */
           int      nd,   /* total number degrees of freedom of element */
           double  *fie); /* nodal forces                               */           
/*----------------------------------------------------------------------*/
/*  w1_call_mat.c                                         ah 06/02      */
/* get density out of material law                                      */
/*----------------------------------------------------------------------*/
void w1_getdensity(MATERIAL  *mat,        /* actual material            */
                   double *density);      /* density of actual material */
/*----------------------------------------------------------------------*/
/*  w1_boplin.c                                           ah 06/02      */
/* evaluation of linear B-operator                                      */
/*----------------------------------------------------------------------*/
void w1_boplin(double    **boplin,        /* Blin                       */
               double    **deriv,     /* derivatives of ansatzfunctions */
               double    **xjm,           /* jacobian matrix            */
               double      det,           /* det of jacobian matrix     */
               int         iel);          /* nodenumber of element      */
/*----------------------------------------------------------------------*/
/*  w1_defgrad.c                                          ah 06/02      */
/*  evaluation of deformation gradient F                                */
/*----------------------------------------------------------------------*/
void w1_defgrad(double    *F,        /* Deformation gradient            */
                double    *strain,   /* Green-Lagrange-strains          */
                double    **xrefe,   /* coordinates in referenz-config. */
                double    **xcure,   /* coordinates in current-config.  */
                double    **boplin,  /* Blin                            */
                int         iel);    /* nodenumber of element           */
/*----------------------------------------------------------------------*
|  w1_call_matgeononl.c                                     ah 06/02     |
|  select proper material law for large deformations                     |
*-----------------------------------------------------------------------*/
void w1_call_matgeononl(ELEMENT   *ele,   /* actual element             */
                        MATERIAL  *mat,   /* material of actual element */
                        WALL_TYPE wtype,/*Info about actual wall element*/
                        double **boplin,  /* linear B-operator          */
                        double **xjm,     /* jacobian matrix            */
                        int ip,           /* integration point Id       */
                        double *strain,  /*strain at act. gaussian point*/
                        double **stress, /*stress at act. gaussian point*/
                        double **d,       /* material tangent           */
                        int istore,
                        int numeps);     /* number of strain components */
/*----------------------------------------------------------------------*
|  w1_mat_linelgeonon.c                                     ah 06/02     |
|  linear elastic material law for large deformations                    |
*-----------------------------------------------------------------------*/
void w1_mat_linelgeonon(double ym,       /* Young's modulus             */
                        double pv,       /* poisson's ratio             */
                        WALL_TYPE wtype,/*Info about actual wall element*/
                        double *strain,  /*strain at act. gaussian point*/
                        double **d,      /* material tangent            */
                        double **stress, /*stress at act. gaussian point*/
                        int numeps);     /* number of strain components */
/*-----------------------------------------------------------------------*
|  w1_cal_kg.c                                              ah 06/02     |
|  evaluation of geometric part of stiffness matrix in geononl. case     |
*-----------------------------------------------------------------------*/
void w1_kg(double  **kg,                 /* geometric stiffness matrix  */
           double  **boplin,             /* linear B-operator           */
           double  **stress,             /*stress at act. gaussian point*/
           double    fac,                /* integration factor          */
           int       nd,                 /* dof's of element            */
           int       neps);              /* number of strain components */
/*-----------------------------------------------------------------------*
|  w1_cal_keu.c                                              ah 06/02    |
|  evaluation of elast+init. disp. part of stiffness in geononl. case    |
*-----------------------------------------------------------------------*/
void w1_keu(double  **keu,    /* elastic + initial deformation stiffness*/
            double  **boplin,            /* linear B-operator           */
            double  **D,                 /* material tangente           */
            double   *F,                 /* deformation gradient        */
            double    fac,               /* integration factor          */
            int       nd,                /* dof's of element            */
            int       neps);             /* number of strain components */
 /*----------------------------------------------------------------------*
 |  w1_cal_fint.c                                              ah 06/02   |
 | evaluate internal element forces for large def (total Lagr)           |
 *----------------------------------------------------------------------*/
void w1_fint( double **stress,           /* 2.PK stresses               */ 
              double  *F,                /* Deformation gradient        */ 
              double **boplin,           /* B-lin-operator              */ 
              double  *fint,             /* internal forces             */ 
              double   fac,              /* detJ*wr*ws*thickness        */ 
              int      nd);              /* Element-DOF                 */
/*----------------------------------------------------------------------*
 | Transform stress and strain local-global                  fh 7/02    |
 | Local 3-direction is zero                                            |
 *----------------------------------------------------------------------*/ 
void w1_lss(double    *a,    /* vector to be transformed                */
            double    **G,   /* elements of transformation matrix       */
            double    **GI,  /* inverse of G                            */
	    int         it); /* flag for local/global strains/stresses  */                                                                                                                                              
/*----------------------------------------------------------------------*
 | Set Transformation Matrices G and GI                      fh 7/02    |
 *----------------------------------------------------------------------*/
void w1_sett(double   **A,   /* matrix of direction cosines             */
            double    **B,   /* G-Matrix     (L=1)                      */
            double    **C);  /* inverse of G (L=2)                      */                                                         
/*----------------------------------------------------------------------*
 | Calculate Transformation Matrices G and G(Inv)            fh 7/02    |
 *----------------------------------------------------------------------*/
void w1_tram(double   **xjm, /* Elements of Jacobian Matrix             */
            double    **G,   /* G-Matrix                                */
            double    **GI,  /* inverse of G                            */
	    double    **dum);/* matrix of direction cosines             */                                                                     
/*----------------------------------------------------------------------*
 | w1_cal_fext.c                                             ah 07/02   |
 | element line and surface loads                                       |
 *----------------------------------------------------------------------*/
void w1_eleload(ELEMENT  *ele,     /* actual element                  */
                W1_DATA  *data,    /* wall1- Data                     */
                double	 *loadvec, /* global element load vector fext */
                int	  init,    /* flag if init or calculation     */
                int       imyrank);        
void w1_fsiload(ELEMENT  *ele,     /* actual element                  */
                W1_DATA  *data,    /* wall1- Data                     */
                double	 *loadvec, /* global element load vector fext */
                int	  init,    /* flag if init or calculation     */
                int       imyrank);
void w1_iedg(int *iegnod, ELEMENT *ele, int line, int init);		
/*----------------------------------------------------------------------*
 | w1_cal_fext.c                                             ah 07/02   |
 | integration of element loads                                         |
 *----------------------------------------------------------------------*/
void w1_fextsurf(ELEMENT    *ele,        /* actuell element             */
                 double    **eload,      /* element load vector         */
                 double     *funct,      /* ansatz-functions            */
                 double      fac,        /* integration factor          */
                 int         iel);       /* element DOF                 */
/*----------------------------------------------------------------------*
 | w1_bop.c                                                  ah 9/02    |
 | additional operator matrix gop at r,s for incomp modes               |
 *----------------------------------------------------------------------*/
void w1_gop(double    *gop,            /* operator matrix G             */
            double    **xjm0,          /* jacobian matrix at r,s=0      */
            double      det0,          /* det J at r,s=0                */
            double      det,           /* det J at actual r,s           */
            double      e1,            /* actual GP coordinate r        */
            double      e2);           /* actual GP coordinate s        */
/*----------------------------------------------------------------------*
 | w1_incompmode.c                                              ah 9/02 |
 | stiffness matrix due to incompatible modes                           |
 *----------------------------------------------------------------------*/
void w1_knc(double  **knc,             /* stiffness knc= BT C G         */
            double  **bop,             /* operator matrix               */
            double   *gop,             /* additional opperator matrix   */
            double  **d,               /* constitutive matrix           */
            double    fac);            /* integration factor            */
/*----------------------------------------------------------------------*
 | w1_incompmode.c                                              ah 9/02 |
 | stiffness matrix  knn      (4*4)                                     |
 *----------------------------------------------------------------------*/
void w1_knn(double  **knn,             /* stiffness knn= GT C G         */
            double   *gop,             /* additional opperator matrix   */
            double  **d,               /* constitutive matrix           */
            double    fac);            /* integration factor            */
/*----------------------------------------------------------------------*
 | w1_incompmode.c                                              ah 9/02 |
 | evaluate internal element forces due to incomp. modes                |
 *----------------------------------------------------------------------*/
void w1_fintn(double  *F,              /* stress                        */
              double   fac,            /* integration factor            */
              double  *gop,            /* additional opperator matrix   */
              double  *fintn);          /* int forces due to inc modes  */
/*----------------------------------------------------------------------*
 | w1_math.c                                                    ah 9/02 |
 | calculate inverse of arbitrary NxN matrix                            |
 *----------------------------------------------------------------------*/
void w1_inverse_matrix
 (int N,          /* I: size of matrix: N x N                           */
  double **A,     /* I: Original Matrix                                 */
  double **Y);    /* O: Inverse Matrix                                  */
/*----------------------------------------------------------------------*
 | w1_math.c                                                    ah 9/02 |
 | LU-Decomposition of  NxN - matrix ->"numerical recipes in C"         |
 *----------------------------------------------------------------------*/
void w1_ludcmp(double **a,  /* I: Original Matrix (will be rearanged)   */
               int n,       /* I: size of matrix: N x N                 */
               int *indx,   /* O:row premutation pf pivoting            */
               double *d);/* O:+1/-1,depending on even or odd row permut*/
/*----------------------------------------------------------------------*
 | w1_math.c                                                    ah 9/02 |
 | solving AX=B, A is LU decomposed  ->"numerical recipes in C"         |
 *----------------------------------------------------------------------*/
void w1_lubksb(double **a,   /* I: LU-decomposed matrix                 */
               int n,        /* I: size of matrix: N x N                */
               int *indx,    /* I:row premutation pf pivoting           */
               double *b);   /* In:B will be changed-> Out:Solution X   */
/*----------------------------------------------------------------------*
 | w1_math.c                                                    ah 9/02 |
 |   allocate a double vector with subscript range v[nl..nh]            |
 *----------------------------------------------------------------------*/
double *w1_vector(int nl, int nh);
/*----------------------------------------------------------------------*
 | w1_math.c                                                    ah 9/02 |
 |   free a double vector allocated with w1_vector()                    |
 *----------------------------------------------------------------------*/
void w1_free_vector(double *v, int nl, int nh);
/*----------------------------------------------------------------------*
 | w1_incompmode.c                                              ah 9/02 |
 | static kondensation:     Kele =   K  - knc*inverse(knn)*kcn,         |
 |                      Fint ele = Fint - knc*inverse(knn)*fintn,       |
 *----------------------------------------------------------------------*/
void  w1_stat_cond(double **knninv,/*I:stiffness inverse(knn) (4x4)     */
                   double **knc,   /*I: knc = BT C G (8x4)              */
                   double **deltak,/*O: knc*inverse(knn)*kcn            */ 
                   double  *fintn, /*I: fintn=GT*singma                 */
                   double  *deltaf,/*O: knc*inverse(knn)*fintn          */ 
                   ELEMENT *ele); /*I: actual element                  */
/*----------------------------------------------------------------------*
 | w1_incompmode.c                                              ah 9/02 |
 | update of internal dof alpha                                         |
 *----------------------------------------------------------------------*/
void  w1_updalpha(double  *alpha,  /*O: internal dof for incomp modes   */
                  ELEMENT *ele,    /*I:actual element                   */
                  double **knc,    /*I: mixed stiffness                 */
                  double **knninv, /*I: inverse of incomp. stiffness    */
                  double  *fintn,  /*I: int. forces of inc. modes       */
                  int      istore);/*I: update after loadstep->istore=1 */
/*----------------------------------------------------------------------*
 |  integration routine for WALL1 element                      al 6/01  |
 *----------------------------------------------------------------------*/
void w1_oint(
             ELEMENT   *ele, 
             W1_DATA   *data, 
             MATERIAL  *mat,
             double    *retval,  /* return value */
             int        init     /* ==2 calc.strain energy */
             );
/*-----------------------------------------------------------------------*
|  w1_mat_plast_mises3D.c                                      sh 8/02   |
|  constitutive matrix - forces - linear elastic- von Mises - 3D         |
|  uses a 3D-Material formulation; constitutive matrix and stresses      |
|  are condensed back to 2D conditions (plane stress, plane strain)      |
|  (rotational symmetry is not implemented)                              |
*-----------------------------------------------------------------------*/
void w1_mat_plast_mises_3D(double    ym,     /* young's modulus              */
                           double    pv,     /* poisson's ratio              */
                           double    ALFAT,  /* temperature expansion factor */
                           double    sigy,   /* yield stresse                */
                           double    hard,   /* hardening modulus            */
                           double    gf,     /* fracture energy              */
                           double    betah,
                           ELEMENT  *ele,    /* actual element               */
                           WALL_TYPE wtype,  /* plane stress/strain...       */
                           double  **bop,    /* derivative operator          */
                           double   *gop,    /* add. derivative operator     */
                           double   *alpha,  /* internal dof                 */
                           int       ip,     /* integration point Id         */
                           double   *stress, /* vector of stresses           */
                           double  **d,      /* constitutive matrix          */
                           int       istore, /* controls storing of stresses */
                           int       newval);/* controls eval. of stresses   */
/*----------------------------------------------------------------------*
 | w1_mat_plast_mises3D.c                                               |
 | constitutive matrix - forces - linear elastic- von Mises - 3D sh 7/02|
 |                Test for 3D Materialmodell in Wall-Element            |
 *----------------------------------------------------------------------*/
void mat_plast_mises_3D(double ym,       /* young's modulus               */
                        double pv,       /* poisson's ration              */
                        double ALFAT,    /* temperature expansion factor  */
                        double sigy,     /* uniaxial yield stress         */
                        double hard,     /* hardening modulus             */
                        double gf,       /* fracture energy               */
                        double betah,
                        double *stress,  /* ele stress (-resultant) vector*/      
                        double **d,      /* material matrix               */
                        int   *iupd,     /* controls update of new stresses to WA */
                                         /* additional parameters         */
                        int    *yip,     /* stress state  1=el;2=pl       */
                        double *epstn,   /* uniaxial plastic strai        */
                        double strain[6],/* current strains from displacements */
                        double sig[6],   /* stresse from last update -> WA     */
                        double eps[6],   /* strains from last update -> WA     */
                        double qn[6],    /* backstress vector from last update -> WA */
                        double dia);     /* internal length parameter from WA  */
/*----------------------------------------------------------------------*
 | w1_mat_plast_mises3D.c                                               |
 | constitutive matrix - linear elastic - 3D [11,22,33,12,23,13] sh 7/02|
 *----------------------------------------------------------------------*/
void mat_linel3D(double ym,       /* young's modulus                    */
                 double pv,       /* poisson's ration                   */
                 double **d);     /* material matrix 3D                 */
/*----------------------------------------------------------------------*
 | w1_mat_plast_mises3D.c                                               |
 | yield criterion for mises plasticity - 3D                     sh 7/02|
 *----------------------------------------------------------------------*/
void yilcr_mises_lin(double E,    /* young's modulus                    */
                     double Eh,   /* hardening modulus                  */
                     double betah,
                     double sigy, /* uniaxial yield stress              */
                     double epstn,/* equivalent uniaxial plastic strain */
                     int    isoft,/* ??                                 */
                     double dia,  /* internal length parameter          */
                     double *tau, /* current stress (lokal)             */
                     double *ft); /* yield function                     */
/*----------------------------------------------------------------------*
 | w1_mat_plast_mises3D.c                                               |
 | radial return for mises plasticity with combined linear              |
 |               isotropic/kinematic hardening law - 3D          sh 7/02|
 *----------------------------------------------------------------------*/
void radi_mises_lin(double e,     /* young's modulus                    */
                    double eh,    /* hardening modulus                  */
                    double betah,
                    double sigy,  /* uniaxial yield stress              */
                    double vnu,   /* poisson's ration                   */
                    double dia,   /* internal length parameter          */
                    double *sigma,/* stress to be projected             */
                    double *qn,   /* backstress vector                  */
                    int    isoft, /* ??                                 */
                    double *epstn,/* uniaxial plastic strain            */
                    double *dlam);/* plastic multiplier                 */
/*----------------------------------------------------------------------*
 | w1_mat_plast_mises3D.c                                               |
 | forms the consistent material tangent for mises plasticity with      |
 | combined linear isotropic/kinematic hardening law - 3D        sh 7/02|
 *----------------------------------------------------------------------*/
void mapl_mises_lin(double e,     /* young's modulus                    */
                    double eh,    /* hardening modulus                  */
                    double betah,
                    double sigy,  /* uniaxial yield stress              */
                    double vnu,   /* poisson's ration                   */
                    double dia,   /* internal length parameter          */
                    double *tau,  /* stresses: predictor (plane_strain) */
                                  /*           projected (plane_stress) */ 
                    int    isoft, /* ??                                 */
                    double *epstn,/* uniaxial plastic strain            */
                    double dlam,  /* plastic multiplier                 */
                    double **d);  /* material matrix 3D                 */
/*-----------------------------------------------------------------------|
| w1_mat_trans.c                                                         |
|      topic: blowing up plane stress/strain conditions           sh 7/02|
|             to 3D --> 3D-Material Law                                  |
|-----------------------------------------------------------------------*/
void w1mat_trans_up (double     ym,        /* Young's modulus               */
                     double     pv,        /* poissons ratio                */
                     ELEMENT   *ele,       /* current element               */                                 
                     WALL_TYPE  wtype,     /* Info about actual wall element*/
                     double   **bop,       /* B-operator                    */
                     double    *gop,
                     double    *alpha,
                     int        ip,        /* integration point ID          */ 
                     double    *stress,    /* current stress condensed      */
                     double    *stress3D,  /* current stress condensed [6]  */
                     double    *strain3D,  /* strains to be calculated [6]  */
                     double    *sig3D,     /* stresses from last update [6] */ 
                     double    *eps3D,     /* strains from last update [6]  */
                     double    *qn3D,      /* backstress vektor [6]         */
                     int        newval);   /* controls evaluation of new stresses */           
/*-----------------------------------------------------------------------|
| w1_mat_trans.c                                                         |
|      topic: kondense 3D conditions                              sh 7/02|
|             to plane stress/strain conditions                          |
|-----------------------------------------------------------------------*/
void w1mat_trans_down (double   **d,        /*current material matrix 3D -> 2D */
                       ELEMENT   *ele,      /* current element                 */ 
                       WALL_TYPE  wtype,    /* Info about actual wall element  */
                       int        ip,       /* integration point ID            */
                       int        yipc,     /* stress state  1=el;2=pl -> copy of yip for correct saving */
                       double    *stressc,  /* condensed stresses              */
                       double    *sig,      /* stresses from last update       */
                       double    *eps,      /* strains from last update        */
                       double    *stress,   /* current stresses                */
                       double    *strain,   /* current strains                 */
                       double    *qn);      /* backstress vector               */
/*-----------------------------------------------------------------------*
| w1_mat_trans.c                                                         |
|     changes 2 rows of a vector                               sh 8/02   |
|  vec[x,x,a,x,b,x,...] -> vec[x,x,b,x,a,x,...]                          |
*-----------------------------------------------------------------------*/
void w1_vec_switch(double *vec,      /* vector do be modified           */
                   int a,            /* row to be changed to b          */
                   int b);           /* row to be changed to a          */

/*-----------------------------------------------------------------------*
| w1_mat_trans.c                                                         |
|     changes 2 rows & columns of a square matrix              sh 8/02   |
|     [-,a,-,b,-]        [-,b,-,a,-]                                     |
|     [a,a,a,c,a]        [b,b,b,c,b]                                     |
|  mat[-,a,-,b,-] ->  mat[-,b,-,a,-]                                     |
|     [b,c,b,b,b]        [a,c,a,a,a]                                     |
|     [-,a,-,b,-]        [-,b,-,a,-]                                     |
*-----------------------------------------------------------------------*/
void w1_matrix_switch(double **mat,   /* matrix do be modified          */
                      int a,          /* row & colum to be changed to b */     
                      int b,          /* row & colum to be changed to a */
                      int l);         /* length of row/column of matrix */
/*----------------------------------------------------------------------*/
#endif /*D_WALL1*/
/*! @} (documentation module close)*/
                                                                        
