/*!----------------------------------------------------------------------
\file
\brief contains all prototypes for gradient enhance wall element

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0771 - 685-6122
</pre>

*----------------------------------------------------------------------*/
#ifdef D_WALLGE

/*! 
\addtogroup WALLGE
*//*! @{ (documentation module open)*/

/*-----------------------------------------------------------------------*
|  wge_setdof.c                                                 mn 05/03  |
*-----------------------------------------------------------------------*/
void wge_setdof(); 
/*-----------------------------------------------------------------------*
|  wge_init.c                                                 mn 05/03  |
*-----------------------------------------------------------------------*/
void wgeinit(PARTITION *actpart,MATERIAL *mat);
/*-----------------------------------------------------------------------*
|  wge_init.c                                                 mn 05/03  |
*-----------------------------------------------------------------------*/
void wgeintg(ELEMENT       *ele,     
             WALLGE_DATA   *data,
             INT            option);
/*-----------------------------------------------------------------------*
|  wge_static_ke.c                                            mn 05/03  |
*-----------------------------------------------------------------------*/
void wgestatic_ke(ELEMENT       *ele, 
                  WALLGE_DATA   *data, 
                  MATERIAL      *mat,
                  ARRAY         *estif_global,
                  ARRAY         *emass_global,
                  DOUBLE        *force,
                  INT            init);
/*-----------------------------------------------------------------------*
|  wge_bope.c                                            mn 05/03  |
*-----------------------------------------------------------------------*/
void wge_bope(DOUBLE    **bope,          
              DOUBLE    **derive,        
              DOUBLE    **xjm,          
              DOUBLE      det); 
/*-----------------------------------------------------------------------*
|  wge_call_mat.c                                            mn 05/03  |
*-----------------------------------------------------------------------*/
void wge_call_mat(ELEMENT   *ele,     /* actual element                */
                  MATERIAL  *mat,     /* material                      */
                  DOUBLE   **bopd,    /* B-operator for displacements  */
                  DOUBLE    *functe,  /* Ansatz-funct.for equiv.strain */
                  DOUBLE   **bope,    /* B-operator for equiv.strain   */
                  INT        ip,      /* ID of actual Gauss point      */
                  DOUBLE    *stress,  /* stresses at GP                */
                  DOUBLE    *eps_vl,  /* local equiv. strain at GP     */
                  DOUBLE    *eps_vnl, /* nonlocal equiv. strain at GP  */
                  DOUBLE    *grad_eps_vnl, /* grad of nonl.equi.strain */
                  DOUBLE   **D,       /* tangent d sig/d eps           */
                  DOUBLE    *E,       /* tangent d sig/d eps_v         */
                  DOUBLE    *F,       /* tangent d eps_v/d eps         */
                  INT        istore,  /* flag: istore=1 -> update      */
                  INT        newval); /* flag: newval=1 -> just stress */ 
/*-----------------------------------------------------------------------*
|  wge_call_mat.c                                            mn 05/03  |
*-----------------------------------------------------------------------*/
void wge_mat_damage(INT      equival,  /* flag for equivalent strains     */
                    INT      damtyp,   /* flag for Damage-Typ             */
                    DOUBLE   youngs,   /* young's modulus                 */
                    DOUBLE   nue,      /* poisson's ratio                 */
                    DOUBLE   kappa_0,  /* initial damage equivalent strain*/
                    DOUBLE   kappa_m,  /* factor for damage-law           */
                    DOUBLE   alpha,    /* factor for exp. damage-law      */
                    DOUBLE   beta,     /* factor for exp. damage-law      */
                    DOUBLE   k_fac,    /* de Vree                         */
                    ELEMENT *ele,      /* actual element                  */
                    DOUBLE **bopd,     /* B-operator for displacements    */
                    DOUBLE  *functe,   /* Ansatz-funct. for equiv. strain */
                    DOUBLE **bope,     /* B-operator for equiv. strain    */
                    INT      ip,       /* integration point Id            */
                    DOUBLE  *stress,   /* stress vector                   */
                    DOUBLE  *eps_vl,   /* local equivalent strain         */
                    DOUBLE  *eps_vnl,  /* nonlocal equivalent strain      */
                    DOUBLE  *grad_eps_vnl,/* grad of nonlocal equi.strain */
                    DOUBLE **D,        /* 1. Material tangent             */
                    DOUBLE  *E,        /* 2. Material tangent             */
                    DOUBLE  *F,        /* 3. Material tangent             */
                    INT      istore,
                    INT      newval);       
/*-----------------------------------------------------------------------*
|  wge_strain.c                                            mn 05/03  |
*-----------------------------------------------------------------------*/
void wge_strain(ELEMENT   *ele,     /* actual element                  */
                DOUBLE   **bopd,    /* B-operator for displacements    */
                DOUBLE    *functe,  /* Ansatz-funct. for equiv. strain */
                DOUBLE   **bope,    /* B-operator for equiv. strain    */
                DOUBLE    *strain,  /* stress vector                   */
                DOUBLE    *eps_vnl, /* nonlocal equivalent strain      */
                DOUBLE    *grad_eps_vnl);/* grad nonl. equiv. strain   */
/*-----------------------------------------------------------------------*
|  wge_damage_service.c                                      mn 05/03  |
*-----------------------------------------------------------------------*/
void wge_damvar(INT     damtyp,        /*  Damage law                   */
                DOUBLE  kappa,         /*  actual nonl. eqvival. strain */
                DOUBLE  kappa_0,       /* parameter of exp. damage law  */
                DOUBLE *eps_vnl,       /* nonlocal equiv. strain        */
                DOUBLE  alpha,         /* parameter of exp. damage law  */
                DOUBLE  beta,          /* parameter of exp. damage law  */
                DOUBLE *damage,        /* actual damage variable        */
                DOUBLE *dam_deriv);    /* (partial D)/(partial eps_equi)*/       
/*-----------------------------------------------------------------------*
|  wge_damage_service.c                                      mn 05/03  |
*-----------------------------------------------------------------------*/
void wge_epsequiv(INT     equival,       /*definition of equiv. strains*/
                  DOUBLE  epsilon[3][3],    /*  strain tensor          */
                  DOUBLE  delta[3][3],      /*  kronerker delta        */
                  DOUBLE  k_f,              /*  parameter for de Vree  */
                  DOUBLE  nue,              /*  poisson rate           */
                  DOUBLE *eps_vl,           /*  local equiv. strains   */
                  DOUBLE  F_ed[3][3]);      /* equ. strain derivative  */         
/*-----------------------------------------------------------------------*
|  wge_damage_service.c                                      mn 05/03  |
*-----------------------------------------------------------------------*/
void wge_stress(DOUBLE   youngs,           /*  youngs modulus          */
                DOUBLE   nue,              /*  poisson ratio           */
                DOUBLE   delta[3][3],      /*  kroneker delta          */
                DOUBLE   epsilon[3][3],    /*  strain tensor           */
                DOUBLE   damage,           /*  actual damage variable  */
                DOUBLE   sigma_el[3][3],   /*  elastic stress tensor   */ 
                DOUBLE   sigma[3][3]);     /*  stress tensor           */      
/*-----------------------------------------------------------------------*
|  wge_damage_service.c                                      mn 05/03  |
*-----------------------------------------------------------------------*/
void wge_tangent(DOUBLE   youngs,           /*  youngs modulus         */
                 DOUBLE   nue,              /*  poisson ratio          */
                 DOUBLE   delta[3][3],      /*  kroneker delta         */
                 DOUBLE   damage,           /*  damaage variable       */
                 DOUBLE   C_ed[3][3][3][3]);/*  elasto-damage-tangent  */         
/*-----------------------------------------------------------------------*
|  wge_damage_service.c                                      mn 05/03  |
*-----------------------------------------------------------------------*/
void wge_condense(DOUBLE    **D,     /* 1.elasto-damage-tangent-matrix */
                  DOUBLE     *E,     /* 2.elasto-damage-tangent-matrix */
                  DOUBLE     *F,     /* 3.elasto-damage-tangent-matrix */
                  WALLGE_TYPE wtype, /* plane-stress or plane strain   */
                  DOUBLE      nue);  /* poisson ratio                  */         
/*-----------------------------------------------------------------------*
|  wge_stiff.c                                      mn 05/03  |
*-----------------------------------------------------------------------*/
void wge_stiff_de(DOUBLE  **Kde,      /* stiffnes displ - equiv,strain */
                  DOUBLE  **bopd,     /* B-operator for displacements  */
                  DOUBLE   *E,        /* Material tangent vector       */
                  DOUBLE   *functe,   /* Anatzfunc. for equiv.strain   */
                  DOUBLE    fac,      /* integration factor            */
                  INT       numdfd,   /* element displacement DOF      */
                  INT       iele,     /* element equiv.strain DOF      */
                  INT       neps);    /* number of strain components   */       
/*-----------------------------------------------------------------------*
|  wge_stiff.c                                      mn 05/03  |
*-----------------------------------------------------------------------*/
void wge_stiff_ed(DOUBLE  **Ked,      /* stiffnes equiv,strain - disp  */
                  DOUBLE  **bopd,     /* B-operator for displacements  */
                  DOUBLE   *F,        /* Material tangent vector       */
                  DOUBLE   *functe,   /* Anatzfunc. for equiv.strain   */
                  DOUBLE    fac,      /* integration factor            */
                  INT       numdfd,   /* element displacement DOF      */
                  INT       iele,     /* element equiv.strain DOF      */
                  INT       neps);    /* number of strain components   */       
/*-----------------------------------------------------------------------*
|  wge_stiff.c                                      mn 05/03  |
*-----------------------------------------------------------------------*/
void wge_stiff_ee(DOUBLE  **Kee,      /* stiffnes equiv,strain - disp  */
                  DOUBLE  **bope,     /* B-operator for displacements  */
                  DOUBLE    crad,     /* influenceradius nonloc.strain */
                  DOUBLE   *functe,   /* Anatzfunc. for equiv.strain   */
                  DOUBLE    fac,      /* integration factor            */
                  INT       iele);     /* element equiv.strain DOF      */      
/*-----------------------------------------------------------------------*
|  wge_stiff.c                                      mn 05/03  |
*-----------------------------------------------------------------------*/
void wge_fintd(DOUBLE   *stress,      /* stresses                      */ 
               DOUBLE    fac,         /* integration factor            */
               DOUBLE   **bopd,       /* B-operator for displacements  */
               INT        numdfd,     /* element-displacement DOF      */
               INT        neps,       /* number of strain components   */
               DOUBLE    *fintd);     /* int. force(displacementDOF)   */      
/*-----------------------------------------------------------------------*
|  wge_stiff.c                                      mn 05/03  |
*-----------------------------------------------------------------------*/
void wge_finte(DOUBLE     eps_vl,       /* local equiv. strain         */
               DOUBLE     eps_vnl,      /* nonlocal equiv. strain      */
               DOUBLE    *grad_eps_vnl, /* grad nonlocal equiv. strain */
               DOUBLE     crad,         /* material parameter          */
               DOUBLE     fac,          /* integration factor          */
               DOUBLE   **bope,         /* B-operator for equiv.strain */
               DOUBLE    *functe,       /* Ansatz-fun for equiv.strain */
               INT        iele,         /* eleDOF for equiv.strain     */
               DOUBLE    *finte);       /* int. force(equiv.strainDOF) */      
/*-----------------------------------------------------------------------*
|  wge_stiff.c                                      mn 05/03  |
*-----------------------------------------------------------------------*/
void wge_permstiff(DOUBLE    **Kdd,     /* 1.stiffness part            */
                   DOUBLE    **Kde,     /* 2.stiffness part            */
                   DOUBLE    **Ked,     /* 3.stiffness part            */
                   DOUBLE    **Kee,     /* 4.stiffness part            */
                   INT         iele,    /* node-number of equiv.strain */
                   INT         ield,    /* node-number of displacement */
                   INT         numdf,   /* total DOF                   */
                   DOUBLE    **estif);  /* "mixed" element stiffness   */      
/*-----------------------------------------------------------------------*
|  wge_stiff.c                                      mn 05/03  |
*-----------------------------------------------------------------------*/
void wge_permforce(DOUBLE     *fintd,   /* 1.part int. force           */
                   DOUBLE     *finte,   /* 2.part int. force           */
                   INT         iele,    /* num.of equiv.strain nodes   */
                   INT         ield,    /* num.of displacement nodes   */
                   INT         numdf,   /* total element DOF           */
                   DOUBLE     *force);  /* "mixed" element int. force  */      
/*-----------------------------------------------------------------------*
|  wge_stress.c                                             mn 05/03  |
*-----------------------------------------------------------------------*/
void wge_cal_stress(ELEMENT   *ele); 
/*-----------------------------------------------------------------------*
|  wge_fext.c                                             mn 05/03  |
*-----------------------------------------------------------------------*/
void wge_eleload(ELEMENT     *ele,       /* actual element              */
                 WALLGE_DATA *data,      /* element integration data    */
                 DOUBLE      *loadvec,   /* external element forces     */
                 INT          init); 

/*----------------------------------------------------------------------*/
#endif /*D_WALLGE*/
/*! @} (documentation module close)*/
