/*!----------------------------------------------------------------------
\file
\brief contains all prototypes for gradient enhance wall element

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0711 - 685-6122
</pre>

*----------------------------------------------------------------------*/
#ifdef D_MLSTRUCT

/*! 
\addtogroup ML_struct
*//*! @{ (documentation module open)*/

/*-----------------------------------------------------------------------*
|  s2ml_static_ke.c                                           ah 03/04  |
*-----------------------------------------------------------------------*/
void s2ml_static_ke(ELEMENT       *actmaele, 
                    W1_DATA       *data, 
                    MATERIAL      *mat,
                    ARRAY         *estif_global,
                    ARRAY         *emass_global,
                    DOUBLE        *force,  
                    INT            init);
/*-----------------------------------------------------------------------*
|  s2ml_cal_smelm.c                                           ah 03/04  |
*-----------------------------------------------------------------------*/
void s2ml_init(ELEMENT  *actmaele,       /*  actual macro element       */
               INT       numsmnodes,     /*  number of submeshnodes     */
               INT       numsmele);      /*  number of submeshelements  */
/*-----------------------------------------------------------------------*
|  s2ml_service.c                                           ah 03/04  |
*-----------------------------------------------------------------------*/
void s2ml_bopstraintonode(FIELD       *actsmfield,/*  submesh field     */
                          ELEMENT     *ele, /* actual makroelement      */
                          DOUBLE       nue, /* poisson rat.of macro mat.*/
                          INT          init); 
/*-----------------------------------------------------------------------*
|  s2ml_cal_smelm.c                                           ah 03/04  |
*-----------------------------------------------------------------------*/
void s2ml_cal_smelm(ELEMENT     *actmaele,   /*  actual macro element    */
                    FIELD       *actsmfield, /*  submesh field           */
                    SOLVAR      *actsmsolv,  /*  submesh SOLVAR->sm assembled stif mi_mi*/
                    INTRA       *actsmintra, /* the sm pseudo intra-communicator */
                    DOUBLE      *smintforcemi, /*  sm assembled int force micro */
                    DBCSR       *smstiffmima_csr, /* sm assebled stiffness_mi_ma*/
                    DBCSR       *smstiffmami_csr, /* sm assebled stiffness_ma_mi*/
                    DOUBLE      *smintforcema, /*  sm "assembled" int force macro */
                    DOUBLE     **smstiffmama, /*  sm "assembled" stiffness macro-macro */
                    INT          istore);    /*  is it update step?      */
/*-----------------------------------------------------------------------*
|  s2ml_stiff_wall.c                                           ah 03/04  |
*-----------------------------------------------------------------------*/
void s2ml_stiff_wall(MATERIAL  *actsmmat,    /* actual submesh material*/
                     ELEMENT   *actmaele,    /* actual macro element   */
                     ELEMENT   *actsmele,    /* actual submesh element */
                     ARRAY     *estif_ma_ma, /* stiffness macro-macro  */
                     ARRAY     *estif_ma_mi, /* stiffness macro-micro  */
                     ARRAY     *estif_mi_ma, /* stiffness micro-macro  */
                     ARRAY     *estif_mi_mi, /* stiffness micro-micro  */
                     ARRAY     *intforce_ma, /* internal force macro   */
                     ARRAY     *intforce_mi, /* internal force micro   */
                     INT        istore,      /* update?                */
                     INT        init);       /* allocate arrays?       */
/*-----------------------------------------------------------------------*
|  s2ml_strain.c                                               ah 03/04  |
*-----------------------------------------------------------------------*/
void s2ml_aequistrain(ELEMENT     *actmaele, /* actual Macro-element   */
                      WALL_TYPE    wtype,    /* plane stress/strain... */         
                      DOUBLE     **bopma,    /* Macro-B-operator       */  
                      DOUBLE       nue,      /* poisson-ratio          */
                      DOUBLE      *eps_equiv);/* equivalent strain      */
/*-----------------------------------------------------------------------*
|  s2ml_strain.c                                               ah 03/04  |
*-----------------------------------------------------------------------*/
void s2ml_bopstrainma(DOUBLE     **bopma,     /* macro B-operator at GP      */
                      DOUBLE      *strainma,  /* macro strain at GP          */
                      DOUBLE      *functmi,   /* micro Ansatzfunctions at GP */
                      ELEMENT     *actsmele); /* actual submesh element      */
/*-----------------------------------------------------------------------*
|  s2ml_strain.c                                               ah 03/04  |
*-----------------------------------------------------------------------*/
void s2ml_strainmi(DOUBLE     **bopmi,     /* micro B-operator at GP    */
                   DOUBLE      *strainmi,  /* micro strain at GP        */
                   ELEMENT     *actsmele,  /* actual submesh element    */
                   ELEMENT     *actmaele,  /* actual submesh element    */
                   DOUBLE       nue);      /* poisson ratio             */
/*-----------------------------------------------------------------------*
|  s2ml_strain.c                                               ah 04/04  |
*-----------------------------------------------------------------------*/
void s2ml_ifbopmajumpu(ELEMENT  *actsmele,    /* actual submesh element    */
                       ELEMENT  *actmaele,    /* actual macroe element     */
                       DOUBLE  **bopmi,       /* micro B-operator at GP    */
                       DOUBLE  **bopma,       /* macro B-operator at GP    */
                       DOUBLE   *jumpuma,     /* displacement jump (macro) */
                       DOUBLE   *DELTAjumpuma,/* incre disp jump (macro)   */
                       DOUBLE   *jumpumi,     /* displacement jump (micro) */
                       DOUBLE   *DELTAjumpumi);/* incre disp jump (micro)   */
/*-----------------------------------------------------------------------*
|  s2ml_serive.c                                             ah 03/04  |
*-----------------------------------------------------------------------*/
void s2ml_callmat(ELEMENT     *actmaele,   /*  actual macro element   */
                  ELEMENT     *actsmele,   /*  actual submeshelement  */
                  MATERIAL    *actsmmat,   /*  actual sm_material     */
                  DOUBLE      *strain,     /* total strain at GP      */
                  INT          ip,         /* ip of GP                */
                  DOUBLE      *stress,     /*  stress                 */
                  DOUBLE     **D,          /*  tangent                */
                  INT          istore);    /*  is it update step?     */
/*-----------------------------------------------------------------------*
|  s2ml_service.c                                           ah 03/04  |
*-----------------------------------------------------------------------*/
void s2ml_cal_stiff(DOUBLE   **stiffmatrix, /* element stiffnes matrix */
                    DOUBLE   **left_bop,    /* left B-operator (to be transp) */
                    DOUBLE   **right_bop,   /* right B-operator        */
                    DOUBLE   **D,           /* material tangent        */
                    DOUBLE   fac,           /* integration factor      */
                    INT      left_nd,       /* number ele-dof left vec */
                    INT      right_nd,      /* number ele-dof right vec*/
                    INT      numeps);       /* number of strains-compon*/
/*-----------------------------------------------------------------------*
|  s2ml_stiff_interf.c                                        ah 04/04  |
*-----------------------------------------------------------------------*/
void s2ml_stiff_interf(MATERIAL  *actsmmat,    /* actual submesh material*/
                       ELEMENT   *actmaele,    /* actual macro element   */
                       ELEMENT   *actsmele,    /* actual submesh element */
                       ARRAY     *estif_ma_ma, /* stiffness macro-macro  */
                       ARRAY     *estif_ma_mi, /* stiffness macro-micro  */
                       ARRAY     *estif_mi_ma, /* stiffness micro-macro  */
                       ARRAY     *estif_mi_mi, /* stiffness micro-micro  */
                       ARRAY     *intforce_ma, /* internal force macro   */
                       ARRAY     *intforce_mi, /* internal force micro   */
                       INT        istore,      /* update?                */
                       INT        init);       /* allocate arrays?       */
/*----------------------------------------------------------------------*/

/*! @} (documentation module close)*/

#endif /* D_MLSTRUCT */
