/*!----------------------------------------------------------------------
\file
\brief headerfile for multilayer shell element, containing structures and
       prototypes

*----------------------------------------------------------------------*/
#ifdef D_SHELL9

/*! 
\addtogroup SHELL9 
*//*! @{ (documentation module open)*/

/*----------------------------------------------------------------------*
 | shell9             modified from shell8              by     sh 10/02 |
 *----------------------------------------------------------------------*/
typedef struct _SHELL9
{
/*----------------------------- general variables needed by the element */
struct _ARRAY    a3ref;
struct _ARRAY    thick_node;
DOUBLE           thick;
DOUBLE           sdc;
INT              nGP[3];
INT              nGP_tri;
INT              num_klay;      /* number of kinematic layers */
DOUBLE          *klayhgt;       /* hgt of a kinematic layer in % of total thickness of the shell */
struct _KINLAY  *kinlay;        /* one kinematic layer -> struct is defined in materials.h */
INT              numdf;         /* number of dofs per node of the s9 element */
/*---------------------------------------------- array of forces at GPs */
enum
    {
    s9_xyz,
    s9_rst,
    s9_rst_ortho
    }            forcetyp;    /*-> forces and stresses in same coordinate system*/
struct _ARRAY4D  forces;      /*stress resultants at GPs -> mid surface */
struct _ARRAY4D  stresses;    /*physical stresses at GPs*/
/*struct _ARRAY    energy;*/
/*-------------------------------------------- array of internal forces */
struct _ARRAY    intforce;
/*-------------------------------------------- variables needed for eas */
INT              oldkstep;                
INT              eas[5];
INT              nhyb; 
struct _ARRAY    alfa;
struct _ARRAY    Dtildinv[MAXKLAY_SHELL9];
struct _ARRAY    Lt[MAXKLAY_SHELL9];
struct _ARRAY    Rtilde[MAXKLAY_SHELL9];
/*-------------------------------------------- variables needed for ans */
INT              ans;
} SHELL9;


/*----------------------------------------------------------------------*
 | shell9 data                                            m.gee 6/01    |
 *----------------------------------------------------------------------*/
typedef struct _S9_DATA
{
DOUBLE        xgpr[3];
DOUBLE        wgtr[3];

DOUBLE        xgps[3];
DOUBLE        wgts[3];

DOUBLE        xgpt[3];
DOUBLE        wgtt[3];
} S9_DATA;


/*----------------------------------------------------------------------*
 |  PROTOTYPES OF ALL SHELL9 ROUTINES      (shell8: m.gee 11/01)        |
 |                modified for shell9                      by sh 10/02  |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 |  s9_a3.c                                              m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s9a3(ELEMENT   *ele);
void s9a3ref_extern(DOUBLE   *funct,
                    DOUBLE  **deriv,
                    DOUBLE  **a3ref, 
                    ELEMENT  *ele);
void s9averdir(DOUBLE **dir_list, INT numa3, DOUBLE *a3, DOUBLE h2);
/*----------------------------------------------------------------------*
 |  s9_ans.c                                            m.gee 03/02    |
 *----------------------------------------------------------------------*/
void s9_ans_colloqpoints(INT       nsansq,     INT       iel,
                         INT       ans,        DIS_TYP   distyp,
                         DOUBLE    xr1[],      DOUBLE    xs1[],
                         DOUBLE    xr2[],      DOUBLE    xs2[],
                         DOUBLE   *funct1q[],  DOUBLE  **deriv1q[], 
                         DOUBLE   *funct2q[],  DOUBLE  **deriv2q[],
                         DOUBLE  **xrefe,      DOUBLE ***a3r,
                         DOUBLE  **xcure,      DOUBLE ***a3c,
                         DOUBLE ***akovr1q[],  DOUBLE ***akonr1q[],
                         DOUBLE ***amkovr1q[], DOUBLE ***amkonr1q[],
                         DOUBLE ***a3kvpr1q[],
                         DOUBLE ***akovc1q[] , DOUBLE ***akonc1q[],
                         DOUBLE ***amkovc1q[], DOUBLE ***amkonc1q[],
                         DOUBLE ***a3kvpc1q[],
                         DOUBLE ***akovr2q[] , DOUBLE ***akonr2q[],
                         DOUBLE ***amkovr2q[], DOUBLE ***amkonr2q[],
                         DOUBLE ***a3kvpr2q[],
                         DOUBLE ***akovc2q[] , DOUBLE ***akonc2q[],
                         DOUBLE ***amkovc2q[], DOUBLE ***amkonc2q[],
                         DOUBLE ***a3kvpc2q[],
                         DOUBLE  **akovh,
                         DOUBLE  **akonh,
                         DOUBLE  **amkovh,
                         DOUBLE  **amkonh,
                         INT       num_klay,   INT       numdf);
void s9_ans_colloqcoords(DOUBLE xqr1[], DOUBLE xqs1[],
                         DOUBLE xqr2[], DOUBLE xqs2[],
                         INT    iel,    INT    ans);
void s9_ansq_funct(DOUBLE frq[], DOUBLE fsq[], DOUBLE r, DOUBLE s,
                   INT    iel,   INT    nsansq);
void s9_ans_bbar_q(DOUBLE  **bop,        DOUBLE    frq[],      DOUBLE fsq[],
                   DOUBLE   *funct1q[],  DOUBLE   *funct2q[],
                   DOUBLE  **deriv1q[],  DOUBLE  **deriv2q[],
                   DOUBLE ***akovc1q[],  DOUBLE ***akovc2q[],
                   DOUBLE ***a3kvpc1q[], DOUBLE ***a3kvpc2q[],
                   INT       iel,        INT       numdf, 
                   INT       num_klay,   INT       klay,
                   DOUBLE    condfac,    INT       nsansq);
void s9_ans_tvhe_q(DOUBLE  **gmkovr,     DOUBLE  **gmkovc,    
                   DOUBLE  **gmkonr,     DOUBLE  **gmkonc,
                   DOUBLE  **gkovr,      DOUBLE  **gkovc,     
                   DOUBLE ***amkovc,     DOUBLE ***amkovr,
                   DOUBLE ***akovc,      DOUBLE ***akovr,     
                   DOUBLE ***a3kvpc,     DOUBLE ***a3kvpr,
                   DOUBLE   *detr,       DOUBLE   *detc,
                   DOUBLE ***amkovr1q[], DOUBLE ***amkovc1q[], 
                   DOUBLE ***akovr1q[] , DOUBLE ***akovc1q[] ,
                   DOUBLE ***a3kvpr1q[], DOUBLE ***a3kvpc1q[],
                   DOUBLE ***amkovr2q[], DOUBLE ***amkovc2q[], 
                   DOUBLE ***akovr2q[] , DOUBLE ***akovc2q[] ,
                   DOUBLE ***a3kvpr2q[], DOUBLE ***a3kvpc2q[],
                   DOUBLE    frq[],      DOUBLE    fsq[], 
                   DOUBLE    e3,         INT       nansq, 
                   INT       iel,       
                   DOUBLE    h,           /* total thickness of this element */
                   DOUBLE   *klayhgt,     /* hight of kin layer in % of total thickness of shell */
                   DOUBLE   *mlayhgt,     /* hight of mat layer in % of adjacent kin layer */
                   INT       num_klay,    /* number of kin layers to this element */
                   INT       num_mlay,    /* number of mat layers to this kin layer */
                   INT       klay,        /* actual kin layer */
                   INT       mlay,        /* actual mat layer of this kin layer */
                   DOUBLE    condfac);
void s9_ans_tvkg(DOUBLE **estif,     DOUBLE  *stress_r,
                 DOUBLE  *funct,     DOUBLE **deriv,
                 INT      numdf,     INT      iel,         
                 DOUBLE   weight,
                 DOUBLE   frq[],     DOUBLE   fsq[],
                 DOUBLE  *funct1q[], DOUBLE  *funct2q[],
                 DOUBLE **deriv1q[], DOUBLE **deriv2q[],
                 INT      ansq,      INT       nsansq,
                 INT      klay,                       /* actual kin layer */
                 INT      num_klay);                  /* number of kin layers to this element */ 
/*----------------------------------------------------------------------*
 |  s9_btdb.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s9_BtDB(DOUBLE **estif, DOUBLE **bop,    DOUBLE **D,   INT iel,
             INT      numdf, DOUBLE   weight, DOUBLE **work);
/*----------------------------------------------------------------------*
 |  s9_cal_dyn.c                                         m.gee 02/02    |
 *----------------------------------------------------------------------*/
void s9_tmas(DOUBLE *funct, DOUBLE *thick, DOUBLE **emass, INT iel, INT numdf,
             DOUBLE facv,   DOUBLE facw,   DOUBLE facvw);
/*----------------------------------------------------------------------*
 |  s9_call_mat.c                                           sh 02/03    |
 *----------------------------------------------------------------------*/
void s9_call_mat(ELEMENT    *ele,
                 MULTIMAT   *multimat,    /* material of actual material layer */
                 DOUBLE     *stress,
                 DOUBLE     *strain,
                 DOUBLE    **C,
                 DOUBLE    **gmkovc,
                 DOUBLE    **gmkonc,
                 DOUBLE    **gmkovr,
                 DOUBLE    **gmkonr,
                 DOUBLE    **gkovc,
                 DOUBLE    **gkonc,
                 DOUBLE    **gkovr,
                 DOUBLE    **gkonr,
                 INT         rot_axis,
                 DOUBLE      phi);
void s9_getdensity(MATERIAL   *mat, DOUBLE *density);
/*----------------------------------------------------------------------*
 |  s9_eas.c                                             m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s9_eas(const INT    nhyb,
            const DOUBLE e1,
            const DOUBLE e2,
            const INT    iel, 
            const INT   *eas, 
            DOUBLE     **P);
void s9_transeas(DOUBLE      **P, 
                 DOUBLE      **transP,
                 DOUBLE      **T,
                 DOUBLE     ***akovr,
                 DOUBLE      **akonr0,
                 DOUBLE        detr,
                 DOUBLE        detr0,
                 INT           nhyb,
                 INT           klay);
/*----------------------------------------------------------------------*
 |  s9_eps.c                                             m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s9_eps(DOUBLE *strain,DOUBLE **gmkovc, DOUBLE **gmkovr);
/*----------------------------------------------------------------------*
 |  s9_funcderiv.c                                       m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s9_funct_deriv(DOUBLE     *funct, 
                    DOUBLE    **deriv, 
                    DOUBLE      r, 
                    DOUBLE      s,
                    DIS_TYP     typ,
                    INT         option);
/*----------------------------------------------------------------------*
 |  s9_init.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s9init(FIELD *actfield);
/*----------------------------------------------------------------------*
 |  s9_inpele.c                                          m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s9inp(ELEMENT *ele);
/*----------------------------------------------------------------------*
 |  s9_intforce.c                                        m.gee 02/02    |
 *----------------------------------------------------------------------*/
void s9_intforce(DOUBLE *intforce, DOUBLE *stress_r, DOUBLE **bop,
                 INT iel, INT numdf, INT nstress_r, DOUBLE weight);
/*----------------------------------------------------------------------*
 |  s9_intg.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s9intg(const ELEMENT   *ele,
            S9_DATA         *data,
            INT              option);
/*----------------------------------------------------------------------*
 |  s9_jaco.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s9jaco(DOUBLE    *funct,
            DOUBLE   **deriv,
            DOUBLE   **x,
            DOUBLE   **xjm,
            DOUBLE    *hte,
            DOUBLE  ***a3r,
            DOUBLE     e3,
            INT        iel,
            DOUBLE    *deta,
            INT        init,
            INT        num_klay,
            DOUBLE    *klayhgt,     /* hight of kin layer in % of total thickness of shell */
            DOUBLE    *mlayhgt);    /* hight of mat layer in % of adjacent kin layer */
/*----------------------------------------------------------------------*
 |  s9_load1.c                                           m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s9eleload(ELEMENT  *ele,
               S9_DATA  *data,
               DOUBLE	*loadvec,
               INT	 init);
void s9loadGP(ELEMENT    *ele,
              DOUBLE    **eload1,
              DOUBLE      wgt,
              DOUBLE    **xjm,
              DOUBLE     *funct,
              DOUBLE    **deriv,
              INT         iel,
              DOUBLE      xi,
              DOUBLE      yi,
              DOUBLE      zi,
              DOUBLE      shift);    /*shell shifter if top,bot */
/*----------------------------------------------------------------------*
 |  s9_load_surf.c                                          sh 12/02    |
 *----------------------------------------------------------------------*/
void s9_surf(DOUBLE **eload, INT nsurf, INT numklay, INT iel);
void s9_surf_P(DOUBLE *pload, INT nsurf, INT numklay);
void s9_surf_onoff(INT *pload, INT nsurf, INT numklay);
/*----------------------------------------------------------------------*
 |  s9_loccoordnode.c                                    m.gee 11/01    |
 *----------------------------------------------------------------------*/
DOUBLE s9_local_coord_node(INT node, INT flag, enum _DIS_TYP typ);
/*----------------------------------------------------------------------*
 |  s9_main.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void shell9(FIELD       *actfield,
            PARTITION   *actpart,
            INTRA       *actintra,
            ELEMENT     *ele,
            ARRAY       *estif_global,
            ARRAY       *emass_global,
            ARRAY       *intforce_global,
            CALC_ACTION *action,
            CONTAINER   *container);    /* contains variables defined in container.h */
/*----------------------------------------------------------------------*
 |  s9_mat_linel.c                                       m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s9_mat_linel(STVENANT *mat, DOUBLE **g, DOUBLE **CC);
void s9_mat_stress1(DOUBLE *stress, DOUBLE *strain, DOUBLE **C);
/*----------------------------------------------------------------------*
 |  s9_mat_linel3D.c                                        sh 02/03    |
 *----------------------------------------------------------------------*/
void s9_mat_linel3D(DOUBLE youngs,
                    DOUBLE possionratio,
                    DOUBLE **d);
/*----------------------------------------------------------------------*
 |  s9_mat_orth3D.c                                         sh 02/03    |
 *----------------------------------------------------------------------*/
void s9_mat_orth3D(ORTHOTROPIC *mat, DOUBLE **d);
/*----------------------------------------------------------------------*
 |  s9_math.c                                               sh 02/03    |
 *----------------------------------------------------------------------*/
void s9_matrix_sort(DOUBLE **M, char flag1[], char flag2[]);
void s9tcuca (DOUBLE *str, DOUBLE **gkov, char flag1[],char flag2[]);
void s9T4sym (DOUBLE **A, DOUBLE **C);
void s9shear_cor(DOUBLE **d, DOUBLE xsi); 
void s9_Tmaca (DOUBLE **akov, DOUBLE phi, INT rot_axis, DOUBLE **C);
void s9_Tcacu (DOUBLE **akov, DOUBLE **C);
/*----------------------------------------------------------------------*
 |  s9_out.c                                                sh 12/02    |
 *----------------------------------------------------------------------*/
void s9out_dof(NODE *actnode,INT j);
void s9out_nodal_dis(NODE *actnode,INT place);
void s9_out_gid_allcoords(FILE *out,INT type);
void s9_out_gid_eletop(FILE *out,INT type,ELEMENT *actele);
void s9_out_gid_sol_dis(FILE *out,FIELD *actfield,INT place);
void s9_out_gid_sol_str(FILE *out,FIELD *actfield,INT place);
/*----------------------------------------------------------------------*
 |  s9_restart.c                                         m.gee 05/02    |
 *----------------------------------------------------------------------*/
void s9_write_restart(ELEMENT *actele, INT nhandle, long int *handles);
void s9_read_restart(ELEMENT *actele, INT nhandle, long int *handles);
/*----------------------------------------------------------------------*
 |  s9_static_keug.c                                       m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s9static_keug(ELEMENT   *ele,                        /* the element structure */
                   S9_DATA   *data,                       /* element integration data */
                   MATERIAL  *mat,                        /* the material structure */
                   ARRAY     *estif_global,               /* element stiffness matrix (NOT initialized!) */
                   ARRAY     *emass_global,               /* element mass matrix (NOT initialized!) */
                   INT        kintyp,                     /* typ of kinematic formulation */
                   DOUBLE    *force,                      /* global vector for internal forces (initialized!) */
                   INT        kstep,                      /* actual step in nonlinear analysis */
                   INT        init);                      /* init=1 -> init phase / init=0 -> calc. phase / init=-1 -> uninit phase */
/*----------------------------------------------------------------------*
 |  s9_stress.c                                          m.gee 02/02    |
 *----------------------------------------------------------------------*/
void s9_stress(ELEMENT      *ele,
               S9_DATA      *data,
               MATERIAL     *mat,
               INT           kintyp,   /* typ of kinematic formulation */
               INT           kstep,
               INT           init);
void s9_stress_reduce(FIELD     *actfield,
                      PARTITION *actpart,
                      INTRA     *actintra,
                      INT        kstep);
/*----------------------------------------------------------------------*
 |  s9_tforce.c                                          m.gee 12/01    |
 *----------------------------------------------------------------------*/
void s9_tforce(DOUBLE   **force,
               INT        ngauss,
               DOUBLE  ***akov,
               INT        klay,    /* actual kin layer */
               ELEMENT   *ele);
void s9_tstress(DOUBLE   **force,
                DOUBLE    *stress,
                INT        ngauss,
                DOUBLE   **akov,
                ELEMENT   *ele);
void s9_tettr(DOUBLE x[3][3], DOUBLE **a, DOUBLE b[3][3]);
/*----------------------------------------------------------------------*
 |  s9_tfte.c         m.gee 02/02     modified for shell9   by sh02/03  |
 *----------------------------------------------------------------------*/
void s9_tfte(DOUBLE  **force,
             INT       ngauss,
             DOUBLE   *stress,
             DOUBLE  **gkov,
             DOUBLE ***akon,
             DOUBLE ***amkov,
             DOUBLE    h,
             DOUBLE    e3,
             DOUBLE   *klayhgt,   /* hight of kin layer in % of total thickness of shell */
             DOUBLE   *mlayhgt,   /* hight of mat layer in % of adjacent kin layer */
             INT       klay,      /* actual kin layer */
             INT       mlay,      /* actual mat layer of this kin layer */
             DOUBLE    fact,
             DOUBLE    detsm,
             DOUBLE    detsr);
/*----------------------------------------------------------------------*
 |  s9_tmtr.c                                               sh 10/02    |
 *----------------------------------------------------------------------*/
void s9_tmtr(DOUBLE   **x,
             DOUBLE  ***a3,
             DOUBLE     e3,
             DOUBLE   **gkov,
             DOUBLE   **gkon,
             DOUBLE   **gmkov,
             DOUBLE   **gmkon,
             DOUBLE    *det,
             DOUBLE    *funct,
             DOUBLE   **deriv,
             INT        iel,
             DOUBLE  ***akov,
             DOUBLE  ***a3kvp,
             DOUBLE     h,           /* total thickness of this element */
             DOUBLE    *klayhgt,     /* hight of kin layer in % of total thickness of shell */
             DOUBLE    *mlayhgt,     /* hight of mat layer in % of adjacent kin layer */
             INT        num_klay,    /* number of kin layers to this element */  
             INT        num_mlay,    /* number of mat layers to this kin layer */
             INT        klay,        /* actual kin layer */
             INT        mlay,        /* actual mat layer of this kin layer */
             DOUBLE     condfac);
/*----------------------------------------------------------------------*
 |  s9_tvbo.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s9_tvbo(DOUBLE      e1,
             DOUBLE      e2,
             DOUBLE    **bop,
             DOUBLE     *funct,
             DOUBLE    **deriv,
             INT         iel,
             INT         numdf,
             DOUBLE   ***akov,
             DOUBLE   ***a3kvp,
             INT         num_klay,    /* number of kin layers to this element */  
             INT         klay,        /* actual kin layer */
             DOUBLE      condfac,
             INT         nsansq);
/*----------------------------------------------------------------------*
 |  s9_tvhe.c          m.gee 02/02     modified for shell9     sh 10/02 |
 *----------------------------------------------------------------------*/
void s9_tvhe(DOUBLE  **gmkovr,
             DOUBLE  **gmkovc,
             DOUBLE  **gmkonr,
             DOUBLE  **gmkonc,
             DOUBLE  **gkovr,
             DOUBLE  **gkovc,
             DOUBLE   *detr,
             DOUBLE   *detc,
             DOUBLE ***amkovc,
             DOUBLE ***amkovr,
             DOUBLE ***akovc,
             DOUBLE ***akovr,
             DOUBLE ***a3kvpc,
             DOUBLE ***a3kvpr,
             DOUBLE    e3,
             INT       kintyp,    /* typ of kinematic formulation */
             DOUBLE    h,         /* total thickness of this element */
             DOUBLE   *klayhgt,   /* hight of kin layer in % of total thickness of shell */
             DOUBLE   *mlayhgt,   /* hight of mat layer in % of adjacent kin layer */
             INT       num_klay,  /* number of kin layers to this element */
             INT       num_mlay,  /* number of mat layers to this kin layer */
             INT       klay,      /* actual kin layer */
             INT       mlay,      /* actual mat layer of this kin layer */
             DOUBLE    condfac);
void s9_tvhe_lin(DOUBLE  **gmkovr,
                 DOUBLE  **gmkovc,
                 DOUBLE  **gmkonr,
                 DOUBLE  **gmkonc,
                 DOUBLE  **gkovr,
                 DOUBLE  **gkovc,
                 DOUBLE   *detr,
                 DOUBLE   *detc);
/*----------------------------------------------------------------------*
 |  s9_tvma.c       m.gee 11/01    modified for shell9   by   sh 10/02  |
 *----------------------------------------------------------------------*/
void s9_tvma(DOUBLE   **D,
             DOUBLE   **C,
             DOUBLE    *stress,
             DOUBLE    *stress_r,
             DOUBLE     e3,
             DOUBLE     fact,
             DOUBLE     h,           /* total thickness of this element */
             DOUBLE    *klayhgt,     /* hight of kin layer in % of total thickness of shell */
             DOUBLE    *mlayhgt,     /* hight of mat layer in % of adjacent kin layer */
             INT        num_klay,    /* number of kin layers to this element */  
             INT        num_mlay,    /* number of mat layers to this kin layer */
             INT        klay,        /* actual kin layer */
             INT        mlay,        /* actual mat layer of this kin layer */
             DOUBLE     condfac);
/*----------------------------------------------------------------------*
 |  s9_tvmr.c                                                sh 10/02   |
 |                                                                      |
 | callculates all the metrics at the middle surface of each kinematic  |
 | layer                                                                |
 *----------------------------------------------------------------------*/
void s9_tvmr(DOUBLE    **x,
             DOUBLE   ***a3,
             DOUBLE   ***akov,
             DOUBLE   ***akon,
             DOUBLE   ***amkov,
             DOUBLE   ***amkon,
             DOUBLE    **akovh,
             DOUBLE    **akonh,
             DOUBLE    **amkovh,
             DOUBLE    **amkonh,
             DOUBLE     *det,
             DOUBLE     *funct,
             DOUBLE    **deriv,
             INT         iel,
             DOUBLE   ***a3kvp,
             INT         num_klay);
/*----------------------------------------------------------------------*
 |  s9_vthv.c          m.gee 02/02     modified for shell9     sh 02/03 |
 *----------------------------------------------------------------------*/
void s9_vthv(DOUBLE **gmkovc,
             DOUBLE **gmkonc,
             DOUBLE  *epsh,
             DOUBLE  *detc,
             DOUBLE   e3,
             DOUBLE   h,         /* total thickness of this element */
             DOUBLE  *klayhgt,   /* hight of kin layer in % of total thickness of shell */
             DOUBLE  *mlayhgt,   /* hight of mat layer in % of adjacent kin layer */
             INT      num_klay,  /* number of kin layers to this element */
             INT      num_mlay,  /* number of mat layers to this kin layer */
             INT      klay,      /* actual kin layer */
             INT      mlay,      /* actual mat layer of this kin layer */
             DOUBLE   condfac);
/*----------------------------------------------------------------------*
 |  s9_tvkg.c         m.gee 02/02     modified for shell9   by sh02/03  |
 *----------------------------------------------------------------------*/
void s9_tvkg(DOUBLE **estif,
             DOUBLE  *stress_r,
             DOUBLE  *funct,
             DOUBLE **deriv,
             INT      numdf,
             INT      iel,
             DOUBLE   weight,
             INT      klay,        /* actual kin layer */
             INT      num_klay);   /* number of kin layers to this element */  
/*----------------------------------------------------------------------*
 |  s9_xint.c                                            m.gee 02/02    |
 *----------------------------------------------------------------------*/
void s9_xint(DOUBLE *result, DOUBLE *values, DOUBLE *funct, INT iel);
DOUBLE s9con(DOUBLE    e3,       /* zeta_kl */
             INT       numlay,   /* number of kinematic layers */
             INT       ilay,     /* kin layer where point is located */
             INT       jlay,     /* kin layer over which is looped */   
             DOUBLE    condfac);
DOUBLE s9notr(INT       numlay,   /* number of kinematic layers */
              INT       ilay,     /* kin layer where point is located */
              INT       jlay);    /* kin layer over which is looped */   
DOUBLE s9ksi(INT       numlay,   /* number of kinematic layers */
             INT       ilay,     /* kin layer where point is located */
             INT       jlay,     /* kin layer over which is looped */  
             DOUBLE    condfac); 
/*----------------------------------------------------------------------*/
#endif /*D_SHELL9*/
/*! @} (documentation module close)*/
