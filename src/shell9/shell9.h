/*!----------------------------------------------------------------------
\file
\brief headerfile for multilayer shell element, containing structures and
       prototypes

*----------------------------------------------------------------------*/
#ifdef D_SHELL9

/*! 
\addtogroup SHELL9 
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief shell9 data

<pre>                                                            sh 03/03
This structure contains the coordinates and weights for numerical
integration of a shell9 element
</pre>

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

/*!----------------------------------------------------------------------
\brief shell9 gausspoint working array

<pre>                                                            sh 03/03
This structure contains the working array for g.p. values 
</pre>

*----------------------------------------------------------------------*/
typedef struct _S9_IP_WA
{
/* mises ...*/
DOUBLE      sig[6]; /*!<  global stresses                             */
DOUBLE      eps[6]; /*!<  global strains                              */
DOUBLE       qn[6]; /*!< 'backstress' vec                             */
DOUBLE       epstn; /*!<  equivalent strain                           */
INT          yip;   /*!<  stress state: 1=elastic 2=plastic           */
/* hoff ...*/
DOUBLE   dkappa[6]; /*!< */
DOUBLE   gamma[6];  /*!< */
DOUBLE   rkappa[9]; /*!< */
DOUBLE   dhard;     /*!< */
} S9_IP_WA;


/*!----------------------------------------------------------------------
\brief shell9 working array

<pre>                                                            sh 03/03
This structure contains the working array for element values 
</pre>

*----------------------------------------------------------------------*/
typedef struct _S9_ELE_WA
{
S9_IP_WA      *ipwa;    /*!< working array for integration points       */
} S9_ELE_WA;




/*!----------------------------------------------------------------------
\brief shell9 element

<pre>                                                           sh 10/02
This structure contains all specific information for a shell9 element
</pre>

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
struct _ARRAY    Dtildinv;
struct _ARRAY    L;
struct _ARRAY    Rtilde;
/*-------------------------------------------- variables needed for ans */
INT              ans;
/*--------------------------------- variables needed nonlinear material */
S9_ELE_WA       *elewa;                      /*!< element working array */
} SHELL9;




/*----------------------------------------------------------------------*
 |  PROTOTYPES OF ALL SHELL9 ROUTINES      (shell8: m.gee 11/01)        |
 |                modified for shell9                      by sh 10/02  |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 |  s9_Ccacu_sym.c                                          sh 09/03    |
 *----------------------------------------------------------------------*/
void s9_Ccacu_sym (DOUBLE **A, DOUBLE **C);
/*----------------------------------------------------------------------*
 |  s9_Ccacu_unsym.c                                        sh 09/03    |
 *----------------------------------------------------------------------*/
void s9_Ccacu_unsym (DOUBLE **A, DOUBLE **C);
/*----------------------------------------------------------------------*
 |  s9_EStrans.c                                            sh 09/03    |
 *----------------------------------------------------------------------*/
void s9_Ecacu (DOUBLE E[6], DOUBLE **gkov);
void s9_Ecuca (DOUBLE E[6], DOUBLE **gkon);
void s9_Scacu (DOUBLE S[6], DOUBLE **gkon);
void s9_Scuca (DOUBLE S[6], DOUBLE **gkov);
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
                         INT       num_klay);
void s9_ans_colloqcoords(DOUBLE xqr1[], DOUBLE xqs1[],
                         DOUBLE xqr2[], DOUBLE xqs2[],
                         INT    iel,    INT    ans);
void s9_ansq_funct(DOUBLE frq[], DOUBLE fsq[], DOUBLE r, DOUBLE s,
                   INT    iel);
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
                   DOUBLE ***amkovc,     DOUBLE ***amkovr,
                   DOUBLE ***akovc,      DOUBLE ***akovr,     
                   DOUBLE ***a3kvpc,     DOUBLE ***a3kvpr,
                   DOUBLE   *detr,       DOUBLE   *detc,
                   DOUBLE ***amkovr1q[], DOUBLE ***amkovc1q[], 
                   DOUBLE ***amkovr2q[], DOUBLE ***amkovc2q[], 
                   DOUBLE    frq[],      DOUBLE    fsq[], 
                   DOUBLE    e3,         INT       nansq, 
                   DOUBLE    h,           /* total thickness of this element */
                   DOUBLE   *klayhgt,     /* hight of kin layer in % of total thickness of shell */
                   DOUBLE   *mlayhgt,     /* hight of mat layer in % of adjacent kin layer */
                   INT       num_klay,    /* number of kin layers to this element */
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
                 DOUBLE      stress[6],
                 DOUBLE      strain[6],
                 DOUBLE    **C,
                 DOUBLE    **gmkovc,
                 DOUBLE    **gmkovr,
                 DOUBLE    **gmkonr,
                 DOUBLE    **gkovr,
                 DOUBLE    **gkonr,
                 INT         rot_axis,
                 DOUBLE      phi,
                 INT         ip,
                 INT         actlay,
                 INT         istore, /* controls storing of new stresses to wa */
                 INT         newval);/* controls evaluation of new stresses    */
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
void s9_eps(DOUBLE strain[6],DOUBLE **gmkovc, DOUBLE **gmkovr);
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
void s9intg_str(S9_DATA         *data,
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
            DOUBLE    *klayhgt);     /* hight of kin layer in % of total thickness of shell */
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
              INT         iel,
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
void s9_mat_stress1(DOUBLE stress[6], DOUBLE strain[6], DOUBLE **C);
/*----------------------------------------------------------------------*
 |  s9_mat_plast_hoff.c                                     sh 03/03    |
 *----------------------------------------------------------------------*/
void s9_mat_plast_hoff(
                 PL_HOFF   *mat,       /*!< material properties          */
                 ELEMENT   *ele,       /*!< actual element               */
                 INT        ip,        /*!< integration point Id         */
                 INT        actlay,    /*!< actual layer                 */
                 DOUBLE     stress[6], /*!< vector of stresses [11,22,33,12,23,13]  */
                 DOUBLE     strain[6], /*!< vector of strains  [11,22,33,12,23,13]  */
                 DOUBLE   **d,         /*!< constitutive matrix          */
                 INT        istore,    /*!< controls storing of stresses */
                 INT        newval);   /*!< controls eval. of stresses   */
/*----------------------------------------------------------------------*
 |  s9_mat_plast_mises.c                                    sh 03/03    |
 *----------------------------------------------------------------------*/
void s9_mat_plast_mises(
                 DOUBLE     ym,        /*!< young's modulus              */
                 DOUBLE     pv,        /*!< poisson's ratio              */
                 DOUBLE     sigy,      /*!< yield stress                 */
                 DOUBLE     hard,      /*!< hardening modulus            */
                 DOUBLE     gf,        /*!< fracture energy              */
                 DOUBLE     betah,     /*!< controls the iso/kin hard.   */
                 ELEMENT   *ele,       /*!< actual element               */
                 INT        ip,        /*!< integration point Id         */
                 INT        actlay,    /*!< actual layer                 */
                 DOUBLE     stress[6], /*!< vector of stresses [11,22,33,12,23,13]  */
                 DOUBLE     strain[6], /*!< vector of strains  [11,22,33,12,23,13]  */
                 DOUBLE   **d,         /*!< constitutive matrix          */
                 INT        istore,    /*!< controls storing of stresses */
                 INT        newval);   /*!< controls eval. of stresses   */
/*----------------------------------------------------------------------*
 |  s9_math.c                                               sh 02/03    |
 *----------------------------------------------------------------------*/
void s9_Msort_bs9(DOUBLE **M);
void s9_Msort_s9b(DOUBLE **M);
void s9_Vsort_bs9(DOUBLE vec[6]);
void s9_Vsort_s9b(DOUBLE vec[6]);
void s9shear_cor(DOUBLE **C, DOUBLE xsi); 
void s9_rot(DOUBLE phi,   INT    rot_axis, DOUBLE **gkov, 
            DOUBLE g1[3], DOUBLE g2[3],    DOUBLE g3[3]);
void s9_Teps(DOUBLE g1[3], DOUBLE g2[3], DOUBLE g3[3],DOUBLE T[6][6]);
void s9_Ecama(DOUBLE E[6], DOUBLE T[6][6]);
void s9_Smaca(DOUBLE S[6], DOUBLE T[6][6]);
void s9_Cmaca(DOUBLE **C, DOUBLE T[6][6]);
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
 |  s9_out_unsmo.c                                          sh 07/03    |
 *----------------------------------------------------------------------*/
void s9_out_gid_allcoords_unsmo(FILE *out,INT type);
void s9_out_gid_eletop_unsmo(FILE *out,INT type,INT j);
void s9_out_gid_sol_dis_unsmo(FILE *out,FIELD *actfield,INT place);
void s9_out_gid_sol_str_unsmo(FILE *out,FIELD *actfield,INT place);
/*----------------------------------------------------------------------*
 |  s9_restart.c                                         m.gee 05/02    |
 *----------------------------------------------------------------------*/
void s9_write_restart(ELEMENT *actele, INT nhandle, long int *handles, INT init);
void s9_read_restart(ELEMENT *actele, long int *handles, INT init);
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
void s9_tmtr(DOUBLE     e3,
             DOUBLE   **gkov,
             DOUBLE   **gkon,
             DOUBLE   **gmkov,
             DOUBLE   **gmkon,
             DOUBLE    *det,
             DOUBLE  ***akov,
             DOUBLE  ***a3kvp,
             DOUBLE     h,           /* total thickness of this element */
             DOUBLE    *klayhgt,     /* hight of kin layer in % of total thickness of shell */
             DOUBLE    *mlayhgt,     /* hight of mat layer in % of adjacent kin layer */
             INT        num_klay,    /* number of kin layers to this element */  
             INT        klay,        /* actual kin layer */
             INT        mlay,        /* actual mat layer of this kin layer */
             DOUBLE     condfac);
/*----------------------------------------------------------------------*
 |  s9_tvbo.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s9_tvbo(DOUBLE    **bop,
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
