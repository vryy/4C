/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Michael Gee
            gee@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/gee/
            0711 - 685-6572
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SHELL8
/*----------------------------------------------------------------------*
 | shell8                                                 m.gee 8/00    |
 *----------------------------------------------------------------------*/
typedef struct _SHELL8
{
/*----------------------------- general variables needed by the element */
struct _ARRAY    a3ref;
struct _ARRAY    thick_node;
DOUBLE           thick;
DOUBLE           sdc;
INT              nGP[3];
INT              nGP_tri;
/*---------------------------------------------- array of forces at GPs */
enum
    {
    s8_xyz,
    s8_rst,
    s8_rst_ortho
    }            forcetyp;
struct _ARRAY4D  forces;
/*struct _ARRAY    energy;*/
/*-------------------------------------------- array of internal forces */
struct _ARRAY    intforce;
/*-------------------------------------------- variables needed for eas */
INT              eas[5];
INT              nhyb;
struct _ARRAY    alfa;
struct _ARRAY    oldalfa;
struct _ARRAY    Dtildinv;
struct _ARRAY    Lt;
struct _ARRAY    Rtilde;
/*-------------------------------------------- variables needed for ans */
INT              ans;
/*------------------------------- variables needed for material history */
struct _ARRAY4D *his1;
struct _ARRAY4D *his2;
} SHELL8;

/*----------------------------------------------------------------------*
 | shell8 data                                            m.gee 6/01    |
 *----------------------------------------------------------------------*/
typedef struct _S8_DATA
{
DOUBLE        xgpr[3];
DOUBLE        wgtr[3];

DOUBLE        xgps[3];
DOUBLE        wgts[3];

DOUBLE        xgpt[3];
DOUBLE        wgtt[3];
} S8_DATA;






/*----------------------------------------------------------------------*
 |  PROTOTYPES OF ALL SHELL8 ROUTINES                    m.gee 11/01    |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 |  s8_a3.c                                              m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8a3(ELEMENT   *ele,
          S8_DATA   *data,/* is' hier ueberfluessig, rausschmeissen!*/
          INT        option);
void s8a3ref_extern(DOUBLE   *funct,
                       DOUBLE  **deriv,
                       DOUBLE   *thick,
                       DOUBLE  **a3ref,
                       ELEMENT  *ele);
void s8averdir(DOUBLE **dir_list, INT numa3, DOUBLE *a3);
/*----------------------------------------------------------------------*
 |  s8_btdb.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8_BtDB(DOUBLE **estif, DOUBLE **bop, DOUBLE **D, INT iel,
                INT numdf, DOUBLE weight, DOUBLE **work);
/*----------------------------------------------------------------------*
 |  s8_eas.c                                             m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8_eas(const INT    nhyb,
               const DOUBLE e1,
               const DOUBLE e2,
               const INT    iel,
               const INT   *eas,
               DOUBLE     **P);
void s8_transeas(DOUBLE      **P,
                    DOUBLE      **transP,
                    DOUBLE      **T,
                    DOUBLE      **akovr,
                    DOUBLE      **akonr0,
                    DOUBLE        detr,
                    DOUBLE        detr0,
                    INT           nhyb);
/*----------------------------------------------------------------------*
 |  s8_eps.c                                             m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8_eps(DOUBLE *strain,DOUBLE **gmkovc, DOUBLE **gmkovr);
/*----------------------------------------------------------------------*
 |  s8_funcderiv.c                                       m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8_funct_deriv(DOUBLE     *funct,
                       DOUBLE    **deriv,
                       DOUBLE      r,
                       DOUBLE      s,
                       DIS_TYP     typ,
                       INT         option);
/*----------------------------------------------------------------------*
 |  s8_init.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8init(FIELD *actfield, PARTITION *actpart);
/*----------------------------------------------------------------------*
 |  s8_inpele.c                                          m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8inp(ELEMENT *ele);
/*----------------------------------------------------------------------*
 |  s8_intg.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8intg(const ELEMENT   *ele,
               S8_DATA         *data,
               INT              option);
/*----------------------------------------------------------------------*
 |  s8_jaco.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8jaco(DOUBLE    *funct,
               DOUBLE   **deriv,
               DOUBLE   **x,
               DOUBLE   **xjm,
               DOUBLE    *hte,
               DOUBLE   **a3ref,
               DOUBLE     e3,
               INT        iel,
               DOUBLE    *det,
               DOUBLE    *deta,
               INT        init);
/*----------------------------------------------------------------------*
 |  s8_load1.c                                           m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8eleload(ELEMENT  *ele,
                  S8_DATA  *data,
                  MATERIAL *mat,
                  DOUBLE   *global_vec,
                  INT       init);
void s8loadGP(ELEMENT    *ele,
                DOUBLE    **eload,
                DOUBLE      hhi,
                DOUBLE      wgt,
                DOUBLE    **xjm,
                DOUBLE     *funct,
                DOUBLE    **deriv,
                INT         iel,
                DOUBLE      xi,
                DOUBLE      yi,
                DOUBLE      zi);
void s8fsiload(ELEMENT  *ele,
               S8_DATA  *data,
               MATERIAL *mat,
               DOUBLE	*loadvec,
               INT	 init);
/*----------------------------------------------------------------------*
 |  s8_loccoordnode.c                                    m.gee 11/01    |
 *----------------------------------------------------------------------*/
DOUBLE s8_local_coord_node(INT node, INT flag, enum _DIS_TYP typ);
/*----------------------------------------------------------------------*
 |  s8_main.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void shell8(FIELD      *actfield,
            PARTITION  *actpart,
            INTRA      *actintra,
            ELEMENT    *ele,
            ARRAY      *estif_global,
            ARRAY      *emass_global,
            ARRAY      *intforce_global,
            CALC_ACTION *action,
            CONTAINER  *container);
/*----------------------------------------------------------------------*
 |  s8_mat_linel.c                                       m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8_mat_linel(STVENANT *mat, DOUBLE **g, DOUBLE **CC);
void s8_mat_stress1(DOUBLE *stress, DOUBLE *strain, DOUBLE **C);
/*----------------------------------------------------------------------*
 |  s8_mtr.c                                             m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8mtr(DOUBLE   **x,
              DOUBLE   **a3,
              DOUBLE     e3,
              DOUBLE   **gkov,
              DOUBLE   **gkon,
              DOUBLE   **gmkov,
              DOUBLE   **gmkon,
              DOUBLE    *det,
              DOUBLE    *funct,
              DOUBLE   **deriv,
              DOUBLE    *hte,
              INT        iel,
              DOUBLE     condfac,
              char       string[]);
/*----------------------------------------------------------------------*
 |  s8_static_ke.c                                       m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8static_ke(ELEMENT   *ele,
                    S8_DATA   *data,
                    MATERIAL  *mat,
                    ARRAY     *estif_global,
                    DOUBLE    *force,
                    INT        iforce,
                    INT        kstep,
                    INT        init);
/*----------------------------------------------------------------------*
 |  s8_static_keug.c                                       m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8static_keug(ELEMENT   *ele,                         /* the element structure */
                   S8_DATA   *data,                       /* element integration data */
                   MATERIAL  *mat,                        /* the material structure */
                   ARRAY     *estif_global,               /* element stiffness matrix (NOT initialized!) */
                   ARRAY     *emass_global,               /* element mass matrix      (NOT initialized!) */
                   DOUBLE    *force,                      /* global vector for internal forces (initialized!) */
                   INT        kstep,                      /* actual step in nonlinear analysis */
                   INT        init);                      /* init=1 -> init phase / init=0 -> calc. phase / init=-1 -> uninit phase */
/*----------------------------------------------------------------------*
 |  s8_tmat.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8_tmat(ELEMENT    *ele,
             MATERIAL   *mat,
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
             DOUBLE      detc,
             DOUBLE      detr,
             DOUBLE      e3,
             INT         option,
             INT         ngauss);
void s8_getdensity(MATERIAL   *mat, DOUBLE *density);
/*----------------------------------------------------------------------*
 |  s8_tmtr.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8_tmtr(DOUBLE   **x,
                DOUBLE   **a3,
                DOUBLE     e3,
                DOUBLE   **gkov,
                DOUBLE   **gkon,
                DOUBLE   **gmkov,
                DOUBLE   **gmkon,
                DOUBLE    *det,
                DOUBLE    *funct,
                DOUBLE   **deriv,
                INT        iel,
                DOUBLE     condfac,
                INT        flag);
/*----------------------------------------------------------------------*
 |  s8_tvbo.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8_tvbo(DOUBLE      e1,
                DOUBLE      e2,
                DOUBLE    **bop,
                DOUBLE     *funct,
                DOUBLE    **deriv,
                INT         iel,
                INT         numdf,
                DOUBLE    **akov,
                DOUBLE    **a3kvp,
                INT         nansq);
/*----------------------------------------------------------------------*
 |  s8_tvhe_linear.c                                     m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8_tvhe_linear(DOUBLE **gmkovr,
                       DOUBLE **gmkovc,
                       DOUBLE **gmkonr,
                       DOUBLE **gmkonc,
                       DOUBLE **gkovr,
                       DOUBLE **gkovc,
                       DOUBLE  *detr,
                       DOUBLE  *detc);
/*----------------------------------------------------------------------*
 |  s8_tvma.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8_tvma(DOUBLE **D, DOUBLE **C, DOUBLE *stress, DOUBLE *stress_r,
                DOUBLE e3, DOUBLE fact, DOUBLE condfac);
/*----------------------------------------------------------------------*
 |  s8_tvmr.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8_tvmr(DOUBLE   **x,
                DOUBLE   **a3,
                DOUBLE   **akov,
                DOUBLE   **akon,
                DOUBLE   **amkov,
                DOUBLE   **amkon,
                DOUBLE    *det,
                DOUBLE    *funct,
                DOUBLE   **deriv,
                INT        iel,
                DOUBLE   **a3kvp,
                INT        flag);
/*----------------------------------------------------------------------*
 |  s8_tforce.c                                          m.gee 12/01    |
 *----------------------------------------------------------------------*/
void s8_tforce(DOUBLE **force,
             INT      ngauss,
             DOUBLE **akov,
             DOUBLE **akon,
             ELEMENT *ele);
void s8_tettr(DOUBLE x[3][3], DOUBLE **a, DOUBLE b[3][3]);
/*----------------------------------------------------------------------*
 |  s8_stress.c                                          m.gee 02/02    |
 *----------------------------------------------------------------------*/
void s8_stress(ELEMENT      *ele,
               S8_DATA      *data,
               MATERIAL     *mat,
               INT           kstep,
               INT           init);
void s8_stress_reduce(FIELD     *actfield,
                      PARTITION *actpart,
                      INTRA     *actintra,
                      INT        kstep);
void s8_strains_res(DOUBLE **akovr, DOUBLE **akovc, DOUBLE **a3kvpr, DOUBLE **a3kvpc,
                   DOUBLE hh, DOUBLE *strains);
/*----------------------------------------------------------------------*
 |  s8_vthv.c                                            m.gee 02/02    |
 *----------------------------------------------------------------------*/
void s8_vthv(DOUBLE **gmkovc,
             DOUBLE **gmkonc,
             DOUBLE  *epsh,
             DOUBLE  *detc,
             DOUBLE   e3,
             DOUBLE   condfac);
/*----------------------------------------------------------------------*
 |  s8_tvhe.c                                            m.gee 02/02    |
 *----------------------------------------------------------------------*/
void s8_tvhe(DOUBLE **gmkovr,
             DOUBLE **gmkovc,
             DOUBLE **gmkonr,
             DOUBLE **gmkonc,
             DOUBLE **gkovr,
             DOUBLE **gkovc,
             DOUBLE  *detr,
             DOUBLE  *detc,
             DOUBLE **amkovc,
             DOUBLE **amkovr,
             DOUBLE **akovc,
             DOUBLE **akovr,
             DOUBLE **a3kvpc,
             DOUBLE **a3kvpr,
             DOUBLE   e3,
             DOUBLE   condfac);
void s8_tvhe_lin(DOUBLE **gmkovr,
                 DOUBLE **gmkovc,
                 DOUBLE **gmkonr,
                 DOUBLE **gmkonc,
                 DOUBLE **gkovr,
                 DOUBLE **gkovc,
                 DOUBLE  *detr,
                 DOUBLE  *detc,
                 DOUBLE **amkovc,
                 DOUBLE **amkovr,
                 DOUBLE **akovc,
                 DOUBLE **akovr,
                 DOUBLE **a3kvpc,
                 DOUBLE **a3kvpr,
                 DOUBLE   e3);
/*----------------------------------------------------------------------*
 |  s8_tvkg.c                                            m.gee 02/02    |
 *----------------------------------------------------------------------*/
void s8_tvkg(DOUBLE **estif,
             DOUBLE  *stress_r,
             DOUBLE  *funct,
             DOUBLE **deriv,
             INT      numdf,
             INT      iel,
             DOUBLE   weight,
             DOUBLE   e1,
             DOUBLE   e2);
/*----------------------------------------------------------------------*
 |  s8_intforce.c                                        m.gee 02/02    |
 *----------------------------------------------------------------------*/
void s8_intforce(DOUBLE *intforce, DOUBLE *stress_r, DOUBLE **bop,
                 INT iel, INT numdf, INT nstress_r, DOUBLE weight);
/*----------------------------------------------------------------------*
 |  s8_cal_dyn.c                                         m.gee 02/02    |
 *----------------------------------------------------------------------*/
void s8_tmas(DOUBLE *funct, DOUBLE *thick, DOUBLE **emass, INT iel, INT numdf,
             DOUBLE facv, DOUBLE facw, DOUBLE facvw);
/*----------------------------------------------------------------------*
 |  s8_xint.c                                            m.gee 02/02    |
 *----------------------------------------------------------------------*/
void s8_xint(DOUBLE *result, DOUBLE *values, DOUBLE *funct, INT iel);
/*----------------------------------------------------------------------*
 |  s8_tfte.c                                            m.gee 02/02    |
 *----------------------------------------------------------------------*/
void s8_tfte(DOUBLE **force,
             INT      ngauss,
             DOUBLE  *stress,
             DOUBLE **gkov,
             DOUBLE **akon,
             DOUBLE **gmkov,
             DOUBLE **gmkon,
             DOUBLE **amkov,
             DOUBLE **amkon,
             DOUBLE   h,
             DOUBLE   e3,
             DOUBLE   fact,
             DOUBLE   detsm,
             DOUBLE   detsr);
/*----------------------------------------------------------------------*
 |  s8_ans.c                                            m.gee 03/02    |
 *----------------------------------------------------------------------*/
void s8_ans_colloqpoints(INT nsansq,INT iel,INT ans,DIS_TYP distyp,
                        DOUBLE xr1[],DOUBLE xs1[],DOUBLE xr2[],DOUBLE xs2[],
                        DOUBLE *funct1q[],DOUBLE **deriv1q[],
                        DOUBLE *funct2q[],DOUBLE **deriv2q[],
                        DOUBLE **xrefe,DOUBLE **a3r,DOUBLE **xcure,DOUBLE **a3c,
                        DOUBLE **akovr1q[] ,DOUBLE **akonr1q[],
                        DOUBLE **amkovr1q[],DOUBLE **amkonr1q[],
                        DOUBLE **a3kvpr1q[],
                        DOUBLE **akovc1q[] ,DOUBLE **akonc1q[],
                        DOUBLE **amkovc1q[],DOUBLE **amkonc1q[],
                        DOUBLE **a3kvpc1q[],
                        DOUBLE **akovr2q[] ,DOUBLE **akonr2q[],
                        DOUBLE **amkovr2q[],DOUBLE **amkonr2q[],
                        DOUBLE **a3kvpr2q[],
                        DOUBLE **akovc2q[] ,DOUBLE **akonc2q[],
                        DOUBLE **amkovc2q[],DOUBLE **amkonc2q[],
                        DOUBLE **a3kvpc2q[],
                        DOUBLE *detr, DOUBLE *detc);
void s8_ans_colloqcoords(DOUBLE xqr1[], DOUBLE xqs1[],
                        DOUBLE xqr2[], DOUBLE xqs2[],
                        INT iel, INT ans);
void s8_ansq_funct(DOUBLE frq[], DOUBLE fsq[], DOUBLE r, DOUBLE s,
                  INT iel, INT nsansq);
void s8_ans_bbar_q(DOUBLE **bop, DOUBLE frq[], DOUBLE fsq[],
                  DOUBLE  *funct1q[],  DOUBLE  *funct2q[],
                  DOUBLE **deriv1q[],  DOUBLE **deriv2q[],
                  DOUBLE **akovc1q[],  DOUBLE **akovc2q[],
                  DOUBLE **a3kvpc1q[], DOUBLE **a3kvpc2q[],
                  INT iel, INT numdf, INT nsansq);
void s8_ans_tvhe_q(DOUBLE **gmkovr,DOUBLE **gmkovc,DOUBLE **gmkonr,DOUBLE **gmkonc,
                  DOUBLE **gkovr,DOUBLE **gkovc,DOUBLE **amkovc,DOUBLE **amkovr,
                  DOUBLE **akovc,DOUBLE **akovr,DOUBLE **a3kvpc,DOUBLE **a3kvpr,
                  DOUBLE *detr,   DOUBLE *detc,
                  DOUBLE **amkovr1q[], DOUBLE **amkovc1q[],
                  DOUBLE **akovr1q[] , DOUBLE **akovc1q[] ,
                  DOUBLE **a3kvpr1q[], DOUBLE **a3kvpc1q[],
                  DOUBLE **amkovr2q[], DOUBLE **amkovc2q[],
                  DOUBLE **akovr2q[] , DOUBLE **akovc2q[] ,
                  DOUBLE **a3kvpr2q[], DOUBLE **a3kvpc2q[],
                  DOUBLE frq[], DOUBLE fsq[], DOUBLE e3, INT nansq, INT iel,
                  DOUBLE condfac);
void s8_ans_tvkg(DOUBLE **estif,DOUBLE *stress_r,DOUBLE *funct,DOUBLE **deriv,
                INT numdf,INT iel,DOUBLE weight,DOUBLE e1,DOUBLE e2,
                DOUBLE frq[], DOUBLE fsq[],DOUBLE *funct1q[],DOUBLE  *funct2q[],
                DOUBLE **deriv1q[], DOUBLE **deriv2q[],INT ansq, INT nsansq);
/*----------------------------------------------------------------------*
 |  s8_restart.c                                         m.gee 05/02    |
 *----------------------------------------------------------------------*/
void s8_write_restart(ELEMENT *actele, INT nhandle, long int *handles);
void s8_read_restart(ELEMENT *actele, INT nhandle, long int *handles);
/*----------------------------------------------------------------------*
 |  s8_mat_ogden.c                                       m.gee 06/03    |
 *----------------------------------------------------------------------*/
void s8_mat_ogden_coupled(COMPOGDEN *mat, DOUBLE *stress_cart, DOUBLE C_cart[3][3][3][3],
                          DOUBLE **gkonr,DOUBLE **gmkovc);
void s8_ogden_Ccart(DOUBLE C[3][3][3][3], DOUBLE C_cart[3][3][3][3], DOUBLE N[3][3]);
void s8_ogden_cartPK2(DOUBLE PK2[3][3], DOUBLE PK2main[3], DOUBLE N[3][3]);
void s8_ogden_principal_CG(DOUBLE CG[3][3], DOUBLE lambda[3], DOUBLE N[3][3]);
void s8_kov_CGcuca(DOUBLE T[3][3], const DOUBLE **gkon);
/*----------------------------------------------------------------------*
 |  s8_mat_ogden2.c                                      m.gee 06/03    |
 *----------------------------------------------------------------------*/
void s8_mat_ogden_uncoupled(COMPOGDEN *mat, DOUBLE *stress_cart, DOUBLE *strain, DOUBLE C_cart[3][3][3][3],
                            DOUBLE **gkovr, DOUBLE **gkonr, DOUBLE **gkovc, DOUBLE **gkonc,
                            DOUBLE **gmkovr,DOUBLE **gmkonr, DOUBLE **gmkovc, DOUBLE **gmkonc);
/*----------------------------------------------------------------------*
 |  s8_mat_ogden3.c                                      m.gee 06/03    |
 *----------------------------------------------------------------------*/
void s8_mat_ogden_uncoupled2(COMPOGDEN *mat, DOUBLE *stress_cart, DOUBLE C_cart[3][3][3][3],
                            DOUBLE **gkonr, DOUBLE **gmkovc);
/*----------------------------------------------------------------------*
 |  s8_mattransform.c                                    m.gee 06/03    |
 *----------------------------------------------------------------------*/
void s8_kov_cuca(DOUBLE *t, const DOUBLE **gkon);
void s8_kon_cacu(DOUBLE *t, const DOUBLE **gkon);
void s8_kov_cacu(DOUBLE *t, const DOUBLE **gkov);
void s8_4kon_cacu(DOUBLE Ccart[][3][3][3], const DOUBLE **gkon);
void s8_c4_to_C2(const DOUBLE C[][3][3][3],DOUBLE **CC);
void s8_mat_linel_cart(STVENANT *mat,DOUBLE C[][3][3][3],DOUBLE **CC,DOUBLE *strain);
/*----------------------------------------------------------------------*
 |  s8_static_mass.c                                     m.gee 06/03    |
 *----------------------------------------------------------------------*/
void s8static_mass(ELEMENT   *ele,                        /* the element structure */
                   S8_DATA   *data,                       /* element integration data */
                   MATERIAL  *mat,                        /* the material structure */
                   ARRAY     *estif_global,               /* element stiffness matrix (NOT initialized!) */
                   ARRAY     *emass_global,               /* element mass matrix      (NOT initialized!) */
                   DOUBLE    *force,                      /* for internal forces (initialized!) */
                   INT        kstep,                      /* actual step in nonlinear analysis */
                   INT        init);                      /* init=1 -> init phase / init=0 -> calc. phase / init=-1 -> uninit phase */
/*----------------------------------------------------------------------*
 |  s8_mat_linel.c                                     m.gee 06/03    |
 *----------------------------------------------------------------------*/
void s8_mat_neohooke(NEO_HOOKE *mat,
                     DOUBLE    *stress,
                     DOUBLE   **CC,
                     DOUBLE   **gmkonr,
                     DOUBLE   **gmkonc,
                     DOUBLE     detr,
                     DOUBLE     detc);
/*----------------------------------------------------------------------*
 |  s8_fserv.f                                           m.gee 06/03    |
 *----------------------------------------------------------------------*/
void s8jacb(DOUBLE *A, DOUBLE *V);




#endif
