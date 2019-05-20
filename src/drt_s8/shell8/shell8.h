/*----------------------------------------------------------------------------*/
/*!
\brief shell8

\level 1

\maintainer Christoph Meier

*/
/*---------------------------------------------------------------------------*/
#include "../../headers/am.h"
#include "../../headers/materials.h"

#ifdef D_SHELL8
/*----------------------------------------------------------------------*
 | shell8                                                 m.gee 8/00    |
 *----------------------------------------------------------------------*/
typedef struct _SHELL8
{
  /*----------------------------- general variables needed by the element */
  struct _ARRAY a3ref;
  struct _ARRAY thick_node;
  double thick;
  double sdc;
  int volc; /* volume constraint toggle */
  int nGP[3];
  int nGP_tri;
  /*---------------------------------------------- array of forces at GPs */
  enum
  {
    s8_xyz,
    s8_rst,
    s8_rst_ortho
  } forcetyp;
#define S8_FORCETYPE                         \
  {                                          \
    "s8_xyz", "s8_rst", "s8_rst_ortho", NULL \
  }
  struct _ARRAY4D forces;
  /*struct _ARRAY    energy;*/
  /*-------------------------------------------- array of internal forces */
  struct _ARRAY intforce;
  /*-------------------------------------------- variables needed for eas */
  int eas[5];
  int nhyb;
  struct _ARRAY alfa;
  struct _ARRAY oldalfa;
  struct _ARRAY Dtildinv;
  struct _ARRAY Lt;
  struct _ARRAY Rtilde;
  /*-------------------------------------------- variables needed for ans */
  int ans;
  /*------------------------------- variables needed for material history */
  struct _ARRAY4D *his1;
  struct _ARRAY4D *his2;
} SHELL8;

/*----------------------------------------------------------------------*
 | shell8 data                                            m.gee 6/01    |
 *----------------------------------------------------------------------*/
typedef struct _S8_DATA
{
  double xgpr[3];
  double wgtr[3];

  double xgps[3];
  double wgts[3];

  double xgpt[3];
  double wgtt[3];
} S8_DATA;


#if 0

/*----------------------------------------------------------------------*
 |  PROTOTYPES OF ALL SHELL8 ROUTINES                    m.gee 11/01    |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 |  s8_a3.c                                              m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8a3(ELEMENT   *ele,
          S8_DATA   *data,/* is' hier ueberfluessig, rausschmeissen!*/
          int        option);
void s8a3ref_extern(double   *funct,
                       double  **deriv,
                       double   *thick,
                       double  **a3ref,
                       ELEMENT  *ele);
void s8averdir(double **dir_list, int numa3, double *a3);
/*----------------------------------------------------------------------*
 |  s8_btdb.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8_BtDB(double **estif, double **bop, double **D, int iel,
                int numdf, double weight, double **work);
/*----------------------------------------------------------------------*
 |  s8_eas.c                                             m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8_eas(const int    nhyb,
               const double e1,
               const double e2,
               const int    iel,
               const int   *eas,
               double     **P);
void s8_transeas(double      **P,
                    double      **transP,
                    double      **T,
                    double      **akovr,
                    double      **akonr0,
                    double        detr,
                    double        detr0,
                    int           nhyb);
/*----------------------------------------------------------------------*
 |  s8_eps.c                                             m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8_eps(double *strain,double **gmkovc, double **gmkovr);
/*----------------------------------------------------------------------*
 |  s8_funcderiv.c                                       m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8_funct_deriv(double     *funct,
                       double    **deriv,
                       double      r,
                       double      s,
                       DIS_TYP     typ,
                       int         option);
/*----------------------------------------------------------------------*
 |  s8_init.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8init(
    FIELD              *actfield,
    PARTITION          *actpart,
    int                 disnum);


/*----------------------------------------------------------------------*
 |  s8_inpele.c                                          m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8inp(ELEMENT *ele);
/*----------------------------------------------------------------------*
 |  s8_intg.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8intg(const ELEMENT   *ele,
               S8_DATA         *data,
               int              option);
/*----------------------------------------------------------------------*
 |  s8_jaco.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8jaco(double    *funct,
               double   **deriv,
               double   **x,
               double   **xjm,
               double    *hte,
               double   **a3ref,
               double     e3,
               int        iel,
               double    *det,
               double    *deta,
               int        init);
/*----------------------------------------------------------------------*
 |  s8_load1.c                                           m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8eleload(ELEMENT  *ele,
                  S8_DATA  *data,
                  MATERIAL *mat,
                  double   *global_vec,
                  int       init);
void s8loadGP(ELEMENT    *ele,
                double    **eload,
                double      hhi,
                double      wgt,
                double    **xjm,
                double     *funct,
                double    **deriv,
                int         iel,
                double      xi,
                double      yi,
                double      zi);
void s8fsiload(ELEMENT  *ele,
               S8_DATA  *data,
               MATERIAL *mat,
               double	*loadvec,
               int	 init);
/*----------------------------------------------------------------------*
 |  s8_loccoordnode.c                                    m.gee 11/01    |
 *----------------------------------------------------------------------*/
double s8_local_coord_node(int node, int flag, enum _DIS_TYP typ);
/*----------------------------------------------------------------------*
 |  s8_main.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void shell8(FIELD      *actfield,
            PARTITION  *actpart,
            intRA      *actintra,
            ELEMENT    *ele,
            ARRAY      *estif_global,
            ARRAY      *emass_global,
            ARRAY      *intforce_global,
            CALC_ACTION *action,
            CONTAINER  *container);
/*----------------------------------------------------------------------*
 |  s8_mat_linel.c                                       m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8_mat_linel(STVENANT *mat, double **g, double **CC);
void s8_mat_stress1(double *stress, double *strain, double **C);
/*----------------------------------------------------------------------*
 |  s8_mtr.c                                             m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8mtr(double   **x,
              double   **a3,
              double     e3,
              double   **gkov,
              double   **gkon,
              double   **gmkov,
              double   **gmkon,
              double    *det,
              double    *funct,
              double   **deriv,
              double    *hte,
              int        iel,
              double     condfac,
              char       string[]);
/*----------------------------------------------------------------------*
 |  s8_static_ke.c                                       m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8static_ke(ELEMENT   *ele,
                    S8_DATA   *data,
                    MATERIAL  *mat,
                    ARRAY     *estif_global,
                    double    *force,
                    int        iforce,
                    int        kstep,
                    int        init);
/*----------------------------------------------------------------------*
 |  s8_static_keug.c                                       m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8static_keug(ELEMENT   *ele,                         /* the element structure */
                   S8_DATA   *data,                       /* element integration data */
                   MATERIAL  *mat,                        /* the material structure */
                   ARRAY     *estif_global,               /* element stiffness matrix (NOT initialized!) */
                   ARRAY     *emass_global,               /* element mass matrix      (NOT initialized!) */
                   double    *force,                      /* global vector for internal forces (initialized!) */
                   int        kstep,                      /* actual step in nonlinear analysis */
                   int        init);                      /* init=1 -> init phase / init=0 -> calc. phase / init=-1 -> uninit phase */
/*----------------------------------------------------------------------*
 |  s8_tmat.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8_tmat(ELEMENT    *ele,
             MATERIAL   *mat,
             double     *stress,
             double     *strain,
             double    **C,
             double    **gmkovc,
             double    **gmkonc,
             double    **gmkovr,
             double    **gmkonr,
             double    **gkovc,
             double    **gkonc,
             double    **gkovr,
             double    **gkonr,
             double      detc,
             double      detr,
             double      e3,
             int         option,
             int         ngauss);
void s8_getdensity(MATERIAL   *mat, double *density);
/*----------------------------------------------------------------------*
 |  s8_tmtr.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8_tmtr(double   **x,
                double   **a3,
                double     e3,
                double   **gkov,
                double   **gkon,
                double   **gmkov,
                double   **gmkon,
                double    *det,
                double    *funct,
                double   **deriv,
                int        iel,
                double     condfac,
                int        flag);
/*----------------------------------------------------------------------*
 |  s8_tvbo.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8_tvbo(double      e1,
                double      e2,
                double    **bop,
                double     *funct,
                double    **deriv,
                int         iel,
                int         numdf,
                double    **akov,
                double    **a3kvp,
                int         nansq);
/*----------------------------------------------------------------------*
 |  s8_tvhe_linear.c                                     m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8_tvhe_linear(double **gmkovr,
                       double **gmkovc,
                       double **gmkonr,
                       double **gmkonc,
                       double **gkovr,
                       double **gkovc,
                       double  *detr,
                       double  *detc);
/*----------------------------------------------------------------------*
 |  s8_tvma.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8_tvma(double **D, double **C, double *stress, double *stress_r,
                double e3, double fact, double condfac);
/*----------------------------------------------------------------------*
 |  s8_tvmr.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8_tvmr(double   **x,
                double   **a3,
                double   **akov,
                double   **akon,
                double   **amkov,
                double   **amkon,
                double    *det,
                double    *funct,
                double   **deriv,
                int        iel,
                double   **a3kvp,
                int        flag);
/*----------------------------------------------------------------------*
 |  s8_tforce.c                                          m.gee 12/01    |
 *----------------------------------------------------------------------*/
void s8_tforce(double **force,
             int      ngauss,
             double **akov,
             double **akon,
             ELEMENT *ele);
void s8_tettr(double x[3][3], double **a, double b[3][3]);
/*----------------------------------------------------------------------*
 |  s8_stress.c                                          m.gee 02/02    |
 *----------------------------------------------------------------------*/
void s8_stress(ELEMENT      *ele,
               S8_DATA      *data,
               MATERIAL     *mat,
               int           kstep,
               int           init);
void s8_stress_reduce(FIELD     *actfield,
                      PARTITION *actpart,
                      intRA     *actintra,
                      int        kstep);
void s8_strains_res(double **akovr, double **akovc, double **a3kvpr, double **a3kvpc,
                   double hh, double *strains);
/*----------------------------------------------------------------------*
 |  s8_vthv.c                                            m.gee 02/02    |
 *----------------------------------------------------------------------*/
void s8_vthv(double **gmkovc,
             double **gmkonc,
             double  *epsh,
             double  *detc,
             double   e3,
             double   condfac);
/*----------------------------------------------------------------------*
 |  s8_tvhe.c                                            m.gee 02/02    |
 *----------------------------------------------------------------------*/
void s8_tvhe(double **gmkovr,
             double **gmkovc,
             double **gmkonr,
             double **gmkonc,
             double **gkovr,
             double **gkovc,
             double  *detr,
             double  *detc,
             double **amkovc,
             double **amkovr,
             double **akovc,
             double **akovr,
             double **a3kvpc,
             double **a3kvpr,
             double   e3,
             double   condfac);
void s8_tvhe_lin(double **gmkovr,
                 double **gmkovc,
                 double **gmkonr,
                 double **gmkonc,
                 double **gkovr,
                 double **gkovc,
                 double  *detr,
                 double  *detc,
                 double **amkovc,
                 double **amkovr,
                 double **akovc,
                 double **akovr,
                 double **a3kvpc,
                 double **a3kvpr,
                 double   e3);
/*----------------------------------------------------------------------*
 |  s8_tvkg.c                                            m.gee 02/02    |
 *----------------------------------------------------------------------*/
void s8_tvkg(double **estif,
             double  *stress_r,
             double  *funct,
             double **deriv,
             int      numdf,
             int      iel,
             double   weight,
             double   e1,
             double   e2);
/*----------------------------------------------------------------------*
 |  s8_intforce.c                                        m.gee 02/02    |
 *----------------------------------------------------------------------*/
void s8_intforce(double *intforce, double *stress_r, double **bop,
                 int iel, int numdf, int nstress_r, double weight);
/*----------------------------------------------------------------------*
 |  s8_cal_dyn.c                                         m.gee 02/02    |
 *----------------------------------------------------------------------*/
void s8_tmas(double *funct, double *thick, double **emass, int iel, int numdf,
             double facv, double facw, double facvw);
/*----------------------------------------------------------------------*
 |  s8_xint.c                                            m.gee 02/02    |
 *----------------------------------------------------------------------*/
void s8_xint(double *result, double *values, double *funct, int iel);
/*----------------------------------------------------------------------*
 |  s8_tfte.c                                            m.gee 02/02    |
 *----------------------------------------------------------------------*/
void s8_tfte(double **force,
             int      ngauss,
             double  *stress,
             double **gkov,
             double **akon,
             double **gmkov,
             double **gmkon,
             double **amkov,
             double **amkon,
             double   h,
             double   e3,
             double   fact,
             double   detsm,
             double   detsr);
/*----------------------------------------------------------------------*
 |  s8_ans.c                                            m.gee 03/02    |
 *----------------------------------------------------------------------*/
void s8_ans_colloqpoints(int nsansq,int iel,int ans,DIS_TYP distyp,
                        double xr1[],double xs1[],double xr2[],double xs2[],
                        double *funct1q[],double **deriv1q[],
                        double *funct2q[],double **deriv2q[],
                        double **xrefe,double **a3r,double **xcure,double **a3c,
                        double **akovr1q[] ,double **akonr1q[],
                        double **amkovr1q[],double **amkonr1q[],
                        double **a3kvpr1q[],
                        double **akovc1q[] ,double **akonc1q[],
                        double **amkovc1q[],double **amkonc1q[],
                        double **a3kvpc1q[],
                        double **akovr2q[] ,double **akonr2q[],
                        double **amkovr2q[],double **amkonr2q[],
                        double **a3kvpr2q[],
                        double **akovc2q[] ,double **akonc2q[],
                        double **amkovc2q[],double **amkonc2q[],
                        double **a3kvpc2q[],
                        double *detr, double *detc);
void s8_ans_colloqcoords(double xqr1[], double xqs1[],
                        double xqr2[], double xqs2[],
                        int iel, int ans);
void s8_ansq_funct(double frq[], double fsq[], double r, double s,
                  int iel, int nsansq);
void s8_ans_bbar_q(double **bop, double frq[], double fsq[],
                  double  *funct1q[],  double  *funct2q[],
                  double **deriv1q[],  double **deriv2q[],
                  double **akovc1q[],  double **akovc2q[],
                  double **a3kvpc1q[], double **a3kvpc2q[],
                  int iel, int numdf, int nsansq);
void s8_ans_tvhe_q(double **gmkovr,double **gmkovc,double **gmkonr,double **gmkonc,
                  double **gkovr,double **gkovc,double **amkovc,double **amkovr,
                  double **akovc,double **akovr,double **a3kvpc,double **a3kvpr,
                  double *detr,   double *detc,
                  double **amkovr1q[], double **amkovc1q[],
                  double **akovr1q[] , double **akovc1q[] ,
                  double **a3kvpr1q[], double **a3kvpc1q[],
                  double **amkovr2q[], double **amkovc2q[],
                  double **akovr2q[] , double **akovc2q[] ,
                  double **a3kvpr2q[], double **a3kvpc2q[],
                  double frq[], double fsq[], double e3, int nansq, int iel,
                  double condfac);
void s8_ans_tvkg(double **estif,double *stress_r,double *funct,double **deriv,
                int numdf,int iel,double weight,double e1,double e2,
                double frq[], double fsq[],double *funct1q[],double  *funct2q[],
                double **deriv1q[], double **deriv2q[],int ansq, int nsansq);
/*----------------------------------------------------------------------*
 |  s8_restart.c                                         m.gee 05/02    |
 *----------------------------------------------------------------------*/
void s8_write_restart(ELEMENT *actele, int nhandle, long int *handles);
void s8_read_restart(ELEMENT *actele, int nhandle, long int *handles);
/*----------------------------------------------------------------------*
 |  s8_mat_ogden.c                                       m.gee 06/03    |
 *----------------------------------------------------------------------*/
void s8_mat_ogden_coupled(
    COMPOGDEN          *mat,
    double             *stress_cart,
    double              C_cart[3][3][3][3],
    double            **gkonr,
    double            **gmkovc
    );


void s8_ogden_Ccart(double C[3][3][3][3], double C_cart[3][3][3][3], double N[3][3]);
void s8_ogden_cartPK2(double PK2[3][3], double PK2main[3], double N[3][3]);
void s8_ogden_principal_CG(double CG[3][3], double lambda[3], double N[3][3]);
void s8_kov_CGcuca(double T[3][3], double **gkon);
/*----------------------------------------------------------------------*
 |  s8_mat_ogden2.c                                      m.gee 06/03    |
 *----------------------------------------------------------------------*/
void s8_mat_ogden_uncoupled(COMPOGDEN *mat, double *stress_cart, double *strain, double C_cart[3][3][3][3],
                            double **gkovr, double **gkonr, double **gkovc, double **gkonc,
                            double **gmkovr,double **gmkonr, double **gmkovc, double **gmkonc);
/*----------------------------------------------------------------------*
 |  s8_mat_ogden3.c                                      m.gee 06/03    |
 *----------------------------------------------------------------------*/
void s8_mat_ogden_uncoupled2(COMPOGDEN *mat, double *stress_cart, double C_cart[3][3][3][3],
                            double **gkonr, double **gmkovc);
/*----------------------------------------------------------------------*
 |  s8_mattransform.c                                    m.gee 06/03    |
 *----------------------------------------------------------------------*/
void s8_kov_cuca(double *t, const double **gkon);
void s8_kon_cacu(double *t, double **gkon);
void s8_kov_cacu(double *t, const double **gkov);
void s8_4kon_cacu(double Ccart[][3][3][3], double **gkon);
void s8_c4_to_C2(double C[][3][3][3],double **CC);
void s8_mat_linel_cart(STVENANT *mat,double C[][3][3][3],double **CC,double *strain);
/*----------------------------------------------------------------------*
 |  s8_static_mass.c                                     m.gee 06/03    |
 *----------------------------------------------------------------------*/
void s8static_mass(ELEMENT   *ele,                        /* the element structure */
                   S8_DATA   *data,                       /* element integration data */
                   MATERIAL  *mat,                        /* the material structure */
                   ARRAY     *estif_global,               /* element stiffness matrix (NOT initialized!) */
                   ARRAY     *emass_global,               /* element mass matrix      (NOT initialized!) */
                   double    *force,                      /* for internal forces (initialized!) */
                   int        kstep,                      /* actual step in nonlinear analysis */
                   int        init);                      /* init=1 -> init phase / init=0 -> calc. phase / init=-1 -> uninit phase */
/*----------------------------------------------------------------------*
 |  s8_mat_linel.c                                     m.gee 06/03    |
 *----------------------------------------------------------------------*/
void s8_mat_neohooke(NEO_HOOKE *mat,
                     double    *stress,
                     double   **CC,
                     double   **gmkonr,
                     double   **gmkonc,
                     double     detr,
                     double     detc);
/*----------------------------------------------------------------------*
 |  s8_fserv.f                                           m.gee 06/03    |
 *----------------------------------------------------------------------*/
void s8jacb(double *A, double *V);



/*----------------------------------------------------------------------*
 | s8_mat_viscoushyper.c                                 m.gee 06/03    |
 *----------------------------------------------------------------------*/
void s8_mat_ogden_viscous(
    ELEMENT            *ele,
    VISCOHYPER         *mat,
    double             *stress_cart,
    double              C_cart[3][3][3][3],
    double            **gkonr,
    double            **gmkovc,
    int                 gp
    );

#else

#endif


/* The DRT shell8 element uses some functions of the old shell8 element. A
 * hack! */

/*----------------------------------------------------------------------*
 |  s8_jaco.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8jaco(double *funct, double **deriv, double **x, double **xjm, double *hte, double **a3ref,
    double e3, int iel, double *det, double *deta, int init);
void s8_getdensity(MATERIAL *mat, double *density);
/*----------------------------------------------------------------------*
 |  s8_mat_linel.c                                       m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8_mat_linel(STVENANT *mat, double **g, double **CC);
void s8_mat_stress1(double *stress, double *strain, double **C);
/*----------------------------------------------------------------------*
 |  s8_mtr.c                                             m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8mtr(double **x, double **a3, double e3, double **gkov, double **gkon, double **gmkov,
    double **gmkon, double *det, double *funct, double **deriv, double *hte, int iel,
    double condfac, char string[]);
/*----------------------------------------------------------------------*
 |  s8_mat_linel.c                                     m.gee 06/03    |
 *----------------------------------------------------------------------*/
void s8_mat_neohooke(NEO_HOOKE *mat, double *stress, double **CC, double **gmkonr, double **gmkonc,
    double detr, double detc);
/*----------------------------------------------------------------------*
 |  s8_mat_ogden.c                                       m.gee 06/03    |
 *----------------------------------------------------------------------*/
void s8_mat_ogden_coupled(COMPOGDEN *mat, double *stress_cart, double C_cart[3][3][3][3],
    double **gkonr, double **gmkovc);
void s8_ogden_Ccart(double C[3][3][3][3], double C_cart[3][3][3][3], double N[3][3]);
void s8_ogden_cartPK2(double PK2[3][3], double PK2main[3], double N[3][3]);
void s8_ogden_principal_CG(double CG[3][3], double lambda[3], double N[3][3]);
void s8_kov_CGcuca(double T[3][3], double **gkon);
/*----------------------------------------------------------------------*
 |  s8_mat_ogden2.c                                      m.gee 06/03    |
 *----------------------------------------------------------------------*/
void s8_mat_ogden_uncoupled(COMPOGDEN *mat, double *stress_cart, double *strain,
    double C_cart[3][3][3][3], double **gkovr, double **gkonr, double **gkovc, double **gkonc,
    double **gmkovr, double **gmkonr, double **gmkovc, double **gmkonc);
/*----------------------------------------------------------------------*
 |  s8_mat_ogden3.c                                      m.gee 06/03    |
 *----------------------------------------------------------------------*/
void s8_mat_ogden_uncoupled2(COMPOGDEN *mat, double *stress_cart, double C_cart[3][3][3][3],
    double **gkonr, double **gmkovc);

/*----------------------------------------------------------------------*
 |  s8_mattransform.c                                    m.gee 06/03    |
 *----------------------------------------------------------------------*/
void s8_kov_cuca(double *t, const double **gkon);
void s8_kon_cacu(double *t, double **gkon);
void s8_kov_cacu(double *t, const double **gkov);
void s8_4kon_cacu(double Ccart[][3][3][3], double **gkon);
void s8_c4_to_C2(double C[][3][3][3], double **CC);
void s8_mat_linel_cart(STVENANT *mat, double C[][3][3][3], double **CC, double *strain);

#endif
