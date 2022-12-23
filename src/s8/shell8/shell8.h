/*----------------------------------------------------------------------------*/
/*! \file
\brief shell8

\level 1


*/
/*---------------------------------------------------------------------------*/
#include "am.h"
#include "materials.h"

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
