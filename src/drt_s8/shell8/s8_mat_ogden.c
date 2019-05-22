/*----------------------------------------------------------------------------*/
/*!
\brief shell8

\level 1

\maintainer Christoph Meier

*/
/*---------------------------------------------------------------------------*/
#ifdef D_SHELL8
#include "../../headers/standardtypes.h"
#include "../../headers/compiler_definitions.h" /* access to fortran routines */
#include "shell8.h"
#include "../../drt_lib/drt_dserror.H"
/*----------------------------------------------------------------------*
 | compressible ogden-material                            m.gee 6/03    |
 | no split in volumetric and deviatoric strains                        |
 *----------------------------------------------------------------------*/
void s8_mat_ogden_coupled(
    COMPOGDEN *mat, double *stress_cart, double C_cart[3][3][3][3], double **gkonr, double **gmkovc)

{
  int i, j, k, l, p;
  double mu;
  double E;
  double beta, mbeta;
  double lame1;
  double *mup;
  double *alfap;
  double J;
  double Jpowmbeta;
#if 0
      double      psi;
#endif

  double CG[3][3];
  double N[3][3];
  double lam2[3];
  double lam[3];
  double lampowalfap[3][3];
  double PK2[3];
  double PK2cart[3][3];

  double C[3][3][3][3];
  double C0000;
  double C0011;
  double C0022;
  double C1111;
  double C1122;
  double C2222;
  double C0101 = 0.0;
  double C0202 = 0.0;
  double C1212 = 0.0;

  double scal;
  double Ncross[3];

  /*----------------------------------------------------------------------*/
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      for (k = 0; k < 3; k++)
        for (l = 0; l < 3; l++)
        {
          C[i][j][k][l] = 0.0;
        }
  /*------------------------------------------------- init some constants */
  if (!(mat->init))
  {
    /* make mu = 2.0 * shear modulus */
    mu = 0.0;
    for (i = 0; i < 3; i++) mu += mat->alfap[i] * mat->mup[i];
    /* make Young's modulus */
    E = mu * (1.0 + mat->nue);
    /* make shear modulus */
    mu /= 2.0;
    /* make bulk modulus */
    mat->kappa = E / (3.0 * (1 - 2.0 * mat->nue));
    /* make lame constant no. 1 */
    mat->lambda = mat->kappa - (2.0 / 3.0) * mu;
    /* set init flag */
    mat->init = 1;
  }
  /*-------------------------------------------------- get some constants */
  beta = mat->beta;
  mbeta = -beta;
  alfap = mat->alfap;
  mup = mat->mup;
  lame1 = mat->lambda;
  /*---------------------------- make right Cauchy-Green strain tensor CG */
  /*
  CG = Ft * F = gmkovc_ij
  */
  CG[0][0] = gmkovc[0][0];
  CG[0][1] = gmkovc[0][1];
  CG[0][2] = gmkovc[0][2];
  CG[1][0] = gmkovc[1][0];
  CG[1][1] = gmkovc[1][1];
  CG[1][2] = gmkovc[1][2];
  CG[2][0] = gmkovc[2][0];
  CG[2][1] = gmkovc[2][1];
  CG[2][2] = gmkovc[2][2];
  /*
  CG is in covariant components in contravariant material bases, transform to cartesian
  */
  s8_kov_CGcuca(CG, gkonr);
  /*----------------------------------------------------------------------*/
  /* make spectral decomposition and principal axes of right Cauchy Green */
  /*
  (CG - lambda_a*I)*PHI_a = 0
  */
  s8_ogden_principal_CG(CG, lam2, N);
  /*---------------------------------------------- make principal strains */
  dsassert(lam2[0] > 0.0 && lam2[1] > 0.0 && lam2[2] > 0.0, "Principal strains smaller zero");
  lam[0] = sqrt(lam2[0]);
  lam[1] = sqrt(lam2[1]);
  lam[2] = sqrt(lam2[2]);
#if 1
  for (i = 0; i < 3; i++) mat->l[i] = lam[i];
#endif
  /*------------------------------------------- make 3. invariant == detF */
  J = lam[0] * lam[1] * lam[2];
  dsassert(J > 0.0, "detF <= 0.0 in Ogden material");
  Jpowmbeta = pow(J, mbeta);
  /*----------------------------------------- make powers lam[i]^alfap[p] */
  for (i = 0; i < 3; i++)
    for (p = 0; p < 3; p++) lampowalfap[i][p] = pow(lam[i], alfap[p]);
/*---------------- test orthogonality and unit length of eigenvectors N */
#if 1
  /*N0 * N1 = 0*/
  scal = N[0][0] * N[1][0] + N[0][1] * N[1][1] + N[0][2] * N[1][2];
  dsassert(fabs(scal) < EPS10, "eigenvectors N0,N1 not orthogonal");
  /*N0 * N2 = 0*/
  scal = N[0][0] * N[2][0] + N[0][1] * N[2][1] + N[0][2] * N[2][2];
  dsassert(fabs(scal) < EPS10, "eigenvectors N0,N2 not orthogonal");
  /*N1 * N2 = 0*/
  scal = N[1][0] * N[2][0] + N[1][1] * N[2][1] + N[1][2] * N[2][2];
  dsassert(fabs(scal) < EPS10, "eigenvectors N1,N2 not orthogonal");
  /*--------------------------- test proper orientation of eigenvectors N */
  /*N2 = N0 x N1*/
  Ncross[0] = N[0][1] * N[1][2] - N[0][2] * N[1][1];
  Ncross[1] = N[0][2] * N[1][0] - N[0][0] * N[1][2];
  Ncross[2] = N[0][0] * N[1][1] - N[0][1] * N[1][0];
  /*N2 * Ncross = 1.0*/
  scal = Ncross[0] * N[2][0] + Ncross[1] * N[2][1] + Ncross[2] * N[2][2];
  dsassert(fabs((scal - 1.0)) < EPS10, "eigenvectors do not form proper othogonal system");
#endif
/*----------------------------------------------------------------------*/
/*--------------------------------------------------------- make energy */
#if 0
psi = 0.0;
for (p=0; p<3; p++)
{
   psi += (mup[p]/alfap[p])*(lampowalfap[0][p]+lampowalfap[1][p]+lampowalfap[2][p]-3.0);
   psi -= mup[p]*log(J);
}
psi += (lame1/(beta*beta))*(Jpowmbeta-1.0+beta*log(J));
printf("  coupled PSI  %20.10f\n\n",psi);fflush(stdout);
#endif
  /*--------------------------- calculate the 2.PK stresses in eigenbases */
  PK2[0] = PK2[1] = PK2[2] = (lame1 / beta) * (1.0 - Jpowmbeta);
  for (p = 0; p < 3; p++)
    for (i = 0; i < 3; i++) PK2[i] += mup[p] * (lampowalfap[i][p] - 1.0);
  for (i = 0; i < 3; i++) PK2[i] /= (lam2[i]);
/*----------------------------------------------------------------------*/
#if 0
printf("PK2        [0] %14.8f PK2        [1] %14.8f PK2        [2] %14.8f\n\n",PK2[0],PK2[1],PK2[2]);
#endif
  /*----------------------- calculate the PK2 stresses in cartesian bases */
  s8_ogden_cartPK2(PK2cart, PK2, N);
  /*------------------ sort cartesian stresses to the vector shell8-style */
  stress_cart[0] = PK2cart[0][0];
  stress_cart[1] = PK2cart[0][1];
  stress_cart[2] = PK2cart[0][2];
  stress_cart[3] = PK2cart[1][1];
  stress_cart[4] = PK2cart[1][2];
  stress_cart[5] = PK2cart[2][2];
  /*---------------------- make deviatoric material tangent in eigenspace */
  /*================== components C_aaaa */
  C0000 = C1111 = C2222 = lame1 * ((2 / beta + 1.0) * Jpowmbeta - 2 / beta);
  for (p = 0; p < 3; p++)
  {
    C0000 += mup[p] * (2.0 + (alfap[p] - 2.0) * lampowalfap[0][p]);
    C1111 += mup[p] * (2.0 + (alfap[p] - 2.0) * lampowalfap[1][p]);
    C2222 += mup[p] * (2.0 + (alfap[p] - 2.0) * lampowalfap[2][p]);
  }
  C0000 /= (lam2[0] * lam2[0]);
  C1111 /= (lam2[1] * lam2[1]);
  C2222 /= (lam2[2] * lam2[2]);
  /*================== components C_aabb */
  C0011 = C0022 = C1122 = lame1 * Jpowmbeta;
  C0011 /= (lam2[0] * lam2[1]);
  C0022 /= (lam2[0] * lam2[2]);
  C1122 /= (lam2[1] * lam2[2]);
  /*================== components C_abab */
  if (fabs(lam2[0] - lam2[1]) > EPS12)
    C0101 = (PK2[0] - PK2[1]) / (lam2[0] - lam2[1]);
  else
    C0101 = 0.5 * (C0000 - C0011);

  if (fabs(lam2[0] - lam2[2]) > EPS12)
    C0202 = (PK2[0] - PK2[2]) / (lam2[0] - lam2[2]);
  else
    C0202 = 0.5 * (C0000 - C0022);

  if (fabs(lam2[1] - lam2[2]) > EPS12)
    C1212 = (PK2[1] - PK2[2]) / (lam2[1] - lam2[2]);
  else
    C1212 = 0.5 * (C1111 - C1122);
  /*--------------------------------------------- put everything together */
  C[0][0][0][0] = C0000;
  C[1][1][1][1] = C1111;
  C[2][2][2][2] = C2222;

  C[1][1][0][0] = C[0][0][1][1] = C0011;
  C[2][2][0][0] = C[0][0][2][2] = C0022;
  C[2][2][1][1] = C[1][1][2][2] = C1122;

  C[1][0][1][0] = C[0][1][0][1] = C0101;
  C[2][0][2][0] = C[0][2][0][2] = C0202;
  C[2][1][2][1] = C[1][2][1][2] = C1212;
  /*--------------------------------- calculate C_cart in cartesian basis */
  s8_ogden_Ccart(C, C_cart, N);
  /*----------------------------------------------------------------------*/
  return;
} /* end of s8_mat_ogden_coupled */



/*----------------------------------------------------------------------*
 | transform C in cartesian bases                         PK2 m.gee 6/03|
 *----------------------------------------------------------------------*/
void s8_ogden_Ccart(double C[3][3][3][3], double C_cart[3][3][3][3], double N[3][3])
{
  int i, j, k, l;
  /*----------------------------------------------------------------------*/
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      for (k = 0; k < 3; k++)
        for (l = 0; l < 3; l++)
        {
          /* a = 0  b = 0 */
          C_cart[i][j][k][l] = C[0][0][0][0] * N[i][0] * N[j][0] * N[k][0] * N[l][0];
          /* a = 0  b = 1 */
          C_cart[i][j][k][l] += C[0][0][1][1] * N[i][0] * N[j][0] * N[k][1] * N[l][1];
          /* a = 0  b = 2 */
          C_cart[i][j][k][l] += C[0][0][2][2] * N[i][0] * N[j][0] * N[k][2] * N[l][2];

          /* a = 1  b = 0 */
          C_cart[i][j][k][l] += C[1][1][0][0] * N[i][1] * N[j][1] * N[k][0] * N[l][0];
          /* a = 1  b = 1 */
          C_cart[i][j][k][l] += C[1][1][1][1] * N[i][1] * N[j][1] * N[k][1] * N[l][1];
          /* a = 1  b = 2 */
          C_cart[i][j][k][l] += C[1][1][2][2] * N[i][1] * N[j][1] * N[k][2] * N[l][2];

          /* a = 2  b = 0 */
          C_cart[i][j][k][l] += C[2][2][0][0] * N[i][2] * N[j][2] * N[k][0] * N[l][0];
          /* a = 2  b = 1 */
          C_cart[i][j][k][l] += C[2][2][1][1] * N[i][2] * N[j][2] * N[k][1] * N[l][1];
          /* a = 2  b = 2 */
          C_cart[i][j][k][l] += C[2][2][2][2] * N[i][2] * N[j][2] * N[k][2] * N[l][2];

          /* a = 0  b = 1 */
          C_cart[i][j][k][l] += C[0][1][0][1] * (N[i][0] * N[j][1] * N[k][0] * N[l][1] +
                                                    N[i][0] * N[j][1] * N[k][1] * N[l][0]);
          /* a = 0  b = 2 */
          C_cart[i][j][k][l] += C[0][2][0][2] * (N[i][0] * N[j][2] * N[k][0] * N[l][2] +
                                                    N[i][0] * N[j][2] * N[k][2] * N[l][0]);

          /* a = 1  b = 0 */
          C_cart[i][j][k][l] += C[1][0][1][0] * (N[i][1] * N[j][0] * N[k][1] * N[l][0] +
                                                    N[i][1] * N[j][0] * N[k][0] * N[l][1]);
          /* a = 1  b = 2 */
          C_cart[i][j][k][l] += C[1][2][1][2] * (N[i][1] * N[j][2] * N[k][1] * N[l][2] +
                                                    N[i][1] * N[j][2] * N[k][2] * N[l][1]);

          /* a = 2  b = 0 */
          C_cart[i][j][k][l] += C[2][0][2][0] * (N[i][2] * N[j][0] * N[k][2] * N[l][0] +
                                                    N[i][2] * N[j][0] * N[k][0] * N[l][2]);
          /* a = 2  b = 1 */
          C_cart[i][j][k][l] += C[2][1][2][1] * (N[i][2] * N[j][1] * N[k][2] * N[l][1] +
                                                    N[i][2] * N[j][1] * N[k][1] * N[l][2]);
        }
  /*----------------------------------------------------------------------*/
  return;
} /* end of s8_ogden_Ccart */


/*----------------------------------------------------------------------*
 | transform principal streeses PK2 in cartesian stresses PK2 m.gee 6/03|
 | PK2 = PK2main_a * N[][a] dyad N[][a] (sum over a)                    |
 *----------------------------------------------------------------------*/
void s8_ogden_cartPK2(double PK2[3][3], double PK2main[3], double N[3][3])
{
  int i, j;
  double dyad0[3][3], dyad1[3][3], dyad2[3][3];
  /*----------------------------------------------------------------------*/
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
    {
      dyad0[i][j] = N[i][0] * N[j][0];
      dyad1[i][j] = N[i][1] * N[j][1];
      dyad2[i][j] = N[i][2] * N[j][2];
    }
  PK2[0][0] = PK2main[0] * dyad0[0][0];
  PK2[0][1] = PK2main[0] * dyad0[0][1];
  PK2[0][2] = PK2main[0] * dyad0[0][2];
  PK2[1][0] = PK2main[0] * dyad0[1][0];
  PK2[1][1] = PK2main[0] * dyad0[1][1];
  PK2[1][2] = PK2main[0] * dyad0[1][2];
  PK2[2][0] = PK2main[0] * dyad0[2][0];
  PK2[2][1] = PK2main[0] * dyad0[2][1];
  PK2[2][2] = PK2main[0] * dyad0[2][2];

  PK2[0][0] += PK2main[1] * dyad1[0][0];
  PK2[0][1] += PK2main[1] * dyad1[0][1];
  PK2[0][2] += PK2main[1] * dyad1[0][2];
  PK2[1][0] += PK2main[1] * dyad1[1][0];
  PK2[1][1] += PK2main[1] * dyad1[1][1];
  PK2[1][2] += PK2main[1] * dyad1[1][2];
  PK2[2][0] += PK2main[1] * dyad1[2][0];
  PK2[2][1] += PK2main[1] * dyad1[2][1];
  PK2[2][2] += PK2main[1] * dyad1[2][2];

  PK2[0][0] += PK2main[2] * dyad2[0][0];
  PK2[0][1] += PK2main[2] * dyad2[0][1];
  PK2[0][2] += PK2main[2] * dyad2[0][2];
  PK2[1][0] += PK2main[2] * dyad2[1][0];
  PK2[1][1] += PK2main[2] * dyad2[1][1];
  PK2[1][2] += PK2main[2] * dyad2[1][2];
  PK2[2][0] += PK2main[2] * dyad2[2][0];
  PK2[2][1] += PK2main[2] * dyad2[2][1];
  PK2[2][2] += PK2main[2] * dyad2[2][2];


  /*
  PK2[0][0] = PK2main[0] * N[0][0]*N[0][0];
  PK2[0][1] = PK2main[0] * N[0][0]*N[1][0];
  PK2[0][2] = PK2main[0] * N[0][0]*N[2][0];
  PK2[1][0] = PK2main[0] * N[1][0]*N[0][0];
  PK2[1][1] = PK2main[0] * N[1][0]*N[1][0];
  PK2[1][2] = PK2main[0] * N[1][0]*N[2][0];
  PK2[2][0] = PK2main[0] * N[2][0]*N[0][0];
  PK2[2][1] = PK2main[0] * N[2][0]*N[1][0];
  PK2[2][2] = PK2main[0] * N[2][0]*N[2][0];

  PK2[0][0] += PK2main[1] * N[0][1]*N[0][1];
  PK2[0][1] += PK2main[1] * N[0][1]*N[1][1];
  PK2[0][2] += PK2main[1] * N[0][1]*N[2][1];
  PK2[1][0] += PK2main[1] * N[1][1]*N[0][1];
  PK2[1][1] += PK2main[1] * N[1][1]*N[1][1];
  PK2[1][2] += PK2main[1] * N[1][1]*N[2][1];
  PK2[2][0] += PK2main[1] * N[2][1]*N[0][1];
  PK2[2][1] += PK2main[1] * N[2][1]*N[1][1];
  PK2[2][2] += PK2main[1] * N[2][1]*N[2][1];

  PK2[0][0] += PK2main[2] * N[0][2]*N[0][2];
  PK2[0][1] += PK2main[2] * N[0][2]*N[1][2];
  PK2[0][2] += PK2main[2] * N[0][2]*N[2][2];
  PK2[1][0] += PK2main[2] * N[1][2]*N[0][2];
  PK2[1][1] += PK2main[2] * N[1][2]*N[1][2];
  PK2[1][2] += PK2main[2] * N[1][2]*N[2][2];
  PK2[2][0] += PK2main[2] * N[2][2]*N[0][2];
  PK2[2][1] += PK2main[2] * N[2][2]*N[1][2];
  PK2[2][2] += PK2main[2] * N[2][2]*N[2][2];
  */
  /* make symmetry */
  /*
  PK2[1][0] = PK2[0][1];
  PK2[2][0] = PK2[0][2];
  PK2[2][1] = PK2[1][2];
  */

  /*----------------------------------------------------------------------*/
  return;
} /* end of s8_ogden_cartPK2 */



void s8jacb(double *, double *);
/*----------------------------------------------------------------------*
 | make eigenvalue decomposition of Cauchy-Green strains  m.gee 6/03    |
 *----------------------------------------------------------------------*/
void s8_ogden_principal_CG(double CG[3][3], double lambda[3], double N[3][3])
{
  int i;
  double fstrain[9];
  double fn[9];
  /*----------------------------------------------------------------------*/
  for (i = 0; i < 9; i++) fn[i] = 0.0;

  fstrain[0] = CG[0][0];
  fstrain[1] = CG[1][0];
  fstrain[2] = CG[2][0];

  fstrain[3] = CG[0][1];
  fstrain[4] = CG[1][1];
  fstrain[5] = CG[2][1];

  fstrain[6] = CG[0][2];
  fstrain[7] = CG[1][2];
  fstrain[8] = CG[2][2];

  s8jacb(fstrain, fn);

  lambda[0] = fstrain[0];
  lambda[1] = fstrain[4];
  lambda[2] = fstrain[8];

  N[0][0] = fn[0];
  N[1][0] = fn[1];
  N[2][0] = fn[2];

  N[0][1] = fn[3];
  N[1][1] = fn[4];
  N[2][1] = fn[5];

  N[0][2] = fn[6];
  N[1][2] = fn[7];
  N[2][2] = fn[8];
  /*----------------------------------------------------------------------*/
  return;
} /* end of s8_ogden_principal_CG */



/*----------------------------------------------------------------------*
 | st.venant-kirchhoff-material                           m.gee 6/03    |
 *----------------------------------------------------------------------*/
void s8_mat_lineltmp(double E, double nue, double **g, double **CC)
{
  int i, j, k, l;

  /*----- shear correction coefficient not yet introduced */
  /* double xsi=1.0; */

  double C[3][3][3][3]; /*--------------------------- constitutive tensor */
  double l1, l2;        /*----------------------------------------- lame constants */
  double emod;          /*--------------------------------------- mat constants */
  /*----------------------------------------------------------------------*/
  emod = E;
  l1 = (emod * nue) / ((1.0 + nue) * (1.0 - 2.0 * nue));
  l2 = emod / (2.0 * (1.0 + nue));
  /*---------this is not very fast, but corresponds nicely with theory... */
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      for (k = 0; k < 3; k++)
        for (l = 0; l < 3; l++)
          C[i][j][k][l] = l1 * g[i][j] * g[k][l] + l2 * (g[i][k] * g[j][l] + g[i][l] * g[k][j]);
  /*----------------------------------------------------------------------*/
  CC[0][0] = C[0][0][0][0];
  CC[0][1] = C[0][0][1][0];
  CC[0][2] = C[0][0][2][0];
  CC[0][3] = C[0][0][1][1];
  CC[0][4] = C[0][0][2][1];
  CC[0][5] = C[0][0][2][2];

  CC[1][0] = C[1][0][0][0];
  CC[1][1] = C[1][0][1][0];
  CC[1][2] = C[1][0][2][0];
  CC[1][3] = C[1][0][1][1];
  CC[1][4] = C[1][0][2][1];
  CC[1][5] = C[1][0][2][2];

  CC[2][0] = C[2][0][0][0];
  CC[2][1] = C[2][0][1][0];
  CC[2][2] = C[2][0][2][0] /*/xsi*/;
  CC[2][3] = C[2][0][1][1];
  CC[2][4] = C[2][0][2][1] /*/xsi*/;
  CC[2][5] = C[2][0][2][2];

  CC[3][0] = C[1][1][0][0];
  CC[3][1] = C[1][1][1][0];
  CC[3][2] = C[1][1][2][0];
  CC[3][3] = C[1][1][1][1];
  CC[3][4] = C[1][1][2][1];
  CC[3][5] = C[1][1][2][2];

  CC[4][0] = C[2][1][0][0];
  CC[4][1] = C[2][1][1][0];
  CC[4][2] = C[2][1][2][0] /*/xsi*/;
  CC[4][3] = C[2][1][1][1];
  CC[4][4] = C[2][1][2][1] /*/xsi*/;
  CC[4][5] = C[2][1][2][2];

  CC[5][0] = C[2][2][0][0];
  CC[5][1] = C[2][2][1][0];
  CC[5][2] = C[2][2][2][0];
  CC[5][3] = C[2][2][1][1];
  CC[5][4] = C[2][2][2][1];
  CC[5][5] = C[2][2][2][2];
  /*----------------------------------------------------------------------*/
  return;
} /* end of s8_mat_lineltmp */
/*----------------------------------------------------------------------*
 |                                                        m.gee 5/03    |
 *----------------------------------------------------------------------*/
void s8_mat_linel_carttmp(double emod, double nue, double C[][3][3][3])
{
  /*
int i,j,k,l;
*/

  double l1, l2, ll2;
  /*
  double e[3][3];
  */
  /*----------------------------------------------------------------------*/
  l1 = (emod * nue) / ((1.0 + nue) * (1.0 - 2.0 * nue));
  l2 = emod / (2.0 * (1.0 + nue));
  ll2 = 2.0 * l2;
  /*----------------------------------------------------------------------*/
  /*
  e[0][0] = 1.0;
  e[1][0] = 0.0;
  e[2][0] = 0.0;
  e[0][1] = 0.0;
  e[1][1] = 1.0;
  e[2][1] = 0.0;
  e[0][2] = 0.0;
  e[1][2] = 0.0;
  e[2][2] = 1.0;
  */
  /*----------------------------------------------------------------------*/
  /*
  for (i=0; i<3; i++)
  for (j=0; j<3; j++)
  for (k=0; k<3; k++)
  for (l=0; l<3; l++)
  C[i][j][k][l] = l1*e[i][j]*e[k][l] + l2*( e[i][k]*e[j][l]+e[i][l]*e[k][j] );
  */
  /*----------------------------------------------------------------------*/
  C[0][0][0][0] = l1 + ll2;
  C[1][0][0][0] = 0.0;
  C[2][0][0][0] = 0.0;
  C[0][0][0][1] = 0.0;
  C[1][0][0][1] = l2;
  C[2][0][0][1] = 0.0;
  C[0][0][0][2] = 0.0;
  C[1][0][0][2] = 0.0;
  C[2][0][0][2] = l2;
  C[0][0][1][0] = 0.0;
  C[1][0][1][0] = l2;
  C[2][0][1][0] = 0.0;
  C[0][0][1][1] = l1;
  C[1][0][1][1] = 0.0;
  C[2][0][1][1] = 0.0;
  C[0][0][1][2] = 0.0;
  C[1][0][1][2] = 0.0;
  C[2][0][1][2] = 0.0;
  C[0][0][2][0] = 0.0;
  C[1][0][2][0] = 0.0;
  C[2][0][2][0] = l2;
  C[0][0][2][1] = 0.0;
  C[1][0][2][1] = 0.0;
  C[2][0][2][1] = 0.0;
  C[0][0][2][2] = l1;
  C[1][0][2][2] = 0.0;
  C[2][0][2][2] = 0.0;
  C[0][1][0][0] = 0.0;
  C[1][1][0][0] = l1;
  C[2][1][0][0] = 0.0;
  C[0][1][0][1] = l2;
  C[1][1][0][1] = 0.0;
  C[2][1][0][1] = 0.0;
  C[0][1][0][2] = 0.0;
  C[1][1][0][2] = 0.0;
  C[2][1][0][2] = 0.0;
  C[0][1][1][0] = l2;
  C[1][1][1][0] = 0.0;
  C[2][1][1][0] = 0.0;
  C[0][1][1][1] = 0.0;
  C[1][1][1][1] = l1 + ll2;
  C[2][1][1][1] = 0.0;
  C[0][1][1][2] = 0.0;
  C[1][1][1][2] = 0.0;
  C[2][1][1][2] = l2;
  C[0][1][2][0] = 0.0;
  C[1][1][2][0] = 0.0;
  C[2][1][2][0] = 0.0;
  C[0][1][2][1] = 0.0;
  C[1][1][2][1] = 0.0;
  C[2][1][2][1] = l2;
  C[0][1][2][2] = 0.0;
  C[1][1][2][2] = l1;
  C[2][1][2][2] = 0.0;
  C[0][2][0][0] = 0.0;
  C[1][2][0][0] = 0.0;
  C[2][2][0][0] = l1;
  C[0][2][0][1] = 0.0;
  C[1][2][0][1] = 0.0;
  C[2][2][0][1] = 0.0;
  C[0][2][0][2] = l2;
  C[1][2][0][2] = 0.0;
  C[2][2][0][2] = 0.0;
  C[0][2][1][0] = 0.0;
  C[1][2][1][0] = 0.0;
  C[2][2][1][0] = 0.0;
  C[0][2][1][1] = 0.0;
  C[1][2][1][1] = 0.0;
  C[2][2][1][1] = l1;
  C[0][2][1][2] = 0.0;
  C[1][2][1][2] = l2;
  C[2][2][1][2] = 0.0;
  C[0][2][2][0] = l2;
  C[1][2][2][0] = 0.0;
  C[2][2][2][0] = 0.0;
  C[0][2][2][1] = 0.0;
  C[1][2][2][1] = l2;
  C[2][2][2][1] = 0.0;
  C[0][2][2][2] = 0.0;
  C[1][2][2][2] = 0.0;
  C[2][2][2][2] = l1 + ll2;
  /*----------------------------------------------------------------------*/
  return;
} /* end of s8_mat_linel_cart */


/*----------------------------------------------------------------------*
 |                                                        m.gee 4/03    |
 | transform covariant components of a 2. Tensor from                   |
 | curvilinear to cartesian                                             |
 | storage mode is t[e11 e12 e13 e22 e23 e33]                           |
 | Must be called with contravariant base vectors !                     |
 | Tensor must be symmetric!                                            |
 *----------------------------------------------------------------------*/
void s8_kov_CGcuca(double T[3][3], double **gkon)
{
  int i, j;
  double Tcart[3][3];
  /*
  double c[3][3];
  */
  /*----------------------------------------------------------------------*/
  /* theory:
  for (k=0; k<3; k++)
  for (l=0; l<3; l++)
  {
     Tcart[k][l] = 0.0;
     for (i=0; i<3; i++)
     for (j=0; j<3; j++)
        Tcart[k][l] += gkon[k][i]*gkon[l][j]*T[i][j];
  }
  */
  /*----------------------------------------------------------------------*/
  Tcart[0][0] = 0.0;
  Tcart[0][1] = 0.0;
  Tcart[0][2] = 0.0;
  Tcart[1][1] = 0.0;
  Tcart[1][2] = 0.0;
  Tcart[2][2] = 0.0;
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
    {
      Tcart[0][0] += gkon[0][i] * gkon[0][j] * T[i][j];
      Tcart[0][1] += gkon[0][i] * gkon[1][j] * T[i][j];
      Tcart[0][2] += gkon[0][i] * gkon[2][j] * T[i][j];
      Tcart[1][1] += gkon[1][i] * gkon[1][j] * T[i][j];
      Tcart[1][2] += gkon[1][i] * gkon[2][j] * T[i][j];
      Tcart[2][2] += gkon[2][i] * gkon[2][j] * T[i][j];
    }
  /*----------------------------------------------------------------------*/
  T[0][0] = Tcart[0][0];
  T[1][0] = T[0][1] = Tcart[0][1];
  T[2][0] = T[0][2] = Tcart[0][2];
  T[1][1] = Tcart[1][1];
  T[2][1] = T[1][2] = Tcart[1][2];
  T[2][2] = Tcart[2][2];
  /*----------------------------------------------------------------------*/
  return;
} /* end of s8_kov_CGcuca */

#endif
