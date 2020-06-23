/*----------------------------------------------------------------------------*/
/*! \file
\brief shell8

\level 1

\maintainer Christoph Meier

*/
/*---------------------------------------------------------------------------*/

#include "../../headers/compiler_definitions.h" /* access to fortran routines */
#include "shell8.h"
#include "../../headers/standardtypes.h"
#include "../../drt_lib/drt_dserror.H"

/* declaration of fortran routine*/
void fortranpow(double *, double *, double *);

/*----------------------------------------------------------------------*
 | compressible ogden-material                            m.gee 6/03    |
 | split in volumetric and deviatoric strains                           |
 *----------------------------------------------------------------------*/
void s8_mat_ogden_uncoupled(COMPOGDEN *mat, DOUBLE *stress_cart, DOUBLE *strain,
    DOUBLE C_cart[3][3][3][3], DOUBLE **gkovr, DOUBLE **gkonr, DOUBLE **gkovc, DOUBLE **gkonc,
    DOUBLE **gmkovr, DOUBLE **gmkonr, DOUBLE **gmkovc, DOUBLE **gmkonc)
{
  INT i, j, k, l, p, a, c;
  DOUBLE beta, minusbeta;
  static DOUBLE monethird;
  DOUBLE *alfap;
  DOUBLE *mup;
  DOUBLE kappa;
  DOUBLE mu, E; /* shear and Young's modulus */
  DOUBLE work, work2, sum1, sum2;
  DOUBLE J;                 /* detF = third invariant */
  DOUBLE CG[3][3];          /* right-Cauchy Green strains */
  DOUBLE CGlambda2[3];      /* eigenvalues of CG */
  DOUBLE lambda[3];         /* principal stretches */
  DOUBLE lambda_bar[3];     /* deviatoric principal stretches */
  DOUBLE N[3][3];           /* eigenvectors of CG - principal stretch directions */
  DOUBLE PK2[3][3];         /* 2.PK stress tensor in cartesian bases */
  DOUBLE PK2main_dev[3];    /* deviatoric 2.PK stresses in principal directions (main stresses) */
  DOUBLE PK2main_vol[3];    /* volumetric 2.PK stresses in principal directions (main stresses) */
  DOUBLE PK2main[3];        /* 2.PK stresses in principal directions (main stresses) */
  DOUBLE C[3][3][3][3];     /* components of material tangent in principal directions */
  DOUBLE C_vol[3][3][3][3]; /* volumetric components of material tangent in principal directions */
  DOUBLE C_cartvol[3][3][3][3];
  /*
  DOUBLE              psi,psi1,psi2;
  */

  DOUBLE Ntest[3];
#if 0
/* for testing against holzapfel fortran routine*/
DOUBLE              df[28],bpr[3],tautil[3],PK2tautil[3],atilp[6][6];
DOUBLE              Ef;
#endif
  /*----------------------------------------------------------------------*/
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      for (k = 0; k < 3; k++)
        for (l = 0; l < 3; l++)
        {
          C[i][j][k][l] = 0.0;
          C_vol[i][j][k][l] = 0.0;

          C_cart[i][j][k][l] = 0.0;
          C_cartvol[i][j][k][l] = 0.0;
        }
  /*------------------------------------------------- init some constants */
  /*if (!(mat->init))
  {*/
  monethird = -1.0 / 3.0;
  /* make mu = 2.0 * shear modulus */
  mu = 0.0;
  for (i = 0; i < 3; i++) mu += mat->alfap[i] * mat->mup[i];
  /* make Young's modulus */
  E = mu * (1.0 + mat->nue);
  /* make shear modulus */
  mu /= 2.0;
  /* make bulk modulus */
  mat->kappa = E / (3.0 * (1 - 2.0 * (mat->nue)));
  /* make lame constant no. 1 */
  mat->lambda = mat->kappa - (2.0 / 3.0) * mu;
  /* set init flag */
  mat->init = 1;
  /*}*/
  /*-------------------------------------------------- get some constants */
  beta = mat->beta;
  minusbeta = -beta;
  alfap = mat->alfap;
  mup = mat->mup;
  kappa = mat->kappa;
  /*------------------------------------------- make deformation gradient */
  /*
  F = gkovc diad gkonr (defined in mixed components and in mixed shell base vectors)
  F = Fi_^j gkovc_i dyad gkonr_^j

  As we are not interested in F itself but in the right Cauchy Green strain tensor we need
  CG = F^T * F = cg_ij gkonr_^i dyad gkonr_^j (covariant components in kontravariant bases)

  To build F^T * F is too much work, but from the definition of the Green-Lagrange strains we know
  that E = 0.5*(F^T*F - I) = 0.5*(gmkovc_ij - gmkovr_ij) gkonr_^i dyad gkonr_^j.

  The unit tensor I = gmkovr_ij gkonr_^i dyad gkonr_^j, so the part F^T*F must be
  CG = F^T*F = gmkovc_ij gkonr_^i dyad gkonr_^j.

  or the components cg_ij = gmkovc_ij. Nice, isn't it?

  */
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
  /* make spectral decomposition and principal axes of right Cauchy Green */
  /*
  (CG - lambda_a*I)*PHI_a = 0
  */
  s8_ogden_principal_CG(CG, CGlambda2, N);
  /*---------------------------------------------- make principal strains */
  dsassert(CGlambda2[0] > 0.0 && CGlambda2[1] > 0.0 && CGlambda2[2] > 0.0,
      "Principal stretches smaller zero");
  lambda[0] = sqrt(CGlambda2[0]);
  lambda[1] = sqrt(CGlambda2[1]);
  lambda[2] = sqrt(CGlambda2[2]);
  /* test for equality of two or three of them */
  if (fabs(lambda[0] - lambda[1]) < EPS9) printf("lambda[0]==lambda[1]\n");
  if (fabs(lambda[0] - lambda[2]) < EPS9) printf("lambda[0]==lambda[2]\n");
  if (fabs(lambda[1] - lambda[2]) < EPS9) printf("lambda[1]==lambda[2]\n");
  /* test orthogonality of eigenvectors */
  /* N0 * N1 = 0.0 */
  work = N[0][0] * N[1][0] + N[0][1] * N[1][1] + N[0][2] * N[1][2];
  if (fabs(work) > EPS10) printf("eigenvectors N0 * N1 = %f\n", work);
  /* N0 * N2 = 0.0 */
  work = N[0][0] * N[2][0] + N[0][1] * N[2][1] + N[0][2] * N[2][2];
  if (fabs(work) > EPS10) printf("eigenvectors N0 * N2 = %f\n", work);
  /* N1 * N2 = 0.0 */
  work = N[1][0] * N[2][0] + N[1][1] * N[2][1] + N[1][2] * N[2][2];
  if (fabs(work) > EPS10) printf("eigenvectors N1 * N2 = %f\n", work);
  /* test length==1 of eigenvectors */
  work = N[0][0] * N[0][0] + N[0][1] * N[0][1] + N[0][2] * N[0][2];
  work = sqrt(work);
  if (fabs(work - 1.0) > EPS10) printf("eigenvectors N0 not unit length %f\n", work);
  work = N[1][0] * N[1][0] + N[1][1] * N[1][1] + N[1][2] * N[1][2];
  work = sqrt(work);
  if (fabs(work - 1.0) > EPS10) printf("eigenvectors N1 not unit length %f\n", work);
  work = N[2][0] * N[2][0] + N[2][1] * N[2][1] + N[2][2] * N[2][2];
  work = sqrt(work);
  if (fabs(work - 1.0) > EPS10) printf("eigenvectors N2 not unit length %f\n", work);
  /* test proper right hand system, build N3 by cross product */
  Ntest[0] = N[0][1] * N[1][2] - N[0][2] * N[1][1];
  Ntest[1] = N[0][2] * N[1][0] - N[0][0] * N[1][2];
  Ntest[2] = N[0][0] * N[1][1] - N[0][1] * N[1][0];
  work = N[2][0] * Ntest[0] + N[2][1] * Ntest[1] + N[2][2] * Ntest[2];
  if (fabs(work - 1.0) > EPS10)
    printf("eigenvectors are not right hand system  - cross product * N3 = %f\n", work);
/*==================call holzapfel routine */
#if 0
bpr[0] = CGlambda2[0];
bpr[1] = CGlambda2[1];
bpr[2] = CGlambda2[2];
df[21] = mup[0];
df[22] = alfap[0];
df[23] = mup[1];
df[24] = alfap[1];
df[25] = mup[2];
df[26] = alfap[2];
df[27] = 3.0;
s8wder3f(df,bpr,tautil,atilp,&Ef);
PK2tautil[0] = tautil[0]/CGlambda2[0];
PK2tautil[1] = tautil[1]/CGlambda2[1];
PK2tautil[2] = tautil[2]/CGlambda2[2];
#endif
  /*=========================================*/
  /*------------------------------------------- make 3. invariant == detF */
  J = lambda[0] * lambda[1] * lambda[2];
  dsassert(J > 0.0, "detF <= 0.0 in Ogden material");
  /*------ make modified principal strains lambda_bar = J^-0.333 * lambda */
  fortranpow(&J, &work, &monethird);
  lambda_bar[0] = work * lambda[0];
  lambda_bar[1] = work * lambda[1];
  lambda_bar[2] = work * lambda[2];
  /*--------------------------------------------------------- make energy */
  /*
  psi1 = 0.0;
  psi2 = 0.0;
  for (p=0; p<3; p++)
  {
     for (a=0; a<3; a++)
     {
         fortranpow(&(lambda_bar[a]),&work,&(alfap[p]));
         psi1 += (mup[p]/alfap[p])*(work-1.0);
     }
  }
  fortranpow(&J,&work,&minusbeta);
  psi2 = (kappa/(beta*beta))*(work-1.0+beta*log(J));
  psi = psi1+psi2;
  printf("uncoupled PSI1 %20.10f PSI2 %20.10f PSI %20.10f\n",psi1,psi2,psi);fflush(stdout);
  */
  /*-------------------------------- do deviatoric principal PK2 stresses */
  for (a = 0; a < 3; a++)
  {
    PK2main_dev[a] = 0.0;
    for (p = 0; p < 3; p++)
    {
      work2 = 0.0;
      for (c = 0; c < 3; c++)
      {
        fortranpow(&(lambda_bar[c]), &work, &(alfap[p]));
        work2 += work;
      }
      fortranpow(&(lambda_bar[a]), &work, &(alfap[p]));
      PK2main_dev[a] += mup[p] * (work + monethird * work2);
    }
    PK2main_dev[a] /= CGlambda2[a];
    PK2main[a] = PK2main_dev[a];
  }
  /*------------------------------ make volumetric principal PK2 stresses */
  fortranpow(&J, &work, &minusbeta);
  for (a = 0; a < 3; a++)
  {
    PK2main_vol[a] = (kappa / (beta * CGlambda2[a])) * (1.0 - work);
    PK2main[a] += PK2main_vol[a];
  }
  /*=========================================*/
  /*
  printf("     tautil[0] %14.8f      tautil[1] %14.8f      tautil[2]
  %14.8f\n",tautil[0],tautil[1],tautil[2]); printf("  PK2tautil[0] %14.8f   PK2tautil[1] %14.8f
  PK2tautil[2] %14.8f\n",PK2tautil[0],PK2tautil[1],PK2tautil[2]); printf("PK2main_dev[0] %14.8f
  PK2main_dev[1] %14.8f PK2main_dev[2] %14.8f\n",PK2main_dev[0],PK2main_dev[1],PK2main_dev[2]);
  printf("PK2main_vol[0] %14.8f PK2main_vol[1] %14.8f PK2main_vol[2]
  %14.8f\n",PK2main_vol[0],PK2main_vol[1],PK2main_vol[2]);
  */
  /*=========================================*/
  /*----------------------- calculate the PK2 stresses in cartesian bases */
  /*
  PK2 = PK2main_a * N_a dyad N_a   (sum over a )
  */
  s8_ogden_cartPK2(PK2, PK2main, N);
  /*------------------ sort cartesian stresses to the vector shell8-style */
  stress_cart[0] = PK2[0][0];
  stress_cart[1] = PK2[0][1];
  stress_cart[2] = PK2[0][2];
  stress_cart[3] = PK2[1][1];
  stress_cart[4] = PK2[1][2];
  stress_cart[5] = PK2[2][2];
  /*================ make components of C[][][][] in principal directions */
  /*-------------------------------------------------- make C[a][a][a][a] */
  /*--------------------------------------------- deviatoric contribution */
  for (a = 0; a < 3; a++)
  {
    C[a][a][a][a] = 0.0;
    for (p = 0; p < 3; p++)
    {
      work2 = 0.0;
      for (c = 0; c < 3; c++)
      {
        fortranpow(&(lambda_bar[c]), &work, &(alfap[p]));
        work2 += work;
      }
      fortranpow(&(lambda_bar[a]), &work, &(alfap[p]));
      C[a][a][a][a] += mup[p] * alfap[p] * (work / 3.0 + work2 / 9.0);
    }
    C[a][a][a][a] /= (CGlambda2[a] * CGlambda2[a]);
  }
  /*-------------------------------------- make C[a][a][b][b] with a != b */
  /*--------------------------------------------- deviatoric contribution */
  C[0][0][1][1] = 0.0;
  C[0][0][2][2] = 0.0;
  C[1][1][2][2] = 0.0;
  for (p = 0; p < 3; p++)
  {
    work2 = 0.0;
    for (c = 0; c < 3; c++)
    {
      fortranpow(&(lambda_bar[c]), &work, &(alfap[p]));
      work2 += work;
    }
    work2 /= 9.0;

    /* a = 0 b = 1 */
    fortranpow(&(lambda_bar[0]), &sum1, &(alfap[p]));
    fortranpow(&(lambda_bar[1]), &sum2, &(alfap[p]));
    C[0][0][1][1] += mup[p] * alfap[p] * (-sum1 / 3.0 - sum2 / 3.0 + work2);
    /* a = 0 b = 2 */
    fortranpow(&(lambda_bar[0]), &sum1, &(alfap[p]));
    fortranpow(&(lambda_bar[2]), &sum2, &(alfap[p]));
    C[0][0][2][2] += mup[p] * alfap[p] * (-sum1 / 3.0 - sum2 / 3.0 + work2);
    /* a = 1 b = 2 */
    fortranpow(&(lambda_bar[1]), &sum1, &(alfap[p]));
    fortranpow(&(lambda_bar[2]), &sum2, &(alfap[p]));
    C[1][1][2][2] += mup[p] * alfap[p] * (-sum1 / 3.0 - sum2 / 3.0 + work2);
  }
  C[0][0][1][1] /= (CGlambda2[0] * CGlambda2[1]);
  C[0][0][2][2] /= (CGlambda2[0] * CGlambda2[2]);
  C[1][1][2][2] /= (CGlambda2[1] * CGlambda2[2]);
  C[1][1][0][0] = C[0][0][1][1];
  C[2][2][0][0] = C[0][0][2][2];
  C[2][2][1][1] = C[1][1][2][2];
  /*--------------------------------------------- make C[a][b][a][b] a!=b */
  /*--------------------------------------------- deviatoric contribution */
  if (fabs(CGlambda2[0] - CGlambda2[1]) > EPS12)
    C[0][1][0][1] = (PK2main_dev[0] - PK2main_dev[1]) / (CGlambda2[0] - CGlambda2[1]);
  else
    C[0][1][0][1] = 0.5 * (C[0][0][0][0] - C[0][0][1][1]);


  if (fabs(CGlambda2[0] - CGlambda2[2]) > EPS12)
    C[0][2][0][2] = (PK2main_dev[0] - PK2main_dev[2]) / (CGlambda2[0] - CGlambda2[2]);
  else
    C[0][2][0][2] = 0.5 * (C[0][0][0][0] - C[0][0][2][2]);


  if (fabs(CGlambda2[1] - CGlambda2[2]) > EPS12)
    C[1][2][1][2] = (PK2main_dev[1] - PK2main_dev[2]) / (CGlambda2[1] - CGlambda2[2]);
  else
    C[1][2][1][2] = 0.5 * (C[1][1][1][1] - C[1][1][2][2]);

  C[1][0][1][0] = C[0][1][0][1];
  C[2][0][2][0] = C[0][2][0][2];
  C[2][1][2][1] = C[1][2][1][2];
  /*---------------------------------------------- make C_vol[a][a][a][a] */
  /*--------------------------------------------- volumetric contribution */
  fortranpow(&J, &work, &minusbeta);
  C_vol[0][0][0][0] = kappa * (((2.0 / beta) + 1.0) * work - (2.0 / beta));
  C_vol[1][1][1][1] = C_vol[2][2][2][2] = C_vol[0][0][0][0];

  C_vol[0][0][0][0] /= (CGlambda2[0] * CGlambda2[0]);
  C_vol[1][1][1][1] /= (CGlambda2[1] * CGlambda2[1]);
  C_vol[2][2][2][2] /= (CGlambda2[2] * CGlambda2[2]);
  /*-------------------------------------- make C[a][a][b][b] with a != b */
  /*fortranpow(&J,&work,&minusbeta);*/
  C_vol[0][0][1][1] = (kappa / (CGlambda2[0] * CGlambda2[1])) * work;
  C_vol[0][0][2][2] = (kappa / (CGlambda2[0] * CGlambda2[2])) * work;
  C_vol[1][1][2][2] = (kappa / (CGlambda2[1] * CGlambda2[2])) * work;
  C_vol[1][1][0][0] = C_vol[0][0][1][1];
  C_vol[2][2][0][0] = C_vol[0][0][2][2];
  C_vol[2][2][1][1] = C_vol[1][1][2][2];
  /*--------------------------------------------- make C[a][b][a][b] a!=b */
  /*--------------------------------------------- volumetric contribution */
  if (fabs(CGlambda2[0] - CGlambda2[1]) > EPS12)
    C_vol[0][1][0][1] = (PK2main_vol[0] - PK2main_vol[1]) / (CGlambda2[0] - CGlambda2[1]);
  else
    C_vol[0][1][0][1] = 0.5 * (C_vol[0][0][0][0] - C_vol[0][0][1][1]);


  if (fabs(CGlambda2[0] - CGlambda2[2]) > EPS12)
    C_vol[0][2][0][2] = (PK2main_vol[0] - PK2main_vol[2]) / (CGlambda2[0] - CGlambda2[2]);
  else
    C_vol[0][2][0][2] = 0.5 * (C_vol[0][0][0][0] - C_vol[0][0][2][2]);


  if (fabs(CGlambda2[1] - CGlambda2[2]) > EPS12)
    C_vol[1][2][1][2] = (PK2main_vol[1] - PK2main_vol[2]) / (CGlambda2[1] - CGlambda2[2]);
  else
    C_vol[1][2][1][2] = 0.5 * (C_vol[1][1][1][1] - C_vol[1][1][2][2]);

  C_vol[1][0][1][0] = C_vol[0][1][0][1];
  C_vol[2][0][2][0] = C_vol[0][2][0][2];
  C_vol[2][1][2][1] = C_vol[1][2][1][2];
  /*--------------------------------- calculate C_cart in cartesian basis */
  s8_ogden_Ccart(C, C_cart, N);
  s8_ogden_Ccart(C_vol, C_cartvol, N);
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      for (k = 0; k < 3; k++)
        for (l = 0; l < 3; l++) C_cart[i][j][k][l] += C_cartvol[i][j][k][l];
  /*----------------------------------------------------------------------*/
  return;
} /* end of s8_mat_ogden_uncoupled */
