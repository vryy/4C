#ifdef D_SHELL8
#include "../headers/standardtypes.h"
#include "shell8.h"
/*----------------------------------------------------------------------*
 | compressible ogden-material                            m.gee 6/03    |
 | no split in volumetric and deviatoric strains                        |
 *----------------------------------------------------------------------*/
void s8_mat_ogden_coupled(COMPOGDEN *mat, DOUBLE *stress_cart, DOUBLE *strain, DOUBLE C_cart[][3][3][3],
                          DOUBLE **gkovr, DOUBLE **gkonr, DOUBLE **gkovc, DOUBLE **gkonc,
                          DOUBLE **gmkovr,DOUBLE **gmkonr, DOUBLE **gmkovc, DOUBLE **gmkonc)
{
INT                 i,j,k,l,p,a,b,c,d;
DOUBLE              nue;
DOUBLE              beta,minusbeta;
DOUBLE             *alfap;
DOUBLE             *mup;
DOUBLE              lame1;
DOUBLE              kappa;
DOUBLE              mu,E;                  /* shear and Young's modulus */
DOUBLE              work;
DOUBLE              F[3][3];               /* deformation gradient */
DOUBLE              J;                     /* detF = third invariant */
DOUBLE              CG[3][3];              /* right-Cauchy Green strains */
/*DOUBLE              GL[3][3];*/              /* Green-Lagrange strains for testing */
DOUBLE              CGlambda2[3];          /* eigenvalues of CG */
DOUBLE              lambda[3];             /* principal stretches */
DOUBLE              N[3][3];               /* eigenvectors of CG - principal stretch directions */
DOUBLE              PK2[3][3];             /* 2.PK stress tensor in cartesian bases */
DOUBLE              PK2main[3]/*,PK2venant[3]*/;/* 2.PK stresses in principal directions (main stresses) */
DOUBLE              C[3][3][3][3];         /* components of material tangent in principal directions */
DOUBLE              psi,psi1,psi2;
#ifdef DEBUG 
dstrc_enter("s8_mat_ogden_coupled");
#endif
/*----------------------------------------------------------------------*/
for (i=0; i<3; i++)
for (j=0; j<3; j++)
for (k=0; k<3; k++)
for (l=0; l<3; l++)
{
   C[i][j][k][l]     = 0.0;
}
/*------------------------------------------------- init some constants */
/*if (!(mat->init))
{*/
   /* make mu = 2.0 * shear modulus */
   mu = 0.0;
   for (i=0; i<3; i++) mu += mat->alfap[i] * mat->mup[i];
   /* make Young's modulus */
   E  = mu*(1.0+mat->nue);
   /* make shear modulus */
   mu /= 2.0;
   /* make bulk modulus */
   mat->kappa = E / (3.0*(1-2.0*mat->nue));
   /* make lame constant no. 1 */
   mat->lambda = mat->kappa - (2.0/3.0)*mu;
   /* set init flag */
   mat->init=1;
/*}*/
/*-------------------------------------------------- get some constants */
nue       = mat->nue;
beta      = mat->beta;
minusbeta = -beta;
alfap     = mat->alfap;
mup       = mat->mup;
lame1     = mat->lambda;
kappa     = mat->kappa;
/*------------------------------------------- make deformation gradient */
/*
F = gkovc diad gkonr (defined in mixed components and in mixed shell base vectors)
F = Fi_^j gkovc_i dyad gkonr_^j

As we are not interested in F itself but in the right Cauchy Green strain tensor we need
CG = F^T * F = cg_ij gkonr_^i dyad gkonr_^j (covariant components in kontravariant bases)

To build F^T * F is too much work, but from the definition of the Green-Lagrange strains we know that
E = 0.5*(F^T*F - I) = 0.5*(gmkovc_ij - gmkovr_ij) gkonr_^i dyad gkonr_^j.

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
s8_kov_CGcuca(CG,gkonr); 
/*----------------------------- make Green-Lagrange strains for testing */
/*
GL[0][0] = 0.5*(CG[0][0]-1.0);
GL[0][1] = 0.5*(CG[0][1]);
GL[0][2] = 0.5*(CG[0][2]);

GL[1][0] = 0.5*(CG[1][0]);
GL[1][1] = 0.5*(CG[1][1]-1.0);
GL[1][2] = 0.5*(CG[1][2]);

GL[2][0] = 0.5*(CG[2][0]);
GL[2][1] = 0.5*(CG[2][1]);
GL[2][2] = 0.5*(CG[2][2]-1.0);

strain[0] = GL[0][0];   
strain[1] = GL[0][1];   
strain[2] = GL[0][2];   
strain[3] = GL[1][1];   
strain[4] = GL[1][2];   
strain[5] = GL[2][2]; 
/*----------------------------------------------------------------------*/
/* make spectral decomposition and principal axes of right Cauchy Green */
/*
(CG - lambda_a*I)*PHI_a = 0
*/
s8_ogden_principal_CG(CG,CGlambda2,N);
/*---------------------------------------------- make principal strains */
dsassert(CGlambda2[0]>0.0 && CGlambda2[1]>0.0 && CGlambda2[2]>0.0,"Principal stretches smaller zero");
lambda[0] = sqrt(CGlambda2[0]);
lambda[1] = sqrt(CGlambda2[1]);
lambda[2] = sqrt(CGlambda2[2]);
/*------------------------------------------- make 3. invariant == detF */
J         = lambda[0]*lambda[1]*lambda[2];
dsassert(J>0.0,"detF <= 0.0 in Ogden material");
/*--------------------------------------------------------- make energy */
/*
psi1 = 0.0;
psi2 = 0.0;
for (p=0; p<3; p++)
{
   for (a=0; a<3; a++)
   {
      fortranpow(&(lambda[a]),&work,&(alfap[p]));
      psi1 += (mup[p]/alfap[p])*work;
   }
   psi1 -= (mup[p]/alfap[p])*3.0;
   psi1 -= mup[p]*log(J);
}
fortranpow(&J,&work,&minusbeta);
psi2 = (lame1/(beta*beta))*(work-1.0+beta*log(J));
psi = psi1+psi2;
printf("  coupled PSI1 %20.10f PSI2 %20.10f PSI %20.10f\n",psi1,psi2,psi);fflush(stdout);
*/
/*------------------------- calculate the 2.PK stresses (contravariant) */
for (a=0; a<3; a++)
{
   PK2main[a] = 0.0;
   for (p=0; p<3; p++) 
   {
      fortranpow(&(lambda[a]),&work,&(alfap[p]));
      PK2main[a] += (mup[p])*(work-1.0);
   }
   fortranpow(&J,&work,&minusbeta);

   PK2main[a] += (lame1/beta)*(1-work);

   PK2main[a] /= CGlambda2[a];
}
/*----------------------- calculate the PK2 stresses in cartesian bases */
/*
PK2 = PK2main_a * N_a dyad N_a   (sum over a )
*/
s8_ogden_cartPK2(PK2,PK2main,N);
/* sort cartesian stresses to the vector shell8-style */
stress_cart[0] = PK2[0][0];   
stress_cart[1] = PK2[0][1];   
stress_cart[2] = PK2[0][2];   
stress_cart[3] = PK2[1][1];   
stress_cart[4] = PK2[1][2];   
stress_cart[5] = PK2[2][2];   
/*------------------------------------------------------ make PK2venant */
/*
work = 0.0;
for (b=0; b<3; b++)
   work += (CGlambda2[b]-1);
for (a=0; a<3; a++)
   PK2venant[a] = (lame1/2.0)*work + mu*(lambda[a]*lambda[a]-1.0);
s8_ogden_cartPK2(PK2,PK2venant,N);
stress_cart[0] = PK2[0][0];   
stress_cart[1] = PK2[0][1];   
stress_cart[2] = PK2[0][2];   
stress_cart[3] = PK2[1][1];   
stress_cart[4] = PK2[1][2];   
stress_cart[5] = PK2[2][2];   
*/
/*================ make components of C[][][][] in principal directions */
/*-------------------------------------------------- make C[a][a][a][a] */
for (a=0; a<3; a++)
{
   C[a][a][a][a] = 0.0;
   for (p=0; p<3; p++) 
   {
      fortranpow(&(lambda[a]),&work,&(alfap[p]));
      C[a][a][a][a] += mup[p]*(2.0+(alfap[p]-2.0)*work);
   }

   fortranpow(&J,&work,&minusbeta);
   C[a][a][a][a] += lame1*(((2.0/beta)+1.0)*work-(2.0/beta));
   
   work = CGlambda2[a]*CGlambda2[a];
   C[a][a][a][a] /= work;
}
/*-------------------------------------- make C[a][a][b][b] with a != b */
fortranpow(&J,&work,&minusbeta);
C[0][0][1][1] = (lame1/(CGlambda2[0]*CGlambda2[1]))*work;
C[0][0][2][2] = (lame1/(CGlambda2[0]*CGlambda2[2]))*work;
C[1][1][2][2] = (lame1/(CGlambda2[1]*CGlambda2[2]))*work;
C[1][1][0][0] = C[0][0][1][1];
C[2][2][0][0] = C[0][0][2][2];
C[2][2][1][1] = C[1][1][2][2];
/*--------------------------------------------- make C[a][b][a][b] a!=b */
if (FABS(CGlambda2[0]-CGlambda2[1])>EPS12)
C[0][1][0][1] = (PK2main[0]-PK2main[1])/(CGlambda2[0]-CGlambda2[1]);
else
C[0][1][0][1] = 0.5*(C[0][0][0][0]-C[0][0][1][1]);


if (FABS(CGlambda2[0]-CGlambda2[2])>EPS12)
C[0][2][0][2] = (PK2main[0]-PK2main[2])/(CGlambda2[0]-CGlambda2[2]);
else
C[0][2][0][2] = 0.5*(C[0][0][0][0]-C[0][0][2][2]);


if (FABS(CGlambda2[1]-CGlambda2[2])>EPS12)
C[1][2][1][2] = (PK2main[1]-PK2main[2])/(CGlambda2[1]-CGlambda2[2]);
else
C[1][2][1][2] = 0.5*(C[1][1][1][1]-C[1][1][2][2]);

C[1][0][1][0] = C[0][1][0][1];
C[2][0][2][0] = C[0][2][0][2];
C[2][1][2][1] = C[1][2][1][2];
/*-------------------------------- calculate C_cart in cartesian basis */
s8_ogden_Ccart(C,C_cart,N);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s8_mat_ogden_coupled */



/*----------------------------------------------------------------------*
 | transform C in cartesian bases                         PK2 m.gee 6/03|
 *----------------------------------------------------------------------*/
void s8_ogden_Ccart(DOUBLE C[3][3][3][3], DOUBLE C_cart[3][3][3][3], DOUBLE N[3][3])
{
INT i,j,k,l;
#ifdef DEBUG 
dstrc_enter("s8_ogden_Ccart");
#endif
/*----------------------------------------------------------------------*/
for (i=0; i<3; i++)
for (j=0; j<3; j++)
for (k=0; k<3; k++)
for (l=0; l<3; l++)
{
   /* a = 0  b = 0 */
   C_cart[i][j][k][l]  = C[0][0][0][0] * N[i][0]*N[j][0]*N[k][0]*N[l][0];
   /* a = 0  b = 1 */
   C_cart[i][j][k][l] += C[0][0][1][1] * N[i][0]*N[j][0]*N[k][1]*N[l][1];
   /* a = 0  b = 2 */
   C_cart[i][j][k][l] += C[0][0][2][2] * N[i][0]*N[j][0]*N[k][2]*N[l][2];

   /* a = 1  b = 0 */
   C_cart[i][j][k][l] += C[1][1][0][0] * N[i][1]*N[j][1]*N[k][0]*N[l][0];
   /* a = 1  b = 1 */
   C_cart[i][j][k][l] += C[1][1][1][1] * N[i][1]*N[j][1]*N[k][1]*N[l][1];
   /* a = 1  b = 2 */
   C_cart[i][j][k][l] += C[1][1][2][2] * N[i][1]*N[j][1]*N[k][2]*N[l][2];

   /* a = 2  b = 0 */
   C_cart[i][j][k][l] += C[2][2][0][0] * N[i][2]*N[j][2]*N[k][0]*N[l][0];
   /* a = 2  b = 1 */
   C_cart[i][j][k][l] += C[2][2][1][1] * N[i][2]*N[j][2]*N[k][1]*N[l][1];
   /* a = 2  b = 2 */
   C_cart[i][j][k][l] += C[2][2][2][2] * N[i][2]*N[j][2]*N[k][2]*N[l][2];
   
   /* a = 0  b = 1 */
   C_cart[i][j][k][l] += C[0][1][0][1] * (N[i][0]*N[j][1]*N[k][0]*N[l][1] + N[i][0]*N[j][1]*N[k][1]*N[l][0]);
   /* a = 0  b = 2 */
   C_cart[i][j][k][l] += C[0][2][0][2] * (N[i][0]*N[j][2]*N[k][0]*N[l][2] + N[i][0]*N[j][2]*N[k][2]*N[l][0]);

   /* a = 1  b = 0 */
   C_cart[i][j][k][l] += C[1][0][1][0] * (N[i][1]*N[j][0]*N[k][1]*N[l][0] + N[i][1]*N[j][0]*N[k][0]*N[l][1]);
   /* a = 1  b = 2 */
   C_cart[i][j][k][l] += C[1][2][1][2] * (N[i][1]*N[j][2]*N[k][1]*N[l][2] + N[i][1]*N[j][2]*N[k][2]*N[l][1]);

   /* a = 2  b = 0 */
   C_cart[i][j][k][l] += C[2][0][2][0] * (N[i][2]*N[j][0]*N[k][2]*N[l][0] + N[i][2]*N[j][0]*N[k][0]*N[l][2]);
   /* a = 2  b = 1 */
   C_cart[i][j][k][l] += C[2][1][2][1] * (N[i][2]*N[j][1]*N[k][2]*N[l][1] + N[i][2]*N[j][1]*N[k][1]*N[l][2]);
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s8_ogden_Ccart */


/*----------------------------------------------------------------------*
 | transform principal streeses PK2 in cartesian stresses PK2 m.gee 6/03|
 | PK2 = PK2main_a * N[][a] dyad N[][a] (sum over a)                    |
 *----------------------------------------------------------------------*/
void s8_ogden_cartPK2(DOUBLE PK2[3][3], DOUBLE PK2main[3], DOUBLE N[3][3])
{
INT i,j;
DOUBLE dyad0[3][3],dyad1[3][3],dyad2[3][3];
#ifdef DEBUG 
dstrc_enter("s8_ogden_cartPK2");
#endif
/*----------------------------------------------------------------------*/
for (i=0; i<3; i++)
for (j=0; j<3; j++)
{
   dyad0[i][j] = N[i][0]*N[j][0];
   dyad1[i][j] = N[i][1]*N[j][1];
   dyad2[i][j] = N[i][2]*N[j][2];
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
/* make symmetry 
PK2[1][0] = PK2[0][1];
PK2[2][0] = PK2[0][2]; 
PK2[2][1] = PK2[1][2];*/
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s8_ogden_cartPK2 */




/*----------------------------------------------------------------------*
 | make eigenvalue decomposition of Cauchy-Green strains  m.gee 6/03    |
 *----------------------------------------------------------------------*/
void s8_ogden_principal_CG(DOUBLE CG[3][3], DOUBLE lambda[3], DOUBLE N[3][3])
{
INT          i;
DOUBLE       fstrain[9];
DOUBLE       fn[9];
#ifdef DEBUG 
dstrc_enter("s8_ogden_principal_CG");
#endif
/*----------------------------------------------------------------------*/
for (i=0; i<9; i++) 
   fn[i] = 0.0;

fstrain[0] = CG[0][0];
fstrain[1] = CG[1][0];
fstrain[2] = CG[2][0];

fstrain[3] = CG[0][1];
fstrain[4] = CG[1][1];
fstrain[5] = CG[2][1];

fstrain[6] = CG[0][2];
fstrain[7] = CG[1][2];
fstrain[8] = CG[2][2];

s8jacb(fstrain,fn);

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
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s8_ogden_principal_CG */













/*----------------------------------------------------------------------*
 | st.venant-kirchhoff-material                           m.gee 6/03    |
 *----------------------------------------------------------------------*/
INT s8_mat_lineltmp(DOUBLE E, DOUBLE nue, DOUBLE **g, DOUBLE **CC)
{
INT i,j,k,l;
DOUBLE xsi=1.0; /*----- shear correction coefficient not yet introduced */
DOUBLE C[3][3][3][3]; /*--------------------------- constitutive tensor */
DOUBLE l1,l2;/*----------------------------------------- lame constants */
DOUBLE emod;/*--------------------------------------- mat constants */
#ifdef DEBUG 
dstrc_enter("s8_mat_lineltmp");
#endif
/*----------------------------------------------------------------------*/
emod = E;
l1 = (emod*nue) / ((1.0+nue)*(1.0-2.0*nue));
l2 = emod/ (2.0*(1.0+nue));
/*---------this is not very fast, but corresponds nicely with theory... */
for (i=0; i<3; i++)
for (j=0; j<3; j++)
for (k=0; k<3; k++)
for (l=0; l<3; l++)
C[i][j][k][l] = l1*g[i][j]*g[k][l] + l2*( g[i][k]*g[j][l]+g[i][l]*g[k][j] );
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
CC[2][2] = C[2][0][2][0]/*/xsi*/;
CC[2][3] = C[2][0][1][1];
CC[2][4] = C[2][0][2][1]/*/xsi*/;
CC[2][5] = C[2][0][2][2];

CC[3][0] = C[1][1][0][0];
CC[3][1] = C[1][1][1][0];
CC[3][2] = C[1][1][2][0];
CC[3][3] = C[1][1][1][1];
CC[3][4] = C[1][1][2][1];
CC[3][5] = C[1][1][2][2];

CC[4][0] = C[2][1][0][0];
CC[4][1] = C[2][1][1][0];
CC[4][2] = C[2][1][2][0]/*/xsi*/;
CC[4][3] = C[2][1][1][1];
CC[4][4] = C[2][1][2][1]/*/xsi*/;
CC[4][5] = C[2][1][2][2];

CC[5][0] = C[2][2][0][0];
CC[5][1] = C[2][2][1][0];
CC[5][2] = C[2][2][2][0];
CC[5][3] = C[2][2][1][1];
CC[5][4] = C[2][2][2][1];
CC[5][5] = C[2][2][2][2];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s8_mat_lineltmp */
/*----------------------------------------------------------------------*
 |                                                        m.gee 5/03    |
 *----------------------------------------------------------------------*/
INT s8_mat_linel_carttmp(DOUBLE emod, DOUBLE nue, 
                         DOUBLE C[][3][3][3])
{
INT i,j,k,l;
DOUBLE l1,l2,ll2;
/*
DOUBLE e[3][3];
*/
#ifdef DEBUG 
dstrc_enter("s8_mat_linel_cart"); 
#endif
/*----------------------------------------------------------------------*/
l1   = (emod*nue) / ((1.0+nue)*(1.0-2.0*nue));
l2   = emod/ (2.0*(1.0+nue));
ll2  = 2.0*l2;
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
C[0][0][0][0]= l1+ll2;      C[1][0][0][0]= 0.0;          C[2][0][0][0]= 0.0;
C[0][0][0][1]= 0.0;         C[1][0][0][1]= l2 ;          C[2][0][0][1]= 0.0;
C[0][0][0][2]= 0.0;         C[1][0][0][2]= 0.0;          C[2][0][0][2]= l2;
C[0][0][1][0]= 0.0;         C[1][0][1][0]= l2 ;          C[2][0][1][0]= 0.0;
C[0][0][1][1]= l1;          C[1][0][1][1]= 0.0;          C[2][0][1][1]= 0.0;
C[0][0][1][2]= 0.0;         C[1][0][1][2]= 0.0;          C[2][0][1][2]= 0.0;
C[0][0][2][0]= 0.0;         C[1][0][2][0]= 0.0;          C[2][0][2][0]= l2;
C[0][0][2][1]= 0.0;         C[1][0][2][1]= 0.0;          C[2][0][2][1]= 0.0;
C[0][0][2][2]= l1;          C[1][0][2][2]= 0.0;          C[2][0][2][2]= 0.0;
C[0][1][0][0]= 0.0;         C[1][1][0][0]= l1;           C[2][1][0][0]= 0.0;
C[0][1][0][1]= l2 ;         C[1][1][0][1]= 0.0;          C[2][1][0][1]= 0.0;
C[0][1][0][2]= 0.0;         C[1][1][0][2]= 0.0;          C[2][1][0][2]= 0.0;
C[0][1][1][0]= l2 ;         C[1][1][1][0]= 0.0;          C[2][1][1][0]= 0.0;
C[0][1][1][1]= 0.0;         C[1][1][1][1]= l1+ll2;       C[2][1][1][1]= 0.0;
C[0][1][1][2]= 0.0;         C[1][1][1][2]= 0.0;          C[2][1][1][2]= l2;
C[0][1][2][0]= 0.0;         C[1][1][2][0]= 0.0;          C[2][1][2][0]= 0.0;
C[0][1][2][1]= 0.0;         C[1][1][2][1]= 0.0;          C[2][1][2][1]= l2;
C[0][1][2][2]= 0.0;         C[1][1][2][2]= l1;           C[2][1][2][2]= 0.0;
C[0][2][0][0]= 0.0;         C[1][2][0][0]= 0.0;          C[2][2][0][0]= l1;
C[0][2][0][1]= 0.0;         C[1][2][0][1]= 0.0;          C[2][2][0][1]= 0.0;
C[0][2][0][2]= l2 ;         C[1][2][0][2]= 0.0;          C[2][2][0][2]= 0.0;
C[0][2][1][0]= 0.0;         C[1][2][1][0]= 0.0;          C[2][2][1][0]= 0.0;
C[0][2][1][1]= 0.0;         C[1][2][1][1]= 0.0;          C[2][2][1][1]= l1;
C[0][2][1][2]= 0.0;         C[1][2][1][2]= l2 ;          C[2][2][1][2]= 0.0;
C[0][2][2][0]= l2 ;         C[1][2][2][0]= 0.0;          C[2][2][2][0]= 0.0;
C[0][2][2][1]= 0.0;         C[1][2][2][1]= l2 ;          C[2][2][2][1]= 0.0;
C[0][2][2][2]= 0.0;         C[1][2][2][2]= 0.0;          C[2][2][2][2]= l1+ll2;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
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
void s8_kov_CGcuca(DOUBLE T[3][3], const DOUBLE **gkon)
{
INT i,j,k,l;
DOUBLE Tcart[3][3];
/*
DOUBLE c[3][3];
*/
#ifdef DEBUG 
dstrc_enter("s8_kov_CGcuca"); 
#endif
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
for (i=0; i<3; i++)
for (j=0; j<3; j++)
{
   Tcart[0][0] += gkon[0][i]*gkon[0][j]*T[i][j];
   Tcart[0][1] += gkon[0][i]*gkon[1][j]*T[i][j];
   Tcart[0][2] += gkon[0][i]*gkon[2][j]*T[i][j];
   Tcart[1][1] += gkon[1][i]*gkon[1][j]*T[i][j];
   Tcart[1][2] += gkon[1][i]*gkon[2][j]*T[i][j];
   Tcart[2][2] += gkon[2][i]*gkon[2][j]*T[i][j];
}
/*----------------------------------------------------------------------*/
T[0][0] =           Tcart[0][0];   
T[1][0] = T[0][1] = Tcart[0][1];   
T[2][0] = T[0][2] = Tcart[0][2];   
T[1][1] =           Tcart[1][1];   
T[2][1] = T[1][2] = Tcart[1][2];   
T[2][2] =           Tcart[2][2];   
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s8_kov_CGcuca */

#endif
