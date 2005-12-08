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
#include "../headers/standardtypes.h"
#include "shell8.h"
/*----------------------------------------------------------------------*
 | acttime and deltat                                     m.gee 9/03    |
 | global variables of time to be seen from element routines            |
 *----------------------------------------------------------------------*/
extern DOUBLE acttime;
extern DOUBLE deltat;
/*----------------------------------------------------------------------*
 | compressible ogden-material                            m.gee 6/03    |
 | split in volumetric and deviatoric strains                           |
 *----------------------------------------------------------------------*/
void s8_mat_ogden_viscous(ELEMENT    *ele,
                          VISCOHYPER *mat,
                          DOUBLE     *stress_cart,
                          DOUBLE      C_cart[3][3][3][3],
                          DOUBLE    **gkonr,
                          DOUBLE    **gmkovc,
                          INT         gp)
{
      INT         i,j,k,l,p;
const DOUBLE      mthird = -0.3333333333333333333333333333;
const DOUBLE      third  =  0.3333333333333333333333333333;
const DOUBLE      ninth  =  0.1111111111111111111111111111;

      DOUBLE      mu;
      DOUBLE      E;
      DOUBLE      kappa;
      DOUBLE      beta,mbeta;
      DOUBLE      lame1;
      DOUBLE      nue;
      DOUBLE     *mup;
      DOUBLE     *alfap;
      DOUBLE      J;

      /*
      DOUBLE      psi1,psi2,psi;
      */

      DOUBLE      CG[3][3];
      DOUBLE      N[3][3];
      DOUBLE      lam2[3];
      DOUBLE      lam[3];
      DOUBLE      lamdev[3];
      DOUBLE      lamdevpowalfap[3][3];

      DOUBLE      PK2dev[3];
      DOUBLE      PK2vol[3];
      DOUBLE      PK2[3];
      DOUBLE      PK2cart[3][3];

      DOUBLE      scal;
      DOUBLE      Ncross[3];

      DOUBLE      Cdev0000=0.0;
      DOUBLE      Cdev0011=0.0;
      DOUBLE      Cdev0022=0.0;
      DOUBLE      Cdev1111=0.0;
      DOUBLE      Cdev1122=0.0;
      DOUBLE      Cdev2222=0.0;
      DOUBLE      Cdev0101=0.0;
      DOUBLE      Cdev0202=0.0;
      DOUBLE      Cdev1212=0.0;

      DOUBLE      Cvol0000=0.0;
      DOUBLE      Cvol0011=0.0;
      DOUBLE      Cvol0022=0.0;
      DOUBLE      Cvol1111=0.0;
      DOUBLE      Cvol1122=0.0;
      DOUBLE      Cvol2222=0.0;
      DOUBLE      Cvol0101=0.0;
      DOUBLE      Cvol0202=0.0;
      DOUBLE      Cvol1212=0.0;

      DOUBLE      C[3][3][3][3];

/* variables for the viscous effects */

      DOUBLE      nmaxw;
      DOUBLE     *tau;
      DOUBLE     *betas;
      DOUBLE   ***hisold;
      DOUBLE   ***hisnew;
      DOUBLE      xi[4];
      DOUBLE      exi[4];
      DOUBLE      e2xi[4];
      DOUBLE      Qnew[4][3][3];
      DOUBLE      PK2devcart[3][3];
      DOUBLE      onepdelta;
#ifdef DEBUG
dstrc_enter("s8_mat_ogden_viscous");
#endif
/*----------------------------------------------------------------------*/
/*-------------------------------------- init some local arrays to zero */
for (i=0; i<3; i++)
{
   PK2dev[i] = 0.0;
   for (j=0; j<3; j++)
   for (k=0; k<3; k++)
   for (l=0; l<3; l++)
      C[i][j][k][l]      = 0.0;
}
/*------------------------------------------------- init some constants */
if (!(mat->init))
{
   /* make mu = 2.0 * shear modulus */
   mu = 0.0;
   for (i=0; i<3; i++)
   mu += mat->alfap[i] * mat->mup[i];
   /* make Young's modulus E=2mu(1+nue)*/
   E  = mu*(1.0+mat->nue);
   /* make shear modulus */
   mu /= 2.0;
   /* make bulk modulus */
   kappa = mat->kappa = E / (3.0*(1-2.0*(mat->nue)));
   /* make lame constant no. 1 */
   mat->lambda = mat->kappa - (2.0/3.0)*mu;
   /* set init flag */
   mat->init=1;
}
/*-------------------------------------------------- get some constants */
nue        = mat->nue;
beta       = mat->beta;
mbeta      = -beta;
alfap      = mat->alfap;
mup        = mat->mup;
lame1      = mat->lambda;
kappa      = mat->kappa;
nmaxw      = mat->nmaxw; dsassert(nmaxw<=4,"Only 4 Maxwell elements allowed in viscohyper material");
tau        = mat->tau;
betas      = mat->betas;
/*
hisold[0..ngauss-1][0..nmaxw][3][3] converged PK2cartdev and Qalfa of time t-dt
hisnew[0..ngauss-1][0..nmaxw][3][3] current   PK2cartdev and Qalfa of time t
*/
hisold     = ele->e.s8->his1->a.d4[gp];
hisnew     = ele->e.s8->his2->a.d4[gp];
/* scaling factor for deviatoric tangent */
onepdelta  = 1.0;
/* some useful variables */
for (i=0; i<nmaxw; i++)
{
   xi[i]      = -deltat/(2.0*tau[i]);
   exi[i]     = exp(xi[i]);
   e2xi[i]    = exp(2.0*xi[i]);
   onepdelta += betas[i]*exi[i];
}
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
/*------------------------------------- transform CG to cartesian basis */
s8_kov_CGcuca(CG,gkonr);
/*------------------------------------ make spectral decoposition of CG */
s8_ogden_principal_CG(CG,lam2,N);
dsassert(lam2[0]>0.0 && lam2[1]>0.0 && lam2[2]>0.0,"Principal stretches smaller zero");
/*-------------------------------------------- make principal stretches */
lam[0] = sqrt(lam2[0]);
lam[1] = sqrt(lam2[1]);
lam[2] = sqrt(lam2[2]);
/*----------------------------------------------------------------------*/
#if 1
/*---------------- test orthogonality and unit length of eigenvectors N */
/*N0 * N1 = 0*/
scal = N[0][0]*N[1][0] + N[0][1]*N[1][1] + N[0][2]*N[1][2];
dsassert(FABS(scal)<EPS10,"eigenvectors N0,N1 not orthogonal");
/*N0 * N2 = 0*/
scal = N[0][0]*N[2][0] + N[0][1]*N[2][1] + N[0][2]*N[2][2];
dsassert(FABS(scal)<EPS10,"eigenvectors N0,N2 not orthogonal");
/*N1 * N2 = 0*/
scal = N[1][0]*N[2][0] + N[1][1]*N[2][1] + N[1][2]*N[2][2];
dsassert(FABS(scal)<EPS10,"eigenvectors N1,N2 not orthogonal");
/*--------------------------- test proper orientation of eigenvectors N */
/*N2 = N0 x N1*/
Ncross[0] = N[0][1]*N[1][2] - N[0][2]*N[1][1];
Ncross[1] = N[0][2]*N[1][0] - N[0][0]*N[1][2];
Ncross[2] = N[0][0]*N[1][1] - N[0][1]*N[1][0];
/*N2 * Ncross = 1.0*/
scal = Ncross[0]*N[2][0] + Ncross[1]*N[2][1] + Ncross[2]*N[2][2];
dsassert(FABS((scal-1.0))<EPS10,"eigenvectors do not form proper othogonal system");
/*----------------------------------------------------------------------*/
#endif
/*----------------------------------------------------- make J = det(F) */
J = lam[0]*lam[1]*lam[2];
dsassert(J>0.0,"detF <= 0.0 in Ogden material");
/* make pure deviatoric principal stretches by multiplicative split of F */
scal = pow(J,mthird);
lamdev[0] = scal*lam[0];
lamdev[1] = scal*lam[1];
lamdev[2] = scal*lam[2];
/*----------------------------------- make powers lamdev[i] ** alpfa[p] */
for (i=0; i<3; i++)
for (p=0; p<3; p++)
   lamdevpowalfap[i][p] = pow(lamdev[i],alfap[p]);
/*--------------------------------------------------------- make energy */
#if 0
psi1 = 0.0;
for (p=0; p<3; p++)
{
   psi1 += (mup[p]/alfap[p]) *
           (lamdevpowalfap[0][p] +
           lamdevpowalfap[1][p] +
           lamdevpowalfap[2][p] - 3.0);
}
psi2 = (kappa/(beta*beta)) *
       (beta*log(J) + pow(J,mbeta)-1.0);
psi = psi1+psi2;
printf("uncoupled PSI1 %20.10f PSI2 %20.10f PSI %20.10f\n\n",psi1,psi2,psi);fflush(stdout);
#endif
/*-------------------------------------------- make deviatoric stresses */
for (p=0; p<3; p++)
{
   scal = mthird * (lamdevpowalfap[0][p]+lamdevpowalfap[1][p]+ lamdevpowalfap[2][p]);
   for (i=0; i<3; i++)
      PK2dev[i] += mup[p]*(lamdevpowalfap[i][p]+scal);
}
for (i=0; i<3; i++)
   PK2dev[i] /= lam2[i];
/*-------------------------------------------- make volumetric stresses */
scal = (kappa/(beta))*(1.0-pow(J,mbeta));
for (i=0; i<3; i++)
   PK2vol[i] = scal/lam2[i];
/*---------------- make total stresses and transform to cartesian bases */
for (i=0; i<3; i++) PK2[i] = PK2dev[i] + PK2vol[i];
s8_ogden_cartPK2(PK2cart,PK2,N);
/*---------------------------------- make cartesian deviatoric stresses */
s8_ogden_cartPK2(PK2devcart,PK2dev,N);
/*-------------------- put the cartesian deviatoric stresses to storage */
for (i=0; i<3; i++)
for (j=0; j<3; j++)
   hisnew[0][i][j] = PK2devcart[i][j];
/*--------------------------------- make non-equillibrium stresses Qnew */
/*                                this is all done in cartesian bases ! */
for (i=0; i<nmaxw; i++)
{
   for (j=0; j<3; j++)
   for (k=0; k<3; k++)
      Qnew[i][j][k] = betas[i]*exi[i]*PK2devcart[j][k] +
                      e2xi[i] * hisold[i+1][j][k]      -
                      betas[i]*exi[i]*hisold[0][j][k];
}
/*-------------------------------------------- put new Qalfa to storage */
for (i=0; i<nmaxw; i++)
for (j=0; j<3; j++)
for (k=0; k<3; k++)
   hisnew[i+1][j][k] = Qnew[i][j][k];
/*----------------------------------------------------------------------*/
#if 0
printf("PK2main_dev[0] %14.8f PK2main_dev[1] %14.8f PK2main_dev[2] %14.8f\n",PK2dev[0],PK2dev[1],PK2dev[2]);
printf("PK2main_vol[0] %14.8f PK2main_vol[1] %14.8f PK2main_vol[2] %14.8f\n",PK2vol[0],PK2vol[1],PK2vol[2]);
printf("PK2        [0] %14.8f PK2        [1] %14.8f PK2        [2] %14.8f\n\n",PK2[0],PK2[1],PK2[2]);
#endif
/*---------------------- add the viscous stresses to the total stresses */
for (i=0; i<nmaxw; i++)
for (j=0; j<3; j++)
for (k=0; k<3; k++)
   PK2cart[j][k] += Qnew[i][j][k];
/*--------------------------------- sort PK2cart to vector shell8 style */
/*--- which then will be tranformed to shell bases outside this routine */
stress_cart[0] = PK2cart[0][0];
stress_cart[1] = PK2cart[0][1];
stress_cart[2] = PK2cart[0][2];
stress_cart[3] = PK2cart[1][1];
stress_cart[4] = PK2cart[1][2];
stress_cart[5] = PK2cart[2][2];
/*---------------------- make deviatoric material tangent in eigenspace */
/*================== components Cdev_aaaa */
for (p=0; p<3; p++)
{
   scal = ninth*(lamdevpowalfap[0][p]+lamdevpowalfap[1][p]+ lamdevpowalfap[2][p]);

   Cdev0000 += alfap[p]*mup[p]*(third*lamdevpowalfap[0][p]+scal);
   Cdev1111 += alfap[p]*mup[p]*(third*lamdevpowalfap[1][p]+scal);
   Cdev2222 += alfap[p]*mup[p]*(third*lamdevpowalfap[2][p]+scal);
}
Cdev0000 /= (lam2[0]*lam2[0]);
Cdev1111 /= (lam2[1]*lam2[1]);
Cdev2222 /= (lam2[2]*lam2[2]);
/*--------- part missing in holzapfel book 'Nonlinear Solid Mechanics' */
Cdev0000 += (-2.0)*PK2dev[0]/lam2[0];
Cdev1111 += (-2.0)*PK2dev[1]/lam2[1];
Cdev2222 += (-2.0)*PK2dev[2]/lam2[2];
/*================== components Cdev_aabb */
for (p=0; p<3; p++)
{
   scal = ninth*(lamdevpowalfap[0][p]+lamdevpowalfap[1][p]+ lamdevpowalfap[2][p]);

   Cdev0011 += alfap[p]*mup[p]*(mthird*lamdevpowalfap[0][p]+mthird*lamdevpowalfap[1][p]+scal);
   Cdev0022 += alfap[p]*mup[p]*(mthird*lamdevpowalfap[0][p]+mthird*lamdevpowalfap[2][p]+scal);
   Cdev1122 += alfap[p]*mup[p]*(mthird*lamdevpowalfap[1][p]+mthird*lamdevpowalfap[2][p]+scal);
}
Cdev0011 /= (lam2[0]*lam2[1]);
Cdev0022 /= (lam2[0]*lam2[2]);
Cdev1122 /= (lam2[1]*lam2[2]);
/*================== components Cdev_abab */
if (FABS(lam2[0]-lam2[1])>EPS12)
Cdev0101 = (PK2dev[0]-PK2dev[1])/(lam2[0]-lam2[1]);
else
Cdev0101 = 0.5*(Cdev0000-Cdev0011);

if (FABS(lam2[0]-lam2[2])>EPS12)
Cdev0202 = (PK2dev[0]-PK2dev[2])/(lam2[0]-lam2[2]);
else
Cdev0202 = 0.5*(Cdev0000-Cdev0022);

if (FABS(lam2[1]-lam2[2])>EPS12)
Cdev1212 = (PK2dev[1]-PK2dev[2])/(lam2[1]-lam2[2]);
else
Cdev1212 = 0.5*(Cdev1111-Cdev1122);
/*---------------------- make volumetric material tangent in eigenspace */
/*================== components Cvol_aaaa */
scal = kappa*((2.0/beta+1)*pow(J,mbeta)-2.0/beta);
Cvol0000 = scal / (lam2[0]*lam2[0]);
Cvol1111 = scal / (lam2[1]*lam2[1]);
Cvol2222 = scal / (lam2[2]*lam2[2]);

/*================== components Cvol_aabb */
scal = kappa*pow(J,mbeta);
Cvol0011 =  scal / (lam2[0]*lam2[1]);
Cvol0022 =  scal / (lam2[0]*lam2[2]);
Cvol1122 =  scal / (lam2[1]*lam2[2]);

/*================== components Cvol_abab */
if (FABS(lam2[0]-lam2[1])>EPS12)
Cvol0101 = (PK2vol[0]-PK2vol[1])/(lam2[0]-lam2[1]);
else
Cvol0101 = 0.5*(Cvol0000-Cvol0011);

if (FABS(lam2[0]-lam2[2])>EPS12)
Cvol0202 = (PK2vol[0]-PK2vol[2])/(lam2[0]-lam2[2]);
else
Cvol0202 = 0.5*(Cvol0000-Cvol0022);

if (FABS(lam2[1]-lam2[2])>EPS12)
Cvol1212 = (PK2vol[1]-PK2vol[2])/(lam2[1]-lam2[2]);
else
Cvol1212 = 0.5*(Cvol1111-Cvol1122);
/*------------ put everything together and transform to cartesian bases */
C[0][0][0][0] = onepdelta*Cdev0000 + Cvol0000;
C[1][1][1][1] = onepdelta*Cdev1111 + Cvol1111;
C[2][2][2][2] = onepdelta*Cdev2222 + Cvol2222;

C[1][1][0][0] = C[0][0][1][1] = onepdelta*Cdev0011 + Cvol0011;
C[2][2][0][0] = C[0][0][2][2] = onepdelta*Cdev0022 + Cvol0022;
C[2][2][1][1] = C[1][1][2][2] = onepdelta*Cdev1122 + Cvol1122;

C[1][0][1][0] = C[0][1][0][1] = onepdelta*Cdev0101 + Cvol0101;
C[2][0][2][0] = C[0][2][0][2] = onepdelta*Cdev0202 + Cvol0202;
C[2][1][2][1] = C[1][2][1][2] = onepdelta*Cdev1212 + Cvol1212;

s8_ogden_Ccart(C,C_cart,N);

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s8_mat_ogden_viscous */



#endif
