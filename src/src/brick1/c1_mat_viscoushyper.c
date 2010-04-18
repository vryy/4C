/*!-----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer:	Robert Metzke
		metzke@lnm.mw.tum.de
		http://www.lnm.mw.tum.de/Members/metzke
		089 - 289 15244
</pre>
------------------------------------------------------------------------*/

#ifndef CCADISCRET

#include "../headers/standardtypes.h"
#include "brick1.h"
#include "brick1_prototypes.h"
#ifdef D_MAT
        #include "../materials/mat_prototypes.h"
#endif
/*----------------------------------------------------------------------*
 | acttime and deltat                                          rm 01/07 |
 | global variables of time to be seen from element routines            |
 *----------------------------------------------------------------------*/
extern DOUBLE acttime;
extern DOUBLE deltat;
/*----------------------------------------------------------------------*
 | compressible ogden-material                            m.gee 6/03    |
 | split in volumetric and deviatoric strains                           |
 | adapted from shell 8 to brick                               rm 03/07 |
 *----------------------------------------------------------------------*/
void c1_mat_ogden_viscous(  ELEMENT    *ele,
                            INT ip,
                            VISCOHYPER *mat,
                            DOUBLE *disd,
                            DOUBLE *stress,
                            DOUBLE **d)
{
#ifdef D_BRICK1

      INT i,j,k,l,p;
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

      DOUBLE	  FT[3][3];
      DOUBLE      lam2[3];
      DOUBLE      lam[3];
      DOUBLE      lamdev[3];
      DOUBLE      lamdevpowalfap[3][3];

      DOUBLE      PK2dev[3];
      DOUBLE      PK2vol[3];
      DOUBLE      PK2[3];

      DOUBLE      scal;

      static DOUBLE **CG;
      static ARRAY CG_a;
      static DOUBLE **N;
      static ARRAY N_a;

      static DOUBLE **PK2cart;
      static ARRAY PK2cart_a;
      static DOUBLE **PK2devcart;
      static ARRAY PK2devcart_a;

      if (CG==NULL)
      {
	      CG = amdef("CG",&CG_a,3,3,"DA");
	      N = amdef("N",&N_a,3,3,"DA");

	      PK2cart = amdef("PK2cart",&PK2cart_a,3,3,"DA");
	      PK2devcart = amdef("PK2devcart",&PK2devcart_a,3,3,"DA");
      }

      amzero(&CG_a);
      amzero(&N_a);
      amzero(&PK2cart_a);
      amzero(&PK2devcart_a);


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

      DOUBLE      Ceigen[3][3][3][3];
      DOUBLE      Ccart[3][3][3][3];
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
      DOUBLE      onepdelta;

#ifdef DEBUG
	dstrc_enter("c1_mat_ogden_viscous");
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

/*hisold     = ele->e.c1->his1->a.d4[ip];*/
hisold = ele->e.c1->his1->a.d4[ip];
hisnew = ele->e.c1->his2->a.d4[ip];
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

/*--------------------------------------Deformation Gradient, transposed*/
FT[0][0]=disd[0]+1.0;
FT[1][1]=disd[1]+1.0;
FT[2][2]=disd[2]+1.0;
FT[0][1]=disd[3];/*---------------------disd[4];*/
FT[1][0]=disd[4];/*---------------------disd[3];*/
FT[1][2]=disd[5];/*---------------------disd[6];*/
FT[2][1]=disd[6];/*---------------------disd[5];*/
FT[0][2]=disd[7];/*---------------------disd[8];*/
FT[2][0]=disd[8];/*---------------------disd[7];*/
/*--------------------------------------------Right Cauchy Green Tensor*/
for (k=0; k<3; k++) {
	for (i=0; i<3; i++) {
		for (j=0; j<3; j++) {
			CG[i][k]+=(FT[i][j]*FT[k][j]);
		}
	}
}
/*------------------------------------ make spectral decoposition of CG */
/*c1_ogden_principal_CG(CG,lam2,N);*/
c1_calc_eigenval_eigenvec_jacobi(CG,lam2,N);
dsassert(lam2[0]>0.0 && lam2[1]>0.0 && lam2[2]>0.0,"Principal stretches smaller zero");
/*-------------------------------------------- make principal stretches */
lam[0] = sqrt(lam2[0]);
lam[1] = sqrt(lam2[1]);
lam[2] = sqrt(lam2[2]);
/*----------------------------------------------------------------------*/
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
	for (i=0; i<3; i++) PK2dev[i] += mup[p]*(lamdevpowalfap[i][p]+scal);
}
for (i=0; i<3; i++) PK2dev[i] /= lam2[i];
/*-------------------------------------------- make volumetric stresses */
scal = (kappa/(beta))*(1.0-pow(J,mbeta));
for (i=0; i<3; i++) PK2vol[i] = scal/lam2[i];
/*------------------------------------------------- make total stresses */
for (i=0; i<3; i++) PK2[i] = PK2dev[i] + PK2vol[i];
#ifdef D_MAT
mat_ogden_cartPK2(PK2cart,PK2,N);
/*---------------------------------- make cartesian deviatoric stresses */
mat_ogden_cartPK2(PK2devcart,PK2dev,N);
#else
dserror("Please use D_MAT in your configure file, mate!");
#endif
/*-------------------- put the cartesian deviatoric stresses to storage */
for (i=0; i<3; i++)
for (j=0; j<3; j++)
   hisnew[0][i][j] = PK2devcart[i][j];
/*--------------------------------- make non-equillibrium stresses Qnew */
for (i=0; i<nmaxw; i++)
{
	for (j=0; j<3; j++)
	for (k=0; k<3; k++)
	Qnew[i][j][k] = betas[i]*exi[i]*PK2devcart[j][k] +
		e2xi[i] * hisold[i+1][j][k] -
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
/*----------------------------------- sort PK2cart to vector brick style */
stress[0] = PK2cart[0][0];
stress[1] = PK2cart[1][1];
stress[2] = PK2cart[2][2];
stress[3] = PK2cart[0][1];
stress[4] = PK2cart[1][2];
stress[5] = PK2cart[0][2];
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
/*-------------------------------------------- put everything together */
Ceigen[0][0][0][0] = onepdelta*Cdev0000 + Cvol0000;
Ceigen[1][1][1][1] = onepdelta*Cdev1111 + Cvol1111;
Ceigen[2][2][2][2] = onepdelta*Cdev2222 + Cvol2222;

Ceigen[1][1][0][0] = Ceigen[0][0][1][1] = onepdelta*Cdev0011 + Cvol0011;
Ceigen[2][2][0][0] = Ceigen[0][0][2][2] = onepdelta*Cdev0022 + Cvol0022;
Ceigen[2][2][1][1] = Ceigen[1][1][2][2] = onepdelta*Cdev1122 + Cvol1122;

Ceigen[1][0][1][0] = Ceigen[0][1][0][1] = onepdelta*Cdev0101 + Cvol0101;
Ceigen[2][0][2][0] = Ceigen[0][2][0][2] = onepdelta*Cdev0202 + Cvol0202;
Ceigen[2][1][2][1] = Ceigen[1][2][1][2] = onepdelta*Cdev1212 + Cvol1212;
#ifdef D_MAT
mat_ogden_Ccart(Ceigen,Ccart,N);
#else
dserror("Please use D_MAT in your configure file, mate!");
#endif

for (i=0; i<3; i++)
for (j=0; j<3; j++)
for (k=0; k<3; k++)
for (l=0; l<3; l++)
{
C[i][j][k][l] = Ccart[i][j][k][l];
}
/*------------------------------------------------------- write C in d */
d[0][0]=C[0][0][0][0];
d[0][1]=C[0][0][1][1];
d[0][2]=C[0][0][2][2];
d[0][3]=C[0][0][1][0];
d[0][4]=C[0][0][2][1];
d[0][5]=C[0][0][2][0];

d[1][0]=C[1][1][0][0];
d[1][1]=C[1][1][1][1];
d[1][2]=C[1][1][2][2];
d[1][3]=C[1][1][1][0];
d[1][4]=C[1][1][2][1];
d[1][5]=C[1][1][2][0];

d[2][0]=C[2][2][0][0];
d[2][1]=C[2][2][1][1];
d[2][2]=C[2][2][2][2];
d[2][3]=C[2][2][1][0];
d[2][4]=C[2][2][2][1];
d[2][5]=C[2][2][2][0];

d[3][0]=C[1][0][0][0];
d[3][1]=C[1][0][1][1];
d[3][2]=C[1][0][2][2];
d[3][3]=C[1][0][1][0];
d[3][4]=C[1][0][2][1];
d[3][5]=C[1][0][2][0];

d[4][0]=C[2][1][0][0];
d[4][1]=C[2][1][1][1];
d[4][2]=C[2][1][2][2];
d[4][3]=C[2][1][1][0];
d[4][4]=C[2][1][2][1];
d[4][5]=C[2][1][2][0];

d[5][0]=C[2][0][0][0];
d[5][1]=C[2][0][1][1];
d[5][2]=C[2][0][2][2];
d[5][3]=C[2][0][1][0];
d[5][4]=C[2][0][2][1];
d[5][5]=C[2][0][2][0];


/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
#endif /* end of ifdef BRICK1 */
return;
} /* end of c1_mat_ogden_viscous */

#endif
