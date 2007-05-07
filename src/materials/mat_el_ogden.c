/*---------------------------------------------------------------------                                    
\file
\brief

<pre>
Maintainer:     Robert Metzke
                metzke@lnm.mw.tum.de
                http://www.lnm.mw.tum.de/Members/metzke
                089 - 289 15244
</pre>
*---------------------------------------------------------------------*/
#ifdef D_MAT
#include "../headers/standardtypes.h"
#include "mat_prototypes.h"

void mat_el_ogden_uncoupled (
		COMPOGDEN *mat,
		DOUBLE lam[3],
		DOUBLE N[3][3],
		DOUBLE stress[3][3],
		DOUBLE C[3][3][3][3])
{
	INT i,j,k,l,p;
	const   DOUBLE      mthird = -0.3333333333333333333333333333;
	const   DOUBLE      third  =  0.3333333333333333333333333333;
	const   DOUBLE      ninth  =  0.1111111111111111111111111111;
		
	DOUBLE  mu;
	DOUBLE  E;
	DOUBLE  kappa;
	DOUBLE  beta,mbeta;
	DOUBLE  lame1;
	DOUBLE  nue;
	DOUBLE  *mup;
	DOUBLE  *alfap;
	DOUBLE  J;
	DOUBLE  scal;
	
	DOUBLE  Neigen[3][3];
	DOUBLE  Ncross[3];
	DOUBLE  lam2[3];
	DOUBLE  lamdev[3];
	DOUBLE  lamdevpowalfap[3][3];

        /*
	DOUBLE  psi1,psi2,psi;
	*/
	
	DOUBLE      PK2dev[3];
	DOUBLE      PK2vol[3];
	DOUBLE      PK2[3];
	DOUBLE      PK2cart[3][3];

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

#ifdef DEBUG
	dstrc_enter("c1_mat_ogden_uncoupled");
#endif
/*-------------------------------------- init some local arrays to zero */
	for (i=0; i<3; i++) {
	Ncross[i] = 0.0;
	lam2[i] = 0.0;
	lamdev[i] = 0.0;
	PK2dev[i] = 0.0;
	PK2vol[i] = 0.0;
	PK2[i] = 0.0;
	for (j=0; j<3; j++) {
		Neigen[i][j] = 0.0;
		lamdevpowalfap[i][j] = 0.0;
		PK2cart[i][j] = 0.0;
		stress[i][j] = 0.0;
		for (k=0; k<3; k++)
			for (l=0; l<3; l++) {
				Ceigen[i][j][k][l] = 0.0;
				C[i][j][k][l] = 0.0;
			}
	}
	}
/*------------------------------------------------- init some constants */
	if (!(mat->init))
	{
		/* make mu = 2.0 * shear modulus */
		mu = 0.0;
		for (i=0; i<3; i++) mu += mat->alfap[i] * mat->mup[i];
		/* make Young's modulus */
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
	nue       = mat->nue;
	beta      = mat->beta;
	mbeta     = -beta;
	alfap     = mat->alfap;
	mup       = mat->mup;
	lame1     = mat->lambda;
	kappa     = mat->kappa;
/*-------------------------------------------- test principal stretches */
	dsassert(lam[0]>0.0 && lam[1]>0.0 && lam[2]>0.0,"Principal stretches smaller zero");
/*----------------------------------------- make principal stretchesi^2 */
	lam2[0] = pow(lam[0],2);
	lam2[1] = pow(lam[1],2);
	lam2[2] = pow(lam[2],2);
/*---------------- test orthogonality and unit length of eigenvectors N */
#if 1
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
#endif
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
	psi1 += (mup[p]/alfap[p]) * (lamdevpowalfap[0][p] + lamdevpowalfap[1][p] + lamdevpowalfap[2][p] - 3.0);
	}
	psi2 = (kappa/(beta*beta)) * (beta*log(J) + pow(J,mbeta)-1.0);
	psi = psi1+psi2;
	printf("uncoupled PSI1 %20.10f PSI2 %20.10f PSI %20.10f\n\n",psi1,psi2,psi);fflush(stdout);
#endif
/*----------------------------------------------------------------------*/
/*-------------------------------------------------------------- stress */
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
/*------------------------------------------------- make total stresses */
	for (i=0; i<3; i++) PK2[i] = PK2dev[i] + PK2vol[i];
	mat_ogden_cartPK2(PK2cart,PK2,N);
/*----------------------------------------------------------------------*/
#if 0
printf("PK2main_dev[0] %14.8f PK2main_dev[1] %14.8f PK2main_dev[2] %14.8f\n",PK2dev[0],PK2dev[1],PK2dev[2]);
printf("PK2main_vol[0] %14.8f PK2main_vol[1] %14.8f PK2main_vol[2] %14.8f\n",PK2vol[0],PK2vol[1],PK2vol[2]);
printf("PK2        [0] %14.8f PK2        [1] %14.8f PK2        [2] %14.8f\n\n",PK2[0],PK2[1],PK2[2]);
#endif
/*------------------------------------------------------- sort stresses */
	stress[0][0] = PK2cart[0][0];
	stress[1][1] = PK2cart[1][1];
	stress[2][2] = PK2cart[2][2];
	stress[0][1] = PK2cart[0][1];
	stress[1][2] = PK2cart[1][2];
	stress[0][2] = PK2cart[0][2];
/*----------------------------------------------------------------------*/
/*--------------------------------------------------- elasticity tensor */
/*---------------------- make deviatoric material tangent in eigenspace */
/*================== components Cdev_aaaa */
	for(p=0; p<3; p++)
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
	Ceigen[0][0][0][0] = Cdev0000+Cvol0000;
	Ceigen[1][1][1][1] = Cdev1111+Cvol1111;
	Ceigen[2][2][2][2] = Cdev2222+Cvol2222;

	Ceigen[1][1][0][0] = Ceigen[0][0][1][1] = Cdev0011+Cvol0011;
	Ceigen[2][2][0][0] = Ceigen[0][0][2][2] = Cdev0022+Cvol0022;
	Ceigen[2][2][1][1] = Ceigen[1][1][2][2] = Cdev1122+Cvol1122;

	Ceigen[1][0][1][0] = Ceigen[0][1][0][1] = Cdev0101+Cvol0101;
	Ceigen[2][0][2][0] = Ceigen[0][2][0][2] = Cdev0202+Cvol0202;
	Ceigen[2][1][2][1] = Ceigen[1][2][1][2] = Cdev1212+Cvol1212;

	mat_ogden_Ccart(Ceigen,Ccart,N);

	for (i=0; i<3; i++)
	for (j=0; j<3; j++)
	for (k=0; k<3; k++)
	for (l=0; l<3; l++)
	{
	C[i][j][k][l] = Ccart[i][j][k][l];
	}
	
#ifdef DEBUG
	dstrc_exit();
#endif
	return;
} /* end of mat_el_ogden_uncoupled */



/*----------------------------------------------------------------------*
 * transform principal streeses PK2 in cartesian stresses PK2 m.gee 6/03|
 * PK2 = PK2main_a * N[][a] dyad N[][a] (sum over a)                    |
 *----------------------------------------------------------------------*/
void mat_ogden_cartPK2(DOUBLE PK2[3][3], DOUBLE PK2main[3], DOUBLE N[3][3])
{
	INT i,j;
	DOUBLE dyad0[3][3],dyad1[3][3],dyad2[3][3];
#ifdef DEBUG
	dstrc_enter("c1_ogden_cartPK2");
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
	
#ifdef DEBUG
	dstrc_exit();
#endif
	return;
} /* end of mat_ogden_cartPK2 */



/*----------------------------------------------------------------------*
 * transform C in cartesian bases                         PK2 m.gee 6/03|
 *----------------------------------------------------------------------*/
void mat_ogden_Ccart(DOUBLE C[3][3][3][3], DOUBLE C_cart[3][3][3][3], DOUBLE N[3][3])
{
	INT i,j,k,l;
#ifdef DEBUG
	dstrc_enter("c1_ogden_Ccart");
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
} /* end of mat_ogden_Ccart */

					
	

#endif /*D_MAT*/
