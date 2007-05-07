/*---------------------------------------------------------------------                                    
\file
\brief contains the routine
 - c1_mat_ogden_uncoupled:	calculates the principal stretches and
 				the eigenvectors of right Cauchy Green
				and calls the decoupled ogden material
				law at the materials folder and gives
				back the stress and constitutive tensor
 - c1_ogden_principal_CG:	prepares the right Cauchy Green tensor
 				for spectral decomposition
<pre>
Maintainer:		Robert Metzke
 			metzke@lnm.mw.tum.de
			http://www.lnm.mw.tum.de/Members/metzke
			089 - 289 15244
</pre>

*-------------------------------------------------------------------*/
#ifdef D_BRICK1
#include "../headers/standardtypes.h"
#include "brick1.h"
#include "brick1_prototypes.h"
/*----------------------------------------------------------------------*
 * compressible ogden-material                               m.gee 6/03 |
 * split in volumetric and deviatoric strains                           |
 * adapted from shell 8 to brick                               rm 03/07 |
 *----------------------------------------------------------------------*/
void c1_mat_ogden_uncoupled(
		 COMPOGDEN *mat,
		 DOUBLE *disd,
		 DOUBLE *stress,
		 DOUBLE **d)
{
	INT i,j,k;
	
	DOUBLE  FT[3][3];
	DOUBLE  CG[3][3];
	DOUBLE  N[3][3];
	DOUBLE  lam[3];
	DOUBLE  lam2[3];

	DOUBLE  PK2[3][3];
	DOUBLE  C[3][3][3][3];

#ifdef DEBUG
	dstrc_enter("c1_mat_ogden_uncoupled");
#endif
/*-------------------------------------- init some local arrays to zero */
	for (i=0; i<3; i++) {
		lam[i] = 0.0;
		lam2[i] = 0.0;
		for (j=0; j<3; j++) {
			N[i][j] = 0.0;
			FT[i][j] = 0.0;
			CG[i][j] = 0.0;
		}
	}
/*----------------------------------- Deformation Gradient, transposed */
	FT[0][0]=disd[0]+1.0;
	FT[1][1]=disd[1]+1.0;
	FT[2][2]=disd[2]+1.0;
	FT[0][1]=disd[3];/*---------------------disd[4];*/
	FT[1][0]=disd[4];/*---------------------disd[3];*/
	FT[1][2]=disd[5];/*---------------------disd[6];*/
	FT[2][1]=disd[6];/*---------------------disd[5];*/
	FT[0][2]=disd[7];/*---------------------disd[8];*/
	FT[2][0]=disd[8];/*---------------------disd[7];*/
/*------------------------------------------- Right Cauchy Green Tensor*/
	for (k=0; k<3; k++) {
		for (i=0; i<3; i++) {
			for (j=0; j<3; j++) {
				CG[i][k]+=(FT[i][j]*FT[k][j]); } } }
/*------------------------------------ make spectral decoposition of CG */
	c1_ogden_principal_CG(CG,lam2,N);
/*-------------------------------------------- make principal stretches */
	lam[0] = sqrt(lam2[0]);
	lam[1] = sqrt(lam2[1]);
	lam[2] = sqrt(lam2[2]);
#if 1
	for (i=0; i<3; i++) mat->l[i] = lam[i];
#endif
/*------------------------------ call ogden material in material folder */
	mat_el_ogden_uncoupled(mat,lam,N,PK2,C);
/*------------------------------------------------------------ Stresses */
	stress[0] = PK2[0][0];
	stress[1] = PK2[1][1];
	stress[2] = PK2[2][2];
	stress[3] = PK2[0][1];
	stress[4] = PK2[1][2];
	stress[5] = PK2[0][2];
/*---------------------------------------------------Constitutive Matrix*/
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
					    
#ifdef DEBUG
	dstrc_exit();
#endif
	return;
} /* end of c1_mat_ogden_uncoupled */




/*----------------------------------------------------------------------*
 * make eigenvalue decomposition of Cauchy-Green strains  m.gee 6/03    |
 * adapted to brick1                                      rm 03.07      |
 *----------------------------------------------------------------------*/
void c1_ogden_principal_CG(DOUBLE CG[3][3], DOUBLE lambda[3], DOUBLE N[3][3])
{
	INT          i;
	DOUBLE       fstrain[9];
	DOUBLE       fn[9];
#ifdef DEBUG
	dstrc_enter("c1_ogden_principal_CG");
#endif
/*----------------------------------------------------------------------*/
	for (i=0; i<9; i++) fn[i] = 0.0;

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
} /* end of c1_ogden_principal_CG */

#endif

