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
#ifdef D_MAT
	#include "../materials/mat_prototypes.h"
#endif
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
#ifdef D_MAT
	    mat_el_ogden_uncoupled(mat,lam,N,PK2,C);
#else
            dserror("Please use D_MAT, mate!");
#endif 
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
	INT          i,j;
	DOUBLE       fstrain[9];
	DOUBLE       fn[9];
	DOUBLE       lam2[3];
	static DOUBLE **CGt;
 	static ARRAY CGt_a;
	static DOUBLE **Nt;
	static ARRAY Nt_a;
	
	CGt = amdef("CGt",&CGt_a,3,3,"DA");
	Nt = amdef("Nt",&Nt_a,3,3,"DA");
 
#ifdef DEBUG
	dstrc_enter("c1_ogden_principal_CG");
#endif
/*----------------------------------------------------------------------*/

	for(i=0; i<3; i++) {
		for(j=0; j<3; j++) {
			CGt[i][j] = CG[i][j];
		}
	}
	
	c1_calc_eigenval_eigenvec_jacobi(CGt,lam2,Nt);

	lambda[0] = lam2[0];
	lambda[1] = lam2[1];
	lambda[2] = lam2[2];

	for(i=0; i<3; i++) {
		for(j=0; j<3; j++) {
			N[i][j] = Nt[i][j];
		}
	}
	
/*
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
*/
/*----------------------------------------------------------------------*/
#ifdef DEBUG
	dstrc_exit();
#endif
return;
} /* end of c1_ogden_principal_CG */


/*----------------------------------------------------------------------*
 * eigenvalues and eigenvectors of a matrix, applying the jacobi method |
 * taken from Numerical Recipes and specialized to 3x3 case    rm 03.07 |
 *----------------------------------------------------------------------*/
void c1_calc_eigenval_eigenvec_jacobi (DOUBLE **C, DOUBLE *d, DOUBLE **V)
{
#ifdef DEBUG
	dstrc_enter("c1_calc_eigenval_eigenvec_jacobi");
#endif
	INT n=3;
	INT j,iq,ip,i;
	DOUBLE tresh,theta,tau,t,sm,s,h,g,c;
	/*DOUBLE C_loc[3][3];*/
	DOUBLE **C_loc;
	ARRAY C_loc_a;
	
	C_loc = amdef("C_loc",&C_loc_a,3,3,"DA");
	DOUBLE b[n];
	DOUBLE z[n];
/*------------------------------------------------------ Einheitsmatrix */
	for (i=0; i<3; i++) {
		for (j=0; j<3; j++) {
			C_loc[i][j]=C[i][j];
		}
	}
	for (ip=0; ip<n; ip++) {
		for (iq=0; iq<n; iq++) {
			V[ip][iq]=0.0;
			V[ip][ip]=1.0;
		}
	}
/*---------------------- vectors b,d = diagonal elements of A, z = NULL */
	for (ip=0;ip<n;ip++) {
		b[ip]=d[ip]=C_loc[ip][ip];
		z[ip]=0.0;
	}
/*------------------------------------------------------- 50 Iterationen*/
	for (i=1;i<=50;i++) {
		sm=0.0;
		for (ip=0;ip<n-1;ip++) {
			for (iq=ip+1;iq<n;iq++) {
				sm += fabs(C_loc[ip][iq]);
			}
		}
		if (sm == 0.0) {
			return;
		}
		if (i < 4) 
			tresh=0.2*sm/(n*n);
		else
			tresh=0.0;
		for (ip=0;ip<n-1;ip++) {
			for (iq=ip+1;iq<n;iq++) {
				g=100.0 * fabs(C_loc[ip][iq]);
				if (i > 4 && (fabs(d[ip])+g) == fabs(d[ip]) && (fabs(d[iq])+g) == fabs(d[iq]))
					C_loc[ip][iq]=0.0;
				else if (fabs(C_loc[ip][iq]) > tresh){
					h=d[iq]-d[ip];
					if ((fabs(h)+g) == fabs(h))
						t=(C_loc[ip][iq])/h;
					else {
						theta=0.5f*h/(C_loc[ip][iq]);
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0f)
							t = -t;
					}
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*C[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					C_loc[ip][iq]=0.0;
					for (j=0;j<=ip-1;j++) {
						c1_rotation(C_loc,j,ip,j,iq,tau,s);
					}
					for (j=ip+1;j<=iq-1;j++) {
						c1_rotation(C_loc,ip,j,j,iq,tau,s);
					}
					for (j=iq+1;j<n;j++) {
						c1_rotation(C_loc,ip,j,iq,j,tau,s);
					}
					for (j=0;j<n;j++) {
						c1_rotation(V,j,ip,j,iq,tau,s);
					}
				}
			}
		}
		for (ip=0;ip<n;ip++) {
			b[ip] += z[ip];
			d[ip]=b[ip];
			z[ip]=0.0;
		}
	}
}


void c1_rotation (DOUBLE **C, INT i, INT j, INT k, INT l, DOUBLE tau, DOUBLE s)
{
#ifdef DEBUG
	dstrc_enter("c1_rotation");
#endif
	DOUBLE g,h;
	g=C[i][j];
	h=C[k][l];
	C[i][j]=g-s*(h+g*tau);
	C[k][l]=h+s*(g-h*tau);
#ifdef DEBUG
	dstrc_exit();
#endif
	return;
}



#endif

