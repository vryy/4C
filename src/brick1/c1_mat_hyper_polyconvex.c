/*-----------------------------------------------------------------------------------------------------*/
/*													|
\file													|
\brief contains the routine 'c1_mat_hyper_polyconvex' to establish local material law,			|
       stress-strain relationship for hyperelastic, anisotropic material for a 3D hex element		|
       used for biological, soft, collagenous tissues.							|
       													|
<pre>													|
Maintainer:												|
													|
</pre>													|
													|
*/
/*-----------------------------------------------------------------------------------------------------*/
#ifdef D_BRICK1

#include "../headers/standardtypes.h"
#include "brick1.h"
#include "brick1_prototypes.h"
/*-----------------------------------------------------------------------------------------------------*/
/*                                            								|
\brief establish local material										|
													|
<pre>                                                              br 09/06				|
This routine establishes a local material law,								|
stress-strain relationship for hyperelastic, anisotropic material for a 3D-hex-element.			|
Used for biological, soft, collagenous tissues.								|
													|
Documented in the 'Diplomarbeit' by Barbara Roehrnbauer at LNM.						|
Based on Holzapfel [1], Ogden [2] and Balzani, Schroeder, Neff [3].					|
													|
[1] G.A.Holzapfel, R.W.Ogden, A New Consitutive Framework for Arterial Wall Mechanics and		|
	a Comparative Study of Material Models, Journal of Elasticity 61, 1-48, 2000.			|
[2] R.W.Ogden, Anisotropy and Nonlinear Elasticity in Arterial Wall Mechanics,				|
	CISM Course on Biomechanical Modeling, Lectures 2,3, 2006.					|
[3] D.Balzani, P.Neff, J.Schroeder, G.A.Holzapfel, A Polyconvex Framework for Soft Biological Tissues	|
	Adjustment to Experimental Data, Report-Preprint No. 22, 2005.					|
													|
													|
</pre>													|
\param      DOUBLE    c			(i)   parameter for ground substance				|
\param 		DOUBLE    k1		(i)   parameter for fiber potential				|
\param      DOUBLE	  k2		(i)   parameter for fiber potential				|
\param      DOUBLE	  gamma		(i)   penalty parameter						|
\param      DOUBLE	  epsilon	(i)   penalty parameter						|
\param		DOUBLE	  *disd		(i)	  displacement derivatives				|
\param      DOUBLE    *stress	(o)   ele stress (-resultant) vector					|
\param     	DOUBLE	  **d  		(o)   constitutive matrix					|
													|
*/
/*------------------------------------------------------------------------------------------------------*/

void c1_mat_hyper_polyconvex (	DOUBLE c,
				DOUBLE k1,
				DOUBLE k2,
				DOUBLE gamma,
				DOUBLE epsilon,
				DOUBLE *disd,
				DOUBLE *stress,
				DOUBLE **d) {


	INT i,j,k,l;/*---------------------------------------------------------------------------Counter Variables*/
	DOUBLE W,W1,W2,W3,J;/*--------------------------------------------------------------Strain-Energy function*/
	DOUBLE drittel;/*-----------------------------------------------------------------------Auxiliary Variable*/
	DOUBLE kappa;/*----------------------------------------------------------------------Dispersions Parameter*/
	DOUBLE phi,theta,K;/*-------------------------------Angles and Invariant for Anisotropic Fiber Orientation*/
	
	static DOUBLE *Inv;
	static ARRAY Inv_a;
	static DOUBLE *deltags;
	static ARRAY deltags_a;
	static DOUBLE *deltapen;
	static ARRAY deltapen_a;
	static DOUBLE *deltafib;
	static ARRAY deltafib_a;
	static DOUBLE *sum;
	static ARRAY sum_a;
	static DOUBLE *ad;
	static ARRAY ad_a;
	
	
	static DOUBLE **FT;
	static ARRAY FT_a;
	static DOUBLE **C;
	static ARRAY C_a;
	static DOUBLE **Cinv;
	static ARRAY Cinv_a;
	static DOUBLE **HC;
	static ARRAY HC_a;
	static DOUBLE **S;
	static ARRAY S_a;
	static DOUBLE **S1;
	static ARRAY S1_a;
	static DOUBLE **S2;
	static ARRAY S2_a;
	static DOUBLE **S3;
	static ARRAY S3_a;
	static DOUBLE **Sigma;
	static ARRAY Sigma_a;
	static DOUBLE **Sigma1;
	static ARRAY Sigma1_a;
	static DOUBLE **Sigma2;
	static ARRAY Sigma2_a;
	static DOUBLE **I;
	static ARRAY I_a;
	static DOUBLE **Celast;
	static ARRAY Celast_a;
	static DOUBLE **I9;
	static ARRAY I9_a;
	static DOUBLE **IC;
	static ARRAY IC_a;
	static DOUBLE **CI;
	static ARRAY CI_a;
	static DOUBLE **ICinv;
	static ARRAY ICinv_a;
	static DOUBLE **CinvI;
	static ARRAY CinvI_a;
	static DOUBLE **CC;
	static ARRAY CC_a;
	static DOUBLE **CCinv;
	static ARRAY CCinv_a;
	static DOUBLE **CinvC;
	static ARRAY CinvC_a;
	static DOUBLE **CinvCinv;
	static ARRAY CinvCinv_a;
	static DOUBLE **CinvoCinv;
	static ARRAY CinvoCinv_a;
	static DOUBLE **II;
	static ARRAY II_a;
	static DOUBLE **H_H;
	static ARRAY H_H_a;
	static DOUBLE **H_C;
	static ARRAY H_C_a;
	static DOUBLE **C_H;
	static ARRAY C_H_a;
	static DOUBLE **H_Cinv;
	static ARRAY H_Cinv_a;
	static DOUBLE **Cinv_H;
	static ARRAY Cinv_H_a;
	static DOUBLE **M;
	static ARRAY M_a;
	static DOUBLE **H;
	static ARRAY H_a;
	
	
	if (FT==NULL)
	{
		Inv = amdef("Inv",&Inv_a,3,1,"DV");
		deltags = amdef("deltags",&deltags_a,8,1,"DV");
		deltapen = amdef("deltapen",&deltapen_a,8,1,"DV");
		deltafib = amdef("deltafib",&deltafib_a,8,1,"DV");
		sum = amdef("sum",&sum_a,8,1,"DV");
		ad = amdef("ad",&ad_a,3,1,"DV");
				
		FT = amdef("FT",&FT_a,3,3,"DA");
		C = amdef("C",&C_a,3,3,"DA");
		Cinv = amdef("Cinv",&Cinv_a,3,3,"DA");
		HC = amdef("HC",&HC_a,3,3,"DA");
		S = amdef("S",&S_a,3,3,"DA");
		S1 = amdef("S1",&S1_a,3,3,"DA");
		S2 = amdef("S2",&S2_a,3,3,"DA");
		S3 = amdef("S3",&S3_a,3,3,"DA");
		Sigma = amdef("Sigma",&Sigma_a,3,3,"DA");
		Sigma1 = amdef("Sigma1",&Sigma1_a,3,3,"DA");
		Sigma2 = amdef("Sigma2",&Sigma2_a,3,3,"DA");
		I = amdef("I",&I_a,3,3,"DA");
		Celast = amdef("Celast",&Celast_a,9,9,"DA");
		I9 = amdef("I9",&I9_a,9,9,"DA");
		IC = amdef("IC",&IC_a,9,9,"DA");
		CI = amdef("CI",&CI_a,9,9,"DA");
		ICinv = amdef("ICinv",&ICinv_a,9,9,"DA");
		CinvI = amdef("CinvI",&CinvI_a,9,9,"DA");
		CC = amdef("CC",&CC_a,9,9,"DA");
		CCinv = amdef("CCinv",&CCinv_a,9,9,"DA");
		CinvC = amdef("CinvC",&CinvC_a,9,9,"DA");
		CinvCinv = amdef("CinvCinv",&CinvCinv_a,9,9,"DA");
		CinvoCinv = amdef("CinvoCinv",&CinvoCinv_a,9,9,"DA");
		II = amdef("II",&II_a,9,9,"DA");
		H_H = amdef("H_H",&H_H_a,9,9,"DA");
		H_C = amdef("H_C",&H_C_a,9,9,"DA");
		C_H = amdef("C_H",&C_H_a,9,9,"DA");
		H_Cinv = amdef("H_Cinv",&H_Cinv_a,9,9,"DA");
		Cinv_H = amdef("Cinv_H",&Cinv_H_a,9,9,"DA");
		M = amdef("M",&M_a,3,3,"DA");
		H = amdef("H",&H_a,3,3,"DA");
	}
                
    amzero(&Inv_a);
    amzero(&deltags_a);
    amzero(&deltapen_a);
    amzero(&deltafib_a);
    amzero(&sum_a);
    amzero(&ad_a);
    
    amzero(&FT_a);
    amzero(&C_a);
    amzero(&Cinv_a);
    amzero(&HC_a);
    amzero(&S_a);
    amzero(&S1_a);
    amzero(&S2_a);
    amzero(&S3_a);
    amzero(&Sigma_a);
    amzero(&Sigma1_a);
    amzero(&Sigma2_a);
    amzero(&I_a);
    amzero(&Celast_a);
    amzero(&I9_a);
    amzero(&IC_a);
    amzero(&CI_a);
    amzero(&ICinv_a);
    amzero(&CinvI_a);
    amzero(&CC_a);
    amzero(&CCinv_a);
    amzero(&CinvC_a);
    amzero(&CinvCinv_a);
    amzero(&CinvoCinv_a);
    amzero(&II_a);
    amzero(&H_H_a);
    amzero(&H_C_a);
    amzero(&C_H_a);
    amzero(&H_Cinv_a);
    amzero(&Cinv_H_a);
    amzero(&M_a);
    amzero(&H_a);

	#ifdef DEBUG
	dstrc_enter("c1_mat_hyper_polyconvex");
	#endif

	drittel=1.0/3.0;
	/*--------------------------------------------------------------------------Define kappa*/
	kappa=drittel;
	/*---------------------------------------------------------------------Fiber Orientation*/
	phi=0.0;
	theta=0.0;
	/*----------------------------------------------------------Vector of Preferred Direction*/
	ad[0]=1;
	ad[1]=0;
	ad[2]=0;
	/*---------------------------------------------------------------------Orientation Tensor*/
	for (i=0; i<3; i++) {
		for (j=0; j<3; j++) {
			M[i][j]=ad[i]*ad[j]; } }
	/*---------------------------------------------------------------------------Unity Matrix*/
	for (i=0; i<3; i++) {
		for (j=0; j<3; j++) {
			if (i==j) {
				I[i][j]=1;}
			else {
				I[i][j]=0;} } }	
	/*----------------------------------------------------------------------Structural Tensor*/
	for (i=0; i<3; i++) {
		for (j=0; j<3; j++) {
			/*H [i][j] = kappa*I[i][j] + (1-(3*kappa))*M[i][j]; } }	*/
			H [i][j] = kappa*I[i][j]; } }
	/*--------------------------------------------------------Deformation Gradient, transposed*/
	
	FT[0][0]=disd[0]+1.0;
	FT[1][1]=disd[1]+1.0;
	FT[2][2]=disd[2]+1.0;
	FT[0][1]=disd[3];/*---------------------disd[4];*/
	FT[1][0]=disd[4];/*---------------------disd[3];*/
	FT[1][2]=disd[5];/*---------------------disd[6];*/
	FT[2][1]=disd[6];/*---------------------disd[5];*/
	FT[0][2]=disd[7];/*---------------------disd[8];*/
	FT[2][0]=disd[8];/*---------------------disd[7];*/
	
	/*---------------------------------------------------------------------Jacobian Determinant*/
	J=(FT[0][0]*FT[1][1]*FT[2][2]) + (FT[0][1]*FT[1][2]*FT[2][0]) + (FT[0][2]*FT[1][0]*FT[2][1]) - (FT[2][0]*FT[1][1]*FT[0][2]) - (FT[2][1]*FT[1][2]*FT[0][0]) - (FT[2][2]*FT[1][0]*FT[0][1]);
	
	/*---------------------------------------------------------------Right Cauchy Green Tensor*/
	for (k=0; k<3; k++) {
		for (i=0; i<3; i++) {	
			for (j=0; j<3; j++) {
				C[i][k]+=(FT[i][j]*FT[k][j]); } } }	
	/*---------------------------------------------------------------------------Invariants of C*/
	c1_calc_invariants(C, Inv);
	/*---------------------------------------------------------Inverse Right Cauchy Green Tensor*/
	c1_calc_inverse (C,Cinv,Inv);
	/*----------------------------------------------------------------------------------------HC*/	
	for (i=0; i<3; i++) {
		for (j=0; j<3; j++) {
			for (k=0; k<3; k++) {
				HC[i][j]+=H[i][k]*C[k][j];} } }
	/*--------------------------------------------------------------------Anisotropic Invariant K*/
	K = HC[0][0]+HC[1][1]+HC[2][2];
	/*---------------------------------------------------------------------Strain-Energy Function*/
	/*-----------------------------------------------------------------------Ground Substance SEF*/
	/*W1=c*(Inv[0]-3.0);}*/
	W1=(c*((Inv[0]/pow(Inv[2],drittel))-3.0));	
	/*----------------------------------------------------------------------------------Fiber SEF*/
	if (K<1)
		W2=0.0;
	else if (K>=1) 
		W2=(k1/(2.0*k2))*(exp(k2*pow((K-1.0),2)-1.0));
	else { exit(1); }
	/*----------------------------------------------------------------------------Penalty Function*/
	W3 = (epsilon*(pow(Inv[2],gamma)+pow(Inv[2],(-gamma))-2));
	/*-----------------------------------------------------------------------------------------SEF*/
	W= W1 + W2 + W3;	
	/*-------------------------------------------------------------------------------2nd PK Stress*/
	for (i=0; i<3; i++)	{
		for (j=0; j<3; j++)	{
	/*------------------------------------------------------------------------------------------GS*/
			S1[i][j]=2.0*c*(pow(Inv[2],(-drittel)))*(I[i][j]-drittel*Inv[0]*Cinv[j][i]);
	/*---------------------------------------------------------------------------------------Fiber*/
			if (K<1)
				S2[i][j]=0.0;
			else if (K>=1) 
				S2[i][j] = 2.0*k1*exp(k2*pow((K-1.0),2))*(K-1.0)*H[i][j];
			else {exit(1);}
	/*-------------------------------------------------------------------------------------Penalty*/
			S3[i][j]=2*epsilon*gamma*Cinv[j][i]*(pow(Inv[2],gamma)-pow(Inv[2],(-gamma)));
	/*------------------------------------------------------------------------2nd PK Stress Tensor*/
			S[i][j]= S1[i][j] + S2[i][j] + S3[i][j];
		} }
	/*-----------------------------------------------------------------------Element Stress Vector*/
	stress[0]=S[0][0];
	stress[1]=S[1][1];
	stress[2]=S[2][2];
	stress[3]=S[0][1];
	stress[4]=S[1][2];
	stress[5]=S[0][2];
	/*-------------------------------------------------------------------------Cauchy Stress Tensor*/
	for (i=0; i<3; i++)	{
		for (k=0; k<3; k++)	{
			for (j=0; j<3; j++)	{
				Sigma1[i][k]+=S[i][j]*FT[j][k]; } } }
					
	for (i=0; i<3; i++)	{
		for (k=0; k<3; k++)	{
			for (j=0; j<3; j++)	{
				Sigma2[i][k]+=FT[j][i]*Sigma1[j][k]; } } }
		
	for (i=0; i<3; i++)	{
		for (j=0; j<3; j++)	{
			Sigma[i][j]=(1/J)*Sigma2[i][j];} }
	/*---------------------------------------------------------------------Tangent Elasticity Tensor*/
	/*---------------------------------------------------------------------------------------Delta's*/
	deltags[0] = 0.0;
	deltags[1] = 0.0;
	deltags[2] = -c*4.0*drittel*pow(Inv[2],-drittel);
	deltags[3] = 0.0;
	deltags[4] = 0.0;
	deltags[5] = ((4.0/9.0)*c*Inv[0]*pow(Inv[2],-drittel));
	deltags[6] = (4.0*drittel*c*Inv[0]*pow(Inv[2],-drittel));
	deltags[7] = 0.0;
	
	if (K<1)
		for (i=0;i<8;i++) { 
			deltafib[i]=0.0;	}		
	else if (K>=1) {
			deltafib[0] = 4.0*k1*exp(k2*pow((K-1.0),2))*((2.0*k2*pow((K-1.0),2))+1.0);
			deltafib[1] = 0.0;
			deltafib[2] = 0.0;
			deltafib[3] = 0.0;
			deltafib[4] = 0.0;
			deltafib[5] = 0.0;
			deltafib[6] = 0.0;
			deltafib[7] = 0.0; }
	else {exit(1);}
	
	deltapen[0]=0.0;
	deltapen[1]=0.0;
	deltapen[2]=0.0;
	deltapen[3]=0.0;
	deltapen[4]=0.0;
	deltapen[5]=4*epsilon*pow(gamma,2)*(pow(Inv[2],gamma)+pow(Inv[2],-gamma));
	deltapen[6]=-4*epsilon*gamma*(pow(Inv[2],gamma)-pow(Inv[2],-gamma));
	deltapen[7]=0.0;
	/*--------------------------------------------------------------------------------Tensor Products*/
	c1_calc_tensorproduct(I,I,I9);
	c1_calc_tensorproduct(I,C,IC);
	c1_calc_tensorproduct(C,I,CI);
	c1_calc_tensorproduct(I,Cinv,ICinv);
	c1_calc_tensorproduct(Cinv,I,CinvI);
	c1_calc_tensorproduct(C,C,CC);
	c1_calc_tensorproduct(C,Cinv,CCinv);
	c1_calc_tensorproduct(Cinv,C,CinvC);
	c1_calc_tensorproduct(Cinv,Cinv,CinvCinv);
	
	for (k=0; k<9; k+=3) {
		for (l=0; l<9; l+=3) {
			for (i=0; i<3; i++) {
				for (j=0; j<3; j++) {
					CinvoCinv[i+k][j+l]= 0.5*+(Cinv[k/3][i]*Cinv[l/3][j]+Cinv[k/3][j]*Cinv[l/3][i]);}}}}

	for (k=0; k<9; k+=3) {
		for (l=0; l<9; l+=3) {
			for (i=0; i<3; i++) {
				for (j=0; j<3; j++) {
					if (i==k/3 && j==l/3){
						II[i+k][j+l]=1;}
					else			{
						II[i+k][j+l]=0;} } } }	}

	c1_calc_tensorproduct(H,H,H_H);
	c1_calc_tensorproduct(H,C,H_C);
	c1_calc_tensorproduct(C,H,C_H);
	c1_calc_tensorproduct(H,Cinv,H_Cinv);
	c1_calc_tensorproduct(Cinv,H,Cinv_H);	
	
	for (k=0; k<9; k+=3) {
	  for (l=0; l<9; l+=3) {
		for (i=0; i<3; i++) {
		    for (j=0; j<3; j++) {
			sum[0]=(deltags[0]+deltapen[0])*I9[i+k][j+l]+deltafib[0]*H_H[i+k][j+l];
			sum[1]=(deltags[1]+deltapen[1])*(IC[i+k][j+l]+CI[i+k][j+l])+deltafib[1]*(H_C[i+k][j+l]+C_H[i+k][j+l]);
			sum[2]=(deltags[2]+deltapen[2])*(ICinv[i+k][j+l]+CinvI[i+k][j+l])+deltafib[2]*(H_Cinv[i+k][j+l]+Cinv_H[i+k][j+l]);
			sum[3]=(deltags[3]+deltapen[3]+deltafib[3])*CC[i+k][j+l];
			sum[4]=(deltags[4]+deltapen[4]+deltafib[4])*(CCinv[i+k][j+l]+CinvC[i+k][j+l]);
			sum[5]=(deltags[5]+deltapen[5]+deltafib[5])*CinvCinv[i+k][j+l];
			sum[6]=(deltags[6]+deltapen[6]+deltafib[6])*CinvoCinv[i+k][j+l];
			sum[7]=(deltags[7]+deltapen[7]+deltafib[7]);
			Celast[i+k][j+l] = sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7];
			} } } }	
	/*--------------------------------------------------------------------------------Constitutive Matrix*/
			d[0][0]=Celast[0][0];
			d[0][1]=Celast[1][1];
			d[0][2]=Celast[2][2];
			d[0][3]=Celast[1][0];
			d[0][4]=Celast[2][1];
			d[0][5]=Celast[2][0];
			
			d[1][0]=Celast[3][3];
			d[1][1]=Celast[4][4];
			d[1][2]=Celast[5][5];
			d[1][3]=Celast[4][3];
			d[1][4]=Celast[5][4];
			d[1][5]=Celast[5][3];
			
			d[2][0]=Celast[6][6];
			d[2][1]=Celast[7][7];
			d[2][2]=Celast[8][8];
			d[2][3]=Celast[7][6];
			d[2][4]=Celast[8][7];
			d[2][5]=Celast[8][6];
		
			d[3][0]=Celast[3][0];
			d[3][1]=Celast[4][1];
			d[3][2]=Celast[5][2];
			d[3][3]=Celast[4][0];
			d[3][4]=Celast[5][1];
			d[3][5]=Celast[5][0];
			
			d[4][0]=Celast[6][3];
			d[4][1]=Celast[7][4];
			d[4][2]=Celast[8][5];
			d[4][3]=Celast[7][3];
			d[4][4]=Celast[8][4];
			d[4][5]=Celast[8][3];
			
			d[5][0]=Celast[6][0];
			d[5][1]=Celast[7][1];
			d[5][2]=Celast[8][2];
			d[5][3]=Celast[7][0];
			d[5][4]=Celast[8][1];
			d[5][5]=Celast[8][0];
	/*-----------------------------------------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
/*exit(0);*/
return;
} /*end of c1_mat_hyper_polyconvex*/

/*---------------------------------------------------------------------------------------------------------*/
/*																											|
\brief calculation of invariants of a 3x3 tensor															|
																											|
<pre>                                                              br 09/06									|
																											|
*/
/*---------------------------------------------------------------------------------------------------------*/

void c1_calc_invariants (DOUBLE **M, DOUBLE *Inv)
{
	INT i,j;
	DOUBLE m1, m2;
	
	#ifdef DEBUG
	dstrc_enter("c1_calc_invariants");
	#endif
	/*-----------------------------------------------------------------------------------------------Inv1*/
	Inv[0] = M[0][0]+M[1][1]+M[2][2];
	/*-----------------------------------------------------------------------------------------------Inv2*/
	m1=0.0;
	for (i=0; i<3; i++) {
		for (j=0; j<3; j++) {
			m1+=(M[i][i] * M[j][j]);}}
	m2=0.0;		
	for (i=0; i<3; i++) {
		for (j=0; j<3; j++) {
			 	m2+=M[j][i]*M[i][j]; }}
			
	Inv[1]=0.5*(m1-m2);	
	/*------------------------------------------------------------------------------------------------Inv3*/
	Inv[2]=(M[0][0]*M[1][1]*M[2][2]) + (M[0][1]*M[1][2]*M[2][0]) + (M[0][2]*M[1][0]*M[2][1]) - (M[2][0]*M[1][1]*M[0][2]) - (M[2][1]*M[1][2]*M[0][0]) - (M[2][2]*M[1][0]*M[0][1]);
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /*end of c1_calc_invariants*/

/*---------------------------------------------------------------------------------------------------------*/
/*																											|
\brief calculation the inverse of a 3x3 tensor																|
																											|
<pre>                                                              br 09/06									|
																											|
*/
/*---------------------------------------------------------------------------------------------------------*/

void c1_calc_inverse (DOUBLE **M, DOUBLE **Minv, DOUBLE *Inv)
{
	#ifdef DEBUG
	dstrc_enter("c1_calc_inverse");
	#endif
	
	if (Inv[2]==0) {
		printf("Matrix nicht invertierbar!\n");
	}
	else {
	Minv[0][0]= 1/Inv[2] * (M[1][1]*M[2][2] - M[2][1]*M[1][2]);
	Minv[1][0]=-1/Inv[2] * (M[0][1]*M[2][2] - M[2][1]*M[0][2]);
	Minv[2][0]= 1/Inv[2] * (M[0][1]*M[1][2] - M[1][1]*M[0][2]);
	Minv[0][1]=-1/Inv[2] * (M[1][0]*M[2][2] - M[2][0]*M[1][2]);
	Minv[1][1]= 1/Inv[2] * (M[0][0]*M[2][2] - M[2][0]*M[0][2]);
	Minv[2][1]=-1/Inv[2] * (M[0][0]*M[1][2] - M[1][0]*M[0][2]);
	Minv[0][2]= 1/Inv[2] * (M[1][0]*M[2][1] - M[2][0]*M[1][1]);
	Minv[1][2]=-1/Inv[2] * (M[0][0]*M[2][1] - M[2][0]*M[0][1]);
	Minv[2][2]= 1/Inv[2] * (M[0][0]*M[1][1] - M[1][0]*M[0][1]);
	}
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /*end of c1_calc_inverse*/

/*---------------------------------------------------------------------------------------------------------*/
/*																											|
\brief calculation of the tensor product of two 3x3 tensors													|
																											|
<pre>                                                              br 09/06									|
																											|
*/
/*---------------------------------------------------------------------------------------------------------*/

void c1_calc_tensorproduct (DOUBLE **A, DOUBLE **B, DOUBLE **AB)
{
	#ifdef DEBUG
	dstrc_enter("c1_calc_tensorproduct");
	#endif
	INT i,j,k,l;
	
	for (k=0; k<9; k+=3) {
		for (l=0; l<9; l+=3) {
			for (i=0; i<3; i++) {
				for (j=0; j<3; j++) {
					AB[i+k][j+l]= A[k/3][l/3]*B[i][j]; } } } }
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /*end of: c1_calc_tensorproduct*/

#endif
