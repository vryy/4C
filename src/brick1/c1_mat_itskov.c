

/*----------------------------------------------------------------------*/
#ifdef D_BRICK1

#include "../headers/standardtypes.h"
#include "brick1.h"
#include "brick1_prototypes.h"
void c1_mat_itskov( COMPOGDEN *mat, DOUBLE *disd, DOUBLE *stress, DOUBLE **d, INT ip);

extern ALLDYNA      *alldyn;

void c1_mat_itskov(
		 COMPOGDEN *mat,
		 DOUBLE *disd,
		 DOUBLE *stress,
		 DOUBLE **d,
		 INT ip)
{
	INT i,j,k,l;
	DOUBLE alpha, beta, gamma, m;    /*Parameter der strain-energy-function*/
	DOUBLE q, epsilonPen, gammaPen;  /*Parameter der Straffunktionen*/
	DOUBLE energy, energyconstr, Kr, K;
	static DOUBLE **FT;
	static ARRAY FT_a;
	static DOUBLE **C;
	static ARRAY C_a;
	static DOUBLE *Inv;
	static ARRAY Inv_a;
	static DOUBLE *J;
	static ARRAY J_a;
	static DOUBLE **I;
	static ARRAY I_a;
	static DOUBLE **SPK;
	static ARRAY SPK_a;
	static DOUBLE **Cinv;
	static ARRAY Cinv_a;
	static DOUBLE **SPKconstr;
	static ARRAY SPKconstr_a;
	static DOUBLE *delta;
	static ARRAY delta_a;
	static DOUBLE **IxI;
	static ARRAY IxI_a;
	static DOUBLE **LxL;
	static ARRAY LxL_a;
	static DOUBLE **CinvxCinv;
	static ARRAY CinvxCinv_a;
	static DOUBLE **CinvLCinvxCinvLCinv;
	static ARRAY CinvLCinvxCinvLCinv_a;
	static DOUBLE **CinvoCinv;
	static ARRAY CinvoCinv_a;
	static DOUBLE **CinvoCinvLCinv;
	static ARRAY CinvoCinvLCinv_a;
	static DOUBLE **CinvLCinvoCinv;
	static ARRAY CinvLCinvoCinv_a;	
	static DOUBLE **II;
	static ARRAY II_a;
	static DOUBLE **Celast;
	static ARRAY Celast_a;
	static DOUBLE **Celastconstraint;
	static ARRAY Celastconstraint_a;
	static DOUBLE **CinvL;
	static ARRAY CinvL_a;
	static DOUBLE **L;
	static ARRAY L_a;
	static DOUBLE **CinvLCinv;
	static ARRAY CinvLCinv_a;
	static DOUBLE **CL;
	static ARRAY CL_a;
	
	/*nur fuer Ableitung nach Holzapfel*/
	/*static DOUBLE **IxC;
	static ARRAY IxC_a;
	static DOUBLE **CxI;
	static ARRAY CxI_a;
	static DOUBLE **IxCinv;
	static ARRAY IxCinv_a;
	static DOUBLE **CinvxI;
	static ARRAY CinvxI_a;
	static DOUBLE **CxC;
	static ARRAY CxC_a;
	static DOUBLE **CxCinv;
	static ARRAY CxCinv_a;
	static DOUBLE **CinvxC;
	static ARRAY CinvxC_a;*/


	if (FT==NULL)
	{
		FT = amdef("FT",&FT_a,3,3,"DA");
		C = amdef("C",&C_a,3,3,"DA");
		Inv = amdef("Inv",&Inv_a,3,1,"DV");
		J = amdef("J",&J_a,3,1,"DV");
		I = amdef("I",&I_a,3,3,"DA");
		SPK = amdef("SPK",&SPK_a,3,3,"DA");
		Cinv = amdef("Cinv",&Cinv_a,3,3,"DA");
		SPKconstr = amdef("SPKconstr",&SPKconstr_a,3,3,"DA");
		delta = amdef("delta",&delta_a,3,1,"DV");
		IxI = amdef("IxI",&IxI_a,9,9,"DA");
		LxL = amdef("LxL",&LxL_a,9,9,"DA");
		CinvxCinv = amdef("CinvxCinv",&CinvxCinv_a,9,9,"DA");
		CinvLCinvxCinvLCinv = amdef("CinvLCinvxCinvLCinv",&CinvLCinvxCinvLCinv_a,9,9,"DA");
		CinvoCinv = amdef("CinvoCinv",&CinvoCinv_a,9,9,"DA");
		CinvoCinvLCinv = amdef("CinvoCinvLCinv",&CinvoCinvLCinv_a,9,9,"DA");
		CinvLCinvoCinv = amdef("CinvLCinvoCinv",&CinvLCinvoCinv_a,9,9,"DA");
		II = amdef("II",&II_a,9,9,"DA");
		Celast = amdef("Celast",&Celast_a,9,9,"DA");
		Celastconstraint = amdef("Celastconstraint",&Celastconstraint_a,9,9,"DA");
		CinvL = amdef("CinvL",&CinvL_a,3,3,"DA");
		L = amdef("L",&L_a,3,3,"DA");
		CinvLCinv = amdef("CinvLCinv",&CinvLCinv_a,3,3,"DA");
		CL = amdef("CL",&CL_a,3,3,"DA");
		
		/*nur fuer Ableitung nach Holzapfel*/
		/*IxC = amdef("IxC",&IxC_a,9,9,"DA");
		CxI = amdef("CxI",&CxI_a,9,9,"DA");
		IxCinv = amdef("IxCinv",&IxCinv_a,9,9,"DA");
		CinvxI = amdef("CinvxI",&CinvxI_a,9,9,"DA");
		CxC = amdef("CxC",&CxC_a,9,9,"DA");
		CxCinv = amdef("CxCinv",&CxCinv_a,9,9,"DA");
		CinvxC = amdef("CinvxC",&CinvxC_a,9,9,"DA");*/		
	}
	
	amzero(&FT_a);
	amzero(&C_a);
	amzero(&Inv_a);
	amzero(&J_a);
	amzero(&I_a);
	amzero(&SPK_a);
	amzero(&Cinv_a);
	amzero(&SPKconstr_a);
	amzero(&delta_a);
	amzero(&IxI_a);
	amzero(&LxL_a);   
    amzero(&CinvxCinv_a);
    amzero(&CinvLCinvxCinvLCinv_a);
    amzero(&CinvoCinv_a);
    amzero(&CinvoCinvLCinv_a);
    amzero(&CinvLCinvoCinv_a);
    amzero(&II_a);
    amzero(&Celast_a);
    amzero(&Celastconstraint_a);
    amzero(&CinvL_a);
    amzero(&L_a);
    amzero(&CinvLCinv_a);
    amzero(&CL_a);
    
    /*nur fuer Ableitung nach Holzapfel*/
    /*amzero(&IxC_a);
    amzero(&CxI_a);
    amzero(&IxCinv_a);
    amzero(&CinvxI_a);
    amzero(&CxC_a);
    amzero(&CxCinv_a);
    amzero(&CinvxC_a);*/
	
	
	#ifdef DEBUG
	dstrc_enter("c1_mat_itskov");
	#endif

    exit(0);

	/*Parameter der strain-energy function*/
	alpha=30.0;
	beta=25.0;
	gamma=1.0;
	m=1.3563;
	
	/*Parameter der Straffunktion nach Balzani*/
	epsilonPen=50.0; 
	gammaPen=50.0;
	
	/*Parameter der Straffunktion im Itskov-Paper*/
	/*q=15.0;*/
	
	
/*---------------------------------------------------------Definition structural tensor L (siehe Itskov-Paper Gleichung (37))*/
	for (k=0; k<3; k++) {
		for (i=0; i<3; i++) {
			if (i==k)
				L[i][k]=1.0/3.0;
			else
				L[i][k]=0.0; } }	

/*--------------------------------------------------um den aktuellen Zeitschritt mit auszugeben*/
STRUCT_DYNAMIC *sdyn;   
INT disnum;
DOUBLE time;
disnum = 0;     
sdyn = alldyn[disnum].sdyn;
time = sdyn->time;


/*---------------------------------------------------------Definition Einheitstensor*/	
	for (k=0; k<3; k++) {
		for (i=0; i<3; i++) {
			if (i==k)
				I[i][k]=1;
			else
				I[i][k]=0; } }
										
/*---------------------------------------------------------Tensoren mit Nulleintraegen fuellen*/	
					
	for (k=0; k<3; k++) {
		for (i=0; i<3; i++) {
					CinvL[i][k]=0;} }	
					
	for (k=0; k<3; k++) {
		for (i=0; i<3; i++) {
					CinvLCinv[i][k]=0;} }	
					
	for (k=0; k<3; k++) {
		for (i=0; i<3; i++) {
					CL[i][k]=0;} }		
	
/*---------------------------------------------------------Deformation Gradient, transposed */
	FT[0][0]=disd[0]+1.0;
	FT[1][1]=disd[1]+1.0;
	FT[2][2]=disd[2]+1.0;
	FT[0][1]=disd[3];
	FT[1][0]=disd[4];
	FT[1][2]=disd[5];
	FT[2][1]=disd[6];
	FT[0][2]=disd[7];
	FT[2][0]=disd[8];
	
	/*nur zum testen*/
	/*double lam11 = 1.5;
	FT[0][0] = lam11;
	FT[1][1] = 1./sqrt(lam11);
	FT[2][2] = 1./sqrt(lam11);
	FT[1][0] = 0.0;
	FT[0][1] = 0.0;
	FT[2][1] = 0.0;
	FT[1][2] = 0.0;
	FT[2][0] = 0.0;
	FT[0][2] = 0.0;*/	
/*---------------------------------------------------------Right Cauchy Green Tensor*/
	for (k=0; k<3; k++) {
		for (i=0; i<3; i++) {
			for (j=0; j<3; j++) {
				C[i][k]+=(FT[i][j]*FT[k][j]); } } }

/*---------------------------------------------------------Bestimmung von q (damit Beitrag der Straffunktion im undeformiertren zustand 0 ist)*/
	if (C[1][1]==1 && C[2][2]==1 && C[0][0]==1)
		q=0;
/*---------------------------------------------------------Invariants I*/
	c1_calc_invariants(C, Inv); 
/*---------------------------------------------------------Inverse Right Cauchy Green Tensor*/				
	c1_calc_inverse (C,Cinv,Inv);
/*---------------------------------------------------------Cinv*L*Cinv*/
	for (i=0; i<3; i++) {
		for (k=0; k<3; k++) {
			for (j=0; j<3; j++) {
				CinvL[i][k]+=Cinv[i][j]*L[j][k];} } }	
					
	for (i=0; i<3; i++) {
		for (k=0; k<3; k++) {
			for (j=0; j<3; j++) {
				CinvLCinv[i][k]+=CinvL[i][j]*Cinv[j][k];} } }	

/*---------------------------------------------------------Invarianten J mit structural tensor*/
/*---------------------------------------------------------1. Invariante: tr(CL)*/
   	for (i=0; i<3; i++) {
		for (k=0; k<3; k++) {
			for (j=0; j<3; j++) {
					CL[i][k]+=C[i][j]*L[j][k];} } }	
   	   	
   	J[0]= CL[0][0]+CL[1][1]+CL[2][2];	
/*---------------------------------------------------------3. Invariante: detC*/
	J[2]=Inv[2];
/*---------------------------------------------------------2. Invariante: tr[(cofC)L]*/
	J[1]=J[2]*(CinvL[0][0]+CinvL[1][1]+CinvL[2][2]);
/*---------------------------------------------------------2. Invariante für inkompressibel: tr[Cinv*L]*/	
	Kr=(CinvL[0][0]+CinvL[1][1]+CinvL[2][2]); 
/*------------------------------------------------------------------------------------------*/

	
/*---------------------------------------------------------constraint function*/
    /*Straffunktion nach Balzani*/
    energyconstr= (epsilonPen*(pow(Inv[2],gammaPen)+pow(Inv[2],(-gammaPen))-2));
       
/*---------------------------------------------------------strain-energy-function*/	

	/*kompressibel*/
	/*energy=m/4*((exp(alpha*(J[0]-1.0))-1.0)/alpha+(exp(beta*(J[1]-1.0))-1.0)/beta+(pow(J[2],-gamma)-1.0)/gamma)+energyconstr;*/
	
	/*inkompressibel (also I3=1 eingesetzt)*/
	energy=m/4.0*((exp(alpha*(J[0]-1.0))-1.0)/alpha+(exp(beta*(Kr-1.0))-1.0)/beta)+energyconstr;
	
	
/*--------------------------------------------------------------------------------*/
/*--------------------------2nd Piola-Kirchhoff Stress ---------------------------*/
/*--------------------------------------------------------------------------------*/	
/*---------------------------------------------------------constraint function*/
for (i=0; i<3; i++) {
	for (k=0; k<3; k++) {
	
	/*Straffunktion nach Itskov*/			
	/*SPKconstr[i][k]=q/3.0*Cinv[i][k];*/
	
	/*Straffunktion nach Balzani*/
	SPKconstr[i][k]=2.0*epsilonPen*gammaPen*Cinv[k][i]*(pow(Inv[2],gammaPen)-pow(Inv[2],(-gammaPen)));
	
	 }}
	
/*---------------------------------------------------------2nd Piola-Kirchhoff Stress tensor*/
for (k=0; k<3; k++) {
		for (i=0; i<3; i++) {
			
						
			/* --> Siehe Itskov Paper --- Inkompressibel*/
			SPK[i][k]=m/2.0*(exp(alpha*(J[0]-1.0))*L[i][k]
						-exp(beta*(Kr-1.0))*CinvLCinv[i][k])
						+ SPKconstr[i][k];
			
			/* --> Ableitung Siehe Holzapfel-Buch */		
			/*SPK[i][k]=m/2.0*(exp(alpha*(Inv[0]/3.0-1.0))/3.0*I[i][k]
						+ exp(beta*(Inv[1]/3.0-1.0))/3.0*(Inv[0]*I[i][k]-C[i][k])
						- pow(Inv[2],-gamma)*Cinv[i][k])
						+ SPKconstr[i][k];*/						
			} }			
						

/*---------------------------------------------------------Element Stress Vector*/
	stress[0]=SPK[0][0];
	stress[1]=SPK[1][1];
	stress[2]=SPK[2][2];
	stress[3]=SPK[0][1];
	stress[4]=SPK[1][2];
	stress[5]=SPK[0][2];
	
	
/*----------------------------------------------------------------------------------------------*/	
/*-------------------------------------tangent elasticity tensor--------------------------------*/
/*----------------------------------------------------------------------------------------------*/

/*---------------------------------------------------------Deltas für Ableitung nach Holzapfel-Buch*/	
	delta[0] = m * ( alpha/9. * exp(alpha*(Inv[0]/3.-1.))
	               + 1./3. * exp(beta*(Inv[1]/3.-1.))
	               + pow(Inv[0],2.) * beta/9. * exp(beta*(Inv[1]/3.-1.)));
	delta[1] = -m*Inv[0]*beta/9.*exp(beta*(Inv[1]/3.-1.));
	delta[2] = 0.0;
	delta[3] = m*beta/9.*exp(beta*(Inv[1]/3.-1.));
	delta[4] = 0.0;
	delta[5] = m*gamma*pow(Inv[2],-gamma);
	delta[6] = m*pow(Inv[2],-gamma);
	delta[7] = -m/3.*exp(beta*(Inv[1]/3.-1.));
/*---------------------------------------------------------Tensor Products*/
	
	c1_calc_tensorproduct(L,L,LxL);
	c1_calc_tensorproduct(Cinv,Cinv,CinvxCinv);
	c1_calc_tensorproduct(CinvLCinv,CinvLCinv,CinvLCinvxCinvLCinv);
	
	/*nur fuer Ableitung nach Holzapfel-Buch*/
	/* c1_calc_tensorproduct(I,I,IxI);
	c1_calc_tensorproduct(I,C,IxC);
	c1_calc_tensorproduct(C,I,CxI);
	c1_calc_tensorproduct(I,Cinv,IxCinv);
	c1_calc_tensorproduct(Cinv,I,CinvxI);
	c1_calc_tensorproduct(C,C,CxC);
	c1_calc_tensorproduct(C,Cinv,CxCinv);
	c1_calc_tensorproduct(Cinv,C,CinvxC);*/
	
			
	for (k=0; k<9; k+=3) {
		for (l=0; l<9; l+=3) {
			for (i=0; i<3; i++) {
				for (j=0; j<3; j++) {
					CinvoCinv[i+k][j+l]= 0.5*+(Cinv[k/3][i]*Cinv[l/3][j]+Cinv[k/3][j]*Cinv[l/3][i]);}}}}

	for (k=0; k<9; k+=3) {
		for (l=0; l<9; l+=3) {
			for (i=0; i<3; i++) {
				for (j=0; j<3; j++) {
					CinvoCinvLCinv[i+k][j+l]= 0.5*+(Cinv[k/3][i]*CinvLCinv[l/3][j]+Cinv[k/3][j]*CinvLCinv[l/3][i]);}}}}

	for (k=0; k<9; k+=3) {
		for (l=0; l<9; l+=3) {
			for (i=0; i<3; i++) {
				for (j=0; j<3; j++) {
					CinvLCinvoCinv[i+k][j+l]= 0.5*+(CinvLCinv[k/3][i]*Cinv[l/3][j]+CinvLCinv[k/3][j]*Cinv[l/3][i]);}}}}
	
					
	/*II nach Definition im Holzapfelbuch IIijkl=dik*djl*/
	/*for (k=0; k<9; k+=3) {
		for (l=0; l<9; l+=3) {
			for (i=0; i<3; i++) {
				for (j=0; j<3; j++) {
					if (i==k/3 && j==l/3){
						II[i+k][j+l]=1;}
					else			{
						II[i+k][j+l]=0;} } } }	}*/
						
	/*II nach andrer Definition (zB siehe NiLiFEM-Skript) IIijkl=0.5(dik*djl+dil*djk)*/					
	/*for (k=0; k<9; k+=3) {
		for (l=0; l<9; l+=3) {
			for (i=0; i<3; i++) {
				for (j=0; j<3; j++) {
					II[i+k][j+l]= 0.5*+(I[k/3][i]*I[l/3][j]+I[k/3][j]*I[l/3][i]);}}}}*/		
	
/*---------------------------------------------------------constraint function*/
for (k=0; k<9; k+=3) {
	  for (l=0; l<9; l+=3) {
		for (i=0; i<3; i++) {
		    for (j=0; j<3; j++) {		
			
			/*Straffunktion nach Itskov*/
			/*Celastconstraint[i+k][j+l] = 2.0/3.0*q*(1.0/3.0*CinvxCinv[i+k][j+l] - CinvoCinv[i+k][j+l]);*/
						
			/*Straffunktion nach Balzani*/
			Celastconstraint[i+k][j+l] = 4.0*epsilonPen*pow(gammaPen,2)*(pow(Inv[2],gammaPen)+pow(Inv[2],-gammaPen))*CinvxCinv[i+k][j+l]
										-4.0*epsilonPen*gammaPen*(pow(Inv[2],gammaPen)-pow(Inv[2],-gammaPen))*CinvoCinv[i+k][j+l];
				
			}}}}	
				
/*---------------------------------------------------------tangent elasticity tensor*/

	for (k=0; k<9; k+=3) {
	  for (l=0; l<9; l+=3) {
		for (i=0; i<3; i++) {
		  for (j=0; j<3; j++) {	
 						   
	/* --> Siehe Itskov Paper --- Inkompressibel */	
	Celast[i+k][j+l] = m*  (alpha*exp(alpha*(J[0]-1.0))*LxL[i+k][j+l]
						   +beta*exp(beta*(Kr-1.0))*CinvLCinvxCinvLCinv[i+k][j+l]
						   +exp(beta*(Kr-1.0))*(CinvoCinvLCinv[i+k][j+l]+CinvLCinvoCinv[i+k][j+l]))
 						   + Celastconstraint[i+k][j+l];
 					
	/* --> Siehe Holzapfel-Buch*/	
 	/*Celast[i+k][j+l] = delta[0]*IxI[i+k][j+l]
						+ delta[1]*(IxC[i+k][j+l]+CxI[i+k][j+l])
						+ delta[2]*(IxCinv[i+k][j+l]+CinvxI[i+k][j+l])
						+ delta[3]*CxC[i+k][j+l]
						+ delta[4]*(CxCinv[i+k][j+l]+CinvxC[i+k][j+l])
						+ delta[5]*CinvxCinv[i+k][j+l]
						+ delta[6]*CinvoCinv[i+k][j+l]
						+ delta[7]*II[i+k][j+l]
						+ Celastconstraint[i+k][j+l];*/
								
			}}}}
/*---------------------------------------------------------Constitutive Matrix*/
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
	
/*---------------------------------------------------------Ausgaben zum testen*/
if (ip==0)
	{
	printf("time:%f      detC:%f\n", time, Inv[2]);
	
	printf("\nRight cauchy green\n");
	for (k=0; k<3; k++) {
		for (i=0; i<3; i++) {
				printf("%f   ", C[i][k]);}
				printf("\n");}
	
	printf("\nSecond Piola Kirchhhoff\n");
	for (k=0; k<3; k++) {
		for (i=0; i<3; i++) {
				printf("%f   ", SPK[i][k]);}
				printf("\n");}
				
	/*printf("\nCelast\n");				
	for (k=0; k<9; k++) {
	  for (l=0; l<9; l++) {
			printf("%f  ", Celast[k][l]);}
			printf("\n");}*/
	
	/*K = d[0][0]*sqrt(C[0][0]) - 0.5*d[0][1]/C[0][0] - 0.5*d[0][2]/C[0][0];
	printf("\nK: %lf \n",K);*/
	
	FILE *pDatei = fopen("test_I3.txt", "at");
	fprintf(pDatei, "time:%f      detC:%f          SPK11:%f          SPKconstr11:%f\n", time, Inv[2], SPK[0][0], SPKconstr[0][0]);
	fclose(pDatei);
	
	}
	
/*exit(0);*/

#ifdef DEBUG
	dstrc_exit();
#endif
return;
}

#endif /*D_BRICK1 */
