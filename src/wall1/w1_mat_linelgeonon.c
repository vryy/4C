#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"
/*----------------------------------------------------------------------*
 | constitutive matrix - linear elastic - 2D              ah   06/02    |
 | plane stress, plane strain,                                          |
 *----------------------------------------------------------------------*/
void w1_mat_linelgeonon(double ym, 
                        double pv, 
                        WALL_TYPE wtype,
                        double *strain, 
                        double **d,
                        double **stress,
                        int numeps)
{
double e1, e2, e3, a1, b1, c1;
double svector[3];
int i,k;
#ifdef DEBUG 
dstrc_enter("w1_mat_linelgeonon");
#endif
/*---------------------------------material-tangente-- plane stress ---*/
  switch(wtype)
  {
  case plane_stress:
    e1=ym/(1. - pv*pv);
    e2=pv*e1;
    e3=e1*(1. - pv)/2.;

    d[0][0]=e1;
    d[0][1]=e2;
    d[0][2]=0.;
    d[0][3]=0.;
    
    d[1][0]=e2;
    d[1][1]=e1;
    d[1][2]=0.;
    d[1][3]=0.;    
    
    d[2][0]=0.;
    d[2][1]=0.;
    d[2][2]=e3;
    d[2][3]=e3;
    
    d[3][0]=0.;
    d[3][1]=0.;
    d[3][2]=e3;
    d[3][3]=e3;
  break;
  default:
/*----------- material-tangente - plane strain, rotational symmetry ---*/
    c1=ym/(1.0+pv);
    b1=c1*pv/(1.0-2.0*pv);
    a1=b1+c1;

    d[0][0]=a1;
    d[0][1]=b1;
    d[0][2]=0.;
    d[0][3]=0.;
    
    d[1][0]=b1;
    d[1][1]=a1;
    d[1][2]=0.;
    d[1][3]=0.;

    d[2][0]=0.;
    d[2][1]=0.;
    d[2][2]=c1/2.;
    d[2][3]=c1/2.;
    
    d[3][0]=0.;
    d[3][1]=0.;
    d[3][2]=c1/2;
    d[3][3]=c1/2;
  break;
  }
/*-------------------------- evaluate 2.PK-stresses -------------------*/ 
/*------------------ Summenschleife -> += (2.PK stored as vecor) ------*/ 
for (k=0; k<3; k++)
{ 
  svector[k]=0;
  for (i=0; i<numeps; i++)
  { 
     svector[k] += d[k][i] * strain[i];                     
  }
}
/*------------------ 2.PK stored as matrix -----------------------------*/
stress[0][0]=svector[0]; 
stress[0][2]=svector[2]; 
stress[1][1]=svector[1]; 
stress[1][3]=svector[2]; 
stress[2][0]=svector[2]; 
stress[2][2]=svector[1]; 
stress[3][1]=svector[2]; 
stress[3][3]=svector[0]; 
                  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of w1_mat_linelgeonon */

