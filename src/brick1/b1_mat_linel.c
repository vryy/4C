#include "../headers/standardtypes.h"
/*----------------------------------------------------------------------*
 | program to establish local material law for brick element    al 9/01 |
 | stress-strain law for isotropic material.                            |
 *----------------------------------------------------------------------*/
void b1_mat_linel(STVENANT *mat, double **d)
{
double d1,d2,d3;
double ym,pv;/*------------------------------------------ mat constants */
#ifdef DEBUG 
dstrc_enter("b1_mat_linel");
#endif
/*----------------------------------------------------------------------*/
ym  = mat->youngs;
pv  = mat->possionratio;
/*----------------------------------- evaluate basic material values ---*/
d1=ym*(1.0 - pv)/((1.0 + pv)*(1.0 - 2.0*pv));
d2=ym*pv/((1.0 + pv)*(1.0 - 2.0*pv));
d3=ym/((1.0 + pv)*2.0);
/*------------------------------------ set values in material-matrix ---*/
d[0][0]=d1;
d[0][1]=d2;
d[0][2]=d2;
d[0][3]=0.0;
d[0][4]=0.0;
d[0][5]=0.0;

d[1][0]=d2;
d[1][1]=d1;
d[1][2]=d2;
d[1][3]=0.0;
d[1][4]=0.0;
d[1][5]=0.0;

d[2][0]=d2;
d[2][1]=d2;
d[2][2]=d1;
d[2][3]=0.0;
d[2][4]=0.0;
d[2][5]=0.0;

d[3][0]=0.0;
d[3][1]=0.0;
d[3][2]=0.0;
d[3][3]=d3;
d[3][4]=0.0;
d[3][5]=0.0;

d[4][0]=0.0;
d[4][1]=0.0;
d[4][2]=0.0;
d[4][3]=0.0;
d[4][4]=d3;
d[4][5]=0.0;

d[5][0]=0.0;
d[5][1]=0.0;
d[5][2]=0.0;
d[5][3]=0.0;
d[5][4]=0.0;
d[5][5]=d3;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of b1_mat_linel */
