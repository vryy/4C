/*!----------------------------------------------------------------------
\file
\brief contains the routine 
 - s9_mat_orth3D: which calculates the constitutive matrix for a linear
                  elastic, otrthotropic material 

*----------------------------------------------------------------------*/
#ifdef D_SHELL9
#include "../headers/standardtypes.h"
#include "shell9.h"

/*! 
\addtogroup SHELL9 
*//*! @{ (documentation module open)*/ 

/*!----------------------------------------------------------------------
\brief calculates the constitutive matrix for a linear elastic, 
orthotropic material.

<pre>                                                             sh 02/03
This routine to establish local material law stress-strain law for linear
elastic orthotropic material (3D-formulation). 
-> sorting [11,22,33,12,23,13]
</pre>
\param   *mat  ORTHOTROPIC  (i)  Material proberties
\param  **d    DOUBLE       (o)  constitutive matrix

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: s9_call_mat()  [s9_call_mat.c]

*----------------------------------------------------------------------*/
void s9_mat_orth3D(ORTHOTROPIC *mat, DOUBLE **d)
{
DOUBLE emod1,emod2,emod3;
DOUBLE gmod12,gmod13,gmod23;
DOUBLE xnue12,xnue13,xnue23;
DOUBLE xnue21,xnue32,xnue31;
DOUBLE emod;
DOUBLE delta;
#ifdef DEBUG 
dstrc_enter("s9_mat_orth3D");
#endif
/*----------------------------------------------------------------------*/
emod1  = mat->emod1;
emod2  = mat->emod2;
emod3  = mat->emod3;
gmod12 = mat->gmod12;
gmod13 = mat->gmod13;
gmod23 = mat->gmod23;
xnue12 = mat->xnue12;
xnue13 = mat->xnue13;
xnue23 = mat->xnue23;
/*----------------------------------- 3D-MATERIALTENSOR (ORTHOTROPIC) --*/
xnue31 = xnue13 * emod3/emod1;
xnue32 = xnue23 * emod3/emod2;
xnue21 = xnue12 * emod2/emod1;
emod   = emod1 * emod2 * emod3;
delta = 1. - xnue13*xnue31 - xnue23*xnue32 - xnue12*xnue21;
delta = (delta - 2.*xnue31*xnue23*xnue12)/emod;
/*---------------------------------------------------------------------*/
d[0][0]=(1. - xnue23*xnue32)/(emod2*emod3*delta);
d[0][1]=(xnue12+xnue13*xnue32)/(emod1*emod3*delta);
d[0][2]=(xnue13+xnue12*xnue23)/(emod1*emod2*delta);
d[0][3]=0.0;
d[0][4]=0.0;
d[0][5]=0.0;

d[1][0]=d[0][1];
d[1][1]=(1. - xnue13*xnue31)/(emod1*emod3*delta);
d[1][2]=(xnue23+xnue21*xnue13)/(emod1*emod2*delta);
d[1][3]=0.0;
d[1][4]=0.0;
d[1][5]=0.0;

d[2][0]=d[0][2];
d[2][1]=d[1][2];
d[2][2]=(1. - xnue12*xnue21)/(emod1*emod2*delta);
d[2][3]=0.0;
d[2][4]=0.0;
d[2][5]=0.0;

d[3][0]=0.0;
d[3][1]=0.0;
d[3][2]=0.0;
d[3][3]=gmod12;
d[3][4]=0.0;
d[3][5]=0.0;

d[4][0]=0.0;
d[4][1]=0.0;
d[4][2]=0.0;
d[4][3]=0.0;
d[4][4]=gmod23;
d[4][5]=0.0;

d[5][0]=0.0;
d[5][1]=0.0;
d[5][2]=0.0;
d[5][3]=0.0;
d[5][4]=0.0;
d[5][5]=gmod13;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s9_mat_orth3D */

/*----------------------------------------------------------------------*/
#endif /*D_SHELL9*/
/*! @} (documentation module close)*/
