/*!----------------------------------------------------------------------
\file
\brief contains the routine 
 - s9_mat_linel3D: which calculates the constitutive matrix for a linear
                   elastic, isotropic material 

*----------------------------------------------------------------------*/
#ifdef D_SHELL9
#include "../headers/standardtypes.h"
#include "shell9.h"

/*! 
\addtogroup SHELL9 
*//*! @{ (documentation module open)*/ 

/*!----------------------------------------------------------------------
\brief establish local material; equal to 'c1_mat_linel()'

<pre>                                                              al 06/02
This routine to establish local material law stress-strain law for linear
elastic isotropic material (3D-formulation). -> sorting [11,22,33,12,23,13]

</pre>
\param       youngs   DOUBLE  (i)   young's modulus 
\param possionratio   DOUBLE  (i)   poisson's ratio 
\param          **d   DOUBLE  (i)   constitutive matrix

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: s9_call_mat()  [s9_call_mat.c]

*----------------------------------------------------------------------*/
void s9_mat_linel3D(DOUBLE youngs,
                    DOUBLE possionratio,
                    DOUBLE **d)
{
DOUBLE d1,d2,d3;
DOUBLE ym,pv;/*------------------------------------------ mat constants */
#ifdef DEBUG 
dstrc_enter("s9_mat_linel3D");
#endif
/*----------------------------------------------------------------------*/
ym = youngs;
pv = possionratio;
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
} /* end of s9_mat_linel3D */

/*----------------------------------------------------------------------*/
#endif /*D_SHELL9*/
/*! @} (documentation module close)*/
