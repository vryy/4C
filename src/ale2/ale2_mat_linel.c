/*!----------------------------------------------------------------------
\file
\brief contains the routine 'ale2_mat_linel' which calculates the
linear elastic constitutive matrix for a 2D ale element

*----------------------------------------------------------------------*/
#ifdef D_ALE
#include "../headers/standardtypes.h"
#include "ale2.h"

/*! 
\addtogroup Ale 
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief  calculates the linear elastic constitutive matrix

<pre>                                                              mn 06/02 
This routine calculates the linear elastic isotropic constitutive 
matrix for a 2D ale element.

</pre>
\param *mat   STVENANT   (i)   my material
\param *d     DOUBLE     (o)   the constitutive matrix

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: ale2_static_ke()

*----------------------------------------------------------------------*/
void ale2_mat_linel(STVENANT *mat, DOUBLE **d)
{
DOUBLE a1, b1, c1;
DOUBLE ym,pv;
#ifdef DEBUG 
dstrc_enter("ale2_mat_linel");
#endif
/*----------------------------------------------------------------------*/
ym  = mat->youngs;
pv  = mat->possionratio;
/*----------------------------------------------------- plane stress ---*/
/*  switch(wtype)
  {
  case plane_stress:
    e1=ym/(1. - pv*pv);
    e2=pv*e1;
    e3=e1*(1. - pv)/2.;

    d[0][0]=e1;
    d[0][1]=e2;
    d[0][2]=0.;
    d[1][0]=e2;
    d[1][1]=e1;
    d[1][2]=0.;
    d[2][0]=0.;
    d[2][1]=0.;
    d[2][2]=e3;
  break;
  default:*/
/*-------------------------------- plane strain, rotational symmetry ---*/
    c1=ym/(1.0+pv);
    b1=c1*pv/(1.0-2.0*pv);
    a1=b1+c1;

    d[0][0]=a1;
    d[0][1]=b1;
    d[0][2]=0.;
    d[0][3]=b1;
    
    d[1][0]=b1;
    d[1][1]=a1;
    d[1][2]=0.;
    d[1][3]=b1;

    d[2][0]=0.;
    d[2][1]=0.;
    d[2][2]=c1/2.;
    d[2][3]=0.;
    
    d[3][0]=b1;
    d[3][1]=b1;
    d[3][2]=0.;
    d[3][3]=a1;
/*  break;
  }*/
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of ale2_mat_linel */
#endif
/*! @} (documentation module close)*/
