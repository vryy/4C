/*!----------------------------------------------------------------------
\file
\brief contains the routine
 - mat_el_orth: which calculates the constitutive matrix for a linear
                elastic, orthotropic material formulated in cartesian
                coordinate system,
                general 3D with the sorting [11,22,33,12,23,13]

<pre>
Maintainer: Stefan Hartmann
            hartmann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hartmann/
            0711 - 685-6120
</pre>

*----------------------------------------------------------------------*/
#ifdef D_MAT
#include "../headers/standardtypes.h"
#include "mat_prototypes.h"

/*!
\addtogroup MAT
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief consitutive matrix for a linear elastic orthotropic material

<pre>                                                            sh 03/03
This routine calculates the constitutive matrix for a linear elastic,
orthotropic material
Within this routine, everything is done in a cartesian coordinate system
with the following sorting of stresses and strains:
"brick" [11,22,33,12,23,13]
</pre>
\param  DOUBLE    emod1,..   (i)  youngs modulus in different directions
\param  DOUBLE    xnue23,..  (i)  poisons ratio ..
\param  DOUBLE    gmod12,..  (i)  shear modulus ..
\param  DOUBLE  **d          (o)  constitutive matrix

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by:

*----------------------------------------------------------------------*/
void mat_el_orth(DOUBLE    emod1,
                 DOUBLE    emod2,
                 DOUBLE    emod3,
                 DOUBLE    xnue12,
                 DOUBLE    xnue23,
                 DOUBLE    xnue13,
                 DOUBLE    gmod12,
                 DOUBLE    gmod23,
                 DOUBLE    gmod13,
                 DOUBLE  **d)
{
DOUBLE xnue21,xnue32,xnue31;
DOUBLE emod;
DOUBLE delta;
/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("mat_el_orth");
#endif
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
} /* end of mat_el_orth */
/*----------------------------------------------------------------------*/
#endif /*D_MAT*/
/*! @} (documentation module close)*/
