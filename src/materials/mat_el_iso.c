/*!----------------------------------------------------------------------
\file
\brief contains the routine 
 - mat_el_iso: which calculates the constitutive matrix for a linear 
               elastic, isotropic material (St. Venant-Kirchhoff)
               formulated in cartesian coordinate system, 
               general 3D with the sorting [11,22,33,12,23,13]
*----------------------------------------------------------------------*/
#ifdef D_MAT
#include "../headers/standardtypes.h"
#include "mat_prototypes.h"

/*! 
\addtogroup MAT 
*//*! @{ (documentation module open)*/ 

/*!----------------------------------------------------------------------
\brief consitutive matrix for a linear elastic isotropic material                                    

<pre>                                                            sh 03/03
This routine calculates the constitutive matrix for a linear elastic,
isotropic material
Within this routine, everything is done in a cartesian coordinate system
with the following sorting of stresses and strains:
"brick" [11,22,33,12,23,13]
</pre>
\param       youngs   DOUBLE  (i)   young's modulus 
\param possionratio   DOUBLE  (i)   poisson's ratio 
\param      **d       DOUBLE  (i)   constitutive matrix

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: 

*----------------------------------------------------------------------*/
void mat_el_iso(DOUBLE   youngs,
                DOUBLE   possionratio,
                DOUBLE **d)
{
DOUBLE d1,d2,d3;
DOUBLE ym,pv;/*------------------------------------------ mat constants */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mat_el_iso");
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
} /* end of mat_el_iso */
/*----------------------------------------------------------------------*/
#endif /*D_MAT*/
/*! @} (documentation module close)*/
