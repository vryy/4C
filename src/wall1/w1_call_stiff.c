/*!----------------------------------------------------------------------
\file
\brief contains the routine 'w1_keku' which calculates the usual
       stiffness matrix for a wall element (total Lagr)

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0771 - 685-6122
</pre>

*----------------------------------------------------------------------*/
#ifdef D_WALL1
#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"

/*! 
\addtogroup WALL1
*//*! @{ (documentation module open)*/

/*----------------------------------------------------------------------*
 | usual stiffness matrix total lagrangian formulation           al 9/01|
 |----------------------------------------------------------------------|
 | S       -->  ELEMENT STIFFNESS MATRIX                                |
 | BS      -->  DERIVATIVE OPERATOR                                     |
 |              CONTENTS DEPENDS ON CALLING PROGRAM                     |
 |              IF T.L BS INCLUDES INITIAL DISPLACEMENT PARTS           |
 | D       -->  CONSTITUTIVE MATRIX                                     |
 | FAC     -->  INTEGRATION FACTOR                                      |
 | ND      -->  TOTAL NUMBER DEGREES OF FREEDOM OF ELEMENT              |
 | NEPS    -->  ACTUAL NUMBER OF STRAIN COMPONENTS                      |
 |              =3 PLANE STRESS/PLANE STRAIN                            |
 |              =4 AXISYMMETRIC                                         |
 |                                                                      |
 *----------------------------------------------------------------------*/
void w1_keku(DOUBLE  **s, 
             DOUBLE  **bs, 
             DOUBLE  **d, 
             DOUBLE    fac, 
             INT       nd,
             INT       neps)
{
INT            i, j, k, l, m;
DOUBLE         dum;
DOUBLE         db[4];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("w1_keku");
#endif
/*----------------------------------------------------------------------*/
   for (j=0; j<nd; j++)
   {
     for (k=0; k<neps; k++)
     {
      db[k] = 0.0 ;                                                              
       for (l=0; l<neps; l++)
       {
       db[k] = db[k] + d[k][l]*bs[l][j]*fac ;
       }
     }
     for (i=0; i<nd; i++)
     {
       dum = 0.0 ;                                                                
       for (m=0; m<neps; m++)
       {
        dum = dum + bs[m][i]*db[m] ;
       }
        s[i][j] = s[i][j] + dum ;
     }
   }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of w1_keku */
/*----------------------------------------------------------------------*/
#endif /*D_WALL1*/
/*! @} (documentation module close)*/



