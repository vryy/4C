/*!----------------------------------------------------------------------
\file
\brief contains the routine 'w1_kg' which calculates the geometric
       stiffness part (total Lagr) for a wall element

*----------------------------------------------------------------------*/
#ifdef D_WALL1
#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"

/*! 
\addtogroup WALL1
*//*! @{ (documentation module open)*/

/*----------------------------------------------------------------------*
 | geometric stiffness part (total lagrange)                    ah 06/02|
 |----------------------------------------------------------------------|
 | kg      -->  ELEMENT STIFFNESS MATRIX kg                             |
 | Boplin  -->  DERIVATIVE OPERATOR                                     |
 | Stress  -->  STRESS MATRIX                                           |
 | fac     -->  INTEGRATION FACTOR: THICK*DET J * WR * WS               |
 | ND      -->  TOTAL NUMBER DEGREES OF FREEDOM OF ELEMENT              |
 | NEPS    -->  ACTUAL NUMBER OF STRAIN COMPONENTS   =4                 |
 |                                                                      |
 *----------------------------------------------------------------------*/
void w1_kg(double  **kg, 
           double  **boplin, 
           double  **stress,
           double    fac, 
           int       nd,
           int       neps)
{
int            i, j, k, l, m;
double         bsb;
double         sb[4];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("w1_kg");
#endif
/*----------------------------------------------------------------------*/
   for (j=0; j<nd; j++)
   {
     for (k=0; k<neps; k++)
     {
      sb[k] = 0.0 ;                                                              
       for (l=0; l<neps; l++)
       {
       sb[k] += stress[k][l]*boplin[l][j]*fac ;
       }
     }
     for (i=0; i<nd; i++)
     {
       bsb = 0.0 ;                                                                
       for (m=0; m<neps; m++)
       {
        bsb += boplin[m][i]*sb[m] ;
       }
       kg[i][j] += bsb ;
     }
   }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of w1_kg */
/*----------------------------------------------------------------------*/
#endif /*D_WALL1*/
/*! @} (documentation module close)*/
