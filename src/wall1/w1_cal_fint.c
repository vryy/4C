/*!----------------------------------------------------------------------
\file
\brief contains the routine 'w1_fint' which evaluates the internal element
       forces for large def (total Lagr) for a wall element

*----------------------------------------------------------------------*/
#ifdef D_WALL1
#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"
#include "wall1_prototypes.h"

/*! 
\addtogroup WALL1
*//*! @{ (documentation module open)*/

/*----------------------------------------------------------------------*
 | evaluate internal element forces for large def (total Lagr) ah 06/02  |
 *----------------------------------------------------------------------*/
void w1_fint(DOUBLE **stress,     /* 2.PK stresses        */ 
             DOUBLE  *F,          /* Deformation gradient */ 
             DOUBLE **boplin,     /* B-lin-operator       */ 
             DOUBLE  *fint,       /* internal forces      */ 
             DOUBLE   fac,        /* detJ*wr*ws*thickness */ 
             INT      nd)         /* Element-DOF          */
{
/*----------------------------------------------------------------------*/
INT i,j;
DOUBLE fs[4];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("w1_fint");
#endif
/*--------------------------------------------- F*S = Vektor fs[i] -----*/

fs[0] = fac * (stress[0][0] * F[0] + stress[0][2] * F[2]);
fs[1] = fac * (stress[1][1] * F[1] + stress[0][2] * F[3]);
fs[2] = fac * (stress[1][1] * F[2] + stress[0][2] * F[0]);
fs[3] = fac * (stress[0][0] * F[3] + stress[0][2] * F[1]);

/*---------------------------------------------------------------------*/
for(i=0; i<nd; i++)
{
  j=i+1;
  fint[i] += boplin[0][i]*fs[0]+boplin[2][i]*fs[2];    /* even          */
  fint[j] += boplin[1][j]*fs[1]+boplin[3][j]*fs[3];    /* odd           */
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of w1_fint */
/*----------------------------------------------------------------------*/
#endif /*D_WALL1*/
/*! @} (documentation module close)*/
