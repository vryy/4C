/*!----------------------------------------------------------------------
\file
\brief contains the routine 'w1_defgrad' which calculates the deformation
       gradient and GL strain at gausian point r,s for a wall element

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
 | calc deformation gradient and GL strain  at gausian point r,s  ah06/02|
 *----------------------------------------------------------------------*/
void w1_defgrad(DOUBLE    *F,
                DOUBLE    *strain,
                DOUBLE    **xrefe,
                DOUBLE    **xcure,
                DOUBLE    **boplin,
                INT         iel)
{
/*----------------------------------------------------------------------*/
INT inode;
#ifdef DEBUG 
dstrc_enter("w1_defgrad");
#endif
/*------------------calculate defgrad --------- (Summenschleife->+=) ---*/
F[0]=1;
F[1]=1;
for (inode=0; inode<iel; inode++)
{ 
   F[0] += boplin[0][2*inode]   * (xcure[0][inode] - xrefe[0][inode]);
   F[1] += boplin[1][2*inode+1] * (xcure[1][inode] - xrefe[1][inode]);
   F[2] += boplin[2][2*inode]   * (xcure[0][inode] - xrefe[0][inode]);
   F[3] += boplin[3][2*inode+1] * (xcure[1][inode] - xrefe[1][inode]);
} /* end of loop over nodes */
/*-----------------------calculate Green-Lagrange strain ---------------*/
strain[0]=0.5 * (F[0] * F[0] + F[3] * F[3] - 1);
strain[1]=0.5 * (F[2] * F[2] + F[1] * F[1] - 1);
strain[2]=0.5 * (F[0] * F[2] + F[3] * F[1]);
strain[3]=strain[2];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of w1_defgrad */
/*----------------------------------------------------------------------*/
#endif /*D_WALL1*/
/*! @} (documentation module close)*/
