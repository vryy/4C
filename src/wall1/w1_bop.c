/*!----------------------------------------------------------------------
\file
\brief contains the routine 'w1_bop' which calculates the operator
       matrix for a wall element

*----------------------------------------------------------------------*/
#ifdef D_WALL1
#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"

/*! 
\addtogroup WALL1 
*//*! @{ (documentation module open)*/

/*----------------------------------------------------------------------*
 | calculate operator matrix at point r,s                    al 9/01    |
 *----------------------------------------------------------------------*/
void w1_bop(DOUBLE    **bop,          
            DOUBLE    **deriv,        
            DOUBLE    **xjm,          
            DOUBLE      det,          
            INT         iel)          
{
/*----------------------------------------------------------------------*/
INT inode, node_start;
DOUBLE dum;
DOUBLE xji[2][2];
#ifdef DEBUG 
dstrc_enter("w1_bop");
#endif
/*---------------------------------------------- inverse of jacobian ---*/
dum = 1.0/det;                                                             
xji[0][0] = xjm[1][1]* dum;
xji[0][1] =-xjm[0][1]* dum;
xji[1][0] =-xjm[1][0]* dum;
xji[1][1] = xjm[0][0]* dum;
/*----------------------------- get operator b of global derivatives ---*/
for (inode=0; inode<iel; inode++)
{
  node_start = inode*2;
  
  bop[0][node_start+0] += xji[0][0] * deriv[0][inode] 
                       +  xji[0][1] * deriv[1][inode];                             
  bop[1][node_start+1] += xji[1][0] * deriv[0][inode]                             
                       +  xji[1][1] * deriv[1][inode];                             
  bop[2][node_start+1] = bop[0][node_start+0];                                                     
  bop[2][node_start+0] = bop[1][node_start+1];                                                     
} /* end of loop over nodes */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of w1_bop */
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | additional operator matrix gop at r,s for incomp modes    ah 9/02    |
 *----------------------------------------------------------------------*/
void w1_gop(DOUBLE    *gop,            /* operator matrix G             */
            DOUBLE    **xjm0,          /* jacobian mtrix at r,s=0       */
            DOUBLE      det0,          /* det J at r,s=0                */
            DOUBLE      det,           /* det J at r,s                  */
            DOUBLE      e1,            /* actual GP coordinate r        */
            DOUBLE      e2)            /* actual GP coordinate s        */
{
/*----------------------------------------------------------------------*/
INT i;
DOUBLE dum;
DOUBLE xji[2][2];
#ifdef DEBUG 
dstrc_enter("w1_gop");
#endif
/*---------------------------------------------- inverse of jacobian ---*/
dum = 1.0/det0;                                                             
xji[0][0] = xjm0[1][1]* dum;
xji[0][1] =-xjm0[0][1]* dum;
xji[1][0] =-xjm0[1][0]* dum;
xji[1][1] = xjm0[0][0]* dum;
/*----------------------------------- calculate additional operator  ---*/
  gop[0] = - xji[0][0] * 2 * e1 * det0 /  det ; 
  gop[1] = - xji[1][0] * 2 * e1 * det0 /  det ; 
  gop[2] = - xji[0][1] * 2 * e2 * det0 /  det ;
  gop[3] = - xji[1][1] * 2 * e2 * det0 /  det ;
/*----------------------------------------------------------------------*/
#ifdef DEBUG             
dstrc_exit();
#endif
return;
} /* end of w1_gop */
/*----------------------------------------------------------------------*/
#endif /*D_WALL1*/
/*! @} (documentation module close)*/
