#ifdef D_WALL1
#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"

#if 1
/*----------------------------------------------------------------------*
 | calculate linear operator matrix at point r,s            ah 06/02    |
 *----------------------------------------------------------------------*/
void w1_boplin(double    **boplin,
               double    **deriv,
               double    **xjm,
               double      det,
               int         iel)
{
/*----------------------------------------------------------------------*/
int inode, dnode;
double dum;
double xji[2][2];
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
  dnode = inode*2;
  
  boplin[0][dnode+0] = deriv[0][inode]*xji[0][0] + deriv[1][inode]*xji[0][1];                             
  boplin[1][dnode+1] = deriv[0][inode]*xji[1][0] + deriv[1][inode]*xji[1][1];                           
  boplin[2][dnode+0] = boplin[1][dnode+1];
  boplin[3][dnode+1] = boplin[0][dnode+0];
} /* end of loop over nodes */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of w1_boplin */
#endif

#endif
