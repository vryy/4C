#ifdef D_ALE
#include "../headers/standardtypes.h"
#include "ale2.h"
/*----------------------------------------------------------------------*
 | calculate operator matrix at point r,s                    al 9/01    |
 *----------------------------------------------------------------------*/
void ale2_bop(double    **bop,
              double    **deriv,
              double    **xjm,
              double      det,
              int         iel)
{
/*----------------------------------------------------------------------*/
int inode, node_start;
double dum;
double xji[2][2];
#ifdef DEBUG 
dstrc_enter("ale2_bop");
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
} /* end of ale2_bop */
/*----------------------------------------------------------------------*/
#endif
