/*!----------------------------------------------------------------------
\file
\brief contains the routine 'w1_disd' which computes the displacement
       derivatives for a wall element

*----------------------------------------------------------------------*/
#ifdef D_WALL1
#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"

/*! 
\addtogroup WALL1 
*//*! @{ (documentation module open)*/

/*----------------------------------------------------------------------*
 | compute displacement derivatives                          al 9/01    |
 *----------------------------------------------------------------------*/
void w1_disd(ELEMENT   *ele,
             double   **bop, 
             double    *gop, 
             double    *alpha, 
             WALL_TYPE wtype,
             double   *disd)
{
/*----------------------------------------------------------------------*/
int i,inode, node_start;
double dum;
double xji[2][2];
#ifdef DEBUG 
dstrc_enter("w1_disd");
#endif
/*----------------------------------------------------------------------*/
for (i=0; i<5; i++) disd[i] = 0.0;
for (inode=0; inode<ele->numnp; inode++)
{
  node_start = inode*2;
  disd[0] += bop[0][node_start+0]*ele->node[inode]->sol.a.da[0][0];
  disd[1] += bop[1][node_start+1]*ele->node[inode]->sol.a.da[0][1];
  disd[2] += bop[2][node_start+0]*ele->node[inode]->sol.a.da[0][0];
  disd[3] += bop[2][node_start+1]*ele->node[inode]->sol.a.da[0][1];
} /* end of loop over nodes */
/*------------------ additional strain componentes due to inc. modes ---*/
if(ele->e.w1->modeltype == incomp_mode)
{
  disd[0] += gop[0] * alpha[0] + gop[2] * alpha[2];
  disd[1] += gop[1] * alpha[1] + gop[3] * alpha[3];
  disd[2] += gop[1] * alpha[0] + gop[3] * alpha[2];
  disd[3] += gop[0] * alpha[1] + gop[2] * alpha[3];
}
/*----------------------------------------------------------------------*/
if(wtype==rotat_symmet)
{
  for (inode=0; inode<ele->numnp; inode++)
  {
    node_start = inode*2;
    disd[4] += bop[3][node_start+0]*ele->node[inode]->sol.a.da[0][1];
  } /* end of loop over nodes */
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of w1_disd */
/*----------------------------------------------------------------------*/
#endif /*D_WALL1*/
/*! @} (documentation module close)*/
