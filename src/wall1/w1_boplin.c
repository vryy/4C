/*!----------------------------------------------------------------------
\file
\brief contains the routine 'w1_boplin' which calculates the
       boplin operator matrix for a wall element at point r,s

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

#if 1

/*!
\addtogroup WALL1
*//*! @{ (documentation module open)*/

/*----------------------------------------------------------------------*
 | calculate boplin operator matrix at point r,s                         |
 *----------------------------------------------------------------------*/
void w1_boplin(DOUBLE    **boplin,
               DOUBLE    **deriv,
               DOUBLE    **xjm,
               DOUBLE      det,
               INT         iel)
{
/*----------------------------------------------------------------------*/
INT inode, dnode;
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
/*----------------------------- get operator boplin of global derivatives -*/
/*-------------- some comments, so that even fluid people are able to
   understand this quickly :-)
   the Boplin looks like
       | Nk,x    0   |
       |   0    Nk,y |
       | Nk,y    0   |
       |  0     Nk,x |
*/
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
/*----------------------------------------------------------------------*/
#endif /*D_WALL1*/
/*! @} (documentation module close)*/
