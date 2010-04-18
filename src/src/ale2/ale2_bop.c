/*!----------------------------------------------------------------------
\file
\brief contains the routine 'ale2_bop' which calclates the operator
       matrix for a 2d ale element

<pre>
Maintainer: Christiane Foerster
            foerster@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/foerster/
            0711 - 685-6572
</pre>

*----------------------------------------------------------------------*/
#ifdef D_ALE
#include "../headers/standardtypes.h"
#include "ale2.h"

/*!
\addtogroup Ale
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief calculate operator matrix at point r,s

<pre>                                                              mn 06/02
This routine calcuates the operator matrix b at the given point r,s
for an 2D-ale-element.
</pre>
\param **bop     DOUBLE  (o)   the calculated operator matrix
\param **deriv   DOUBLE  (i)   the derivatives of the shape functions
\param **xjm     DOUBLE  (i)   the Jacobian matrix
\param det       DOUBLE  (i)   the determinant of the Jacobian matrix
\param iel       INT     (i)   number of nodes per element

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: ale2_static_ke()

*----------------------------------------------------------------------*/
void ale2_bop(DOUBLE    **bop,
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
/*! @} (documentation module close)*/
