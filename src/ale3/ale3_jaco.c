/*!----------------------------------------------------------------------
\file
\brief contains the routine 'ale3_jaco' which calculates the Jacobian 
matrix for a 3d ale element

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

*----------------------------------------------------------------------*/
#ifdef D_ALE
#include "../headers/standardtypes.h"
#include "ale3.h"

/*! 
\addtogroup Ale 
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief  calculate the Jacobian matrix  

<pre>                                                              mn 06/02 
This routine calculates the Jacobian matrix  at a point r,s,t for 
a 3D ale element.

</pre>
\param deriv    DOUBLE[][] (i)   derivatives of the shape functions
\param xjm      DOUBLE[][] (o)   the Jacobian matrix
\param *det     DOUBLE     (i)   determinant of the Jacobian matrix
\param *ele     ELEMENT    (i)   the element
\param iel      INT        (i)   number of nodes of the element

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: ale3_static_ke

*----------------------------------------------------------------------*/
void ale3_jaco(
    DOUBLE      deriv[3][MAXNOD_BRICK1],
    DOUBLE      xjm[3][3],
    DOUBLE     *det,
    ELEMENT    *ele,
    INT         iel
    )
{

  INT i,j,l;
  DOUBLE dum;

#ifdef DEBUG 
  dstrc_enter("ale3_jaco");
#endif

  /* determine jacobian at point r,s,t */       
  for (i=0; i<3; i++)
  {
    for (j=0; j<3; j++)
    {
      dum=0.0;
      for (l=0; l<iel; l++)
      {
        dum += deriv[i][l]*ele->node[l]->x[j];
      }
      xjm[i][j]=dum;
    }
  }

  /* determinant of jacobian */        
  *det = xjm[0][0]*xjm[1][1]*xjm[2][2]+
    xjm[0][1]*xjm[1][2]*xjm[2][0]+
    xjm[0][2]*xjm[1][0]*xjm[2][1]-
    xjm[0][2]*xjm[1][1]*xjm[2][0]-
    xjm[0][0]*xjm[1][2]*xjm[2][1]-
    xjm[0][1]*xjm[1][0]*xjm[2][2];

  if (*det<0.0) dserror("NEGATIVE JACOBIAN DETERMINANT");         

#ifdef DEBUG 
  dstrc_exit();
#endif

  return;
} /* end of ale3_jaco */

#endif
/*! @} (documentation module close)*/
