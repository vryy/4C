/*!----------------------------------------------------------------------
\file
\brief contains the routine 'c1_jaco' which calculates the Jacobian
matrix for a 3D hex element

<pre>
Maintainer: Andreas Lipka
            lipka@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/lipka/
            0771 - 685-6575
</pre>

*----------------------------------------------------------------------*/
#ifdef D_BRICK1

#include "../headers/standardtypes.h"
#include "brick1.h"
#include "brick1_prototypes.h"

/*!
\addtogroup BRICK1
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief calculate the Jacobian matrix

<pre>                                                              al 06/02
This routine calcuates the Jacobian matrix for a 3D-hex-element.

</pre>
\param deriv    DOUBLE*     (i)   derivatives of the shape functions
\param xjm      DOUBLE**    (o)   the Jacobian matrix
\param det      DOUBLE*     (i)   determinant of the Jacobian matrix
\param ele      ELEMENT*    (i)   the element
\param iel      INT         (i)   number of nodes of the element

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | calculate operator matrix at point r,s,t                  al 9/01    |
 *----------------------------------------------------------------------*/
void c1_jaco(DOUBLE  **deriv,
             DOUBLE    **xjm,
             DOUBLE     *det,
             DOUBLE    *xyze,
             INT         iel)
{
/*----------------------------------------------------------------------*/
INT i,j,l,pc;
DOUBLE dum;
#ifdef DEBUG
dstrc_enter("c1_jaco");
#endif
/*-------------------------------- determine jacobian at point r,s,t ---*/
  for (i=0; i<3; i++)
  {
    for (j=0; j<3; j++)
    {
      dum = 0.;
      pc  = j;
      for (l=0; l<iel; l++)
      {
        dum += deriv[i][l] * xyze[pc];
        pc+=3;
      }
      xjm[i][j]=dum;
    }
  }
/*------------------------------------------ determinant of jacobian ---*/
  *det = xjm[0][0]*xjm[1][1]*xjm[2][2]+
         xjm[0][1]*xjm[1][2]*xjm[2][0]+
         xjm[0][2]*xjm[1][0]*xjm[2][1]-
         xjm[0][2]*xjm[1][1]*xjm[2][0]-
         xjm[0][0]*xjm[1][2]*xjm[2][1]-
         xjm[0][1]*xjm[1][0]*xjm[2][2];

  if (*det<0.0)
  {
     dserror("NEGATIVE JACOBIAN DETERMINANT");
  }
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of c1_jaco */
/*----------------------------------------------------------------------*/
#endif
/*! @} (documentation module close)*/
