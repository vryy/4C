/*!-----------------------------------------------------------------------
\file
\brief contains the routine 'saxi_B' which calculates the linear operator
matrix for a symmetric shell element at point xsi=s/l

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

*-----------------------------------------------------------------------*/
#ifdef D_AXISHELL

#include "../headers/standardtypes.h"
#include "axishell.h"
#include "axishell_prototypes.h"

/*!
\addtogroup AXISHELL
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief calculate operator matrix at point xsi=s/l

<pre>                                                             mn 05/03
This routine calcuates the operator matrix B at the given point xsi=s/l
for an axisymmetric shell element.

</pre>
\param   B       DOUBLE[][]  (o)   the calculated operator matrix
\param   xsi     DOUBLE  (i)   the integration point where B is computed
\param   r       DOUBLE  (i)   the radius of the integration point
\param   dl      DOUBLE  (i)   the length of the current element
\param   cosa    DOUBLE  (i)   the cosine of the angle of the current element
\param   sina    DOUBLE  (i)   the sine of the angle of the current element

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: saxi_static_ke(), saxi_cal_stress()

*----------------------------------------------------------------------*/
void saxi_B(
    DOUBLE      B[4][7],
    DOUBLE      xsi,
    DOUBLE      r,
    DOUBLE      dl,
    DOUBLE      cosa,
    DOUBLE      sina
    )
{

  INT i,j;

#ifdef DEBUG
  dstrc_enter("saxi_B");
#endif

  for (i=0; i<4; i++)
  {
    for (j=0; j<7; j++)
    {
      B[i][j] = 0.0;
    }
  }

  /* compute operator B  */
  /* compare with Schalen-Skriptum S.11.8  */
  B[0][0] = (-3.0+4.0*xsi)/dl;
  B[0][3] = (-1.0+4.0*xsi)/dl;
  B[1][0] = (1.0-3.0*xsi+2.0*xsi*xsi)*cosa/r;
  B[1][1] = (1.0-3.0*xsi*xsi+2.0*xsi*xsi*xsi)*sina/r;
  B[1][2] = (1.0-2.0*xsi+xsi*xsi)*sina*dl*xsi/r;
  B[1][3] = (-xsi+2.0*xsi*xsi)*cosa/r;
  B[1][4] = (3.0*xsi*xsi-2.0*xsi*xsi*xsi)*sina/r;
  B[1][5] = (-1.0+xsi)*xsi*xsi*sina*dl/r;
  B[2][1] = (-6.0+12.0*xsi)/dl/dl;
  B[2][2] = (-4.0+6.0*xsi)/dl;
  B[2][4] = -B[2][1];
  B[2][5] = (-2.0+6.0*xsi)/dl;
  B[3][1] = 6.0*xsi*(-1.0+xsi)*cosa/(dl*r);
  B[3][2] = (1.0-4.0*xsi+3.0*xsi*xsi)*cosa/r;
  B[3][4] = -B[3][1];
  B[3][5] = xsi*(-2.0+3.0*xsi)*cosa/r;
  B[0][6] = (4.0-8.0*xsi)/dl;
  B[1][6] = 4.0*xsi*(1.0-xsi)*cosa/r;

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of saxi_B */

/*! @} (documentation module close)*/
#endif /*D_AXISHELL*/
