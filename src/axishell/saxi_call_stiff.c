/*!-----------------------------------------------------------------------
\file
\brief contains the routine 'saxi_keku' which calculates the usual
stiffness matrix at one integration point for a axisymmetric shell element

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
\brief calculates usual stiffness matrix in total lagrangian formulation

<pre>                                                              mn 05/03
This routine calculates usual stiffness matrix in total lagrangian
formulation.

</pre>
\param   s       DOUBLE[][]    (o)  element stiffness matrix
\param   b       DOUBLE[][]    (i)  derivative operator
\param   d       DOUBLE[][]    (i)  constitutive matrix
\param   work1   DOUBLE[][]    (i)  working array (4,7)
\param   work2   DOUBLE[][]    (i)  working array (7,7)
\param   fac     DOUBLE    (i)  integration factor
\param   r       DOUBLE    (i)  radius of current integration point
\param   dl      DOUBLE    (i)  length of current element

\warning There is nothing special to this routine
\return void
\sa calling:  ---;
    caled by: saxi_static_ke();

*----------------------------------------------------------------------*/
void saxi_keku(
    DOUBLE    s[7][7],
    DOUBLE    b[4][7],
    DOUBLE    d[4][4],
    DOUBLE    work1[4][7],
    DOUBLE    work2[7][7],
    DOUBLE    fac,
    DOUBLE    r,
    DOUBLE    dl
    )
{

  INT            i, j, k;

#ifdef DEBUG
  dstrc_enter("saxi_keku");
#endif

  /*  make multiplication: work1 = D * B */
  for (k=0; k<7; k++)
  {
    for (i=0; i<4; i++)
    {
      work1[i][k] = 0.0;
      for (j=0; j<4; j++)
      {
        work1[i][k] += d[i][j] * b[j][k];
      }
    }
  }

  /*  make multiplication: work2 = B^T * work1 */
  for (k=0; k<7; k++)
  {
    for (i=0; i<7; i++)
    {
      work2[i][k] = 0.0;
      for (j=0; j<4; j++)
      {
        work2[i][k] += b[j][i] * work1[j][k];
      }
    }
  }

  /* compute: estiff[i][j] += fac*r*dl/2*work2[i][j] */
  for (i=0; i<7; i++)
  {
    for (j=0; j<7; j++)
    {
      s[i][j] += fac*r*dl*work2[i][j]/1.0;
      /* im fortran code wird durch 2 dividiert, dieser faktor ist bei
         mir schon in den Gewichten der numerischen integration */
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of saxi_keku */

/*! @} (documentation module close)*/
#endif /*D_AXISHELL*/
