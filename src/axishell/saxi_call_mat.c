/*!-----------------------------------------------------------------------
\file
\brief contains the routine 'saxi_mat' which is the material law 
 for the axisymmetric shell element

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
\brief  contains the material law for the axisymmetric shell element

<pre>                                                              mn 05/03 
This routine computes the constitutive  matrix of the axisymmetric shell element

</pre>
\param *mat          STVENANT   (i)   my material model
\param  D            DOUBLE[][] (i/o) the constitutive matrix (4/4)
\param thick         DOUBLE     (i)   thickness of the element

\warning There is nothing special to this routine
\return void                                               
\sa calling:   ---; 
    called by: saxi_static_ke(), saxi_cal_stress();

*----------------------------------------------------------------------*/
void saxi_mat(
    STVENANT  *mat, 
    DOUBLE     D[4][4],
    DOUBLE     thick
    )
{

  INT i,j;
  DOUBLE E, mu, fac;

#ifdef DEBUG 
  dstrc_enter("saxi_mat");
#endif

  E  = mat->youngs;      
  mu = mat->possionratio;
  fac= E*thick/(1.0-mu*mu);

  /* Zero out all entries of the 4/4 constitutive matrix */
  for (i=0; i<4; i++)
  {
    for (j=0; j<4; j++)
    {
      D[i][j] = 0.0;
    }
  }

  D[0][0] = fac;
  D[0][1] = fac*mu;
  D[1][0] = D[0][1];
  D[1][1] = D[0][0];
  D[2][2] = fac*thick*thick/12.0;
  D[2][3] = fac*mu*thick*thick/12.0;
  D[3][2] = D[2][3];
  D[3][3] = D[2][2];

#ifdef DEBUG 
  dstrc_exit();
#endif

  return; 
} /* end of saxi_mat */

/*! @} (documentation module close)*/
#endif /*D_AXISHELL*/
