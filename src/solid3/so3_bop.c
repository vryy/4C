/*======================================================================*/
/*!
\file
\brief Calculates the (linear) B-operator matrix for SOLID3 element

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/frenzel
            089-289-15240
</pre>
*/
#ifdef D_SOLID3


/*----------------------------------------------------------------------*/
/* headers */
#include "../headers/standardtypes.h"
#include "solid3.h"


/*!
\addtogroup WALL1
*//*! @{ (documentation module open)*/


/*======================================================================*/
/*!
\brief Calculate B-operator matrix at point (r,s,t) (e.g. Gauss point)

\param   enod       INT     (i)  number of element nodes
\param   **deriv    DOUBLE  (i)  natural derivatives of shape functions
\param   **xji      DOUBLE  (i)  Inverse Jacobi matrix at (r,s,t)
\param   **bop      DOUBLE  (o)  B-operator contribution at (r,s,t)

\return void

\author mf
\date 10/06
*/
void so3_lin_bop(const INT        enod,
                 const DOUBLE   **deriv,
                 const DOUBLE   **xji,
                       DOUBLE   **bop)
{
  /*--------------------------------------------------------------------*/
  INT i, j;  /* dimension indices */
  INT inod;  /* current node index */
  DOUBLE bopcomp;  /* temporary B-operator component */
  DOUBLE xji[2][2];

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_lin_bop");
#endif

  /*--------------------------------------------------------------------*/
  /* construct operator B with global derivatives */
  /* loop over element nodes */
  for (inod=0; inod<enod; inod++)
  {
    for (i=0; i<NDIM_SOLID3; i++)
    {
      bopcomp = 0.0;
      for (j=0; j<NDIM_SOLID3; j++)
      {
        bopcomp = bopcomp + xji[i][j] * deriv[inod][j];
      }
      bop[i][inod] = bopcomp;
    }
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of so3_lin_bop */


/*======================================================================*/
#endif /* end of #ifdef D_SOLID3 */
/*! @} (documentation module close)*/
