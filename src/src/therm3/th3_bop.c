/*======================================================================*/
/*!
\file
\brief Calculates the (linear) B-operator matrix for THERM3 element

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>
*/
#ifndef CCADISCRET
#ifdef D_THERM3


/*----------------------------------------------------------------------*/
/* headers */
#include "../headers/standardtypes.h"
#include "therm3.h"


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

\author bborn
\date 09/06
*/
void th3_bop(INT        enod,
             DOUBLE     deriv[MAXNOD_THERM3][NDIM_THERM3],
             DOUBLE     xji[NDIM_THERM3][NDIM_THERM3],
             DOUBLE     bop[NUMTMGR_THERM3][NUMDOF_THERM3*MAXNOD_THERM3])
{
  /*--------------------------------------------------------------------*/
  INT i, j;  /* dimension indices */
  INT inod;  /* current node index */
  DOUBLE bopcomp;  /* temporary B-operator component */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th3_bop");
#endif

  /*--------------------------------------------------------------------*/
  /* construct operator B with global derivatives */
  /* loop over element nodes */
  for (inod=0; inod<enod; inod++)
  {
    for (i=0; i<NUMTMGR_THERM3; i++)
    {
      bopcomp = 0.0;
      for (j=0; j<NDIM_THERM3; j++)
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
} /* end of th3_bop */


/*======================================================================*/
#endif /* end of #ifdef D_THERM3 */
/*! @} (documentation module close)*/
#endif
