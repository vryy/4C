/*======================================================================*/
/*!
\file
\brief Calculates the B-operator matrix for THERM2 element

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>
*/
#ifdef D_THERM2


/*----------------------------------------------------------------------*/
/* headers */
#include "../headers/standardtypes.h"
#include "therm2.h"


/*!
\addtogroup WALL1
*//*! @{ (documentation module open)*/


/*======================================================================*/
/*!
\brief Calculate B-operator matrix at point r,s

\param   **bop      DOUBLE  (o)  B-operator at (r,s)
\param   **deriv    DOUBLE  (i)  natural derivatives of shapge functions
                                 at the point (r,s)
\param   **xjm      DOUBLE  (i)  Jacobi matrix at (r,s)
\param   det        DOUBLE  (i)  Jacobi determinant at (r,s)
\param   iel        INT     (i)  number of element nodes
\return void

\author bborn
\date 03/06
*/
void th2_bop(DOUBLE    **bop,
             DOUBLE    **deriv,
             DOUBLE    **xjm,
             DOUBLE      det,
             INT         iel)
{
  /*--------------------------------------------------------------------*/
  INT inode;
  DOUBLE xji[2][2];

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th2_bop");
#endif

  /*--------------------------------------------------------------------*/
  /* inverse of Jacobian */
  xji[0][0] = xjm[1][1] / det;
  xji[0][1] =-xjm[0][1] / det;
  xji[1][0] =-xjm[1][0] / det;
  xji[1][1] = xjm[0][0] / det;

  /*--------------------------------------------------------------------*/
  /* construct operator B with global derivatives */
  /* loop over element nodes */
  for (inode=0; inode<iel; inode++)
  {
    bop[0][inode] 
      += xji[0][0] * deriv[0][inode]
      +  xji[0][1] * deriv[1][inode];
    bop[1][inode] 
      += xji[1][0] * deriv[0][inode]
      +  xji[1][1] * deriv[1][inode];
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of th2_bop */


/*======================================================================*/
#endif /* end of #ifdef D_THERM2 */
/*! @} (documentation module close)*/
