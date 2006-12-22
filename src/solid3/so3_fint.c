/*======================================================================*/
/*!
\file
\brief Internal element force vector

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
\addtogroup SOLID3
*//*! @{ (documentation module open)*/


/*======================================================================*/
/*!
\brief Calculate internal force vector 
       contribution at Gauss point

\param   enod       INT     (i)  number of element nodes
\param   **deriv    DOUBLE  (i)  natural derivatives of shape functions
\param   **xji      DOUBLE  (i)  Inverse Jacobi matrix at (r,s,t)
\param   *grdisv    DOUBLE  (o)  displacement gradient at (r,s,t)

\return void

\author bborn
\date 12/06
*/
void so3_fint(INT enod,
              DOUBLE bopl[NUMSTR_SOLID3][MAXNOD_SOLID3*NUMDOF_SOLID3],
              DOUBLE defgrdm[NUMDFGR_SOLID3][NUMSTR_SOLID3],
              DOUBLE stress[NUMSTRS_SOLID3],
              DOUBLE fint[MAXNOD_SOLID3*NUMDOF_SOLID3])
{
  DOUBLE tmp, tmp2;  /* intermediate sums occ. in loops */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_fint");
#endif

  /*--------------------------------------------------------------------*/
  /* f_int = Bl^T . Fm . Sv */
  for (idof=0; idof<enod*NUMDOF_SOLID3; idof++)
  {
    tmp = 0.0
    for (igrd=0; igrd<NUMDFGR_SOLID3; igrd++)
    {
      tmp2 = 0.0;
      for (ists=0; ists<NUMSTR_SOLID3; ists++)
      {
        tmp2 += defgrdm[ists][igrd] * stress[ists];
      }
      tmp += bopl[igrd][idof] * tmp2;
    }
    fint[idof] = tmp;
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of so3_fint */


/*======================================================================*/
#endif /* end of #ifdef D_SOLID3 */
/*! @} (documentation module close)*/
