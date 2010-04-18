/*======================================================================*/
/*!
\file
\brief Select proper material law

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>

\author bborn
\date 03/06
*/
#ifndef CCADISCRET
#ifdef D_THERM2


/*----------------------------------------------------------------------*/
/* headers */
#include "../headers/standardtypes.h"
#include "therm2.h"

/*!
\addtogroup THERM2
*//*! @{ (documentation module open)*/


/*======================================================================*/
/*!
\brief Select proper material law

\param *ele       ELEMENT   (i)   pointer to current element
\param *mat       MATERIAL  (i)   pointer to current material
\param **bop      DOUBLE    (i)   B-operator
\param ip         INT       (i)   current Gauss point index
\param *heatflux  DOUBLE    (o)   heat flux
\param **cmat     DOUBLE    (o)   constitutive matrix
\return void

\author bborn
\date 03/06
*/
void th2_mat_sel(ELEMENT   *ele,
                 MATERIAL  *mat,
                 DOUBLE **bop,
                 INT ip,
                 DOUBLE *heatflux,
                 DOUBLE **cmat)
{
  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th2_mat_sel");
#endif

  /*--------------------------------------------------------------------*/
  /* the material law (it's a material world!) */
  switch (mat->mattyp)
  {
    case m_th_fourier_iso:
      th2_matlin_iso(mat->m.th_fourier_iso->conduct,
                     ele,
                     bop,
                     heatflux,
                     cmat);
      break;
    default:
      dserror("Type of material law is not applicable");
      break;
  }  /* end of switch (mat->mattyp) */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of th2_mat_sel */


/*======================================================================*/
#endif  /* end of #ifdef D_THERM2 */
#endif
