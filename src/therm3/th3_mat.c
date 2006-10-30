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
\date 09/06
*/
#ifdef D_THERM3


/*----------------------------------------------------------------------*/
/* headers */
#include "../headers/standardtypes.h"
#include "therm3.h"

/*!
\addtogroup THERM3
*//*! @{ (documentation module open)*/


/*======================================================================*/
/*!
\brief Select proper material law

\param *cont      CONTAIMER (i)   container data
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
void th3_mat_sel(CONTAINER *cont,
                 ELEMENT *ele,
                 MATERIAL *mat,
                 DOUBLE bop[NDIM_THERM3][NUMDOF_THERM3*MAXNOD_THERM3],
                 INT ip,
                 DOUBLE heatflux[NUMHFLX_THERM3],
                 DOUBLE cmat[NUMHFLX_THERM3][NUMTMGR_THERM3])
{
  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th3_mat_sel");
#endif

  /*--------------------------------------------------------------------*/
  /* the material law (it's a material world!) */
  switch (mat->mattyp)
  {
    case m_th_fourier_iso:
      th3_matlin_iso(cont,
                     mat->m.th_fourier_iso->conduct,
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
}  /* end of th3_mat_sel */


/*======================================================================*/
#endif  /* end of #ifdef D_THERM3 */
/*! @} (documentation module close) */
