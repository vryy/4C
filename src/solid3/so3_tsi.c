/*======================================================================*/
/*!
\file
\brief Extra routines for thermo-structure interaction (TSI)

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>

\author bborn
\date 08/06
*/
#ifdef D_SOLID3
#ifdef D_TSI

/*----------------------------------------------------------------------*/
/* headers */
#include "../headers/standardtypes.h"
#include "solid3.h"
#include "../therm3/therm3.h"

/*======================================================================*/
/*!
\brief Compute strain due to temperature elevation

\param  container CONTAINER*  (i)    container
\param  ele       ELEMENT*    (i)    pointer to current (wall) element
\param  r         DOUBLE      (i)    r-coord of (r,s,t)
\param  s         DOUBLE      (i)    s-coord of (r,s,t)
\param  t         DOUBLE      (i)    s-coord of (r,s,t)
\param  temper    DOUBLE*     (o)    temperarture at (r,s,t)

\return void

\author bborn
\date 03/07
*/
void so3_tsi_temper(const CONTAINER *container,
                    const ELEMENT *ele,
                    const DOUBLE r, 
                    const DOUBLE s,
                    const DOUBLE t,
                    DOUBLE *temper)
{

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_tsi_temper");
#endif

  /*--------------------------------------------------------------------*/
  /* temperature at Gauss point (r,s,t) */
  switch (ele->e.so3->therm_ele->eltyp)
  {
#ifdef D_THERM3
    case el_therm3:
      th3_temper_cal(container, ele->e.so3->therm_ele, r, s, t, temper);
      break;
#endif
    default:
      dserror("Associated thermal element does not exist!");
      break;
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}


/*======================================================================*/
#endif  /* end of #ifdef D_TSI */
#endif  /* end of #ifdef D_SOLID3 */
