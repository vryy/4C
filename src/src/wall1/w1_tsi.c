/*======================================================================*/
/*!
\file
\brief Extra routines for thermo-structure-interaction (TSI)

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>

\author bborn
\date 08/06
*/
#ifndef CCADISCRET
#ifdef D_WALL1
#ifdef D_TSI

/*----------------------------------------------------------------------*/
/* headers */
#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"
#include "../therm2/therm2.h"

/*======================================================================*/
/*!
\brief Compute strain due to temperature elevation

\param  *ele     ELEMENT     (i)    pointer to current (wall) element
\param  *mat     MATERIAL    (i)    material type
\param  r        DOUBLE      (i)    r-coord of (r,s)
\param  s        DOUBLE      (i)    s-coord of (r,s)
\param  numstr   INT         (i)    dimension of strain vector 'strain'
\param  *strain  DOUBLE      (i/o)  strain vector

\return void

\author bborn
\date 08/06
*/
void w1_tsi_thstrain(CONTAINER *container,
                     ELEMENT *ele,
                     MATERIAL *mat,
                     DOUBLE r,
                     DOUBLE s,
                     INT numstr,
                     DOUBLE *strain)
{
  DOUBLE tem;  /* temperature at (r,s) */
  DOUBLE thermexpans;  /* coeff/ of thermal expansion */
  ARRAY thstrain_a;  /* thermal strain */
  DOUBLE *thstrain;
  INT i;  /* loop index */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("w1_tsi_thstrain");
#endif

  /*--------------------------------------------------------------------*/
  /* allocate strain vector due to temperature */
  thstrain = amdef("thstrain", &thstrain_a, numstr, 1, "DV");
  amzero(&thstrain_a);

  /*--------------------------------------------------------------------*/
  /* temperature at Gauss point (r,s) */
  th2_temper_cal(container, ele->e.w1->therm_ele, r, s, &tem);

  /*--------------------------------------------------------------------*/
  /* type of material */
  switch (mat->mattyp)
  {
      case m_stvenant:
        /* coefficient of linear thermal expansion */
        thermexpans = mat->m.stvenant->thermexpans;
        /* thermal strain vector */
        thstrain[0] = thermexpans * tem;  /* E_xx */
        thstrain[1] = thermexpans * tem;  /* E_yy */
        thstrain[2] = 0.0;                /* E_xy */
        thstrain[3] = 0.0;                /* E_yx */
        /* add thermal strain to kinematic strain */
        for (i=0; i<numstr; i++)
        {
          strain[i] = strain[i] - thstrain[i];
        }
        break;
      default:
        dserror("Unknown type of material law");
        break;
  }  /* end of switch */

  /*--------------------------------------------------------------------*/
  /* deallocate thermal strain vector */
  amdel(&thstrain_a);

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of w1_tsi_thstrain */


/*======================================================================*/
#endif  /* end of #ifdef D_TSI */
#endif  /* end of #ifdef D_WALL1 */
#endif
