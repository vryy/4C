/*======================================================================*/
/*!
\file
\brief Main routine for THERM2 element

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
	    http://www.lnm.mw.tum.de/Members/bornemann
	    089-289-15237
</pre>

\author bborn
\date 03/06
*/
/*======================================================================*/


#ifndef CCADISCRET
/*----------------------------------------------------------------------*/
/* only for THERM2 */
#ifdef D_THERM2


/*----------------------------------------------------------------------*/
/* header files */
#include "../headers/standardtypes.h"
#include "therm2.h"

/*!
\addtogroup THERM2
*//*! @{ (documentation module open)*/


/*----------------------------------------------------------------------*/
/*!
\brief General problem data

struct _GENPROB       genprob; defined in global_control.c

\author bborn
\date 03/06
*/
/*----------------------------------------------------------------------*/
/* extern struct _GENPROB genprob; */


/*----------------------------------------------------------------------*/
/*!
\brief Vector of material laws

defined in global_control.c

\author bborn
\date 03/06
*/
/*----------------------------------------------------------------------*/
extern struct _MATERIAL *mat;


/*----------------------------------------------------------------------*/
/*!
\brief global variable th2_data

Contains Gauss point coordinates and weights

\author bborn
\date 03/06
*/
/*----------------------------------------------------------------------*/
static TH2_DATA th2_data;


/*----------------------------------------------------------------------*/
/*!
\brief Main THERM2 routine

This routine addresses a manifold of THERM2 routines/functions. These
functions are primarily accessed by the appropriate ACTION keyword.
The ACTION keywords are defined in headers/enum.h.

\author bborn
\date 03/06
*/
/*----------------------------------------------------------------------*/
void therm2(PARTITION *actpart,
	    INTRA *actintra,
	    ELEMENT *ele,
	    ARRAY *estif_global,
	    ARRAY *emass_global,
	    ARRAY *intforce_global,
	    CALC_ACTION *action,
	    CONTAINER *container)   /* contains variables defined
				     * in container.h */
{

  MATERIAL *actmat;

  INT imyrank;
  DOUBLE *intforce;


  /*====================================================================*/
#ifdef DEBUG
  dstrc_enter("therm2");
#endif

  if (intforce_global)
  {
    intforce = intforce_global->a.dv;
  }

  /*--------------------------------------------------------------------*/
  /* switch according to ACTION */
  switch (*action)
  {
    case calc_therm_init:
      /* these are only called once! */
      th2_intg_init(&(th2_data));
      th2_cfg_init();
      th2_lin_init();
      th2_load_init();
      th2_hflux_init(actpart);
/*      th2_write_restart(NULL, NULL, 0, NULL, 1); */
/*      th2_read_restart(NULL, NULL, NULL, 1); */
#ifdef D_TSI
      th2_temper_init();
#endif
      break;
    case calc_therm_tang_stat:
      actmat = &(mat[ele->mat-1]);
      th2_lin_stiff(ele, &(th2_data), actmat, estif_global, NULL, NULL);
      break;
    case calc_therm_heatload:
      imyrank = actintra->intra_rank;
      actmat = &(mat[ele->mat-1]);
      th2_load_heat(ele, &(th2_data), intforce, imyrank);
      break;
    case calc_therm_heatflux:
      imyrank = actintra->intra_rank;
      actmat = &(mat[ele->mat-1]);
      th2_hflux_cal(ele, &(th2_data), actmat, 0);
      break;
    case calc_therm_final:
      /* these are only called once */
      th2_lin_final();
      th2_load_final();
      th2_hflux_final(actpart);
#ifdef D_TSI
      th2_temper_final();
#endif
      break;
    default:
      dserror("action unknown");
      break;
  }  /* end of switch (*action) */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
}  /* end of therm2() */
/*! @} (documentation module close)*/


/*----------------------------------------------------------------------*/
#endif  /* end of #ifdef D_THERM2 */
#endif
