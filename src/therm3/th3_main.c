/*======================================================================*/
/*!
\file
\brief Main routine for THERM3 element

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>

\author bborn
\date 09/06
*/
#ifndef CCADISCRET
#ifdef D_THERM3


/*----------------------------------------------------------------------*/
/* header files */
#include "../headers/standardtypes.h"
#include "therm3.h"

/*!
\addtogroup THERM3
*//*! @{ (documentation module open)*/


/*----------------------------------------------------------------------*/
/*!
\brief General problem data

struct _GENPROB       genprob; defined in global_control.c

\author bborn
\date 09/06
*/
/*----------------------------------------------------------------------*/
/* extern GENPROB genprob; */


/*----------------------------------------------------------------------*/
/*!
\brief Vector of material laws

defined in global_control.c

\author bborn
\date 09/06
*/
/*----------------------------------------------------------------------*/
extern MATERIAL *mat;


/*----------------------------------------------------------------------*/
/*!
\brief Main THERM3 routine

This routine addresses a manifold of THERM3 routines/functions. These
functions are primarily accessed by the appropriate ACTION keyword.
The ACTION keywords are defined in headers/enum.h.

\author bborn
\date 09/06
*/
/*----------------------------------------------------------------------*/
void therm3(PARTITION *actpart,
            INTRA *actintra,
            ELEMENT *ele,
            ARRAY *estif_global,
            ARRAY *emass_global,
            ARRAY *eforc_global,
            CALC_ACTION *action,
            CONTAINER *container)   /* contains variables defined
                                     * in container.h */
{

  static TH3_DATA th3_data; /* global variable th3_data */

  MATERIAL* actmat;

  INT imyrank;

  /*====================================================================*/
#ifdef DEBUG
  dstrc_enter("therm3");
#endif

  /*--------------------------------------------------------------------*/
  /* switch according to ACTION */
  switch (*action)
  {
    /*------------------------------------------------------------------*/
    /* initialise generally element type */
    case calc_therm_init:
      /* these are only called once! */
      th3_cfg_chkdef();
      th3_intg_init(&(th3_data));
      th3_cfg_init(&(th3_data));
#ifdef TEST_THERM3
      th3_cfg_test(&(th3_data));
#endif
       th3_hflux_init(actpart);
/*      th2_write_restart(NULL, NULL, 0, NULL, 1); */
/*      th2_read_restart(NULL, NULL, NULL, 1); */
#ifdef D_TSI
      th3_temper_init();
#endif
      break;
    /*------------------------------------------------------------------*/
    /* stationary tangent */
    case calc_therm_tang_stat:
      actmat = &(mat[ele->mat-1]);
      th3_lin_tang(container, ele, &(th3_data), actmat,
                   estif_global, NULL, NULL);
      break;
    /*------------------------------------------------------------------*/
    /* in-stationary tangent */
    case calc_therm_tang_instat:
      actmat = &(mat[ele->mat-1]);
      th3_lin_tang(container, ele, &(th3_data), actmat,
                   estif_global, emass_global, eforc_global);
      break;
    /*------------------------------------------------------------------*/
    /* external element nodal heat flux vector */
    case calc_therm_heatload:
      imyrank = actintra->intra_rank;
      actmat = &(mat[ele->mat-1]);
      th3_load_heat(container, ele, &(th3_data), imyrank, eforc_global);
      break;
    /*------------------------------------------------------------------*/
    /* heat flux at Gauss points */
    case calc_therm_heatflux:
      imyrank = actintra->intra_rank;
      actmat = &(mat[ele->mat-1]);
      th3_hflux_cal(container, ele, &(th3_data), actmat);
      break;
    /*------------------------------------------------------------------*/
    /* finalise */
    case calc_therm_final:
      th3_hflux_final(actpart);
#ifdef D_TSI
      th3_temper_final();
#endif
      break;
    /*------------------------------------------------------------------*/
    /* catch */
    default:
      dserror("action unknown");
      break;
  }  /* end of switch (*action) */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of therm3() */

/*======================================================================*/
#endif  /* end of #ifdef D_THERM3 */
/*! @} (documentation module close)*/
#endif
