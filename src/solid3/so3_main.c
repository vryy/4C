/*======================================================================*/
/*!
\file
\brief Main routine for SOLID3 element

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
	    http://www.lnm.mw.tum.de/Members/frenzel
	    089-289-15240
</pre>

\author mf
\date 10/06
*/
#ifdef D_SOLID3


/*----------------------------------------------------------------------*/
/* header files */
#include "../headers/standardtypes.h"
#include "solid3.h"

/*!
\addtogroup SOLID3
*//*! @{ (documentation module open)*/


/*----------------------------------------------------------------------*/
/*!
\brief General problem data

struct _GENPROB       genprob; defined in global_control.c

\author mf
\date 10/06
*/
/*----------------------------------------------------------------------*/
/* extern struct _GENPROB genprob; */


/*----------------------------------------------------------------------*/
/*!
\brief Vector of material laws

defined in global_control.c

\author mf
\date 10/06
*/
/*----------------------------------------------------------------------*/
extern struct _MATERIAL *mat;


/*----------------------------------------------------------------------*/
/*!
\brief global variable so3_data

Contains Gauss point coordinates and weights

\author mf
\date 10/06
*/
/*----------------------------------------------------------------------*/
static SOLID3_DATA so3_data;


/*----------------------------------------------------------------------*/
/*!
\brief Main SOLID3 routine

This routine addresses a manifold of SOLID3 routines/functions. These
functions are primarily accessed by the appropriate ACTION keyword.
The ACTION keywords are defined in headers/enum.h.

\author mf
\date 10/06
*/
/*----------------------------------------------------------------------*/
void solid3(PARTITION *actpart,
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
  dstrc_enter("solid3");
#endif

  if (intforce_global)
  {
    intforce = intforce_global->a.dv;
  }

  /*--------------------------------------------------------------------*/
  /* switch according to ACTION */
  switch (*action)
  {
/* initialize   --------------------------------------------------------*/
    case calc_struct_init:
      /* these are only called once! */
      so3_intg_init(&(so3_data));
      so3_cfg_init(&(so3_data));
/*       th2_lin_init(); */
/*       th2_load_init(); */
/*       th2_hflux_init(actpart); */
/*      th2_write_restart(NULL, NULL, 0, NULL, 1); */
/*      th2_read_restart(NULL, NULL, NULL, 1); */
      break;
/* linear stiffness matrix ----------------------------------------------*/
    case calc_struct_linstiff:
      actmat = &(mat[ele->mat-1]);
      so3_stiff(ele, &(so3_data), actmat, estif_global, NULL, NULL);
      break;
/* nonlinear stiffness matrix -------------------------------------------*/
    case calc_struct_nlnstiff:
      actmat = &(mat[ele->mat-1]);
/*       so3_stiff(ele, &(so3_data), actmat, estif_global, NULL, NULL); */
      break;
/* nonlinear stiffness and mass matrix ----------------------------------*/
    case calc_struct_nlnstiffmass:
      actmat = &(mat[ele->mat-1]);
/*       so3_stiff(ele, &(so3_data), actmat, estif_global, NULL, NULL); */
      break;
/* load vector of element loads -----------------------------------------*/
    case calc_struct_eleload:
      imyrank = actintra->intra_rank;
      actmat = &(mat[ele->mat-1]);
/*       so3_load(ele, &(so3_data), intforce, imyrank); */
      break;
/* calculate stresses --------------------------------------------------*/
    case calc_struct_stress:
      imyrank = actintra->intra_rank;
      actmat = &(mat[ele->mat-1]);
/*       so3_stress(ele, &(so3_data), actmat, 0); */
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
}  /* end of therm3() */

/*======================================================================*/
#endif  /* end of #ifdef D_SOLID3 */
/*! @} (documentation module close)*/
