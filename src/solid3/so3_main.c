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
/* extern struct _GENPROB genprob; */


/*----------------------------------------------------------------------*/
/*!
\brief Vector of material laws

defined in global_control.c

\author mf
\date 10/06
*/
extern MATERIAL *mat;


/*----------------------------------------------------------------------*/
/*!
\brief global variable so3_data

Contains constant Gauss point coordinates and weights

\author mf
\date 10/06
*/
static SO3_DATA so3_data;


/*----------------------------------------------------------------------*/
/*!
\brief All Gauss point coordinates and weights, shape functions and their 
       parametric derivates readily evaluated at Gauss point.
       
These data will only change if a different discretisation/Gauss point
combination occurs of the element.

\author bborn
\date 12/06
*/
static SO3_GPSHAPEDERIV so3_gpshade;


/*======================================================================*/
/*!
\brief Main SOLID3 routine

This routine addresses a manifold of SOLID3 routines/functions. These
functions are primarily accessed by the appropriate ACTION keyword.
The ACTION keywords are defined in headers/enum.h.

\param actpart         PARTITION*   (i)   my partition
\param actintra        INTRA*       (i)   my intra-communicator
\param ele             ELEMENT*     (i)   my element
\param estif_global    ARRAY*       (i)   global stiffness matrix
\param emass_global    ARRAY*       (i)   global mass matrix
\param eforc_global    ARRAY*       (i)   global force vector
                                          (internal/external)
\param action          CALC_ACTION* (i)   option passed to element
\param container       CONTAINER*   (i)   contains variables defined 
                                          in container.h
\return void

\author mf
\date 10/06
*/
void solid3(PARTITION *actpart,
	    INTRA *actintra,
	    ELEMENT *ele,
	    ARRAY *estif_global,
	    ARRAY *emass_global,
	    ARRAY *eforc_global,
	    CALC_ACTION *action,
	    CONTAINER *container)
{

  MATERIAL *actmat;
  INT imyrank;

  /*====================================================================*/
#ifdef DEBUG
  dstrc_enter("solid3");
#endif

  /*--------------------------------------------------------------------*/
  /* switch according to ACTION */
  switch (*action)
  {
    /*------------------------------------------------------------------*/
    /* initialize --- these are only called once! */
    case calc_struct_init:
      so3_intg_init(&(so3_data));
      so3_cfg_init(&(so3_data));
#ifdef TEST_SOLID3
      so3_cfg_test(&(so3_data));
      so3_shape_test(&(so3_data));
#endif
      so3_shape_gpshade_init(&(so3_gpshade));
      so3_mat_init(actpart, mat);
      so3_stress_init(actpart);
/*      th2_write_restart(NULL, NULL, 0, NULL, 1); */
/*      th2_read_restart(NULL, NULL, NULL, 1); */
      break;
    /*------------------------------------------------------------------*/
    /* linear stiffness matrix */
    case calc_struct_linstiff:
      actmat = &(mat[ele->mat-1]);
      so3_shape_gpshade(ele, &(so3_data), &(so3_gpshade));
      so3_int_fintstifmass(container, ele, &(so3_gpshade), actmat, 
                           NULL, estif_global, NULL);
      break;
    /*------------------------------------------------------------------*/
    /* nonlinear stiffness matrix and internal force vector */
    case calc_struct_nlnstiff:
      actmat = &(mat[ele->mat-1]);
      so3_shape_gpshade(ele, &(so3_data), &(so3_gpshade));
      so3_int_fintstifmass(container, ele, &(so3_gpshade), actmat, 
                           eforc_global, estif_global, NULL);
      break;
    /*------------------------------------------------------------------*/
    /* linear stiffness and consistent mass matrix */
    case calc_struct_linstiffmass:
      actmat = &(mat[ele->mat-1]);
      so3_shape_gpshade(ele, &(so3_data), &(so3_gpshade));
      so3_int_fintstifmass(container, ele, &(so3_gpshade), actmat, 
                           eforc_global, estif_global, emass_global);
      break;
    /*------------------------------------------------------------------*/
    /* calculate linear stiffness and lumped mass matrix */
    case calc_struct_linstifflmass:
      dserror("Lumped mass matrix is not implemented");
      actmat = &(mat[ele->mat-1]);
      so3_shape_gpshade(ele, &(so3_data), &(so3_gpshade));
      so3_int_fintstifmass(container, ele, &(so3_gpshade), actmat, 
                           eforc_global, estif_global, emass_global);
      break;
    /*------------------------------------------------------------------*/
    /* internal forces, nonlinear stiffness and mass matrix */
    case calc_struct_nlnstiffmass:
      actmat = &(mat[ele->mat-1]);
      so3_shape_gpshade(ele, &(so3_data), &(so3_gpshade));
      so3_int_fintstifmass(container, ele, &(so3_gpshade), actmat, 
                           eforc_global, estif_global, emass_global);
      break;
    /*------------------------------------------------------------------*/
    /* internal forces */
    case calc_struct_internalforce:
      actmat = &(mat[ele->mat-1]);
      so3_shape_gpshade(ele, &(so3_data), &(so3_gpshade));
      so3_int_fintstifmass(container, ele, &(so3_gpshade), actmat, 
                           eforc_global, NULL, NULL);
      break;
    /*------------------------------------------------------------------*/
    /* load vector of element loads */
    case calc_struct_eleload:
      imyrank = actintra->intra_rank;
      actmat = &(mat[ele->mat-1]);
      so3_shape_gpshade(ele, &(so3_data), &(so3_gpshade));
      so3_load(ele, &(so3_data), &(so3_gpshade), imyrank, eforc_global);
      break;
    /*------------------------------------------------------------------*/
    /* incremental update of internal variables */
    case calc_struct_update_istep:
      actmat = &(mat[ele->mat-1]);
      /* The update of internal variables is necessary for
       * advanced material laws, e.g. plastic, visco-elastic or
       * visco-plastic materials, or anti-locking means, e.g.
       * enhanced assumed strain. */
      so3_iv_upd(container, ele, actmat);
      break;
    /*------------------------------------------------------------------*/
    /* calculate stresses */
    case calc_struct_stress:
      imyrank = actintra->intra_rank;
      actmat = &(mat[ele->mat-1]);
      so3_shape_gpshade(ele, &(so3_data), &(so3_gpshade));
      so3_stress(container, ele, &(so3_data), &(so3_gpshade), imyrank, actmat);
      break;
    /*------------------------------------------------------------------*/
    /* finalise */
    case calc_struct_final:
      so3_mat_final(actpart, mat);
      so3_stress_final(actpart);
      break;
    /*------------------------------------------------------------------*/
    /* catch errors */
    default:
      dserror("SOLID3: action unknown");
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
