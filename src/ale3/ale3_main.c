/*!----------------------------------------------------------------------
\file
\brief contains the routine 'ale3', the main routine of the 3d ale element

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "ale3.h"
/*----------------------------------------------------------------------*
 |                                                         mn 06/02     |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

/*!
\addtogroup Ale
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief  control routine for the 3d ale element

<pre>                                                              mn 06/02
This routine controles the calculation of the element stiffness, acts
according to the action.

</pre>
\param *actpart      PARTITION   (i)   my partition
\param *actintra     INTRA       (i)   my intra-communicator
\param *ele          ELEMENT     (i)   my element
\param *estif_global ARRAY       (i)   global stiffness matrix
\param *action       CALC_ACTION (i)   option passed to element
\param *container    CONTAINER   (i/o) contains variables defined in container.h

\warning There is nothing special to this routine
\return void
\sa calling: ale3_static_ke; called by: ale_calelm(), ale_rhs()

*----------------------------------------------------------------------*/
void ale3(
    PARTITION   *actpart,
    INTRA       *actintra,
    ELEMENT     *ele,
    ARRAY       *estif_global,
    CALC_ACTION *action,
    CONTAINER   *container
    )
{

#ifdef D_ALE
  ALE3_DATA     actdata;
  MATERIAL    *actmat;

#ifdef DEBUG
  dstrc_enter("ale1");
#endif

  /* switch to do option */
  switch (*action)
  {
    /* init the element routines */
    case calc_ale_init:
      ale3_static_ke(NULL,NULL,NULL,NULL,1);
      break;

    /* calculate linear stiffness matrix */
    case calc_ale_stiff:
      actmat = &(mat[ele->mat-1]);

      if (ele->distyp == hex20 && ele->e.ale3->hex20_red == 1)
        ale3_static_ke_red(ele,&actdata,actmat,estif_global,0);
      else
        ale3_static_ke(ele,&actdata,actmat,estif_global,0);

      break;

    /* do nothig */
    case calc_ale_rhs:
      break;

    /* default */
    default:
      dserror("action unknown");
      break;
  }

#ifdef DEBUG
  dstrc_exit();
#endif

#endif
  return;
} /* end of ale3 */


/*! @} (documentation module close)*/
