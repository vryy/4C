/*!----------------------------------------------------------------------
\file
\brief contains the routine 'axishell' the main routine for the 
 axisymmetric shell element

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

*----------------------------------------------------------------------*/
#ifdef D_AXISHELL

#include "../headers/standardtypes.h"
#include "axishell.h"
#include "axishell_prototypes.h"
/*----------------------------------------------------------------------*
 |                                                          mn 06/02    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

/*! 
\addtogroup AXISHELL 
*//*! @{ (documentation module open)*/


/*!----------------------------------------------------------------------
\brief  control routine for the axisymmetric shell element

<pre>                                                              mn 05/03 
This routine acts according to the action and either initializes the element
or computes the linear stiffness matrix, the stresses or the right hand side 
vector.

</pre>
\param *actpart      PARTITION   (i)   my partition
\param *actintra     INTRA       (i)   my intra-communicator 
\param *ele          ELEMENT     (i)   my element
\param *estif_global ARRAY       (i)   global stiffness matrix
\param *action       CALC_ACTION (i)   option passed to element
\param *container    CONTAINER   (i/o) contains variables defined in container.h

\warning *container is not needed here
\return void                                               
\sa calling:   saxiinit, saxistatic_ke, saxi_cal_stress, saxi_eleload; 
    called by: global_calelm();

*----------------------------------------------------------------------*/
void axishell(
    PARTITION   *actpart,
    INTRA       *actintra,
    ELEMENT     *ele,
    ARRAY       *estif_global,
    ARRAY       *intforce_global,
    CALC_ACTION *action,
    CONTAINER   *container
    )
{

  SAXI_DATA    actdata;
  MATERIAL    *actmat;

  INT          imyrank;
  DOUBLE      *intforce = NULL;

#ifdef DEBUG 
  dstrc_enter("axishell");
#endif

  if (intforce_global)
    intforce = intforce_global->a.dv;

  /* switch to do option */
  switch (*action)
  {

    /* init the element routines */
    case calc_struct_init:
      saxiinit(actpart);
      saxistatic_ke(NULL,NULL,NULL,NULL,1);
      saxi_cal_stress(NULL,NULL,NULL,1); 
      saxi_eleload(ele,&actdata,intforce,1);         
      break;

    /* calculate linear stiffness matrix */
    case calc_struct_linstiff:
      actmat = &(mat[ele->mat-1]);
      saxistatic_ke(ele,&actdata,actmat,estif_global,0);
      break;

    /* calculate stresses in a certain step */
    case calc_struct_stress:
      actmat = &(mat[ele->mat-1]);
      saxi_cal_stress(ele,&actdata,actmat,0);
      break;

    /* calculate load vector of element loads */
    case calc_struct_eleload:
      imyrank = actintra->intra_rank;
      if (imyrank==ele->proc) 
      {
        actmat = &(mat[ele->mat-1]);
        saxi_eleload(ele,&actdata,intforce,0);
      }
      break;

    /* do nothig */
    case calc_struct_update_istep:
      break;

    /* do nothig */
    case write_restart:
      break;

    /* do nothig */
    case  read_restart:
      break;

#ifdef D_OPTIM                   /* include optimization code to ccarat */
      /* do nothig */
    case calc_struct_opt_init:
      break;

    /* do nothig */
    case calc_struct_ste:
      break;

    /* do nothig */
    case calc_struct_stm:
      break;

    /* do nothig */
    case calc_struct_dee:
      break;

    /* do nothig */
    case calc_struct_dmc:
      break;

    /* do nothig */
    case update_struct_odens:
      break;
#endif                   /* stop including optimization code to ccarat :*/

    /* default */
    default:
      dserror("action unknown");
      break;
  }

#ifdef DEBUG 
  dstrc_exit();
#endif

  return; 
} /* end of axishell */

/*! @} (documentation module close)*/
#endif /*D_AXISHELL*/
