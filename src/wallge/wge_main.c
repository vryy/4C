/*!----------------------------------------------------------------------
\file
\brief contains the routine 'wallge' the main routine for the
 gradient enhanced wall element

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0771 - 685-6122
</pre>
*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "wallge.h"
#include "wallge_prototypes.h"
/*----------------------------------------------------------------------*
 |                                                          mn 06/02    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

/*!
\addtogroup WALLGE
*//*! @{ (documentation module open)*/


/*!----------------------------------------------------------------------
\brief  control routine for the gradient enhanced wall element

<pre>                                                              ah 05/03
This routine acts according to the action and either initializes the element
or computes the linear stiffness matrix, internal forces,the stresses or the
right hand side vector.
</pre>

\param *actpart         PARTITION   (i)   my partition
\param *actintra        INTRA       (i)   my intra-communicator
\param *ele             ELEMENT     (i)   my element
\param *estif_global    ARRAY       (i)   global stiffness matrix
\param *emass_global    ARRAY       (i)   global mass matrix
\param *intforce_global ARRAY       (i)   global internal force vector
\param *action          CALC_ACTION (i)   option passed to element
\param *container       CONTAINER   (i/o) contains variables defined in container.h

\return void
\sa calling:   ifinit, ifstatic_ke, if_cal_stress, if_eleload;
    called by: global_calelm();

*----------------------------------------------------------------------*/
void wallge(PARTITION   *actpart,
            INTRA       *actintra,
            ELEMENT     *ele,
            ARRAY       *estif_global,
            ARRAY       *emass_global,
            ARRAY       *intforce_global,
            CALC_ACTION *action,
            CONTAINER   *container)
{
#ifdef D_WALLGE
WALLGE_DATA  actdata;
MATERIAL    *actmat;

INT          imyrank;
DOUBLE      *intforce;

#ifdef DEBUG
dstrc_enter("wallge");
#endif
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
if (intforce_global)
   intforce = intforce_global->a.dv;
/*------------------------------------------------- switch to do option */
switch (*action)
{
/*------------------------------------------- init the element routines */
case calc_struct_init:
    wgeinit(actpart,mat);
    wgestatic_ke(NULL,NULL,NULL,NULL,NULL,NULL,1);
    wge_eleload(ele,&actdata,intforce,1);
# if 0
    wge_stress(NULL,NULL,NULL,1);
# endif
break;/*----------------------------------------------------------------*/
/*----------------------------------- calculate linear stiffness matrix */
case calc_struct_nlnstiff:
      actmat = &(mat[ele->mat-1]);
      wgestatic_ke(ele,&actdata,actmat,estif_global,NULL,intforce,0);
break;/*----------------------------------------------------------------*/
/*-------------------------------- calculate stresses in a certain step */
case calc_struct_stress:
   imyrank = actintra->intra_rank;
   if (imyrank==ele->proc)
   {
     wge_cal_stress(ele);
   }
break;/*----------------------------------------------------------------*/
/*--------------------- calculate external load vector of element loads */
case calc_struct_eleload:
   imyrank = actintra->intra_rank;
   if (imyrank==ele->proc)
   {
     wge_eleload(ele,&actdata,intforce,0);
   }
break;/*----------------------------------------------------------------*/
/*----------------------- calculate nonlinear stiffness and mass matrix */
case calc_struct_nlnstiffmass:
break;/*----------------------------------------------------------------*/
/*------------------------------------------- update after load step ---*/
case calc_struct_update_istep:
      actmat = &(mat[ele->mat-1]);
      wgestatic_ke(ele,&actdata,actmat,estif_global,NULL,intforce,2);
  break;/*--------------------------------------------------------------*/
/*----------------------------------------------------------- do nothig */
case write_restart:
  break;/*--------------------------------------------------------------*/
/*----------------------------------------------------------- do nothig */
case  read_restart:
  break;/*--------------------------------------------------------------*/
/*------------------------------------------------------------- default */
default:
   dserror("action unknown");
break;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

#endif /*D_WALLGE*/
return;
} /* end of wallge */
/*----------------------------------------------------------------------*/
/*! @} (documentation module close)*/
