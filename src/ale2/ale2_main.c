/*!----------------------------------------------------------------------
\file
\brief contains the routine 'ale2', the main routine of the 2d ale element 

<pre>
Maintainer: Christiane Foerster
            foerster@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/foerster/
            0711 - 685-6127
</pre>

*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "ale2.h"
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
\brief  control routine for the 2d ale element

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
\sa calling: ale2_static_ke(); called by: ale_calelm(), ale_rhs()

*----------------------------------------------------------------------*/
void ale2(     PARTITION   *actpart,
               INTRA       *actintra,
               ELEMENT     *ele,
               ARRAY       *estif_global,
               CALC_ACTION *action,
               CONTAINER   *container)
{
#ifdef D_ALE
ALE2_DATA     actdata;
MATERIAL     *actmat;

#ifdef DEBUG 
dstrc_enter("ale2");
#endif
/*----------------------------------------------------------------------*/
/*------------------------------------------------- switch to do option */
switch (*action)
{
/*------------------------------------------- init the element routines */
case calc_ale_init:
   ale2_static_ke(NULL,NULL,NULL,NULL,1);
break;
/*----------------------------------- calculate linear stiffness matrix */
case calc_ale_stiff:
   actmat = &(mat[ele->mat-1]);
   ale2_static_ke(ele,&actdata,actmat,estif_global,0);
break;
/*----------------------------------------------------------- do nothig */
case calc_ale_rhs:
break;  
/*---------------------------------- init nonlinear element routines ---*/
case calc_ale_init_nln:
   ale2_static_ke_stiff(NULL,NULL,NULL,NULL,1,container->quality);
   ale2_static_ke_prestress(NULL,NULL,NULL,NULL,1,NULL,0,container->quality);
break;
/*------------------------------------ calculate nonlinear stiffness ---*/
case calc_ale_stiff_nln:  
   actmat = &(mat[ele->mat-1]);
   ale2_static_ke_stiff(ele,&actdata,actmat,estif_global,0,
			container->quality); 
break;
/*------------- calculate stiffness with prestress for initial steps ---*/
case calc_ale_stiff_stress:  
   actmat = &(mat[ele->mat-1]);
   ale2_static_ke_prestress(ele,&actdata,actmat,estif_global,0,
                            container->dirich,container->global_numeq
			    ,container->quality); 
break;
/*-------------------------- init element routines for two step calc ---*/
case calc_ale_init_step2:  /* init element routines for both steps */
   ale2_static_ke(NULL,NULL,NULL,NULL,1);
   ale2_static_ke_step2(NULL,NULL,NULL,NULL,1,0,NULL,NULL,
                        NULL,NULL); 
break;
/*---------------------------- calculate stiffness for two step calc ---*/
case calc_ale_stiff_step2:
   actmat = &(mat[ele->mat-1]);
   ale2_static_ke_step2(ele,&actdata,actmat,estif_global,0,
                       container->quality,&container->min_stiff,
		       &container->max_stiff,&container->min,
		       &container->max); 
break;
/*--------------------- init element routines for spring stiffnesses ---*/
case calc_ale_init_spring:
   ale2_static_ke_spring(NULL,NULL,0,1);
break;
/*------------------------------------- calculate spring stiffnesses ---*/
case calc_ale_stiff_spring:
   actmat = &(mat[ele->mat-1]);
   ale2_static_ke_spring(ele,estif_global,container->quality,0);
break;
/*-------------------- init element routines for laplace stiffnesses ---*/
case calc_ale_init_laplace:
   ale2_static_ke_laplace(NULL,NULL,NULL,1,
			container->quality);
break;
/*------------------------------------ calculate laplace stiffnesses ---*/
case calc_ale_stiff_laplace:
   ale2_static_ke_laplace(ele,&actdata,estif_global,0,
			container->quality);
break;
/*------------------------------------------------------------- default */
default:
   dserror("action unknown");
break;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
#endif
return; 
} /* end of ale3 */
/*! @} (documentation module close)*/
