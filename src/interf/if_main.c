/*!----------------------------------------------------------------------
\file
\brief contains the routine 'interf' the main routine for the 
 interface element

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0771 - 685-6122
</pre>

*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "interf.h"
#include "interf_prototypes.h"
/*----------------------------------------------------------------------*
 |                                                          ah 06/02    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

/*! 
\addtogroup INTERF 
*//*! @{ (documentation module open)*/


/*!----------------------------------------------------------------------
\brief  control routine for the 1D interface element

<pre>                                                              mn 05/03 
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
void interf(PARTITION   *actpart,
            INTRA       *actintra,
            ELEMENT     *ele,
            ARRAY       *estif_global,
            ARRAY       *emass_global,
            ARRAY       *intforce_global,
            CALC_ACTION *action,
            CONTAINER   *container) 
{
#ifdef D_INTERF
INTERF_DATA  actdata;
MATERIAL    *actmat;

DOUBLE      *intforce;

#ifdef DEBUG 
dstrc_enter("interf");
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
    ifinit(actpart);
    ifstatic_ke(NULL,NULL,NULL,NULL,NULL,NULL,1);  
    if_stress(NULL,NULL,NULL,1);       
    if_write_restart(NULL,NULL,0,NULL,1);
    if_read_restart(NULL,NULL,NULL,1);
 /*   if_eleload(ele,&actdata,intforce,1);   */      
break;/*----------------------------------------------------------------*/
/*----------------------------------- calculate linear stiffness matrix */
case calc_struct_nlnstiff:
   actmat = &(mat[ele->mat-1]);                       
   ifstatic_ke(ele,&actdata,actmat,estif_global,NULL,intforce,0);  
break;/*----------------------------------------------------------------*/
/*-------------------------------- calculate stresses in a certain step */
case calc_struct_stress:
   actmat = &(mat[ele->mat-1]);                   
   if_stress(ele,&actdata,actmat,0);    
break;/*----------------------------------------------------------------*/
/*------------------------------ calculate load vector of element loads */
case calc_struct_eleload:
 /*  imyrank = actintra->intra_rank;               */
 /*   if (imyrank==ele->proc)                      */
 /*   {                                            */
 /*      actmat = &(mat[ele->mat-1]);              */
 /*      if_eleload(ele,&actdata,intforce,0);    */
  /*   }                                           */
break;/*----------------------------------------------------------------*/
/*----------------------- calculate nonlinear stiffness and mass matrix */
case calc_struct_nlnstiffmass:
   actmat = &(mat[ele->mat-1]);
   ifstatic_ke(ele,&actdata,actmat,estif_global,emass_global,intforce,0);  
   if (intforce && container->isdyn && ele->proc == actintra->intra_rank)
     solserv_sol_localassemble(actintra,ele,intforce,1,2);
break;/*----------------------------------------------------------------*/
/*------------------------------------------- update after load step ---*/
case calc_struct_update_istep:
     actmat = &(mat[ele->mat-1]);                       
     if (emass_global) 
     {
       ifstatic_ke(ele,&actdata,actmat,estif_global,emass_global,intforce,2);
     } 
     else
     {
       ifstatic_ke(ele,&actdata,actmat,estif_global,NULL,intforce,2);
     }
break;/*----------------------------------------------------------------*/
case write_restart:
   actmat = &(mat[ele->mat-1]);
   if_write_restart(ele,actmat,container->handsize,container->handles,0);
break;/*----------------------------------------------------------------*/
case  read_restart:
   actmat = &(mat[ele->mat-1]);
   if_read_restart(ele,actmat,container->handles,0);
break;/*----------------------------------------------------------------*/
/*------------------------------------------------------------- default */
default:
   dserror("action unknown");
break;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
#endif /*D_INTERF*/
return; 
} /* end of interf */
/*----------------------------------------------------------------------*/
/*! @} (documentation module close)*/
