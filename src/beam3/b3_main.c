/*!----------------------------------------------------------------------
\file
\brief contains the routine 'beam3' which is the main beam3 control 
routine

*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "beam3.h"
#include "beam3_prototypes.h"


/*!----------------------------------------------------------------------
\brief contains all material data sets

<pre>                                                              fh 09/02
declared in src/headers/materials.h

</pre>
*----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;


/*! 
\addtogroup BEAM3
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief controls beam3 element

<pre>                                                              fh 09/02
This routine is the main beam3 control routine

</pre>
\param *actfield        FIELD        (i)  actual field
\param *actpart         PARTITION    (i)  actual partition
\param *actintra        INTRA        (i)  actual intra
\param *actele          ELEMENT     (i/o) actual element
\param *estif_global    ARRAY       (i/o) global element stiffness vector
\param *emass_global    ARRAY       (i/o) global element mass vector
\param *intforce_global ARRAY       (i/o) global internal force vector
\param *action          CALC_ACTION (i/o)  actual action
\param *container       CONTAINER   (i/o) container

\warning There is nothing special in this routine
\return void                                               
\sa calling:   b3_init() , b3_cal_ele() , b3_boplin3D, b3_trans_stf() , 
               b3_load() , b3_loadlin() , b3_cal_force() , 
	       b3_mat_plast_mises() , b3_setdirich() , b3_write_restart() ,
	       b3_read_restart()
    called by: calinit() , calelm() 

*----------------------------------------------------------------------*/
void beam3(FIELD       *actfield,
           PARTITION   *actpart,
           INTRA       *actintra,
           ELEMENT     *ele,
           ARRAY       *estif_global,
           ARRAY       *emass_global,
           ARRAY       *intforce_global,
           CALC_ACTION *action,
	   CONTAINER   *container)
{
#ifdef D_BEAM3
B3_DATA      actdata;  /* actual data */
MATERIAL    *actmat;   /* actual material */

INT          imyrank;  
DOUBLE      *intforce; /* internal force vector */

#ifdef DEBUG 
dstrc_enter("beam3");
#endif
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
if (intforce_global)
intforce = intforce_global->a.dv;
/*--- switch to do option ----------------------------------------------*/
switch (*action)
{
/*--- init the element routines ----------------------------------------*/
case calc_struct_init:
   b3_init(actpart, mat);
   b3_cal_ele(NULL,NULL,NULL,NULL,NULL,NULL,1);
   b3_boplin3D(NULL,NULL,NULL,NULL,NULL,NULL,ZERO,ZERO,ZERO,ZERO,0,1);
   b3_trans_stf(NULL,NULL,NULL,NULL,1);
   b3_load(NULL,NULL,NULL,NULL,NULL,1);
   b3_loadlin(NULL,NULL,NULL,ZERO,ZERO,NULL,1);
   b3_cal_force(NULL,NULL,NULL,1);
   b3_mat_plast_mises(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,NULL,NULL,0,NULL,NULL,0,0,1);
   b3_setdirich(actfield);

   
break;/*----------------------------------------------------------------*/
/*----------- calculate linear stiffness matrix ------------------------*/
case calc_struct_linstiff:
   actmat = &(mat[ele->mat-1]);
   b3_cal_ele(ele,&actdata,actmat,estif_global,NULL,action,0);
break;
/*--------------- calculate global element load vector -----------------*/
case calc_struct_eleload:
   imyrank = actintra->intra_rank;
   if (imyrank==ele->proc) 
   {
      actmat = &(mat[ele->mat-1]);
      b3_cal_ele(ele,&actdata,actmat,estif_global,intforce,action,0);
   }
break;
/*--------------- calculate element internal force vector --------------*/
case calc_struct_stress:
   actmat = &(mat[ele->mat-1]);
   b3_cal_ele(ele,&actdata,actmat,estif_global,intforce,action,0);
break;
/*----------- calculate nonlinear stiffness matrix ---------------------*/
case calc_struct_nlnstiff:
   actmat = &(mat[ele->mat-1]);
   b3_cal_ele(ele,&actdata,actmat,estif_global,intforce,action,0);
break;
/*----------- update stresses from last step ---------------------------*/
case calc_struct_update_istep:
   actmat = &(mat[ele->mat-1]);
   b3_cal_ele(ele,&actdata,actmat,estif_global,intforce,action,0);   
break;
/*------------------------------------------------------- write restart */
case write_restart:
   b3_write_restart(ele,container->handsize,container->handles);
break;
/*------------------------------------------------------- read restart  */
case  read_restart:
   b3_read_restart(ele,container->handsize,container->handles);
break;
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
} /* end of beam3 */
/*! @} (documentation module close)*/
