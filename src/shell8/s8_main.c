#include "../headers/standardtypes.h"
#include "shell8.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

/*----------------------------------------------------------------------*
 | main shell8 control routine                            m.gee 6/01    |
 *----------------------------------------------------------------------*/
void shell8(      FIELD     *actfield,
                  PARTITION *actpart,
                  INTRA     *actintra,
                  ELEMENT   *ele,
                  ARRAY     *estif_global,
                  ARRAY     *emass_global,
                  double    *global_vec,
                  int        global_numeq,
                  int        kstep,
            const int        option
            )
{
int          i;
int          imyrank;
int          inprocs;
S8_DATA      actdata;
MATERIAL    *actmat;

#ifdef DEBUG 
dstrc_enter("shell8");
#endif
/*----------------------------------------------------------------------*/
/*------------------------------------------------- switch to do option */
switch (option)
{
case 0:/*---------------------- init the element routines and directors */
   s8init(actfield);
   s8static_ke(NULL,NULL,NULL,NULL,NULL,0,0,1);
   s8static_keug(NULL,NULL,NULL,NULL,NULL,0,0,1);
   s8eleload(NULL,NULL,NULL,NULL,0,1);
   s8jaco(NULL,NULL,NULL,NULL,NULL,NULL,0.0,0,NULL,NULL,1);
   s8_stress(NULL,NULL,NULL,0,1);
break;/*----------------------------------------------------------------*/
case 1:/*---------------------------- calculate linear stiffness matrix */
   actmat = &(mat[ele->mat-1]);
   s8static_ke(ele,&actdata,actmat,estif_global,NULL,0,0,0);
break;/*----------------------------------------------------------------*/
case 2:/*--------------------------calculate nonlinear stiffness matrix */
   actmat = &(mat[ele->mat-1]);
   s8static_keug(ele,&actdata,actmat,estif_global,global_vec,global_numeq,kstep,0);
break;/*----------------------------------------------------------------*/
case 3:/*------------------- calculate linear stiffness and mass matrix */
break;/*----------------------------------------------------------------*/
case 4:/*---------------- calculate nonlinear stiffness and mass matrix */
break;/*----------------------------------------------------------------*/
case 5:/*------------------------- calculate stresses in a certain step */
   imyrank = actintra->intra_rank;
   if (imyrank==ele->proc) 
   {
      actmat = &(mat[ele->mat-1]);
      s8_stress(ele,&actdata,actmat,kstep,0);
   }
break;/*----------------------------------------------------------------*/
case 6:/*----------------------- calculate load vector of element loads */
   imyrank = actintra->intra_rank;
   if (imyrank==ele->proc) 
   {
      actmat = &(mat[ele->mat-1]);
      s8eleload(ele,&actdata,actmat,global_vec,global_numeq,0);
   }
break;/*----------------------------------------------------------------*/
case 7:/*--------------------------------- reduce stresses to all procs */
   /*------------------------------------- not necessary in sequentiell */
   if (actintra->intra_nprocs==1) goto end;
   s8_stress_reduce(actfield,actpart,actintra,kstep);      
break;/*----------------------------------------------------------------*/
default:
   dserror("action unknown");
break;
}
/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of shell8 */
