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
void shell8(FIELD      *actfield,
            PARTITION  *actpart,
            INTRA      *actintra,
            ELEMENT    *ele,
            ARRAY      *estif_global,
            ARRAY      *emass_global,
            ARRAY      *intforce_global,
            int         kstep,
            CALC_ACTION *action)
{
int          i;
int          imyrank;
int          inprocs;

double      *intforce;

S8_DATA      actdata;
MATERIAL    *actmat;

#ifdef DEBUG 
dstrc_enter("shell8");
#endif
/*----------------------------------------------------------------------*/
if (intforce_global)
intforce = intforce_global->a.dv;
/*------------------------------------------------- switch to do option */
switch (*action)
{
/*----------------------------- init the element routines and directors */
case calc_struct_init:
   s8init(actfield);
   s8static_ke(NULL,NULL,NULL,NULL,NULL,0,0,1);
   s8static_keug(NULL,NULL,NULL,NULL,NULL,NULL,0,1);
   s8eleload(NULL,NULL,NULL,NULL,1);
   s8jaco(NULL,NULL,NULL,NULL,NULL,NULL,0.0,0,NULL,NULL,1);
   s8_stress(NULL,NULL,NULL,0,1);
break;/*----------------------------------------------------------------*/
/*----------------------------------- calculate linear stiffness matrix */
case calc_struct_linstiff:
   actmat = &(mat[ele->mat-1]);
   s8static_ke(ele,
               &actdata,
               actmat,
               estif_global,
               NULL,
               0,
               0,
               0);
break;/*----------------------------------------------------------------*/
/*---------------------------------calculate nonlinear stiffness matrix */
case calc_struct_nlnstiff:
   actmat = &(mat[ele->mat-1]);
   s8static_keug(ele,
                 &actdata,
                 actmat,
                 estif_global,
                 NULL,
                 intforce,
                 kstep,
                 0);
break;/*----------------------------------------------------------------*/
/*-------------------------- calculate linear stiffness and mass matrix */
case calc_struct_linstiffmass:
break;/*----------------------------------------------------------------*/
/*----------------------- calculate nonlinear stiffness and mass matrix */
case calc_struct_nlnstiffmass:
   actmat = &(mat[ele->mat-1]);
   s8static_keug(ele,
                 &actdata,
                 actmat,
                 estif_global,
                 emass_global,
                 intforce,
                 kstep,
                 0);
break;/*----------------------------------------------------------------*/
/*-------------------------------- calculate stresses in a certain step */
case calc_struct_stress:
   imyrank = actintra->intra_rank;
   if (imyrank==ele->proc) 
   {
      actmat = &(mat[ele->mat-1]);
      s8_stress(ele,&actdata,actmat,kstep,0);
   }
break;/*----------------------------------------------------------------*/
/*------------------------------ calculate load vector of element loads */
case calc_struct_eleload:
   imyrank = actintra->intra_rank;
   if (imyrank==ele->proc) 
   {
      actmat = &(mat[ele->mat-1]);
      s8eleload(ele,&actdata,actmat,intforce,0);
   }
break;/*----------------------------------------------------------------*/
/*---------------------------------------- reduce stresses to all procs */
case calc_struct_stressreduce:
   /*------------------------------------- not necessary in sequentiell */
   if (actintra->intra_nprocs==1) goto end;
   s8_stress_reduce(actfield,actpart,actintra,kstep);      
break;/*----------------------------------------------------------------*/
/*-----------------------------------------------------update variables */
case calc_struct_update_istep:
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
