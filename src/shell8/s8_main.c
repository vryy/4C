/*!----------------------------------------------------------------------
\file
\brief 

<pre>
Maintainer: Michael Gee
            gee@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/gee/
            0771 - 685-6572
</pre>

*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "shell8.h"

/* 
prototypes from prototypes_sol.h which are necessary here, but prototypes_sol.h
shall not be included here
*/
void solserv_sol_localassemble(INTRA *actintra, ELEMENT *actele, DOUBLE *localvec, INT arraynum,
                              INT place);
void dyn_ekin_local(ELEMENT *actele,ARRAY *emass, CONTAINER  *container);

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
            CALC_ACTION *action,
            CONTAINER  *container)    /* contains variables defined in container.h */
{
#ifdef D_SHELL8
INT          i;
INT          imyrank;

DOUBLE      *intforce;

S8_DATA      actdata;
MATERIAL    *actmat;

#ifdef DEBUG 
dstrc_enter("shell8");
#endif
/*----------------------------------------------------------------------*/
if (intforce_global)
intforce = intforce_global->a.dv;
else
intforce = NULL;
/*------------------------------------------------- switch to do option */
switch (*action)
{
/*----------------------------- init the element routines and directors */
case calc_struct_init:
   s8init(actfield);
   s8static_ke(NULL,NULL,NULL,NULL,NULL,0,0,1);
   s8static_keug(NULL,NULL,NULL,NULL,NULL,NULL,0,1);
   s8static_mass(NULL,NULL,NULL,NULL,NULL,NULL,0,1);
   s8eleload(NULL,NULL,NULL,NULL,1);
   s8jaco(NULL,NULL,NULL,NULL,NULL,NULL,0.0,0,NULL,NULL,1);
   s8_stress(NULL,NULL,NULL,0,1);
break;/*----------------------------------------------------------------*/
/*----------------------------------- calculate linear stiffness matrix */
case calc_struct_linstiff:
   actmat = &(mat[ele->mat-1]);
   s8static_keug(ele,
                 &actdata,
                 actmat,
                 estif_global,
                 NULL,
                 intforce,
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
                 container->kstep,
                 0);
break;/*----------------------------------------------------------------*/
/*---------------------------------calculate nonlinear stiffness matrix */
case calc_struct_internalforce:
   actmat = &(mat[ele->mat-1]);
   s8static_keug(ele,
                 &actdata,
                 actmat,
                 estif_global,
                 NULL,
                 intforce,
                 container->kstep,
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
                 container->kstep,
                 0);
   if (intforce && container->isdyn)
      solserv_sol_localassemble(actintra,ele,intforce,1,2);
break;/*----------------------------------------------------------------*/
/*-------------------------------- calculate stresses in a certain step */
case calc_struct_stress:
   imyrank = actintra->intra_rank;
   if (imyrank==ele->proc) 
   {
      actmat = &(mat[ele->mat-1]);
      s8_stress(ele,&actdata,actmat,container->kstep,0);
   }
break;/*----------------------------------------------------------------*/
/*------------------------------ calculate load vector of element loads */
case calc_struct_eleload:
   imyrank = actintra->intra_rank;
/*   if (imyrank==ele->proc) AL
   {*/
      actmat = &(mat[ele->mat-1]);
      s8eleload(ele,&actdata,actmat,intforce,0);
/*   }*/
break;/*----------------------------------------------------------------*/
/*---------------------------------------- reduce stresses to all procs */
case calc_struct_stressreduce:
   /*------------------------------------- not necessary in sequentiell */
   if (actintra->intra_nprocs==1) goto end;
   s8_stress_reduce(actfield,actpart,actintra,container->kstep);      
break;/*----------------------------------------------------------------*/
/*-----------------------------------------------------update variables */
case calc_struct_update_istep:
   if (container->isdyn && ele->proc == actintra->intra_rank)
   {
   actmat = &(mat[ele->mat-1]);
   s8static_mass(ele,&actdata,actmat,estif_global,emass_global,NULL,container->kstep,0);
   dyn_ekin_local(ele,emass_global,container);
   }
   for (i=0; i<ele->e.s8->nhyb; i++)
   ele->e.s8->oldalfa.a.da[0][i] = ele->e.s8->alfa.a.da[0][i];
break;/*----------------------------------------------------------------*/
/*--------------------------------------set variables back to last step */
case calc_struct_update_stepback:
   for (i=0; i<ele->e.s8->nhyb; i++)
   ele->e.s8->alfa.a.da[0][i] = ele->e.s8->oldalfa.a.da[0][i];
break;/*----------------------------------------------------------------*/
/*--------------------------------------------------------write restart */
case write_restart:
   s8_write_restart(ele,container->handsize,container->handles);
break;/*----------------------------------------------------------------*/
/*---------------------------------------------------------read restart */
case read_restart:
   s8_read_restart(ele,container->handsize,container->handles);
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
#endif
return; 
} /* end of shell8 */
