/*!----------------------------------------------------------------------
\file
\brief contains the routine 'wall1' the main routine for the wall element

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0771 - 685-6122
</pre>

*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"

/*! 
\addtogroup WALL1 
*//*! @{ (documentation module open)*/

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;
/*----------------------------------------------------------------------*
 | main wall1  control routine                               al 6/01    |
 *----------------------------------------------------------------------*/
void wall1(PARTITION   *actpart,
           INTRA       *actintra,
           ELEMENT     *ele,
           ARRAY       *estif_global,
           ARRAY       *emass_global,
           ARRAY       *intforce_global,
           CALC_ACTION *action,
           CONTAINER   *container)   /* contains variables defined in container.h */
{
#ifdef D_WALL1
INT  i, iloc;
W1_DATA      actdata;
MATERIAL    *actmat;

INT          imyrank;
DOUBLE      *intforce;
DOUBLE       getval;

#ifdef DEBUG 
dstrc_enter("wall1");
#endif
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
if (intforce_global)
intforce = intforce_global->a.dv;
/*------------------------------------------------- switch to do option */
switch (*action)
{
/*-init the element routines(geolin and geononlin,no matter if I need both)*/
case calc_struct_init:
   w1init(actpart, mat);
   w1static_ke(NULL,NULL,NULL,NULL,NULL,NULL,1);
   w1static_keug(NULL,NULL,NULL,NULL,NULL,NULL,1);
   w1_cal_stress(NULL,NULL,NULL,NULL,NULL,0,1);
   w1_eleload(ele,&actdata,intforce,1,NULL);
   w1_iedg(NULL,NULL,0,1);
break;/*----------------------------------------------------------------*/
/*----------------------------------- calculate linear stiffness matrix */
case calc_struct_linstiff:
   actmat = &(mat[ele->mat-1]);
   w1static_ke(ele,&actdata,actmat,estif_global,NULL,NULL,0);
break;/*----------------------------------------------------------------*/
/*---------------------------------calculate nonlinear stiffness matrix */
case calc_struct_nlnstiff:
   actmat = &(mat[ele->mat-1]);
   if(ele->e.w1->kintype==total_lagr)
   {
      w1static_keug(ele,&actdata,actmat,estif_global,NULL,intforce,0);
   }
   else if(ele->e.w1->kintype==geo_lin)
   {
      w1static_ke(ele,&actdata,actmat,estif_global,NULL,intforce,0);
   }
   else
   {
      dserror("action unknown");
      break;
   }
break;/*----------------------------------------------------------------*/
/*-------------------------- calculate linear stiffness and mass matrix */
case calc_struct_linstiffmass:
      actmat = &(mat[ele->mat-1]);
      w1static_ke(ele,&actdata,actmat,estif_global,emass_global,intforce,0);
break;/*----------------------------------------------------------------*/
/*----------------------- calculate nonlinear stiffness and mass matrix */
case calc_struct_nlnstiffmass:
     actmat = &(mat[ele->mat-1]);
     w1static_keug(ele,&actdata,actmat,estif_global,emass_global,intforce,0);
     if (intforce && container->isdyn && ele->proc == actintra->intra_rank)
      solserv_sol_localassemble(actintra,ele,intforce,1,2);
break;/*----------------------------------------------------------------*/
/*-------------------------------- calculate stresses in a certain step */
case calc_struct_stress:
   imyrank = actintra->intra_rank;
   if (imyrank==ele->proc)
   { 
    actmat = &(mat[ele->mat-1]);
    w1_cal_stress(ele,&actdata,actmat,estif_global,intforce,0,0);
   }
break;/*----------------------------------------------------------------*/
/*------------------------------ calculate load vector of element loads */
case calc_struct_eleload:
   imyrank = actintra->intra_rank;
   actmat = &(mat[ele->mat-1]);
   w1_eleload(ele,&actdata,intforce,0,imyrank);
break;/*----------------------------------------------------------------*/
case calc_struct_fsiload:
   imyrank = actintra->intra_rank;
   actmat = &(mat[ele->mat-1]);
   w1_fsiload(ele,&actdata,intforce,0,imyrank);
break;/*----------------------------------------------------------------*/
/*--------------------------------------- update after incremental step */
case calc_struct_update_istep:
   actmat = &(mat[ele->mat-1]);
   if(ele->e.w1->kintype == geo_lin)
   {
   w1static_ke(ele,&actdata,actmat,estif_global,NULL,intforce,2);
   }
   else
   {
   w1static_keug(ele,&actdata,actmat,estif_global,NULL,intforce,2);
   }
   /* calculate mass matrix and kinetic energy (m.gee 4/03) */
   /* only mass matrix needed, it would be efficient to have a separate */
   /* routine w1static_mass */
   if (container->isdyn && ele->proc == actintra->intra_rank)
   {
   actmat = &(mat[ele->mat-1]);
   w1static_keug(ele,&actdata,actmat,estif_global,emass_global,NULL,0);
   dyn_ekin_local(ele,emass_global,container);
   }
#ifdef GEMM  
   if (container->isdyn) w1_update_history(ele,&actdata,actmat);
#endif
break;/*----------------------------------------------------------------*/
/*------------------------------------------------------- write restart */
case write_restart:
#if 0
  momentan nicht notwendig, da keine notwendige Info im Element vgl. shell8
      w1_write_restart(ele,container->handsize,container->handles);
#endif
break;/*----------------------------------------------------------------*/
/*------------------------------------------------------- write restart */
case  read_restart:
#if 0
  momentan nicht notwendig, da keine notwendige Info im Element vgl. shell8
      w1_read_restart(ele,container->handsize,container->handles);
#endif
break;/*----------------------------------------------------------------*/
#ifdef D_OPTIM                   /* include optimization code to ccarat */
/*-------------- init the element integration routines for optimization */
case calc_struct_opt_init:
  w1_oint(NULL,NULL,NULL,NULL,1);
break;/*----------------------------------------------------------------*/
/*---------------------------------------------- evaluate strain energy */
case calc_struct_ste:
#ifdef PARALLEL 
   if(ele->proc!=actintra->intra_rank) break;
#endif
   actmat = &(mat[ele->mat-1]);
   w1_oint(ele,
           &actdata,
           actmat,
           &getval,
           2 /* flag for strain energy */
           );
   (*container).getvalue += getval;
break;/*----------------------------------------------------------------*/
/*------------------------------------------------------- evaluate mass */
case calc_struct_stm:
#ifdef PARALLEL 
   if(ele->proc!=actintra->intra_rank) break;
#endif
   actmat = &(mat[ele->mat-1]);
   w1_oint(ele,
           &actdata,
           actmat,
           &getval,
           3 /* flag for mass */
           );
   (*container).getvalue += getval;
break;/*----------------------------------------------------------------*/
/*-------------------------------- evaluate derivative of strain energy */
case calc_struct_dee:
#ifdef PARALLEL 
   if(ele->proc!=actintra->intra_rank) break;
#endif
   if(ele->optdata==NULL) break; /* element does not take part in opt. */
   if(ele->optdata[0]==0) break; /* position in variable vector        */
   iloc = ele->optdata[0];
   
   actmat = &(mat[ele->mat-1]);
   w1_oint(ele,
           &actdata,
           actmat,
           &getval,
           4 /* flag for derivative of strain energy */
           );
   (*container).getvector[iloc-1] += getval;
break;/*----------------------------------------------------------------*/
/*------------------------------ evaluate derivative of mass constraint */
case calc_struct_dmc:
#ifdef PARALLEL 
   if(ele->proc!=actintra->intra_rank) break;
#endif
   if(ele->optdata==NULL) break; /* element does not take part in opt. */
   if(ele->optdata[0]==0) break; /* position in variable vector        */
   iloc = ele->optdata[0];
   
   actmat = &(mat[ele->mat-1]);
   w1_oint(ele,
           &actdata,
           actmat,
           &getval,
           5 /* flag for volume */
           );
   (*container).getvector[iloc-1] += getval;
break;/*----------------------------------------------------------------*/
/*---------------------------------------------------- update densities */
case update_struct_odens:
   if(ele->optdata==NULL) break; /* element does not take part in opt. */
   if(ele->optdata[0]==0) break; /* position in variable vector        */
   iloc = ele->optdata[0];

   getval = (*container).getvector[iloc-1];
   ele->e.w1->elewa->matdata[0]   =  getval;
break;/*----------------------------------------------------------------*/
#endif /* stop including optimization code to ccarat :*/
default:
   dserror("action unknown");
break;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
#endif /*D_WALL1*/
return; 
} /* end of wall1 */
/*----------------------------------------------------------------------*/
/*! @} (documentation module close)*/
