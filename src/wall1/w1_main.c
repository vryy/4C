/*!----------------------------------------------------------------------
\file
\brief contains the routine 'wall1' the main routine for the wall element

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0711 - 685-6122
</pre>

*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"
#include "../struct2_ml/s2ml_prototypes.h"

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
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

#ifdef D_MLSTRUCT
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | struct _GENPROB       genprob; defined in global_control.c           |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
#endif /* D_MLSTRUCT */

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
W1_DATA      actdata;
MATERIAL    *actmat;

#ifdef D_SSI
INT          i;
#endif  /*D_SSI*/
INT          imyrank;
DOUBLE      *intforce;

#ifdef D_OPTIM                   /* include optimization code to ccarat */
DOUBLE       getval;
INT          iloc;
#endif /* D_OPTIM*/
ARRAY_POSITION* ipos;

#ifdef DEBUG
dstrc_enter("wall1");
#endif

if (container!=NULL)
  ipos = &(field[genprob.numsf].dis[container->disnum].ipos);
else
  ipos = NULL;

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
if (intforce_global)
intforce = intforce_global->a.dv;
/*------------------------------------------------- switch to do option */
switch (*action)
{
/*-init the element routines(geolin and geononlin,no matter if I need both)*/
case calc_struct_init:
#ifdef D_MLSTRUCT
   if (genprob.multisc_struct == 1)
   {
     printf("structural Multiscale only for BILINEAR, RECTENGULAR Makroelements and submesh-left-down-corner-node must have coordinates 0.0/0.0\n");
     printf("For more than 3000 submesh-DOF's the size of the fields ele.e.w1.fint_mi,..for rueckrechnung in static condensation has to be enlarged\n");
     w1init(actpart, mat);
     s2ml_static_ke(NULL,NULL,NULL,NULL,NULL,NULL,1);
     w1_eleload(NULL,NULL,NULL,1,NULL);
     s2ml_bopstraintonode(NULL,NULL,0.0,1);
     s2ml_stiff_wall(NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,0,1);
     s2ml_stiff_interf(NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,0,1);
   }
   else
#endif /* D_MLSTRUCT */
   {
     w1init(actpart, mat);
     w1static_ke(NULL,NULL,NULL,NULL,NULL,NULL,1);
     w1static_keug(NULL,NULL,NULL,NULL,NULL,NULL,1,NULL);
     w1_cal_stress(NULL,NULL,NULL,0,1);
     w1_eleload(ele,&actdata,intforce,1,-1);
     w1_iedg(NULL,NULL,0,1);
     w1_write_restart(NULL,NULL,0,NULL,1);
     w1_read_restart(NULL,NULL,NULL,1);
   }/* end of else: if (genprob.multisc_struct == 1) */
#ifndef LOCALSYSTEMS_ST
   dsassert(ele->locsys==locsys_no,
            "locsys not compiled for WALL1 element, define LOCALSYSTEMS_ST\n");
#endif
break;/*----------------------------------------------------------------*/
/*----------------------------------- calculate linear stiffness matrix */
case calc_struct_linstiff:
   actmat = &(mat[ele->mat-1]);
   w1static_ke(ele,&actdata,actmat,estif_global,NULL,NULL,0);
   /* rotate components at Dirichlet nodes into local system */
#ifdef LOCALSYSTEMS_ST
   if (ele->locsys == locsys_yes)
   {
     locsys_trans_equant_dirich(ele, estif_global, NULL, NULL, NULL);
   }
#endif
break;/*----------------------------------------------------------------*/
/*---------------------------------calculate nonlinear stiffness matrix */
case calc_struct_nlnstiff:
#ifdef D_MLSTRUCT
   if (genprob.multisc_struct == 1)
   {
     actmat = &(mat[ele->mat-1]);
     s2ml_static_ke(ele,&actdata,actmat,estif_global,NULL,intforce,0);
   }
   else
#endif /* D_MLSTRUCT */
   {
     actmat = &(mat[ele->mat-1]);
     if(ele->e.w1->kintype==total_lagr)
     {
       w1static_keug(ele,&actdata,actmat,estif_global,NULL,intforce,0,
                     container);
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
   }/* end of else: if (genprob.multisc_struct == 1) */
   /* rotate components at Dirichlet nodes into local system */
#ifdef LOCALSYSTEMS_ST
   if (ele->locsys == locsys_yes)
   {
     locsys_trans_equant_dirich(ele, estif_global, NULL, NULL, intforce);
   }
#endif
break;/*----------------------------------------------------------------*/
/*---------------------------------- calculate nonlinear internal force */
case calc_struct_internalforce:
   actmat = &(mat[ele->mat-1]);
   if(ele->e.w1->kintype==total_lagr)
   {
     w1static_keug(ele,&actdata,actmat,estif_global,NULL,intforce,0,
                   container);
   }
   else 
   { 
     if(ele->e.w1->kintype==geo_lin)
     {
       w1static_ke(ele,&actdata,actmat,estif_global,NULL,intforce,0);
     }
     else
     {
       dserror("action unknown");
       break;
     }
   }
   /* rotate components at Dirichlet nodes into local system */
#ifdef LOCALSYSTEMS_ST
   if (ele->locsys == locsys_yes)
   {
     locsys_trans_equant_dirich(ele, estif_global, NULL, NULL, intforce);
   }
#endif
break;/*----------------------------------------------------------------*/
/*-------------------------- calculate linear stiffness and mass matrix */
case calc_struct_linstiffmass:
      actmat = &(mat[ele->mat-1]);
      w1static_ke(ele,&actdata,actmat,estif_global,emass_global,intforce,0);
break;/*----------------------------------------------------------------*/
/*----------------------- calculate nonlinear stiffness and mass matrix */
case calc_struct_nlnstiffmass:
   actmat = &(mat[ele->mat-1]);
   if(ele->e.w1->kintype==total_lagr)
   {
     w1static_keug(ele,&actdata,actmat,estif_global,emass_global,intforce,0,
                   container);
   }
   else if(ele->e.w1->kintype==geo_lin)
   {
     w1static_ke(ele,&actdata,actmat,estif_global,emass_global,intforce,0);
   }
   if (intforce && container->isdyn && ele->proc == actintra->intra_rank)
      solserv_sol_localassemble(actintra,ele,intforce,1,2);
   /* rotate components at Dirichlet nodes into local system */
#ifdef LOCALSYSTEMS_ST
   if (ele->locsys == locsys_yes)
   {
     locsys_trans_equant_dirich(ele, estif_global, emass_global, NULL, intforce);
   }
#endif
break;/*----------------------------------------------------------------*/
/*-------------------------------- calculate stresses in a certain step */
case calc_struct_stress:
#ifdef D_MLSTRUCT
   if (genprob.multisc_struct == 1)
   {
     actmat = &(mat[ele->mat-1]);
     s2ml_static_ke(ele,&actdata,actmat,estif_global,NULL,intforce,3);
   }
   else
#endif /* D_MLSTRUCT */
   {
     imyrank = actintra->intra_rank;
     if (imyrank==ele->proc)
     {
      actmat = &(mat[ele->mat-1]);
      w1_cal_stress(ele,&actdata,actmat,0,0);
     }
   }/* end of else: if (genprob.multisc_struct == 1) */
break;/*----------------------------------------------------------------*/
/*------------------------------ calculate load vector of element loads */
case calc_struct_eleload:
   imyrank = actintra->intra_rank;
   actmat = &(mat[ele->mat-1]);
   w1_eleload(ele,&actdata,intforce,0,imyrank);
   /* rotate components at Dirichlet nodes into local system */
#ifdef LOCALSYSTEMS_ST
   if (ele->locsys == locsys_yes)
   {
     locsys_trans_equant_dirich(ele, NULL, NULL, NULL, intforce);
   }
#endif
break;/*----------------------------------------------------------------*/
case calc_struct_fsiload:
   imyrank = actintra->intra_rank;
   actmat = &(mat[ele->mat-1]);
   w1_fsiload(ele,&actdata,intforce,0,ipos,imyrank);
break;/*----------------------------------------------------------------*/
/*------------------------------ iterative update of internal variables */
case calc_struct_update_iterstep:
/* do nothing */
/* The iterative update of internal variables is used in the SOLID3
 * element (May 2007). This action is occurs here only for compatiblity. */
break;/*----------------------------------------------------------------*/
/*--------------------------------------- update after incremental step */
case calc_struct_update_istep:
#ifdef D_MLSTRUCT
   if (genprob.multisc_struct == 1)
   {
     actmat = &(mat[ele->mat-1]);
     s2ml_static_ke(ele,&actdata,actmat,estif_global,NULL,intforce,2);
   }
   else
#endif /* D_MLSTRUCT */
   {
     actmat = &(mat[ele->mat-1]);
     if(ele->e.w1->kintype == geo_lin)
     {
       w1static_ke(ele,&actdata,actmat,estif_global,NULL,intforce,2);
     }
     else
     {
       w1static_keug(ele,&actdata,actmat,estif_global,NULL,intforce,2,
                     container);
     }
     /* calculate mass matrix and kinetic energy (m.gee 4/03) */
     /* only mass matrix needed, it would be efficient to have a separate */
     /* routine w1static_mass */
     if (container->isdyn && ele->proc == actintra->intra_rank)
     {
       actmat = &(mat[ele->mat-1]);
       w1static_keug(ele,&actdata,actmat,estif_global,emass_global,NULL,0,
                     container);
       dyn_ekin_local(ele,emass_global,container);
     }
#ifdef GEMM
     if (container->isdyn) w1_update_history(ele,&actdata,actmat);
#endif
   }/* end of else: if (genprob.multisc_struct == 1) */
break;/*----------------------------------------------------------------*/
/*------------------------------------------------------- write restart */
case write_restart:
   actmat = &(mat[ele->mat-1]);
   w1_write_restart(ele,actmat,container->handsize,container->handles,0);
break;/*----------------------------------------------------------------*/
/*------------------------------------------------------- write restart */
case  read_restart:
   actmat = &(mat[ele->mat-1]);
   w1_read_restart(ele,actmat,container->handles,0);
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
#ifdef D_SSI
/*------- calculate coupling forces for ssi problems, firl / genk 10/03 */
case calc_struct_ssi_coup_force:
   for(i=0; i<ele->numnp; i++)
   {
     if (ele->node[i]->gnode->ssicouple == NULL) continue;
     if (ele->node[i]->gnode->ssicouple->ssi_couptyp == ssi_slave)
     {
       actmat = &(mat[ele->mat-1]);
       w1static_keug(ele,&actdata,actmat,estif_global,emass_global,intforce,0,
                     container);
       break;
     }
   }
break;/*----------------------------------------------------------------*/
#endif
/*---------------------------------------------------- finalise element */
case calc_struct_final:
   /* In essence, dynamically allocated variables are deallocated */
   w1_final(actpart, mat);
   w1static_ke(NULL,NULL,NULL,NULL,NULL,NULL,-1);
   w1static_keug(NULL,NULL,NULL,NULL,NULL,NULL,-1,NULL);
   w1_cal_stress(NULL,NULL,NULL,0,-1);
   w1_eleload(NULL,NULL,NULL,-1,-1);
   w1_write_restart(NULL,NULL,0,NULL,-1);
   w1_read_restart(NULL,NULL,NULL,-1);
break;
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
