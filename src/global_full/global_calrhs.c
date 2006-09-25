/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../solver/solver.h"

#ifdef D_SSI
#include "../ssi_full/ssi_prototypes.h"
#endif

#ifdef D_SHELL9
#include "../shell9/shell9.h"
#endif /*D_SHELL9*/

/*----------------------------------------------------------------------*
 |  routine to call rhs-routines                         m.gee 10/01    |
 |  in here, only step dependent loads are calculated, which means      |
 |  neumann conditions which may change from time/load step to step.    |
 *----------------------------------------------------------------------*/
void calrhs(FIELD        *actfield,     /* the active field */
            SOLVAR       *actsolv,      /* the active SOLVAR */
            PARTITION    *actpart,      /* my partition of this field */
            INTRA        *actintra,     /* the field's intra-communicator */
            INT           actsysarray,  /* the active sparse array */
            DIST_VECTOR  *rhs1,         /* 1 dist. vectors for rhs */
            CALC_ACTION  *action,       /* action to be passed to element routines */
            CONTAINER    *container)     /*!< contains variables defined in container.h */
{
INT i;
static ARRAY rhs_a;
static DOUBLE *rhs;
#ifdef PARALLEL
static ARRAY rhsrecv_a;
static DOUBLE *rhsrecv;
#endif
SPARSE_TYP   *sysarraytyp;
SPARSE_ARRAY *sysarray;
#ifdef DEBUG
dstrc_enter("calrhs");
#endif
/*----------------------------------------------------------------------*/
sysarraytyp = &(actsolv->sysarray_typ[actsysarray]);
sysarray    = &(actsolv->sysarray[actsysarray]);
/*-------------------- create a temporary vector of redundant full size */
if (rhs_a.Typ != cca_DV)
{
   rhs = amdef("tmprhs",&rhs_a,rhs1->numeq_total,1,"DV");
}
if (rhs_a.fdim < rhs1->numeq_total)
{
   amdel(&rhs_a);
   rhs = amdef("tmprhs",&rhs_a,rhs1->numeq_total,1,"DV");
}
amzero(&rhs_a);
#ifdef PARALLEL
if (rhsrecv_a.Typ != cca_DV)
{
   rhsrecv = amdef("tmprhs",&rhsrecv_a,rhs1->numeq_total,1,"DV");
}
if (rhsrecv_a.fdim < rhs1->numeq_total)
{
   amdel(&rhsrecv_a);
   rhsrecv = amdef("tmprhs",&rhsrecv_a,rhs1->numeq_total,1,"DV");
}
amzero(&rhsrecv_a);
#endif
/*--------- inherit the Neumann conditions from design to discretization */
if (container->inherit>0)
for (i=0; i<actfield->ndis; i++) inherit_design_dis_neum(&(actfield->dis[i]));
/*---------------------------------- calculate point neumann conditions */
if (container->point_neum>0) rhs_point_neum(rhs,rhs1->numeq_total,actpart);
/*-------------------------------- assemble rhs for ssi coupling forces */
#ifdef D_SSI
if (*action == calc_struct_ssiload) 
{  
  ssiserv_rhs_point_neum(rhs,rhs1->numeq_total,actpart);
  *action = calc_struct_eleload;
}
#endif
#ifdef D_MORTAR
#ifdef D_FSI
if (*action == calc_struct_fsiload_mtr) 
{  
  fsiserv_rhs_point_neum(rhs,rhs1->numeq_total,actpart);
  *action = calc_struct_eleload;
}
#endif
#endif
/*--- line/surface/volume loads are a matter of element integration and */
/*                                            different in each element */
#if 0
/* changed: action is now set before calling calrhs!!!!                 */
*action = calc_struct_eleload;
#endif
container->dvec         = rhs;
container->dirich       = NULL;
container->global_numeq = rhs1->numeq_total;
calelm(actfield,actsolv,actpart,actintra,actsysarray,-1,container,action);
/*-------------------------------------------- allreduce the vector rhs */
#ifdef PARALLEL
MPI_Allreduce(rhs,rhsrecv,rhs1->numeq_total,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
/*----------------- assemble rhs to rhs1, which is a distributed vector */
assemble_vec(actintra,sysarraytyp,sysarray,rhs1,rhsrecv,1.0);
#else
assemble_vec(actintra,sysarraytyp,sysarray,rhs1,rhs,1.0);
#endif
/*
m.gee
testen, ob linienlasten ueber prozessorengrenzen hinweg richtig assembliert werden!
*/
/* genk 05/03
GLINEs haben nun einen Prozessor als eindeutigen Besitzer gline->proc.
Integrale ueber die Elementkanten werden nur vom besitzenden Prozessor
ausgewertet.
Umsetzung fuer
 - WALL1 - Linienlasten
 - WALL1 - FSI-Lasten
*/
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of calrhs */

/*----------------------------------------------------------------------*
 |  point neumann conditions                             m.gee 3/02     |
 | Attention! This assembly of nodal forces works only correct with     |
 | partitioning in the "Cut_Elements" style !!!!
 *----------------------------------------------------------------------*/
void rhs_point_neum(DOUBLE *rhs, INT dimrhs, PARTITION *actpart)
{
INT             i,j;
INT             dof;
NODE           *actnode;
NEUM_CONDITION *actneum;
#ifdef D_SHELL9
INT             nsurf   = 0;     /* 1=MID; 2=TOP; 3=BOT */
INT             numklay = 0;   /* number of kinematic layers if shell9 */
#endif
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("rhs_point_neum");
#endif
/*----------------------------------------------------------------------*/
for (i=0; i<actpart->pdis[0].numnp; i++)
{
   /*----------------------- check presence of nodal neumann conditions */
   if (actpart->pdis[0].node[i]->gnode->neum == NULL) continue;
   /*-------------------------------------------------- set active node */
   actnode = actpart->pdis[0].node[i];
   /*------------------------------------- set active neumann condition */
   actneum = actnode->gnode->neum;

#ifdef D_SHELL9
     if (actnode->element[0]->eltyp == el_shell9)
     {
       numklay = (actnode->numdf-3)/3;
       /* modify the loadvector if load is applied on the surface of shell elements */
       /* */
       switch(actneum->neum_surf)
       {
       case mid:
          nsurf = 1;
       break;
       case top:
          nsurf = 2;
       break;
       case bot:
          nsurf = 3;
       break;
       default:
          dserror("Unknown type of neum_surf");
       break;
       }/*end of switch(actneum->neum_surf*/

       /* modify the nodal load vector due to nsurf */
       s9_surf_P(actneum->neum_val.a.dv, nsurf, numklay);

       /*switch the neum_onoff to 1 due to nsurf */
       s9_surf_onoff(actneum->neum_onoff.a.iv, nsurf, numklay);
     }
#endif /*D_SHELL9*/

   /*------------------------------------------------ loop dofs of node */
   for (j=0; j<actnode->numdf; j++)
   {
      /*-------------------------- check flag whether a dof has a value */
      if (actneum->neum_onoff.a.iv[j]==0) continue;
      /*----------------------------------- get dof which has the value */
      dof = actnode->dof[j];
      /*---------- check whether this dof is inside free dofs (<dimrhs) */
#if defined(SOLVE_DIRICH) || defined(SOLVE_DIRICH2)
      if (actnode->gnode->dirich!=NULL &&
          actnode->gnode->dirich->dirich_onoff.a.iv[j]!=0)
        continue;
#else
      if (dof >= dimrhs) continue;
#endif
      /*-------------------------------------------- assemble the value */
      rhs[dof] += actneum->neum_val.a.dv[j];
   }
   /* unset the neumann condition from the discretization, to make sure */
   /* it is not assembled again */
   actnode->gnode->neum=NULL;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of rhs_point_neum */
