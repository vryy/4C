#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
/*----------------------------------------------------------------------*
 |  routine to call rhs-routines                         m.gee 10/01    |
 |  in here, only step dependent loads are calculated, which means      |
 |  neumann conditions which may change from time/load step to step.    |
 *----------------------------------------------------------------------*/
void calrhs(FIELD        *actfield,     /* the active field */
            SOLVAR       *actsolv,      /* the active SOLVAR */
            PARTITION    *actpart,      /* my partition of this field */
            INTRA        *actintra,     /* the field's intra-communicator */
            int           actsysarray,  /* the active sparse array */
            DIST_VECTOR  *rhs1,         /* 1 dist. vectors for rhs */
            int           kstep,
            CALC_ACTION  *action)       /* action to be passed to element routines */
{
int i;
static ARRAY rhs_a;
static double *rhs;
#ifdef PARALLEL 
static ARRAY rhsrecv_a;
static double *rhsrecv;
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
/*--------- inherit the neuman conditions from design to discretization */
for (i=0; i<actfield->ndis; i++) inherit_design_dis_neum(&(actfield->dis[i]));
/*---------------------------------- calculate point neumann conditions */
rhs_point_neum(rhs,rhs1->numeq_total,actpart);
/*--- line/surface/volume loads are a matter of element integration and */
/*                                            different in each element */
*action = calc_struct_eleload;
calelm(actfield,actsolv,actpart,actintra,actsysarray,-1,rhs,NULL,rhs1->numeq_total,kstep,action);
/*-------------------------------------------- allreduce the vector rhs */
#ifdef PARALLEL 
MPI_Allreduce(rhs,rhsrecv,rhs1->numeq_total,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
/*----------------- assemble rhs to rhs1, which is a distributed vector */
assemble_vec(actintra,sysarraytyp,sysarray,rhs1,rhsrecv,1.0);
#else
assemble_vec(actintra,sysarraytyp,sysarray,rhs1,rhs,1.0);
#endif
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
void rhs_point_neum(double *rhs, int dimrhs, PARTITION *actpart)     
{
int             i,j;
int             dof;
NODE           *actnode;
NEUM_CONDITION *actneum;
#ifdef DEBUG
dstrc_enter("rhs_point_neum");
#endif
/*----------------------------------------------------------------------*/
for (i=0; i<actpart->pdis[0].numnp; i++)
{
   /*------------------------ check presence of nodal eumann conditions */
   if (actpart->pdis[0].node[i]->gnode->neum == NULL) continue;
   /*-------------------------------------------------- set active node */
   actnode = actpart->pdis[0].node[i];
   /*------------------------------------- set active neumann condition */
   actneum = actnode->gnode->neum;
   /*------------------------------------------------ loop dofs of node */
   for (j=0; j<actnode->numdf; j++)
   {
      /*-------------------------- check flag whether a dof has a value */
      if (actneum->neum_onoff.a.iv[j]==0) continue;
      /*----------------------------------- get dof which has the value */
      dof = actnode->dof[j];
      /*---------- check whether this dof is inside free dofs (<dimrhs) */
      if (dof >= dimrhs) continue;
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

 


