#include "../headers/standardtypes.h"
#include "../headers/solution.h"
int cmp_int(const void *a, const void *b );
double cmp_double(const void *a, const void *b );


/*----------------------------------------------------------------------*
 |  calculate the mask of an vbr matrix                  m.gee 5/02     |
 *----------------------------------------------------------------------*/
void mask_vbr(FIELD         *actfield, 
              PARTITION     *actpart, 
              SOLVAR        *actsolv,
              INTRA         *actintra, 
              AZ_ARRAY_VBR  *vbr)
{
int       i,j,k,l;
int       numeq;
#ifdef DEBUG 
dstrc_enter("mask_vbr");
#endif
/*----------------------------------------------------------------------*/
/*------------------------------------------- put total size of problem */
vbr->numeq_total = actfield->dis[0].numeq;
/* count number of eqns on proc and build processor-global couplingdof 
                                                                 matrix */
msr_numeq(actfield,actpart,actsolv,actintra,&numeq);
vbr->numeq = numeq;
/*---------------------------------------------- allocate vector update */
amdef("update",&(vbr->update),numeq,1,"IV");
amzero(&(vbr->update));
/*--------------------------------put dofs in update in ascending order */
vbr_update(actfield,actpart,actsolv,actintra,vbr);
/*--------------------------------------------- make the block topology */
vbr_block_connect(actfield,actpart,actsolv,actintra,vbr);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mask_vbr */


/*----------------------------------------------------------------------*
 |  make the block connectivity                             m.gee 5/02  |
 *----------------------------------------------------------------------*/
int vbr_block_connect(FIELD         *actfield, 
                      PARTITION     *actpart, 
                      SOLVAR        *actsolv,
                      INTRA         *actintra,
                      AZ_ARRAY_VBR  *vbr)
{
int       i,j,k,l;
int       counter;
int       imyrank;
int       inprocs;
int       nodestack[300];
int       stacksize;
NODE     *actnode;
ELEMENT  *actele;
int       actnodeId;

int      *bupdate;
int      *rpntr;

#ifdef DEBUG 
dstrc_enter("vbr_block_connect");
#endif
/*----------------------------------------------------------------------*/
imyrank = actintra->intra_rank;
inprocs = actintra->intra_nprocs;
/*----------------------------------------------------------------------*/
vbr->block_consize = actfield->dis[0].numnp;
vbr->block_connect = (int**)CALLOC(vbr->block_consize,sizeof(int*));
if (!(vbr->block_connect)) dserror("Allocation of memory failed");
/*--------------------------------------------- loop nodes on partition */
for (i=0; i<actfield->dis[0].numnp; i++)
{
   /* make a node patch */
   actnode   = actpart->pdis[0].node[i];
   actnodeId = actpart->pdis[0].node[i]->Id_loc;
   counter=0;
   for (j=0; j<actnode->numele; j++)
   {
      actele = actnode->element[j];
      for (k=0; k<actele->numnp; k++)
      {
         nodestack[counter++] = actele->node[k]->Id_loc;
      }
   }
   /* delete doubles on patch */
   stacksize=counter;
   for (j=0; j<stacksize; j++)
   {
      k = nodestack[j];
      if (nodestack[j]==actnodeId) nodestack[j]=-1;
      if (nodestack[j]==-1) continue;
      for (l=j+1; l<stacksize; l++)
         if (nodestack[l] == k ||
             nodestack[l] == actnodeId)
             nodestack[l]=-1;
   }
   /* count true stacksize */
   counter=0;
   for (j=0; j<stacksize; j++)
      if (nodestack[j] != -1) counter++;
   /* allocate vector of block connectivity */
   /* block_connect[actnodeId][0] = my own block size */
   /* block_connect[actnodeId][1] = number off-diagonal blocks in this row */
   /* block_connect[actnodeId][2..] = id_loc of off-diagonal blocks */
   vbr->block_connect[actnodeId] = (int*)CALLOC(2+counter,sizeof(int));
   if (!(vbr->block_connect[actnodeId]))
   dserror("Allocation of memory failed");
   vbr->block_connect[actnodeId][1] = counter;
   /* move the off-diagonals to block_connect */
   counter=0;
   for (j=0; j<stacksize; j++)
      if (nodestack[j] != -1)
         vbr->block_connect[actnodeId][2+counter++] = nodestack[j];
   /* sort them */
   qsort((int*)(&(vbr->block_connect[actnodeId][2])),
                vbr->block_connect[actnodeId][1],
                sizeof(int),cmp_int);
   /* calc my own block size */
   counter=0;
   for (j=0; j<actnode->numdf; j++)
   if (actnode->dof[j] < vbr->numeq_total) 
   counter++;
   vbr->block_connect[actnodeId][0] = counter;
}/* end of for (i=0; i<actpart->pdis[0].numnp; i++)*/
/*---------------------- make list of nodal blocks updated on this proc */
bupdate = amdef("bupdate",&(vbr->bupdate),actpart->pdis[0].numnp,1,"IV");
for (i=0; i<actpart->pdis[0].numnp; i++)
   bupdate[i] = actpart->pdis[0].node[i]->Id_loc;
/*--------------------------- now loop only my own nodes and make rpntr */
rpntr    = amdef("rpntr",&(vbr->rpntr),actpart->pdis[0].numnp+1,1,"IV");
counter  = 1;
rpntr[0] = 0;
for (i=0; i<actpart->pdis[0].numnp; i++)
{
   actnodeId = actpart->pdis[0].node[i]->Id_loc;
   rpntr[counter] = vbr->block_connect[actnodeId][0] + rpntr[counter-1];
   counter++;
}
/*-------------------------- now loop all nodes and make column pointer */





/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of vbr_block_connect */



/*----------------------------------------------------------------------*
 |  allocate update put dofs in update in ascending order   m.gee 5/01  |
 *----------------------------------------------------------------------*/
int vbr_update(FIELD         *actfield, 
                PARTITION     *actpart, 
                SOLVAR        *actsolv,
                INTRA         *actintra,
                AZ_ARRAY_VBR  *vbr)
{
int       i,j,k,l;
int       counter;
int      *update;
int       dof;
int       foundit;
int       imyrank;
int       inprocs;
NODE     *actnode;
ARRAY     coupledofs;
#ifdef DEBUG 
dstrc_enter("vbr_update");
#endif
/*----------------------------------------------------------------------*/
imyrank = actintra->intra_rank;
inprocs = actintra->intra_nprocs;
/*------------------ make a local copy of the array actpart->coupledofs */
am_alloc_copy(&(actpart->pdis[0].coupledofs),&coupledofs);
/*----------------------------------------------------------------------*/
update = vbr->update.a.iv;
counter=0;
/*------------------------------------- loop the nodes on the partition */
for (i=0; i<actpart->pdis[0].numnp; i++)
{
   actnode = actpart->pdis[0].node[i];
   for (l=0; l<actnode->numdf; l++)
   {
      dof = actnode->dof[l];
      /* dirichlet condition on dof */
      if (dof >= actfield->dis[0].numeq) continue;
      /* no coupling on dof */
      if (actnode->gnode->couple==NULL)
      {
         update[counter] = dof;
         counter++;
         continue;
      }
      else /* coupling on node */
      {
         foundit=0;
         /* find dof in coupledofs */
         for (k=0; k<coupledofs.fdim; k++)
         {
            if (dof == coupledofs.a.ia[k][0])
            {
               /* am I owner of this dof or not */
               if (coupledofs.a.ia[k][imyrank+1]==2) 
               foundit=2;
               else if (coupledofs.a.ia[k][imyrank+1]==1)                                     
               foundit=1;
               break;
            }
         }
         /* dof found in coupledofs */
         if (foundit==2)/* I am master owner of this coupled dof */
         {
            update[counter] = dof;
            counter++;
            coupledofs.a.ia[k][imyrank+1]=1;
            continue;
         }
         else if (foundit==1)/* I am slave owner of this coupled dof */
         {
           /* do nothing, this dof doesn't exist for me (no more)*/
         }
         else /* this dof is not a coupled one */
         {
            update[counter] = dof;
            counter++;
            continue;
         }
      }
      
   }
}
/*----------- check whether the correct number of dofs has been counted */
if (counter != vbr->numeq) dserror("Number of dofs in VBR-vector update wrong");
/*---------------------------- sort the vector update just to make sure */
qsort((int*) update, counter, sizeof(int), cmp_int);
/*----------------------------------------------------------------------*/
amdel(&coupledofs);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of vbr_update */
