#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
int cmp_int(const void *a, const void *b );
double cmp_double(const void *a, const void *b );


/*----------------------------------------------------------------------*
 |  calculate the mask of an rc_ptr matrix               m.gee 1/02     |
 *----------------------------------------------------------------------*/
void mask_skyline(FIELD         *actfield, 
                  PARTITION     *actpart, 
                  SOLVAR        *actsolv,
                  INTRA         *actintra, 
                  SKYMATRIX     *sky)
{
int       i;
int       numeq;
int     **dof_connect;
ARRAY     red_dof_connect;
#ifdef DEBUG 
dstrc_enter("mask_skyline");
#endif
/*----------------------------------------------------------------------*/
/* remember some facts:
   PARTITION is different on every proc.
   AZ_ARRAY_MSR will be different on every proc
   FIELD is the same everywhere
   In this routine, the vectors update and bindx and val are determined
   in size and allocated, the contents of the vectors update and bindx 
   are calculated
*/
/*------------------------------------------- put total size of problem */
sky->numeq_total = actfield->dis[0].numeq;
/* count number of eqns on proc and build processor-global couplingdof 
                                                                 matrix */
msr_numeq(actfield,actpart,actsolv,actintra,&numeq);
sky->numeq = numeq;
/*---------------------------------------------- allocate vector update */
amdef("update",&(sky->update),numeq,1,"IV");
amzero(&(sky->update));
/*--------------------------------put dofs in update in ascending order */
skyline_update(actfield,actpart,actsolv,actintra,sky);
/*------------------------ count number of nonzero entries on partition 
                                    and calculate dof connectivity list */
   /*
      dof_connect[i][0] = lenght of dof_connect[i]
      dof_connect[i][1] = iscoupled ( 1 or 2 ) 
      dof_connect[i][2] = dof
      dof_connect[i][ 2..dof_connect[i][0]-1 ] = connected dofs exluding itself 
   */
dof_connect = (int**)CCACALLOC(sky->numeq_total,sizeof(int*));
if (!dof_connect) dserror("Allocation of dof_connect failed");
skyline_nnz_topology(actfield,actpart,actsolv,actintra,sky,dof_connect);
/*------------------------------------------------------ make nnz_total */
sky->nnz_total=sky->nnz;
/*------------------------------make dof_connect redundant on all procs */
skyline_make_red_dof_connect(actfield,
                             actpart,
                             actsolv,
                             actintra,
                             sky,
                             dof_connect,
                             &red_dof_connect);
/*---------------------------------------- make arrays from dof_connect */
skyline_make_sparsity(sky,&red_dof_connect);
/*---------------------------------------- delete the array dof_connect */
for (i=0; i<sky->numeq_total; i++)
{
   if (dof_connect[i]) CCAFREE(dof_connect[i]);
}
CCAFREE(dof_connect);
amdel(&red_dof_connect);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mask_skyline */



/*----------------------------------------------------------------------*
 |  allocate update put dofs in update in ascending order   m.gee 1/02  |
 *----------------------------------------------------------------------*/
void  skyline_update(FIELD         *actfield, 
                    PARTITION     *actpart, 
                    SOLVAR        *actsolv,
                    INTRA         *actintra,
                    SKYMATRIX     *sky)
{
int       i,k,l;
int       counter;
int      *update;
int       dof;
int       foundit;
int       imyrank;
int       inprocs;
NODE     *actnode;
ARRAY     coupledofs;
#ifdef DEBUG 
dstrc_enter("skyline_update");
#endif
/*----------------------------------------------------------------------*/
imyrank = actintra->intra_rank;
inprocs = actintra->intra_nprocs;
/*------------------ make a local copy of the array actpart->coupledofs */
am_alloc_copy(&(actpart->pdis[0].coupledofs),&coupledofs);
/*------------------------------------- loop the nodes on the partition */
update = sky->update.a.iv;
counter=0;
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
/*---------- check whether the correct number of dofs have been counted */
if (counter != sky->numeq) dserror("Number of dofs in update wrong");
/*---------------------------- sort the vector update just to make sure */
qsort((int*) update, counter, sizeof(int), cmp_int);
/*----------------------------------------------------------------------*/
amdel(&coupledofs);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of skyline_update */



/*----------------------------------------------------------------------*
 |  calculate number of nonzero entries and dof topology    m.gee 1/02  |
 *----------------------------------------------------------------------*/
void  skyline_nnz_topology(FIELD      *actfield, 
                         PARTITION    *actpart, 
                         SOLVAR       *actsolv,
                         INTRA        *actintra,
                         SKYMATRIX    *sky,
                         int         **dof_connect)
{
int        i,j,k,l,m,n;
int        counter,counter2;
int        dof;
int        nnz;
int        iscoupled;
int       *update;
int        numeq;
int        actdof;
int        dofmaster;
int        dofslave;
int        sendlenght,recvlenght;
int        recvflag;
NODE      *centernode;
NODE      *actnode;
ELEMENT   *actele;
ARRAY      dofpatch;
ARRAY     *coupledofs;
int        imyrank;
int        inprocs;
int        actndis;

#ifdef PARALLEL 
MPI_Status status;
#endif

#ifdef DEBUG 
dstrc_enter("skyline_nnz_topology");
#endif
/*----------------------------------------------------------------------*/
imyrank = actintra->intra_rank;
inprocs = actintra->intra_nprocs;
/*------------------------------------------- set actual discretisation */
actndis=0;
/*----------------------------------------------------------- shortcuts */
sky->nnz=0;
/*numeq  = sky->numeq;*/
numeq  = sky->numeq_total;
update = sky->update.a.iv;
for (i=0; i<sky->numeq_total; i++) dof_connect[i]=NULL;
amdef("tmp",&dofpatch,1000,1,"IV");
amzero(&dofpatch);
/*----------------------------------------------------------------------*/
for (i=0; i<numeq; i++)
{
   dof = i;
   /*dof = update[i];*/
   /*------------------------------ check whether this is a coupled dof */
   iscoupled=0;
   dof_in_coupledofs(dof,actpart,&iscoupled);
   if (iscoupled==1) continue;
   /*--------------------------------- find the centernode for this dof */
   centernode=NULL;
   for (j=0; j<actfield->dis[actndis].numnp; j++)
   {
      for (k=0; k<actfield->dis[actndis].node[j].numdf; k++)
      {
         if (actfield->dis[actndis].node[j].dof[k] == dof)
         {
            centernode = &actfield->dis[actndis].node[j];
          goto nodefound1;
         }
      }
   }
   nodefound1:
   dsassert(centernode!=NULL,"Cannot find centernode of a patch");
   /*--------------------------------- make dof patch around centernode */
   counter=0;
   for (j=0; j<centernode->numele; j++)
   {
      actele = centernode->element[j];
      for (k=0; k<actele->numnp; k++)
      {
         actnode = actele->node[k];
         for (l=0; l<actnode->numdf; l++)
         {
            if (actnode->dof[l] < actfield->dis[0].numeq)
            {
               if (counter>=dofpatch.fdim) amredef(&dofpatch,dofpatch.fdim+500,1,"IV");
               dofpatch.a.iv[counter] = actnode->dof[l];
               counter++;
            }
         }
      }
   }
   /*----------------------------------------- delete doubles on patch */
   /*------------------------------- also delete dof itself from patch */
   for (j=0; j<counter; j++)
   {
      actdof = dofpatch.a.iv[j];
      if (dofpatch.a.iv[j]==dof) dofpatch.a.iv[j]=-1;
      if (actdof==-1) continue;
      for (k=j+1; k<counter; k++)
      {
         if (dofpatch.a.iv[k] == actdof ||
             dofpatch.a.iv[k] == dof      ) dofpatch.a.iv[k]=-1;
      }
   }
   /*----------------------------------- count number of dofs on patch */
   counter2=0;
   for (j=0; j<counter; j++)
   {
      if (dofpatch.a.iv[j] != -1) counter2++;
   }
   /*-------------- allocate the dof_connect vector and put dofs in it */
   dof_connect[dof] = (int*)CCACALLOC(counter2+3,sizeof(int));
   if (!dof_connect[dof]) dserror("Allocation of dof connect list failed");
   dof_connect[dof][0] = counter2+3;
   dof_connect[dof][1] = 0; 
   dof_connect[dof][2] = dof;
   /*
      dof_connect[i][0] = lenght of dof_connect[i]
      dof_connect[i][1] = iscoupled ( 1 or 2 ) done later on 
      dof_connect[i][2] = dof
      dof_connect[i][ 2..dof_connect[i][0]-1 ] = connected dofs exluding itself 
   */
   counter2=0;
   for (j=0; j<counter; j++)
   {
      if (dofpatch.a.iv[j] != -1) 
      {
         dof_connect[dof][counter2+3] = dofpatch.a.iv[j];
         counter2++;
      }
   }
}  /* end of loop over numeq */ 
/*--------------------------------------------- now do the coupled dofs */
coupledofs = &(actpart->pdis[0].coupledofs);
for (i=0; i<coupledofs->fdim; i++)
{
   dof = coupledofs->a.ia[i][0];
   /*------------------------------------- find all patches to this dof */
   counter=0;
   for (j=0; j<actfield->dis[actndis].numnp; j++)
   {
      centernode=NULL;
      /* */
      for (l=0; l<actfield->dis[actndis].node[j].numdf; l++)
      {
         if (dof == actfield->dis[actndis].node[j].dof[l])
         {
            centernode = &actfield->dis[actndis].node[j];
            break;
         }
      }
      if (centernode !=NULL)
      {
         /*--------------------------- make dof patch around centernode */
         for (k=0; k<centernode->numele; k++)
         {
            actele = centernode->element[k];
            for (m=0; m<actele->numnp; m++)
            {
               actnode = actele->node[m];
               for (l=0; l<actnode->numdf; l++)
               {
                  if (actnode->dof[l] < actfield->dis[0].numeq)
                  {
                     if (counter>=dofpatch.fdim) amredef(&dofpatch,dofpatch.fdim+500,1,"IV");
                     dofpatch.a.iv[counter] = actnode->dof[l];
                     counter++;
                  }
               }
            }
         }
      }
   }/* end of making dofpatch */
   /*----------------------------------------- delete doubles on patch */
   for (j=0; j<counter; j++)
   {
      actdof = dofpatch.a.iv[j];
      if (actdof==-1) continue;
      if (actdof==dof) dofpatch.a.iv[j]=-1;
      for (k=j+1; k<counter; k++)
      {
         if (dofpatch.a.iv[k] == actdof ||
             dofpatch.a.iv[k] == dof      ) dofpatch.a.iv[k]=-1;
      }
   }
   /*----------------------------------- count number of dofs on patch */
   counter2=0;
   for (j=0; j<counter; j++)
   {
      if (dofpatch.a.iv[j] != -1) counter2++;
   }
   /*-------------- allocate the dof_connect vector and put dofs in it */
   dof_connect[dof] = (int*)CCACALLOC(counter2+3,sizeof(int));
   if (!dof_connect[dof]) dserror("Allocation of dof connect list failed");
   dof_connect[dof][0] = counter2+3;
   dof_connect[dof][1] = 0;
   dof_connect[dof][2] = dof;
   /*-------------------------- put the patch to the dof_connect array */
   counter2=0;
   for (j=0; j<counter; j++)
   {
      if (dofpatch.a.iv[j] != -1) 
      {
         dof_connect[dof][counter2+3] = dofpatch.a.iv[j];
         counter2++;
      }
   }
} /* end of loop over coupled dofs */
/*---------------------------- now go through dof_connect and count nnz */
nnz=0;
for (i=0; i<numeq; i++)
{
   nnz += (dof_connect[i][0]-2);
}
sky->nnz=nnz;
/*--------- last thing to do is to order dof_connect in ascending order */
for (i=0; i<numeq; i++)
{
   qsort((int*)(&(dof_connect[i][3])), dof_connect[i][0]-3, sizeof(int), cmp_int);
}
/*----------------------------------------------------------------------*/
amdel(&dofpatch);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of skyline_nnz_topology */


/*----------------------------------------------------------------------*
 |  make the dof_connect list redundant                    m.gee 01/02  |
 *----------------------------------------------------------------------*/
void   skyline_make_red_dof_connect(FIELD         *actfield, 
                                   PARTITION     *actpart, 
                                   SOLVAR        *actsolv,
                                   INTRA         *actintra,
                                   SKYMATRIX     *sky,
                                   int          **dof_connect,
                                   ARRAY         *red_dof_connect)
{
int        i,j;

int        imyrank;
int        inprocs;

int      **reddof;

int        max_dof_connect;
/*----------------------------------------------------------------------*/
/* communicate coupled dofs */
/*----------------------------------------------------------------------*/

#ifdef DEBUG 
dstrc_enter("skyline_make_red_dof_connect");
#endif
/*----------------------------- check for largest row in my dof_connect */
max_dof_connect=0;
for (i=0; i<sky->numeq_total; i++)
{
   if (dof_connect[i])
   if (dof_connect[i][0]>max_dof_connect)
   max_dof_connect=dof_connect[i][0];
}
/*---------------- allocate temporary array to hold global connectivity */
reddof = amdef("tmp",red_dof_connect,sky->numeq_total,max_dof_connect,"IA");
/*----------------------------- put my own dof_connect values to reddof */
for (i=0; i<sky->numeq_total; i++)
{
   if (dof_connect[i])
   {
      for (j=0; j<dof_connect[i][0]; j++) 
         reddof[i][j] = dof_connect[i][j];
   }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of skyline_make_red_dof_connect */







/*----------------------------------------------------------------------*
 |  make sparsity mask for skyline matrix                   m.gee 1/02  |
 *----------------------------------------------------------------------*/
void  skyline_make_sparsity(SKYMATRIX  *sky, ARRAY *red_dof_connect)
{
int        i,j,k,l;
int      **reddof;
int       *maxa;
int        actdof;
int        lenght;
int        counter=0;
int        mindof;
#ifdef DEBUG 
dstrc_enter("skyline_make_sparsity");
#endif
/*----------------------------------------------------------------------*/
   /*
      reddof[i][0] = lenght of reddof[i]
      reddof[i][1] = iscoupled ( 1 or 2 ) done later on 
      reddof[i][2] = dof
      reddof[i][ 2..reddof[i][0]-1 ] = connected dofs exluding itself 
   */
reddof = red_dof_connect->a.ia;
/*------------------------------------------------------- allocate maxa */
maxa = amdef("maxa",&(sky->maxa),sky->numeq_total+1,1,"IV");
/*------------------------------------------------------- loop the dofs */       
for (i=0; i<sky->numeq_total; i++)
{
   actdof = reddof[i][2];
   if (actdof != i) printf("Warning: skyline format mixed up!\n");
   
   /*-------------- search reddof[i][2..reddof[i][0] ] for smallest dof */
   mindof=10000000;
   for (j=2; j<reddof[i][0]; j++)
   {
      if (mindof>reddof[i][j]) 
      mindof = reddof[i][j];
   }
   /*-------------------------- check distance between actdof and mindof */
   lenght = actdof-mindof+1;
   /*--------------------------- maxa[i] holds start of column of actdof */
   maxa[i]=counter;
   counter+=lenght;
}
maxa[i]=counter;
/*-----------------------------------------------------------allocate A */
amdef("A",&(sky->A),maxa[i],1,"DV");
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of skyline_make_sparsity */
