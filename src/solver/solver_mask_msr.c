#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
int cmp_int(const void *a, const void *b );
double cmp_double(const void *a, const void *b );

static int kk;
/*----------------------------------------------------------------------*
 |  calculate the mask of an msr matrix                  m.gee 5/01     |
 *----------------------------------------------------------------------*/
void mask_msr(FIELD         *actfield, 
              PARTITION     *actpart, 
              SOLVAR        *actsolv,
              INTRA         *actintra, 
              AZ_ARRAY_MSR  *msr,
	      int            actndis)
{
int       i,j,k,l;
int       numeq;
int     **dof_connect;
#ifdef DEBUG 
dstrc_enter("mask_msr");
#endif

/*------------------------------------------- set actual discretisation */
kk=actndis;
/*----------------------------------------------------------------------*/
/* remember some facts:
   PARTITION is different on every proc.
   AZ_ARRAY_MSR will be different on every proc
   FIELD is the same everywhere
   In this routine, the vectors update and bindx and val are determined
   in size and allocated, the contents of the vectors update and bindx 
   are calculated
/*------------------------------------------- put total size of problem */
msr->numeq_total = actfield->dis[kk].numeq;
/* count number of eqns on proc and build processor-global couplingdof 
                                                                 matrix */
msr_numeq(actfield,actpart,actsolv,actintra,&numeq);
msr->numeq = numeq;
/*---------------------------------------------- allocate vector update */
amdef("update",&(msr->update),numeq,1,"IV");
amzero(&(msr->update));
/*--------------------------------put dofs in update in ascending order */
msr_update(actfield,actpart,actsolv,actintra,msr);
/*------------------------ count number of nonzero entries on partition 
                                    and calculate dof connectivity list */
   /*
      dof_connect[i][0] = lenght of dof_connect[i] 
      dof_connect[i][1] = iscoupled ( 1 or 2 ) 
      dof_connect[i][2] = dof
      dof_connect[i][ 2..dof_connect[i][0]-1 ] = connected dofs exluding itself 
   */
dof_connect = (int**)CCACALLOC(msr->numeq_total,sizeof(int*));
if (!dof_connect) dserror("Allocation of dof_connect failed");
msr_nnz_topology(actfield,actpart,actsolv,actintra,msr,dof_connect);
/*---------------------------------------------- allocate bindx and val */
amdef("bindx",&(msr->bindx),(msr->nnz+1),1,"IV");
amdef("val"  ,&(msr->val)  ,(msr->nnz+1),1,"DV");
/*---------------------------------------------------------- make bindx */
msr_make_bindx(actfield,actpart,actsolv,msr,dof_connect);
/*---------------------------------------- delete the array dof_connect */
for (i=0; i<msr->numeq_total; i++)
{
   if (dof_connect[i]) CCAFREE(dof_connect[i]);
}
CCAFREE(dof_connect);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mask_msr */






/*----------------------------------------------------------------------*
 |  count processor local and global number of equations    m.gee 5/01  |
 *----------------------------------------------------------------------*/
void msr_numeq(FIELD         *actfield, 
               PARTITION    *actpart, 
               SOLVAR       *actsolv,
               INTRA        *actintra,
               int          *numeq)
{
int       i,j,k,l;
int       counter;
int       dof;
int       iscoupled;
int      *sendbuff,*recvbuff, sendsize;
int      *tmp;
int       inter_proc;
long int  min;
int       proc;
int       inprocs;
int       imyrank;
NODE     *actnode;
#ifdef DEBUG 
dstrc_enter("msr_numeq");
#endif
/*----------------------------------------------------------------------*/
imyrank = actintra->intra_rank;
inprocs = actintra->intra_nprocs;
/*------------------------- first make a list of dofs which are coupled */
/*----------------------------------- estimate size of coupdofs to 5000 */
amdef("coupledofs",&(actpart->pdis[kk].coupledofs),5000,1,"IV");
amzero(&(actpart->pdis[kk].coupledofs));
counter=0;
/*-------------------------------- loop all nodes and find coupled dofs */
for (i=0; i<actfield->dis[kk].numnp; i++)
{
   actnode = &(actfield->dis[0].node[i]);
   if (actnode->gnode->couple==NULL) continue;
   for (l=0; l<actnode->numdf; l++)
   {
      if (actnode->dof[l]>=actfield->dis[kk].numeq) continue;
      /* there is coupling on this dof */
      if (actnode->gnode->couple->couple.a.ia[l][0] != 0 ||
          actnode->gnode->couple->couple.a.ia[l][1] != 0 )
      {
         if (counter>=actpart->pdis[kk].coupledofs.fdim) 
         amredef(&(actpart->pdis[kk].coupledofs),(actpart->pdis[kk].coupledofs.fdim+5000),1,"IV");
         /* the coupled dof could be dirichlet conditioned */
         if (actnode->dof[l]<actfield->dis[kk].numeq)
         {
            actpart->pdis[kk].coupledofs.a.iv[counter] = actnode->dof[l];
            counter++;
         }
      }
   }
}
amredef(&(actpart->pdis[kk].coupledofs),counter,1,"IV");
/*---------------------------------- delete the doubles in coupledofs */
for (i=0; i<actpart->pdis[kk].coupledofs.fdim; i++)
{
   if (actpart->pdis[kk].coupledofs.a.iv[i]==-1) continue;
   dof = actpart->pdis[kk].coupledofs.a.iv[i];
   for (j=i+1; j<actpart->pdis[kk].coupledofs.fdim; j++)
   {
      if (actpart->pdis[kk].coupledofs.a.iv[j]==dof) actpart->pdis[kk].coupledofs.a.iv[j]=-1;
   }
}
/*--------- move all remaining coupdofs to the front and redefine again */
counter=0;
for (i=0; i<actpart->pdis[kk].coupledofs.fdim; i++)
{
   if (actpart->pdis[kk].coupledofs.a.iv[i]!=-1)
   {
      actpart->pdis[kk].coupledofs.a.iv[counter] = actpart->pdis[kk].coupledofs.a.iv[i];
      counter++;
   }
}
amredef(&(actpart->pdis[kk].coupledofs),counter,inprocs+1,"IA");
/*------------------- the newly allocated columns have to be initialized */
for (i=1; i<actpart->pdis[kk].coupledofs.sdim; i++)
for (j=0; j<actpart->pdis[kk].coupledofs.fdim; j++) 
actpart->pdis[kk].coupledofs.a.ia[j][i]=0;

/* processor looks on his own domain whether he has some of these coupdofs, 
   puts this information in the array coupledofs in the column myrank+1, so it 
   can be allreduced 

   The matrix has the following style (after allreduce on all procs the same):
   
               ----------------------
               | 12 | 1 | 0 | 1 | 0 |
               | 40 | 1 | 0 | 0 | 0 |
               | 41 | 1 | 1 | 1 | 1 |
               | 76 | 0 | 1 | 1 | 0 |
               ----------------------
               
               column 0                : number of the coupled equation
               column 1 - inprocs+1 : proc has coupled equation or not
               
*/
if (inprocs==1) /*--------------------------------- sequentiell version */
{
   for (k=0; k<actpart->pdis[kk].coupledofs.fdim; k++)
   {
      actpart->pdis[kk].coupledofs.a.ia[k][imyrank+1]=2;
   }
}
else /*----------------------------------------------- parallel version */
{
/*
   actpart->node[i] really loops only nodes with dofs updated on this proc
*/
   for (i=0; i<actpart->pdis[kk].numnp; i++) /* now loop only my nodes */
   {
      for (l=0; l<actpart->pdis[kk].node[i]->numdf; l++)
      {
         dof = actpart->pdis[kk].node[i]->dof[l];
         for (k=0; k<actpart->pdis[kk].coupledofs.fdim; k++)
         {
            if (actpart->pdis[kk].coupledofs.a.ia[k][0]==dof)
            {
               actpart->pdis[kk].coupledofs.a.ia[k][imyrank+1]=1;
               break;
            }
         }
      }
   }
}
/* ----- Allreduce the whole array, so every proc knows about where all 
                                                         coupledofs are */
#ifdef PARALLEL
sendsize = (actpart->pdis[kk].coupledofs.fdim)*(inprocs);
sendbuff = (int*)CCACALLOC(sendsize,sizeof(int));
recvbuff = (int*)CCACALLOC(sendsize,sizeof(int));
if (sendbuff==NULL || recvbuff==NULL) dserror("Allocation of temporary memory failed");
counter=0;
for (i=0; i<actpart->pdis[kk].coupledofs.fdim; i++)
{
   for (j=0; j<inprocs; j++)
   {
      sendbuff[counter] = actpart->pdis[kk].coupledofs.a.ia[i][j+1];
      counter++;
   }
}
MPI_Allreduce(sendbuff,
              recvbuff,
              sendsize,
              MPI_INT,
              MPI_SUM,
              actintra->MPI_INTRA_COMM);
counter=0;
for (i=0; i<actpart->pdis[kk].coupledofs.fdim; i++)
{
   for (j=0; j<inprocs; j++)
   {
      actpart->pdis[kk].coupledofs.a.ia[i][j+1] = recvbuff[counter];
      counter++;
   }
}
CCAFREE(sendbuff);CCAFREE(recvbuff);
#endif
/*------- count number of equations on partition including coupled dofs */
/*---------------------------------------- count the coupled ones first */
counter=0;
for (i=0; i<actpart->pdis[kk].coupledofs.fdim; i++)
{
   if (actpart->pdis[kk].coupledofs.a.ia[i][imyrank+1]!=0) counter++;
}
/*-------------------------------- count all dofs which are not coupled */
for (i=0; i<actpart->pdis[kk].numnp; i++)
{
   actnode = actpart->pdis[kk].node[i];
   for (l=0; l<actnode->numdf; l++)
   {
      dof = actnode->dof[l];
      iscoupled=0;
      for (k=0; k<actpart->pdis[kk].coupledofs.fdim; k++)
      {
         if (dof == actpart->pdis[kk].coupledofs.a.ia[k][0]) 
         {
            iscoupled=1;
            break;
         }
      }
      if (iscoupled==0) 
      {
         if (dof < actfield->dis[kk].numeq)
         counter++;
      }
   }
}
/*--- number of equations on this partition including the coupled ones */
*numeq = counter;
/* 
   An inter-proc coupled equation produces communications calculating the 
   sparsity mask of the matrix
   An inter-proc coupled equation produces communications adding element
   matrices to the system matrix
   An inter-proc coupled equation ruins the bandwith locally
   ->
   Now one processor has to be owner of the coupled equation. 
   Try to distribute the coupled equations equally over the processors

   The matrix has the following style (after allreduce on all procs the same):
   
               ----------------------
               | 12 | 2 | 0 | 1 | 0 |
               | 40 | 2 | 0 | 0 | 0 |
               | 41 | 1 | 2 | 1 | 1 |
               | 76 | 0 | 1 | 2 | 0 |
               ----------------------
               
               column 0                : number of the coupled equation
               column 1 - inprocs+1 : proc has coupled equation or not
                                         2 indicates owner of equation
*/
if (inprocs > 1)
{
   tmp = (int*)CCACALLOC(inprocs,sizeof(int));
   if (!tmp) dserror("Allocation of temporary memory failed");
   for (i=0; i<actpart->pdis[kk].coupledofs.fdim; i++)/*  loop coupled eqns */
   {
   /*--------------------------------- check whether its inter-proc eqn */
      inter_proc=0;
      for (j=0; j<inprocs; j++) inter_proc += actpart->pdis[kk].coupledofs.a.ia[i][j+1];
      if (inter_proc==1)/*----------------- no inter-processor coupling */
      {
         for (j=0; j<inprocs; j++)
         {
            if (actpart->pdis[kk].coupledofs.a.ia[i][j+1]==1) 
            {
               actpart->pdis[kk].coupledofs.a.ia[i][j+1]=2;
               break;
            }
         }
      }
      else/*----------------------------- eqn is an inter-proc equation */
      {
         /* there won't be more than a million procs in the near future....*/
         min=1000000;
         proc=-1;
         for (j=0; j<inprocs; j++)
         {
            if (actpart->pdis[kk].coupledofs.a.ia[i][j+1]==1) 
            {
               if (tmp[j]<=min)
               {
                  min = tmp[j];
                  proc = j;
               }
            }
         }
         actpart->pdis[kk].coupledofs.a.ia[i][proc+1]=2;
         tmp[proc] += 1;
      }
   }/* end loop over coupling eqns */
   CCAFREE(tmp);
}
/* procs who have not become owner of a coupling equation have to reduce there
   number of equations */
if (inprocs > 1)
{
   for (i=0; i<actpart->pdis[kk].coupledofs.fdim; i++)/* loop coupled eqns */
   {
      /* ------Yes, I am slave owner of an inter_proc coupling equation */
      if (actpart->pdis[kk].coupledofs.a.ia[i][imyrank+1]==1)
      {
         (*numeq) = (*numeq)-1;
      }
      /* master owner of equation do nothing, 'cause the equation has been
        counted already */
   }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of msr_numeq */






/*----------------------------------------------------------------------*
 |  allocate update put dofs in update in ascending order   m.gee 5/01  |
 *----------------------------------------------------------------------*/
void msr_update(FIELD         *actfield, 
                PARTITION     *actpart, 
                SOLVAR        *actsolv,
                INTRA         *actintra,
                AZ_ARRAY_MSR  *msr)
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
dstrc_enter("msr_update");
#endif
/*----------------------------------------------------------------------*/
imyrank = actintra->intra_rank;
inprocs = actintra->intra_nprocs;
/*------------------ make a local copy of the array actpart->coupledofs */
am_alloc_copy(&(actpart->pdis[kk].coupledofs),&coupledofs);
/*----------------------------------------------------------------------*/
update = msr->update.a.iv;
counter=0;
/*------------------------------------- loop the nodes on the partition */
for (i=0; i<actpart->pdis[kk].numnp; i++)
{
   actnode = actpart->pdis[kk].node[i];
   for (l=0; l<actnode->numdf; l++)
   {
      dof = actnode->dof[l];
      /* dirichlet condition on dof */
      if (dof >= actfield->dis[kk].numeq) continue;
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
if (counter != msr->numeq) dserror("Number of dofs in MSR-vector update wrong");
/*---------------------------- sort the vector update just to make sure */
qsort((int*) update, counter, sizeof(int), cmp_int);
/*----------------------------------------------------------------------*/
amdel(&coupledofs);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of msr_update */




/*----------------------------------------------------------------------*
 |  calculate number of nonzero entries and dof topology    m.gee 6/01  |
 *----------------------------------------------------------------------*/
void msr_nnz_topology(FIELD         *actfield, 
                      PARTITION     *actpart, 
                      SOLVAR        *actsolv,
                      INTRA         *actintra,
                      AZ_ARRAY_MSR  *msr,
                      int          **dof_connect)
{
int        i,j,k,l,m,n;
int        counter,counter2;
int        dof;
int        nnz;
int        iscoupled;
int       *update;
int        numeq;
int        actdof;
int        dofflag;
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

#ifdef PARALLEL 
MPI_Status status;
#endif

#ifdef DEBUG 
dstrc_enter("msr_nnz_topology");
#endif
/*----------------------------------------------------------------------*/
imyrank = actintra->intra_rank;
inprocs = actintra->intra_nprocs;
/*----------------------------------------------------------- shortcuts */
msr->nnz=0;
numeq  = msr->numeq;
update = msr->update.a.iv;
for (i=0; i<msr->numeq_total; i++) dof_connect[i]=NULL;
amdef("tmp",&dofpatch,MAX_NNZPERROW,1,"IV");
amzero(&dofpatch);
/*----------------------------------------------------------------------*/
for (i=0; i<numeq; i++)
{
   dof = update[i];
   /*------------------------------ check whether this is a coupled dof */
   iscoupled=0;
   dof_in_coupledofs(dof,actpart,&iscoupled);
   if (iscoupled==1) continue;
   /*--------------------------------- find the centernode for this dof */
   centernode=NULL;
   dof_find_centernode(dof,actpart,&centernode);
   dsassert(centernode!=NULL,"Cannot make sparsity pattern for Aztec");
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
            if (actnode->dof[l] < actfield->dis[kk].numeq)
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
coupledofs = &(actpart->pdis[kk].coupledofs);
for (i=0; i<coupledofs->fdim; i++)
{
   dof = coupledofs->a.ia[i][0];
   /*--------------------------- check for my own ownership of this dof */
   dofflag = coupledofs->a.ia[i][imyrank+1];
   /*----------- if dofflag is zero this dof has nothing to do with me */
   if (dofflag==0) continue;
   /*------------------------------------- find all patches to this dof */
   counter=0;
   for (j=0; j<actpart->pdis[kk].numnp; j++)
   {
      centernode=NULL;
      for (l=0; l<actpart->pdis[kk].node[j]->numdf; l++)
      {
         if (dof == actpart->pdis[kk].node[j]->dof[l])
         {
            centernode = actpart->pdis[kk].node[j];
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
                  if (actnode->dof[l] < actfield->dis[kk].numeq)
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
   dof_connect[dof][1] = dofflag;
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
/* make the who-has-to-send-whom-how-much-and-what-arrays and communicate */
#ifdef PARALLEL 
counter=0;
for (i=0; i<coupledofs->fdim; i++)
{
   dof = coupledofs->a.ia[i][0];
   /*-------------------------------------- find the master of this dof */
   for (j=1; j<coupledofs->sdim; j++)
   {
      if (coupledofs->a.ia[i][j]==2) 
      {
         dofmaster = j-1;
         break;
      }
   }
   /*-------------------------------------- find the slaves of this dof */
   for (j=1; j<coupledofs->sdim; j++)
   {
      if (coupledofs->a.ia[i][j]==1)
      {
         dofslave = j-1;
         /*----------------------------------- if I am master I receive */
         if (imyrank==dofmaster)
         {
            /* note:
               This is a nice example to do individual communication
               between two procs without communicating the size
               of the message in advance
            */
            /*--------------------------------- get envelope of message */
            MPI_Probe(dofslave,counter,actintra->MPI_INTRA_COMM,&status);
            /*----------------------------------- get lenght of message */
            MPI_Get_count(&status,MPI_INT,&recvlenght);
            /*--------------------------------------- realloc the array */
            dof_connect[dof] = (int*)CCAREALLOC(dof_connect[dof],
                                             (dof_connect[dof][0]+recvlenght)*
                                             sizeof(int));
            if (!dof_connect[dof]) dserror("Reallocation of dof_connect failed");
            /*----------------------------------------- receive message */
            MPI_Recv(&(dof_connect[dof][ dof_connect[dof][0] ]),recvlenght,MPI_INT,
                     dofslave,counter,actintra->MPI_INTRA_COMM,&status);
            /*--------------------------------- put new lenght to array */
            dof_connect[dof][0] += recvlenght;
            /*-------------------------------- delete the doubles again */
            for (m=2; m<dof_connect[dof][0]; m++)
            {
               actdof = dof_connect[dof][m];
               if (actdof==-1) continue;
               for (k=m+1; k<dof_connect[dof][0]; k++)
               {
                  if (dof_connect[dof][k] == actdof) 
                  dof_connect[dof][k] = -1;
               }
            }
            /*-------------------- move all remaining dofs to the front */
            counter2=2;
            for (m=2; m<dof_connect[dof][0]; m++)
            {
               if (dof_connect[dof][m]!=-1)
               {
                  dof_connect[dof][counter2] = dof_connect[dof][m];
                  counter2++;
               }
            }
            /*--------------------------------------- realloc the array */
            dof_connect[dof] = (int*)CCAREALLOC(dof_connect[dof],
                                             counter2*sizeof(int));
            if (!dof_connect[dof]) dserror("Reallocation of dof_connect failed");
            dof_connect[dof][0] = counter2;
         }
         if (imyrank==dofslave)
         {
            MPI_Send(
                     &(dof_connect[dof][3]),
                     (dof_connect[dof][0]-3),
                     MPI_INT,
                     dofmaster,
                     counter,
                     actintra->MPI_INTRA_COMM
                    );
         }
         counter++;
      }
   }
}
#endif
/*--------------------------------- now go through update and count nnz */
nnz=0;
for (i=0; i<msr->update.fdim; i++)
{
   dof = msr->update.a.iv[i];
   nnz += (dof_connect[dof][0]-2);
}
msr->nnz=nnz;
/*--------- last thing to do is to order dof_connect in ascending order */
for (i=0; i<numeq; i++)
{
   dof = update[i];
   qsort((int*)(&(dof_connect[dof][3])), dof_connect[dof][0]-3, sizeof(int), cmp_int);
}
/*----------------------------------------------------------------------*/
amdel(&dofpatch);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of msr_nnz_topology */







/*----------------------------------------------------------------------*
 |  check whether this dof is in coupledofs              m.gee 6/01     |
 *----------------------------------------------------------------------*/
void dof_in_coupledofs(int dof, PARTITION *actpart, int *iscoupled)
{
int       i;
#ifdef DEBUG 
dstrc_enter("dof_in_coupledofs");
#endif
/*----------------------------------------------------------------------*/
   for (i=0; i<actpart->pdis[kk].coupledofs.fdim; i++)
   {
      if (dof==actpart->pdis[kk].coupledofs.a.ia[i][0])
      {
         *iscoupled = 1;
         break;
      }
   }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of dof_in_coupledofs */






/*----------------------------------------------------------------------*
 |  find the node to this dof in partition               m.gee 6/01     |
 *----------------------------------------------------------------------*/
void dof_find_centernode(int dof, PARTITION *actpart, NODE **centernode)
{
int       j,k;
#ifdef DEBUG 
dstrc_enter("dof_find_centernode");
#endif
/*----------------------------------------------------------------------*/
for (j=0; j<actpart->pdis[kk].numnp; j++)
{
   for (k=0; k<actpart->pdis[kk].node[j]->numdf; k++)
   {
      if (actpart->pdis[kk].node[j]->dof[k] == dof)
      {
         *centernode = actpart->pdis[kk].node[j];
         goto nodefound1;
      }
   }
}
nodefound1:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of dof_find_centernode */





/*----------------------------------------------------------------------*
 |  make the DMSR vector bindx                              m.gee 6/01  |
 | for format see Aztec manual                                          |
 *----------------------------------------------------------------------*/
void msr_make_bindx(FIELD         *actfield, 
                       PARTITION     *actpart, 
                       SOLVAR        *actsolv,
                       AZ_ARRAY_MSR  *msr,
                       int          **dof_connect)
{
int        i,j,k,l;
int        count1,count2;
int        dof;

#ifdef DEBUG 
dstrc_enter("msr_make_bindx");
#endif
/*----------------------------------------------------------------------*/
/*-------------------------------------------------------------do bindx */
count1=0;
count2=msr->numeq+1;
for (i=0; i<msr->update.fdim; i++)
{
   dof = msr->update.a.iv[i];
   msr->bindx.a.iv[count1] = count2;
   count1++;
   for (j=3; j<dof_connect[dof][0]; j++)
   {
      msr->bindx.a.iv[count2] = dof_connect[dof][j];
      count2++;
   }   
}
msr->bindx.a.iv[msr->numeq] = msr->nnz+1;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of msr_make_bindx */







