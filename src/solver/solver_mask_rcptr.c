#include "../headers/standardtypes.h"
#include "../headers/solution.h"
int cmp_int(const void *a, const void *b );
double cmp_double(const void *a, const void *b );


/*----------------------------------------------------------------------*
 |  calculate the mask of an rc_ptr matrix               m.gee 1/02     |
 *----------------------------------------------------------------------*/
void mask_rc_ptr(FIELD         *actfield, 
                 PARTITION     *actpart, 
                 SOLVAR        *actsolv,
                 INTRA         *actintra, 
                 RC_PTR        *rc_ptr)
{
int       i,j,k,l;
int       numeq;
int     **dof_connect;
#ifdef DEBUG 
dstrc_enter("mask_rc_ptr");
#endif
/*----------------------------------------------------------------------*/
/* remember some facts:
   PARTITION is different on every proc.
   AZ_ARRAY_MSR will be different on every proc
   FIELD is the same everywhere
   In this routine, the vectors update and bindx and val are determined
   in size and allocated, the contents of the vectors update and bindx 
   are calculated
/*------------------------------------------- put total size of problem */
rc_ptr->numeq_total = actfield->numeq;
/* count number of eqns on proc and build processor-global couplingdof 
                                                                 matrix */
msr_numeq(actfield,actpart,actsolv,actintra,&numeq);
rc_ptr->numeq = numeq;
/*---------------------------------------------- allocate vector update */
amdef("update",&(rc_ptr->update),numeq,1,"IV");
amzero(&(rc_ptr->update));
/*--------------------------------put dofs in update in ascending order */
rc_ptr_update(actfield,actpart,actsolv,actintra,rc_ptr);
/*------------------------ count number of nonzero entries on partition 
                                    and calculate dof connectivity list */
   /*
      dof_connect[i][0] = lenght of dof_connect[i]
      dof_connect[i][1] = iscoupled ( 1 or 2 ) 
      dof_connect[i][2] = dof
      dof_connect[i][ 2..dof_connect[i][0]-1 ] = connected dofs exluding itself 
   */
dof_connect = (int**)calloc(rc_ptr->numeq_total,sizeof(int*));
if (!dof_connect) dserror("Allocation of dof_connect failed");
rc_ptr_nnz_topology(actfield,actpart,actsolv,actintra,rc_ptr,dof_connect);
/*----------------------------------------------------- allocate arrays */
/*                                                     see MUMPS manual */
amdef("rowptr" ,&(rc_ptr->rowptr) ,rc_ptr->numeq    ,1,"IV");
amdef("irn_loc",&(rc_ptr->irn_loc),rc_ptr->nnz      ,1,"IV");
amdef("jcn_loc",&(rc_ptr->jcn_loc),rc_ptr->nnz      ,1,"IV");
amdef("A"      ,&(rc_ptr->A_loc)  ,rc_ptr->nnz      ,1,"DV");
/*------------------------------------------------------ allocate bindx */
amdef("bindx",&(rc_ptr->bindx),(rc_ptr->nnz+1),1,"IV");
/*---------------------------------------------------------- make bindx */
rc_ptr_make_bindx(actfield,actpart,actsolv,rc_ptr,dof_connect);
/*----------------- make rowptr, irn_loc, jcn_loc from bindx and update */
rc_ptr_make_sparsity(rc_ptr);
/*------------------------------------- make irn, jcn, irn_loc, jcn_loc */

/*---------------------------------------- delete the array dof_connect */
for (i=0; i<rc_ptr->numeq_total; i++)
{
   if (!dof_connect[i]) free(dof_connect[i]);
}
free(dof_connect);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mask_rc_ptr */



/*----------------------------------------------------------------------*
 |  allocate update put dofs in update in ascending order   m.gee 1/02  |
 *----------------------------------------------------------------------*/
int  rc_ptr_update(FIELD         *actfield, 
                   PARTITION     *actpart, 
                   SOLVAR        *actsolv,
                   INTRA         *actintra,
                   RC_PTR        *rc_ptr)
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
dstrc_enter("rc_ptr_update");
#endif
/*----------------------------------------------------------------------*/
imyrank = actintra->intra_rank;
inprocs = actintra->intra_nprocs;
/*------------------ make a local copy of the array actpart->coupledofs */
am_alloc_copy(&(actpart->coupledofs),&coupledofs);
/*------------------------------------- loop the nodes on the partition */
update = rc_ptr->update.a.iv;
counter=0;
for (i=0; i<actpart->numnp; i++)
{
   actnode = actpart->node[i];
   for (l=0; l<actnode->numdf; l++)
   {
      dof = actnode->dof[l];
      /* dirichlet condition on dof */
      if (dof >= actfield->numeq) continue;
      /* no condition on dof */
      if (actnode->c==NULL)
      {
         update[counter] = dof;
         counter++;
         continue;
      }
      /* no coupling on dof */
      if (actnode->c->iscoupled==0)
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
if (counter != rc_ptr->numeq) dserror("Number of dofs in update wrong");
/*---------------------------- sort the vector update just to make sure */
qsort((int*) update, counter, sizeof(int), cmp_int);
/*----------------------------------------------------------------------*/
amdel(&coupledofs);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of rc_ptr_update */



/*----------------------------------------------------------------------*
 |  calculate number of nonzero entries and dof topology    m.gee 1/02  |
 *----------------------------------------------------------------------*/
int  rc_ptr_nnz_topology(FIELD         *actfield, 
                         PARTITION    *actpart, 
                         SOLVAR       *actsolv,
                         INTRA        *actintra,
                         RC_PTR       *rc_ptr,
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
dstrc_enter("rc_ptr_nnz_topology");
#endif
/*----------------------------------------------------------------------*/
imyrank = actintra->intra_rank;
inprocs = actintra->intra_nprocs;
/*----------------------------------------------------------- shortcuts */
rc_ptr->nnz=0;
numeq  = rc_ptr->numeq;
update = rc_ptr->update.a.iv;
for (i=0; i<rc_ptr->numeq_total; i++) dof_connect[i]=NULL;
amdef("tmp",&dofpatch,1000,1,"IV");
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
   dof_find_centernode(dof,actpart,&centernode);
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
            if (actnode->dof[l] < actfield->numeq)
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
   dof_connect[dof] = (int*)calloc(counter2+3,sizeof(int));
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
coupledofs = &(actpart->coupledofs);
for (i=0; i<coupledofs->fdim; i++)
{
   dof = coupledofs->a.ia[i][0];
   /*--------------------------- check for my own ownership of this dof */
   dofflag = coupledofs->a.ia[i][imyrank+1];
   /*----------- if dofflag is zero this dof has nothing to do with me */
   if (dofflag==0) continue;
   /*------------------------------------- find all patches to this dof */
   counter=0;
   for (j=0; j<actpart->numnp; j++)
   {
      centernode=NULL;
      for (l=0; l<actpart->node[j]->numdf; l++)
      {
         if (dof == actpart->node[j]->dof[l])
         {
            centernode = actpart->node[j];
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
                  if (actnode->dof[l] < actfield->numeq)
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
   dof_connect[dof] = (int*)calloc(counter2+3,sizeof(int));
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
            dof_connect[dof] = (int*)realloc(dof_connect[dof],
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
            dof_connect[dof] = (int*)realloc(dof_connect[dof],
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
for (i=0; i<rc_ptr->update.fdim; i++)
{
   dof = rc_ptr->update.a.iv[i];
   nnz += (dof_connect[dof][0]-2);
}
rc_ptr->nnz=nnz;
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
} /* end of rc_ptr_nnz_topology */


/*----------------------------------------------------------------------*
 |  make the DMSR vector bindx                              m.gee 1/02  |
 | for format see Aztec manual                                          |
 *----------------------------------------------------------------------*/
int  rc_ptr_make_bindx(FIELD         *actfield, 
                       PARTITION     *actpart, 
                       SOLVAR        *actsolv,
                       RC_PTR        *rc_ptr,
                       int          **dof_connect)
{
int        i,j,k,l;
int        count1,count2;
int        dof;

#ifdef DEBUG 
dstrc_enter("rc_ptr_make_bindx");
#endif
/*----------------------------------------------------------------------*/
/*-------------------------------------------------------------do bindx */
count1=0;
count2=rc_ptr->numeq+1;
for (i=0; i<rc_ptr->update.fdim; i++)
{
   dof = rc_ptr->update.a.iv[i];
   rc_ptr->bindx.a.iv[count1] = count2;
   count1++;
   for (j=3; j<dof_connect[dof][0]; j++)
   {
      rc_ptr->bindx.a.iv[count2] = dof_connect[dof][j];
      count2++;
   }   
}
rc_ptr->bindx.a.iv[rc_ptr->numeq] = rc_ptr->nnz+1;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of rc_ptr_make_bindx */
/*----------------------------------------------------------------------*
 |  make the vectors                                        m.gee 1/02  |
 | irn_loc, jcn_loc, rowptr from update and bindx                       |
 *----------------------------------------------------------------------*/
int  rc_ptr_make_sparsity(RC_PTR        *rc_ptr)
{
int        i,j,k,l;
int        start,end;
int        actdof;
int        numeq;
int        numeq_total;
int        nnz;
int       *update;
int       *bindx;
int       *irn;
int       *jcn;
int       *rptr;

#ifdef DEBUG 
dstrc_enter("rc_ptr_make_sparsity");
#endif
/*----------------------------------------------------------------------*/
numeq       = rc_ptr->numeq;
numeq_total = rc_ptr->numeq_total;
nnz         = rc_ptr->nnz;
update      = rc_ptr->update.a.iv;
bindx       = rc_ptr->bindx.a.iv;
irn         = rc_ptr->irn_loc.a.iv;
jcn         = rc_ptr->jcn_loc.a.iv;
rptr        = rc_ptr->rowptr.a.iv;
/*------------------------------------------ loop all dofs on this proc */
for (i=0; i<numeq; i++)
{
   actdof = update[i];
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of rc_ptr_make_sparsity */
