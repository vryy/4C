/*!---------------------------------------------------------------------
\file
\brief contains memory managing functions

---------------------------------------------------------------------*/

/*!----------------------------------------------------------------------
\brief the header of everything
*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"


int cmp_int(const void *a, const void *b );
double cmp_double(const void *a, const void *b );

/*! 
\addtogroup MLPCG 
*//*! @{ (documentation module open)*/

 
/*!---------------------------------------------------------------------
\brief calculate the mask of a dbcsr matrix                                              

<pre>                                                        m.gee 6/02 
calculates the mask of a block distributed csr matrix
in structure DBCSR.
DBCSR contains a rectangular part of a distributed square
compressed sparse row matrix with additional information about
the dofs belonging to one node
It is supposed to serve as the finest grid matrix in an mulitlevel environment
This routine also renumbers the dofs in the nodes in a sense, that partition's
internal dofs are numbers first followed by the boundary dofs coupled to other 
processor's equations.
In the one processor or sequentiell case this is not done
</pre>
\param actfield   FIELD*      (i)   catual physical field (structure)                                
\param actpart    PARTITION*  (i)   this processors partition
\param actsolv    SOLVAR      (i)   general structure of solver informations                   
\param actintra   INTRA       (i)   the intra-communicator of this field                  
\param bdcsr      DBCSR       (o)   the empty dbcsr matrix                 
\warning this routine renumbers the dofs!
\return void                                               

------------------------------------------------------------------------*/
void mask_bdcsr(FIELD         *actfield, 
                PARTITION     *actpart, 
                SOLVAR        *actsolv,
                INTRA         *actintra, 
                DBCSR         *bdcsr)
{
int            i,j,counter;
int            imyrank;
int            inproc;
int            numeq;
int            numdf;
int          **dof_connect;
int            actdof; 
int          **blocks;
PARTDISCRET   *actpdiscret;
NODE          *actnode;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mask_bdcsr");
#endif
/*----------------------------------------------------------------------*/
bdcsr->numeq_total = actfield->dis[0].numeq;
bdcsr_numeq(actfield,actpart,actsolv,actintra,&numeq);
bdcsr->numeq       = numeq;
numdf              = actfield->dis[0].node[0].numdf;
/*----------------------------------------------------------------------*/
#ifdef PARALLEL
actpdiscret = actpart->pdis;
imyrank     = actintra->intra_rank;
inproc      = actintra->intra_nprocs;
if (inproc>1) mlpcg_renumberdofs(imyrank,inproc,actfield,actpdiscret,actintra,bdcsr);
#endif
/*---------------------------------------------- allocate vector update */
amdef("update",&(bdcsr->update),numeq,1,"IV");
amzero(&(bdcsr->update));
/*--------------------------------put dofs in update in ascending order */
bdcsr_update(actfield,actpart,actsolv,actintra,bdcsr);
/*------------------------ count number of nonzero entries on partition 
                                    and calculate dof connectivity list */
   /*
      dof_connect[i][0] = lenght of dof_connect[i]
      dof_connect[i][1] = iscoupled ( 1 or 2 ) 
      dof_connect[i][2] = dof
      dof_connect[i][ 3..dof_connect[i][0]-1 ] = connected dofs exluding itself 
   */
dof_connect = (int**)CCACALLOC(bdcsr->numeq_total,sizeof(int*));
if (!dof_connect) dserror("Allocation of dof_connect failed");
bdcsr_nnz_topology(actfield,actpart,actsolv,actintra,bdcsr,dof_connect);
/*--------------------------------------------------- allocate a and ja */
amdef("a"  ,&(bdcsr->a)  ,(bdcsr->nnz+1),1,"DV");
amdef("ja" ,&(bdcsr->ja) ,(bdcsr->nnz+1),1,"IV");
amdef("ia" ,&(bdcsr->ia) ,(bdcsr->numeq+1),1,"IV");
bdcsr_make_csr(actfield,actpart,actsolv,bdcsr,dof_connect);
/*--------------------------------------------------- make nodal blocks */
blocks=amdef("blocks",&(bdcsr->blocks),actpart->pdis[0].numnp,numdf+1,"IA");
for (i=0; i<actpart->pdis[0].numnp; i++)
{
   actnode = actpart->pdis[0].node[i];
   counter = 0;
   for (j=0; j<actnode->numdf; j++)
   {
      actdof = actnode->dof[j];
      if (actdof>=bdcsr->numeq_total) continue;
/*
      actdofindex = find_index(actdof,bdcsr->update.a.iv,bdcsr->update.fdim);
      dsassert(actdofindex!=-1,"error in local dof numbering");
      blocks[i][1+counter++] = actdofindex;
*/
      blocks[i][1+counter++] = actdof;
   }
   blocks[i][0] = counter;
}
/*---------------------------------------- delete the array dof_connect */
for (i=0; i<bdcsr->numeq_total; i++)
{
   if (dof_connect[i]) CCAFREE(dof_connect[i]);
}
CCAFREE(dof_connect);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mask_bdcsr */





/*!---------------------------------------------------------------------
\brief make csr matrix                                              

<pre>                                                        m.gee 6/02 
make csr matrix (see book of Y.Saad)
</pre>
\param actfield    FIELD*      (i)   catual physical field (structure)                                
\param actpart     PARTITION*  (i)   this processors partition
\param actsolv     SOLVAR      (i)   general structure of solver informations                   
\param bdcsr       DBCSR_ROOT  (i/o) the empty dbcsr matrix                 
\param dof_connect int**       (i)   the connectivity
\return void                                               

------------------------------------------------------------------------*/
void bdcsr_make_csr(FIELD         *actfield, 
                    PARTITION     *actpart, 
                    SOLVAR        *actsolv,
                    DBCSR         *bdcsr,
                    int          **dof_connect)
{
int        i,j;
int        count;
int        dof;
int        *ja,*ia,*update;
#ifdef DEBUG 
dstrc_enter("bdcsr_make_csr");
#endif
/*----------------------------------------------------------------------*/
ja     = bdcsr->ja.a.iv;
ia     = bdcsr->ia.a.iv;
update = bdcsr->update.a.iv;
/*----------------------------------------------------------------------*/
count=0;
ia[0] =0;
for (i=0; i<bdcsr->update.fdim; i++)
{
   dof = update[i];
   qsort((int*)(&(dof_connect[dof][2])), dof_connect[dof][0]-2, sizeof(int), cmp_int);
   for (j=2; j<dof_connect[dof][0]; j++)
   {
      ja[count] = dof_connect[dof][j];
      count++;
   }
   ia[i+1] = count;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of bdcsr_make_csr */





/*!---------------------------------------------------------------------
\brief count number of nonzeros on this proc                                              

<pre>                                                        m.gee 6/02 
count number of nonzeros on this proc and make a dof connecticvity list

dof_connect[i][0] = lenght of dof_connect[i]
dof_connect[i][1] = iscoupled ( 1 or 2 ) 
dof_connect[i][2] = dof
dof_connect[i][ 3..dof_connect[i][0]-1 ] = connected dofs exluding itself 

</pre>
\param actfield    FIELD*      (i)   catual physical field (structure)                                
\param actpart     PARTITION*  (i)   this processors partition
\param actsolv     SOLVAR      (i)   general structure of solver informations                   
\param actintra    INTRA       (i)   the intra-communicator of this field                  
\param bdcsr       DBCSR_ROOT  (i/o) the empty dbcsr matrix                 
\param dof_connect int**       (o)   the connectivity
\warning this routine does not support coupling!
\return void                                               

------------------------------------------------------------------------*/
void  bdcsr_nnz_topology(FIELD         *actfield, 
                         PARTITION     *actpart, 
                         SOLVAR        *actsolv,
                         INTRA         *actintra,
                         DBCSR         *bdcsr,
                         int          **dof_connect)
{
int        i,j,k,l,m;
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
int        recvlenght;
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
dstrc_enter("bdcsr_nnz_topology");
#endif
/*----------------------------------------------------------------------*/
imyrank = actintra->intra_rank;
inprocs = actintra->intra_nprocs;
/*----------------------------------------------------------- shortcuts */
bdcsr->nnz=0;
numeq  = bdcsr->numeq;
update = bdcsr->update.a.iv;
for (i=0; i<bdcsr->numeq_total; i++) dof_connect[i]=NULL;
amdef("tmp",&dofpatch,MAX_NNZPERROW,1,"IV");
amzero(&dofpatch);
/*----------------------------------------------------------------------*/
for (i=0; i<numeq; i++)
{
   dof = update[i];
   /*------------------------------ check whether this is a coupled dof */
   iscoupled=0;
   dof_in_coupledofs(dof,actpart,&iscoupled);
   if (iscoupled==1) dserror("BDCSR does not support dofcoupling");
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
      dof_connect[i][ 3..dof_connect[i][0]-1 ] = connected dofs exluding itself 
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
   /*--------------------------- check for my own ownership of this dof */
   dofflag = coupledofs->a.ia[i][imyrank+1];
   /*----------- if dofflag is zero this dof has nothing to do with me */
   if (dofflag==0) continue;
   /*------------------------------------- find all patches to this dof */
   counter=0;
   for (j=0; j<actpart->pdis[0].numnp; j++)
   {
      centernode=NULL;
      for (l=0; l<actpart->pdis[0].node[j]->numdf; l++)
      {
         if (dof == actpart->pdis[0].node[j]->dof[l])
         {
            centernode = actpart->pdis[0].node[j];
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
for (i=0; i<bdcsr->update.fdim; i++)
{
   dof = bdcsr->update.a.iv[i];
   nnz += (dof_connect[dof][0]-2);
}
bdcsr->nnz=nnz;
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
} /* end of bdcsr_nnz_topology */






/*!---------------------------------------------------------------------
\brief make list of dofs on each proc (sorted)                                              

<pre>                                                        m.gee 6/02 
make list of dofs on each proc (sorted in ascending order)
</pre>
\param actfield   FIELD*      (i)   catual physical field (structure)                                
\param actpart    PARTITION*  (i)   this processors partition
\param actsolv    SOLVAR      (i)   general structure of solver informations                   
\param actintra   INTRA       (i)   the intra-communicator of this field                  
\param bdcsr      DBCSR_ROOT  (o)   the empty dbcsr matrix                 
\warning this routine does not support coupling!
\return void                                               

------------------------------------------------------------------------*/
void bdcsr_update(FIELD          *actfield, 
                  PARTITION      *actpart, 
                  SOLVAR         *actsolv,
                  INTRA          *actintra,
                  DBCSR          *bdcsr)
{
int       i,l;
int       counter;
int      *update;
int       dof;
int       imyrank;
int       inprocs;
NODE     *actnode;
#ifdef DEBUG 
dstrc_enter("bdcsr_update");
#endif
/*----------------------------------------------------------------------*/
imyrank = actintra->intra_rank;
inprocs = actintra->intra_nprocs;
/*----------------------------------------------------------------------*/
update = bdcsr->update.a.iv;
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
          dserror("BDCSR does not support dofcoupling");
      }
      
   }
}
/*----------- check whether the correct number of dofs has been counted */
if (counter != bdcsr->numeq) dserror("Number of dofs in update-vector update wrong");
/*---------------------------- sort the vector update just to make sure */
qsort((int*) update, counter, sizeof(int), cmp_int);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of bdcsr_update */





/*!---------------------------------------------------------------------
\brief renumber dofs                                              

<pre>                                                        m.gee 6/02 
This routine renumbers the dofs in the nodes in a sense, that partition's
internal dofs are numbers first followed by the boundary dofs coupled to other 
processor's equations.
In the one processor or sequentiell case this is not done
</pre>
\param myrank     int         (i)   this processors intra rank
\param nproc      int         (i)   number of processors in this intra-communicator
\param actfield   FIELD*      (i)   catual physical field (structure)                                
\param actpart    PARTITION*  (i)   this processors partition
\param actintra   INTRA       (i)   the intra-communicator of this field                  
\warning <b>this routine renumbers the dofs and does NOT support coupling!!!!</b>
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_renumberdofs(int            myrank,
                        int            nproc,
                        FIELD         *actfield, 
                        PARTDISCRET   *actpdiscret, 
                        INTRA         *actintra,
                        DBCSR         *bdcsr)
{
#ifdef PARALLEL
int            i,j,k;
ELEMENT       *actele;
NODE          *actnode;
ARRAY          dofflag_a;
int           *dofflag;
ARRAY          newdof_a;
int           *newdof;
ARRAY          newdofrecv_a;
int           *newdofrecv;
int            locsize_send[MAXPROC];
int            locsize[MAXPROC];
int            startdof;
int            savedof;
int            actdof;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_renumberdofs");
#endif
/*----------------------------------------------------------------------*/
dofflag = amdef("tmp",&dofflag_a,bdcsr->numeq_total,1,"IV");
amzero(&dofflag_a);
newdof = amdef("tmp",&newdof_a,bdcsr->numeq_total,1,"IV");
amzero(&newdof_a);
newdofrecv = amdef("tmp",&newdofrecv_a,bdcsr->numeq_total,1,"IV");
amzero(&newdofrecv_a);
/*----------------------------------------------------------------------*/
for (i=0; i<nproc; i++) locsize_send[i]=0;
locsize_send[myrank]=bdcsr->numeq;
MPI_Allreduce(locsize_send,locsize,nproc,MPI_INT,MPI_SUM,actintra->MPI_INTRA_COMM);
startdof=0;
for (i=0; i<myrank; i++) startdof += locsize[i];
savedof = startdof;
/*---------------------------------------- make flags for boundary dofs */
for (i=0; i<actpdiscret->bou_numele; i++)
{
   actele = actpdiscret->bou_element[i];
   for (j=0; j<actele->numnp; j++)
   {
      actnode = actele->node[j];
      if (actnode->proc != myrank) continue;
      for (k=0; k<actnode->numdf; k++)
      {
         actdof = actnode->dof[k];
         if (actdof>=bdcsr->numeq_total) continue;
         dofflag[actdof]=1;
      }
   }
}
/*-------- give ascending numbers to the dofs, start with internal dofs */
for (i=0; i<actpdiscret->numnp; i++)
{
   actnode = actpdiscret->node[i];
   dsassert(actnode->proc==myrank,"Partitioning got mixed up");
   for (j=0; j<actnode->numdf; j++)
   {
      actdof = actnode->dof[j];
      if (actdof>=bdcsr->numeq_total) continue;
      if (dofflag[actdof]==1) continue;
      newdof[actdof] = startdof;
      startdof++;
   }
}
/* save the first dof which is interproc-coupled in the DBCSR_ROOT matrix */
bdcsr->firstcoupledof = startdof;
/*----- now number the dofs, which have interproc off-diagonal entries */
for (i=0; i<actpdiscret->numnp; i++)
{
   actnode = actpdiscret->node[i];
   for (j=0; j<actnode->numdf; j++)
   {
      actdof = actnode->dof[j];
      if (actdof>=bdcsr->numeq_total) continue;
      if (dofflag[actdof]!=1) continue;
      newdof[actdof] = startdof;
      startdof++;
   }
}
/*------------------------------------------------- check local dof sum */
if (startdof-savedof != bdcsr->numeq)
dserror("Local nuber of equations wrong");
/*--------------------------------------- allreduce the new dof numbers */
MPI_Allreduce(newdof,newdofrecv,bdcsr->numeq_total,MPI_INT,MPI_SUM,actintra->MPI_INTRA_COMM);
/*----------------------- now put new dofnumbers to the node structures */
for (i=0; i<actfield->dis[0].numnp; i++)
{
   actnode = &(actfield->dis[0].node[i]);
   for (j=0; j<actnode->numdf; j++)
   {
      actdof = actnode->dof[j];
      if (actdof>=bdcsr->numeq_total) continue;
      actnode->dof[j] = newdofrecv[actdof];
   }
}
/*----------------------------------------------------------------------*/
amdel(&dofflag_a);
amdel(&newdof_a);
amdel(&newdofrecv_a);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
#endif
return;
} /* end of mlpcg_renumberdofs */






/*!---------------------------------------------------------------------
\brief count processors number of equations                                              

<pre>                                                        m.gee 6/02 
counts processors number of equations, including coupling conditions
</pre>
\param actfield   FIELD*      (i)   catual physical field (structure)                                
\param actpart    PARTITION*  (i)   this processors partition
\param actsolv    SOLVAR      (i)   general structure of solver informations                   
\param actintra   INTRA       (i)   the intra-communicator of this field                  
\param numeq      int*        (o)   local number of dofs
\return void                                               

------------------------------------------------------------------------*/
void bdcsr_numeq(FIELD         *actfield, 
                 PARTITION     *actpart, 
                 SOLVAR        *actsolv,
                 INTRA         *actintra,
                 int           *numeq)
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
dstrc_enter("bdcsr_numeq");
#endif
/*----------------------------------------------------------------------*/
imyrank = actintra->intra_rank;
inprocs = actintra->intra_nprocs;
/*------------------------- first make a list of dofs which are coupled */
/*----------------------------------- estimate size of coupdofs to 5000 */
amdef("coupledofs",&(actpart->pdis[0].coupledofs),5000,1,"IV");
amzero(&(actpart->pdis[0].coupledofs));
counter=0;
/*-------------------------------- loop all nodes and find coupled dofs */
for (i=0; i<actfield->dis[0].numnp; i++)
{
   actnode = &(actfield->dis[0].node[i]);
   if (actnode->gnode->couple==NULL && actnode->gnode->dirich==NULL) continue;
   if (actnode->gnode->couple==NULL) continue;
   for (l=0; l<actnode->numdf; l++)
   {
      if (actnode->dof[l]>=actfield->dis[0].numeq) continue;
      /* there is coupling on this dof */
      if (actnode->gnode->couple->couple.a.ia[l][0] != 0 ||
          actnode->gnode->couple->couple.a.ia[l][1] != 0 )
      {
         if (counter>=actpart->pdis[0].coupledofs.fdim) 
         amredef(&(actpart->pdis[0].coupledofs),(actpart->pdis[0].coupledofs.fdim+5000),1,"IV");
         /* the coupled dof could be dirichlet conditioned */
         if (actnode->dof[l]<actfield->dis[0].numeq)
         {
            actpart->pdis[0].coupledofs.a.iv[counter] = actnode->dof[l];
            counter++;
         }
      }
   }
}
amredef(&(actpart->pdis[0].coupledofs),counter,1,"IV");
/*---------------------------------- delete the doubles in coupledofs */
for (i=0; i<actpart->pdis[0].coupledofs.fdim; i++)
{
   if (actpart->pdis[0].coupledofs.a.iv[i]==-1) continue;
   dof = actpart->pdis[0].coupledofs.a.iv[i];
   for (j=i+1; j<actpart->pdis[0].coupledofs.fdim; j++)
   {
      if (actpart->pdis[0].coupledofs.a.iv[j]==dof) actpart->pdis[0].coupledofs.a.iv[j]=-1;
   }
}
/*--------- move all remaining coupdofs to the front and redefine again */
counter=0;
for (i=0; i<actpart->pdis[0].coupledofs.fdim; i++)
{
   if (actpart->pdis[0].coupledofs.a.iv[i]!=-1)
   {
      actpart->pdis[0].coupledofs.a.iv[counter] = actpart->pdis[0].coupledofs.a.iv[i];
      counter++;
   }
}
amredef(&(actpart->pdis[0].coupledofs),counter,inprocs+1,"IA");
/*------------------- the newly allocated columns have to be initialized */
for (i=1; i<actpart->pdis[0].coupledofs.sdim; i++)
for (j=0; j<actpart->pdis[0].coupledofs.fdim; j++) 
actpart->pdis[0].coupledofs.a.ia[j][i]=0;

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
   for (k=0; k<actpart->pdis[0].coupledofs.fdim; k++)
   {
      actpart->pdis[0].coupledofs.a.ia[k][imyrank+1]=2;
   }
}
else /*----------------------------------------------- parallel version */
{
/*
   actpart->node[i] really loops only nodes with dofs updated on this proc
*/
   for (i=0; i<actpart->pdis[0].numnp; i++) /* now loop only my nodes */
   {
      for (l=0; l<actpart->pdis[0].node[i]->numdf; l++)
      {
         dof = actpart->pdis[0].node[i]->dof[l];
         for (k=0; k<actpart->pdis[0].coupledofs.fdim; k++)
         {
            if (actpart->pdis[0].coupledofs.a.ia[k][0]==dof)
            {
               actpart->pdis[0].coupledofs.a.ia[k][imyrank+1]=1;
               break;
            }
         }
      }
   }
}
/* ----- Allreduce the whole array, so every proc knows about where all 
                                                         coupledofs are */
#ifdef PARALLEL
sendsize = (actpart->pdis[0].coupledofs.fdim)*(inprocs);
sendbuff = (int*)CCACALLOC(sendsize,sizeof(int));
recvbuff = (int*)CCACALLOC(sendsize,sizeof(int));
if (sendbuff==NULL || recvbuff==NULL) dserror("Allocation of temporary memory failed");
counter=0;
for (i=0; i<actpart->pdis[0].coupledofs.fdim; i++)
{
   for (j=0; j<inprocs; j++)
   {
      sendbuff[counter] = actpart->pdis[0].coupledofs.a.ia[i][j+1];
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
for (i=0; i<actpart->pdis[0].coupledofs.fdim; i++)
{
   for (j=0; j<inprocs; j++)
   {
      actpart->pdis[0].coupledofs.a.ia[i][j+1] = recvbuff[counter];
      counter++;
   }
}
CCAFREE(sendbuff);CCAFREE(recvbuff);
#endif
/*------- count number of equations on partition including coupled dofs */
/*---------------------------------------- count the coupled ones first */
counter=0;
for (i=0; i<actpart->pdis[0].coupledofs.fdim; i++)
{
   if (actpart->pdis[0].coupledofs.a.ia[i][imyrank+1]!=0) counter++;
}
/*-------------------------------- count all dofs which are not coupled */
for (i=0; i<actpart->pdis[0].numnp; i++)
{
   actnode = actpart->pdis[0].node[i];
   for (l=0; l<actnode->numdf; l++)
   {
      dof = actnode->dof[l];
      iscoupled=0;
      for (k=0; k<actpart->pdis[0].coupledofs.fdim; k++)
      {
         if (dof == actpart->pdis[0].coupledofs.a.ia[k][0]) 
         {
            iscoupled=1;
            break;
         }
      }
      if (iscoupled==0) 
      {
         if (dof < actfield->dis[0].numeq)
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
   for (i=0; i<actpart->pdis[0].coupledofs.fdim; i++)/*  loop coupled eqns */
   {
   /*--------------------------------- check whether its inter-proc eqn */
      inter_proc=0;
      for (j=0; j<inprocs; j++) inter_proc += actpart->pdis[0].coupledofs.a.ia[i][j+1];
      if (inter_proc==1)/*----------------- no inter-processor coupling */
      {
         for (j=0; j<inprocs; j++)
         {
            if (actpart->pdis[0].coupledofs.a.ia[i][j+1]==1) 
            {
               actpart->pdis[0].coupledofs.a.ia[i][j+1]=2;
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
            if (actpart->pdis[0].coupledofs.a.ia[i][j+1]==1) 
            {
               if (tmp[j]<=min)
               {
                  min = tmp[j];
                  proc = j;
               }
            }
         }
         actpart->pdis[0].coupledofs.a.ia[i][proc+1]=2;
         tmp[proc] += 1;
      }
   }/* end loop over coupling eqns */
   CCAFREE(tmp);
}
/* procs who have not become owner of a coupling equation have to reduce there
   number of equations */
if (inprocs > 1)
{
   for (i=0; i<actpart->pdis[0].coupledofs.fdim; i++)/* loop coupled eqns */
   {
      /* ------Yes, I am slave owner of an inter_proc coupling equation */
      if (actpart->pdis[0].coupledofs.a.ia[i][imyrank+1]==1)
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
} /* end of bdcsr_numeq */


/*! @} (documentation module close)*/
