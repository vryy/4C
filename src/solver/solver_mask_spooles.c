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


#ifdef SPOOLES_PACKAGE


#include "../headers/standardtypes.h"
#include "../solver/solver.h"
INT cmp_int(const void *a, const void *b );
DOUBLE cmp_double(const void *a, const void *b );


/*----------------------------------------------------------------------*
 |  calculate the mask of an spooles matrix              m.gee 5/01     |
 *----------------------------------------------------------------------*/
void mask_spooles(FIELD         *actfield, 
                  PARTITION     *actpart, 
                  SOLVAR        *actsolv,
                  INTRA         *actintra, 
                  SPOOLMAT      *spo)
{
INT       i;
INT       numeq;
INT     **dof_connect;
ARRAY     bindx_a;
INT      *bindx;

ELEMENT  *actele;

#ifdef DEBUG 
dstrc_enter("mask_spooles");
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
spo->numeq_total = actfield->dis[0].numeq;
/* count number of eqns on proc and build processor-global couplingdof 
                                                                 matrix */
mask_numeq(actfield,actpart,actsolv,actintra,&numeq,0);
spo->numeq = numeq;
/*---------------------------------------------- allocate vector update */
amdef("update",&(spo->update),numeq,1,"IV");
amzero(&(spo->update));
/*--------------------------------put dofs in update in ascending order */
spo_update(actfield,actpart,actsolv,actintra,spo);
/*------------------------ count number of nonzero entries on partition 
                                    and calculate dof connectivity list */
   /*
      dof_connect[i][0] = lenght of dof_connect[i]
      dof_connect[i][1] = iscoupled ( 1 or 2 ) 
      dof_connect[i][2] = dof
      dof_connect[i][ 2..dof_connect[i][0]-1 ] = connected dofs exluding itself 
   */
dof_connect = (INT**)CCACALLOC(spo->numeq_total,sizeof(INT*));
if (!dof_connect) dserror("Allocation of dof_connect failed");
/*---------------------- make the dof_connect list locally on each proc */
spo_nnz_topology(actfield,actpart,actsolv,actintra,spo,dof_connect);
/*----------------------------------------------------- allocate arrays */
amdef("rowptr" ,&(spo->rowptr) ,spo->numeq+1  ,1,"IV");
amdef("irn_loc",&(spo->irn_loc),spo->nnz      ,1,"IV");
amdef("jcn_loc",&(spo->jcn_loc),spo->nnz      ,1,"IV");
amdef("A"      ,&(spo->A_loc)  ,spo->nnz      ,1,"DV");
/*------------------------------------------------------ allocate bindx */
bindx = amdef("bindx",&(bindx_a),(spo->nnz+1),1,"IV");
/*---------------------------------------------------------- make bindx */
spo_make_bindx(actfield,actpart,actsolv,spo,dof_connect,bindx);
/*----------------- make rowptr, irn_loc, jcn_loc from bindx and update */
spo_make_sparsity(spo,bindx);
/*---------------------------------------- delete the array dof_connect */
for (i=0; i<spo->numeq_total; i++)
{
   if (dof_connect[i]) CCAFREE(dof_connect[i]);
}
CCAFREE(dof_connect);
amdel(&bindx_a);


/* make the index vector for faster assembling */
  for (i=0; i<actpart->pdis[0].numele; i++)
  {
    actele = actpart->pdis[0].element[i];
    spo_make_index(actfield,actpart,actintra,actele,spo);
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mask_spooles */


/*----------------------------------------------------------------------*
 |  make the vectors                                        m.gee 1/02  |
 | irn_loc, jcn_loc, rowptr from update and bindx                       |
 *----------------------------------------------------------------------*/
void     spo_make_sparsity(SPOOLMAT        *spo,
                           INT           *bindx)
{
INT        i,j;
INT        start,end;
INT        counter;
INT        actdof;
INT        numeq;
INT        numeq_total;
INT        nnz;
INT       *update;
INT       *irn;
INT       *jcn;
INT       *rptr;

#ifdef DEBUG 
dstrc_enter("spo_make_sparsity");
#endif
/*----------------------------------------------------------------------*/
numeq       = spo->numeq;
numeq_total = spo->numeq_total;
nnz         = spo->nnz;
update      = spo->update.a.iv;
irn         = spo->irn_loc.a.iv;
jcn         = spo->jcn_loc.a.iv;
rptr        = spo->rowptr.a.iv;
/*------------------------------------------ loop all dofs on this proc */
counter=0;
for (i=0; i<numeq; i++)
{
   actdof   = update[i];
   start    = bindx[i];
   end      = bindx[i+1];
   rptr[i]  = counter;
   j=start;
   while (j<end && bindx[j]<actdof)/* dofs lower then actdof */ 
   {
      irn[counter]=actdof;
      jcn[counter]=bindx[j];
      counter++;
      j++;
   }
   /*------------------------------- main diagonal of actdof */
   irn[counter]=actdof;
   jcn[counter]=actdof;
   counter++;
   while(j<end)/*------------------- dofs higher then actdof */
   {
      irn[counter]=actdof;
      jcn[counter]=bindx[j];
      counter++;
      j++;
   }
}
rptr[i]=counter;
if (counter != nnz) dserror("sparsity mask failure");
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of spo_make_sparsity */



/*----------------------------------------------------------------------*
 |  make the DMSR vector bindx                              m.gee 1/02  |
 | for format see Aztec manual                                          |
 *----------------------------------------------------------------------*/
void    spo_make_bindx(FIELD         *actfield, 
                       PARTITION     *actpart, 
                       SOLVAR        *actsolv,
                       SPOOLMAT      *spo,
                       INT          **dof_connect,
                       INT           *bindx)
{
INT        i,j;
INT        count1,count2;
INT        dof;

#ifdef DEBUG 
dstrc_enter("spo_make_bindx");
#endif
/*----------------------------------------------------------------------*/
/*-------------------------------------------------------------do bindx */
count1=0;
count2=spo->numeq+1;
for (i=0; i<spo->update.fdim; i++)
{
   dof = spo->update.a.iv[i];
   bindx[count1] = count2;
   count1++;
   for (j=3; j<dof_connect[dof][0]; j++)
   {
      bindx[count2] = dof_connect[dof][j];
      count2++;
   }   
}
bindx[spo->numeq] = spo->nnz+1;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of spo_make_bindx */


/*----------------------------------------------------------------------*
 |  calculate number of nonzero entries and dof topology    m.gee 1/02  |
 *----------------------------------------------------------------------*/
void  spo_nnz_topology(FIELD         *actfield, 
                       PARTITION    *actpart, 
                       SOLVAR       *actsolv,
                       INTRA        *actintra,
                       SPOOLMAT     *spo,
                       INT         **dof_connect)
{
INT        i,j,k,l,m;
INT        counter,counter2;
INT        dof;
INT        nnz;
INT        iscoupled;
INT       *update;
INT        numeq;
INT        actdof;
INT        dofflag;
#ifdef PARALLEL
INT        dofmaster;
INT        dofslave;
INT        recvlenght;
#endif
NODE      *centernode;
NODE      *actnode;
ELEMENT   *actele;
ARRAY      dofpatch;
ARRAY     *coupledofs;
INT        imyrank;
INT        inprocs;

NODE     **node_dof;
INT        max_dof,dof_id;

#ifdef PARALLEL 
MPI_Status status;
#endif

#ifdef DEBUG 
dstrc_enter("spo_nnz_topology");
#endif
/*----------------------------------------------------------------------*/
imyrank = actintra->intra_rank;
inprocs = actintra->intra_nprocs;
/*----------------------------------------------------------- shortcuts */
spo->nnz=0;
numeq  = spo->numeq;
update = spo->update.a.iv;
for (i=0; i<spo->numeq_total; i++) dof_connect[i]=NULL;
amdef("tmp",&dofpatch,1000,1,"IV");
amzero(&dofpatch);

/* create node pointers for all dofs in this partition */
  max_dof = 0;
  for (j=0; j<actpart->pdis[0].numnp; j++)
  {
    for (k=0; k<actpart->pdis[0].node[j]->numdf; k++)
    {
      if (actpart->pdis[0].node[j]->dof[k] >= max_dof)
      {
        max_dof = actpart->pdis[0].node[j]->dof[k];
      }
    }
  }
  /* allocate pointer vector to the nodes */
  node_dof = (NODE**)CCACALLOC(max_dof+1,sizeof(NODE*));
  /* store pointers to nodes in node_dof at position accord. to dof */
  for (j=0; j<actpart->pdis[0].numnp; j++)
  {
    for (k=0; k<actpart->pdis[0].node[j]->numdf; k++)
    {
      dof_id = actpart->pdis[0].node[j]->dof[k];
      dsassert(dof_id <= max_dof,"zu kleiner node_dof Vector");
      node_dof[dof_id] = actpart->pdis[0].node[j];
    }
  }

/*----------------------------------------------------------------------*/
for (i=0; i<numeq; i++)
{
   dof = update[i];
   /*------------------------------ check whether this is a coupled dof */
   iscoupled=0;
   dof_in_coupledofs(dof,actpart,&iscoupled);
   if (iscoupled==1) continue;
   /*--------------------------------- find the centernode for this dof */
   centernode = node_dof[dof];
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
   dof_connect[dof] = (INT*)CCACALLOC(counter2+3,sizeof(INT));
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
   dof_connect[dof] = (INT*)CCACALLOC(counter2+3,sizeof(INT));
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
            dof_connect[dof] = (INT*)CCAREALLOC(dof_connect[dof],
                                             (dof_connect[dof][0]+recvlenght)*
                                             sizeof(INT));
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
            dof_connect[dof] = (INT*)CCAREALLOC(dof_connect[dof],
                                             counter2*sizeof(INT));
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
for (i=0; i<spo->update.fdim; i++)
{
   dof = spo->update.a.iv[i];
   nnz += (dof_connect[dof][0]-2);
}
spo->nnz=nnz;
/*--------- last thing to do is to order dof_connect in ascending order */
for (i=0; i<numeq; i++)
{
   dof = update[i];
   qsort((INT*)(&(dof_connect[dof][3])), dof_connect[dof][0]-3, sizeof(INT), cmp_int);
}
/*----------------------------------------------------------------------*/
CCAFREE(node_dof);
amdel(&dofpatch);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of spo_nnz_topology */







/*----------------------------------------------------------------------*
 |  allocate update put dofs in update in ascending order   m.gee 5/01  |
 *----------------------------------------------------------------------*/
void spo_update(FIELD         *actfield, 
                PARTITION     *actpart, 
                SOLVAR        *actsolv,
                INTRA         *actintra,
                SPOOLMAT      *spo)
{
INT       i,k,l;
INT       counter;
INT      *update;
INT       dof;
INT       foundit;
INT       imyrank;
INT       inprocs;
NODE     *actnode;
ARRAY     coupledofs;
#ifdef DEBUG 
dstrc_enter("spo_update");
#endif
/*----------------------------------------------------------------------*/
imyrank = actintra->intra_rank;
inprocs = actintra->intra_nprocs;
/*------------------ make a local copy of the array actpart->coupledofs */
am_alloc_copy(&(actpart->pdis[0].coupledofs),&coupledofs);
/*----------------------------------------------------------------------*/
update = spo->update.a.iv;
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
if (counter != spo->numeq) dserror("Number of dofs in spooles-vector update wrong");
/*---------------------------- sort the vector update just to make sure */
qsort((INT*) update, counter, sizeof(INT), cmp_int);
/*----------------------------------------------------------------------*/
amdel(&coupledofs);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of spo_update */



/*----------------------------------------------------------------------*/
/*!
 \brief 

 This routine determines the location vactor for the actele and stores it
 in the element structure.  Furthermore for each component [i][j] in the 
 element stiffness matrix the position in the 1d sparse matrix is 
 calculated and stored in actele->index[i][j]. These can be used later on 
 for the assembling procedure.
  
 \param actfield  *FIELD        (i)  the field we are working on
 \param actpart   *PARTITION    (i)  the partition we are working on
 \param actintra  *INTRA        (i)  the intra-communicator we do not need
 \param actele    *ELEMENT      (i)  the element we would like to work with
 \param spo1      *SPOOLMAT     (i)  the sparse matrix we will assemble into

 \author mn
 \date 07/04

 */
/*----------------------------------------------------------------------*/
void spo_make_index(
    FIELD                 *actfield, 
    PARTITION             *actpart,
    INTRA                 *actintra,
    ELEMENT               *actele,
    struct _SPOOLMAT      *spo1
    )
{

  INT         i,j,k,l,counter;          /* some counter variables */
  INT         istwo=0;
  INT         start,index,lenght;       /* some more special-purpose counters */
  INT         ii,jj;                    /* counter variables for system matrix */
  INT         ii_iscouple;              /* flag whether ii is a coupled dof */
  INT         ii_owner;                 /* who is owner of dof ii -> procnumber */
  INT         ii_index;                 /* place of ii in dmsr format */
  INT         jj_index;                 /* place of jj in dmsr format */
  INT         nd,ndnd;                  /* size of estif */
  INT         nnz;                      /* number of nonzeros in sparse system matrix */
  INT         numeq_total;              /* total number of equations */
  INT         numeq;                    /* number of equations on this proc */
  INT         lm[MAXDOFPERELE];         /* location vector for this element */
  INT         owner[MAXDOFPERELE];      /* the owner of every dof */
  INT         myrank;                   /* my intra-proc number */
  INT         nprocs;                   /* my intra- number of processes */
  DOUBLE    **emass;                    /* element matrix to be added to system matrix */
  INT        *update;                   /* vector update see AZTEC manual */
  DOUBLE     *A_loc;                    /*    "       A_loc see MUMPS manual */
  DOUBLE     *B_loc;                    /*    "       A_loc see MUMPS manual */
  INT        *irn;                      /*    "       irn see MUMPS manual */
  INT        *jcn;                      /*    "       jcn see MUMPS manual */
  INT        *rowptr;                   /*    "       rowptr see rc_ptr structure */
  INT       **cdofs;                    /* list of coupled dofs and there owners, see init_assembly */
  INT         ncdofs;                   /* total number of coupled dofs */
  INT       **isend1;                   /* pointer to sendbuffer to communicate coupling conditions */
  DOUBLE    **dsend1;                   /* pointer to sendbuffer to communicate coupling conditions */
  INT       **isend2;                   /* pointer to sendbuffer to communicate coupling conditions */
  DOUBLE    **dsend2;                   /* pointer to sendbuffer to communicate coupling conditions */
  INT         nsend;

  struct _ARRAY ele_index;
  struct _ARRAY ele_locm;
#ifdef PARALLEL
  struct _ARRAY ele_owner;
#endif


#ifdef DEBUG 
  dstrc_enter("spo_make_index");
#endif

  /* set some pointers and variables */
  myrank     = actintra->intra_rank;
  nprocs     = actintra->intra_nprocs;
  nnz        = spo1->nnz;
  numeq_total= spo1->numeq_total;
  numeq      = spo1->numeq;
  update     = spo1->update.a.iv;
  A_loc      = spo1->A_loc.a.dv;
  irn        = spo1->irn_loc.a.iv;
  jcn        = spo1->jcn_loc.a.iv;
  rowptr     = spo1->rowptr.a.iv;
  cdofs      = actpart->pdis[0].coupledofs.a.ia;
  ncdofs     = actpart->pdis[0].coupledofs.fdim;

  /* put pointers to sendbuffers if any */
#ifdef PARALLEL 
  if (spo1->couple_i_send) 
  {
    isend1 = spo1->couple_i_send->a.ia;
    dsend1 = spo1->couple_d_send->a.da;
    nsend  = spo1->couple_i_send->fdim;
  }
#endif


  /* determine the size of estiff */
  counter=0;
  for (i=0; i<actele->numnp; i++)
  {
    for (j=0; j<actele->node[i]->numdf; j++)
    {
      counter++;
    }
  }
  /* end of loop over element nodes */
  nd = counter;
  actele->nd = counter;


  /* allocate locm, index and owner */
  actele->locm  = amdef("locm" ,&ele_locm ,nd, 1,"IV");
  actele->index = amdef("index",&ele_index,nd,nd,"IA");
#ifdef PARALLEL
  actele->owner = amdef("owner",&ele_owner,nd, 1,"IV");
#endif


  /* make location vector locm */
  counter=0;
  for (i=0; i<actele->numnp; i++)
  {
    for (j=0; j<actele->node[i]->numdf; j++)
    {
      actele->locm[counter]    = actele->node[i]->dof[j];
#ifdef PARALLEL 
      actele->owner[counter]   = actele->node[i]->proc;
#endif
      counter++;
    }
  }
  /* end of loop over element nodes */



  /* now start looping the dofs */
  /* loop over i (the element row) */
  ii_iscouple = 0;
  ii_owner    = myrank;

  for (i=0; i<nd; i++)
  {
    ii = actele->locm[i];

    /* loop only my own rows */
#ifdef PARALLEL 
    if (actele->owner[i]!=myrank)
    {
      for (j=0; j<nd; j++) actele->index[i][j] = -1;
      continue;
    }
#endif

    /* check for boundary condition */
    if (ii>=numeq_total)
    {
      for (j=0; j<nd; j++) actele->index[i][j] = -1;
      continue;
    }

    /* check for coupling condition */
#ifdef PARALLEL 
    if (ncdofs)
    {
      ii_iscouple = 0;
      ii_owner    = -1;
      add_msr_checkcouple(ii,cdofs,ncdofs,&ii_iscouple,&ii_owner,nprocs);
    }
#endif

    /* ii is not a coupled dofs or I am master owner */
    ii_index      = find_index(ii,update,numeq);

    if (!ii_iscouple || ii_owner==myrank)
    {

      if (ii_index==-1) dserror("dof %4i not found on this proc",ii);
      start         = rowptr[ii_index];
      lenght        = rowptr[ii_index+1]-rowptr[ii_index];

    }
    /* loop over j (the element column) */
    for (j=0; j<nd; j++)
    {
      jj = actele->locm[j];

      /* check for boundary condition */
      if (jj>=numeq_total)
      {
        actele->index[i][j] = -1;
        continue;
      }

      /* do main-diagonal entry */
      /* (either not a coupled dof or I am master owner) */
      if (!ii_iscouple || ii_owner==myrank)
      {
        index         = find_index(jj,&(jcn[start]),lenght);
        if (index==-1) dserror("dof jj not found in this row ii");
        index        += start;
        actele->index[i][j] = index;
      }

      /* do main-diagonal entry */
      /* (a coupled dof and I am slave owner) */
      else
      {
        actele->index[i][j] = -2;
      }

    } /* end loop over j */
  }/* end loop over i */


#ifdef DEBUG 
  dstrc_exit();
#endif

  return;
} /* end of spo_make_index */


#endif /* ifdef SPOOLES_PACKAGE */

