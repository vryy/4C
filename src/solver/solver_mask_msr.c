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

#ifdef AZTEC_PACKAGE

#include "../headers/standardtypes.h"
#include "../solver/solver.h"
INT cmp_int(const void *a, const void *b );
DOUBLE cmp_double(const void *a, const void *b );

#if 0
#ifdef DEBUG

static void out_ivector(INTRA* actintra, INT* vector, INT size, CHAR* vname)
{
  INT i;
  FILE* f;
  CHAR name[256];
/*   static INT count = 0; */

  sprintf(name, "ivector.%s.%d", vname, actintra->intra_rank);
  f = fopen(name, "w");
  for (i=0; i<size; ++i) {
    fprintf(f, "%d: %d\n", i, vector[i]);
  }
  fclose(f);
  
/*   count++; */
}

#endif
#endif

static INT kk;
/*----------------------------------------------------------------------*
 |  calculate the mask of an msr matrix                  m.gee 5/01     |
 *----------------------------------------------------------------------*/
void mask_msr(FIELD         *actfield, 
              PARTITION     *actpart, 
              SOLVAR        *actsolv,
              INTRA         *actintra, 
              AZ_ARRAY_MSR  *msr,
	      INT            actndis)
{
INT       i;
INT       numeq;
INT     **dof_connect;

ELEMENT  *actele;

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
*/
/*------------------------------------------- put total size of problem */
msr->numeq_total = actfield->dis[kk].numeq;
/* count number of eqns on proc and build processor-global couplingdof 
                                                                 matrix */
mask_numeq(actfield,actpart,actsolv,actintra,&numeq,kk);
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
dof_connect = (INT**)CCACALLOC(msr->numeq_total,sizeof(INT*));
if (!dof_connect) dserror("Allocation of dof_connect failed");
msr_nnz_topology(actfield,actpart,actsolv,actintra,msr,dof_connect);
/*---------------------------------------------- allocate bindx and val */
amdef("bindx",&(msr->bindx),(msr->nnz+1),1,"IV");
amdef("val"  ,&(msr->val)  ,(msr->nnz+1),1,"DV");
msr->bindx_backup.Typ = cca_XX;
/*---------------------------------------------------------- make bindx */
msr_make_bindx(actfield,actpart,actsolv,msr,dof_connect);
/*---------------------------------------- delete the array dof_connect */
for (i=0; i<msr->numeq_total; i++)
{
   if (dof_connect[i]) CCAFREE(dof_connect[i]);
}
CCAFREE(dof_connect);

/* make the index vector for faster assembling */
  for (i=0; i<actpart->pdis[0].numele; i++)
  {
    actele = actpart->pdis[0].element[i];
    msr_make_index(actfield,actpart,actintra,actele,msr);
  }

#if 0
#ifdef DEBUG 
  out_ivector(actintra, msr->bindx.a.iv, msr->bindx.fdim, "bindx_global");
  out_ivector(actintra, msr->update.a.iv, msr->update.fdim, "update");
#endif
#endif

#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mask_msr */




/*----------------------------------------------------------------------*
 |  allocate update put dofs in update in ascending order   m.gee 5/01  |
 *----------------------------------------------------------------------*/
void msr_update(FIELD         *actfield, 
                PARTITION     *actpart, 
                SOLVAR        *actsolv,
                INTRA         *actintra,
                AZ_ARRAY_MSR  *msr)
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
qsort((INT*) update, counter, sizeof(INT), cmp_int);
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
                      INT          **dof_connect)
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

  /* create node pointers for all dofs in this partition */
  max_dof = 0;
  for (j=0; j<actpart->pdis[kk].numnp; j++)
  {
    for (k=0; k<actpart->pdis[kk].node[j]->numdf; k++)
    {
      if (actpart->pdis[kk].node[j]->dof[k] >= max_dof)
      {
        max_dof = actpart->pdis[kk].node[j]->dof[k];
      }
    }
  }
  /* allocate pointer vector to the nodes */
  node_dof = (NODE**)CCACALLOC(max_dof+1,sizeof(NODE*));
  /* store pointers to nodes in node_dof at position accord. to dof */
  for (j=0; j<actpart->pdis[kk].numnp; j++)
  {
    for (k=0; k<actpart->pdis[kk].node[j]->numdf; k++)
    {
      dof_id = actpart->pdis[kk].node[j]->dof[k];
      dsassert(dof_id <= max_dof,"zu kleiner node_dof Vector");
      node_dof[dof_id] = actpart->pdis[kk].node[j];
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
   centernode=NULL;
   centernode = node_dof[dof];
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
   qsort((INT*)(&(dof_connect[dof][3])), dof_connect[dof][0]-3, sizeof(INT), cmp_int);
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
 |  make the DMSR vector bindx                              m.gee 6/01  |
 | for format see Aztec manual                                          |
 *----------------------------------------------------------------------*/
void msr_make_bindx(FIELD         *actfield, 
                       PARTITION     *actpart, 
                       SOLVAR        *actsolv,
                       AZ_ARRAY_MSR  *msr,
                       INT          **dof_connect)
{
INT        i,j;
INT        count1,count2;
INT        dof;

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
 \param msr1      *AZ_ARRAY_MSR (i)  the sparse matrix we will assemble into

 \author mn
 \date 07/04

 */
/*----------------------------------------------------------------------*/
void msr_make_index(
    FIELD                 *actfield, 
    PARTITION             *actpart,
    INTRA                 *actintra,
    ELEMENT               *actele,
    struct _AZ_ARRAY_MSR  *msr1
    )
{

  INT         i,j,counter;          /* some counter variables */
  INT         start,index,lenght;       /* some more special-purpose counters */
  INT         ii,jj;                    /* counter variables for system matrix */
  INT         ii_iscouple;              /* flag whether ii is a coupled dof */
  INT         ii_owner;                 /* who is owner of dof ii -> procnumber */
  INT         ii_index;                 /* place of ii in dmsr format */
  INT         nd;                       /* size of estif */
  INT         numeq_total;              /* total number of equations */
  INT         numeq;                    /* number of equations on this proc */
  INT         myrank;                   /* my intra-proc number */
  INT         nprocs;                   /* my intra- number of processes */
  INT        *update;                   /* msr-vector update see AZTEC manual */
  INT         shift;                    /* variables for aztec quick finding algorithms */
  INT        *bins;
  INT        *bindx;                    /*    "       bindx         "         */
  INT       **cdofs;                    /* list of coupled dofs and there owners, see init_assembly */
  INT         ncdofs;                   /* total number of coupled dofs */

  struct _ARRAY ele_index;
  struct _ARRAY ele_locm;
#ifdef PARALLEL
  struct _ARRAY ele_owner;
#endif

#ifdef DEBUG 
  dstrc_enter("msr_make_index");
#endif

  /* set some pointers and variables */
  myrank     = actintra->intra_rank;
  nprocs     = actintra->intra_nprocs;
  numeq_total= msr1->numeq_total;
  numeq      = msr1->numeq;
  update     = msr1->update.a.iv;
  bindx      = msr1->bindx.a.iv;
  cdofs      = actpart->pdis[0].coupledofs.a.ia;
  ncdofs     = actpart->pdis[0].coupledofs.fdim;


  /* allocate and calculate shifts and bins for quick_find routines */
  if (!(msr1->bins))
  {
    msr1->bins = (INT*)CCACALLOC( ABS(4+numeq/4),sizeof(INT));
    if (!(msr1->bins)) dserror("Allocation of msr->bins failed");
    AZ_init_quick_find(update,numeq,&(msr1->shift),msr1->bins);
  }
  shift      = msr1->shift;
  bins       = msr1->bins;

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
        if (ii==jj)
        {
          ii_index = AZ_quick_find(ii,update,numeq,shift,bins);
          if (ii_index==-1) dserror("dof ii not found on this proc");
          actele->index[i][j] = ii_index;
        } 
        /* do off-diagonal entry in row ii */
        else
        {
          ii_index    = AZ_quick_find(ii,update,numeq,shift,bins);
          if (ii_index==-1) dserror("dof ii not found on this proc");
          start       = bindx[ii_index];
          lenght      = bindx[ii_index+1]-bindx[ii_index];
          index       = AZ_find_index(jj,&(bindx[start]),lenght);
          if (index==-1) dserror("dof jj not found in this row ii");
          index      += start;
          actele->index[i][j] = index;
        }
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
} /* end of msr_make_index */


#endif /* ifdef AZTEC_PACKAGE */

