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


static INT  disnum;


/*----------------------------------------------------------------------*
  |  calculate the mask of an msr matrix                  m.gee 5/01   |
 *----------------------------------------------------------------------*/
void mask_msr(
    FIELD         *actfield,
    PARTITION     *actpart,
    SOLVAR        *actsolv,
    INTRA         *actintra,
    AZ_ARRAY_MSR  *msr,
    INT            disnum_
    )

{
  INT       i;
  INT       numeq;
  INT     **dof_connect;

#ifdef FAST_ASS
  ELEMENT  *actele;
#endif

#ifdef DEBUG
  dstrc_enter("mask_msr");
#endif


  disnum = disnum_;

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
  msr->numeq_total = actfield->dis[disnum].numeq;
  /* count number of eqns on proc and build processor-global couplingdof
     matrix */
  mask_numeq(actfield,actpart,actsolv,actintra,&numeq,disnum);
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

#ifdef FAST_ASS
  /* make the index vector for faster assembling */
  for (i=0; i<actpart->pdis[disnum].numele; i++)
  {
    actele = actpart->pdis[disnum].element[i];
    msr_make_index(actfield,actpart,actintra,actele,msr);
  }
#endif


#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of mask_msr */





/*----------------------------------------------------------------------*
  |  allocate update put dofs in update in ascending order  m.gee 5/01 |
 *----------------------------------------------------------------------*/
void msr_update(
    FIELD         *actfield,
    PARTITION     *actpart,
    SOLVAR        *actsolv,
    INTRA         *actintra,
    AZ_ARRAY_MSR  *msr
    )

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
  memset(&coupledofs, 0, sizeof(ARRAY));
  if (actpart->pdis[disnum].coupledofs.Typ != cca_XX)
    am_alloc_copy(&(actpart->pdis[disnum].coupledofs),&coupledofs);
  /*----------------------------------------------------------------------*/
  update = msr->update.a.iv;
  counter=0;
  /*------------------------------------- loop the nodes on the partition */
  for (i=0; i<actpart->pdis[disnum].numnp; i++)
  {
    actnode = actpart->pdis[disnum].node[i];
    for (l=0; l<actnode->numdf; l++)
    {
      dof = actnode->dof[l];
      /* dirichlet condition on dof */
      if (dof >= actfield->dis[disnum].numeq) continue;
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
  if (coupledofs.fdim > 0)
    amdel(&coupledofs);
  /*----------------------------------------------------------------------*/

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of msr_update */





/*----------------------------------------------------------------------*
  | calculate number of nonzero entries and dof topology   m.gee 6/01  |
 *----------------------------------------------------------------------*/
void msr_nnz_topology(
    FIELD         *actfield,
    PARTITION     *actpart,
    SOLVAR        *actsolv,
    INTRA         *actintra,
    AZ_ARRAY_MSR  *msr,
    INT          **dof_connect
    )

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
  for (j=0; j<actpart->pdis[disnum].numnp; j++)
  {
    for (k=0; k<actpart->pdis[disnum].node[j]->numdf; k++)
    {
      if (actpart->pdis[disnum].node[j]->dof[k] >= max_dof)
      {
        max_dof = actpart->pdis[disnum].node[j]->dof[k];
      }
    }
  }
  /* allocate pointer vector to the nodes */
  node_dof = (NODE**)CCACALLOC(max_dof+1,sizeof(NODE*));
  /* store pointers to nodes in node_dof at position accord. to dof */
  for (j=0; j<actpart->pdis[disnum].numnp; j++)
  {
    for (k=0; k<actpart->pdis[disnum].node[j]->numdf; k++)
    {
      dof_id = actpart->pdis[disnum].node[j]->dof[k];
      dsassert(dof_id <= max_dof,"zu kleiner node_dof Vector");
      node_dof[dof_id] = actpart->pdis[disnum].node[j];
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
          if (actnode->dof[l] < actfield->dis[disnum].numeq)
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

  /* free node_dof */
  CCAFREE(node_dof);

  /*--------------------------------------------- now do the coupled dofs */
  coupledofs = &(actpart->pdis[disnum].coupledofs);
  for (i=0; i<coupledofs->fdim; i++)
  {
    dof = coupledofs->a.ia[i][0];
    /*--------------------------- check for my own ownership of this dof */
    dofflag = coupledofs->a.ia[i][imyrank+1];
    /*----------- if dofflag is zero this dof has nothing to do with me */
    if (dofflag==0) continue;
    /*------------------------------------- find all patches to this dof */
    counter=0;
    for (j=0; j<actpart->pdis[disnum].numnp; j++)
    {
      centernode=NULL;
      for (l=0; l<actpart->pdis[disnum].node[j]->numdf; l++)
      {
        if (dof == actpart->pdis[disnum].node[j]->dof[l])
        {
          centernode = actpart->pdis[disnum].node[j];
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
              if (actnode->dof[l] < actfield->dis[disnum].numeq)
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
  |  make the DMSR vector bindx                            m.gee 6/01  |
  | for format see Aztec manual                                        |
 *----------------------------------------------------------------------*/
void msr_make_bindx(
    FIELD         *actfield,
    PARTITION     *actpart,
    SOLVAR        *actsolv,
    AZ_ARRAY_MSR  *msr,
    INT          **dof_connect
    )

{
  INT        i,j;
  INT        count1,count2;
  INT        dof;

#ifdef DEBUG
  dstrc_enter("msr_make_bindx");
#endif

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




#ifdef FAST_ASS
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

#endif  /* ifdef FAST_ASS */



/*----------------------------------------------------------------------*/
/*!
  \brief assemble into a msr matrix (original version)

  This routine assembles one or two element matrices (estiff_global and
  emass_global) into the global matrices in the msr format.
  It makes extensive use of the search algorithms provided by AZTEC.

  \param actpart   *PARTITION    (i)  the partition we are working on
  \param actsolv   *SOLVAR       (i)  the solver we are using
  \param actintra  *INTRA        (i)  the intra-communicator we do not need
  \param actele    *ELEMENT      (i)  the element we would like to work with
  \param msr1      *AZ_ARRAY_MSR (i)  one sparse matrix we will assemble into
  \param msr2      *AZ_ARRAY_MSR (i)  the other sparse matrix we will assemble into

  \author mn
  \date 07/04

*/
/*----------------------------------------------------------------------*/
void  add_msr(
    struct _PARTITION     *actpart,
    struct _SOLVAR        *actsolv,
    struct _INTRA         *actintra,
    struct _ELEMENT       *actele,
    struct _AZ_ARRAY_MSR  *msr1,
    struct _AZ_ARRAY_MSR  *msr2,
    struct _ARRAY         *elearray1,
    struct _ARRAY         *elearray2
    )

{

  INT         i,j,counter;          /* some counter variables */
  INT         start,index,lenght;   /* some more special-purpose counters */
  INT         ii,jj;                /* counter variables for system matrix */
  INT         ii_iscouple;          /* flag whether ii is a coupled dof */
  INT         ii_owner;             /* who is owner of dof ii -> procnumber */
  INT         ii_index;             /* place of ii in dmsr format */
  INT         nd;                   /* size of estif */
  INT         nnz;                  /* number of nonzeros in sparse system matrix */
  INT         numeq_total;          /* total number of equations */
  INT         numeq;                /* number of equations on this proc */
  INT         lm[MAXDOFPERELE];     /* location vector for this element */
#ifdef PARALLEL
  INT         owner[MAXDOFPERELE];  /* the owner of every dof */
#endif
  INT         myrank;               /* my intra-proc number */
  INT         nprocs;               /* my intra- number of processes */
  DOUBLE    **estif;                /* element matrix to be added to system matrix */
  DOUBLE    **emass;                /* element matrix to be added to system matrix */
  INT        *update;               /* msr-vector update see AZTEC manual */
  INT         shift;                /* variables for aztec quick finding algorithms */
  INT        *bins;
  INT        *bindx;                /*    "       bindx         "         */
  DOUBLE     *val1,*val2;           /*    "       val           "         */
  INT       **cdofs;                /* list of coupled dofs and there owners */
  INT         ncdofs;               /* total number of coupled dofs */
  INT       **isend1 = NULL;        /* p to sendbuffer to communicate coupling cond */
  DOUBLE    **dsend1 = NULL;        /* p to sendbuffer to communicate coupling cond */
  INT       **isend2 = NULL;        /* p to sendbuffer to communicate coupling cond */
  DOUBLE    **dsend2 = NULL;        /* p to sendbuffer to communicate coupling cond */
  INT         nsend =0;


#ifdef DEBUG
  dstrc_enter("add_msr");
#endif

  /* set some pointers and variables */
  myrank     = actintra->intra_rank;
  nprocs     = actintra->intra_nprocs;
  estif      = elearray1->a.da;
  if (msr2) emass = elearray2->a.da;
  else      emass = NULL;
  nd         = actele->numnp * actele->node[0]->numdf;
  nnz        = msr1->nnz;
  numeq_total= msr1->numeq_total;
  numeq      = msr1->numeq;
  update     = msr1->update.a.iv;
  bindx      = msr1->bindx.a.iv;
  val1       = msr1->val.a.dv;
  if (msr2) val2 = msr2->val.a.dv;
  else      val2 = NULL;
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


  /* put pointers to sendbuffers if any */
#ifdef PARALLEL
  if (msr1->couple_i_send)
  {
    isend1 = msr1->couple_i_send->a.ia;
    dsend1 = msr1->couple_d_send->a.da;
    nsend  = msr1->couple_i_send->fdim;
    if (msr2)
    {
      isend2 = msr2->couple_i_send->a.ia;
      dsend2 = msr2->couple_d_send->a.da;
    }
  }
#endif

  /* make location vector lm*/
  counter=0;
  for (i=0; i<actele->numnp; i++)
  {
    for (j=0; j<actele->node[i]->numdf; j++)
    {
      lm[counter]    = actele->node[i]->dof[j];
#ifdef PARALLEL
      owner[counter] = actele->node[i]->proc;
#endif
      counter++;
    }
  }
  /* end of loop over element nodes */
  /* this check is not possible any more for fluid element with implicit
     free surface condition: nd not eqaual numnp*numdf!!! */
#if 0
  if (counter != nd) dserror("assemblage failed due to wrong dof numbering");
#endif
  nd = counter;


  /* now start looping the dofs */
  /* loop over i (the element row) */
  ii_iscouple = 0;
  ii_owner    = myrank;
  for (i=0; i<nd; i++)
  {
    ii = lm[i];
    /* loop only my own rows */
#ifdef PARALLEL
    if (owner[i]!=myrank) continue;
#endif

    /* check for boundary condition */
    if (ii>=numeq_total) continue;

    /* check for coupling condition */
#ifdef PARALLEL
    if (ncdofs)
    {
      ii_iscouple = 0;
      ii_owner    = -1;
      add_msr_checkcouple(ii,cdofs,ncdofs,&ii_iscouple,&ii_owner,nprocs);
    }
#endif


    ii_index = AZ_quick_find(ii,update,numeq,shift,bins);
    if (ii_index==-1) dserror("dof ii not found on this proc");

    start       = bindx[ii_index];
    lenght      = bindx[ii_index+1]-bindx[ii_index];

    /* loop over j (the element column) */
    /* This is the full unsymmetric version ! */
    for (j=0; j<nd; j++)
    {
      jj = lm[j];

      /* check for boundary condition */
      if (jj>=numeq_total) continue;

      /* do main-diagonal entry */
      /* (either not a coupled dof or I am master owner) */
      if (!ii_iscouple || ii_owner==myrank)
      {

        /* if (i==j) AL (coupling of several dofs in the same element) */
        if (ii==jj)
        {
          val1[ii_index] += estif[i][j];
          if (msr2)
            val2[ii_index] += emass[i][j];
        }

        /* do off-diagonal entry in row ii */
        /* (either not a coupled dof or I am master owner) */
        else
        {
          index       = AZ_find_index(jj,&(bindx[start]),lenght);
          if (index==-1) dserror("dof jj not found in this row ii");

          index      += start;

          val1[index] += estif[i][j];
          if (msr2)
            val2[index] += emass[i][j];
        }
      }

      /* do main-diagonal entry */
      /* (a coupled dof and I am slave owner) */
      else
      {
        add_msr_sendbuff(ii,jj,i,j,ii_owner,isend1,dsend1,estif,nsend);
        if (msr2)
          add_msr_sendbuff(ii,jj,i,j,ii_owner,isend2,dsend2,emass,nsend);
      }
    } /* end loop over j */
  }/* end loop over i */


#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of add_msr */



#ifdef FAST_ASS
/*----------------------------------------------------------------------*/
/*!
  \brief assemble into a msr matrix (faster version)

  This routine assembles one or two element matrices (estiff_global and
  emass_global) into the global matrices in the msr format. (same as add_msr)
  It makes use of the information saved in actele->index to determine the
  correct position in the sparse matrix.
  This is faster then searching every time, but consumes a lot of memory!!

  \param actpart   *PARTITION    (i)  the partition we are working on
  \param actsolv   *SOLVAR       (i)  the solver we are using
  \param actintra  *INTRA        (i)  the intra-communicator we do not need
  \param actele    *ELEMENT      (i)  the element we would like to work with
  \param msr1      *AZ_ARRAY_MSR (i)  one sparse matrix we will assemble into
  \param msr2      *AZ_ARRAY_MSR (i)  the other sparse matrix we will assemble into

  \author mn
  \date 07/04

*/
/*----------------------------------------------------------------------*/
void  add_msr_fast(
    struct _PARTITION     *actpart,
    struct _SOLVAR        *actsolv,
    struct _INTRA         *actintra,
    struct _ELEMENT       *actele,
    struct _AZ_ARRAY_MSR  *msr1,
    struct _AZ_ARRAY_MSR  *msr2,
    struct _ARRAY         *elearray1,
    struct _ARRAY         *elearray2
    )

{

  INT         i,j;           /* some counter variables */
  INT         index;         /* some more special-purpose counters */
  INT         ii,jj;         /* counter variables for system matrix */
  INT         ii_owner;      /* who is owner of dof ii -> procnumber */
  INT         nd;            /* size of estif */
  INT         myrank;        /* my intra-proc number */
  INT         nprocs;        /* my intra- number of processes */
  DOUBLE    **estif;         /* element matrix to be added to system matrix */
  DOUBLE    **emass;         /* element matrix to be added to system matrix */
  DOUBLE     *val1,*val2;    /*    "       val           "         */
  INT       **isend1;        /* pointer to sendbuffer to communicate coupling conditions */
  DOUBLE    **dsend1;        /* pointer to sendbuffer to communicate coupling conditions */
  INT       **isend2;        /* pointer to sendbuffer to communicate coupling conditions */
  DOUBLE    **dsend2;        /* pointer to sendbuffer to communicate coupling conditions */
  INT         nsend;

#ifdef DEBUG
  dstrc_enter("add_msr_fast");
#endif

  /* set some pointers and variables */
  myrank     = actintra->intra_rank;
  nprocs     = actintra->intra_nprocs;
  estif      = elearray1->a.da;
  if (msr2) emass = elearray2->a.da;
  else      emass = NULL;
  nd         = actele->nd;
  val1       = msr1->val.a.dv;
  if (msr2) val2 = msr2->val.a.dv;
  else      val2 = NULL;


  /* put pointers to sendbuffers if any */
#ifdef PARALLEL
  if (msr1->couple_i_send)
  {
    isend1 = msr1->couple_i_send->a.ia;
    dsend1 = msr1->couple_d_send->a.da;
    nsend  = msr1->couple_i_send->fdim;
    if (msr2)
    {
      isend2 = msr2->couple_i_send->a.ia;
      dsend2 = msr2->couple_d_send->a.da;
    }
  }
#endif

  /* loop over i (the element row) */
  ii_owner    = myrank;
  for (i=0; i<nd; i++)
  {
    ii = actele->locm[i];

    /* loop over j (the element column) */
    for (j=0; j<nd; j++)
    {
      jj = actele->locm[j];
      index = actele->index[i][j];

      if(index >= 0)  /* normal dof */
      {
        val1[index] += estif[i][j];
        if (msr2)
          val2[index] += emass[i][j];
      }

      if(index == -1)  /* boundary condition dof */
        continue;

      if(index == -2)  /* coupled dof */
      {
        add_msr_sendbuff(ii,jj,i,j,ii_owner,isend1,dsend1,estif,nsend);
        if (msr2)
          add_msr_sendbuff(ii,jj,i,j,ii_owner,isend2,dsend2,emass,nsend);
      }
    } /* end loop over j */
  }/* end loop over i */


#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of add_msr_fast */

#endif /* ifdef FAST_ASS */




#ifdef FAST_ASS2
/*----------------------------------------------------------------------*/
/*!
  \brief assemble into a msr matrix (faster version 2)

  This routine assembles one or two element matrices (estiff_global and
  emass_global) into the global matrices in the msr format. (same as add_msr)
  It makes use of "inversion" of update and parts of bindx to avoid the
  searching.
  This is faster then the original, and consumes only very little memory.

  \param actpart   *PARTITION    (i)  the partition we are working on
  \param actsolv   *SOLVAR       (i)  the solver we are using
  \param actintra  *INTRA        (i)  the intra-communicator we do not need
  \param actele    *ELEMENT      (i)  the element we would like to work with
  \param msr1      *AZ_ARRAY_MSR (i)  one sparse matrix we will assemble into
  \param msr2      *AZ_ARRAY_MSR (i)  the other sparse matrix we will assemble into

  \author mn
  \date 07/04

*/
/*----------------------------------------------------------------------*/
void  add_msr_fast2(
    struct _PARTITION     *actpart,
    struct _SOLVAR        *actsolv,
    struct _INTRA         *actintra,
    struct _ELEMENT       *actele,
    struct _AZ_ARRAY_MSR  *msr1,
    struct _AZ_ARRAY_MSR  *msr2,
    struct _ARRAY         *elearray1,
    struct _ARRAY         *elearray2
    )

{

  INT         i,j,k,counter;        /* some counter variables */
  INT         start,index,lenght;   /* some more special-purpose counters */
  INT         ii,jj;                /* counter variables for system matrix */
  INT         ii_iscouple;          /* flag whether ii is a coupled dof */
  INT         ii_owner;             /* who is owner of dof ii -> procnumber */
  INT         ii_index;             /* place of ii in dmsr format */
  INT         nd;                   /* size of estif */
  INT         numeq_total;          /* total number of equations */
  INT         numeq;                /* number of equations on this proc */
  INT         lm[MAXDOFPERELE];     /* location vector for this element */
#ifdef PARALLEL
  INT         owner[MAXDOFPERELE];  /* the owner of every dof */
#endif
  INT         myrank;               /* my intra-proc number */
  INT         nprocs;               /* my intra- number of processes */
  DOUBLE    **estif;                /* element matrix to be added to system matrix */
  DOUBLE    **emass;                /* element matrix to be added to system matrix */
  INT        *update;               /* msr-vector update see AZTEC manual */
  INT         shift;                /* variables for aztec quick finding algorithms */
  INT        *bins;
  INT        *bindx;                /*    "       bindx         "         */
  DOUBLE     *val1,*val2;           /*    "       val           "         */
  INT       **cdofs;                /* list of coupled dofs and there owners */
  INT         ncdofs;               /* total number of coupled dofs */
  INT       **isend1 = NULL;        /* p to sendbuffer to communicate coupling cond */
  DOUBLE    **dsend1 = NULL;        /* p to sendbuffer to communicate coupling cond */
  INT       **isend2 = NULL;        /* p to sendbuffer to communicate coupling cond */
  DOUBLE    **dsend2 = NULL;        /* p to sendbuffer to communicate coupling cond */
  INT         nsend =0;

  INT        *invupdate = NULL;
  INT        *invbindx  = NULL;
  INT index2;

#ifdef DEBUG
  dstrc_enter("add_msr_fast2");
#endif

  /* set some pointers and variables */
  myrank     = actintra->intra_rank;
  nprocs     = actintra->intra_nprocs;
  estif      = elearray1->a.da;
  if (msr2) emass = elearray2->a.da;
  else      emass = NULL;
  nd         = actele->numnp * actele->node[0]->numdf;
  numeq_total= msr1->numeq_total;
  numeq      = msr1->numeq;
  update     = msr1->update.a.iv;
  bindx      = msr1->bindx.a.iv;
  val1       = msr1->val.a.dv;
  invupdate  = msr1->invupdate;
  invbindx   = msr1->invbindx;
  if (msr2) val2 = msr2->val.a.dv;
  else      val2 = NULL;
  cdofs      = actpart->pdis[0].coupledofs.a.ia;
  ncdofs     = actpart->pdis[0].coupledofs.fdim;


  if (!invupdate)
  {
    /* allocate invupdate */
    invupdate = (INT*)CCACALLOC( numeq_total,sizeof(INT));
    if (!invupdate) dserror("Allocation of invupdate failed");

    /* initialize with minus one */
    for (k=0; k<numeq_total; k++)
    {
      invupdate[k] = -1;
    }

    /* fill invupdate */
    for (k=0; k<numeq; k++)
    {
      invupdate[update[k]] = k;
    }
  }


  if (!invbindx)
  {
    /* allocate invbindx */
    invbindx = (INT*)CCACALLOC( numeq_total,sizeof(INT));
    if (!invbindx) dserror("Allocation of invbindx failed");
  }


  /* put pointers to sendbuffers if any */
#ifdef PARALLEL
  if (msr1->couple_i_send)
  {
    isend1 = msr1->couple_i_send->a.ia;
    dsend1 = msr1->couple_d_send->a.da;
    nsend  = msr1->couple_i_send->fdim;
    if (msr2)
    {
      isend2 = msr2->couple_i_send->a.ia;
      dsend2 = msr2->couple_d_send->a.da;
    }
  }
#endif

  /* make location vector lm*/
  counter=0;
  for (i=0; i<actele->numnp; i++)
  {
    for (j=0; j<actele->node[i]->numdf; j++)
    {
      lm[counter]    = actele->node[i]->dof[j];
#ifdef PARALLEL
      owner[counter] = actele->node[i]->proc;
#endif
      counter++;
    }
  }
  /* end of loop over element nodes */
  nd = counter;


  /* now start looping the dofs */
  /* loop over i (the element row) */
  ii_iscouple = 0;
  ii_owner    = myrank;
  for (i=0; i<nd; i++)
  {
    ii = lm[i];
    /* loop only my own rows */
#ifdef PARALLEL
    if (owner[i]!=myrank) continue;
#endif

    /* check for boundary condition */
    if (ii>=numeq_total) continue;

    /* check for coupling condition */
#ifdef PARALLEL
    if (ncdofs)
    {
      ii_iscouple = 0;
      ii_owner    = -1;
      add_msr_checkcouple(ii,cdofs,ncdofs,&ii_iscouple,&ii_owner,nprocs);
    }
#endif


    ii_index = invupdate[ii];
    if (ii_index==-1) dserror("dof ii not found on this proc");


    /* NO initialization: only values that are set below will be used for this row!! */
    /* fill invbindx */
    for (k=bindx[ii_index]; k<bindx[ii_index+1]; k++)
    {
      invbindx[bindx[k]] = k;
    }

    /* loop over j (the element column) */
    /* This is the full unsymmetric version ! */
    for (j=0; j<nd; j++)
    {
      jj = lm[j];

      /* check for boundary condition */
      if (jj>=numeq_total) continue;

      /* do main-diagonal entry */
      /* (either not a coupled dof or I am master owner) */
      if (!ii_iscouple || ii_owner==myrank)
      {
        /* if (i==j) AL (coupling of several dofs in the same element) */
        if (ii==jj)
        {
          val1[ii_index] += estif[i][j];
          if (msr2)
            val2[ii_index] += emass[i][j];
        }

        /* do off-diagonal entry in row ii */
        /* (either not a coupled dof or I am master owner) */
        else
        {
          index = invbindx[jj];
          val1[index] += estif[i][j];
          if (msr2)
            val2[index] += emass[i][j];
        }
      }

      /* do main-diagonal entry */
      /* (a coupled dof and I am slave owner) */
      else
      {
        add_msr_sendbuff(ii,jj,i,j,ii_owner,isend1,dsend1,estif,nsend);
        if (msr2)
          add_msr_sendbuff(ii,jj,i,j,ii_owner,isend2,dsend2,emass,nsend);
      }
    } /* end loop over j */
  }/* end loop over i */


#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of add_msr_fast2 */

#endif /* ifdef FAST_ASS2 */





/*----------------------------------------------------------------------*
  |  checks coupling for the add_msr routine                 m.gee 9/01|
 *----------------------------------------------------------------------*/
void add_msr_checkcouple(
    INT       ii,
    INT     **cdofs,
    INT       ncdofs,
    INT      *iscouple,
    INT      *isowner,
    INT       nprocs
    )

{
  INT         i,k;

#ifdef DEBUG
  dstrc_enter("add_msr_checkcouple");
#endif

  /*----------------------------------------------------------------------*/
  for (k=0; k<ncdofs; k++)
  {
    if (ii==cdofs[k][0])
    {
      *iscouple=1;
      for (i=1; i<=nprocs; i++)
      {
        if (cdofs[k][i]==2)
        {
          *isowner=i-1;
          break;
        }
      }
    }
  }
  /*----------------------------------------------------------------------*/

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of add_msr_checkcouple */



/*----------------------------------------------------------------------*
  |  fill sendbuffer isend and dsend                         m.gee 9/01|
 *----------------------------------------------------------------------*/
void add_msr_sendbuff(
    INT         ii,
    INT         jj,
    INT         i,
    INT         j,
    INT         ii_owner,
    INT       **isend,
    DOUBLE    **dsend,
    DOUBLE    **estif,
    INT         numsend
    )

{
  INT         k;

#ifdef DEBUG
  dstrc_enter("add_msr_sendbuff");
#endif

  /*----------------------------------------------------------------------*/
  for (k=0; k<numsend; k++)
  {
    if (isend[k][0]==ii) break;
  }
  isend[k][1]  = ii_owner;
  dsend[k][jj]+= estif[i][j];
  /*----------------------------------------------------------------------*/

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of add_msr_sendbuff */





/*----------------------------------------------------------------------*
  |  exchange coupled dofs and add to dmsr matrix            m.gee 9/01|
 *----------------------------------------------------------------------*/
void exchange_coup_msr(
    PARTITION     *actpart,
    SOLVAR        *actsolv,
    INTRA         *actintra,
    AZ_ARRAY_MSR  *msr
    )

{

#ifdef PARALLEL
  INT            i,j;
  INT            ii,ii_index;
  INT            start;
  INT            lenght;
  INT            index;
  INT            tag;
  INT            source;
  INT            numeq,numeq_total;
  INT            numsend;
  INT            numrecv;
  INT            shift;
  INT           *bins;
  INT           *bindx;
  DOUBLE        *val;
  INT           *update;
  INT          **isend;
  DOUBLE       **dsend;
  INT          **irecv;
  DOUBLE       **drecv;
  INT            imyrank;
  INT            inprocs;

  MPI_Status    *irecv_status;
  MPI_Status    *drecv_status;

  MPI_Request   *isendrequest;
  MPI_Request   *dsendrequest;

  MPI_Comm      *ACTCOMM;
#endif

#ifdef DEBUG
  dstrc_enter("exchange_coup_msr");
#endif

  /*----------------------------------------------------------------------*/
#ifdef PARALLEL
  /*----------------------------------------------------------------------*/
  imyrank = actintra->intra_rank;
  inprocs = actintra->intra_nprocs;
  ACTCOMM = &(actintra->MPI_INTRA_COMM);
  /*---------------------------------------- set some pointers and values */
  numsend     = msr->numcoupsend;
  numrecv     = msr->numcouprecv;
  shift       = msr->shift;
  bins        = msr->bins;
  bindx       = msr->bindx.a.iv;
  val         = msr->val.a.dv;
  update      = msr->update.a.iv;
  numeq_total = msr->numeq_total;
  numeq       = msr->numeq;
  if (msr->couple_i_send) isend   = msr->couple_i_send->a.ia;
  if (msr->couple_d_send) dsend   = msr->couple_d_send->a.da;
  if (msr->couple_i_recv) irecv   = msr->couple_i_recv->a.ia;
  if (msr->couple_d_recv) drecv   = msr->couple_d_recv->a.da;
  /*--------------------------------------------- allocate some envelopes */
  if (numrecv)
  {
    irecv_status = (MPI_Status*)CCACALLOC(numrecv,sizeof(MPI_Status));
    drecv_status = (MPI_Status*)CCACALLOC(numrecv,sizeof(MPI_Status));
    if (!irecv_status || !drecv_status) dserror("Allocation of memory failed");
  }
  if (numsend)
  {
    isendrequest = (MPI_Request*)CCACALLOC(numsend,sizeof(MPI_Request));
    dsendrequest = (MPI_Request*)CCACALLOC(numsend,sizeof(MPI_Request));
    if ( !isendrequest || !dsendrequest) dserror("Allocation of memory failed");
  }
  /*-------------------------------------------- loop the dofs to be send */
  /* do all non-blocking sends and don't care about arrival (wird scho' klappe)*/
  /*     the only thing to care for is the order in which things are send */
  for (i=0; i<numsend; i++)
  {
    /*            sendbuffer       lenght    typ        dest.        tag          comm      request-handle */
    MPI_Isend(&(isend[i][0]),          2,MPI_INT   ,isend[i][1],isend[i][0],(*ACTCOMM),&(isendrequest[i]));
    MPI_Isend(&(dsend[i][0]),numeq_total,MPI_DOUBLE,isend[i][1],isend[i][0],(*ACTCOMM),&(dsendrequest[i]));
  }/*------------------------------------------------ end of sending loop */
  /*------------------------------- now loop over the dofs to be received */
  /*
     do blocking receives, 'cause one can't add something to the system
     matrix, which has not yet arrived, easy, isn't it?
     */
  for (i=0; i<numrecv; i++)
  {
    /*--------------------------- use wildcards to receive first to come */
    /*          recv-buf  lenght typ     source           tag       comm              status */
    MPI_Recv(&(irecv[i][0]),2,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,(*ACTCOMM),&(irecv_status[i]));
    if (irecv_status[i].MPI_ERROR) dserror("An error in MPI - communication occured !");

    /*---------------------- the dof number was sent as tag and as entry */
    tag    = irecv_status[i].MPI_TAG;
    if (tag != irecv[i][0]) dserror("MPI messages somehow got mixed up");
    source = irecv_status[i].MPI_SOURCE;

    /* do not use wildcards for second recv, we know now where it should come from */
    MPI_Recv(&(drecv[i][0]),numeq_total,MPI_DOUBLE,source,tag,(*ACTCOMM),&(drecv_status[i]));
    if (drecv_status[i].MPI_ERROR) dserror("An error in MPI - communication occured !");

    /* now add the received data properly to my own piece of sparse matrix */
    ii = tag;
    ii_index = AZ_quick_find(ii,update,numeq,shift,bins);
    if (ii_index==-1) dserror("dof ii not found on this proc");
    /*---------------------------------------------- main diagonal entry */
    val[ii_index] += drecv[i][ii];
    /*--------------------------------------------- off-diagonal entries */
    start  = bindx[ii_index];
    lenght = bindx[ii_index+1]-bindx[ii_index];
    for (j=0; j<lenght; j++)
    {
      index         = bindx[start+j];
      val[start+j] += drecv[i][index];
    }
  }/*---------------------------------------------- end of receiving loop */
  /*-------------------------------------------- free allocated MPI-stuff */
  if (numrecv){CCAFREE(irecv_status);CCAFREE(drecv_status);}
  if (numsend){CCAFREE(isendrequest);CCAFREE(dsendrequest);}
  /*----------------------------------------------------------------------
    do a barrier, because this is the end of the assembly, the msr matrix
    is now ready for solve
    */
  MPI_Barrier(*ACTCOMM);
#endif /*---------------------------------------------- end of PARALLEL */
  /*----------------------------------------------------------------------*/

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of exchange_coup_msr */



#endif /* ifdef AZTEC_PACKAGE */


