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




static INT  disnum;



/*----------------------------------------------------------------------*
  | calculate the mask of an spooles matrix              m.gee 5/01    |
 *----------------------------------------------------------------------*/
void mask_spooles(
    FIELD         *actfield,
    PARTITION     *actpart,
    SOLVAR        *actsolv,
    INTRA         *actintra,
    SPOOLMAT      *spo,
    INT            disnum_
    )

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
  spo->numeq_total = actfield->dis[disnum].numeq;
  /* count number of eqns on proc and build processor-global couplingdof
     matrix */
  mask_numeq(actfield,actpart,actsolv,actintra,&numeq,disnum);
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


#ifdef FAST_ASS
  /* make the index vector for faster assembling */
  for (i=0; i<actpart->pdis[disnum].numele; i++)
  {
    actele = actpart->pdis[disnum].element[i];
    spo_make_index(actfield,actpart,actintra,actele,spo);
  }
#endif


#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of mask_spooles */





/*----------------------------------------------------------------------*
  | make the vectors                                        m.gee 1/02 |
  | irn_loc, jcn_loc, rowptr from update and bindx                     |
 *----------------------------------------------------------------------*/
void spo_make_sparsity(
    SPOOLMAT      *spo,
    INT           *bindx
    )

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
  | make the DMSR vector bindx                              m.gee 1/02 |
  | for format see Aztec manual                                        |
 *----------------------------------------------------------------------*/
void spo_make_bindx(
    FIELD         *actfield,
    PARTITION     *actpart,
    SOLVAR        *actsolv,
    SPOOLMAT      *spo,
    INT          **dof_connect,
    INT           *bindx
    )

{
  INT        i,j;
  INT        count1,count2;
  INT        dof;

#ifdef DEBUG
  dstrc_enter("spo_make_bindx");
#endif

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
  | calculate number of nonzero entries and dof topology    m.gee 1/02 |
 *----------------------------------------------------------------------*/
void  spo_nnz_topology(
    FIELD        *actfield,
    PARTITION    *actpart,
    SOLVAR       *actsolv,
    INTRA        *actintra,
    SPOOLMAT     *spo,
    INT         **dof_connect
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
  | allocate update put dofs in update in ascending order   m.gee 5/01 |
 *----------------------------------------------------------------------*/
void spo_update(
    FIELD         *actfield,
    PARTITION     *actpart,
    SOLVAR        *actsolv,
    INTRA         *actintra,
    SPOOLMAT      *spo
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
  dstrc_enter("spo_update");
#endif

  /*----------------------------------------------------------------------*/
  imyrank = actintra->intra_rank;
  inprocs = actintra->intra_nprocs;
  /*------------------ make a local copy of the array actpart->coupledofs */
  memset(&coupledofs, 0, sizeof(ARRAY));
  if (actpart->pdis[disnum].coupledofs.Typ != cca_XX)
    am_alloc_copy(&(actpart->pdis[disnum].coupledofs),&coupledofs);
  /*----------------------------------------------------------------------*/
  update = spo->update.a.iv;
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
  if (counter != spo->numeq) dserror("Number of dofs in spooles-vector update wrong");
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
} /* end of spo_update */





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
  cdofs      = actpart->pdis[disnum].coupledofs.a.ia;
  ncdofs     = actpart->pdis[disnum].coupledofs.fdim;

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
      add_spo_checkcouple(ii,cdofs,ncdofs,&ii_iscouple,&ii_owner,nprocs);
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

#endif /* ifdef FAST_ASS */






/*----------------------------------------------------------------------*/
/*!
  \brief assemble into a msr matrix (original version)

  This routine assembles one or two element matrices (estiff_global and
  emass_global) into the global matrices in the spooles format.
  It makes extensive use of the searchs provided by the function
  'find_index'.

  \param actpart   *PARTITION    (i)  the partition we are working on
  \param actsolv   *SOLVAR       (i)  the solver we are using
  \param actintra  *INTRA        (i)  the intra-communicator we do not need
  \param actele    *ELEMENT      (i)  the element we would like to work with
  \param spo1      *SPOOLMAT     (i)  one sparse matrix we will assemble into
  \param spo2      *SPOOLMAT     (i)  the other sparse matrix we will assemble into

  \author mn
  \date 07/04

*/
/*----------------------------------------------------------------------*/
void  add_spo(
    struct _PARTITION     *actpart,
    struct _SOLVAR        *actsolv,
    struct _INTRA         *actintra,
    struct _ELEMENT       *actele,
    struct _SPOOLMAT      *spo1,
    struct _SPOOLMAT      *spo2,
    struct _ARRAY         *elearray1,
    struct _ARRAY         *elearray2
    )

{

  INT         i,j,k,l,counter;    /* some counter variables */
  INT         istwo=0;
  INT         start,index,lenght; /* some more special-purpose counters */
  INT         ii,jj;              /* counter variables for system matrix */
  INT         ii_iscouple;        /* flag whether ii is a coupled dof */
  INT         ii_owner;           /* who is owner of dof ii -> procnumber */
  INT         ii_index;           /* place of ii in dmsr format */
  INT         jj_index;           /* place of jj in dmsr format */
  INT         nd,ndnd;            /* size of estif */
  INT         nnz;                /* number of nonzeros in sparse system matrix */
  INT         numeq_total;        /* total number of equations */
  INT         numeq;              /* number of equations on this proc */
  INT         lm[MAXDOFPERELE];   /* location vector for this element */
  INT         owner[MAXDOFPERELE];/* the owner of every dof */
  INT         myrank;             /* my intra-proc number */
  INT         nprocs;             /* my intra- number of processes */
  DOUBLE    **estif;              /* element matrix to be added to system matrix */
  DOUBLE    **emass;              /* element matrix to be added to system matrix */
  INT        *update;             /* vector update see AZTEC manual */
  DOUBLE     *A_loc;              /*    "       A_loc see MUMPS manual */
  DOUBLE     *B_loc;              /*    "       A_loc see MUMPS manual */
  INT        *irn;                /*    "       irn see MUMPS manual */
  INT        *jcn;                /*    "       jcn see MUMPS manual */
  INT        *rowptr;             /*    "       rowptr see rc_ptr structure */
  INT       **cdofs;              /* list of coupled dofs and there owners, see init_assembly */
  INT         ncdofs;             /* total number of coupled dofs */
  INT       **isend1;             /* pointer to sendbuffer to communicate coupling conditions */
  DOUBLE    **dsend1;             /* pointer to sendbuffer to communicate coupling conditions */
  INT       **isend2;             /* pointer to sendbuffer to communicate coupling conditions */
  DOUBLE    **dsend2;             /* pointer to sendbuffer to communicate coupling conditions */
  INT         nsend;


#ifdef DEBUG
  dstrc_enter("add_spo");
#endif

  /* check whether to assemble one or two matrices */
  if (spo2) istwo=1;

  /* set some pointers and variables */
  myrank     = actintra->intra_rank;
  nprocs     = actintra->intra_nprocs;
  estif      = elearray1->a.da;
  if (istwo)
    emass      = elearray2->a.da;
  nd         = actele->numnp * actele->node[0]->numdf;
  ndnd       = nd*nd;
  nnz        = spo1->nnz;
  numeq_total= spo1->numeq_total;
  numeq      = spo1->numeq;
  update     = spo1->update.a.iv;
  A_loc      = spo1->A_loc.a.dv;
  if (istwo)
    B_loc      = spo2->A_loc.a.dv;
  irn        = spo1->irn_loc.a.iv;
  jcn        = spo1->jcn_loc.a.iv;
  rowptr     = spo1->rowptr.a.iv;
  cdofs      = actpart->pdis[disnum].coupledofs.a.ia;
  ncdofs     = actpart->pdis[disnum].coupledofs.fdim;

  /* put pointers to sendbuffers if any */
#ifdef PARALLEL
  if (spo1->couple_i_send)
  {
    isend1 = spo1->couple_i_send->a.ia;
    dsend1 = spo1->couple_d_send->a.da;
    nsend  = spo1->couple_i_send->fdim;
    if (istwo)
    {
      isend2 = spo2->couple_i_send->a.ia;
      dsend2 = spo2->couple_d_send->a.da;
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
  }/* end of loop over element nodes */
  /* end of loop over element nodes */

  /* this check is not possible any more for fluid element with implicit
     free surface condition: nd not eqaual numnp*numdf!!!                    */
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
      add_spo_checkcouple(ii,cdofs,ncdofs,&ii_iscouple,&ii_owner,nprocs);
    }
#endif

    /* ii is not a coupled dofs or I am master owner */
    ii_index      = find_index(ii,update,numeq);
#ifndef D_CONTACT
    if (!ii_iscouple || ii_owner==myrank)
    {

      if (ii_index==-1) dserror("dof ii not found on this proc");
      start         = rowptr[ii_index];
      lenght        = rowptr[ii_index+1]-rowptr[ii_index];

    }
#endif

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
#ifdef D_CONTACT
        add_val_spo(ii,ii_index,jj,spo1,estif[i][j],actintra);
        if (istwo)
          add_val_spo(ii,ii_index,jj,spo2,emass[i][j],actintra);
#else
        index         = find_index(jj,&(jcn[start]),lenght);
        if (index==-1) dserror("dof jj not found in this row ii");
        index        += start;
        A_loc[index] += estif[i][j];
        if (istwo)
          B_loc[index] += emass[i][j];
#endif
      }

      /* do main-diagonal entry */
      /* (a coupled dof and I am slave owner) */
      else
      {
        add_spo_sendbuff(ii,jj,i,j,ii_owner,isend1,dsend1,estif,nsend);
        if (istwo)
          add_spo_sendbuff(ii,jj,i,j,ii_owner,isend2,dsend2,emass,nsend);
      }

    } /* end loop over j */
  }/* end loop over i */


#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of add_spo */






#ifdef FAST_ASS

/*----------------------------------------------------------------------*/
/*!
  \brief assemble into a spooles matrix (faster version)

  This routine assembles one or two element matrices (estiff_global and
  emass_global) into the global matrices in the spooles format.
  It makes use of the information saved in actele->index to determine the
  correct position in the sparse matrix.
  This is faster then searching every time, but consumes a lot of memory!!

  \param actpart   *PARTITION    (i)  the partition we are working on
  \param actsolv   *SOLVAR       (i)  the solver we are using
  \param actintra  *INTRA        (i)  the intra-communicator we do not need
  \param actele    *ELEMENT      (i)  the element we would like to work with
  \param spo1      *SPOOLMAT     (i)  one sparse matrix we will assemble into
  \param spo2      *SPOOLMAT     (i)  the other sparse matrix we will assemble into

  \author mn
  \date 07/04

*/
/*----------------------------------------------------------------------*/
void  add_spo_fast(
    struct _PARTITION     *actpart,
    struct _SOLVAR        *actsolv,
    struct _INTRA         *actintra,
    struct _ELEMENT       *actele,
    struct _SPOOLMAT      *spo1,
    struct _SPOOLMAT      *spo2,
    struct _ARRAY         *elearray1,
    struct _ARRAY         *elearray2
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
  DOUBLE    **estif;                    /* element matrix to be added to system matrix */
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

#ifdef DEBUG
  dstrc_enter("add_spo_fast");
#endif

  /* check whether to assemble one or two matrices */
  if (spo2) istwo=1;

  /* set some pointers and variables */
  myrank     = actintra->intra_rank;
  nprocs     = actintra->intra_nprocs;
  estif      = elearray1->a.da;
  if (istwo)
    emass      = elearray2->a.da;
  nd         = actele->nd;
  ndnd       = nd*nd;
  nnz        = spo1->nnz;
  numeq_total= spo1->numeq_total;
  numeq      = spo1->numeq;
  update     = spo1->update.a.iv;
  A_loc      = spo1->A_loc.a.dv;
  if (istwo)
    B_loc      = spo2->A_loc.a.dv;
  irn        = spo1->irn_loc.a.iv;
  jcn        = spo1->jcn_loc.a.iv;
  rowptr     = spo1->rowptr.a.iv;
  cdofs      = actpart->pdis[disnum].coupledofs.a.ia;
  ncdofs     = actpart->pdis[disnum].coupledofs.fdim;

  /* put pointers to sendbuffers if any */
#ifdef PARALLEL
  if (spo1->couple_i_send)
  {
    isend1 = spo1->couple_i_send->a.ia;
    dsend1 = spo1->couple_d_send->a.da;
    nsend  = spo1->couple_i_send->fdim;
    if (istwo)
    {
      isend2 = spo2->couple_i_send->a.ia;
      dsend2 = spo2->couple_d_send->a.da;
    }
  }
#endif


  /* loop over i (the element row) */
  ii_iscouple = 0;
  ii_owner    = myrank;

  for (i=0; i<nd; i++)
  {
    ii = actele->locm[i];
#ifdef D_CONTACT
    ii_index = find_index(ii,update,numeq);
#endif

    /* loop over j (the element column) */
    for (j=0; j<nd; j++)
    {
      jj = actele->locm[j];
      index = actele->index[i][j];

      if(index >= 0)  /* normal dof */
      {
#ifdef D_CONTACT
        ii_index      = find_index(ii,update,numeq);
        add_val_spo(ii,ii_index,jj,spo1,estif[i][j],actintra);
        if (istwo)
          add_val_spo(ii,ii_index,jj,spo2,emass[i][j],actintra);
#else
        A_loc[index] += estif[i][j];
        if (istwo)
          B_loc[index] += emass[i][j];
#endif
      }


      if(index == -1)  /* boundary condition dof */
        continue;


      if(index == -2)  /* coupled dof */
      {
        add_spo_sendbuff(ii,jj,i,j,ii_owner,isend1,dsend1,estif,nsend);
        if (istwo)
          add_spo_sendbuff(ii,jj,i,j,ii_owner,isend2,dsend2,emass,nsend);
      }

    } /* end loop over j */
  }/* end loop over i */


#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of add_spo_fast */

#endif  /* ifdef FAST_ASS */





/*----------------------------------------------------------------------*
  | add value to spooles matrix (with enlargment)          m.gee 11/02 |
 *----------------------------------------------------------------------*/
void add_val_spo(
    INT                ii,
    INT                index,
    INT                jj,
    struct _SPOOLMAT  *spo,
    DOUBLE             val,
    INTRA             *actintra
    )

{
  INT     i,j,k,l,colstart,colend,foundit;
  INT     counter,hasmoved;
  INT    *irn,*jcn,*update,*rptr,numeq;
  DOUBLE *A;
  INT     nnz,nnz_new;
  INT     rsize;
  INT     move_index;
  INT     sf_col,ef_col,st_col,et_col;

#ifdef DEBUG
  dstrc_enter("add_val_spo");
#endif

  /*----------------------------------------------------------------------*/
  irn    = spo->irn_loc.a.iv;
  jcn    = spo->jcn_loc.a.iv;
  update = spo->update.a.iv;
  numeq  = spo->numeq;
  rptr   = spo->rowptr.a.iv;
  A      = spo->A_loc.a.dv;
  /*----------------------------------------------------------------------*/
  /*index  = find_index(ii,update,numeq);*/
  if (index==-1) dserror("Cannot find local dof");
  colstart = rptr[index];
  colend   = rptr[index+1];
  /* check whether entry jj already exists */
  foundit=find_index(jj,&jcn[colstart],colend-colstart);
  /* found the entry jj */
  if (foundit != -1)
  {
    A[colstart+foundit] += val;
  }
  /* entry does not exists in this row */
  else
  {
    if (jcn[colstart]==-1) /* there is still room in this row */
    {
      irn[colstart] = ii;
      jcn[colstart] = jj;
      A[colstart] = val;
      mg_sort(&(jcn[colstart]),colend-colstart,&(irn[colstart]),&(A[colstart]));
    }
    else /* there is no room in matrix */
    {
startenlarge:
      nnz = spo->A_loc.fdim;
      nnz_new = (INT)(1.5*nnz);
      rsize   = (INT)(nnz_new/numeq);
      nnz_new = rsize*numeq;
      irn = amredef(&(spo->irn_loc),nnz_new,1,"IV");
      jcn = amredef(&(spo->jcn_loc),nnz_new,1,"IV");
      A   = amredef(&(spo->A_loc)  ,nnz_new,1,"DV");
      /* make sure, the new rowsize is larger then the old one */
      for (l=0; l<numeq; l++)
        if (rptr[l+1]-rptr[l] > rsize-1)
          goto startenlarge;
      /* init the new part of irn and jcn */
      for (i=nnz; i<nnz_new; i++)
      {
        irn[i]=-1;
        jcn[i]=-1;
      }
      /* loop row from back to front and move them */
      for (k=numeq-1; k>=0; k--)
      {
        move_index = k;
        /* get old column range */
        sf_col = rptr[move_index];
        ef_col = rptr[move_index+1];
        /* get new column range */
        st_col = move_index * rsize;
        et_col = (move_index+1) *rsize;
        /* loop the old row and copy to new location */
        counter=0;
        hasmoved=0;
        for (l=sf_col; l<ef_col; l++)
        {
          if (jcn[l]==-1) continue;
          hasmoved++;
          jcn[st_col+counter] = jcn[l];
          irn[st_col+counter] = irn[l];
          A[st_col+counter]   = A[l];
          counter++;
        }
        /* make sure, there is no old stuff in the new row */
        for (l=st_col+counter; l<et_col; l++)
        {
          jcn[l] = -1;
          irn[l] = -1;
          A[l]   = 0.0;
        }
        /* transfer complete, now sort the values to the back of new row */
        mg_sort(&(jcn[st_col]),et_col-st_col,&(irn[st_col]),&(A[st_col]));
        /* set new ptr in rptr */
        rptr[move_index+1] = et_col;
      }
      /* now there is room in the row, add value */
      /* the row has moved, so get the range again */
      colstart = rptr[index];
      colend   = rptr[index+1];
      if (jcn[colstart]==-1)
      {
        jcn[colstart] = jj;
        irn[colstart] = ii;
        A[colstart]   = val;
        mg_sort(&(jcn[colstart]),colend-colstart,&(irn[colstart]),&(A[colstart]));
      }
      else
        dserror("Fatal error in redimensioning of system matrix");
    }
  }
  /*----------------------------------------------------------------------*/

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of add_val_spo */






/*----------------------------------------------------------------------*
  | set value to spooles matrix (with enlargment)          m.gee 11/02 |
 *----------------------------------------------------------------------*/
void set_val_spo(
    INT                ii,
    INT                index,
    INT                jj,
    struct _SPOOLMAT  *spo,
    DOUBLE             val,
    INTRA             *actintra
    )

{
  INT     i,j,k,l,colstart,colend,foundit;
  INT     counter,hasmoved;
  INT    *irn,*jcn,*update,*rptr,numeq;
  DOUBLE *A;
  INT     nnz,nnz_new;
  INT     rsize;
  INT     move_index;
  INT     sf_col,ef_col,st_col,et_col;

#ifdef DEBUG
  dstrc_enter("set_val_spo");
#endif

  /*----------------------------------------------------------------------*/
  irn    = spo->irn_loc.a.iv;
  jcn    = spo->jcn_loc.a.iv;
  update = spo->update.a.iv;
  numeq  = spo->numeq;
  rptr   = spo->rowptr.a.iv;
  A      = spo->A_loc.a.dv;
  /*----------------------------------------------------------------------*/
  /*index  = find_index(ii,update,numeq);*/
  if (index==-1) dserror("Cannot find local dof");
  colstart = rptr[index];
  colend   = rptr[index+1];
  /* check whether entry jj already exists */
  foundit=find_index(jj,&jcn[colstart],colend-colstart);
  /* found the entry jj */
  if (foundit != -1)
  {
    A[colstart+foundit] = val;
  }
  /* entry does not exists in this row */
  else
  {
    if (jcn[colstart]==-1) /* there is still room in this row */
    {
      irn[colstart] = ii;
      jcn[colstart] = jj;
      A[colstart] = val;
      mg_sort(&(jcn[colstart]),colend-colstart,&(irn[colstart]),&(A[colstart]));
    }
    else /* there is no room in matrix */
    {
startenlarge:
      printf("MESSAGE: Enlargment of system matrix\n");
      nnz = spo->A_loc.fdim;
      nnz_new = (INT)(1.5*nnz);
      rsize   = (INT)(nnz_new/numeq);
      nnz_new = rsize*numeq;
      irn = amredef(&(spo->irn_loc),nnz_new,1,"IV");
      jcn = amredef(&(spo->jcn_loc),nnz_new,1,"IV");
      A   = amredef(&(spo->A_loc)  ,nnz_new,1,"DV");
      /* make sure, the new rowsize is larger then the old one */
      for (l=0; l<numeq; l++)
        if (rptr[l+1]-rptr[l] > rsize-1)
          goto startenlarge;
      /* init the new part of irn and jcn */
      for (i=nnz; i<nnz_new; i++)
      {
        irn[i]=-1;
        jcn[i]=-1;
      }
      /* loop row from back to front and move them */
      for (k=numeq-1; k>=0; k--)
      {
        move_index = k;
        /* get old column range */
        sf_col = rptr[move_index];
        ef_col = rptr[move_index+1];
        /* get new column range */
        st_col = move_index * rsize;
        et_col = (move_index+1) *rsize;
        /* loop the old row and copy to new location */
        counter=0;
        hasmoved=0;
        for (l=sf_col; l<ef_col; l++)
        {
          if (jcn[l]==-1) continue;
          hasmoved++;
          jcn[st_col+counter] = jcn[l];
          irn[st_col+counter] = irn[l];
          A[st_col+counter]   = A[l];
          counter++;
        }
        /* make sure, there is no old stuff in the new row */
        for (l=st_col+counter; l<et_col; l++)
        {
          jcn[l] = -1;
          irn[l] = -1;
          A[l]   = 0.0;
        }
        /* transfer complete, now sort the values to the back of new row */
        mg_sort(&(jcn[st_col]),et_col-st_col,&(irn[st_col]),&(A[st_col]));
        /* set new ptr in rptr */
        rptr[move_index+1] = et_col;
      }
      /* now there is room in the row, add value */
      /* the row has moved, so get the range again */
      colstart = rptr[index];
      colend   = rptr[index+1];
      if (jcn[colstart]==-1)
      {
        jcn[colstart] = jj;
        irn[colstart] = ii;
        A[colstart]   = val;
        mg_sort(&(jcn[colstart]),colend-colstart,&(irn[colstart]),&(A[colstart]));
      }
      else
        dserror("Fatal error in redimensioning of system matrix");
    }
  }
  /*----------------------------------------------------------------------*/

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of set_val_spo */





/*----------------------------------------------------------------------*
  | finalize the assembly to the spooles matrix            m.gee 11/02 |
 *----------------------------------------------------------------------*/
void close_spooles_matrix(
    struct _SPOOLMAT   *spo,
    INTRA              *actintra
    )

{
  INT     i,j,k,index,colstart,colend,foundit,actrow,offset;
  INT    *irn,*jcn,*update,*rptr,numeq;
  DOUBLE *A;

#ifdef DEBUG
  dstrc_enter("close_spooles_matrix");
#endif

  /*----------------------------------------------------------------------*/
  irn    = spo->irn_loc.a.iv;
  jcn    = spo->jcn_loc.a.iv;
  update = spo->update.a.iv;
  numeq  = spo->numeq;
  rptr   = spo->rowptr.a.iv;
  A      = spo->A_loc.a.dv;
  /*----------------------------------------------------------------------*/
  /*-------------------------- loop the rows and move values to the front */
  for (i=0; i<numeq; i++)
  {
    actrow = i;
    colstart = rptr[actrow];
    colend   = rptr[actrow+1];
    offset   = 0;
    /* count number of unused entries in this row */
    for (j=colstart; j<colend; j++)
    {
      if (jcn[j]==-1) offset++;
      else            break;
    }
    /* do nothing, if there is no offset */
    if (offset==0) continue;
    /* move all values offset to the front */
    for (k=j; k<colend; k++)
    {
      jcn[k-offset] = jcn[k];
      irn[k-offset] = irn[k];
      A[k-offset]   = A[k];
      jcn[k]        = -1;
      irn[k]        = -1;
      A[k]          = 0.0;
    }
    /* resize the row */
    rptr[actrow+1] = colend-offset;
  }
  /* redefine the array */
  if (spo->jcn_loc.fdim > (INT)(1.2 * rptr[numeq]))
  {
    amredef(&(spo->jcn_loc),rptr[numeq],1,"IV");
    amredef(&(spo->irn_loc),rptr[numeq],1,"IV");
    amredef(&(spo->A_loc)  ,rptr[numeq],1,"DV");
  }
  /* set correct number of nonzeros */
  spo->nnz = rptr[numeq];
  /*----------------------------------------------------------------------*/

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of close_spooles_matrix */






/*----------------------------------------------------------------------*
  | add one spooles matrix to another                      m.gee 11/02 |
  | init=0 : to = to + from * factor                                   |
  | init=1:  to =      from * factor                                   |
 *----------------------------------------------------------------------*/
void add_spooles_matrix(
    struct _SPOOLMAT   *to,
    struct _SPOOLMAT   *from,
    DOUBLE              factor,
    INT                 init,
    INTRA              *actintra
    )

{
  INT     i,j,k,index,foundit,actrow,offset;
  INT     colstart_to,colend_to;
  INT     colstart_from,colend_from;
  INT    *irn_to,*jcn_to,*update_to,*rptr_to,numeq_to;
  DOUBLE *A_to;
  INT    *irn_from,*jcn_from,*update_from,*rptr_from,numeq_from;
  DOUBLE *A_from;

#ifdef DEBUG
  dstrc_enter("add_spooles_matrix");
#endif

  /*----------------------------------------------------------------------*/
  irn_to    = to->irn_loc.a.iv;
  jcn_to    = to->jcn_loc.a.iv;
  update_to = to->update.a.iv;
  numeq_to  = to->numeq;
  rptr_to   = to->rowptr.a.iv;
  A_to      = to->A_loc.a.dv;
  irn_from    = from->irn_loc.a.iv;
  jcn_from    = from->jcn_loc.a.iv;
  update_from = from->update.a.iv;
  numeq_from  = from->numeq;
  rptr_from   = from->rowptr.a.iv;
  A_from      = from->A_loc.a.dv;
  /*----------------------------------------------------------------------*/
  if (numeq_to != numeq_from) dserror("Number of rows not equal");
  /*-------------------------- loop the rows and move values to the front */
  for (i=0; i<numeq_to; i++)
  {
    actrow        = i;
    colstart_to   = rptr_to[actrow];
    colend_to     = rptr_to[actrow+1];
    colstart_from = rptr_from[actrow];
    colend_from   = rptr_from[actrow+1];
    for (j=colstart_from; j<colend_from; j++)
    {
      index = find_index(irn_from[j],update_to,numeq_to);
      if (!init)
        add_val_spo(irn_from[j],index,jcn_from[j],to,A_from[j]*factor,actintra);
      else
        set_val_spo(irn_from[j],index,jcn_from[j],to,A_from[j]*factor,actintra);
    }
  }
  close_spooles_matrix(to,actintra);
  /*----------------------------------------------------------------------*/

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of add_spooles_matrix */




/*----------------------------------------------------------------------*
  | checks coupling for the add_spo routine                 m.gee 9/01 |
 *----------------------------------------------------------------------*/
void add_spo_checkcouple(
    INT         ii,
    INT       **cdofs,
    INT         ncdofs,
    INT        *iscouple,
    INT        *isowner,
    INT         nprocs
    )

{

  INT         i,k;

#ifdef DEBUG
  dstrc_enter("add_spo_checkcouple");
#endif

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

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of add_spo_checkcouple */





/*----------------------------------------------------------------------*
  | fill sendbuffer isend and dsend                         m.gee 1/02 |
 *----------------------------------------------------------------------*/
void add_spo_sendbuff(
    INT          ii,
    INT          jj,
    INT          i,
    INT          j,
    INT          ii_owner,
    INT        **isend,
    DOUBLE     **dsend,
    DOUBLE     **estif,
    INT          numsend
    )

{
  INT         k,l;

#ifdef DEBUG
  dstrc_enter("add_spo_sendbuff");
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
} /* end of add_spo_sendbuff */





/*----------------------------------------------------------------------*
  | exchange coupled dofs and add to row/column ptr matrix  m.gee 1/02 |
 *----------------------------------------------------------------------*/
void exchange_coup_spo(
    PARTITION     *actpart,
    SOLVAR        *actsolv,
    INTRA         *actintra,
    SPOOLMAT      *spo
    )

{
  INT            i,j;
  INT            ii,jj,ii_index;
  INT            start;
  INT            lenght;
  INT            index;
  INT            tag;
  INT            source;
  INT            numeq,numeq_total;
  INT            numsend;
  INT            numrecv;
  INT           *update;
  INT          **isend;
  DOUBLE       **dsend;
  INT          **irecv;
  DOUBLE       **drecv;
  INT            imyrank;
  INT            inprocs;
  DOUBLE        *A_loc;                    /*    "       A_loc see MUMPS manual */
  INT           *irn;                      /*    "       irn see MUMPS manual */
  INT           *jcn;                      /*    "       jcn see MUMPS manual */
  INT           *rowptr;                   /*    "       rowptr see rc_ptr structure */

#ifdef PARALLEL
  MPI_Status    *irecv_status;
  MPI_Status    *drecv_status;

  MPI_Request   *isendrequest;
  MPI_Request   *dsendrequest;

  MPI_Comm      *ACTCOMM;
#endif

#ifdef DEBUG
  dstrc_enter("exchange_coup_spo");
#endif

  /*----------------------------------------------------------------------*/
#ifdef PARALLEL
  /*----------------------------------------------------------------------*/
  imyrank = actintra->intra_rank;
  inprocs = actintra->intra_nprocs;
  ACTCOMM = &(actintra->MPI_INTRA_COMM);
  /*---------------------------------------- set some pointers and values */
  numsend     = spo->numcoupsend;
  numrecv     = spo->numcouprecv;
  A_loc      = spo->A_loc.a.dv;
  irn        = spo->irn_loc.a.iv;
  jcn        = spo->jcn_loc.a.iv;
  rowptr     = spo->rowptr.a.iv;
  update      = spo->update.a.iv;
  numeq_total = spo->numeq_total;
  numeq       = spo->numeq;
  if (spo->couple_i_send) isend = spo->couple_i_send->a.ia;
  if (spo->couple_d_send) dsend = spo->couple_d_send->a.da;
  if (spo->couple_i_recv) irecv = spo->couple_i_recv->a.ia;
  if (spo->couple_d_recv) drecv = spo->couple_d_recv->a.da;
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
    ii_index = find_index(ii,update,numeq);
    if (ii_index==-1) dserror("dof ii not found on this proc");
    start         = rowptr[ii_index];
    lenght        = rowptr[ii_index+1]-rowptr[ii_index];
    for (j=0; j<lenght; j++)
    {
      index         = start+j;
      jj            = jcn[index];
      if (jj==-1) continue;
#ifdef D_CONTACT
      add_val_spo(ii,ii_index,jj,spo,drecv[i][jj],actintra);
#else
      A_loc[index] += drecv[i][jj];
#endif
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
} /* end of exchange_coup_spo */



#endif /* ifdef SPOOLES_PACKAGE */



