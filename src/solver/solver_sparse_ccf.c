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


#ifdef UMFPACK


#include "../headers/standardtypes.h"
#include "../solver/solver.h"


INT cmp_int(const void *a, const void *b );
DOUBLE cmp_double(const void *a, const void *b );



static INT  disnum;




/*----------------------------------------------------------------------*
  |  calculate the mask of an ccf matrix           s.offermanns 05/02    |
 *----------------------------------------------------------------------*/
void mask_ccf(
    FIELD         *actfield,
    PARTITION     *actpart,
    SOLVAR        *actsolv,
    INTRA         *actintra,
    CCF           *ccf,
    INT            disnum_
    )

{
  INT       i;
  INT       numeq;
  INT     **dof_connect;
  ARRAY     bindx_a;
  INT      *bindx;
  ARRAY     red_dof_connect;

#ifdef FAST_ASS
  ELEMENT  *actele;
#endif

#ifdef DEBUG
  dstrc_enter("mask_ccf");
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
  /*------------------------------------------- put size of problem */
  ccf->numeq_total = actfield->dis[disnum].numeq;
  /* count number of eqns on proc and build processor-global couplingdof
     matrix */
  mask_numeq(actfield,actpart,actsolv,actintra,&numeq,disnum);
  ccf->numeq = numeq;
  /*---------------------------------------------- allocate vector update */
  amdef("update",&(ccf->update),numeq,1,"IV");
  amzero(&(ccf->update));
  /*--------------------------------put dofs in update in ascending order */
  ccf_update(actfield,actpart,actsolv,actintra,ccf);
  /*------------------------ count number of nonzero entries on partition
    and calculate dof connectivity list */
  /*
     dof_connect[i][0] = lenght of dof_connect[i]
     dof_connect[i][1] = iscoupled ( 1 or 2 )
     dof_connect[i][2] = dof
     dof_connect[i][ 2..dof_connect[i][0]-1 ] = connected dofs exluding itself
     */
  dof_connect = (INT**)CCACALLOC(ccf->numeq_total,sizeof(INT*));
  if (!dof_connect) dserror("Allocation of dof_connect failed");
  /*---------------------- make the dof_connect list locally on each proc */
  ccf_nnz_topology(actfield,actpart,actsolv,actintra,ccf,dof_connect);
  /*------------------------------------------------------ make nnz_total */
#ifdef PARALLEL
  ccf->nnz_total=0;
  MPI_Allreduce(&(ccf->nnz),&(ccf->nnz_total),1,MPI_INT,MPI_SUM,actintra->MPI_INTRA_COMM);
#else
  ccf->nnz_total=ccf->nnz;
#endif
  /*------------------------------------------ make dof_connect redundant */
  ccf_red_dof_connect(actfield,actpart,actsolv,actintra,ccf,dof_connect,&red_dof_connect);
  /*----------------------------------------------------- allocate arrays */
  /*                                                   see UMFPACK manual */
  amdef("Ap_loc",&(ccf->Ap),ccf->numeq_total+1 ,1,"IV");
  amdef("Ai_loc",&(ccf->Ai),ccf->nnz_total     ,1,"IV");
  amdef("Ax",&(ccf->Ax),ccf->nnz_total     ,1,"DV");
  amzero(&(ccf->Ax));
  /*------------------------------------------------------ allocate bindx */
  bindx = amdef("bindx",&(bindx_a),(ccf->nnz_total+1),1,"IV");
  /*---------------------------------------------------------- make bindx */
  ccf_make_bindx(actfield,actpart,actsolv,ccf,bindx,&red_dof_connect);
  /*----------------- make rowptr, irn_loc, jcn_loc from bindx and update */
  ccf_make_sparsity(ccf,bindx);
  /*---------------------------------------- delete the array dof_connect */
  for (i=0; i<ccf->numeq_total; i++)
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
    ccf_make_index(actfield,actpart,actintra,actele,ccf);
  }
#endif

  /*----------------------------------------------------------------------*/

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of mask_ccf */





/*----------------------------------------------------------------------*
  |  allocate update put dofs in update in ascending order   m.gee 1/02  |
 *----------------------------------------------------------------------*/
void  ccf_update(
    FIELD         *actfield,
    PARTITION   *actpart,
    SOLVAR      *actsolv,
    INTRA       *actintra,
    CCF         *ccf
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
  dstrc_enter("ccf_update");
#endif

  /*----------------------------------------------------------------------*/
  imyrank = actintra->intra_rank;
  inprocs = actintra->intra_nprocs;
  /*------------------ make a local copy of the array actpart->coupledofs */
  am_alloc_copy(&(actpart->pdis[disnum].coupledofs),&coupledofs);
  /*------------------------------------- loop the nodes on the partition */
  update = ccf->update.a.iv;
  counter=0;
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
  /*---------- check whether the correct number of dofs have been counted */
  if (counter != ccf->numeq) dserror("Number of dofs in update wrong");
  /*---------------------------- sort the vector update just to make sure */
  qsort((INT*) update, counter, sizeof(INT), cmp_int);
  /*----------------------------------------------------------------------*/
  amdel(&coupledofs);
  /*----------------------------------------------------------------------*/

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of ccf_update */





/*----------------------------------------------------------------------*
  |  calculate number of nonzero entries and dof topology    m.gee 1/02  |
 *----------------------------------------------------------------------*/
void  ccf_nnz_topology(
    FIELD         *actfield,
    PARTITION    *actpart,
    SOLVAR       *actsolv,
    INTRA        *actintra,
    CCF          *ccf,
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

#ifdef PARALLEL
  MPI_Status status;
#endif

#ifdef DEBUG
  dstrc_enter("ccf_nnz_topology");
#endif

  /*----------------------------------------------------------------------*/
  imyrank = actintra->intra_rank;
  inprocs = actintra->intra_nprocs;
  /*----------------------------------------------------------- shortcuts */
  ccf->nnz=0;
  numeq  = ccf->numeq;
  update = ccf->update.a.iv;
  for (i=0; i<ccf->numeq_total; i++) dof_connect[i]=NULL;
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
  for (i=0; i<ccf->update.fdim; i++)
  {
    dof = ccf->update.a.iv[i];
    nnz += (dof_connect[dof][0]-2);
  }
  ccf->nnz=nnz;
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
} /* end of ccf_nnz_topology */





/*----------------------------------------------------------------------*
  |  make redundant dof_connect array                        m.gee 1/02  |
 *----------------------------------------------------------------------*/
void ccf_red_dof_connect(
    FIELD        *actfield,
    PARTITION    *actpart,
    SOLVAR       *actsolv,
    INTRA        *actintra,
    CCF          *ccf,
    INT         **dof_connect,
    ARRAY        *red_dof_connect
    )

{
  INT        i,j;
  INT        max_dof_connect_send;
  INT        max_dof_connect_recv;

#ifdef PARALLEL
  ARRAY      tmps_a;
  INT      **tmps;
#endif
  INT      **reddof;

#ifdef DEBUG
  dstrc_enter("ccf_red_dof_connect");
#endif

  /*----------------------------- check for largest row in my dof_connect */
  max_dof_connect_send=0;
  max_dof_connect_recv=0;
  for (i=0; i<ccf->numeq_total; i++)
  {
    if (dof_connect[i])
      if (dof_connect[i][0]>max_dof_connect_send)
        max_dof_connect_send=dof_connect[i][0];
  }
#ifdef PARALLEL
  MPI_Allreduce(&max_dof_connect_send,&max_dof_connect_recv,1,MPI_INT,MPI_MAX,actintra->MPI_INTRA_COMM);
#else
  max_dof_connect_recv=max_dof_connect_send;
#endif
  /*---------------- allocate temporary array to hold global connectivity */
  reddof = amdef("tmp",red_dof_connect,ccf->numeq_total,max_dof_connect_recv,"IA");
  /*---------------- allocate temporary array to hold global connectivity */
#ifdef PARALLEL
  tmps = amdef("tmp",&tmps_a,ccf->numeq_total,max_dof_connect_recv,"IA");
  amzero(&tmps_a);
  /*-------------------------------- put my own dof_connect values to tmp */
  for (i=0; i<ccf->numeq_total; i++)
  {
    if (dof_connect[i])
    {
      for (j=0; j<dof_connect[i][0]; j++)
        tmps[i][j] = dof_connect[i][j];
    }
  }
#else
  /*----------------------------- put my own dof_connect values to reddof */
  for (i=0; i<ccf->numeq_total; i++)
  {
    if (dof_connect[i])
    {
      for (j=0; j<dof_connect[i][0]; j++)
        reddof[i][j] = dof_connect[i][j];
    }
  }
#endif
  /*--------------------------------------------- allreduce the array tmp */
#ifdef PARALLEL
  MPI_Allreduce(tmps[0],reddof[0],(tmps_a.fdim*tmps_a.sdim),MPI_INT,MPI_SUM,actintra->MPI_INTRA_COMM);
#endif
  /*----------------------------------------------------------------------*/
#ifdef PARALLEL
  amdel(&tmps_a);
#endif
  /*----------------------------------------------------------------------*/

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of ccf_red_dof_connect */




/*----------------------------------------------------------------------*
  |  make the DMSR vector bindx                              m.gee 1/02  |
  | for format see Aztec manual                                          |
 *----------------------------------------------------------------------*/
void  ccf_make_bindx(
    FIELD         *actfield,
    PARTITION     *actpart,
    SOLVAR        *actsolv,
    CCF           *ccf,
    INT           *bindx,
    ARRAY         *red_dof_connect
    )

{
  INT        i,j;
  INT        count1,count2;
  INT      **reddof;

#ifdef DEBUG
  dstrc_enter("ccf_make_bindx");
#endif

  reddof = red_dof_connect->a.ia;
  /*----------------------------------------------------------------------*/
  /*-------------------------------------------------------------do bindx */
  count1=0;
  count2=ccf->numeq_total+1;
  for (i=0; i<ccf->numeq_total; i++)
  {
    if (reddof[i])
    {
      bindx[count1] = count2;
      for (j=3; j<reddof[i][0]; j++)
      {
        bindx[count2] = reddof[i][j];
        count2++;
      }
    }
    count1++;
  }
  bindx[ccf->numeq_total] = ccf->nnz_total+1;
  /*----------------------------------------------------------------------*/

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of ccf_make_bindx */





/*----------------------------------------------------------------------*
  |  make the vectors                                s.offermanns 05/02  |
  | irn_loc, jcn_loc, rowptr from update and bindx                       |
 *----------------------------------------------------------------------*/
void  ccf_make_sparsity(
    CCF     *ccf,
    INT     *bindx
    )

{
  INT        i,j;
  INT        start,end;
  INT        counter;
  INT        numeq;
  INT        numeq_total;
  INT        nnz;
  INT       *update;
  INT       *Ap;
  INT       *Ai;

#ifdef DEBUG
  dstrc_enter("ccf_make_sparsity");
#endif

  /*----------------------------------------------------------------------*/
  numeq       = ccf->numeq;
  numeq_total = ccf->numeq_total;
  nnz         = ccf->nnz;
  update      = ccf->update.a.iv;
  Ap          = ccf->Ap.a.iv;
  Ai          = ccf->Ai.a.iv;
  /*------------------------------------------ loop all dofs on this proc */
  counter=0;
  for (i=0; i<numeq_total; i++)
  {
    start    = bindx[i];
    end      = bindx[i+1];
    Ap[i]    = counter;
    j=start;
    while (j<end && bindx[j]<i)/* dofs lower then actdof */
    {
      Ai[counter]=bindx[j];
      counter++;
      j++;
    }
    /*------------------------------- main diagonal of actdof */
    Ai[counter]=i;
    counter++;
    while(j<end)/*------------------- dofs higher then actdof */
    {
      Ai[counter]=bindx[j];
      counter++;
      j++;
    }
  }
  Ap[i]=counter;

  if (counter != ccf->nnz_total) dserror("sparsity mask failure");
  /*----------------------------------------------------------------------*/

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of ccf_make_sparsity */





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
  \param ccf1      *CCF          (i)  the sparse matrix we will assemble into

  \author mn
  \date 07/04

*/
/*----------------------------------------------------------------------*/
void ccf_make_index(
    FIELD                 *actfield,
    PARTITION             *actpart,
    INTRA                 *actintra,
    ELEMENT               *actele,
    struct _CCF           *ccf1
    )
{

  INT         i,j,counter;              /* some counter variables */
  INT         start,index,lenght;       /* some more special-purpose counters */
  INT         ii,jj;                    /* counter variables for system matrix */
  INT         nd;                       /* size of estif */
  INT         nnz;                      /* number of nonzeros in sparse system matrix */
  INT         numeq_total;              /* total number of equations */
  INT         numeq;                    /* number of equations on this proc */
  INT         myrank;                   /* my intra-proc number */
  INT         nprocs;                   /* my intra- number of processes */
  INT        *Ai;                       /*    "       Ap (column pointer) see UMFPACK manual */
  INT        *Ap;                       /*    "       Ai see UMFPACK manual */
  DOUBLE     *Ax;                       /*    "       Ax see UMFPACK manual */
  INT        *update;                   /* vector update see AZTEC manual */
  INT       **cdofs;                    /* list of coupled dofs and there owners, see init_assembly */
  INT         ncdofs;                   /* total number of coupled dofs */

  struct _ARRAY ele_index;
  struct _ARRAY ele_locm;
#ifdef PARALLEL
  struct _ARRAY ele_owner;
#endif


#ifdef DEBUG
  dstrc_enter("ccf_make_index");
#endif

  /* set some pointers and variables */
  myrank     = actintra->intra_rank;
  nprocs     = actintra->intra_nprocs;
  nnz        = ccf1->nnz;
  numeq_total= ccf1->numeq_total;
  numeq      = ccf1->numeq;
  update     = ccf1->update.a.iv;
  Ax         = ccf1->Ax.a.dv;
  Ai         = ccf1->Ai.a.iv;
  Ap         = ccf1->Ap.a.iv;
  cdofs      = actpart->pdis[disnum].coupledofs.a.ia;
  ncdofs     = actpart->pdis[disnum].coupledofs.fdim;


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
  /* loop over i (the element column) */
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

    /* loop over j (the element row) */
    start         = Ap[ii];
    lenght        = Ap[ii+1]-Ap[ii];
    for (j=0; j<nd; j++)
    {
      jj = actele->locm[j];

      /* check for boundary condition */
      if (jj>=numeq_total)
      {
        actele->index[i][j] = -1;
        continue;
      }

      index         = find_index(jj,&(Ai[start]),lenght);
      if (index==-1) dserror("dof jj not found in this row ii");
      index        += start;
      actele->index[i][j] = index;

    } /* end loop over j */
  }/* end loop over i */


#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of ccf_make_index */

#endif /* ifdef FAST_ASS */


/*----------------------------------------------------------------------*/
/*!
  \brief assemble into a ccf matrix (original version)

  This routine assembles one or two element matrices (estiff_global and
  emass_global) into the global matrices in the ccf format.
  It makes extensive use of the searchs provided by the function
  'find_index'.

  \param actpart   *PARTITION    (i)  the partition we are working on
  \param actsolv   *SOLVAR       (i)  the solver we are using
  \param actintra  *INTRA        (i)  the intra-communicator we do not need
  \param actele    *ELEMENT      (i)  the element we would like to work with
  \param ccf1      *CCF          (i)  one sparse matrix we will assemble into
  \param ccf2      *CCF          (i)  the other sparse matrix we will assemble into

  \author mn
  \date 07/04

*/
/*----------------------------------------------------------------------*/
void  add_ccf(
    struct _PARTITION     *actpart,
    struct _SOLVAR        *actsolv,
    struct _INTRA         *actintra,
    struct _ELEMENT       *actele,
    struct _CCF           *ccf1,
    struct _CCF           *ccf2,
    struct _ARRAY         *elearray1,
    struct _ARRAY         *elearray2)
{

  INT         i,j,counter;        /* some counter variables */
  INT         istwo=0;
  INT         start,index,lenght; /* some more special-purpose counters */
  INT         ii,jj;              /* counter variables for system matrix */
  INT         nd,ndnd;            /* size of estif */
  INT         nnz;                /* number of nonzeros in sparse system matrix */
  INT         numeq_total;        /* total number of equations */
  INT         numeq;              /* number of equations on this proc */
  INT         lm[MAXDOFPERELE];   /* location vector for this element */
  INT         myrank;             /* my intra-proc number */
  INT         nprocs;             /* my intra- number of processes */
  INT        *Ai;                 /*    "       Ap (column pointer) see UMFPACK manual */
  INT        *Ap;                 /*    "       Ai see UMFPACK manual */
  DOUBLE     *Ax;                 /*    "       Ax see UMFPACK manual */
  DOUBLE     *Bx;                 /*    "       Ax see UMFPACK manual */
  DOUBLE    **estif;              /* element matrix to be added to system matrix */
  DOUBLE    **emass;              /* element matrix to be added to system matrix */
  INT        *update;             /* vector update see AZTEC manual */
  INT       **cdofs;              /* list of coupled dofs and there owners, see init_assembly */
  INT         ncdofs;             /* total number of coupled dofs */

#ifdef PARALLEL
  INT         owner[MAXDOFPERELE];/* the owner of every dof */
  INT       **isend1;             /* pointer to sendbuffer to communicate coupling conditions */
  DOUBLE    **dsend1;             /* pointer to sendbuffer to communicate coupling conditions */
  INT       **isend2;             /* pointer to sendbuffer to communicate coupling conditions */
  DOUBLE    **dsend2;             /* pointer to sendbuffer to communicate coupling conditions */
  INT         nsend;
#endif

#ifdef DEBUG
  dstrc_enter("add_ccf");
#endif

  /* check whether to assemble one or two matrices */
  if (ccf2) istwo=1;

  /* set some pointers and variables */
  myrank     = actintra->intra_rank;
  nprocs     = actintra->intra_nprocs;
  estif      = elearray1->a.da;
  if (istwo)
    emass      = elearray2->a.da;
  nd         = actele->numnp * actele->node[0]->numdf;
  ndnd       = nd*nd;
  nnz        = ccf1->nnz;
  numeq_total= ccf1->numeq_total;
  numeq      = ccf1->numeq;
  update     = ccf1->update.a.iv;
  Ax         = ccf1->Ax.a.dv;
  if (istwo)
    Bx       = ccf2->Ax.a.dv;
  else
    Bx       = NULL;
  Ai         = ccf1->Ai.a.iv;
  Ap         = ccf1->Ap.a.iv;
  cdofs      = actpart->pdis[disnum].coupledofs.a.ia;
  ncdofs     = actpart->pdis[disnum].coupledofs.fdim;

  /* put pointers to sendbuffers if any */
#ifdef PARALLEL
  if (ccf1->couple_i_send)
  {
    isend1 = ccf1->couple_i_send->a.ia;
    dsend1 = ccf1->couple_d_send->a.da;
    nsend  = ccf1->couple_i_send->fdim;
    if (istwo)
    {
      isend2 = ccf2->couple_i_send->a.ia;
      dsend2 = ccf2->couple_d_send->a.da;
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
     free surface condition: nd not eqaual numnp*numdf!!! */
#if 0
  if (counter != nd) dserror("assemblage failed due to wrong dof numbering");
#endif
  nd = counter;


  /* now start looping the dofs */
  /* loop over i (the element column) */
  for (i=0; i<nd; i++)
  {
    ii = lm[i];

    /* loop only my own rows */
#ifdef PARALLEL
    if (owner[i]!=myrank) continue;
#endif

    /* check for boundary condition */
    if (ii>=numeq_total) continue;

    /* loop over j (the element row) */
    start         = Ap[ii];
    lenght        = Ap[ii+1]-Ap[ii];
    for (j=0; j<nd; j++)
    {
      jj = lm[j];

      /* check for boundary condition */
      if (jj>=numeq_total) continue;
      index         = find_index(jj,&(Ai[start]),lenght);
      if (index==-1) dserror("dof jj not found in this row ii");
      index        += start;
      Ax[index] += estif[j][i];
      if (istwo)
        Bx[index] += emass[j][i];
    } /* end loop over j */
  }/* end loop over i */


#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}




#ifdef FAST_ASS

/*----------------------------------------------------------------------*/
/*!
  \brief assemble into a ccf matrix (original version)

  This routine assembles one or two element matrices (estiff_global and
  emass_global) into the global matrices in the ccf format.
  It makes use of the information saved in actele->index to determine the
  correct position in the sparse matrix.
  This is faster then searching every time, but consumes a lot of memory!!

  \param actpart   *PARTITION    (i)  the partition we are working on
  \param actsolv   *SOLVAR       (i)  the solver we are using
  \param actintra  *INTRA        (i)  the intra-communicator we do not need
  \param actele    *ELEMENT      (i)  the element we would like to work with
  \param ccf1      *CCF          (i)  one sparse matrix we will assemble into
  \param ccf2      *CCF          (i)  the other sparse matrix we will assemble into

  \author mn
  \date 07/04

*/
/*----------------------------------------------------------------------*/
void  add_ccf_fast(
    struct _PARTITION     *actpart,
    struct _SOLVAR        *actsolv,
    struct _INTRA         *actintra,
    struct _ELEMENT       *actele,
    struct _CCF           *ccf1,
    struct _CCF           *ccf2,
    struct _ARRAY         *elearray1,
    struct _ARRAY         *elearray2)
{

  INT         i,j;                      /* some counter variables */
  INT         istwo=0;
  INT         start,index,lenght;       /* some more special-purpose counters */
  INT         ii,jj;                    /* counter variables for system matrix */
  INT         nd;                       /* size of estif */
  INT         nnz;                      /* number of nonzeros in sparse system matrix */
  INT         numeq_total;              /* total number of equations */
  INT         numeq;                    /* number of equations on this proc */
  INT         myrank;                   /* my intra-proc number */
  INT         nprocs;                   /* my intra- number of processes */
  INT        *Ai;                       /*    "       Ap (column pointer) see UMFPACK manual */
  INT        *Ap;                       /*    "       Ai see UMFPACK manual */
  DOUBLE     *Ax;                       /*    "       Ax see UMFPACK manual */
  DOUBLE     *Bx;                       /*    "       Ax see UMFPACK manual */
  DOUBLE    **estif;                    /* element matrix to be added to system matrix */
  DOUBLE    **emass;                    /* element matrix to be added to system matrix */
  INT        *update;                   /* vector update see AZTEC manual */
  INT       **cdofs;                    /* list of coupled dofs and there owners, see init_assembly */
  INT         ncdofs;                   /* total number of coupled dofs */

#ifdef PARALLEL
  INT       **isend1;                   /* pointer to sendbuffer to communicate coupling conditions */
  DOUBLE    **dsend1;                   /* pointer to sendbuffer to communicate coupling conditions */
  INT       **isend2;                   /* pointer to sendbuffer to communicate coupling conditions */
  DOUBLE    **dsend2;                   /* pointer to sendbuffer to communicate coupling conditions */
  INT         nsend;
#endif

#ifdef DEBUG
  dstrc_enter("add_ccf");
#endif

  /* check whether to assemble one or two matrices */
  if (ccf2) istwo=1;

  /* set some pointers and variables */
  myrank     = actintra->intra_rank;
  nprocs     = actintra->intra_nprocs;
  estif      = elearray1->a.da;
  if (istwo)
    emass      = elearray2->a.da;
  nd         = actele->nd;
  nnz        = ccf1->nnz;
  numeq_total= ccf1->numeq_total;
  numeq      = ccf1->numeq;
  update     = ccf1->update.a.iv;
  Ax         = ccf1->Ax.a.dv;
  if (istwo)
    Bx       = ccf2->Ax.a.dv;
  else
    Bx       = NULL;
  Ai         = ccf1->Ai.a.iv;
  Ap         = ccf1->Ap.a.iv;
  cdofs      = actpart->pdis[disnum].coupledofs.a.ia;
  ncdofs     = actpart->pdis[disnum].coupledofs.fdim;

  /* put pointers to sendbuffers if any */
#ifdef PARALLEL
  if (ccf1->couple_i_send)
  {
    isend1 = ccf1->couple_i_send->a.ia;
    dsend1 = ccf1->couple_d_send->a.da;
    nsend  = ccf1->couple_i_send->fdim;
    if (istwo)
    {
      isend2 = ccf2->couple_i_send->a.ia;
      dsend2 = ccf2->couple_d_send->a.da;
    }
  }
#endif


  /* loop over i (the element column) */
  for (i=0; i<nd; i++)
  {
    ii = actele->locm[i];

    /* loop over j (the element row) */
    start         = Ap[ii];
    lenght        = Ap[ii+1]-Ap[ii];
    for (j=0; j<nd; j++)
    {
      jj = actele->locm[j];
      index = actele->index[i][j];

      if(index >= 0)  /* normal dof */
      {
        Ax[index] += estif[j][i];
        if (istwo)
          Bx[index] += emass[j][i];
      }

      if(index == -1)  /* boundary condition dof */
        continue;

    } /* end loop over j */
  }/* end loop over i */


#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

#endif /* ifdef FAST_ASS */





/*----------------------------------------------------------------------*
 * redundant_ccf()
 *----------------------------------------------------------------------*/
void redundant_ccf(PARTITION *actpart,
    SOLVAR    *actsolv,
    INTRA     *actintra,
    CCF       *ccf1,
    CCF       *ccf2)
{

#ifdef PARALLEL
  INT      imyrank;
  INT      inprocs;

  ARRAY    recv_a;
  DOUBLE  *recv;
#endif

#ifdef DEBUG
  dstrc_enter("redundant_ccf");
#endif


  /*----------------------------------------------------------------------*/
  /*  NOTE:
      In this routine, for a relatively short time the system matrix
      exists 2 times. This takes a lot of memory and may be a
      bottle neck!
      In MPI2 there exists a flag for in-place-Allreduce:

      MPI_Allreduce(MPI_IN_PLACE,
      ucchb->a.a.dv,
      (ucchb->a.fdim)*(ucchb->a.sdim),
      MPI_DOUBLE,
      MPI_SUM,
      actintra->MPI_INTRA_COMM);

      But there is no MPI2 in for HP, yet.
      */
  /*----------------------------------------------------------------------*/
#ifdef PARALLEL
  /*----------------------------------------------------------------------*/
  imyrank = actintra->intra_rank;
  inprocs = actintra->intra_nprocs;
  /*--- very comfortable: the only thing to do is to alreduce the ARRAY a */
  /*                      (all coupling conditions are done then as well) */
  /*--------------------------------------------------- allocate recvbuff */
  recv = amdef("recv_a",&recv_a,ccf1->Ax.fdim,ccf1->Ax.sdim,"DV");
  /*----------------------------------------------------------- Allreduce */
  MPI_Allreduce(ccf1->Ax.a.dv,
      recv,
      (ccf1->Ax.fdim)*(ccf1->Ax.sdim),
      MPI_DOUBLE,
      MPI_SUM,
      actintra->MPI_INTRA_COMM);
  /*----------------------------------------- copy reduced data back to a */
  amcopy(&recv_a,&(ccf1->Ax));
  if (ccf2)
  {
    /*----------------------------------------------------------- Allreduce */
    MPI_Allreduce(ccf2->Ax.a.dv,
        recv,
        (ccf2->Ax.fdim)*(ccf2->Ax.sdim),
        MPI_DOUBLE,
        MPI_SUM,
        actintra->MPI_INTRA_COMM);
    /*----------------------------------------- copy reduced data back to a */
    amcopy(&recv_a,&(ccf2->Ax));
  }
  /*----------------------------------------------------- delete recvbuff */
  amdel(&recv_a);
#endif /*---------------------------------------------- end of PARALLEL */
  /*----------------------------------------------------------------------*/

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of redundant_ccf */



#endif /* ifdef UMFPACK */


