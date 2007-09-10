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

#ifndef CCADISCRET

#ifdef MUMPS_PACKAGE


#include "../headers/standardtypes.h"
#include "../solver/solver.h"


INT cmp_int(const void *a, const void *b );
DOUBLE cmp_double(const void *a, const void *b );



/*----------------------------------------------------------------------*
  |  calculate the mask of an rc_ptr matrix             m.gee 1/02     |
 *----------------------------------------------------------------------*/
void mask_rc_ptr(
    FIELD         *actfield,
    PARTITION     *actpart,
    SOLVAR        *actsolv,
    INTRA         *actintra,
    RC_PTR        *rc_ptr
    )

{
  INT       i;
  INT       numeq;
  INT     **dof_connect;
  ARRAY     bindx_a;
  INT      *bindx;

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
     */
  /*------------------------------------------- put total size of problem */
  rc_ptr->numeq_total = actfield->dis[0].numeq;
  /* count number of eqns on proc and build processor-global couplingdof
     matrix */
  mask_numeq(actfield,actpart,actsolv,actintra,&numeq,0);
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
  dof_connect = (INT**)CCACALLOC(rc_ptr->numeq_total,sizeof(INT*));
  if (!dof_connect) dserror("Allocation of dof_connect failed");
  /*---------------------- make the dof_connect list locally on each proc */
  rc_ptr_nnz_topology(actfield,actpart,actsolv,actintra,rc_ptr,dof_connect);
  /*------------------------------------------ make dof_connect redundant */
  rc_ptr_red_dof_connect(actfield,actpart,actsolv,actintra,rc_ptr,dof_connect);
  /*----------------------------------------------------- allocate arrays */
  /*                                                     see MUMPS manual */
  amdef("rowptr" ,&(rc_ptr->rowptr) ,rc_ptr->numeq+1  ,1,"IV");
  amdef("irn_loc",&(rc_ptr->irn_loc),rc_ptr->nnz      ,1,"IV");
  amdef("jcn_loc",&(rc_ptr->jcn_loc),rc_ptr->nnz      ,1,"IV");
  amdef("A"      ,&(rc_ptr->A_loc)  ,rc_ptr->nnz      ,1,"DV");
  /*------------------------------------------------------ allocate bindx */
  bindx = amdef("bindx",&(bindx_a),(rc_ptr->nnz+1),1,"IV");
  /*---------------------------------------------------------- make bindx */
  rc_ptr_make_bindx(actfield,actpart,actsolv,rc_ptr,dof_connect,bindx);
  /*----------------- make rowptr, irn_loc, jcn_loc from bindx and update */
  rc_ptr_make_sparsity(rc_ptr,bindx);
  /*---------------------------------------- delete the array dof_connect */
  for (i=0; i<rc_ptr->numeq_total; i++)
  {
    if (dof_connect[i]) CCAFREE(dof_connect[i]);
  }
  CCAFREE(dof_connect);
  amdel(&bindx_a);
  /*----------------------------------------------------------------------*/

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of mask_rc_ptr */




/*----------------------------------------------------------------------*
  | make redundant dof_connect array                       m.gee 1/02  |
 *----------------------------------------------------------------------*/
void rc_ptr_red_dof_connect(
    FIELD        *actfield,
    PARTITION    *actpart,
    SOLVAR       *actsolv,
    INTRA        *actintra,
    RC_PTR       *rc_ptr,
    INT         **dof_connect
    )

{
  INT        i,j,counter;
  INT        imyrank;
  INT        inprocs;
  INT        max_dof_connect_send;
  INT        max_dof_connect_recv;
  INT        dof,actdof,lastdof;
  INT        ismaindiag;
  INT        lenght;

  ARRAY      tmps_a;
  INT      **tmps;
  ARRAY      tmpr_a;
  INT      **tmpr;

  INT       *irn;
  INT       *jcn;

#ifdef DEBUG
  dstrc_enter("rc_ptr_red_dof_connect");
#endif

  /*----------------------------------------------------------------------*/
  imyrank = actintra->intra_rank;
  inprocs = actintra->intra_nprocs;
  /*----------------------- Allreduce the total number of nonzero entries */
  rc_ptr->nnz_total=0;
#ifdef PARALLEL
  MPI_Allreduce(&(rc_ptr->nnz),&(rc_ptr->nnz_total),1,MPI_INT,MPI_SUM,actintra->MPI_INTRA_COMM);
#endif
  /*----------------------------- check for largest row in my dof_connect */
  max_dof_connect_send=0;
  max_dof_connect_recv=0;
  for (i=0; i<rc_ptr->numeq_total; i++)
  {
    if (dof_connect[i])
      if (dof_connect[i][0]>max_dof_connect_send)
        max_dof_connect_send=dof_connect[i][0];
  }
#ifdef PARALLEL
  MPI_Allreduce(&max_dof_connect_send,&max_dof_connect_recv,1,MPI_INT,MPI_MAX,actintra->MPI_INTRA_COMM);
#endif
  /*---------------- allocate temporary array to hold global connectivity */
  tmps = amdef("tmp",&tmps_a,rc_ptr->numeq_total,max_dof_connect_recv,"IA");
  amzero(&tmps_a);
  tmpr = amdef("tmp",&tmpr_a,rc_ptr->numeq_total,max_dof_connect_recv,"IA");
  amzero(&tmpr_a);
  /*-------------------------------- put my own dof_connect values to tmp */
  for (i=0; i<rc_ptr->numeq_total; i++)
  {
    if (dof_connect[i])
    {
      for (j=0; j<dof_connect[i][0]; j++)
        tmps[i][j] = dof_connect[i][j];
    }
  }
  /*--------------------------------------------- allreduce the array tmp */
#ifdef PARALLEL
  MPI_Allreduce(tmps[0],tmpr[0],(tmps_a.fdim*tmps_a.sdim),MPI_INT,MPI_SUM,actintra->MPI_INTRA_COMM);
#endif
  /*------------------ allocate the arrays irn_glob jcn_glob on imyrank=0 */
  if (imyrank==0)
  {
    irn = amdef("irn_g",&(rc_ptr->irn_glob),rc_ptr->nnz_total,1,"IV");
    jcn = amdef("jcn_g",&(rc_ptr->jcn_glob),rc_ptr->nnz_total,1,"IV");
    counter=0;
    for (i=0; i<rc_ptr->numeq_total; i++)
    {
      dof        = tmpr[i][2];
      lenght     = tmpr[i][0]-3;
      lastdof    = -1;
      ismaindiag = 0;
      for (j=0; j<lenght; j++)
      {
        actdof = tmpr[i][j+3];
        if (actdof<dof)
        {
          irn[counter] = dof;
          jcn[counter] = actdof;
          lastdof      = actdof;
          counter++;
          continue;
        }
        if (actdof>dof && lastdof<dof)
        {
          irn[counter] = dof;
          jcn[counter] = dof;
          ismaindiag++;
          counter++;
          irn[counter] = dof;
          jcn[counter] = actdof;
          lastdof      = actdof;
          counter++;
          continue;
        }
        if (actdof>dof)
        {
          irn[counter] = dof;
          jcn[counter] = actdof;
          lastdof      = actdof;
          counter++;
          continue;
        }
      }
      if (ismaindiag==0)
      {
        irn[counter] = dof;
        jcn[counter] = dof;
        counter++;
      }

    }
  }/* end of (imyrank==0) */
  /*--------------------------------------------- delete temporary arrays */
  amdel(&tmps_a);
  amdel(&tmpr_a);
  /*----------------------------------------------------------------------*/

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of rc_ptr_red_dof_connect */





/*----------------------------------------------------------------------*
  | allocate update put dofs in update in ascending order   m.gee 1/02 |
 *----------------------------------------------------------------------*/
void  rc_ptr_update(
    FIELD         *actfield,
    PARTITION     *actpart,
    SOLVAR        *actsolv,
    INTRA         *actintra,
    RC_PTR        *rc_ptr
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
  dstrc_enter("rc_ptr_update");
#endif

  /*----------------------------------------------------------------------*/
  imyrank = actintra->intra_rank;
  inprocs = actintra->intra_nprocs;
  /*------------------ make a local copy of the array actpart->coupledofs */
  memset(&coupledofs, 0, sizeof(ARRAY));
  if (actpart->pdis[dis].coupledofs.Typ != cca_XX)
    am_alloc_copy(&(actpart->pdis[0].coupledofs),&coupledofs);
  /*------------------------------------- loop the nodes on the partition */
  update = rc_ptr->update.a.iv;
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
  if (counter != rc_ptr->numeq) dserror("Number of dofs in update wrong");
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
} /* end of rc_ptr_update */




/*----------------------------------------------------------------------*
  | calculate number of nonzero entries and dof topology    m.gee 1/02 |
 *----------------------------------------------------------------------*/
void  rc_ptr_nnz_topology(
    FIELD         *actfield,
    PARTITION    *actpart,
    SOLVAR       *actsolv,
    INTRA        *actintra,
    RC_PTR       *rc_ptr,
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
    qsort((INT*)(&(dof_connect[dof][3])), dof_connect[dof][0]-3, sizeof(INT), cmp_int);
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
  | make the DMSR vector bindx                              m.gee 1/02 |
  | for format see Aztec manual                                        |
 *----------------------------------------------------------------------*/
void  rc_ptr_make_bindx(
    FIELD         *actfield,
    PARTITION     *actpart,
    SOLVAR        *actsolv,
    RC_PTR        *rc_ptr,
    INT          **dof_connect,
    INT           *bindx
    )

{
  INT        i,j;
  INT        count1,count2;
  INT        dof;

#ifdef DEBUG
  dstrc_enter("rc_ptr_make_bindx");
#endif

  /*-------------------------------------------------------------do bindx */
  count1=0;
  count2=rc_ptr->numeq+1;
  for (i=0; i<rc_ptr->update.fdim; i++)
  {
    dof = rc_ptr->update.a.iv[i];
    bindx[count1] = count2;
    count1++;
    for (j=3; j<dof_connect[dof][0]; j++)
    {
      bindx[count2] = dof_connect[dof][j];
      count2++;
    }
  }
  bindx[rc_ptr->numeq] = rc_ptr->nnz+1;
  /*----------------------------------------------------------------------*/

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of rc_ptr_make_bindx */




/*----------------------------------------------------------------------*
  | make the vectors                                        m.gee 1/02 |
  | irn_loc, jcn_loc, rowptr from update and bindx                     |
 *----------------------------------------------------------------------*/
void  rc_ptr_make_sparsity(
    RC_PTR        *rc_ptr,
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
  dstrc_enter("rc_ptr_make_sparsity");
#endif

  /*----------------------------------------------------------------------*/
  numeq       = rc_ptr->numeq;
  numeq_total = rc_ptr->numeq_total;
  nnz         = rc_ptr->nnz;
  update      = rc_ptr->update.a.iv;
  irn         = rc_ptr->irn_loc.a.iv;
  jcn         = rc_ptr->jcn_loc.a.iv;
  rptr        = rc_ptr->rowptr.a.iv;
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
} /* end of rc_ptr_make_sparsity */





/*----------------------------------------------------------------------*
  | routine to assemble element array to global rcptr-matrix           |
  | in parallel,taking care of coupling conditions                     |
  |                                                                    |
  |                                                                    |
  |                                                        m.gee 1/02  |
 *----------------------------------------------------------------------*/
void  add_rc_ptr(
    struct _PARTITION     *actpart,
    struct _SOLVAR        *actsolv,
    struct _INTRA         *actintra,
    struct _ELEMENT       *actele,
    struct _RC_PTR        *rc_ptr1,
    struct _RC_PTR        *rc_ptr2,
    struct _ARRAY         *elearray1,
    struct _ARRAY         *elearray2
    )

{

  INT         i,j,k,l,counter;       /* some counter variables */
  INT         istwo=0;
  INT         start,index,lenght;    /* some more special-purpose counters */
  INT         ii,jj;                 /* counter variables for system matrix */
  INT         ii_iscouple;           /* flag whether ii is a coupled dof */
  INT         ii_owner;              /* who is owner of dof ii -> procnumber */
  INT         ii_index;              /* place of ii in dmsr format */
  INT         jj_index;              /* place of jj in dmsr format */
  INT         nd,ndnd;               /* size of estif */
  INT         nnz;                   /* number of nonzeros in sparse system matrix */
  INT         numeq_total;           /* total number of equations */
  INT         numeq;                 /* number of equations on this proc */
  INT         lm[MAXDOFPERELE];      /* location vector for this element */
  INT         owner[MAXDOFPERELE];   /* the owner of every dof */
  INT         myrank;                /* my intra-proc number */
  INT         nprocs;                /* my intra- number of processes */
  DOUBLE    **estif;                 /* element matrix to be added to system matrix */
  DOUBLE    **emass;                 /* element matrix to be added to system matrix */
  INT        *update;                /* vector update see AZTEC manual */
  DOUBLE     *A_loc;                 /*    "       A_loc see MUMPS manual */
  DOUBLE     *B_loc;                 /*    "       A_loc see MUMPS manual */
  INT        *irn;                   /*    "       irn see MUMPS manual */
  INT        *jcn;                   /*    "       jcn see MUMPS manual */
  INT        *rowptr;                /*    "       rowptr see rc_ptr structure */
  INT       **cdofs;                 /* list of coupled dofs and there owners, see init_assembly */
  INT         ncdofs;                /* total number of coupled dofs */
  INT       **isend1;                /* pointer to sendbuffer to communicate coupling conditions */
  DOUBLE    **dsend1;                /* pointer to sendbuffer to communicate coupling conditions */
  INT       **isend2;                /* pointer to sendbuffer to communicate coupling conditions */
  DOUBLE    **dsend2;                /* pointer to sendbuffer to communicate coupling conditions */
  INT         nsend;

#ifdef DEBUG
  dstrc_enter("add_rc_ptr");
#endif

  /*----------------------- check whether to assemble one or two matrices */
  if (rc_ptr2) istwo=1;
  /*------------------------------------- set some pointers and variables */
  myrank     = actintra->intra_rank;
  nprocs     = actintra->intra_nprocs;
  estif      = elearray1->a.da;
  if (istwo)
    emass      = elearray1->a.da;
  else
    emass      = NULL;
  nd         = actele->numnp * actele->node[0]->numdf;
  ndnd       = nd*nd;
  nnz        = rc_ptr1->nnz;
  numeq_total= rc_ptr1->numeq_total;
  numeq      = rc_ptr1->numeq;
  update     = rc_ptr1->update.a.iv;
  A_loc      = rc_ptr1->A_loc.a.dv;
  if (istwo)
    B_loc      = rc_ptr2->A_loc.a.dv;
  irn        = rc_ptr1->irn_loc.a.iv;
  jcn        = rc_ptr1->jcn_loc.a.iv;
  rowptr     = rc_ptr1->rowptr.a.iv;
  cdofs      = actpart->pdis[0].coupledofs.a.ia;
  ncdofs     = actpart->pdis[0].coupledofs.fdim;
  /*---------------------------------- put pointers to sendbuffers if any */
#ifdef PARALLEL
  if (rc_ptr1->couple_i_send)
  {
    isend1 = rc_ptr1->couple_i_send->a.ia;
    dsend1 = rc_ptr1->couple_d_send->a.da;
    nsend  = rc_ptr1->couple_i_send->fdim;
    if (istwo)
    {
      isend2 = rc_ptr2->couple_i_send->a.ia;
      dsend2 = rc_ptr2->couple_d_send->a.da;
    }
  }
#endif
  /*---------------------------------------------- make location vector lm*/
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
  if (counter != nd) dserror("assemblage failed due to wrong dof numbering");
  /*========================================== now start looping the dofs */
  /*======================================= loop over i (the element row) */
  ii_iscouple = 0;
  ii_owner    = myrank;
  for (i=0; i<nd; i++)
  {
    ii = lm[i];
    /*-------------------------------------------- loop only my own rows */
#ifdef PARALLEL
    if (owner[i]!=myrank) continue;
#endif
    /*------------------------------------- check for boundary condition */
    if (ii>=numeq_total) continue;
    /*------------------------------------- check for coupling condition */
#ifdef PARALLEL
    if (ncdofs)
    {
      ii_iscouple = 0;
      ii_owner    = -1;
      add_rcptr_checkcouple(ii,cdofs,ncdofs,&ii_iscouple,&ii_owner,nprocs);
    }
#endif
    /*-------------------- ii is not a coupled dofs or I am master owner */
    if (!ii_iscouple || ii_owner==myrank)
    {
      ii_index      = find_index(ii,update,numeq);
      if (ii_index==-1) dserror("dof ii not found on this proc");
      start         = rowptr[ii_index];
      lenght        = rowptr[ii_index+1]-rowptr[ii_index];
    }
    /*================================= loop over j (the element column) */
    /*                            This is the full unsymmetric version ! */
    for (j=0; j<nd; j++)
    {
      jj = lm[j];
      /*---------------------------------- check for boundary condition */
      if (jj>=numeq_total) continue;
      /*---------------------------------- check for coupling condition */
      /*
         coupling condition for jj is not checked, because we add to
         row ii here, which must also hold the coupled columns jj
         */
      /*======================================== do main-diagonal entry */
      /*                (either not a coupled dof or I am master owner) */
      if (!ii_iscouple || ii_owner==myrank)
      {
        index         = find_index(jj,&(jcn[start]),lenght);
        if (index==-1) dserror("dof jj not found in this row ii");
        index        += start;
        A_loc[index] += estif[i][j];
        if (istwo)
          B_loc[index] += emass[i][j];
      }
      /*======================================== do main-diagonal entry */
      /*                           (a coupled dof and I am slave owner) */
      else
      {
        add_rcptr_sendbuff(ii,jj,i,j,ii_owner,isend1,dsend1,estif,nsend);
        if (istwo)
          add_rcptr_sendbuff(ii,jj,i,j,ii_owner,isend2,dsend2,emass,nsend);
      }


    } /* end loop over j */
  }/* end loop over i */
  /*----------------------------------------------------------------------*/

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of add_rc_ptr */





/*----------------------------------------------------------------------*
  | checks coupling for the add_rcptr routine               m.gee 9/01 |
 *----------------------------------------------------------------------*/
void add_rcptr_checkcouple(
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
  dstrc_enter("add_rcptr_checkcouple");
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
} /* end of add_rcptr_checkcouple */





/*----------------------------------------------------------------------*
  | fill sendbuffer isend and dsend                         m.gee 1/02 |
 *----------------------------------------------------------------------*/
void add_rcptr_sendbuff(
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
  dstrc_enter("add_rcptr_sendbuff");
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
} /* end of add_rcptr_sendbuff */





/*----------------------------------------------------------------------*
  | exchange coupled dofs and add to row/column ptr matrix  m.gee 1/02 |
 *----------------------------------------------------------------------*/
void exchange_coup_rc_ptr(
    PARTITION     *actpart,
    SOLVAR        *actsolv,
    INTRA         *actintra,
    RC_PTR        *rc_ptr
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
  dstrc_enter("exchange_coup_rc_ptr");
#endif

  /*----------------------------------------------------------------------*/
#ifdef PARALLEL
  /*----------------------------------------------------------------------*/
  imyrank = actintra->intra_rank;
  inprocs = actintra->intra_nprocs;
  ACTCOMM = &(actintra->MPI_INTRA_COMM);
  /*---------------------------------------- set some pointers and values */
  numsend     = rc_ptr->numcoupsend;
  numrecv     = rc_ptr->numcouprecv;
  A_loc      = rc_ptr->A_loc.a.dv;
  irn        = rc_ptr->irn_loc.a.iv;
  jcn        = rc_ptr->jcn_loc.a.iv;
  rowptr     = rc_ptr->rowptr.a.iv;
  update      = rc_ptr->update.a.iv;
  numeq_total = rc_ptr->numeq_total;
  numeq       = rc_ptr->numeq;
  if (rc_ptr->couple_i_send) isend = rc_ptr->couple_i_send->a.ia;
  if (rc_ptr->couple_d_send) dsend = rc_ptr->couple_d_send->a.da;
  if (rc_ptr->couple_i_recv) irecv = rc_ptr->couple_i_recv->a.ia;
  if (rc_ptr->couple_d_recv) drecv = rc_ptr->couple_d_recv->a.da;
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
      A_loc[index] += drecv[i][jj];
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
} /* end of exchange_coup_rc_ptr */




#endif /* ifdef MUMPS_PACKAGE */


#endif
