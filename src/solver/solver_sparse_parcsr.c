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
#ifdef HYPRE_PACKAGE

#include "../headers/standardtypes.h"
#include "../solver/solver.h"


INT cmp_int(const void *a, const void *b );
DOUBLE cmp_double(const void *a, const void *b );




/*----------------------------------------------------------------------*
  |  calculate the mask of an parcsr matrix              m.gee 10/01     |
 *----------------------------------------------------------------------*/
void  mask_parcsr(
    FIELD         *actfield,
    PARTITION     *actpart,
    SOLVAR        *actsolv,
    INTRA         *actintra,
    H_PARCSR      *parcsr
    )

{
  INT       i;
  INT       numeq;
  INT     **dof_connect;

#ifdef DEBUG
  dstrc_enter("mask_parcsr");
#endif

  /*----------------------------------------------------------------------*/
  /* remember some facts:
     PARTITION is different on every proc.
     H_PARCSR will be different on every proc
     FIELD is the same everywhere
     */
  /*------------------------------------------- put total size of problem */
  parcsr->numeq_total = actfield->dis[0].numeq;
  /* count number of eqns on proc and build processor-global couplingdof
     matrix */
  mask_numeq(actfield,actpart,actsolv,actintra,&numeq,0);
  parcsr->numeq = numeq;
  /*---------------------------------------------- allocate vector update */
  amdef("update",&(parcsr->update),numeq,1,"IV");
  amzero(&(parcsr->update));
  /*--------------------------------put dofs in update in ascending order */
  parcsr_update(actfield,actpart,actsolv,actintra,parcsr);
  /*------------------------ count number of nonzero entries on partition
    and calculate dof connectivity list */
  /*
     dof_connect[i][0] = lenght of dof_connect[i]
     dof_connect[i][1] = iscoupled ( 1 or 2 )
     dof_connect[i][2] = dof
     dof_connect[i][ 2..dof_connect[i][0]-1 ] = connected dofs exluding itself
     */
  dof_connect = (INT**)CCACALLOC(parcsr->numeq_total,sizeof(INT*));
  if (!dof_connect) dserror("Allocation of dof_connect failed");
  parcsr_nnz_topology(actfield,actpart,actsolv,actintra,parcsr,dof_connect);
  /*---------------------------------------------- allocate bindx and val */
  amdef("bindx",&(parcsr->bindx),(parcsr->nnz+1),1,"IV");
  /*---------------------------------------------------------- make bindx */
  parcsr_make_bindx(actfield,actpart,actsolv,parcsr,dof_connect);
  /*------------------------------- make the permutation vector of update */
  /*---------------- vector update is redefined to a global array in here */
  parcsr_update_perm(actintra,parcsr);
  /*---------------------------------------- delete the array dof_connect */
  for (i=0; i<parcsr->numeq_total; i++)
  {
    if (dof_connect[i]) CCAFREE(dof_connect[i]);
  }
  CCAFREE(dof_connect);
  /*----------------------------------------------------------------------*/

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of mask_parcsr */




/*----------------------------------------------------------------------*
  |  allocate update put dofs in update in ascending order  m.gee 10/01  |
 *----------------------------------------------------------------------*/
void parcsr_update(
    FIELD         *actfield,
    PARTITION     *actpart,
    SOLVAR        *actsolv,
    INTRA         *actintra,
    H_PARCSR      *parcsr
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
  dstrc_enter("parcsr_update");
#endif

  /*----------------------------------------------------------------------*/
  imyrank = actintra->intra_rank;
  inprocs = actintra->intra_nprocs;
  /*------------------ make a local copy of the array actpart->coupledofs */
  memset(&coupledofs, 0, sizeof(ARRAY));
  if (actpart->pdis[dis].coupledofs.Typ != cca_XX)
    am_alloc_copy(&(actpart->pdis[0].coupledofs),&coupledofs);
  /*------------------------------------- loop the nodes on the partition */
  update = parcsr->update.a.iv;
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
  if (counter != parcsr->numeq)
    dserror("Number of dofs in ParCSR-vector update wrong");
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
} /* end of parcsr_update */




/*----------------------------------------------------------------------*
  |  allocate perm and make permutation of update           m.gee 10/01  |
 *----------------------------------------------------------------------*/
void parcsr_update_perm(
    INTRA         *actintra,
    H_PARCSR      *parcsr
    )

{
  INT       i,j;
  INT       inprocs;
  INT       imyrank;
  INT       prevdofs;
  INT       maxdofs;
  INT       minusone=-1;
  INT      *p_sizes;
  INT     **perm;
  ARRAY     tmp;
  ARRAY     recv_a;

#ifdef PARALLEL
  INT      *recvbuff;
#endif

#ifdef DEBUG
  dstrc_enter("parcsr_update_perm");
#endif

  /*----------------------------------------------------------------------*/
  /*
     The HYPRE solvers want to have the dofs in ascending order and contigous,
     starting with proc 0 (See HYPRE manuals). So it is necesary to perform a
     permutation of the dof numbers especially for the solver.
     The vector update is already ordered in ascending order on each proc,
     but it does not necessariely start with the small numbers on proc 0.
     In addition, update is not contigous.
     */
  /*----------------------------------------------------------------------*/
  inprocs  = actintra->intra_nprocs;
  imyrank  = actintra->intra_rank;
  p_sizes = amdef("p_sizes",&(parcsr->perm_sizes),inprocs,1,"IV");
  amzero(&(parcsr->perm_sizes));
  p_sizes[imyrank] = parcsr->numeq;

#ifdef PARALLEL
  recvbuff = (INT*)CCACALLOC(inprocs,sizeof(INT));
  if (!recvbuff) dserror("Allocation of memory failed");
  MPI_Allreduce(p_sizes,recvbuff,inprocs,MPI_INT,MPI_SUM,actintra->MPI_INTRA_COMM);
  for (i=0; i<inprocs; i++) p_sizes[i] = recvbuff[i];
  CCAFREE(recvbuff);
#endif

  /*--------now every proc knows the sizes of all procs and can calculate */
  /*                                    the range of permuted dof numbers */
  /*---------------------------- find the maximum number of dofs per proc */
  maxdofs=0;
  for (i=0; i<inprocs; i++)
  {
    if (p_sizes[i]>maxdofs) maxdofs = p_sizes[i];
  }
  /*-------------------------------------------- allocate the matrix perm */
  perm = amdef("perm",&(parcsr->perm),inprocs,maxdofs,"IA");
  aminit(&(parcsr->perm),&minusone);
  /*--------- make the permutation vectors for each proc including myself */
  prevdofs=0;
  for (i=0; i<inprocs; i++)
  {
    for (j=0; j<p_sizes[i]; j++)
    {
      perm[i][j]=prevdofs;
      prevdofs++;
    }
  }
  /*---------- now redefine update to hold the update vector of all procs */
  am_alloc_copy(&(parcsr->update),&tmp);
  amdel(&(parcsr->update));
  amdef("update",&(parcsr->update),inprocs,maxdofs,"IA");
  amdef("update",&recv_a          ,inprocs,maxdofs,"IA");
  amzero(&(parcsr->update));
  amzero(&recv_a);
  for (i=0; i<tmp.fdim; i++)
  {
    parcsr->update.a.ia[imyrank][i]=tmp.a.iv[i];
  }
  amdel(&tmp);

#ifdef PARALLEL
  MPI_Allreduce(parcsr->update.a.ia[0],
      recv_a.a.ia[0],
      (parcsr->update.fdim)*(parcsr->update.sdim),
      MPI_INT,
      MPI_SUM,
      actintra->MPI_INTRA_COMM);
  amcopy(&recv_a,&(parcsr->update));
#endif


#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of parcsr_update_perm */




/*----------------------------------------------------------------------*
  |  calculate number of nonzero entries and dof topology   m.gee 10/01  |
 *----------------------------------------------------------------------*/
void  parcsr_nnz_topology(
    FIELD         *actfield,
    PARTITION     *actpart,
    SOLVAR        *actsolv,
    INTRA         *actintra,
    H_PARCSR      *parcsr,
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

#ifdef PARALLEL
  MPI_Status status;
#endif

#ifdef DEBUG
  dstrc_enter("parcsr_nnz_topology");
#endif

  /*----------------------------------------------------------------------*/
  imyrank = actintra->intra_rank;
  inprocs = actintra->intra_nprocs;
  /*----------------------------------------------------------- shortcuts */
  parcsr->nnz=0;
  numeq  = parcsr->numeq;
  update = parcsr->update.a.iv;
  for (i=0; i<parcsr->numeq_total; i++) dof_connect[i]=NULL;
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
    centernode=NULL;
    dof_find_centernode(dof,actpart,&centernode);
    dsassert(centernode!=NULL,"cannot find centernode for patch");
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
            if (counter>=dofpatch.fdim)
              amredef(&dofpatch,dofpatch.fdim+500,1,"IV");
            dofpatch.a.iv[counter] = actnode->dof[l];
            counter++;
          }
        }
      }
    }
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
                if (counter>=dofpatch.fdim)
                  amredef(&dofpatch,dofpatch.fdim+500,1,"IV");
                dofpatch.a.iv[counter] = actnode->dof[l];
                counter++;
              }
            }
          }
        }
      }
    }/* end of making dofpatch */
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
  for (i=0; i<parcsr->update.fdim; i++)
  {
    dof = parcsr->update.a.iv[i];
    nnz += (dof_connect[dof][0]-2);
  }
  parcsr->nnz=nnz;
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
} /* end of parcsr_nnz_topology */




/*----------------------------------------------------------------------*
  |  make the vector bindx                                  m.gee 10/01  |
  | for format see Aztec manual                                          |
 *----------------------------------------------------------------------*/
void parcsr_make_bindx(
    FIELD         *actfield,
    PARTITION     *actpart,
    SOLVAR        *actsolv,
    H_PARCSR      *parcsr,
    INT          **dof_connect
    )

{
  INT        i,j;
  INT        count1,count2;
  INT        dof;

#ifdef DEBUG
  dstrc_enter("parcsr_make_bindx");
#endif

  /*-------------------------------------------------------------do bindx */
  count1=0;
  count2=parcsr->numeq+1;
  for (i=0; i<parcsr->update.fdim; i++)
  {
    dof = parcsr->update.a.iv[i];
    parcsr->bindx.a.iv[count1] = count2;
    count1++;
    for (j=3; j<dof_connect[dof][0]; j++)
    {
      parcsr->bindx.a.iv[count2] = dof_connect[dof][j];
      count2++;
    }
  }
  parcsr->bindx.a.iv[parcsr->numeq] = parcsr->nnz+1;
  /*----------------------------------------------------------------------*/

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of parcsr_make_bindx */





/*----------------------------------------------------------------------*
  |  routine to assemble element array to global PARCSR-matrix           |
  |  in parallel and sequentiell,taking care of coupling conditions      |
  |                                                                      |
  |                                                                      |
  |                                                         m.gee 10/01  |
 *----------------------------------------------------------------------*/
void  add_parcsr(
    struct _PARTITION     *actpart,
    struct _SOLVAR        *actsolv,
    struct _INTRA         *actintra,
    struct _ELEMENT       *actele,
    struct _H_PARCSR      *parcsr,
    struct _ARRAY         *elearray1
    )

{
  INT         i,j,counter;

  INT         ii;
  INT         iiperm;
  INT         ii_index;
  INT         ii_iscouple;              /* flag whether ii is a coupled dof */
  INT         ii_owner;                 /* who is owner of dof ii -> procnumber */
  INT         jj;
  INT         jjperm;
  INT         jj_index;
  INT         jj_iscouple;              /* flag whether ii is a coupled dof */
  INT         jj_owner;                 /* who is owner of dof ii -> procnumber */

  INT         err;
  const INT   nrows=1;

  INT         rows[1];
  INT         ncols[1];
  INT         colcounter;
  INT         cols[MAX_NNZPERROW];
  DOUBLE      values[MAX_NNZPERROW];

  INT         nd;
  INT         numeq_total;
  INT         numeq;
  INT         lm[MAXDOFPERELE];         /* location vector for this element */
  INT         owner[MAXDOFPERELE];      /* the owner of every dof */

  INT         myrank;
  INT         nprocs;

  INT       **isend;
  DOUBLE    **dsend;
  INT         nsend;

  DOUBLE    **estif;
  INT       **cdofs;
  INT         ncdofs;
  INT       **perm;
  INT        *perm_sizes;
  INT       **update;

#ifdef DEBUG
  dstrc_enter("add_parcsr");
#endif

  /*------------------------------------- set some pointers and variables */
  myrank           = actintra->intra_rank;
  nprocs           = actintra->intra_nprocs;

  estif            = elearray1->a.da;
  nd               = actele->numnp * actele->node[0]->numdf;
  numeq_total      = parcsr->numeq_total;
  numeq            = parcsr->numeq;
  cdofs            = actpart->pdis[0].coupledofs.a.ia;
  ncdofs           = actpart->pdis[0].coupledofs.fdim;

  perm             = parcsr->perm.a.ia;
  perm_sizes       = parcsr->perm_sizes.a.iv;
  update           = parcsr->update.a.ia;

  /*---------------------------------- put pointers to sendbuffers if any */
#ifdef PARALLEL
  if (parcsr->couple_i_send)
  {
    isend = parcsr->couple_i_send->a.ia;
    dsend = parcsr->couple_d_send->a.da;
    nsend = parcsr->couple_i_send->fdim;
  }
#else
  isend = NULL;
  dsend = NULL;
  nsend = NULL;
#endif

  /*---------------------------------------------- make location vector lm*/
  counter=0;
  for (i=0; i<actele->numnp; i++)
  {
    for (j=0; j<actele->node[i]->numdf; j++)
    {
      /* location matrix */
      lm[counter]    = actele->node[i]->dof[j];
      /* proc that owns this dof */
      owner[counter] = actele->node[i]->proc;
      counter++;
    }/* end of loop over dofs */
  }/* end of loop over element nodes */
  /*========================================== now start looping the dofs */
  /*======================================= loop over i (the element row) */
  ii_iscouple = 0;
  ii_owner    = myrank;
  for (i=0; i<nd; i++)
  {
    /*------------------------------------------------ set row indize ii */
    ii     = lm[i];
    /*-------------------------------------------- loop only my own rows */
    /*------------------- I am not master owner and I am not slave owner */
    if (owner[i]!=myrank) continue;
    /*------------------------------------- check for boundary condition */
    if (ii>=numeq_total) continue;
    /*------------------------------------- check for coupling condition */
    /*                                                (only in parallel) */
    if (ncdofs)
    {
      ii_iscouple =  0;
      ii_owner    = -1;
      add_parcsr_checkcouple(ii,cdofs,ncdofs,&ii_iscouple,&ii_owner,nprocs);
    }
    /*--------------------------------------- find the permutation of ii */
    if (!ii_iscouple)
    {
      ii_index = find_index(ii,&(update[owner[i]][0]),perm_sizes[owner[i]]);
      if (ii_index==-1) dserror("dof ii not found on expected proc");
      iiperm = perm[owner[i]][ii_index];
    }
    else
    {
      ii_index = find_index(ii,&(update[ii_owner][0]),perm_sizes[ii_owner]);
      if (ii_index==-1) dserror("dof ii not found on expected proc");
      iiperm = perm[ii_owner][ii_index];
    }
    /*------------------------------------------------- prepare assembly */
    rows[0]  = iiperm;
    /*================================= loop over j (the element column) */
    jj_iscouple =  0;
    jj_owner    =  myrank;

    colcounter = 0;
    /*------------------ if ii is not coupled or I am master owner of ii */
    /*                                  this is the normal standard case */
    if ( ii_iscouple==0 || ii_owner==myrank)
    {
      for (j=0; j<nd; j++)
      {
        /*----------------------------------------- check for overflow */
        if (colcounter>=MAX_NNZPERROW)
          dserror("Overflow in element assembly, increase MAX_NNZPERROW in defines.h");
        /*-------------------------------------- set column indizee jj */
        jj = lm[j];
        /*------------------------------- check for boundary condition */
        if (jj>=numeq_total) continue;
        /*------------------------------- check for coupling condition */
        /*                                          (only in parallel) */
        if (ncdofs)
        {
          jj_iscouple = 0;
          jj_owner    = -1;
          add_parcsr_checkcouple(jj,cdofs,ncdofs,&jj_iscouple,&jj_owner,nprocs);
        }
        /*--------------------------------- find the permutation of jj */
        /*                    NOTE: jj is not necessarily on this proc */
        if (!jj_iscouple)
        {
          jj_index = find_index(jj,update[owner[j]],perm_sizes[owner[j]]);
          if (jj_index==-1) dserror("dof jj not found on expected proc");
          jjperm = perm[owner[j]][jj_index];
        }
        else
        {
          jj_index = find_index(jj,update[jj_owner],perm_sizes[owner[j]]);
          if (jj_index==-1) dserror("dof jj not found on expected proc");
          jjperm = perm[jj_owner][jj_index];
        }
        /*------------------------------------------- prepare assembly */
        cols[colcounter]   = jjperm;
        values[colcounter] = estif[i][j];
        colcounter++;
        /*-------------------------------------------------------------*/
      }/* end loop over j */
      /*------------------------- prepare rest of assembly for this row */
      ncols[0]=colcounter;
      mg_sort(cols,colcounter,NULL,values);
      /*-------------------------------------------------- assemble row */
      /* for detailed description of this assembly format see HYPRE manual */
      err=HYPRE_IJMatrixAddToValues(
          parcsr->ij_matrix,
          nrows,
          ncols,
          rows,
          cols,
          values
          );
      if (err) dserror("Error occured adding to ParCSR matrix");
    }/* end of adding myself */
    /*----------------------------- if ii is coupled and I am slave owner */
    /*----- add to the sendbuffer, sendbuffer is initialized in calelm */
    /*            this is the parallel meets coupling conditions exeption */
    else
    {
      for (j=0; j<nd; j++)
      {
        /*-------------------------------------- set column indizee jj */
        jj = lm[j];
        /*------------------------------- check for boundary condition */
        if (jj>=numeq_total) continue;
        /*---------------------------------------add to the sendbuffer */
        add_parcsr_sendbuff(ii,jj,i,j,ii_owner,isend,dsend,estif,nsend);
      }/* end loop over j */
    }/* end of ii is coupled and I am slave owner */
    /*-------------------------------------------------------------------*/
  }/* end loop over i */
  /*----------------------------------------------------------------------*/

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of add_parcsr */




/*----------------------------------------------------------------------*
  |  fill sendbuffer isend and dsend                          m.gee 10/01|
 *----------------------------------------------------------------------*/
void add_parcsr_sendbuff(
    INT       ii,
    INT       jj,
    INT       i,
    INT       j,
    INT       ii_owner,
    INT     **isend,
    DOUBLE  **dsend,
    DOUBLE  **estif,
    INT       numsend
    )

{
  INT         k;

#ifdef DEBUG
  dstrc_enter("add_parcsr_sendbuff");
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
} /* end of add_parcsr_sendbuff */




/*----------------------------------------------------------------------*
  |  checks coupling for the add_msr routine                   m.gee 9/01|
 *----------------------------------------------------------------------*/
void add_parcsr_checkcouple(
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
  dstrc_enter("add_parcsr_checkcouple");
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
} /* end of add_parcsr_checkcouple */





/*----------------------------------------------------------------------*
  |  exchange coupled dofs and add to hypre matrix            m.gee 10/01|
 *----------------------------------------------------------------------*/
void exchange_coup_parcsr(
    PARTITION     *actpart,
    SOLVAR        *actsolv,
    INTRA         *actintra,
    H_PARCSR      *parcsr
    )

{

#ifdef PARALLEL
  INT            i,j,k;
  INT            ii,ii_index;
  INT            jj,jj_index;
  INT            start;
  INT            lenght;
  INT            tag;
  INT            source;
  INT            owner;
  INT            numeq,numeq_total;
  INT            numsend;
  INT            numrecv;
  INT           *bindx;
  INT          **update;
  INT          **perm;
  INT           *perm_sizes;
  INT          **isend;
  DOUBLE       **dsend;
  INT          **irecv;
  DOUBLE       **drecv;
  INT            imyrank;
  INT            inprocs;

  INT            err;

  INT            iiperm;
  INT            jjperm;
  const INT      nrows=1;
  INT            rows[1];
  INT            ncols[1];
  INT            colcounter;
  INT            cols[MAX_NNZPERROW];
  DOUBLE         values[MAX_NNZPERROW];

  MPI_Status    *irecv_status;
  MPI_Status    *drecv_status;

  MPI_Request   *isendrequest;
  MPI_Request   *dsendrequest;

  MPI_Comm      *ACTCOMM;
#endif

#ifdef DEBUG
  dstrc_enter("exchange_coup_parcsr");
#endif

  /*----------------------------------------------------------------------*/
#ifdef PARALLEL
  /*----------------------------------------------------------------------*/
  imyrank = actintra->intra_rank;
  inprocs = actintra->intra_nprocs;
  ACTCOMM = &(actintra->MPI_INTRA_COMM);
  /*---------------------------------------- set some pointers and values */
  numsend     = parcsr->numcoupsend;
  numrecv     = parcsr->numcouprecv;
  bindx       = parcsr->bindx.a.iv;
  update      = parcsr->update.a.ia;
  perm        = parcsr->perm.a.ia;
  perm_sizes  = parcsr->perm_sizes.a.iv;
  numeq_total = parcsr->numeq_total;
  numeq       = parcsr->numeq;
  if (parcsr->couple_i_send) isend   = parcsr->couple_i_send->a.ia;
  if (parcsr->couple_d_send) dsend   = parcsr->couple_d_send->a.da;
  if (parcsr->couple_i_recv) irecv   = parcsr->couple_i_recv->a.ia;
  if (parcsr->couple_d_recv) drecv   = parcsr->couple_d_recv->a.da;
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
    ii_index = find_index(ii,update[imyrank],perm_sizes[imyrank]);
    if (ii_index==-1) dserror("dof ii not found on this proc");
    /*---------------- add to my piece of system matrix in parcsr format */
    iiperm  = perm[imyrank][ii_index];
    rows[0] = iiperm;
    colcounter=0;
    /*------------------------------------- do main diagonal entry first */
    cols[colcounter]   = iiperm;
    values[colcounter] = drecv[i][ii];
    colcounter++;
    /*------------------------------------------ do off-diagonal entries */
    start  = bindx[ii_index];
    lenght = bindx[ii_index+1]-bindx[ii_index];
    for (j=0; j<lenght; j++)
    {
      jj                 = bindx[start+j];
      jj_index           = find_index(jj,update[imyrank],perm_sizes[imyrank]);
      owner              = imyrank;
      /*------------------------ the dof jj is not updated on this proc */
      if (jj_index==-1)
      {
        for (k=0; k<inprocs; k++)
        {
          jj_index = find_index(jj,update[k],perm_sizes[k]);
          if (jj_index != -1)
          {
            owner = k;
            break;
          }
        }
      }
      if (jj_index==-1) dserror("dof jj not found on expected proc");
      jjperm             = perm[owner][jj_index];
      /*----------------------------------------- check for overflow */
      if (colcounter>=MAX_NNZPERROW)
        dserror("Overflow in element assembly, increase MAX_NNZPERROW in defines.h");
      cols[colcounter]   = jjperm;
      values[colcounter] = drecv[i][jj];
      colcounter++;
    }
    ncols[0] = colcounter;
    err=HYPRE_IJMatrixAddToValues(
        parcsr->ij_matrix,
        nrows,
        ncols,
        rows,
        cols,
        values
        );
    if (err) dserror("Error occured adding to ParCSR matrix");
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
} /* end of exchange_coup_parcsr */


#endif /* ifdef HYPRE_PACKAGE */

#endif
