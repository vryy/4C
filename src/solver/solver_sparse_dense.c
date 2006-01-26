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


#include "../headers/standardtypes.h"
#include "../solver/solver.h"



INT cmp_int(const void *a, const void *b );
DOUBLE cmp_double(const void *a, const void *b );




/*----------------------------------------------------------------------*
 |  calculate the mask of an msr matrix                  m.gee 5/01     |
 *----------------------------------------------------------------------*/
void mask_dense(
    FIELD         *actfield,
    PARTITION     *actpart,
    SOLVAR        *actsolv,
    INTRA         *actintra,
    DENSE         *dense
    )

{
  INT       numeq;

#ifdef DEBUG
  dstrc_enter("mask_dense");
#endif

  /*----------------------------------------------------------------------*/
  /* remember some facts:
     PARTITION is different on every proc.
     FIELD is the same everywhere
     In this routine, the vector update is determined
     in size and allocated, the contents of the vector update
     are calculated
     */
  /*------------------------------------------- put total size of problem */
  dense->numeq_total = actfield->dis[0].numeq;
  /* count number of eqns on proc and build processor-global couplingdof
     matrix */
  dense_numeq(actfield,actpart,actsolv,actintra,&numeq);
  dense->numeq = numeq;
  /*---------------------------------------------- allocate vector update */
  amdef("update",&(dense->update),numeq,1,"IV");
  amzero(&(dense->update));
  /*--------------------------------put dofs in update in ascending order */
  dense_update(actfield,actpart,actsolv,actintra,dense);
  /*----------------------------------------------------- allocate matrix */
  amdef("A",&(dense->A),(dense->numeq_total),(dense->numeq_total),"DA");
  /*----------------------------------------------------------------------*/

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of mask_dense */




/*----------------------------------------------------------------------*
  |  count processor local and global number of equations    m.gee 5/01  |
 *----------------------------------------------------------------------*/
void  dense_numeq(
    FIELD         *actfield,
    PARTITION    *actpart,
    SOLVAR       *actsolv,
    INTRA        *actintra,
    INT          *numeq
    )

{
  INT       i,j,k,l;
  INT       counter;
  INT       dof;
  INT       iscoupled;
#ifdef PARALLEL
  INT      *sendbuff,*recvbuff, sendsize;
#endif
  INT      *tmp;
  INT       inter_proc;
  long int  min;
  INT       proc;
  INT       inprocs;
  INT       imyrank;
  NODE     *actnode;

  INT       no_coupling = 0;

#ifdef DEBUG
  dstrc_enter("dense_numeq");
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

  if (counter ==0)
    no_coupling = 1;
  else
    no_coupling = 0;

  amredef(&(actpart->pdis[0].coupledofs),counter,1,"IV");
  /*---------------------------------- delete the doubles in coupledofs */

  if (!no_coupling)
  {

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
      for (j=0; j<actpart->pdis[0].coupledofs.fdim; j++) actpart->pdis[0].coupledofs.a.ia[j][i]=0;


  } /* end of if(!no_coupling) */

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

  if (!no_coupling)
  {

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
    sendbuff = (INT*)CCACALLOC(sendsize,sizeof(INT));
    recvbuff = (INT*)CCACALLOC(sendsize,sizeof(INT));
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

  } /* end of if(!no_coupling) */

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

  if (!no_coupling)
  {

    if (inprocs > 1)
    {
      tmp = (INT*)CCACALLOC(inprocs,sizeof(INT));
      if (!tmp) dserror("Allocation of temporary memory failed");
      for (i=0; i<actpart->pdis[0].coupledofs.fdim; i++)/*------ loop coupled eqns */
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
      for (i=0; i<actpart->pdis[0].coupledofs.fdim; i++)/*------ loop coupled eqns */
      {
        /* ------Yes, I am slave owner of an inter_proc coupling equation */
        if (actpart->pdis[0].coupledofs.a.ia[i][imyrank+1]==1)
        {
          (*numeq) = (*numeq)-1;
        }
        /* master owner of equation do nothing, 'cause the equation has been
           counted anyway */
      }
    }

  } /* end of if(!no_coupling) */

  /*----------------------------------------------------------------------*/

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of dense_numeq */






/*----------------------------------------------------------------------*
  |  allocate update put dofs in update in ascending order   m.gee 5/01  |
 *----------------------------------------------------------------------*/
void  dense_update(
    FIELD         *actfield,
    PARTITION     *actpart,
    SOLVAR        *actsolv,
    INTRA         *actintra,
    DENSE         *dense
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
  dstrc_enter("dense_update");
#endif

  /*----------------------------------------------------------------------*/
  imyrank = actintra->intra_rank;
  inprocs = actintra->intra_nprocs;
  /*------------------ make a local copy of the array actpart->coupledofs */
  memset(&coupledofs, 0, sizeof(ARRAY));
  if (actpart->pdis[0].coupledofs.Typ != cca_XX)
    am_alloc_copy(&(actpart->pdis[0].coupledofs),&coupledofs);
  /*------------------------------------- loop the nodes on the partition */
  update = dense->update.a.iv;
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
  if (counter != dense->numeq) dserror("Number of dofs in DENSE-vector update wrong");
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
} /* end of dense_update */




/*----------------------------------------------------------------------*
  |  routine to assemble element array to global DENSE-matrix            |
  |  in parallel and sequentiell,taking care of coupling conditions      |
  |                                                                      |
  |                                                                      |
  |                                                         m.gee 9/01   |
 *----------------------------------------------------------------------*/
void  add_dense(
    struct _PARTITION     *actpart,
    struct _SOLVAR        *actsolv,
    struct _INTRA         *actintra,
    struct _ELEMENT       *actele,
    struct _DENSE         *dense1,
    struct _DENSE         *dense2,
    struct _ARRAY         *elearray1,
    struct _ARRAY         *elearray2
    )

{
  INT               i,j,counter;           /* some counter variables */
  INT               ii,jj;                 /* counter variables for system matrix */
  INT               nd;                    /* size of estif */
  INT               numeq_total;           /* total number of equations */
  INT               numeq;                 /* number of equations on this proc */
  INT               lm[MAXDOFPERELE];      /* location vector for this element */
#ifdef PARALLEL
  INT               owner[MAXDOFPERELE];   /* the owner of every dof */
#endif
  INT               myrank;                /* my intra-proc number */
  INT               nprocs;                /* my intra- number of processes */
  enum _SOLVER_TYP *solvertyp;
  DOUBLE           **estif;                 /* element matrix to be added to system matrix */
  DOUBLE           **emass;                 /* second element matrix to be assembled */
  DOUBLE           **A;                     /* the dense matrix */
  DOUBLE           **B;                     /* the second dense matrix */
  /*
     INT       **cdofs;
     INT         ncdofs;
     */

#ifdef DEBUG
  dstrc_enter("add_dense");
#endif

  /*------------------------------------- set some pointers and variables */
  myrank     = actintra->intra_rank;
  nprocs     = actintra->intra_nprocs;
  solvertyp  = &(actsolv->solvertyp);
  estif      = elearray1->a.da;
  if (dense2)
    emass      = elearray2->a.da;
  else
    emass      = NULL;
  nd         = actele->numnp * actele->node[0]->numdf;
  numeq_total= dense1->numeq_total;
  numeq      = dense1->numeq;
  A          = dense1->A.a.da;
  /* check for assembly of 2 matrices */
  if (dense2)
    B          = dense2->A.a.da;
  else
    B          = NULL;
  /*
     cdofs      = actpart->coupledofs.a.ia;
     ncdofs     = actpart->coupledofs.fdim;
     */
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
  /*
    NOTE:
    I don't have to care for coupling at all in this case, because
    system matrix is redundant on all procs, every proc adds his part
    (also slave and master owners of a coupled dof) and the system matrix
    is then allreduced. This makes things very comfortable for the moment.
  */
  /*======================================= loop over i (the element row) */
  if (*solvertyp==lapack_nonsym)/*-------------- unsymmetric lapack solve */
  {
    for (i=0; i<nd; i++)
    {
      ii = lm[i];
      /*---------------------------------- check for boundary condition */
      if (ii>=numeq_total) continue;
      /*--------------------------------- check for ownership of row ii */
#ifdef PARALLEL
      if (owner[i]!=myrank) continue;
#endif
      /*============================== loop over j (the element column) */
      /*                         This is the full unsymmetric version ! */
      for (j=0; j<nd; j++)
      {
        jj = lm[j];
        /*------------------------------- check for boundary condition */
        if (jj>=numeq_total) continue;
        /*--------------------------------------- add to system matrix */
        A[jj][ii] += estif[i][j];
        if (B)
          B[jj][ii] += emass[i][j];
      } /* end loop over j */
    }/* end loop over i */
  }
  else/*------------------------------------------ symmetric lapack solve */
  {
    for (i=0; i<nd; i++)
    {
      ii = lm[i];
      /*---------------------------------- check for boundary condition */
      if (ii>=numeq_total) continue;
      /*--------------------------------- check for ownership of row ii */
#ifdef PARALLEL
      if (owner[i]!=myrank) continue;
#endif
      /*============================== loop over j (the element column) */
      /*                         This is the full unsymmetric version ! */
      for (j=i; j<nd; j++)
      {
        jj = lm[j];
        /*------------------------------- check for boundary condition */
        if (jj>=numeq_total) continue;
        if (i==j)/*------------------------------- main diagonal entry */
        {
          A[jj][ii] += estif[i][j];
          if (B)
            B[jj][ii] += emass[i][j];
        }
        else/*------------------------------------- off diagonal entry */
        {
          A[ii][jj] += estif[i][j];
          A[jj][ii] += estif[i][j];
          if (B) {
            B[ii][jj] += emass[i][j];
            B[jj][ii] += emass[i][j];
          }
        }
      } /* end loop over j */
    }/* end loop over i */
  }
  /*----------------------------------------------------------------------*/

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of add_dense */




/*----------------------------------------------------------------------*
  |  make redundant dense matrix on all procs                 m.gee 11/01|
 *----------------------------------------------------------------------*/
void redundant_dense(
    PARTITION     *actpart,
    SOLVAR        *actsolv,
    INTRA         *actintra,
    DENSE         *dense1,
    DENSE         *dense2
    )

{

#ifdef PARALLEL
  INT      imyrank;
  INT      inprocs;

  ARRAY    recv_a;
  DOUBLE **recv;
#endif

#ifdef DEBUG
  dstrc_enter("redundant_dense");
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
  recv = amdef("recv_a",&recv_a,dense1->A.fdim,dense1->A.sdim,"DA");
  /*----------------------------------------------------------- Allreduce */
  MPI_Allreduce(dense1->A.a.da[0],
      recv[0],
      (dense1->A.fdim)*(dense1->A.sdim),
      MPI_DOUBLE,
      MPI_SUM,
      actintra->MPI_INTRA_COMM);
  /*----------------------------------------- copy reduced data back to a */
  amcopy(&recv_a,&(dense1->A));
  /* check for presence of second system matrix and reduce this one as well */
  if (dense2)
  {
    /*-------------------------------------------------------- Allreduce */
    MPI_Allreduce(dense2->A.a.da[0],
        recv[0],
        (dense2->A.fdim)*(dense2->A.sdim),
        MPI_DOUBLE,
        MPI_SUM,
        actintra->MPI_INTRA_COMM);

    /*-------------------------------------- copy reduced data back to a */
    amcopy(&recv_a,&(dense2->A));
  }
  /*----------------------------------------------------- delete recvbuff */
  amdel(&recv_a);
#endif /*---------------------------------------------- end of PARALLEL */
  /*----------------------------------------------------------------------*/

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of redundant_dense */


