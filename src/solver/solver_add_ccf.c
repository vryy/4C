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
/*----------------------------------------------------------------------*
 | global dense matrices for element routines             m.gee 9/01    |
 | (defined in global_calelm.c, so they are extern here)                |
 *----------------------------------------------------------------------*/
extern struct _ARRAY estif_global;
extern struct _ARRAY emass_global;


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
    struct _CCF           *ccf2)
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
  estif      = estif_global.a.da;
  emass      = emass_global.a.da;
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
  cdofs      = actpart->pdis[0].coupledofs.a.ia;
  ncdofs     = actpart->pdis[0].coupledofs.fdim;

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
    struct _CCF           *ccf2)
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
  estif      = estif_global.a.da;
  emass      = emass_global.a.da;
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
  cdofs      = actpart->pdis[0].coupledofs.a.ia;
  ncdofs     = actpart->pdis[0].coupledofs.fdim;

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



#endif /* end of ifdef UMFPACK */


