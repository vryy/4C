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
/*!----------------------------------------------------------------------
\file
\brief contains functions to assemble oll matrices 

*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
/*----------------------------------------------------------------------*
 | global dense matrices for element routines             m.gee 9/01    |
 | (defined in global_calelm.c, so they are extern here)                |                
 *----------------------------------------------------------------------*/
extern struct _ARRAY estif_global;
extern struct _ARRAY emass_global;


/*! 
\addtogroup OLL 
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief fills sendbuffers isend und dsend 

<pre>                                                              mn 02/03
This function fills the sendbuffers isend und dsend for oll matrices.
</pre>
\param    ii            INT    (i)   
\param    jj            INT    (i)   
\param    i             INT    (i)   
\param    j             INT    (i)   
\param    ii_owner      INT    (i)   
\param  **isend         INT    (i)   
\param  **dsend         DOUBLE (i)   
\param  **estiff        DOUBLE (i)   
\param    numsend       INT    (i)   

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: --- 

*----------------------------------------------------------------------*/
void add_oll_sendbuff(
    INT ii,
    INT jj,
    INT i,
    INT j,
    INT ii_owner,
    INT **isend,
    DOUBLE **dsend,
    DOUBLE **estif,
    INT numsend)
{
  INT         k;
#ifdef DEBUG 
  dstrc_enter("add_oll_sendbuff");
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
} /* end of add_oll_sendbuff */




/*!----------------------------------------------------------------------
\brief checks coupling for the add_oll function

<pre>                                                              mn 02/03
This function counts the number of equations on this processor
</pre>
\param   ii             INT    (i)   
\param **cdofs          INT    (i)   
\param   ncdofs         INT    (i)   
\param  *iscoupled      INT    (i)   
\param  *isowner        INT    (i)   
\param   nprocs         INT    (i)   

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: --- 

 *----------------------------------------------------------------------*/
void add_oll_checkcouple(
    INT ii,
    INT **cdofs,
    INT ncdofs,
    INT *iscouple,
    INT *isowner,
    INT nprocs)
{
  INT         i,k;
#ifdef DEBUG 
  dstrc_enter("add_oll_checkcouple");
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
} /* end of add_oll_checkcouple */



/*!----------------------------------------------------------------------
\brief assembles global oll matrix

<pre>                                                              mn 02/03
This function assembles the element matrices to a global oll matrix, 
sequentially or parallel taking care of the couplig conditions
</pre>
\param *actpart        PARTITION (i)   the active partition
\param *actintra       INTRA     (i)   the active communicator
\param *actele         ELEMENT   (i)   the active element
\param *oll1           OLL       (o)   the first global oll matrix
\param *oll2           OLL       (o)   the second global oll matrix

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: --- 

*----------------------------------------------------------------------*/
void  add_oll(
    struct _PARTITION     *actpart,
    struct _INTRA         *actintra,
    struct _ELEMENT       *actele,
    struct _OLL           *oll1,
    struct _OLL           *oll2)
{
  INT         i,j,counter;       /* some counter variables */
  INT         istwo=0;
  INT         ii,jj;                 /* counter variables for system matrix */
  INT         ii_iscouple;           /* flag whether ii is a coupled dof */
  INT         ii_owner;              /* who is owner of dof ii -> procnumber */
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
  INT       **cdofs;                 /* list of coupled dofs and there owners, see init_assembly */
  INT         ncdofs;                /* total number of coupled dofs */
  INT       **isend1;                /* pointer to sendbuffer to communicate coupling conditions */
  DOUBLE    **dsend1;                /* pointer to sendbuffer to communicate coupling conditions */
  INT       **isend2;                /* pointer to sendbuffer to communicate coupling conditions */
  DOUBLE    **dsend2;                /* pointer to sendbuffer to communicate coupling conditions */
  INT         nsend;
#ifdef DEBUG 
  dstrc_enter("add_oll");
#endif
  /*----------------------------------------------------------------------*/
  /*----------------------- check whether to assemble one or two matrices */
  if (oll2) istwo=1;
  /*------------------------------------- set some pointers and variables */
  myrank     = actintra->intra_rank;
  nprocs     = actintra->intra_nprocs;
  estif      = estif_global.a.da;
  emass      = emass_global.a.da;
  nd         = actele->numnp * actele->node[0]->numdf;
  ndnd       = nd*nd;
  nnz        = oll1->nnz;
  numeq_total= oll1->numeq_total;
  numeq      = oll1->numeq;
  update     = oll1->update.a.iv;
  cdofs      = actpart->pdis[0].coupledofs.a.ia;
  ncdofs     = actpart->pdis[0].coupledofs.fdim;
  /*---------------------------------- put pointers to sendbuffers if any */
#ifdef PARALLEL 
  if (oll1->couple_i_send) 
  {
    isend1 = oll1->couple_i_send->a.ia;
    dsend1 = oll1->couple_d_send->a.da;
    nsend  = oll1->couple_i_send->fdim;
    if (istwo)
    {
      isend2 = oll2->couple_i_send->a.ia;
      dsend2 = oll2->couple_d_send->a.da;
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
      add_oll_checkcouple(ii,cdofs,ncdofs,&ii_iscouple,&ii_owner,nprocs);
    }
#endif
    /* -------------- either not a coupled dof or I am master owner */
    if (!ii_iscouple || ii_owner==myrank)
    {
      oll_addrow(oll1, ii, lm, estif[i], nd);
      if (istwo)
        oll_addrow(oll2, ii, lm, emass[i], nd);

    }
    /* ------------------------- a coupled dof and I am slave owner */
    else
    {
      /*=============================== loop over j (the element column) */
      for (j=0; j<nd; j++)
      {
        jj = lm[j];
        /* --------------------------------- check for boundary condition */
        if (jj>=numeq_total) continue;
        /* -------------------------------------------------------------- */
        add_oll_sendbuff(ii,jj,i,j,ii_owner,isend1,dsend1,estif,nsend);
        if (istwo)
          add_oll_sendbuff(ii,jj,i,j,ii_owner,isend2,dsend2,emass,nsend);
      } /* end loop over j */
    }
    /* ------------------------------------------------------------------- */ 
  }/* end loop over i */
  /*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  return;
} /* end of add_oll */




/*!----------------------------------------------------------------------
\brief exchanges the coupled dofs

<pre>                                                              mn 02/03
This function exchanges the coupled dofs for an oll matrix
</pre>
\param *actpart        PARTITION  (i)   the active partition
\param *actintra       INTRA  (i)   the active communicator
\param *oll            OLL    (i)   the oll matrix

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: --- 

 *----------------------------------------------------------------------*/
void exchange_coup_oll(
    PARTITION     *actpart,
    INTRA         *actintra,
    OLL           *oll)
{
  INT            i;
  INT            ii,jj;
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

#ifdef PARALLEL 
  MPI_Status    *irecv_status;
  MPI_Status    *drecv_status;

  MPI_Request   *isendrequest;
  MPI_Request   *dsendrequest;

  MPI_Comm      *ACTCOMM;
#endif

#ifdef DEBUG 
  dstrc_enter("exchange_coup_oll");
#endif
  /*----------------------------------------------------------------------*/
#ifdef PARALLEL 
  /*----------------------------------------------------------------------*/
  imyrank = actintra->intra_rank;
  inprocs = actintra->intra_nprocs;
  ACTCOMM = &(actintra->MPI_INTRA_COMM);
  /*---------------------------------------- set some pointers and values */
  numsend     = oll->numcoupsend;
  numrecv     = oll->numcouprecv;
  update      = oll->update.a.iv;
  numeq_total = oll->numeq_total;
  numeq       = oll->numeq;
  if (oll->couple_i_send) isend = oll->couple_i_send->a.ia;
  if (oll->couple_d_send) dsend = oll->couple_d_send->a.da;
  if (oll->couple_i_recv) irecv = oll->couple_i_recv->a.ia;
  if (oll->couple_d_recv) drecv = oll->couple_d_recv->a.da;
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
    for (jj=0; jj<numeq_total; jj++)
    {
      if(FABS(drecv[i][jj]) >= EPS12)
        oll_addval(oll, ii, jj, drecv[i][jj]);
    }
  }/*---------------------------------------------- end of receiving loop */
  /*-------------------------------------------- free allocated MPI-stuff */
  if (numrecv){CCAFREE(irecv_status);CCAFREE(drecv_status);}
  if (numsend){CCAFREE(isendrequest);CCAFREE(dsendrequest);}
  /*----------------------------------------------------------------------
    do a barrier, because this is the end of the assembly, the oll matrix
    is now ready for solve
    */ 
  MPI_Barrier(*ACTCOMM);
#endif /*---------------------------------------------- end of PARALLEL */ 
  /*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  return;
} /* end of exchange_coup_oll */

/*! @} (documentation module close)*/
