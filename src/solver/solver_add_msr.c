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


/*----------------------------------------------------------------------*
  | global dense matrices for element routines             m.gee 9/01  |
  | (defined in global_calelm.c, so they are extern here)              |                
 *----------------------------------------------------------------------*/
extern struct _ARRAY estif_global;
extern struct _ARRAY emass_global;

/*----------------------------------------------------------------------*
  |  routine to assemble element array to global DMSR-matrix           |
  |  in parallel and sequentiell,taking care of coupling conditions    |
  |                                                                    |
  |                                                                    |
  |                                                         m.gee 9/01 |
 *----------------------------------------------------------------------*/
void  add_msr(
    struct _PARTITION     *actpart,
    struct _SOLVAR        *actsolv,
    struct _INTRA         *actintra,
    struct _ELEMENT       *actele,
    struct _AZ_ARRAY_MSR  *msr1,
    struct _AZ_ARRAY_MSR  *msr2)
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
  dstrc_enter("add_msr");
#endif

  /* set some pointers and variables */
  myrank     = actintra->intra_rank;
  nprocs     = actintra->intra_nprocs;
  estif      = estif_global.a.da;
  if (msr2) emass = emass_global.a.da;
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
} /* end of add_msr */



/*----------------------------------------------------------------------*
 |  checks coupling for the add_msr routine                   m.gee 9/01|
 *----------------------------------------------------------------------*/
void add_msr_checkcouple(INT ii,INT **cdofs,INT ncdofs,INT *iscouple,
                           INT *isowner, INT nprocs)
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
 |  fill sendbuffer isend and dsend                           m.gee 9/01|
 *----------------------------------------------------------------------*/
void add_msr_sendbuff(INT ii,INT jj,INT i,INT j,INT ii_owner,INT **isend,
                    DOUBLE **dsend,DOUBLE **estif, INT numsend)
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
 |  exchange coupled dofs and add to dmsr matrix              m.gee 9/01|
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
