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
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
/*----------------------------------------------------------------------*
 | global dense matrices for element routines             m.gee 9/01    |
 | (defined in global_calelm.c, so they are extern here)                |                
 *----------------------------------------------------------------------*/
extern struct _ARRAY estif_global;
extern struct _ARRAY emass_global;
/*----------------------------------------------------------------------*
 |  routine to assemble element array to global rcptr-matrix            |
 |  in parallel,taking care of coupling conditions                      |
 |                                                                      |
 |                                                                      |
 |                                                         m.gee 1/02   |
 *----------------------------------------------------------------------*/
void  add_rc_ptr(struct _PARTITION     *actpart,
                struct _SOLVAR        *actsolv,
                struct _INTRA         *actintra,
                struct _ELEMENT       *actele,
                struct _RC_PTR        *rc_ptr1,
                struct _RC_PTR        *rc_ptr2)
{
#ifdef MUMPS_PACKAGE
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
dstrc_enter("add_rc_ptr");
#endif
/*----------------------------------------------------------------------*/
/*----------------------- check whether to assemble one or two matrices */
if (rc_ptr2) istwo=1;
/*------------------------------------- set some pointers and variables */
myrank     = actintra->intra_rank;
nprocs     = actintra->intra_nprocs;
estif      = estif_global.a.da;
emass      = emass_global.a.da;
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
      add_msr_checkcouple(ii,cdofs,ncdofs,&ii_iscouple,&ii_owner,nprocs);
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
#endif /* end of ifdef MUMPS_PACKAGE */
return;
} /* end of add_rc_ptr */

/*----------------------------------------------------------------------*
 |  fill sendbuffer isend and dsend                           m.gee 1/02|
 *----------------------------------------------------------------------*/
void add_rcptr_sendbuff(INT ii,INT jj,INT i,INT j,INT ii_owner,INT **isend,
                    DOUBLE **dsend,DOUBLE **estif, INT numsend)
{
#ifdef MUMPS_PACKAGE
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
#endif /* end of ifdef MUMPS_PACKAGE */
return;
} /* end of add_rcptr_sendbuff */


/*----------------------------------------------------------------------*
 |  exchange coupled dofs and add to row/column ptr matrix    m.gee 1/02|
 *----------------------------------------------------------------------*/
void exchange_coup_rc_ptr(
                         PARTITION     *actpart,
                         SOLVAR        *actsolv,
                         INTRA         *actintra,
                         RC_PTR        *rc_ptr
                        )
{
#ifdef MUMPS_PACKAGE
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
#endif /* end of ifdef MUMPS_PACKAGE */
return;
} /* end of exchange_coup_rc_ptr */
