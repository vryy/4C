#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
/*----------------------------------------------------------------------*
 | global dense matrices for element routines             m.gee 10/01   |
 | (defined in global_calelm.c, so they are extern here)                |                
 *----------------------------------------------------------------------*/
extern struct _ARRAY estif_global;
extern struct _ARRAY emass_global;
/*----------------------------------------------------------------------*
 |  routine to assemble element array to global PARCSR-matrix           |
 |  in parallel and sequentiell,taking care of coupling conditions      |
 |                                                                      |
 |                                                                      |
 |                                                         m.gee 10/01  |
 *----------------------------------------------------------------------*/
void  add_parcsr(struct _PARTITION     *actpart,
                   struct _SOLVAR        *actsolv,
                   struct _INTRA         *actintra,
                   struct _ELEMENT       *actele,
                   struct _H_PARCSR      *parcsr)
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
/*----------------------------------------------------------------------*/
/*------------------------------------- set some pointers and variables */
myrank           = actintra->intra_rank;
nprocs           = actintra->intra_nprocs;

estif            = estif_global.a.da;
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
#ifdef HYPRE_PACKAGE
      err=HYPRE_IJMatrixAddToValues(
                                    parcsr->ij_matrix,
                                    nrows,
                                    ncols,
                                    rows,
                                    cols,
                                    values
                                   );
#endif
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
void add_parcsr_sendbuff(INT ii,INT jj,INT i,INT j,INT ii_owner,INT **isend,
                    DOUBLE **dsend,DOUBLE **estif, INT numsend)
{
INT         k,l;
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
void add_parcsr_checkcouple(INT ii,INT **cdofs,INT ncdofs,INT *iscouple,INT *isowner, INT nprocs)
{
INT         i,j,k;
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

#ifdef PARALLEL 
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
#ifdef HYPRE_PACKAGE
   err=HYPRE_IJMatrixAddToValues(
                                 parcsr->ij_matrix,
                                 nrows,
                                 ncols,
                                 rows,
                                 cols,
                                 values
                                );
#endif
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
