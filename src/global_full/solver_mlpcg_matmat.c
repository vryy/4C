/*!---------------------------------------------------------------------
\file
\brief contains the multilevel preconditioner for shells

---------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/prototypes_mlpcg.h"
#include "../shell8/shell8.h"
INT cmp_int(const void *a, const void *b );
DOUBLE cmp_double(const void *a, const void *b );
/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h                                                  
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | struct _GENPROB       genprob; defined in global_control.c           |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*! 
\addtogroup MLPCG 
*//*! @{ (documentation module open)*/
/*!----------------------------------------------------------------------
\brief the multilevel preconditioner main structure

<pre>                                                         m.gee 09/02    
defined in solver_mlpcg.c
</pre>

*----------------------------------------------------------------------*/
extern struct _MLPRECOND mlprecond;




/*!---------------------------------------------------------------------
\brief make the parallel matrix-matrix-matrix product                                         

<pre>                                                        m.gee 10/02 
make the parallel matrix-matrix-matrix product
- outcsr = P(transposed) * incsr * P
- outcsr gets the proc partitioning from P
- P      is a closed matrix with no 'ellbow-space'
- incsr  is a closed matrix with no 'ellbow-space'
- outcsr is an open matrix with 'ellbow-space'
- work   is an open matrix with 'ellbow-space'
- incsr is symmetric !
- P is not symmetric !

After the product 
work = P(transposed) * incsr
is done, the matrix work is closed by call to mlpcg_csr_close
Then the product 
outcsr = work * P 
is done, afterwards the matrix outcsr is closed.

THIS IS NOT AN ALL PUPOSE ROUTINE, IT MAKES USE OF SPECIAL KNOWLEDGE
ABOUT THE PROLONGATOR!
</pre>
\param P           DBCSR*    (i) the matrix P of the product                   
\param incsr       DBCSR*    (i) the matrix incsr of the product                   
\param outcsr      DBCSR*    (o) the matrix outcsr of the product                   
\param work        DBCSR*    (o) the working matrix                   
\param actintra    INTRA*    (i) the communicator                 
\warning the matrix incsr has to be symmetric!!
\return void                                               
------------------------------------------------------------------------*/
void mlpcg_precond_PtKP(DBCSR *P, DBCSR *incsr, DBCSR *outcsr, DBCSR *work,
                        AGG *agg, INT nagg, INTRA *actintra)
{
INT        i,j,k,m,n,counter,foundit,ilength,tag;
INT        nproc,myrank;
DOUBLE     sum;
INT        index;
INT        owner;
INT        ownertmp[MAXPROC][2];
INT        nsend,nrecv;
INT        numeq;

INT        sendtos[MAXPROC][MAXPROC],sendtor[MAXPROC][MAXPROC];
ARRAY      isbuff[MAXPROC];
ARRAY      dsbuff[MAXPROC];
ARRAY      irbuff;
ARRAY      drbuff;
ARRAY      rupdate; INT *updateP;
ARRAY      ria;     INT *iaP;
ARRAY      rja;     INT *jaP;
ARRAY      ra;   DOUBLE *aP;
INT        rcounter[MAXPROC];
INT       *isend,*irecv;
DOUBLE    *dsend,*drecv;
 
INT        actrow,actcol;
INT        colstart,colend;

DOUBLE     block[500][500];
INT        rindex[500];
INT        cindex[500];

DOUBLE     *col,*colP;
INT        *rcol,*rcolP;
INT        scol = 5000, scolP = 5000;
INT        nrow,ncol;

 /* these are tricky column pointers to the incsr matrix */
INT     ***icol;    
INT       *colsize;
DOUBLE  ***dcol;

INT        shift;
INT        bins[1000];

INT        numeq_cscP,*update_cscP,*ia_cscP,*ja_cscP;
DOUBLE    *a_cscP;
INT        numeq_in,*update_in,*ia_in,*ja_in;
DOUBLE    *a_in;
#if 0
DOUBLE     t1,t2;
#endif
INT        min1,max1,min2,max2;

#ifdef PARALLEL 
MPI_Status   status;  
MPI_Request  *request;   
#endif

#if 0 /* for development and debugging */
ARRAY      sdense_in, rdense_in, sdense_P, rdense_P, dense_work, dense_out;
ARRAY      sw,rw;
ARRAY      so,ro;
#endif

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_PtKP");
#endif
/*----------------------------------------------------------------------*/
/*--------------------------- allocate row and column extractor buffers */
rcol  = (INT*)   CCAMALLOC(scol*sizeof(INT));
col   = (DOUBLE*)CCAMALLOC(scol*sizeof(DOUBLE));
rcolP = (INT*)   CCAMALLOC(scolP*sizeof(INT));
colP  = (DOUBLE*)CCAMALLOC(scolP*sizeof(DOUBLE));
/*----------------------------------------------------------------------*/
myrank = actintra->intra_rank;
nproc  = actintra->intra_nprocs;
/*----------------------------------------------------------------------*/
#if 0 /*------------------- debugging, make dense copies of P and incsr */
amdef("tmp",&sdense_in,incsr->numeq_total,incsr->numeq_total,"DA");
amdef("tmp",&rdense_in,incsr->numeq_total,incsr->numeq_total,"DA");
amzero(&sdense_in);
amzero(&rdense_in);
amdef("tmp",&sdense_P,incsr->numeq_total,work->numeq_total,"DA");
amdef("tmp",&rdense_P,incsr->numeq_total,work->numeq_total,"DA");
amzero(&sdense_P);
amzero(&rdense_P);
amdef("tmp",&dense_work,work->numeq_total,incsr->numeq_total,"DA");
amzero(&dense_work);
amdef("tmp",&dense_out,work->numeq_total,work->numeq_total,"DA");
amzero(&dense_out);
/* fill sdense_in */
for (i=0; i<incsr->numeq; i++)
{
   actrow = incsr->update.a.iv[i];
   colstart = incsr->ia.a.iv[i];
   colend   = incsr->ia.a.iv[i+1];
   for (j=colstart; j<colend; j++)
   {
      actcol = incsr->ja.a.iv[j];
      sdense_in.a.da[actrow][actcol] = incsr->a.a.dv[j];
   }
}
#ifdef PARALLEL 
MPI_Allreduce(sdense_in.a.da[0],rdense_in.a.da[0],sdense_in.fdim*sdense_in.sdim,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
#else
for (i=0; i<sdense_in.fdim*sdense_in.sdim; i++) 
   rdense_in.a.da[0][i] = sdense_in.a.da[0][i];
#endif
/* check symmetry */
for (i=0; i<sdense_in.fdim; i++)
{
   for (j=i+1; j<sdense_in.sdim; j++)
     if (FABS(rdense_in.a.da[i][j]-rdense_in.a.da[j][i])>EPS12)
      fprintf(allfiles.out_err,"K not symmetric in ij %d %d %20.10f ji %d %d %20.10f diff %30.20f\n",i,j,rdense_in.a.da[i][j],j,i,rdense_in.a.da[j][i],FABS(rdense_in.a.da[i][j]-rdense_in.a.da[j][i])); 
}
/* print to err */
fprintf(allfiles.out_err,"in-Matrix\n");
for (i=0; i<rdense_in.fdim; i++)
{
   for (j=0; j<rdense_in.sdim; j++)
   {
    /*  fprintf(allfiles.out_err,"%8.2f ",rdense_in.a.da[i][j]);*/
        if (FABS(rdense_in.a.da[i][j])>EPS14)
        fprintf(allfiles.out_err,"X ");
        else
        fprintf(allfiles.out_err,"* ");
   }
   fprintf(allfiles.out_err," \n");
}
/* fill sdense_P */
for (i=0; i<P->numeq; i++)
{
   actrow = P->update.a.iv[i];
   colstart = P->ia.a.iv[i];
   colend   = P->ia.a.iv[i+1];
   for (j=colstart; j<colend; j++)
   {
      actcol = P->ja.a.iv[j];
      sdense_P.a.da[actrow][actcol] = P->a.a.dv[j];
   }
}
#ifdef PARALLEL 
MPI_Allreduce(sdense_P.a.da[0],rdense_P.a.da[0],sdense_P.fdim*sdense_P.sdim,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
#else
for (i=0; i<sdense_P.fdim*sdense_P.sdim; i++)
   rdense_P.a.da[0][i] = sdense_P.a.da[0][i];
#endif
/* print to err */
fprintf(allfiles.out_err,"P-Matrix\n");
for (i=0; i<rdense_P.fdim; i++)
{
   for (j=0; j<rdense_P.sdim; j++)
   {
     fprintf(allfiles.out_err,"%15.10E ",rdense_P.a.da[i][j]);
    /*     if (FABS(rdense_P.a.da[i][j])>EPS14)
        fprintf(allfiles.out_err,"X ");
        else
        fprintf(allfiles.out_err,"* ");*/
   }
   fprintf(allfiles.out_err," \n");
}
/*-------------------- multiply work = Ptransposed * in */
for (i=0; i<dense_work.fdim; i++)/* row of work */
{
   for (j=0; j<dense_work.sdim; j++)/* cols of work */
   {
      sum = 0.0;
      for (k=0; k<rdense_P.fdim; k++)
         sum += rdense_P.a.da[k][i] * rdense_in.a.da[k][j];
      dense_work.a.da[i][j] = sum;
   }
}
/* print to err */
fprintf(allfiles.out_err,"work-Matrix\n");
for (i=0; i<dense_work.fdim; i++)
{
   for (j=0; j<dense_work.sdim; j++)
   {
      fprintf(allfiles.out_err,"%25.8E ",dense_work.a.da[i][j]);
     /*   if (FABS(dense_work.a.da[i][j])>EPS14)
        fprintf(allfiles.out_err,"X ");
        else
        fprintf(allfiles.out_err,"* ");*/
   }
   fprintf(allfiles.out_err," \n");
}
/*------------------------------- multiply out = work * P */
for (i=0; i<dense_out.fdim; i++)/* rows of out */
{
   for (j=0; j<dense_out.sdim; j++)/* columns of work */
   {
      sum = 0.0;
      for (k=0; k<dense_work.sdim; k++)
         sum += dense_work.a.da[i][k] * rdense_P.a.da[k][j];
      dense_out.a.da[i][j] = sum;
   }
}
/* print to err */
fprintf(allfiles.out_err,"out-Matrix\n");
/* check symmetry */
for (i=0; i<dense_out.fdim; i++)
{
   for (j=i+1; j<dense_out.sdim; j++)
     if (FABS(dense_out.a.da[i][j]-dense_out.a.da[j][i])>EPS12)
      fprintf(allfiles.out_err,"K not symmetric in ij %d %d %20.10f ji %d %d %20.10f\n",i,j,dense_out.a.da[i][j],j,i,dense_out.a.da[j][i]); 
}
for (i=0; i<dense_out.fdim; i++)
{
   for (j=0; j<dense_out.sdim; j++)
   {
    /*  fprintf(allfiles.out_err,"%8.2f ",dense_out.a.da[i][j]);*/
        if (FABS(dense_out.a.da[i][j])>EPS14)
        fprintf(allfiles.out_err,"X ");
        else
        fprintf(allfiles.out_err,"* ");
   }
   fprintf(allfiles.out_err," \n");
}
fflush(allfiles.out_err);
#endif
/*------------------------------------------------------end of debugging*/



/*- do a distributed compressed sparse column copy of the prolongator P */
mlpcg_csr_csrtocsc(P,actintra);
/*------------------ construct the column pointers for the incsr matrix */
/*
t1 = ds_cputime();
*/

colsize = (INT*)CCACALLOC(incsr->numeq_total,sizeof(INT));
icol    = (INT***)CCAMALLOC(incsr->numeq_total*sizeof(INT**));
dcol    = (DOUBLE***)CCAMALLOC(incsr->numeq_total*sizeof(DOUBLE**));
mlpcg_extractcollocal_init(incsr,colsize,icol,dcol);
/*
t2 = ds_cputime();
if (myrank==0) printf("collocal_init            : %20.10f\n",t2-t1);
*/
/*--------------------------------------------------- set some pointers */
numeq_cscP  = P->csc->numeq;
update_cscP = P->csc->update.a.iv;
ia_cscP     = P->csc->ia.a.iv;
ja_cscP     = P->csc->ja.a.iv;
a_cscP      = P->csc->a.a.dv;

numeq_in    = incsr->numeq;
update_in   = incsr->update.a.iv;
ia_in       = incsr->ia.a.iv;
ja_in       = incsr->ja.a.iv;
a_in        = incsr->a.a.dv;
/*---------------------------------------- create the ownership of work */
for (n=0; n<nproc; n++) 
{ 
   work->owner[n][0] = 0; 
   work->owner[n][1] = 0;
   rcounter[n]       = 0;
}
for (n=0; n<nproc; n++)
for (m=0; m<nproc; m++) sendtos[n][m] = 0;
work->owner[myrank][0] = work->update.a.iv[0];
work->owner[myrank][1] = work->update.a.iv[work->numeq-1];
#ifdef PARALLEL 
#if 1
for (n=0; n<nproc; n++) 
{
   ownertmp[n][0] = 0;
   ownertmp[n][1] = 0;
}
ownertmp[myrank][0] = work->update.a.iv[0];
ownertmp[myrank][1] = work->update.a.iv[work->numeq-1];
MPI_Allreduce(&(ownertmp[0][0]),&(work->owner[0][0]),MAXPROC*2,MPI_INT,MPI_SUM,actintra->MPI_INTRA_COMM);
#endif
#if 0
for (n=0; n<nproc; n++)
   MPI_Bcast(work->owner[n],2,MPI_INT,n,actintra->MPI_INTRA_COMM);
#endif
/*======================================================================*/
/* do interproc computation work = Ptransposed * incsr   part I         */
/*======================================================================*/
/*
t1 = ds_cputime();
*/
/*--------- loop columns in P which as rows in work do not belong to me */
/*-------------------------- make a list of the owners of these columns */
/*---------------------- count how many rows I have to send to somebody */
nsend=0;
for (i=0; i<numeq_cscP; i++)
{
   actcol = update_cscP[i];
   owner  = mlpcg_getowner(actcol,work->owner,nproc);
   if (owner==myrank) continue;
   sendtos[myrank][owner]++;
   nsend++;
}
MPI_Allreduce(&(sendtos[0][0]),&(sendtor[0][0]),MAXPROC*MAXPROC,MPI_INT,MPI_SUM,actintra->MPI_INTRA_COMM);
/*--------------------------------------------------- count my receives */
nrecv=0;
for (n=0; n<nproc; n++) 
{
   if (n==myrank) continue;
   nrecv += sendtor[n][myrank];
}
/*------------------------------- loop sendtor and allocate sendbuffers */
for (n=0; n<nproc; n++)
{
   if (n==myrank) continue;
   if (sendtor[myrank][n] != 0) /* I have to send to n */
   {
      amdef("isbuff",&(isbuff[n]),sendtor[myrank][n],50,"IA");
      amzero(&(isbuff[n]));
      amdef("dsbuff",&(dsbuff[n]),sendtor[myrank][n],50,"DA");
   }
}
/*--------- loop columns in P which as rows in work do not belong to me */
/*------------------- make remote multiplication and put to sendbuffers */
for (i=0; i<numeq_cscP; i++)/* loop remote columns in P */
{
   actrow = update_cscP[i];
   owner  = mlpcg_getowner(actrow,work->owner,nproc);
   if (owner==myrank) 
      continue;
   dsassert(rcounter[owner]<isbuff[owner].fdim,"sendbuffer not enough rows");
   isend    = &(isbuff[owner].a.ia[rcounter[owner]][0]);   
   dsend    = &(dsbuff[owner].a.da[rcounter[owner]][0]);
   dsend[0] = (DOUBLE)actrow;
   rcounter[owner]++;
   mlpcg_extractcolcsc(actrow,numeq_cscP,update_cscP,ia_cscP,ja_cscP,a_cscP,
                       &colP,&rcolP,&scolP,&ncol);
   if (ncol==0) 
      continue;
   min2 = rcolP[0];
   max2 = rcolP[ncol-1];
   for (j=0; j<incsr->numeq_total; j++) /* loop all columns in incsr */
   {
      actcol = j;
      mlpcg_extractcollocal_fast(incsr,actcol,&col,&rcol,&scol,&nrow,colsize,icol,dcol);
      if (nrow==0)
         continue;
      init_quick_find(rcol,nrow,&shift,bins);
      min1    = rcol[0];
      max1    = rcol[nrow-1];
      if (min2 > max1 || min1 > max2 || IMAX(min1,min2) > IMIN(max2,max2))
         continue;
      foundit = 0;
      sum     = 0.0;
      for (k=0; k<ncol; k++)
      {
         index = quick_find(rcolP[k],rcol,nrow,shift,bins);
         if (index==-1) continue;
         foundit=1;
         sum += colP[k]*col[index];
      }
      if (foundit)
      {
         isend[isend[0]+1] = actcol;
         dsend[isend[0]+1] = sum;
         isend[0]++;
         if (isend[0]+1 >= isbuff[owner].sdim)
         {
            amredef(&(isbuff[owner]),isbuff[owner].fdim,isbuff[owner].sdim+500,"IA");
            amredef(&(dsbuff[owner]),dsbuff[owner].fdim,dsbuff[owner].sdim+500,"DA");
            isend = &(isbuff[owner].a.ia[rcounter[owner]-1][0]);   
            dsend = &(dsbuff[owner].a.da[rcounter[owner]-1][0]);   
         }
      }
   }
}
/*----------------------------------------------- check number of lines */
for (n=0; n<nproc; n++)
{
   if (myrank==n) continue;
   if (sendtor[myrank][n] != rcounter[n]) dserror("Number of lines wrong");
}
/*--------------------------------------------------- allocate requests */
request = (MPI_Request*)CCACALLOC(2*nsend,sizeof(MPI_Request));
/*-------------------------------------------------- now make the sends */
counter=0;
if (2*nsend>=10000) dserror("Processors unique tag range too small");
tag = myrank*10000;
for (n=0; n<nproc; n++)
{
   if (myrank==n) continue;
   for (k=0; k<rcounter[n]; k++)
   {
      isend = &(isbuff[n].a.ia[k][0]);
      dsend = &(dsbuff[n].a.da[k][0]);
      MPI_Isend(isend,isend[0]+1,MPI_INT,n,tag,actintra->MPI_INTRA_COMM,&(request[counter]));
      counter++;tag++;
      MPI_Isend(dsend,isend[0]+1,MPI_DOUBLE,n,tag,actintra->MPI_INTRA_COMM,&(request[counter]));
      counter++;tag++;
   }
}
dsassert(counter==2*nsend,"number of sends wrong");
/*----------------------------------------------------------------------*/
/*
t2 = ds_cputime();
if (myrank==0) printf("work = Pt*incsr inter I  : %20.10f\n",t2-t1);
*/
#endif
/*======================================================================*/
/*                   do local computation of work = Ptransposed * incsr */
/*======================================================================*/
/*
t1 = ds_cputime();
*/
for (i=0; i<incsr->numeq_total; i++)/* loop all columns of work */
{
   actcol = i;
   /* extract the column from incsr */
   mlpcg_extractcollocal_fast(incsr,actcol,&col,&rcol,&scol,&nrow,colsize,icol,dcol);
   if (nrow==0) continue;
   min1 = rcol[0];
   max1 = rcol[nrow-1];
   init_quick_find(rcol,nrow,&shift,bins);
   for (j=0; j<work->numeq; j++) /* loop the local rows */
   {
      actrow = work->update.a.iv[j];
      /* extract the column from P */
      mlpcg_extractcolcsc(actrow,numeq_cscP,update_cscP,ia_cscP,ja_cscP,a_cscP,
                          &colP,&rcolP,&scolP,&ncol);
      if (ncol==0) continue;
      min2    = rcolP[0];
      max2    = rcolP[ncol-1];
      if (min2 > max1 || min1 > max2 || IMAX(min1,min2)>IMIN(max1,max2))
         continue;
      foundit = 0;
      sum     = 0.0;
      for (k=0; k<ncol; k++)
      {
         if (rcolP[k] < min1 || rcolP[k] > max1)
            continue;
         index = quick_find(rcolP[k],rcol,nrow,shift,bins);
         if (index==-1) continue;
         foundit=1;
         sum += colP[k]*col[index];
      }
      if (foundit)
         mlpcg_csr_addentry(work,sum,actrow,actcol,actintra);
   }
}
/*
t2 = ds_cputime();
if (myrank==0) printf("work = Pt*incsr local    : %20.10f\n",t2-t1);
fflush(stdout);
*/
/*======================================================================*/
/* do interproc computation work = Ptransposed * incsr   part II        */
/*======================================================================*/
#ifdef PARALLEL 
/*
t1 = ds_cputime();
*/
irecv = amdef("irbuff",&irbuff,1000,1,"IV");
drecv = amdef("drbuff",&drbuff,1000,1,"DV");
for (n=0; n<nproc; n++)
{
   if (n==myrank)
      continue;
   for (k=0; k<sendtor[n][myrank]; k++)
   {
   MPI_Probe(n,MPI_ANY_TAG,actintra->MPI_INTRA_COMM,&status);
   /*--------------------------------------------------- get the sender */
   if (n != status.MPI_SOURCE)
     dserror("Sender not as expected\n");
   /*-------------------------------------------------------- check tag */
   tag = status.MPI_TAG;
   /*------------------------------------------------- get message size */
   MPI_Get_count(&status,MPI_INT,&ilength);
   /*------------------------------------ check size of receive buffers */
   if (ilength > irbuff.fdim)
   {
      amdel(&irbuff);
      amdel(&drbuff);
      irecv = amdef("irbuff",&irbuff,ilength+10,1,"IV");
      drecv = amdef("drbuff",&drbuff,ilength+10,1,"DV");
   }
   /*-------------------------------------- receive the integer message */
   MPI_Recv(irecv,ilength,MPI_INT,n,tag,actintra->MPI_INTRA_COMM,&status);
   /*---------------------------------- receive matching DOUBLE-message */
   MPI_Recv(drecv,ilength,MPI_DOUBLE,n,tag+1,actintra->MPI_INTRA_COMM,&status);
   /*--------------------------- decrease number of unreceived messages */
   nrecv--;
   /*---------------------------------------------- check for right row */
   actrow = (INT)drecv[0];
   owner = mlpcg_getowner(actrow,work->owner,nproc);
   if (owner != myrank) 
   {
      printf("myrank %d message from %d row %f tag %d myowner from %d to %d!!!!!!!!!!!!!!!!!!!!!\n",
             myrank,n,drecv[0],tag,work->owner[myrank][0],work->owner[myrank][1]);
      dserror("Received message does not fit my piece of matrix");
   }
   /*------------------------------------ put row to my piece of matrix */
   ncol = irecv[0];
   mlpcg_csr_addrow(work,actrow,&(drecv[1]),&(irecv[1]),ncol,actintra);
   }
} /* end of for (n=0; n<nproc; n++) */
if (nrecv != 0) dserror("Numer of receives wrong");
/*------------------------------------------ wait for all sent messages */
for (i=0; i<2*nsend; i++) 
   MPI_Wait(&(request[i]),&status);
/*======================================================================*/
/*                                                       tidy up        */
/*======================================================================*/
for (n=0; n<nproc; n++)
{
   if (n==myrank) continue;
   if (sendtor[myrank][n] != 0) /* I had to send to n */
   {
      amdel(&(isbuff[n]));
      amdel(&(dsbuff[n]));
   }
}
amdel(&irbuff);
amdel(&drbuff);
CCAFREE(request);
/*
t2 = ds_cputime();
if (myrank==0) printf("work = Pt*incsr inter II : %20.10f\n",t2-t1);
*/
#endif
/*======================================================================*/
/*                                         close the work matrix        */
/*======================================================================*/
/*
t1 = ds_cputime();
*/
mlpcg_csr_close(work);
/*
t2 = ds_cputime();
if (myrank==0) printf("work close               : %20.10f\n",t2-t1);
*/
/*------------------------------debugging */
#if 0 /* make dense printout from work */
amdef("sw",&sw,work->numeq_total,incsr->numeq_total,"DA");
amzero(&sw);
amdef("rw",&rw,work->numeq_total,incsr->numeq_total,"DA");
for (i=0; i<work->numeq; i++)
{
   actrow   = work->update.a.iv[i];
   colstart = work->ia.a.iv[i]; 
   colend   = work->ia.a.iv[i+1];
   for (j=colstart; j<colend; j++)
      sw.a.da[actrow][work->ja.a.iv[j]] =  work->a.a.dv[j];
}
#ifdef PARALLEL 
MPI_Allreduce(sw.a.da[0],rw.a.da[0],sw.fdim*sw.sdim,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
#else
for (i=0; i<sw.fdim*sw.sdim; i++)
  rw.a.da[0][i] = sw.a.da[0][i];
#endif
fprintf(allfiles.out_err,"work-Matrix from sparse multiplication\n");
for (i=0; i<rw.fdim; i++)
{
   for (j=0; j<rw.sdim; j++)
   {
      fprintf(allfiles.out_err,"%25.10E ",dense_work.a.da[i][j]);
    /*    if (FABS(rw.a.da[i][j])>EPS14)
        fprintf(allfiles.out_err,"X ");
        else
        fprintf(allfiles.out_err,"* ");*/
   }
   fprintf(allfiles.out_err," \n");
}
fprintf(allfiles.out_err,"work-Matrix from sparse multiplication minus work matrix from seq. dense mult.\n");
for (i=0; i<rw.fdim; i++)
{
   for (j=0; j<rw.sdim; j++)
   {
      fprintf(allfiles.out_err,"%13.4E ",rw.a.da[i][j]-dense_work.a.da[i][j]);
    /*    if (FABS(rw.a.da[i][j]-dense_work.a.da[i][j])>EPS12)
        fprintf(allfiles.out_err,"X ");
        else
        fprintf(allfiles.out_err,"* ");*/
   }
   fprintf(allfiles.out_err," \n");
}
fflush(allfiles.out_err);
#endif
/*-------------------------------debugging */
/*======================================================================*/
/*             do interproc computation out = work * P   part I         */
/*======================================================================*/
#ifdef PARALLEL 
/*
t1 = ds_cputime();
*/
/*------------------------------ send my prolongator to all other procs */
nsend = nproc-1;
nrecv = nproc-1;
/*--------------------------------------------------- allocate requests */
request = (MPI_Request*)CCACALLOC(4*nsend,sizeof(MPI_Request));
counter=0;
for (n=0; n<nproc; n++)
{
   if (n==myrank) continue;
   MPI_Isend(update_cscP,numeq_cscP,MPI_INT,n,counter,actintra->MPI_INTRA_COMM,&(request[counter]));
   counter++;
   MPI_Isend(ia_cscP,numeq_cscP+1,MPI_INT,n,counter,actintra->MPI_INTRA_COMM,&(request[counter]));
   counter++;
   MPI_Isend(ja_cscP,ia_cscP[numeq_cscP],MPI_INT,n,counter,actintra->MPI_INTRA_COMM,&(request[counter]));
   counter++;
   MPI_Isend(a_cscP,ia_cscP[numeq_cscP],MPI_DOUBLE,n,counter,actintra->MPI_INTRA_COMM,&(request[counter]));
   counter++;
}
dsassert(counter==4*nsend,"Number of sends wrong");
/*
t2 = ds_cputime();
if (myrank==0) printf("out  = work * P inter I  : %20.10f\n",t2-t1);
*/
#endif
/*======================================================================*/
/*                               do local computation of out = work * P */
/*======================================================================*/
/*
t1 = ds_cputime();
*/
#ifdef PARALLEL 
MPI_Barrier(actintra->MPI_INTRA_COMM);
#endif
for (i=0; i<work->numeq; i++)/* loop all my rows of work */
{
   actrow = work->update.a.iv[i];
   mlpcg_extractrowcsr(actrow,work->numeq,work->update.a.iv,work->ia.a.iv,
                       work->ja.a.iv,work->a.a.dv,&col,&rcol,&scol,&nrow);
   if (nrow==0) continue;
   init_quick_find(rcol,nrow,&shift,bins);
   min1 = rcol[0];
   max1 = rcol[nrow-1];
   for (j=0; j<numeq_cscP; j++)/* loop all columns in P */
   {
      actcol = update_cscP[j];
      mlpcg_extractcolcsc(actcol,numeq_cscP,update_cscP,ia_cscP,ja_cscP,a_cscP,
                          &colP,&rcolP,&scolP,&ncol);
      if (ncol==0) continue; 
      min2 = rcolP[0];
      max2 = rcolP[ncol-1];
      if (min2 > max1 || min1 > max2 || IMAX(min1,min2)>IMIN(max1,max2))
         continue;
      /* make the multiplication */
      foundit = 0;
      sum     = 0.0;
      for (k=0; k<ncol; k++)
      {
         if (rcolP[k] < min1 || rcolP[k] > max1)
            continue;
         index = quick_find(rcolP[k],rcol,nrow,shift,bins);
         if (index==-1) continue;
         foundit=1;
         sum += colP[k]*col[index];
      }
      if (foundit)
         mlpcg_csr_addentry(outcsr,sum,actrow,actcol,actintra);
   }
}
/*
t2 = ds_cputime();
if (myrank==0) printf("out  = work * P local    : %20.10f\n",t2-t1);
*/
#ifdef PARALLEL 
MPI_Barrier(actintra->MPI_INTRA_COMM);
#endif
/*======================================================================*/
/*             do interproc computation out = work * P   part II        */
/*======================================================================*/
/*-------------------------------------------- allocate receive buffers */
#ifdef PARALLEL 
/*
t1 = ds_cputime();
*/
amdef("upd",&rupdate,1,1,"IV");
amdef("ia" ,&ria    ,1,1,"IV");
amdef("ja" ,&rja    ,1,1,"IV");
amdef("a"  ,&ra     ,1,1,"DV");
while(nrecv!=0)
{
   MPI_Probe(MPI_ANY_SOURCE,MPI_ANY_TAG,actintra->MPI_INTRA_COMM,&status);
   /*--------------------------------------------------- get the sender */
   n = status.MPI_SOURCE;
   /*-------------------------------------------------------- check tag */
   tag = status.MPI_TAG;
   /*------------------------------------------------- get message size */
   MPI_Get_count(&status,MPI_INT,&numeq);
   if (numeq>rupdate.fdim) 
   {
      amdel(&rupdate);
      amdel(&ria);
      updateP = amdef("upd",&rupdate,numeq  ,1,"IV");
      iaP     = amdef("upd",&ria    ,numeq+1,1,"IV");
   }
   /*--------------------------------------------------- receive update */  
   MPI_Recv(updateP,numeq,MPI_INT,n,tag,actintra->MPI_INTRA_COMM,&status);
   /*------------------------------------------------------- receive ia */
   MPI_Recv(iaP    ,numeq+1,MPI_INT,n,tag+1,actintra->MPI_INTRA_COMM,&status);
   /*------------------------------------------- get length of ja and a */
   if (iaP[numeq] > rja.fdim)
   {
      amdel(&rja);
      amdel(&ra);
      jaP = amdef("ja" ,&rja    ,iaP[numeq],1,"IV");
      aP  = amdef("a"  ,&ra     ,iaP[numeq],1,"DV");
   } 
   /*------------------------------------------------------- receive ja */
   MPI_Recv(jaP,iaP[numeq],MPI_INT,n,tag+2,actintra->MPI_INTRA_COMM,&status);
   /*-------------------------------------------------------- receive a */
   MPI_Recv(aP,iaP[numeq],MPI_DOUBLE,n,tag+3,actintra->MPI_INTRA_COMM,&status);
   /*--------------------------- decrease number of unreceived messages */
   nrecv--;
   /*---------------- make multiplication and add to my piece of outcsr */
#if 1
   for (i=0; i<work->numeq; i++)
   {
      actrow = work->update.a.iv[i];
      mlpcg_extractrowcsr(actrow,work->numeq,work->update.a.iv,work->ia.a.iv,
                          work->ja.a.iv,work->a.a.dv,&col,&rcol,&scol,&ncol);
      if (ncol==0) continue;
      min1 = rcol[0];
      max1 = rcol[ncol-1];
      init_quick_find(rcol,ncol,&shift,bins);
      for (j=0; j<numeq; j++) /* loop all columns in received P */
      {
         actcol = updateP[j];
         mlpcg_extractcolcsc(actcol,numeq,updateP,iaP,jaP,aP,&colP,&rcolP,&scolP,&nrow);
         if (nrow==0) continue;
         /* make the multiplication */
         min2    = rcolP[0];
         max2    = rcolP[nrow-1];
         if (min2 > max1 || min1 > max2 || IMAX(min1,min2)>IMIN(max1,max2))
            continue;
         foundit = 0;
         sum     = 0.0;
         for (k=0; k<nrow; k++)
         {
            if (rcolP[k] < min1 || rcolP[k] > max1)
               continue;
            index = quick_find(rcolP[k],rcol,ncol,shift,bins);
            if (index==-1) continue;
            foundit=1;
            sum += col[index]*colP[k]; 
         }
         if (foundit)
            mlpcg_csr_addentry(outcsr,sum,actrow,actcol,actintra);
      } 
   }
#endif
#if 0
   for (i=0; i<work->numeq; i++)
   {
      actrow = work->update.a.iv[i];
      mlpcg_extractrowcsr(actrow,work->numeq,work->update.a.iv,work->ia.a.iv,
                          work->ja.a.iv,work->a.a.dv,&col,&rcol,&scol,&ncol);
      if (ncol==0) continue;
      for (j=0; j<numeq; j++) /* loop all columns in received P */
      {
         actcol = updateP[j];
         mlpcg_extractcolcsc(actcol,numeq,updateP,iaP,jaP,aP,&colP,&rcolP,&scolP,&nrow);
         if (nrow==0) continue;
         /* make the multiplication */
         foundit = 0;
         sum     = 0.0;
         min1    = rcolP[0];
         max1    = rcolP[nrow-1];
         init_quick_find(rcolP,nrow,&shift,bins);
         for (k=0; k<ncol; k++)
         {
            if (rcol[k] < min1 || rcol[k] > max1)
               continue;
            index = quick_find(rcol[k],rcolP,nrow,shift,bins);
            if (index==-1) continue;
            foundit=1;
            sum += col[k]*colP[index]; 
         }
         if (foundit)
            mlpcg_csr_addentry(outcsr,sum,actrow,actcol,actintra);
      } 
   }
#endif
}
/*----------------------------------------------- wait for all messages */
for (n=0; n<(nproc-1)*4; n++) MPI_Wait(&(request[n]),&status);
/*
t2 = ds_cputime();
if (myrank==0) printf("out  = work * P inter II : %20.10f\n",t2-t1);
fflush(stdout);
*/
/*======================================================================*/
/*                                                       tidy up        */
/*======================================================================*/
CCAFREE(request);
amdel(&rupdate);
amdel(&ria);
amdel(&rja);
amdel(&ra);
#endif
/*--------------- uninitialize the column pointers of the incsr matrix */
mlpcg_extractcollocal_uninit(incsr,colsize,icol,dcol);
CCAFREE(colsize);
CCAFREE(icol);
CCAFREE(dcol);

CCAFREE(rcol);
CCAFREE(col);
CCAFREE(rcolP);
CCAFREE(colP);

#ifdef PARALLEL 
MPI_Barrier(actintra->MPI_INTRA_COMM);
#endif
/*----------------------------------------------------------------------*/
/*------------------------------debugging */
#if 0 /* make dense printout from work */
mlpcg_csr_close(outcsr);
amdef("so",&so,outcsr->numeq_total,outcsr->numeq_total,"DA");
amzero(&so);
amdef("rw",&ro,outcsr->numeq_total,outcsr->numeq_total,"DA");
for (i=0; i<outcsr->numeq; i++)
{
   actrow   = outcsr->update.a.iv[i];
   colstart = outcsr->ia.a.iv[i]; 
   colend   = outcsr->ia.a.iv[i+1];
   for (j=colstart; j<colend; j++)
      so.a.da[actrow][outcsr->ja.a.iv[j]] =  outcsr->a.a.dv[j];
}
#ifdef PARALLEL 
MPI_Allreduce(so.a.da[0],ro.a.da[0],so.fdim*so.sdim,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
#else
for (i=0; i<so.fdim*so.sdim; i++)
   ro.a.da[0][i] = so.a.da[0][i];
#endif
fprintf(allfiles.out_err,"out-Matrix from sparse multiplication\n");
for (i=0; i<ro.fdim; i++)
{
   for (j=0; j<ro.sdim; j++)
   {
    /*  fprintf(allfiles.out_err,"%8.2f ",dense_work.a.da[i][j]);*/
        if (FABS(ro.a.da[i][j])>EPS14)
        fprintf(allfiles.out_err,"X ");
        else
        fprintf(allfiles.out_err,"* ");
   }
   fprintf(allfiles.out_err," \n");
}
fprintf(allfiles.out_err,"out-Matrix from sparse multiplication minus out matrix from seq. dense mult.\n");
for (i=0; i<ro.fdim; i++)
{
   for (j=0; j<ro.sdim; j++)
   {
    /*  fprintf(allfiles.out_err,"%8.2f ",dense_work.a.da[i][j]);*/
        if (FABS(ro.a.da[i][j]-dense_out.a.da[i][j])>EPS12)
        fprintf(allfiles.out_err,"X ");
        else
        fprintf(allfiles.out_err,"* ");
   }
   fprintf(allfiles.out_err," \n");
}
fflush(allfiles.out_err);
/* tidy up debugging */
amdel(&sdense_in);
amdel(&rdense_in);
amdel(&sdense_P);
amdel(&rdense_P);
amdel(&dense_work);
amdel(&dense_out);
amdel(&sw);
amdel(&rw);
amdel(&so);
amdel(&ro);
#endif
/*-------------------------------debugging */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_PtKP */



/*!---------------------------------------------------------------------
\brief make a print of the matrix                                              

<pre>                                                        m.gee 11/02 

</pre>
\param mlpcgvars    MLPCGVARS*   (i)   variables needed for mlpcg                   
\param bdcsr        DBCSR*       (i)   the distributed  csr matrix                   
\param actintra     INTRA*       (i)   the intra-communicator of this field                  
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_printfmatrix(DBCSR     *bdcsr, 
                        INTRA     *actintra)
{
INT        i,j;
ARRAY      sdense,rdense;
INT        actrow,actcol,colstart,colend;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_printfmatrix");
#endif
/*----------------------------------------------------------------------*/
#if 1
amdef("tmp",&sdense,bdcsr->numeq_total,bdcsr->numeq_total,"DA");
amdef("tmp",&rdense,bdcsr->numeq_total,bdcsr->numeq_total,"DA");
amzero(&sdense);
amzero(&rdense);
/* fill matrix */
for (i=0; i<bdcsr->numeq; i++)
{
   actrow = bdcsr->update.a.iv[i];
   colstart = bdcsr->ia.a.iv[i];
   colend   = bdcsr->ia.a.iv[i+1];
   for (j=colstart; j<colend; j++)
   {
      actcol = bdcsr->ja.a.iv[j];
      sdense.a.da[actrow][actcol] = bdcsr->a.a.dv[j];
   }
}
#ifdef PARALLEL 
MPI_Allreduce(sdense.a.da[0],rdense.a.da[0],sdense.fdim*sdense.sdim,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
#else
for (i=0; i<sdense.fdim*sdense.sdim; i++)
   rdense.a.da[0][i] = sdense.a.da[0][i];
#endif
amdel(&sdense);
/* check symmetry */
for (i=0; i<rdense.fdim; i++)
{
   for (j=i+1; j<rdense.sdim; j++)
     if (FABS(rdense.a.da[i][j]-rdense.a.da[j][i])>EPS12)
      fprintf(allfiles.out_err,"K not symmetric in ij %d %d %20.10f ji %d %d %20.10f\n",i,j,rdense.a.da[i][j],j,i,rdense.a.da[j][i]); 
}
/* print to err */
fprintf(allfiles.out_err,"Matrix\n");
for (i=0; i<rdense.fdim; i++)
{
   for (j=0; j<rdense.sdim; j++)
   {
      fprintf(allfiles.out_err,"%16.6f ",rdense.a.da[i][j]);
    /*    if (FABS(rdense.a.da[i][j])>EPS14)
        fprintf(allfiles.out_err,"X ");
        else
        fprintf(allfiles.out_err,"* ");*/
   }
   fprintf(allfiles.out_err," \n");
}
fflush(allfiles.out_err);
amdel(&rdense);
#endif
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_printfmatrix */

/*!---------------------------------------------------------------------
\brief make a print of a vector to gid                                             

<pre>                                                        m.gee 11/02 

</pre>
\param z            DOUBLE*      (i)   the distributed vector to print
\param csr          DBCSR*       (i)   the matching csr matrix
\param fielddis     DISCRET*     (i)   the field
\param partdis      PARTDISCRET* (i)   the partition of the field
\param actintra     INTRA*       (i)   the intra-communicator of this field                  
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_printvec(INT          iter,
                   DOUBLE      *z, 
                   DBCSR*       csr,
                   DISCRET     *fielddis,
                   PARTDISCRET *partdis,
                   INTRA       *actintra)
{
INT        i,dof;
FILE      *out = allfiles.out_err;
ARRAY      send_a,recv_a;
DOUBLE    *zs,*zr;
char       sign='"';
NODE      *actnode;
DOUBLE     x[3],a[3],scal,sdc;
INT        nnode;
nnode = genprob.nnode;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_printvec");
#endif
/*----------------------------------------------------------------------*/
zs = amdef("tmp",&send_a,csr->numeq_total,1,"DV");
zr = amdef("tmp",&recv_a,csr->numeq_total,1,"DV");
amzero(&send_a);
for (i=0; i<csr->numeq; i++)
{
   dof = csr->update.a.iv[i];
   zs[dof] = z[i];
}
#ifdef PARALLEL 
MPI_Allreduce(zs,zr,csr->numeq_total,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
#else
for (i=0; i<csr->numeq_total; i++) zr[i] = zs[i];
#endif
amdel(&send_a);
/*----------------------------------------------- make printout of head */
fprintf(out,"#-------------------------------------------------------------------------------\n");
fprintf(out,"# RESULT displacement on FIELD structure\n");
fprintf(out,"#-------------------------------------------------------------------------------\n");
fprintf(out,"RESULT %cdisplacement%c %cpcarat%c %d VECTOR ONNODES\n",sign,sign,sign,sign,iter);
fprintf(out,"RESULTRANGESTABLE %cstandard_structure%c\n",sign,sign);
fprintf(out,"COMPONENTNAMES %cx-displ%c,%cy-displ%c,%cz-displ%c\n",sign,sign,sign,sign,sign,sign);
fprintf(out,"VALUES\n");
#ifdef D_SHELL8
sdc = fielddis->element[0].e.s8->sdc;
#endif
for (i=0; i<fielddis->numnp; i++)
{
actnode = &(fielddis->node[i]);
scal    = 1.0;
/* mid surface */
if (actnode->dof[0] < csr->numeq_total) x[0] = zr[actnode->dof[0]];
else                                    x[0] = 0.0;
if (actnode->dof[1] < csr->numeq_total) x[1] = zr[actnode->dof[1]];
else                                    x[1] = 0.0;
if (actnode->dof[2] < csr->numeq_total) x[2] = zr[actnode->dof[2]];
else                                    x[2] = 0.0;
/* director */
if (actnode->dof[3] < csr->numeq_total) a[0] = zr[actnode->dof[3]];
else                                    a[0] = 0.0;
if (actnode->dof[4] < csr->numeq_total) a[1] = zr[actnode->dof[4]];
else                                    a[1] = 0.0;
if (actnode->dof[5] < csr->numeq_total) a[2] = zr[actnode->dof[5]];
else                                    a[2] = 0.0;
/* lower surface */
fprintf(out," %6d %23.15E %23.15E %23.15E\n",actnode->Id+1,
                                          x[0]-a[0]*scal/sdc,
                                          x[1]-a[1]*scal/sdc,
                                          x[2]-a[2]*scal/sdc);
/* upper surface */
fprintf(out," %6d %23.15E %23.15E %23.15E\n",actnode->Id+1+nnode,
                                          x[0]+a[0]*scal/sdc,
                                          x[1]+a[1]*scal/sdc,
                                          x[2]+a[2]*scal/sdc);
}
fprintf(out,"END VALUES\n");
fflush(out);
amdel(&recv_a);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_printvec */



/*! @} (documentation module close)*/
