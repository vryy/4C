/*!---------------------------------------------------------------------
\file
\brief contains the multilevel preconditioner for shells

---------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/prototypes_mlpcg.h"
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
\brief ilu(n) smoother with overlap                                             

<pre>                                                        m.gee 11/02 

</pre>
\param z            double*      (o)   the solution of the smoothing
\param r            double*      (i)   the right hand side
\param csr          DBCSR*       (i)   the matrix to smooth with
\param nsweep       int          (i)   n in ilu(n)
\param actintra     INTRA*       (i)   the intra-communicator of this field                  
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_precond_smo_ILUn_overlap(double *z, double *r, DBCSR *csr, int nsweep, INTRA *actintra)
{
int            i,j,k,n;
int            dof,index,nrequest;
int            myrank,nproc;
DBCSR         *ilu;
DBCSR         *asm;
ARRAY          levs,w,jw;
int            size;
int            ierr;

double        *rwork;
double        *zwork;

int          **irecv;
ARRAY          irecv_a;
double        *drecv;
ARRAY          drecv_a;
double       **dsend;
ARRAY          dsend_a;

#ifdef PARALLEL
MPI_Request   *request;
MPI_Status     status;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_smo_ILUn_overlap");
#endif
/*----------------------------------------------------------------------*/
myrank = actintra->intra_rank;
nproc  = actintra->intra_nprocs;
/*-------------- do the decomposition of the matrix, if not done before */
if (csr->ilu==NULL)
{
   csr->ilu = (DBCSR*)CCACALLOC(1,sizeof(DBCSR));
   csr->asm = (DBCSR*)CCACALLOC(1,sizeof(DBCSR));
   ilu      = csr->ilu;
   asm      = csr->asm;
   /*------------------------- find the ghost rows due to overlap of ilu */
   mlpcg_csr_overlap(csr,asm,ilu,mlprecond.overlap,actintra);
   /*---------------- change the enumeration to local and fortran style */
   mlpcg_csr_localnumsf_overlap(asm);
   /*---------------------- allocate space for the decomposition in ilu */
   if (nsweep==0) size = (asm->a.fdim+1);
   if (nsweep==1) size = (int)((asm->a.fdim+1)*2.0);
   if (nsweep==2) size = (int)((asm->a.fdim+1)*2.5);
   if (nsweep==3) size = (int)((asm->a.fdim+1)*3.0);
   if (nsweep==4) size = (int)((asm->a.fdim+1)*4.0);
   if (nsweep==5) size = (int)((asm->a.fdim+1)*5.0);
   if (nsweep==6) size = (int)((asm->a.fdim+1)*6.0);
   if (nsweep>=7) size = (int)((asm->a.fdim+1)*7.0);
   tryagain:
   amdef("ilu_val"  ,&(ilu->a) ,size          ,1,"DV");
   amdef("ilu_bindx",&(ilu->ja),size          ,1,"IV");
   amdef("ilu_ia"   ,&(ilu->ia),asm->numeq+1  ,1,"IV");
   amdef("levs"     ,&levs     ,size          ,1,"IV");
   amdef("w"        ,&w        ,asm->numeq    ,1,"DV");
   amdef("jw"       ,&jw       ,3*(asm->numeq),1,"IV");
   /*------------------------------------ call the ilu(k) factorization */
   ierr = 1;
   i    = size;
   iluk(&(asm->numeq),
        asm->a.a.dv,
        asm->ja.a.iv,
        asm->ia.a.iv,
        &nsweep,
        ilu->a.a.dv,
        ilu->ja.a.iv,
        ilu->ia.a.iv,
        levs.a.iv,
        &i,
        w.a.dv,
        jw.a.iv,
        &ierr);
/*
c ierr    = integer. Error message with the following meaning.
c           ierr  = 0    --> successful return.
c           ierr .gt. 0  --> zero pivot encountered at step number ierr.
c           ierr  = -1   --> Error. input matrix may be wrong.
c                            (The elimination process has generated a
c                            row in L or U whose length is .gt.  n.)
c           ierr  = -2   --> The matrix L overflows the array al.
c           ierr  = -3   --> The matrix U overflows the array alu.
c           ierr  = -4   --> Illegal value for lfil.
c           ierr  = -5   --> zero row encountered in A or U.
*/        
   if (ierr != 0)
   {
      if (ierr>0)
      dserror("Zero pivot in ilu(k)");
      if (ierr==-1)
      dserror("Fatal error in ilu(k)");
      if (ierr==-4)
      dserror("Illegal value for fill-in");
      if (ierr==-2 || ierr==-3)
      {
         printf("rank %d: Enlargment of storage for ilu happened\n",myrank);
         size = (int)(size*1.3);
         amdel(&(ilu->a) );
         amdel(&(ilu->ja));
         amdel(&(ilu->ia));
         amdel(&levs     );
         amdel(&w        );
         amdel(&jw       );
         goto tryagain;
      } 
   }
   /*----------------------------------- set flag, that ilu is factored */
   ilu->is_factored = mlprecond.ncall;
   /*---------------------------------------------------------- tidy up */
   amdel(&levs);
   amdel(&w);
   amdel(&jw);
}
else if (csr->ilu->is_factored != mlprecond.ncall && mlprecond.mod==0)
{
   mlpcg_csr_destroy(asm);
   csr->asm = CCAFREE(asm);
   csr->asm = (DBCSR*)CCACALLOC(1,sizeof(DBCSR));
   ilu      = csr->ilu;
   asm      = csr->asm;
   /*------------------------- find the ghost rows due to overlap of ilu */
   mlpcg_csr_overlap(csr,asm,ilu,mlprecond.overlap,actintra);
   /*---------------- change the enumeration to local and fortran style */
   mlpcg_csr_localnumsf_overlap(asm);
   /*---------------------- allocate space for the decomposition in ilu */
   if (nsweep==0) size = (asm->a.fdim+1);
   if (nsweep==1) size = (int)((asm->a.fdim+1)*2.0);
   if (nsweep==2) size = (int)((asm->a.fdim+1)*2.5);
   if (nsweep==3) size = (int)((asm->a.fdim+1)*3.0);
   if (nsweep==4) size = (int)((asm->a.fdim+1)*4.0);
   if (nsweep==5) size = (int)((asm->a.fdim+1)*5.0);
   if (nsweep==6) size = (int)((asm->a.fdim+1)*6.0);
   if (nsweep>=7) size = (int)((asm->a.fdim+1)*7.0);
   tryagain2:
   amdef("levs"     ,&levs     ,size          ,1,"IV");
   amdef("w"        ,&w        ,asm->numeq    ,1,"DV");
   amdef("jw"       ,&jw       ,3*(asm->numeq),1,"IV");
   /*------------------------------------ call the ilu(k) factorization */
   ierr = 1;
   i    = size;
   iluk(&(asm->numeq),
        asm->a.a.dv,
        asm->ja.a.iv,
        asm->ia.a.iv,
        &nsweep,
        ilu->a.a.dv,
        ilu->ja.a.iv,
        ilu->ia.a.iv,
        levs.a.iv,
        &i,
        w.a.dv,
        jw.a.iv,
        &ierr);
/*
c ierr    = integer. Error message with the following meaning.
c           ierr  = 0    --> successful return.
c           ierr .gt. 0  --> zero pivot encountered at step number ierr.
c           ierr  = -1   --> Error. input matrix may be wrong.
c                            (The elimination process has generated a
c                            row in L or U whose length is .gt.  n.)
c           ierr  = -2   --> The matrix L overflows the array al.
c           ierr  = -3   --> The matrix U overflows the array alu.
c           ierr  = -4   --> Illegal value for lfil.
c           ierr  = -5   --> zero row encountered in A or U.
*/        
   if (ierr != 0)
   {
      if (ierr>0)
      dserror("Zero pivot in ilu(k)");
      if (ierr==-1)
      dserror("Fatal error in ilu(k)");
      if (ierr==-4)
      dserror("Illegal value for fill-in");
      if (ierr==-2 || ierr==-3)
      {
         printf("rank %d: Enlargment of storage for ilu happened\n",myrank);
         size = (int)(size*1.3);
         amdel(&(ilu->a) );
         amdel(&(ilu->ja));
         amdef("ilu_val"  ,&(ilu->a) ,size          ,1,"DV");
         amdef("ilu_bindx",&(ilu->ja),size          ,1,"IV");
         amdel(&levs     );
         amdel(&w        );
         amdel(&jw       );
         goto tryagain2;
      } 
   }
   /*----------------------------------- set flag, that ilu is factored */
   ilu->is_factored = mlprecond.ncall;
   /*---------------------------------------------------------- tidy up */
   amdel(&levs);
   amdel(&w);
   amdel(&jw);
}
/*----------------------------------------------------------------------*/
ilu = csr->ilu;
/*----------------------------------------- allocate approbiate z and r */
rwork = (double*)CCACALLOC(ilu->numeq,sizeof(double));
zwork = (double*)CCACALLOC(ilu->numeq,sizeof(double));
/*------------------------------------------ fill my sendbuffers from r */
for (n=0; n<nproc; n++)
{
   if (n==myrank || ilu->computebuff.a.ia[myrank][n]==0) continue;
   ilu->sendbuff.a.da[n][0] = 0.0;
   for (i=0; i<ilu->gdofsend.a.ia[n][0]; i++)
   {
      dof   = ilu->gdofsend.a.ia[n][i+1];
      index = find_index(dof,csr->update.a.iv,csr->numeq);
      if (index==-1) dserror("Cannot find local dof");
      ilu->sendbuff.a.da[n][i+1] = r[index]; 
   }
}
/*----------------------------------------------------------- make send */
nrequest=0;
for (n=0; n<nproc; n++) nrequest += ilu->computebuff.a.ia[myrank][n];
nrequest *= 2;
request = (MPI_Request*)CCAMALLOC(nrequest*sizeof(MPI_Request));
j = 0;
for (n=0; n<nproc; n++)
{
   if (n==myrank || ilu->computebuff.a.ia[myrank][n]==0) continue;
   MPI_Isend(ilu->gdofsend.a.ia[n],ilu->gdofsend.a.ia[n][0]+1,MPI_INT   ,n,0,actintra->MPI_INTRA_COMM,&(request[j]));
   j++;
   MPI_Isend(ilu->sendbuff.a.da[n],ilu->gdofsend.a.ia[n][0]+1,MPI_DOUBLE,n,1,actintra->MPI_INTRA_COMM,&(request[j]));
   j++;
}
dsassert(j==nrequest,"Number of send wrong");
/*-------------------------------------- put my own piece of r to rwork */
/* find the first entry */
index = find_index(csr->update.a.iv[0],ilu->update.a.iv,ilu->numeq);
if (index==-1) dserror("Cannot find local dof");
/* put values from r to rwork in the right place */
for (i=0; i<csr->numeq; i++)
   rwork[index+i] = r[i];
/*---------------------------------- receive the incoming overlap parts */
irecv = amdef("irecv",&irecv_a,nproc,10000,"IA");
drecv = amdef("drecv",&drecv_a,10000,1,"DV");
for (n=0; n<nproc; n++)
{
   if (n==myrank || ilu->computebuff.a.ia[n][myrank]==0) continue;
   MPI_Probe(n,0,actintra->MPI_INTRA_COMM,&status);
   MPI_Get_count(&status,MPI_INT,&k);
   if (k>irecv_a.sdim)
   {
      irecv = amredef(&irecv_a,nproc,k,"IA");
      amdel(&drecv_a);
      drecv = amdef("drecv",&drecv_a,k,1,"DV");
   }
   MPI_Recv(irecv[n],k,MPI_INT   ,n,0,actintra->MPI_INTRA_COMM,&status);
   MPI_Recv(drecv   ,k,MPI_DOUBLE,n,1,actintra->MPI_INTRA_COMM,&status);
   /* put received values to my rwork */
   k--;
   for (i=0; i<k; i++)
   {
      index = find_index(irecv[n][i+1],ilu->update.a.iv,ilu->numeq);
      if (index==-1) dserror("Cannot find dof");
      rwork[index] += drecv[i+1];
   }
}
/*-------------------------------------- make solve with the ilu matrix */
lusol(&(ilu->numeq),rwork,zwork,ilu->a.a.dv,ilu->ja.a.iv,ilu->ia.a.iv);
/*-------------------------------------------- wait for sends to finish */
for (i=0; i<nrequest; i++) MPI_Wait(&(request[i]),&status);
/*-------------------------------------- fill my sendbuffers from zwork */
dsend = amdef("dsend",&dsend_a,nproc,irecv_a.sdim,"DA");
for (n=0; n<nproc; n++)
{
   if (n==myrank || ilu->computebuff.a.ia[n][myrank]==0) continue;
   for (i=0; i<irecv[n][0]; i++)
   {
      dof   = irecv[n][i+1];
      index = find_index(dof,ilu->update.a.iv,ilu->numeq);
      if (index==-1) dserror("Cannot find overlap dof in ilu solution");
      dsend[n][i+1] = zwork[index];
   }
}
/*----------------------------------------------------------- make send */
j=0;
for (n=0; n<nproc; n++)
{
   if (n==myrank || ilu->computebuff.a.ia[n][myrank]==0) continue;
   MPI_Isend(dsend[n],irecv[n][0]+1,MPI_DOUBLE,n,10,actintra->MPI_INTRA_COMM,&(request[j]));
   j++;
}
nrequest = j;
/*---------------------------------------- sort my own values back to z */
/* find the first entry */
index = find_index(csr->update.a.iv[0],ilu->update.a.iv,ilu->numeq);
if (index==-1) dserror("Cannot find local dof");
/* put values from zwork to z */
for (i=0; i<csr->numeq; i++)
   z[i] = zwork[index+i];
/*--------------------------------- receive incoming parts and put to z */
for (n=0; n<nproc; n++)
{
   if (n==myrank || ilu->computebuff.a.ia[myrank][n]==0) continue;
   MPI_Probe(n,10,actintra->MPI_INTRA_COMM,&status);
   MPI_Get_count(&status,MPI_DOUBLE,&k);
   if (k-1 != ilu->gdofsend.a.ia[n][0]) dserror("Message size wrong");
   MPI_Recv(ilu->sendbuff.a.da[myrank],k,MPI_DOUBLE,n,10,actintra->MPI_INTRA_COMM,&status);
   /* put received values to my z */
   for (i=0; i<ilu->gdofsend.a.ia[n][0]; i++)
   {
      dof   = ilu->gdofsend.a.ia[n][i+1];
      index = find_index(dof,csr->update.a.iv,csr->numeq);
      if (index==-1) dserror("Cannot find local dof");
      z[index] += ilu->sendbuff.a.da[myrank][i+1];
   }
}
/*------------------------------------------ wait for request to finish */
for (i=0; i<nrequest; i++) MPI_Wait(&(request[i]),&status);
CCAFREE(request);
/*----------------------------------------------------------------------*/
rwork = CCAFREE(rwork);
zwork = CCAFREE(zwork);
amdel(&irecv_a);
amdel(&drecv_a);
amdel(&dsend_a);
#endif /* end of #ifdef PARALLEL */
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_smo_ILUn_overlap */



/*!---------------------------------------------------------------------
\brief get ghost rows to a given overlap                                          

<pre>                                                        m.gee 12/02 

</pre>
\param csr        DBCSR*          (i)   the csr matrix
\param ocsr       DBCSR*          (o)   the csr matrix which overlaps
\param ilu        DBCSR*          (o)   the csr matrix which overlaps and will hold the ilu 
\param overlap    int*            (i)   degree of overlap
\param oupdate    int**           (o)   adress of the overalping update vector pointer
\param actintra   INTRA*          (i)   the communicator                 
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_csr_overlap(DBCSR *csr, DBCSR *ocsr, DBCSR *ilu, int overlap, INTRA *actintra)
{
int       i,j,k,n,counter;
int       fcd;
int       fcdindex;
int       index;
int       owner;
int       nrequest;
int       size;

int       actrow;
int       actcol;
int       colstart;
int       colend;

int       sendtos[MAXPROC][MAXPROC],sendtor[MAXPROC][MAXPROC];

int       nupdatesend;
int     **updatesend;
ARRAY     updatesend_a;

int       niasend;
int     **iasend;
ARRAY     iasend_a;

int       njasend;
int     **jasend;
ARRAY     jasend_a;

int       nasend;
double  **asend;
ARRAY     asend_a;

int       noupdate;
int      *oupdate;
ARRAY     oupdate_a;

int       nupdaterecv;
int      *updaterecv;
int       niarecv;
int      *iarecv;
int       njarecv;
int      *jarecv;
int       narecv;
double   *arecv;

int       myrank,nproc;

int      *update;
int      *ia;
int      *ja;
double   *a;
int       numeq;

int     *irecv;
ARRAY    irecv_a;

int      rowguess;
int      nnzguess;
int      ione=-1;

#ifdef PARALLEL
MPI_Request   *request;
MPI_Status     status;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_csr_overlap");
#endif
/*----------------------------------------------------------------------*/
ia      = csr->ia.a.iv;
ja      = csr->ja.a.iv;
update  = csr->update.a.iv;
a       = csr->a.a.dv;
numeq   = csr->numeq;
/*----------------------------------------------------------------------*/
myrank = actintra->intra_rank;
nproc  = actintra->intra_nprocs;
/*--------------------------------------- make a neighbourhood relation */
for (i=0; i<nproc; i++)
for (j=0; j<nproc; j++) sendtos[i][j] = 0;
/*---------------------------------------------------- get coupled rows */
fcd      = csr->firstcoupledof;
fcdindex = mlpcg_getindex(fcd,update,numeq);
if (fcdindex==-1) dserror("Cannot find local dof");
/* loop coupled rows */
for (i=ia[fcdindex]; i<ia[numeq]; i++)
{
   owner = mlpcg_getowner(ja[i],csr->owner,nproc);
   if (owner != myrank) sendtos[myrank][owner] = 1;
}
MPI_Allreduce(&(sendtos[0][0]),&(sendtor[0][0]),MAXPROC*MAXPROC,MPI_INT,MPI_SUM,actintra->MPI_INTRA_COMM);
/* count the requests */
nrequest=0;
for (i=0; i<nproc; i++) nrequest += sendtor[myrank][i];
if (nrequest == 0) dserror("Something wrong in overlap of ILU");
/*---------------------------------------------- prepare my send pieces */
/*
NOTE:
- The piece of my stiffness matrix is different for every neighbour
- Overlap is at least 1, so find the rows 
*/
/* make a guess */
nupdatesend = numeq-fcdindex;
updatesend  = amdef("updsend",&updatesend_a,nproc,nupdatesend+1,"IA");
              amzero(&updatesend_a);
for (n=0; n<nproc; n++)
{
   if (n==myrank || sendtor[myrank][n]==0) continue;
   for (i=fcdindex; i<numeq; i++)
   {
      actrow = update[i];
      colstart = ia[i];
      colend   = ia[i+1];
      for (j=colstart; j<colend; j++)
      {
         actcol = ja[j];
         owner  = mlpcg_getowner(actcol,csr->owner,nproc);
         if (owner==n)
         {
            if (updatesend[n][0]>=nupdatesend)
            {
               updatesend = amredef(&updatesend_a,nproc,(int)(1.5*nupdatesend+1),"IA");
               nupdatesend = (int)(1.5*nupdatesend);
            }
            updatesend[n][updatesend[n][0]+1] = actrow;
            updatesend[n][0]++;
            break;
         }
      }
   }
}
/* make more overlap (overlap>1) */
for (i=1; i<overlap; i++)
{
   for (n=0; n<nproc; n++)
   {
      if (n==myrank || sendtor[myrank][n]==0) continue;
      /* loop all dofs in updatesend[k] */
      size = updatesend[n][0];
      for (j=0; j<size; j++)
      {
         actrow = updatesend[n][j+1];
         index  = mlpcg_getindex(actrow,update,numeq);
         if (index==-1) dserror("Cannot find local dof");
         colstart = ia[index];
         colend   = ia[index+1];
         for (k=colstart; k<colend; k++)
         {
            actcol = ja[k];
            index  = find_index(actcol,&(updatesend[n][1]),updatesend[n][0]);
            if (index != -1) continue;
            owner  = mlpcg_getowner(actcol,csr->owner,nproc);
            if (owner != myrank) continue;
            if (updatesend[n][0]>=nupdatesend)
            {
               updatesend = amredef(&updatesend_a,nproc,(int)(3.0*nupdatesend+1),"IA");
               nupdatesend = (int)(3.0*nupdatesend);
            }
            updatesend[n][updatesend[n][0]+1] = actcol;
            updatesend[n][0]++;
         }
         
      }
      /* delete the doubles */
      for (j=0; j<updatesend[n][0]; j++)
      {
         actcol = updatesend[n][j+1];
         for (k=j+1; k<updatesend[n][0]; k++)
            if (updatesend[n][k+1] == actcol)
               updatesend[n][k+1] = -1;
      }
      /* move remaining values to the front */
      counter = 0;
      for (j=0; j<updatesend[n][0]; j++)
      {
         if (updatesend[n][j+1]!=-1)
         {
            updatesend[n][counter+1] = updatesend[n][j+1];
            counter++;
         }
      }
      updatesend[n][0] = counter;
      mg_sort(&(updatesend[n][1]),counter,NULL,NULL);
   }/* end of for (n=0; n<nproc; n++) */
}/* end of for (i=1; i<overlap; i++)*/
/*---------------------------------- now send to those who requested it */
request = (MPI_Request*)CCAMALLOC(nrequest*sizeof(MPI_Request));
k=0;
for (n=0; n<nproc; n++)
{
   if (n==myrank || sendtor[n][myrank]==0) continue;
   MPI_Isend(updatesend[n],updatesend[n][0]+1,MPI_INT,n,myrank,actintra->MPI_INTRA_COMM,&(request[k]));
   k++;
}
if (k != nrequest) dserror("Number of sends wrong");
/*---------------- make a guess how large my overlapping update will be */
noupdate = (int)(numeq*3.0);
oupdate = amdef("oupdate",&oupdate_a,noupdate,1,"IV");
/*------------------------------------------ put my own dofs to oupdate */
counter=0;
for (i=0; i<numeq; i++)
   oupdate[counter++] = update[i];
/*----------------------------------------- receive the messages for me */
irecv = amdef("irecv",&irecv_a,1,1,"IV");
for (n=0; n<nproc; n++)
{
   if (n==myrank || sendtor[myrank][n]==0) continue;
   /* probe for message from n with tag n */
   MPI_Probe(n,n,actintra->MPI_INTRA_COMM,&status);
   /* get length */
   MPI_Get_count(&status,MPI_INT,&k);
   if (k>irecv_a.fdim)
   {
      amdel(&irecv_a);
      irecv = amdef("irecv",&irecv_a,k,1,"IV");
   }
   MPI_Recv(irecv,k,MPI_INT,n,n,actintra->MPI_INTRA_COMM,&status);
   /* put the values to my oupdate */
   for (i=0; i<irecv[0]; i++)
   {
      if (counter>=noupdate)
      {
         noupdate = (int)(numeq*3.0);
         oupdate = amredef(&oupdate_a,noupdate,1,"IV");
      }
      oupdate[counter++] = irecv[i+1];
   }
}
noupdate = counter;
oupdate  = amredef(&oupdate_a,noupdate,1,"IV");
mg_sort(oupdate,noupdate,NULL,NULL);
/*--------------------------------------------------- wait for requests */
for (i=0; i<nrequest; i++) MPI_Wait(&(request[i]),&status);
CCAFREE(request);
/*---------------------------------------------- create the matrix ocsr */
rowguess    = csr->ja.fdim / csr->numeq;
rowguess    = (int)(rowguess*1.2);
nnzguess    = rowguess * noupdate;
ocsr->numeq = noupdate;
amdef("update",&(ocsr->update),noupdate  ,1,"IV");
amdef("ia"    ,&(ocsr->ia)    ,noupdate+1,1,"IV");
amdef("ja"    ,&(ocsr->ja)    ,nnzguess  ,1,"IV");
amdef("a"     ,&(ocsr->a)     ,nnzguess  ,1,"DV");
aminit(&(ocsr->ja),(void*)(&ione));
counter=0;
for (i=0; i<ocsr->numeq; i++) 
{
   ocsr->ia.a.iv[i] = counter;
   counter += rowguess;
}
   ocsr->ia.a.iv[i] = counter;
for (i=0; i<noupdate; i++) 
   ocsr->update.a.iv[i] = oupdate[i];
for (i=0; i<nproc; i++)
{
   ocsr->owner[i][0] = csr->owner[i][0];
   ocsr->owner[i][1] = csr->owner[i][1];
}
/*------------------ send my piece of overlap to everyone who wanted it */
/* find maximum of updatesend and a guess for nnz */
counter=0;
for (n=0; n<nproc; n++)
   if (counter<updatesend[n][0]) counter = updatesend[n][0];
/* allocate some arrays */
iasend = amdef("iasend",&iasend_a,nproc,counter+1,"IA");
jasend = amdef("jasend",&jasend_a,nproc,5000,"IA");
 asend = amdef("asend" ,&asend_a ,nproc,5000,"DA");
for (n=0; n<nproc; n++)
{
   if (n==myrank || sendtor[myrank][n]==0) continue;
   nupdatesend = updatesend[n][0];
   niasend     = nupdatesend+1;
   njasend     = jasend_a.sdim;
   nasend      = asend_a.sdim;
   counter     = 0;
   /* fill the sendbuffer */
   for (i=0; i<nupdatesend; i++)
   {
      actrow       = updatesend[n][i+1];
      index        = find_index(actrow,update,numeq);
      if (index==-1) dserror("Cannot find local dof");
      colstart     = ia[index];
      colend       = ia[index+1];
      iasend[n][i] = counter;
      for (j=colstart; j<colend; j++)
      {
         if (counter>=njasend)
         {
            njasend = (int)(njasend*1.5);
            nasend   = njasend;
            jasend   = amredef(&jasend_a,nproc,njasend,"IA");
            asend    = amredef(&asend_a ,nproc,nasend ,"DA");
         }
         jasend[n][counter] = ja[j];
         asend[n][counter]  = a[j];
         counter++;
      }
   }
   iasend[n][i] = counter;
}/* end of for (n=0; n<nproc; n++) */
/*--------------------------- now send these pieces of stiffness matrix */
nrequest *= 4;
request = (MPI_Request*)CCAMALLOC(nrequest*sizeof(MPI_Request));
counter = 0;
for (n=0; n<nproc; n++)
{
   if (n==myrank || sendtor[myrank][n]==0) continue;
   MPI_Isend(&(updatesend[n][1]),updatesend[n][0]           ,MPI_INT   ,n,0,actintra->MPI_INTRA_COMM,&(request[counter]));
   counter++;
   MPI_Isend(iasend[n]          ,updatesend[n][0]+1         ,MPI_INT   ,n,1,actintra->MPI_INTRA_COMM,&(request[counter]));
   counter++;
   MPI_Isend(jasend[n]          ,iasend[n][updatesend[n][0]],MPI_INT   ,n,2,actintra->MPI_INTRA_COMM,&(request[counter]));
   counter++;
   MPI_Isend(asend[n]           ,iasend[n][updatesend[n][0]],MPI_DOUBLE,n,3,actintra->MPI_INTRA_COMM,&(request[counter]));
   counter++;
}
if (counter != nrequest) dserror("Number of sends wrong");
/*---------------------------------- put my own piece of matrix to ocsr */
for (i=0; i<numeq; i++)
{
   actrow   = update[i];
   colstart = ia[i];
   colend   = ia[i+1];
   for (j=colstart; j<colend; j++)
   {
      actcol = ja[j];
      owner  = mlpcg_getowner(actcol,csr->owner,nproc);
      if (owner != myrank)
      {
         owner = find_index(actcol,ocsr->update.a.iv,ocsr->update.fdim);
         if (owner==-1) continue;
      }
      mlpcg_csr_setentry_overlap(ocsr,a[j],actrow,actcol,actintra);
   }
}
/*--------------------------------- receive the overlap and put to ocsr */
for (n=0; n<nproc; n++)
{
   if (n==myrank || sendtor[n][myrank]==0) continue;
   MPI_Probe(n,0,actintra->MPI_INTRA_COMM,&status);
   MPI_Get_count(&status,MPI_INT,&nupdaterecv);
   niarecv = nupdaterecv+1;
   updaterecv = (int*)CCAMALLOC(nupdaterecv*sizeof(int));
   iarecv     = (int*)CCAMALLOC(niarecv    *sizeof(int));
   MPI_Recv(updaterecv,nupdaterecv,MPI_INT,n,0,actintra->MPI_INTRA_COMM,&status);
   MPI_Recv(iarecv    ,niarecv    ,MPI_INT,n,1,actintra->MPI_INTRA_COMM,&status);
   MPI_Probe(n,2,actintra->MPI_INTRA_COMM,&status);
   MPI_Get_count(&status,MPI_INT,&njarecv);
   if (njarecv != iarecv[nupdaterecv]) dserror("Something wrong with message length");
   narecv = njarecv;
   jarecv = (int*)CCAMALLOC(njarecv*sizeof(int));
   arecv = (double*)CCAMALLOC(narecv*sizeof(double));
   MPI_Recv(jarecv,njarecv,MPI_INT   ,n,2,actintra->MPI_INTRA_COMM,&status);
   MPI_Recv(arecv ,narecv ,MPI_DOUBLE,n,3,actintra->MPI_INTRA_COMM,&status);
   /* put the received stuff to my ocsr matrix */
   for (i=0; i<nupdaterecv; i++)
   {
      actrow = updaterecv[i];
      owner  = find_index(actrow,ocsr->update.a.iv,ocsr->update.fdim);
      if (owner==-1) dserror("Received something that is not in my overlapped ocsr");
      colstart = iarecv[i];
      colend   = iarecv[i+1];
      for (j=colstart; j<colend; j++)
      {
         actcol = jarecv[j];
         owner  = mlpcg_getowner(actcol,csr->owner,nproc);
         if (owner != myrank)
         {
            owner = find_index(actcol,ocsr->update.a.iv,ocsr->update.fdim);
            if (owner==-1) continue;
         }
         mlpcg_csr_setentry_overlap(ocsr,arecv[j],actrow,actcol,actintra);
      }
   }
   updaterecv = CCAFREE(updaterecv);
   iarecv     = CCAFREE(iarecv);
   jarecv     = CCAFREE(jarecv);
   arecv      = CCAFREE(arecv);
}
/*--------------------------------- wait for incomplete sends to finish */
for (i=0; i<nrequest; i++) MPI_Wait(&(request[i]),&status);
CCAFREE(request);
/*---------------------------------------------------- close the matrix */
mlpcg_csr_close(ocsr);
/*--------------- in ilu, create the necessary send and receive buffers */
/*---------------------------- the array updatesend_a is much oversized */
if (ilu->gdofsend.Typ != cca_IA)
{
   counter=0;
   for (n=0; n<nproc; n++)
      if (counter<updatesend[n][0]) 
         counter = updatesend[n][0];
   counter++;
   amredef(&updatesend_a,nproc,counter,"IA");
   am_alloc_copy(&updatesend_a,&(ilu->gdofsend));
}
/*-------- sendbuff will be the corresponding dounble values sendbuffer */
if (ilu->sendbuff.Typ != cca_DA)
   amdef("sbuff",&(ilu->sendbuff),nproc,updatesend_a.sdim,"DA");
/*-------------- computebuff holds the who-sends-who-what array sendtor */
if (ilu->computebuff.Typ != cca_IA)
{
   amdef("senddof",&(ilu->computebuff),nproc,nproc,"IA");
   for (i=0; i<nproc; i++)
   for (j=0; j<nproc; j++)
      ilu->computebuff.a.ia[i][j] = sendtor[i][j];
}
/*------------------------------------------- make update vector in ilu */
ilu->numeq = noupdate;
if (ilu->update.Typ != cca_IV)
   am_alloc_copy(&(ocsr->update),&(ilu->update));
/*----------------------------------------------------------------------*/
amdel(&oupdate_a);
amdel(&irecv_a);
amdel(&updatesend_a);
amdel(&iasend_a);
amdel(&jasend_a);
amdel(&asend_a);
#endif /* end of #ifdef PARALLEL */
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_csr_overlap */

/*!---------------------------------------------------------------------
\brief change global dofs to local dofs in fortran style                                        

<pre>                                                        m.gee 11/02 
</pre>
\param matrix      DBCSR*    (i) the csr to be extracted from
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_csr_localnumsf_overlap(DBCSR *matrix)
{
int        i;
int        numeq,*ia,*ja,*update,index;
int        nnz;
int        shift,*bins;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_csr_localnumsf_overlap");
#endif
/*----------------------------------------------------------------------*/
numeq  = matrix->numeq;
update = matrix->update.a.iv;
ia     = matrix->ia.a.iv;
ja     = matrix->ja.a.iv;
nnz    = ia[numeq];
bins   = (int*)CCAMALLOC(((int)(numeq/4+5))*sizeof(int));
init_quick_find(update,numeq,&shift,bins);
/*----------------------------------------------------------------------*/
for (i=0; i<nnz; i++)
{
    index = quick_find(ja[i],update,numeq,shift,bins);
    if (index==-1)
       dserror("Cannot find local dof");
    ja[i] = index+1;
}

for (i=0; i<=numeq; i++)
   ia[i]++;
/*----------------------------------------------------------------------*/
CCAFREE(bins);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_csr_localnumsf_overlap */

/*!---------------------------------------------------------------------
\brief overlap mutliply with the preconditioner y += A tilde * x * fac                                            

<pre>                                                        m.gee 1/03 

</pre>
\param y            double*      (o)   the solution of y += A*x*fac
\param A            DBCSR*       (i)   the matrix 
\param x            double*      (i)   the right hand side
\param fac          double       (i)   the factor fac
\param init         int          (i)   init=1-> y = A*x*fac init=0->y += A*x*fac
\param actintra     INTRA*       (i)   the intra-communicator of this field                  
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_matvec_asm_overlap(double       *y, 
                              DBCSR        *A,
                              double       *x,
                              double        fac,
                              int           init,
                              INTRA        *actintra)
{
int            i,j,k,n,dof,index;
int            myrank,nproc;
double         sum;
DBCSR         *ilu;
DBCSR         *asm;
int           *update,*ia,*ja,numeq;
double        *a;

int            nrequest;

double        *xwork;
double        *ywork;

int          **irecv;
ARRAY          irecv_a;
double        *drecv;
ARRAY          drecv_a;
double       **dsend;
ARRAY          dsend_a;

#ifdef PARALLEL
MPI_Request   *request;
MPI_Status     status;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_matvec_asm_overlap");
#endif
/*----------------------------------------------------------------------*/
myrank = actintra->intra_rank;
nproc  = actintra->intra_nprocs;
/*----------------------------------------------------------------------*/
asm    = A->asm;
ilu    = A->ilu;
numeq  = ilu->numeq;
update = asm->update.a.iv;
ia     = asm->ia.a.iv;
ja     = asm->ja.a.iv;
a      = asm->a.a.dv;
/*----------------------------------------------------------------------*/
xwork = (double*)CCACALLOC(ilu->numeq,sizeof(double));
ywork = (double*)CCACALLOC(ilu->numeq,sizeof(double));
/*--------------------------------------------- fill sendbuffers from x */
for (n=0; n<nproc; n++)
{
   if (n==myrank || ilu->computebuff.a.ia[myrank][n]==0) continue;
   ilu->sendbuff.a.da[n][0] = 0.0;
   for (i=0; i<ilu->gdofsend.a.ia[n][0]; i++)
   {
      dof   = ilu->gdofsend.a.ia[n][i+1];
      index = find_index(dof,A->update.a.iv,A->numeq);
      if (index==-1) dserror("Cannot find local dof");
      ilu->sendbuff.a.da[n][i+1] = x[index];
   }
}
/*----------------------------------------------------------- make send */
nrequest = 0;
for (n=0; n<nproc; n++) nrequest += ilu->computebuff.a.ia[myrank][n];
nrequest *= 2;
request   = (MPI_Request*)CCAMALLOC(nrequest*sizeof(MPI_Request));
j = 0;
for (n=0; n<nproc; n++)
{
   if (n==myrank || ilu->computebuff.a.ia[myrank][n]==0) continue;
   MPI_Isend(ilu->gdofsend.a.ia[n],ilu->gdofsend.a.ia[n][0]+1,MPI_INT   ,n,0,actintra->MPI_INTRA_COMM,&(request[j]));
   j++;
   MPI_Isend(ilu->sendbuff.a.da[n],ilu->gdofsend.a.ia[n][0]+1,MPI_DOUBLE,n,1,actintra->MPI_INTRA_COMM,&(request[j]));
   j++;
}
/*------------------------------------- put my own piece of x to xwork */
index = find_index(A->update.a.iv[0],ilu->update.a.iv,ilu->numeq);
if (index==-1) dserror("Cannot find local dof");
/*----------------------------------------- put values from x to xwork */
for (i=0; i<A->numeq; i++)
   xwork[index+i] = x[i];
/*------------------------------------------------------- make receive */
irecv = amdef("irecv",&irecv_a,nproc,10000,"IA");
drecv = amdef("drecv",&drecv_a,10000,1,"DV");

for (n=0; n<nproc; n++)
{
   if (n==myrank || ilu->computebuff.a.ia[n][myrank]==0) continue;
   MPI_Probe(n,0,actintra->MPI_INTRA_COMM,&status);
   MPI_Get_count(&status,MPI_INT,&k);
   if (k>irecv_a.sdim)
   {
      irecv = amredef(&irecv_a,nproc,k,"IA");
      amdel(&drecv_a);
      drecv = amdef("drecv",&drecv_a,k,1,"DV");
   }
   MPI_Recv(irecv[n],k,MPI_INT   ,n,0,actintra->MPI_INTRA_COMM,&status);
   MPI_Recv(drecv   ,k,MPI_DOUBLE,n,1,actintra->MPI_INTRA_COMM,&status);
   /*--------------------------------------------- put values to xwork */
   k--;
   for (i=0; i<k; i++)
   {
      index = find_index(irecv[n][i+1],ilu->update.a.iv,ilu->numeq);
      if (index==-1) dserror("Cannot find dof");
      xwork[index] += drecv[i+1];
   }
}
/*----------------------------------------------- make matvec with asm */
/*-------------------------- carefull, the numbering is fortran style! */
for (i=0; i<ilu->numeq; i++)
{
   sum      = 0.0;
   for (j=ia[i]-1; j<ia[i+1]-1; j++)
      sum += a[j] * xwork[ja[j]-1];
   ywork[i] = sum;
}
/*------------------------------------------ wait for isends to finish */
for (i=0; i<nrequest; i++) 
   MPI_Wait(&(request[i]),&status);
/*----------------------------------------- fill sendbuffer from ywork */
dsend = amdef("dsend",&dsend_a,nproc,irecv_a.sdim,"DA");
for (n=0; n<nproc; n++)
{
   if (n==myrank || ilu->computebuff.a.ia[n][myrank]==0) continue;
   for (i=0; i<irecv[n][0]; i++)
   {
      dof   = irecv[n][i+1];
      index = find_index(dof,ilu->update.a.iv,ilu->numeq);
      if (index==-1) dserror("Cannot find overlap dof in ilu solution");
      dsend[n][i+1] = ywork[index];
   }
}
/*---------------------------------------------------------- make send */
j=0;
for (n=0; n<nproc; n++)
{
   if (n==myrank || ilu->computebuff.a.ia[n][myrank]==0) continue;
   MPI_Isend(dsend[n],irecv[n][0]+1,MPI_DOUBLE,n,10,actintra->MPI_INTRA_COMM,&(request[j]));
   j++;
}
nrequest = j;
/*------------------------------- put my own piece of ywork back to y */
if (init) 
   dveczero(y,&(A->numeq));
/*-------------------------------------------------- find first entry */
index = find_index(A->update.a.iv[0],ilu->update.a.iv,ilu->numeq);
if (index==-1) dserror("Cannot find local dof");
/*---------------------------------------- put values from ywork to y */
for (i=0; i<A->numeq; i++)
   y[i] += ywork[index+i]*fac;
/*--------------------------- receive incoming parts and put to ywork */
for (n=0; n<nproc; n++)
{
   if (n==myrank || ilu->computebuff.a.ia[myrank][n]==0) continue;
   MPI_Probe(n,10,actintra->MPI_INTRA_COMM,&status);
   MPI_Get_count(&status,MPI_DOUBLE,&k);
   if (k-1 != ilu->gdofsend.a.ia[n][0]) dserror("Message size wrong");
   MPI_Recv(ilu->sendbuff.a.da[myrank],k,MPI_DOUBLE,n,10,actintra->MPI_INTRA_COMM,&status);
   /*-------------------------------- put received values to my ywork */
   for (i=0; i<ilu->gdofsend.a.ia[n][0]; i++)
   {
      dof   = ilu->gdofsend.a.ia[n][i+1];
      index = find_index(dof,A->update.a.iv,A->numeq);
      if (index==-1) dserror("Cannot find local dof");
      y[index] += ilu->sendbuff.a.da[myrank][i+1]*fac;
   }
}
/*----------------------------------------- wait for isends to finish */
for (i=0; i<nrequest; i++) 
   MPI_Wait(&(request[i]),&status);
/*----------------------------------------------------------- tidy up */
CCAFREE(request);
CCAFREE(xwork);
CCAFREE(ywork);
amdel(&irecv_a);
amdel(&drecv_a);
amdel(&dsend_a);
#ifdef DEBUG 
dstrc_exit();
#endif
#endif /* end of #ifdef PARALLEL */
return;
}



/*! @} (documentation module close)*/
