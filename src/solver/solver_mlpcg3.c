/*!---------------------------------------------------------------------
\file
\brief contains the multilevel solver for shells

---------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
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
/*!----------------------------------------------------------------------
\brief the multilevel preconditioned solver main structure

<pre>                                                         m.gee 09/02    
defined in solver_mlpcg.c
</pre>

*----------------------------------------------------------------------*/
extern struct _MLSOLVER mlsolver;
/*!---------------------------------------------------------------------
\brief create multilevel solver                                              

<pre>                                                        m.gee 9/02 

</pre>
\param bdcsr    DBCSR*         (i)   the distributed  csr matrix                   
\param sol      DIST_VECTOR*   (i)   the distributed  solution vector         
\param rhs      DIST_VECTOR*   (i)   the distributed  rhs vector         
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_solver_create(DBCSR       *bdcsr, 
                         DIST_VECTOR *sol,
                         DIST_VECTOR *rhs,
                         MLPCGVARS   *mlpcgvars)
{
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_solver_create");
#endif
/*----------------------------------------------------------------------*/
mlsolver.tol = mlpcgvars->tol;
mlsolver.maxiter = mlpcgvars->maxiter;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_solver_create */



/*!---------------------------------------------------------------------
\brief multilevel preconditioned cg                                              

<pre>                                                        m.gee 9/02 

</pre>
\param bdcsr    DBCSR_ROOT*    (i)   the distributed csr matrix                   
\param sol      DIST_VECTOR*   (i)   the distributed  solution vector         
\param rhs      DIST_VECTOR*   (i)   the distributed  rhs vector         
\param actintra INTRA*         (i)   the intra-communicator of this field                  
\return void                                               
\sa mlpcg_solver_init mlpcg_matvec_init mlpcg_matvec
------------------------------------------------------------------------*/
void mlpcg_pcg(DBCSR       *bdcsr, 
               DIST_VECTOR *sol,
               DIST_VECTOR *rhs,
               INTRA       *actintra)
{
INT        i;
INT        myrank,nproc;
INT        numeq;
DOUBLE    *r,*z,*p,*q,*x,*b;
DOUBLE     rho1,rho2=1.0,beta,alfa;
DOUBLE     work;
DOUBLE     eps=1000000.0;
DOUBLE     t1,t2;
DOUBLE     done=1.0;
INT        ione=1;
INT        izero=0;
INT        one=1,two;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_pcg");
#endif
/*----------------------------------------------------------------------*/
myrank = actintra->intra_rank;
nproc  = actintra->intra_nprocs;
/*----------------------------------------------------------------------*/
numeq       = bdcsr->numeq;
/*----------------------------------------------------------------------*/
x = sol->vec.a.dv;
b = rhs->vec.a.dv;
r = mlsolver.r.a.dv;
z = mlsolver.z.a.dv; 
p = mlsolver.p.a.dv;
q = mlsolver.q.a.dv;
/*------------------------------------------------------------take time */
t1 = ds_cputime();
/*----------------------------------------------------- make r = b - Ax */
mlpcg_matvec(r,bdcsr,x,-1.0,1,actintra);
mlpcgupdvec(r,b,&done,&izero,&numeq);
/*----------------------------------- loop maximum number of iterations */
for (i=0; i<mlsolver.maxiter; i++)
{
   /*--------------------------------------------- preconditioner solve */
   amzero(&(mlsolver.z));
   /*----- this is the v- or w-cycle algebraic multigrid preconditioner */
   two=2;
   mlpcg_precond_amgVW(0,&(mlsolver.z),&(mlsolver.r),actintra,&one);
   /*----------- this is the f-cycle algebraic multigrid preconditioner */
   /*two = 2;
   mlpcg_precond_amgF(0,&(mlsolver.z),&(mlsolver.r),actintra,&two);
   /*------- this is a simple additive schwarz one level preconditioner */
   /*mlpcg_precond_presmo(z,r,bdcsr,&(mlprecond.level[0]),actintra,0);
   /*--------------------------------- this is no preconditioner at all */
   /*mlpcgupdvec(z,r,&done,&ione,&numeq);
   /*----------------------------------------- inner product rho1 = r*z */
   mlpcg_vecvec(&rho1,r,z,numeq,actintra);
   /*------------------------------- -------- check for first iteration */
   if (i==0)
      mlpcgupdvec(p,z,&done,&ione,&numeq);
   /*---------------------------------------------- not first iteration */
   else
   {
      /*---------------------------------------------- beta = rho1/rho2 */
      beta = rho1/rho2;
      /*---------------------------------------------- p = z + beta * p */
      mlpcgupdupdvec(p,z,&done,p,&beta,&ione,&numeq);
   } 
   /*----------------------------------------------------------- q = Ap */
   mlpcg_matvec(q,bdcsr,p,1.0,1,actintra);
   /*------------------------------------------------- alfa = rho1 / pq */
   mlpcg_vecvec(&work,p,q,numeq,actintra);
   alfa = rho1/work;
   /*------------------------------------------------- x = x + alfa * p */
   mlpcgupdvec(x,p,&alfa,&izero,&numeq);
   /*------------------------------------------------ write x to output */
   /*mlpcg_printvec(i,x,bdcsr,mlprecond.fielddis,mlprecond.partdis,actintra);
   /*------------------------------------------------- r = r - alfa * q */
   work = -1.0 * alfa;
   mlpcgupdvec(r,q,&work,&izero,&numeq);
   /*------------------------------------------------ check convergence */
   /* make l2-norm of r */
   mlpcg_vecvec(&work,r,r,numeq,actintra);
   /* make eps */
   eps = sqrt(work);
   if (i!=0)
   if (i%10==0)
   if (myrank==0)
   {
      printf("%d : %20.15f\n",i,eps);
      fflush(stdout);
   }
   if (eps<=mlsolver.tol) 
   break;
   /*------------------------------------------------------ update rho2 */
   rho2 = rho1;
   /*---------------------------------------------- goto next iteration */
} /* end of for (i=0; i<mlsolver.maxiter; i++) */
/*------------------------------------------------------------take time */
t2 = ds_cputime();
/*-------------------------------------------- print time for iteration */
if (actintra->intra_rank==0 /*&& mlprecond.ncall==0*/)
printf("Time iteration   : %20.10f\n",t2-t1);
/*----------------------------------------------- check reason for exit */
if ((i==mlsolver.maxiter || eps>mlsolver.tol ) && myrank==0)
{
   printf("MLPCG: No convergence in maxiter steps\n");
   fprintf(allfiles.out_err,"MLPCG: No convergence in maxiter steps\n");
   i--;
}
if (myrank==0)
   printf("MLPCG: numiter: %d   eps: %E\n",i+1,eps);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_pcg */



/*!---------------------------------------------------------------------
\brief init multilevel solver                                              

<pre>                                                        m.gee 9/02 

</pre>
\param bdcsr    DBCSR_ROOT*    (i)   the distributed csr matrix                   
\param sol      DIST_VECTOR*   (i)   the distributed  solution vector         
\param rhs      DIST_VECTOR*   (i)   the distributed  rhs vector         
\param actintra   INTRA*       (i)   the intra-communicator of this field                  
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_solver_init(DBCSR       *bdcsr, 
                       DIST_VECTOR *sol,
                       DIST_VECTOR *rhs,
                       INTRA       *actintra)
{
INT        myrank,nproc;
INT        numeq;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_solver_init");
#endif
/*----------------------------------------------------------------------*/
myrank = actintra->intra_rank;
nproc  = actintra->intra_nprocs;
/*----------------------------------------------------------------------*/
numeq       = bdcsr->numeq;
/*-------------------------------------- allocate the iteration vectors */
if (mlsolver.r.Typ != cca_DV)
{
   amdef("r",&(mlsolver.r),numeq,1,"DV");
   amdef("z",&(mlsolver.z),numeq,1,"DV");
   amdef("p",&(mlsolver.p),numeq,1,"DV");
   amdef("q",&(mlsolver.q),numeq,1,"DV");
}
else if (mlsolver.r.fdim != numeq) 
   dserror("Mismatch in dimensions");
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_solver_init */



/*!---------------------------------------------------------------------
\brief the parallel matrix vector product                                              

<pre>                                                        m.gee 9/02 
the parallel matrix vector product y += fac * A * x
Note that this routine is piece of the mos inner computation, so quite
an effort is taken to make very smooth communication using
incomplete send, testing for multiple messages, overlayering communication
with computation 
</pre>
\param y        DOUBLE*        (o)   the result vector of the matrix vector product
\param bdcsr    DBCSR_ROOT*    (i)   the distributed csr matrix                   
\param x        DOUBLE*        (i)   the input vector of the product
\param fac      DOUBLE         (i)   scaling factor
\param init     INT            (i)   init=1: y = fac*A*x / init=0: y += fac*A*x
\param actintra INTRA*         (i)   the intra-communicator of this field                  
\return void   
\warning the incomplete communication used needs mpi-internal buffer space.
         It may be necessary to allocate extra buffer for MPI, see MPI manual.                                            
\sa mlpcg_matvec_init
------------------------------------------------------------------------*/
void mlpcg_matvec(DOUBLE       *y, 
                  DBCSR        *A,
                  DOUBLE       *x,
                  DOUBLE        fac,
                  INT           init,
                  INTRA        *actintra)
{
INT        i,j,counter;
INT        myrank,nproc;
INT        index;

INT        numeq;
INT        numeq_total;

INT       *ia,*ja,*update;
DOUBLE    *a;

INT       *gdofr,ngdofr;
DOUBLE   **rbuff,*cbuff;

INT      **gdofs;
DOUBLE   **sbuff;

INT        recv_size[MAXPROC],nrecv;
INT        own,flag;

INT        fcd,fcdindex;
INT        column;

INT        shift;
INT        bins[5000];

#ifdef PARALLEL 
MPI_Status local_status;
#endif
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_matvec");
#endif
/*----------------------------------------------------------------------*/
myrank = actintra->intra_rank;
nproc  = actintra->intra_nprocs;
/*----------------------------------------------------------------------*/
mlpcg_matvec_init(A,actintra);
if (A->sendbuff.Typ != cca_DA) dserror("Matrix has not been initialized");
/*----------------------------------------------------------------------*/
numeq       = A->numeq;
numeq_total = A->numeq_total;
/*----------------------------------------------------------------------*/
ia          = A->ia.a.iv;
ja          = A->ja.a.iv;
a           = A->a.a.dv;
update      = A->update.a.iv;
gdofr       = A->gdofrecv.a.iv;
ngdofr      = A->gdofrecv.fdim;
rbuff       = A->recvbuff.a.da;
gdofs       = A->gdofsend.a.ia;
sbuff       = A->sendbuff.a.da;
cbuff       = A->computebuff.a.dv;
fcd         = A->firstcoupledof;
fcdindex    = mlpcg_getindex(fcd,update,numeq);
if (fcdindex==-1) dserror("Cannot find local dof in update");
/*----------------------------------------------------------------------*/
/* 
we will shade the communication of the ghost dofs by the computation of
the processor local parts of the product. So we first make incomplete
communication
*/
/*----------------------------------------------------------------------*/
/*------------------------ fill the sendbuffer and make incomplete send */
/*
t1 = ds_cputime();
*/
for (i=0; i<nproc; i++)
{
   if (i==myrank)      continue;/* do not send to myself */
   if (gdofs[i][0]==0) continue;/* do not send messenges of length zero */
   for (j=0; j<gdofs[i][0]; j++)
   {
      index = mlpcg_getindex(gdofs[i][j+1],update,numeq);
      if (index==-1) dserror("local dof not found on proc");
      dsassert(0<=i<A->sendbuff.fdim,"buffer overflow");
      dsassert(0<=j<A->sendbuff.sdim,"buffer overflow");
      sbuff[i][j] = x[index];
   }
   dsassert(j==gdofs[i][0],"size mismatch");
#ifdef PARALLEL 
   /* send the approbiate values of x */          
   MPI_Isend(&(sbuff[i][0]),gdofs[i][0],MPI_DOUBLE,i,myrank,
             actintra->MPI_INTRA_COMM,&(A->request[i]));
#endif
}
/*
t2 = ds_cputime();
if (myrank==0) printf("matvec INTER  I: %20.10f\n",t2-t1);
*/
/*---------------------------------------------init the solution vector */
if (init)
   mlpcgveczero(y,&numeq);
/*----------------------------- make the local piece of the computation */
/*
t1 = ds_cputime();
*/
/* loop indizes of rows, which are strictly local  */
for (i=0; i<fcdindex; i++)
   /* loop column indizes */
   for (j=ia[i]; j<ia[i+1]; j++)
   {
      /* get global column dof number */
      index  = mlpcg_getindex(ja[j],update,numeq);
      dsassert(index != -1,"Cannot find local dof in update");
      /*================================================================*/
      dsassert(0<=i<numeq,"vector overflow");
      y[i] += a[j]*x[index]*fac;
      /*================================================================*/
   }
/* loop indizes of rows, which also contain interproc off-diagonals */
for (i=fcdindex; i<numeq; i++)
   /* loop column indizes */
   for (j=ia[i]; j<ia[i+1]; j++)
   {
      /* get global column dof number */
      index  = mlpcg_getindex(ja[j],update,numeq);
      /* if the column if not proc local do nothing */
      if (index==-1) continue;
      /*================================================================*/
      dsassert(0<=i<numeq,"vector overflow");
      y[i] += a[j]*x[index]*fac;
      /*================================================================*/
   }
/*
t2 = ds_cputime();
if (myrank==0) printf("matvec LOCAL   : %20.10f\n",t2-t1);
*/
/*---------------------------------------------------- make the receive */
/*-------------- check how many and from who I have to receive messages */
/*
t1 = ds_cputime();
*/
for (i=0; i<nproc; i++) recv_size[i]=0;
for (i=0; i<ngdofr; i++)
{
   own = mlpcg_getowner(gdofr[i],A->owner,nproc);
   recv_size[own]++;   
}
/*---------------------- especially receive in the order as they arrive */
startrecv:
flag=0;
for (i=0; i<nproc; i++)
{
   /* do not receive messages of length zero */
   if (recv_size[i]==0) continue;
#ifdef PARALLEL 
   /* check for message from proc i */
   MPI_Iprobe(i,i,actintra->MPI_INTRA_COMM,&flag,&(A->status[i]));
   if (flag==0) continue; /* there is no message from proc i yet */
   else                   /* message from proc i has arrived */
   {
      /* check the sender */
      if (A->status[i].MPI_SOURCE != i || A->status[i].MPI_TAG != i)
      dserror("Messages got mixed up");
      /* check the size of the message */
      MPI_Get_count(&(A->status[i]),MPI_DOUBLE,&nrecv);
      if (nrecv != recv_size[i]) dserror("Size of message is not what was expected");
      /* set recv_size[i] to zero, to indicate, that this message from i is finished */
      recv_size[i]=0;
      /* receive the message */
      dsassert(nrecv<=A->recvbuff.sdim,"Buffer overflow");
      MPI_Recv(&(rbuff[i][0]),nrecv,MPI_DOUBLE,i,i,
               actintra->MPI_INTRA_COMM,&(A->status[i]));
      /* put the values to the computebuff vector */
      counter=0;
      for (j=0; j<ngdofr; j++)
      {
         own = mlpcg_getowner(gdofr[j],A->owner,nproc);
         if (own != i) continue;
         else
         {
            dsassert(0<=j<A->computebuff.fdim,"Buffer overflow");
            dsassert(counter<A->recvbuff.sdim,"Buffer overflow");
            cbuff[j] = rbuff[i][counter];
            counter++;
         }
      }   
   }
#endif
}
/*----------------------- check whether all messages have been received */
counter=0;
for (i=0; i<nproc; i++) counter += recv_size[i];
if (counter != 0) goto startrecv;
/*---------------- make computation of the interproc part of the matrix */
/* loop the indizes of the local rows */
if (nproc>1)
{
if (ngdofr>18000) dserror("static bins too small");
init_quick_find(gdofr,ngdofr,&shift,bins);
for (i=fcdindex; i<numeq; i++) 
   /* loop indizes of colmns */
   for (j=ia[i]; j<ia[i+1]; j++)
   {
      /* get global column number */
      column = ja[j];
      /* check for ownership */
      own    = mlpcg_getowner(column,A->owner,nproc);
      /* if it is a non-interproc offdiagonal or a diagonal entry, it has 
                                                  been computed already */
      if (own==myrank) continue;
      /* it is an interproc column */
      index = quick_find(column,gdofr,ngdofr,shift,bins);
      if (index==-1) dserror("Cannot find interproc coupling column");
      /*========================== yeah! do the interproc computation ! */
      dsassert(0<=i<numeq,"buffer overflow");
      dsassert(0<=index<A->computebuff.fdim,"buffer overflow");
      y[i] += a[j] * cbuff[index] * fac;
      /*================================================================*/
   }   
}
/*--------------------------------- close the incomplete communications */
#ifdef PARALLEL 
for (i=0; i<nproc; i++)
if (A->request[i])
MPI_Wait(&(A->request[i]),&local_status);
#endif
/*
t2 = ds_cputime();
if (myrank==0) printf("matvec INTER II: %20.10f\n",t2-t1);
*/
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_matvec */



/*!---------------------------------------------------------------------
\brief init the matrix vector product of a dbcsr matrix                                              

<pre>                                                        m.gee 9/02 
This routine checks the ghost dofs of the distributed csr matrix, means
it looks for the dofs on other procs, which are locally necessary to
perform a multiplication with this piece of csr matrix. Plausibility tests
of the matrix are made and it is checked whether the given matrix has been 
initialized by this routine before

   For simple computation we demand the following properties of the given
   bdcsr matrix:
   -1.) dof numbering and indexing is C-style (starts with 0)
   0.) the bdcsr matrix is partitioned in a row-wise manner
   1.) the dofs on one proc must be continous and sorted in ascending order
   2.) all dofs on proc-1 must be lower then all dofs on proc
   3.) rows which hold interprocessor coupling off-diagonal entries
       are sorted to the end of the local piece of the matrix
       -> bdcsr->firstcoupledof refers to the first row with interprocessor
                                coupling
       -> no row before bdcsr->firstcoupledof has interproc coupling

</pre>
\param bdcsr      DBCSR*       (i/o) the distributed csr matrix                   
\param actintra   INTRA*       (i)   the intra-communicator of this field                  
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_matvec_init(DBCSR       *bdcsr, 
                       INTRA       *actintra)
{
INT        i,j,k,counter;
INT        myrank,nproc;
INT        numeq,numeq_total;
INT       *update;
INT       *ja;
INT       *ia;
INT       *gdof,**gdofsend;
INT        fcd,index_fcd;
INT        start_index,end_index;
INT        actdof,index_actdof;
INT        ngdof;
ARRAY      gupdatesend_a,gupdaterecv_a;
INT       *gupdatesend,*gupdaterecv;
INT        column;
INT        column_owner;
INT        index;
INT        rbuffsize[MAXPROC];
INT        own;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_matvec_init");
#endif
/*----------------------------------------------------------------------*/
/*--------------------------------- this matrix has already been inited */
if (bdcsr->gdofrecv.Typ==cca_IV) goto end;
/*----------------------------------------------------------------------*/
myrank = actintra->intra_rank;
nproc  = actintra->intra_nprocs;
/*----------------------------------------------------------------------*/
numeq        = bdcsr->numeq;
numeq_total  = bdcsr->numeq_total;
update       = bdcsr->update.a.iv;
ja           = bdcsr->ja.a.iv;
ia           = bdcsr->ia.a.iv;
/*----------------------------------------------------------------------*/
/* 
   make some plausibility checks for the given bdcsr matrix.
   
   For efficient computation we demand the following properties of the given
   bdcsr matrix:
   -1.) dof numbering and indexing is C-style (starts with 0)
   0.) the bdcsr matrix is partitioned in a rowwise manner
   1.) the dofs on one proc must be continous and sorted in ascending order
   2.) all dofs on proc-1 must be lower then all dofs on proc
   3.) rows which hold interprocessor coupling off-diagonal entries
       are sorted to the end of the local piece of the matrix
       -> bdcsr->firstcoupledof refers to the first row with interprocessor coupling
       -> no row before bdcsr->firstcoupledof has interproc coupling
   4.) no row is allowed to be laeger then the define MAX_NNZPERROW
*/
/* check dof 0 is on proc 0 */
if (myrank==0) {if (update[0] != 0) dserror("Distribution of BDCSR matrix not ascending");}
else           {if (update[0] == 0) dserror("Distribution of BDCSR matrix not ascending");}

/* check dofs are continous on proc */
for (i=0; i<numeq-1; i++)
if (update[i+1] != update[i]+1) dserror("BDCSR matrix is not continous");
/* check interproc coupled equations are at the end of matrix */
fcd       = bdcsr->firstcoupledof;
index_fcd = mlpcg_getindex(fcd,update,numeq);
if (index_fcd==-1) dserror("cannot find dof in update");
for (i=0; i<ia[index_fcd]; i++)
{
   j = mlpcg_getindex(ja[i],update,numeq);
   if (j==-1) 
      dserror("Interprocessor coupled rows are not sorted correctly");
}
/* check sizes of rows */
for (i=0; i<numeq; i++)
   if (ia[i+1]-ia[i] >= MAX_NNZPERROW)
      dserror("Define of MAX_NNZPERROW is too small");
/*------- make number of external values needed by this piece of matrix */
/*--------------------------------------- make list of the "ghost dofs" */
gdof = amdef("tmp",&(bdcsr->gdofrecv),10000,1,"IV");
/*--------- loop the coupled dofs and check off-diagonal entries of them */
/*------------------------- start in the row of the fist coupled dof fcd */
start_index = ia[index_fcd];
/*--------------------------------------- end at the end of the last row */
end_index   = ia[numeq];
/*--------------------- count number of dofs, which are not in this proc */
ngdof = 0;
for (i=start_index; i<end_index; i++)
{
   /* active dof column */
   actdof = ja[i];
   /* check for active dof column on this proc */
   index_actdof = mlpcg_getindex(actdof,update,numeq);
   /* if it is on this proc, forget it */
   if (index_actdof != -1) continue;
   /* check whether active dof has been counted before */
   for (j=0; j<ngdof; j++)
   if (actdof==gdof[j]) 
   goto nextdof;
   /* it's a previously unknown ghost-dof */
   gdof[ngdof]=actdof;
   ngdof++;
   if (ngdof==bdcsr->gdofrecv.fdim) 
   amredef(&(bdcsr->gdofrecv),bdcsr->gdofrecv.fdim+2000,1,"IV");
   nextdof:;
}
qsort((INT*)gdof,ngdof,sizeof(INT),cmp_int);
gdof = (INT*)amredef(&(bdcsr->gdofrecv),ngdof,1,"IV");
amdef("cbuff",&(bdcsr->computebuff),ngdof,1,"DV");
/*----------------------------------------------------------------------*/
/* we have to find the owners of the ghost dofs to store them in the 
   second column of gdof */
/*----------------------------------------------------------------------*/
gupdatesend = amdef("tmp",&gupdatesend_a,numeq_total,1,"IV");
#ifdef PARALLEL
gupdaterecv = amdef("tmp",&gupdaterecv_a,numeq_total,1,"IV");
#endif
amzero(&gupdatesend_a);
for (i=0; i<numeq; i++) 
   gupdatesend[update[i]] = myrank;
#ifdef PARALLEL
MPI_Allreduce(gupdatesend,gupdaterecv,numeq_total,MPI_INT,MPI_SUM,actintra->MPI_INTRA_COMM);
#else
gupdaterecv = gupdatesend;
#endif
/*------------------------------- make the ownership array bdcsr->owner */
counter=0;
if (gupdaterecv[0] != 0) dserror("ownership of procs mixed up");
bdcsr->owner[0][0] = 0;
for (i=0; i<numeq_total; i++)
{
   if (gupdaterecv[i] != counter)
   {
      bdcsr->owner[counter][1] = i-1;
      counter++;
      bdcsr->owner[counter][0] = i;
   }
}
bdcsr->owner[counter][1] = numeq_total-1;
dsassert(counter+1==nproc,"number of processors wrong");
/*----------------------------------------------------------------------*/
/* now I(the proc) know in ascending order from which proc I will receive which
   ghost dof. We now have to find out to who we have to send our own dofs. 
   Due to the symmetry of the problem, my coupled dofs are ghost to all 
   which are ghost to me.
   But in the send direction my dof can be ghost to several other procs, 
   so we need this information more detailed
*/
/*----------------------------------------------------------------------*/
gdofsend = (INT**)amdef("gdofsnd",&(bdcsr->gdofsend),nproc,500+1,"IA");
amzero(&(bdcsr->gdofsend));
/*------------------------------------------- loop my coupled equations */
for (i=index_fcd; i<numeq; i++)
{
   /* set the coupled row */
   actdof = update[i];
   /* set columns index range */
   start_index = ia[i];
   end_index   = ia[i+1];
   /* loop the columns in this row */
   for (j=start_index; j<end_index; j++)
   {
      /* the active column */
      column = ja[j];
      /* active column has to belong to ghost dofs */
      index = -1;
      for (k=0; k<ngdof; k++)
      if (gdof[k]==column)
      {
         index = k;
         break;
      }
      /* if it does not belong to ghost dofs continue */
      if (index == -1) continue;
      /* find the owner of this ghost dof */
      column_owner = gupdaterecv[column];
/*      column_owner = mlpcg_getowner(column,bdcsr->owner,nproc);*/
      gdofsend[column_owner][ gdofsend[column_owner][0]+1 ] = actdof;
      gdofsend[column_owner][0]++;
      if (gdofsend[column_owner][0]+1 == bdcsr->gdofsend.sdim)
      gdofsend = amredef(&(bdcsr->gdofsend),nproc,bdcsr->gdofsend.sdim+1000,"IA");
   }
}
/*----------------------- delete the doubles to be send to the same proc */
for (i=0; i<nproc; i++)
{
   if (i == myrank) continue;
   for (j=0; j<gdofsend[i][0]; j++)
   {
      actdof = gdofsend[i][j+1];
      if (actdof == -1) continue;
      for (k=j+2; k<gdofsend[i][0]+1; k++)
      if (gdofsend[i][k] == actdof)
         gdofsend[i][k]=-1;
   }
   counter=0;
   for (j=0; j<gdofsend[i][0]; j++)
   {
      if (gdofsend[i][j+1] != -1)
      {
         gdofsend[i][counter+1] = gdofsend[i][j+1];
         counter++;
      }
   }
   gdofsend[i][0]=counter;
}
/*---------- find the largest number in gdofsend and redefine the array */
j=0;
for (i=0; i<nproc; i++)
{
   if (i==myrank) continue;
   if (gdofsend[i][0] > j) j = gdofsend[i][0];
}
gdofsend = amredef(&(bdcsr->gdofsend),nproc,j+1,"IA");
/*---------------------------------------- sort them in ascending order */
for (i=0; i<nproc; i++)
{
   if (i==myrank) continue;
   qsort((INT*)&(gdofsend[i][1]),gdofsend[i][0],sizeof(INT),cmp_int);
}
/*----------------------------------------------------------------------*/
/* 
   now we have gdofrecv and gdofsend, we can now allocate appropriate 
   DOUBLE precision send- and recvbuffers for the communication
*/
/*---------------------------- make recvbuff size and allocate recvbuff */   
counter=0;
for (i=0; i<MAXPROC; i++) 
rbuffsize[i] = 0;
for (i=0; i<bdcsr->gdofrecv.fdim; i++)
{
   actdof = gdof[i];
   own = gupdaterecv[actdof];
   if (own==myrank) continue;
   rbuffsize[own]++;
}
for (i=0; i<nproc; i++)
   if (rbuffsize[i]>counter) 
      counter = rbuffsize[i];
amdef("rbuff",&(bdcsr->recvbuff),nproc,counter,"DA");
/*--------------------------------------------- make size of sendbuffer */      
counter=0;
for (i=0; i<nproc; i++)
   if (gdofsend[i][0]>counter)
      counter = gdofsend[i][0];
amdef("sbuff",&(bdcsr->sendbuff),nproc,counter,"DA");      
/*----------------------------------------- allocate request and status */
#ifdef PARALLEL
bdcsr->status  = (MPI_Status*)CCACALLOC(nproc,sizeof(MPI_Status));
bdcsr->request = (MPI_Request*)CCACALLOC(nproc,sizeof(MPI_Request));
#endif
/*----------------------------------------------------------------------*/
amdel(&gupdatesend_a);
#ifdef PARALLEL
amdel(&gupdaterecv_a);
#endif
/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_matvec_init */



/*!---------------------------------------------------------------------
\brief uninit the matrix vector product of a dbcsr matrix                                              

<pre>                                                        m.gee 9/02 

</pre>
\param bdcsr    DBCSR*    (i/o)   the distributed csr matrix                   
\warning After a call to this routine, no matrix vector product with 
         the DBCSR matrix can be done
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_matvec_uninit(DBCSR       *bdcsr)
{
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_matvec_uninit");
#endif
/*----------------------------------------------------------------------*/
/* check whether this matrix is inited */
if (bdcsr->gdofrecv.Typ!=cca_IA) goto end;
amdel(&(bdcsr->gdofrecv));
amdel(&(bdcsr->recvbuff));
amdel(&(bdcsr->computebuff));
amdel(&(bdcsr->gdofsend));
amdel(&(bdcsr->sendbuff));
#ifdef PARALLEL
CCAFREE(bdcsr->status);
CCAFREE(bdcsr->request);
#endif
/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_matvec_uninit */



/*!---------------------------------------------------------------------
\brief make the vector product scalar = x * y with distributed vectors                                              

<pre>                                                        m.gee 9/02 
make the vector product scalar = x * y with distributed vectors 
</pre>
\param scalar   DOUBLE*    (o)        the result
\param x        DOUBLE*    (i)        the vector to update y                  
\param y        DOUBLE*    (i)        the updated vector length numeq                   
\param dim      const INT  (i)        the dimension of vectors y and x                   
\param actintra INTRA*     (i)        the corresponding intra-communicator
\return void                                         

------------------------------------------------------------------------*/
void mlpcg_vecvec(DOUBLE *scalar, DOUBLE *x, DOUBLE *y, INT dim, INTRA *actintra)
{
DOUBLE  sbuff=0.0;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_vecvec");
#endif
/*----------------------------------------------------------------------*/
mlpcgvecvec(x,y,&sbuff,&dim);
#ifdef PARALLEL
MPI_Allreduce(&sbuff,scalar,1,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
#else
*scalar = sbuff;
#endif
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_vecvec */









/*!---------------------------------------------------------------------
\brief return index of a given dof from local dof list                                              

<pre>                                                        m.gee 9/02 
return index of a given dof from local dof list. The given vector update
of length length has to be sorted and continous
</pre>
\param dof      INT    (i)           the dof the index is needed for                   
\param update   INT*   (i)           the sorted and continous vector update                   
\param length   INT    (i)           length of update                   
\return the index of dof in update (INT) or -1 if index not in range                                              
\sa mlpcg_matvec_init  mlpcg_getowner                                       

------------------------------------------------------------------------*/
INT mlpcg_getindex(INT dof, INT *update, INT length)
{
INT index;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_getindex");
#endif
/*----------------------------------------------------------------------*/
index = dof-(*update);
if (index<0 || index>=length) 
{
#ifdef DEBUG 
dstrc_exit();
#endif
   return(-1);
}
else
{
#ifdef DEBUG 
dstrc_exit();
#endif
   return(index);
}
/*----------------------------------------------------------------------*/
} /* end of mlpcg_getindex */



/*!---------------------------------------------------------------------
\brief return the owner of a given dof                                              

<pre>                                                        m.gee 9/02 
return the proc a given dof of a BDCSR matrix is held on. This only works
after a function call to mlpcg_matvec_init
</pre>
\param dof      INT      (i)           the dof the index is needed for                   
\param owner    INT[][2] (i)           the owner array of the bdcsr matrix                   
\param nproc    INT      (i)           number of processors                   
\return the owner of a given dof    
\sa mlpcg_matvec_init mlpcg_getindex                                        

------------------------------------------------------------------------*/
INT mlpcg_getowner(INT dof, INT owner[][2], INT nproc)
{
INT i,own;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_getowner");
#endif
/*----------------------------------------------------------------------*/
own=-1;
for (i=0; i<nproc; i++)
{
   if (owner[i][0] <= dof && owner[i][1] >= dof)
   {
      own=i;
      break;
   }
}
if (own==-1) dserror("Cannot find owner of a given dof");
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return(own);
} /* end of mlpcg_getowner */






/*! @} (documentation module close)*/
