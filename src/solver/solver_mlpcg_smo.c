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
\brief Jacobi smoother                                              

<pre>                                                        m.gee 10/02 

</pre>
\param z            double*      (o)   the solution of the smoothing
\param r            double*      (i)   the right hand side
\param csr          DBCSR*       (i)   the matrix to smooth with
\param nsweep       int          (i)   number of smoothing cycles
\param actintra     INTRA*       (i)   the intra-communicator of this field                  
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_precond_smoJacobi(double *z, double *r, DBCSR *csr, int nsweep, INTRA *actintra)
{
int     i,n;
int     numeq;
ARRAY   Dinv_a;
double *Dinv;
ARRAY   work_a;
double *work;
double  done=1.0;
int     ione=1;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_smoJacobi");
#endif
/*----------------------------------------------------------------------*/
if (nsweep==0) goto exit;
/*----------------------------------------------------------------------*/
numeq = csr->numeq;
/*-------------- allocate vector for the inverse of the diagonal of csr */
Dinv = amdef("Dinv",&Dinv_a,numeq,1,"DV");
/*------------------------------------------- allocate a working vector */
work = amdef("work",&work_a,numeq,1,"DV");
/*----------------------- get the inverse of the diagonal of csr matrix */
mlpcg_csr_getdinv(Dinv,csr,numeq);
/*--------------------------------------------------- make z = Dinv * r */
for (i=0; i<numeq; i++) 
   z[i] = Dinv[i] * r[i];
/*--------------------------------------- copy r to working vector work */
mlpcgupdvec(work,r,&done,&ione,&numeq);
/*----------------------------------------- loop about number of sweeps */
for (n=1; n<nsweep; n++)
{
   /* make work = work - A * z */
   mlpcg_matvec(work,csr,z,-1.0,0,actintra);
   /* make z = z + Dinv * work */
   for (i=0; i<numeq; i++)
      z[i] += Dinv[i] * work[i];
   /* copy r to work */
   mlpcgupdvec(work,r,&done,&ione,&numeq);
}
/*------------------------------------------------------------- tidy up */
amdel(&Dinv_a);
amdel(&work_a);
/*----------------------------------------------------------------------*/
exit:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_smoJacobi */


/*!---------------------------------------------------------------------
\brief ilu(n) smoother                                              

<pre>                                                        m.gee 11/02 

</pre>
\param z            double*      (o)   the solution of the smoothing
\param r            double*      (i)   the right hand side
\param csr          DBCSR*       (i)   the matrix to smooth with
\param nsweep       int          (i)   n in ilu(n)
\param actintra     INTRA*       (i)   the intra-communicator of this field                  
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_precond_smo_ILUn(double *z, double *r, DBCSR *csr, int nsweep, INTRA *actintra)
{
int     i,n;
int     myrank,nproc;
int     numeq;
DBCSR  *ilu;
DBCSR  *asm;
int     size,ierr=0;
ARRAY   levs,w,jw;
int    *update,nupdate;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_smo_ILUn");
#endif
/*----------------------------------------------------------------------*/
myrank = actintra->intra_rank;
nproc  = actintra->intra_nprocs;
/*-------------- do the decomposition of the matrix, if not done before */
if (csr->ilu==NULL)
{
   csr->ilu = (DBCSR*)CCACALLOC(1,sizeof(DBCSR));
   csr->asm = (DBCSR*)CCACALLOC(1,sizeof(DBCSR));
   asm        = csr->asm;
   asm->numeq = csr->numeq;
   ilu        = csr->ilu;
   ilu->numeq = csr->numeq;
   /* extract the quadratic matrix on this processor (Additive Schwartz) */
   /*---------------------- the guess of the size is SURELY large enough */
   mlpcg_csr_open(asm,
                  csr->owner[myrank][0],
                  csr->owner[myrank][1],
                  csr->numeq_total,
                  csr->a.fdim,
                  actintra);
   /*--------------------------------- extract the local block from csr */
   mlpcg_csr_extractsubblock(csr,
                             asm,
                             csr->owner[myrank][0],
                             csr->owner[myrank][1],
                             csr->owner[myrank][0],
                             csr->owner[myrank][1],
                             actintra);
   /*--------------------------------------- close the extracted matrix */
   mlpcg_csr_close(asm);
   /*---------------- change the enumeration to local and fortran style */
   mlpcg_csr_localnumsf(asm);
   /*---------------------- allocate space for the decomposition in ilu */
   if (nsweep==0) size = (asm->a.fdim+1);
   if (nsweep==1) size = (asm->a.fdim+1)*2.0;
   if (nsweep==2) size = (asm->a.fdim+1)*2.5;
   if (nsweep==3) size = (asm->a.fdim+1)*3.5;
   if (nsweep==4) size = (asm->a.fdim+1)*4.5;
   if (nsweep==5) size = (asm->a.fdim+1)*5.5;
   if (nsweep==6) size = (asm->a.fdim+1)*6.5;
   if (nsweep>=7) size = (asm->a.fdim+1)*7.5;
   tryagain:
   amdef("ilu_val"  ,&(ilu->a) ,size          ,1,"DV");
   amdef("ilu_bindx",&(ilu->ja),size          ,1,"IV");
   amdef("ilu_ia"   ,&(ilu->ia),asm->numeq    ,1,"IV");
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
         size *= 1.3;
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
   mlpcg_csr_destroy(asm);
   csr->asm   = CCAFREE(asm);
}
else if (csr->ilu->is_factored != mlprecond.ncall && mlprecond.mod==0)
{
   ilu        = csr->ilu;
   csr->asm   = (DBCSR*)CCACALLOC(1,sizeof(DBCSR));
   asm        = csr->asm;
   mlpcg_csr_open(asm,
                  csr->owner[myrank][0],
                  csr->owner[myrank][1],
                  csr->numeq_total,
                  csr->a.fdim,
                  actintra);
   /*--------------------------------- extract the local block from csr */
   mlpcg_csr_extractsubblock(csr,
                             asm,
                             csr->owner[myrank][0],
                             csr->owner[myrank][1],
                             csr->owner[myrank][0],
                             csr->owner[myrank][1],
                             actintra);
   /*--------------------------------------- close the extracted matrix */
   mlpcg_csr_close(asm);
   /*---------------- change the enumeration to local and fortran style */
   mlpcg_csr_localnumsf(asm);
   /*---------------------- allocate space for the decomposition in ilu */
   if (nsweep==0) size = (asm->a.fdim+1);
   if (nsweep==1) size = (asm->a.fdim+1)*2.0;
   if (nsweep==2) size = (asm->a.fdim+1)*2.5;
   if (nsweep==3) size = (asm->a.fdim+1)*3.5;
   if (nsweep==4) size = (asm->a.fdim+1)*4.5;
   if (nsweep==5) size = (asm->a.fdim+1)*5.5;
   if (nsweep==6) size = (asm->a.fdim+1)*6.5;
   if (nsweep>=7) size = (asm->a.fdim+1)*7.5;
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
         size *= 1.3;
         amdel(&(ilu->a) );
         amdel(&(ilu->ja));
         amdel(&(ilu->ia));
         amdef("ilu_val"  ,&(ilu->a) ,size          ,1,"DV");
         amdef("ilu_bindx",&(ilu->ja),size          ,1,"IV");
         amdef("ilu_ia"   ,&(ilu->ia),asm->numeq    ,1,"IV");
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
   mlpcg_csr_destroy(asm);
   csr->asm   = CCAFREE(asm);
}
/*-------------------------------------- make solve with the ilu matrix */
ilu = csr->ilu;
lusol(&(ilu->numeq),r,z,ilu->a.a.dv,ilu->ja.a.iv,ilu->ia.a.iv);
/*----------------------------------------------------------------------*/
exit:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_smo_ILUn */




/*!---------------------------------------------------------------------
\brief lapack solver                                              

<pre>                                                        m.gee 11/02 

</pre>
\param z            double*      (o)   the solution of the solve 
\param r            double*      (i)   the rhs 
\param csr          DBCSR*       (i)   the matrix to be solved with 
\param actintra     INTRA*       (i)   the intra-communicator of this field                  
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_precond_lapacksolve(double *z, double *r, DBCSR *csr, INTRA *actintra)
{
int      i,j,counter,info;
int      actrow,actcol,colstart,colend,index;
int      myrank,nproc;
int      numeq_total,numeq;
int     *update,*ia,*ja;
double  *a;
ARRAY    sdense_a,srhs_a,rrhs_a;
double **sdense,*srhs,*rrhs;
char     trans[1];
int      ione=1;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_lapacksolve");
#endif
/*----------------------------------------------------------------------*/
myrank      = actintra->intra_rank;
nproc       = actintra->intra_nprocs;
numeq_total = csr->numeq_total;
numeq       = csr->numeq;
update      = csr->update.a.iv;
ia          = csr->ia.a.iv;
ja          = csr->ja.a.iv;
a           = csr->a.a.dv;
if (csr->dense==NULL)
{
   /*- allocate temporary dense matrix and right hand side and solution */
   csr->dense = (ARRAY*)CCACALLOC(1,sizeof(ARRAY));
   csr->ipiv  = (ARRAY*)CCACALLOC(1,sizeof(ARRAY));
   amdef("dense",csr->dense,numeq_total,numeq_total,"DA");
   amdef("ipiv",csr->ipiv,numeq_total,1,"IV");

   sdense = amdef("tmp",&sdense_a,numeq_total,numeq_total,"DA");
            amzero(&sdense_a);
#ifndef PARALLEL 
   amzero(csr->dense);
#endif
   /*------------------------------------------------------ fill sdense */
   for (i=0; i<numeq; i++)
      for (j=ia[i]; j<ia[i+1]; j++)
         {
#ifdef PARALLEL 
            sdense[update[i]][ja[j]] = a[j];
#else
            csr->dense->a.da[update[i]][ja[j]] = a[j];
#endif
         }
   /*--------------------------------------------------- allreduce them */
#ifdef PARALLEL 
   MPI_Allreduce(sdense[0],csr->dense->a.da[0],sdense_a.fdim*sdense_a.sdim,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
#endif
   amdel(&sdense_a);
   /*------------------------------------------- now factor using lapack */
   info = 1;
   dgetrf(&numeq_total,
          &numeq_total,
          csr->dense->a.da[0],
          &numeq_total,
          csr->ipiv->a.iv,
          &info);
   if (info != 0) dserror("Lapack dgetrf returned info nonzero");
   csr->is_factored = mlprecond.ncall;
}/* end of if (csr->dense==NULL) */
else if (csr->is_factored != mlprecond.ncall)
{
   sdense = amdef("tmp",&sdense_a,numeq_total,numeq_total,"DA");
            amzero(&sdense_a);
#ifndef PARALLEL 
   amzero(csr->dense);
#endif
   /*------------------------------------------------------ fill sdense */
   for (i=0; i<numeq; i++)
      for (j=ia[i]; j<ia[i+1]; j++)
         {
#ifdef PARALLEL 
            sdense[update[i]][ja[j]] = a[j];
#else
            csr->dense->a.da[update[i]][ja[j]] = a[j];
#endif
         }
   /*--------------------------------------------------- allreduce them */
#ifdef PARALLEL 
   MPI_Allreduce(sdense[0],csr->dense->a.da[0],sdense_a.fdim*sdense_a.sdim,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
#endif
   amdel(&sdense_a);
   /*------------------------------------------- now factor using lapack */
   info = 1;
   dgetrf(&numeq_total,
          &numeq_total,
          csr->dense->a.da[0],
          &numeq_total,
          csr->ipiv->a.iv,
          &info);
   if (info != 0) dserror("Lapack dgetrf returned info nonzero");
   csr->is_factored = mlprecond.ncall;
}/* end of if (csr->is_factored != mlprecond.ncall) */
/*--------------------------------------------------- make rhs redundant */
srhs   = amdef("tmp",&srhs_a,numeq_total,1,"DV");
         amzero(&srhs_a);
rrhs   = amdef("tmp",&rrhs_a,numeq_total,1,"DV");
for (i=0; i<numeq; i++)
#ifdef PARALLEL 
   srhs[update[i]] = r[i];
#else
   rrhs[update[i]] = r[i];
#endif

#ifdef PARALLEL 
MPI_Allreduce(srhs     ,rrhs     ,srhs_a.fdim                ,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
#endif
amdel(&srhs_a);
/*---------------------------------------- do forward and backward solve */
info     = 1;
trans[0] = 'N';
#ifndef AZTEC_PACKAGE
dgetrs(trans,
       &numeq_total,
       &ione,
       csr->dense->a.da[0],
       &numeq_total,
       csr->ipiv->a.iv,
       rrhs,
       &numeq_total,
       &info);
if (info != 0) dserror("Lapack dgetrf returned info nonzero");
#else
   dserror("solver Lapack conflicts with compilation with -DAZTEC_PACKAGE");
#endif
/*------------------ fill the solution back to the distributed vector z */
for (i=0; i<numeq; i++)
   z[i]   = rrhs[update[i]];
/*------------------------------------------------------------- tidy up */
amdel(&rrhs_a);
/*----------------------------------------------------------------------*/
exit:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_lapacksolve */

/*!---------------------------------------------------------------------
\brief Modified Gram-Schmidt Orthonormalization                                              

<pre>                                                        m.gee 11/02 

</pre>
\param P            double**      (i/o)   The prolongator to be orthonormalized
\param R            double**      (o)     the R part of the P = QR factorization
\param nrow         int           (i)     row dimension of P 
\param ncol         int           (i)     column dimension of P
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_precond_gramschmidt(double **P, double **R,const int nrow,const int ncol)
{
int     i,j,k;
double  sum;
int     start=0;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_gramschmidt");
#endif
/*----------------------------------------------------------------------*/
dsassert(nrow<=500,"Local array Q too small");
/*------------------------------------------- make norm of first column */
tryagain:
sum = 0.0;
for (i=0; i<nrow; i++) 
   sum += P[i][start]*P[i][start];
if (FABS(sum) < EPS13)
{ 
   R[start][start]=0.0;
   start++;
   if (start==ncol) dserror("Prolongator is completely zero");
   goto tryagain;
}
sum = sqrt(sum);
/*----------------------------------------------- put first column to Q */
for (i=0; i<nrow; i++) 
   P[i][start] /= sum;
/*----------------------------------------------- put first column in R */
R[start][start] = sum;
/*------------------------------------------ loop all remaining columns */
for (j=start+1; j<ncol; j++)
{
   for (i=0; i<j; i++)
   {
      sum = 0.0;
      for (k=0; k<nrow; k++)
         sum += P[k][j] * P[k][i];

      R[i][j] = sum;

      for (k=0; k<nrow; k++)
         P[k][j] -= P[k][i] * sum;
   }

   sum = 0.0;
   for (k=0; k<nrow; k++)
      sum += P[k][j] * P[k][j];

   if (FABS(sum) > EPS12)
   {
      sum     = sqrt(sum);
      R[j][j] = sum;
      for (k=0; k<nrow; k++)
         P[k][j] /= sum;
   }
   else
      R[j][j] = 0.0;
} /* end of for (j=0; j<ncol; j++) */
/*----------------------------------------------------------------------*/
exit:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_gramschmidt */







/*! @} (documentation module close)*/
