#ifdef MLPCG

/*!---------------------------------------------------------------------
\file
\brief contains the multilevel preconditioner for shells

<pre>
Maintainer: Michael Gee
            gee@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/gee/
            0711 - 685-6572
</pre>

---------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
INT cmp_int(const void *a, const void *b );
DOUBLE cmp_double(const void *a, const void *b );
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
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;



/*!---------------------------------------------------------------------
\brief create the tentative prolongator from actlev+1 to actlev

<pre>                                                        m.gee 9/02
create the tentative prolongator from actlev+1 to actlev
from the aggregation done before , this routine only from 0 to 1
</pre>
\param actlev      MLLEVEL*     (i/o) the active level in the ml-precond.
\param actintra   INTRA*       (i)   the intra-communicator of this field
\return void

------------------------------------------------------------------------*/
void mlpcg_precond_P_fish(MLLEVEL  *actlev, INTRA *actintra)
{
INT          i,j,k,n,counter=0;
INT          myrank,nproc;
DBCSR       *actstiff;
DBCSR       *P;
AGG         *actagg;
INT          nrow,ncol;
DOUBLE       aggblock[1000][500];
INT          rindex[1000],cindex[500];
INT          firstdof=0;
INT          sendbuff[MAXPROC],recvbuff[MAXPROC];
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("mlpcg_precond_P_fish");
#endif
/*----------------------------------------------------------------------*/
myrank = actintra->intra_rank;
nproc  = actintra->intra_nprocs;
/*----------------------------------------------------------------------*/
actstiff    = actlev->csr;
/*----------------------------------------------------------------------*/
/*
loop the aggregates and create the tentative prolongator
*/
for (i=0; i<actlev->nagg; i++)
{
   actagg = &(actlev->agg[i]);
   /* create the tentative prolongator diagonal block of this aggregate */
   mlpcg_precond_oneP_fish(actagg,aggblock,rindex,cindex,&nrow,&ncol,actstiff,actintra);
   if (!(actagg->tentP))
   {
      actagg->tentP = (ARRAY*)CCACALLOC(1,sizeof(ARRAY));
      amdef("tentP",actagg->tentP,nrow,ncol,"DA");
   }
   else if (actagg->tentP->fdim != nrow)
       dserror("Size mismatch in aggregate");
   else if (actagg->tentP_nrow != nrow)
       dserror("Size mismatch in aggregate");
   if (!(actagg->tentP_rindex))
      actagg->tentP_rindex = (INT*)CCAMALLOC(nrow*sizeof(INT));
   /*------------------------------- put rindex to actagg->tentP_rindex */
   if (mlprecond.ncall==0)
      for (j=0; j<nrow; j++) actagg->tentP_rindex[j] = rindex[j];
   /*--------------------------- put the aggblock to actagg->tentP.a.da */
   for (j=0; j<ncol; j++)
   for (k=0; k<nrow; k++)
      actagg->tentP->a.da[k][j] = aggblock[k][j];
   /*------------------ set number of rows in this piece of Prolongator */
   actagg->tentP_nrow = nrow;
   /*------------------------ store the R part to use in the next level */
} /* end of for (i=0; i<actlev->nagg; i++) */
/*------------------------------- set the dof numbers to the aggregates */
if (mlprecond.ncall==0)
{
   counter=0;
   for (i=0; i<actlev->nagg; i++) counter += actlev->agg[i].numdf;
   for (n=0; n<nproc; n++)
   {
      if (n==myrank) sendbuff[n] = counter;
      else           sendbuff[n] = 0;
   }
#ifdef PARALLEL
   MPI_Allreduce(sendbuff,recvbuff,nproc,MPI_INT,MPI_SUM,actintra->MPI_INTRA_COMM);
#endif
   for (n=0; n<myrank; n++) firstdof += recvbuff[n];
   for (i=0; i<actlev->nagg; i++)
   {
      actlev->agg[i].dof = (INT*)CCAMALLOC(actlev->agg[i].numdf * sizeof(INT));
      for (j=0; j<actlev->agg[i].numdf; j++) actlev->agg[i].dof[j] = firstdof++;
   }
/*-------------------------------------------------- now open the matrix */
   if (actlev->P == NULL)
      actlev->P = (DBCSR*)CCACALLOC(1,sizeof(DBCSR));
   P = actlev->P;
   P->firstcoupledof=actstiff->firstcoupledof;
   am_alloc_copy(&(actstiff->blocks),&(P->blocks));
   /* make a guess of the size of the prolongator */
   counter = (INT)(counter*(actstiff->numeq)*4.0);
   /* open the prolongator */
   mlpcg_csr_open(P,actstiff->update.a.iv[0],actstiff->update.a.iv[actstiff->update.fdim-1],
                  actstiff->numeq_total,counter, actintra);
}/* end of if (mlprecond.ncall==0) */
else
   mlpcg_csr_zero(P,actintra);
/*--------------- loop all aggregates again and fill the DBCSR matrix P */
for (i=0; i<actlev->nagg; i++)
{
   for (j=0; j<actlev->agg[i].numdf; j++) /* column loop */
   for (k=0; k<actlev->agg[i].tentP_nrow; k++) /* row loop */
   aggblock[k][j] = actlev->agg[i].tentP->a.da[k][j];
   mlpcg_csr_setblock(P,
                      aggblock,
                      actlev->agg[i].tentP_rindex,
                      actlev->agg[i].dof,
                      actlev->agg[i].tentP_nrow,
                      actlev->agg[i].numdf,
                      actintra);
}
/*------------------------------- make the smoothing of the prolongator */
/*
if (mlprecond.omega>0.0)
   mlpcg_smoothP(P,aggblock,rindex,cindex,&nrow,&ncol,actlev->csr,actlev->agg,actlev->nagg,actintra);
*/
/*------------------------ tent. prolongator is ready, close the matrix */
mlpcg_csr_close(P);
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_P_fish */





/*!---------------------------------------------------------------------
\brief create the tentative prolongator for one aggregate

<pre>                                                        m.gee 12/02

</pre>
\param actagg         AGG*    (i/o) the active aggregate
\param aggblock       DOUBLE[1000][500] (o) the aggregate's block prolongator
\param rindex         INT[1000]         (o) global indizes of aggblock
\param cindex         INT[500]          (o) global indizes of aggblock
\param nrow           INT*              (o) dimension of rindex
\param ncol           INT*              (o) dimension of cindex
\param actstiff       DBCSR*            (i) fine grid stiffness matrix
\param actintra       INTRA*            (i) the intra-communicator of this field
\return void

------------------------------------------------------------------------*/
void mlpcg_precond_oneP_fish(AGG     *actagg,
                             DOUBLE   aggblock[][500],
                             INT      rindex[],
                             INT      cindex[],
                             INT     *nrow,
                             INT     *ncol,
                             DBCSR   *actstiff,
                             INTRA   *actintra)
{
INT           i,j,counter;
DOUBLE      **A;
ARRAY         A_a;
DOUBLE      **Z;
ARRAY         Z_a;

INT           itype=1;
char          jobz[1];
char          range[1];
char          uplo[1];
DOUBLE        vl,vu;
INT           il,iu;
DOUBLE        abstol = EPS14;
DOUBLE        W[500];
INT           lwork = 10000;
DOUBLE        work[10000];
INT           iwork[10000];
INT           ifail[500];
INT           info;
jobz[0]  = 'V';
range[0] = 'I'; /* A gives all eigenvalues and vectors */
uplo[0]  = 'U';
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("mlpcg_precond_oneP_fish");
#endif
/*------------------- get size and global dofs of the prolongator block */
*nrow=0;
for (i=0; i<actagg->nblock; i++)
   (*nrow) += actagg->block[i][0];
if (*nrow >= 1000)
dserror("Local variable aggblock[1000][500] too small");
/*------------------------------------------ get the global row indizes */
counter=0;
for (i=0; i<actagg->nblock; i++)
   for (j=0; j<actagg->block[i][0]; j++)
   {
      rindex[counter] = actagg->block[i][j+1];
      counter++;
   }
dsassert(counter==*nrow,"Number of dofs in prolongator wrong");
qsort((INT*)rindex,*nrow,sizeof(INT),cmp_int);
/*------------------ extract the block, which belongs to this aggregate */
A = amdef("A",&A_a,*nrow,*nrow,"DA");
Z = amdef("Z",&Z_a,*nrow,*nrow,"DA");
amzero(&A_a);
mlpcg_csr_extractsubblock_dense(actstiff,A,rindex,*nrow,actintra);
/*---------------------------------------------------- call eigensolver */
vl = -10.0;
vu = mlprecond.gamma;
il = 1;
iu = IMIN(mlprecond.numdf,*nrow);
mydsyevx(jobz,range,uplo,nrow,A[0],nrow,&vl,&vu,&il,&iu,
         &abstol,ncol,W,Z[0],nrow,work,&lwork,iwork,ifail,&info);
if (info != 0) dserror("Eigenanalysis for prolongator went wrong");
j=0;
for (i=0; i<*ncol; i++) j += ifail[i];
if (j != 0) dserror("Eigenanalysis for prolongator went wrong");

if (*ncol > mlprecond.numdf) *ncol = mlprecond.numdf;

for (i=0; i<*ncol; i++)
{
   if (W[i] > mlprecond.gamma)
   {
      printf("Cutoff of eigenvalue %f appeared\n",W[i]);
      break;
   }
}
*ncol = i;

if (actagg->numdf != 0 && actagg->numdf != *ncol) dserror("Number of dofs changed");

actagg->numdf = *ncol;

if (actagg->numdf==0) dserror("Aggregate with no dofs detected");

for (j=0; j<*ncol; j++)
for (i=0; i<*nrow; i++)
   aggblock[i][j] = Z[j][i]; /* this is the transformation which happens between c and fortran */
/*----------------------------------------------------------------------*/
amdel(&A_a);
amdel(&Z_a);
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_oneP_fish */





/*! @} (documentation module close)*/
#endif
