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
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_smoJacobi");
#endif
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
exit:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_smoJacobi */
/*!---------------------------------------------------------------------
\brief Jacobi smoother                                              

<pre>                                                        m.gee 11/02 

</pre>
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
ARRAY    sdense_a,rdense_a,srhs_a,rrhs_a,ipiv;
double **sdense,**rdense,*srhs,*rrhs;
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
/*---- allocate temporary dense matrix and right hand side and solution */
sdense = amdef("tmp",&sdense_a,numeq_total,numeq_total,"DA");
         amzero(&sdense_a);
rdense = amdef("tmp",&rdense_a,numeq_total,numeq_total,"DA");
srhs   = amdef("tmp",&srhs_a,numeq_total,1,"DV");
         amzero(&srhs_a);
rrhs   = amdef("tmp",&rrhs_a,numeq_total,1,"DV");

         amdef("tmp",&ipiv,numeq_total,1,"IV");
/*------------------------------------------------ fill sdense and srhs */
for (i=0; i<numeq; i++)
{
   actrow = update[i];
   colstart = ia[i];
   colend   = ia[i+1];
   for (j=colstart; j<colend; j++)
   {
      actcol = ja[j];
      sdense[actrow][actcol] = a[j];
   }
   srhs[actrow] = r[i];
}
/*------------------------------------------------------ allreduce them */
#ifdef PARALLEL 
MPI_Allreduce(sdense[0],rdense[0],sdense_a.fdim*sdense_a.sdim,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
MPI_Allreduce(srhs     ,rrhs     ,srhs_a.fdim                ,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
#endif
/*---------------------------------------------- now solve using lapack */
info = 1;
dgetrf(&numeq_total,
       &numeq_total,
       rdense[0],
       &numeq_total,
       ipiv.a.iv,
       &info);
if (info != 0) dserror("Lapack dgetrf returned info nonzero");
info     = 1;
trans[0] = 'N';
#ifndef AZTEC_PACKAGE
dgetrs(trans,
       &numeq_total,
       &ione,
       rdense[0],
       &numeq_total,
       ipiv.a.iv,
       rrhs,
       &numeq_total,
       &info);
if (info != 0) dserror("Lapack dgetrf returned info nonzero");
#else
   dserror("solver Lapack conflicts with compilation with -DAZTEC_PACKAGE");
#endif
/*------------------ fill the solution back to the distributed vector z */
for (i=0; i<numeq; i++)
{
   actrow = update[i];
   z[i]   = rrhs[actrow];
}
/*------------------------------------------------------------- tidy up */
amdel(&sdense_a);
amdel(&rdense_a);
amdel(&srhs_a);
amdel(&rrhs_a);
amdel(&ipiv);
/*----------------------------------------------------------------------*/
exit:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_lapacksolve */

/*! @} (documentation module close)*/
