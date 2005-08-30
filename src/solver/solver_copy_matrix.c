/*!----------------------------------------------------------------------
\file
\brief routines to copy matricex to another storage format

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../solver/solver.h"

/*!---------------------------------------------------------------------
\brief copy matrix csr-format to a other format

<pre>                                                        genk  11/02
                                                             basol 12/02
In this routine a matrix in CSR-format is copied to the actual matrix
format which is used by the solver

</pre>

\param  *amatrix      SPARSE_ARRAY   a  sparse matrix              (o)
\param  *amatrix_typ  SPARSE_TYP     sparse_typ of the matrx       (i)
\param  *amatrix_csr  DBCSR          a matrix in CSR-format        (i)
\param  *actintra     INTRA	     total number of equations     (i)

\return void
\warning this is all in progress
\warning only csr to msr available at the moment

------------------------------------------------------------------------*/
void solver_copy_csr(
                      SPARSE_ARRAY  *amatrix,
                      SPARSE_TYP    *amatrix_typ,
                      DBCSR         *amatrix_csr,
		      int            numeq_total
                     )
{

int     i,j,k;
int     iptr,iptr2;
int     icount;
int    *bindx;
int    *ia,*ja;
int    *iwk;
int     nnz;
double *a;
double *val;
double *wk;
ARRAY   wk_a,iwk_a;
ARRAY   dummy_a;
int *dummy;

#ifdef DEBUG
dstrc_enter("fluid_pmcpamatrix");
#endif

/*--------------------------------------------------- set some pointers */
a  = amatrix_csr->a.a.dv;
ia = amatrix_csr->ia.a.iv;
ja = amatrix_csr->ja.a.iv;
nnz = amatrix_csr->a.fdim;
dummy = amdef("dummy",&dummy_a,amatrix_csr->a.fdim,1,"IV");
amzero(&dummy_a);

#ifdef PARALLEL
switch (*amatrix_typ)
{
case msr:
dserror("copy of matrix failed: parallel version csrmsr not implemented yet!!!\n");
break;
default:
dserror("copy of matrix failed: sysarray_typ not implemented yet!!!\n");
}


#else
switch (*amatrix_typ)
{
case msr: /* copy A-matrix csr to msr format */
   /*------------------------------------------ allocate working arrays */
   wk  = amdef("wk" ,&wk_a ,numeq_total  ,1,"DV");
   iwk = amdef("iwk",&iwk_a,numeq_total+1,1,"IV");
   amzero(&wk_a);
   amzero(&iwk_a);
   /*------------------------------------------------ set some pointers */
   /*------redefine the size of the msr array------------------------*/
   amredef(&(amatrix->msr->val),nnz+1,1,"DV");
   amzero(&(amatrix->msr->val));
   amredef(&(amatrix->msr->bindx),nnz+1,1,"IV");
   amzero(&(amatrix->msr->bindx));
   amatrix->msr->nnz = nnz;
   /*-----------------------------------------------------------------*/
   val   = amatrix->msr->val.a.dv;
   bindx = amatrix->msr->bindx.a.iv;
   /*-------------------------------------------------- copy csr to msr */
   icount=0;
   /*- store away diagonal elements and count nonzero diagonal elements */
   for (i=0;i<numeq_total;i++)
   {
      wk[i]=ZERO;
      iwk[i+1] = ia[i+1]-ia[i];
      for (k=ia[i];k<ia[i+1];k++)
      {
         if (ja[k]==i)
	 {
	    wk[i]=a[k];
	    icount++;
	    iwk[i+1] = iwk[i+1]-1;
	 }
      }
   }
   /*--------------------------------------------- compute total length */
   iptr = numeq_total + ia[numeq_total] - icount + 1;
   iptr2 = numeq_total + ia[numeq_total] - icount;

   /*--------------------------------------------------- copy backwards */
   for (i=numeq_total-1;i>=0;i--)
   {
      for (k=ia[i+1]-1;k>=ia[i];k--)
      {
         j=ja[k];
	 if (j!=i)
	 {
	    iptr--;
	    val[iptr] = a[k];
	    bindx[iptr-1] = j;
         }
      }
   }
   for (i=0;i<iptr2;i++) dummy[i]=bindx[i];
   /*------------------------------- compute pointer values and copy wk */
   bindx[0] = numeq_total+1;
   for (i=0;i<numeq_total;i++)
   {
      val[i] = wk[i];
      bindx[i+1] = bindx[i]+iwk[i+1];
   }
   for (i=numeq_total;i<iptr2+1;i++)
   {
      bindx[i+1] = dummy[i];
   }
   /*-------------------------------------------- delete working arrays */
   amdel(&wk_a);
   amdel(&iwk_a);
break;
default:
dserror("copy of matrix failed: sysarray_typ not implemented yet!!!\n");
}
#endif

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of fluid_pmcpamatrix */
