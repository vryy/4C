#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
/*----------------------------------------------------------------------*
|  control solver Umpfack                            s.offermanns 05/02 |
*----------------------------------------------------------------------*/
void solver_umfpack(struct _SOLVAR         *actsolv,
                    struct _INTRA          *actintra,
                    struct _CCF            *ccf,
                    struct _DIST_VECTOR    *sol,
                    struct _DIST_VECTOR    *rhs,
                    int                     option)
{
#ifdef UMFPACK 
int            i;
int            dof;
int            imyrank;
int            inprocs;

int            job;
int            sym;
int            n;
int            nz;
int            nz_loc;
int           *Ap;
int           *Ai;
double        *Ax;
int           *icntl;

ARRAY          b_a;
double        *b;
ARRAY          x_a;
double        *x;
ARRAY          tmp_a;
double        *tmp;

/* some umfpack variables see umfpack manual */
double        *control = (double *) NULL;
double        *info    = (double *) NULL;
void          *symbolic, *numeric;


#ifdef DEBUG 
dstrc_enter("solver_umfpack");
#endif
/*----------------------------------------------------------------------*/
imyrank      = actintra->intra_rank;
inprocs      = actintra->intra_nprocs;

switch(option)
{
/*----------------------------------------------------------------------*/
/*                                                           init phase */
/*----------------------------------------------------------------------*/
case 1:

   /* set flag, that this matrix has been initialized and is ready for solve */   
   ccf->is_init     = 1;
   ccf->ncall       = 0;
   ccf->is_factored = 0;

   break;
      
case 0:
   /*--------------------------------- allocate rhs and solution vector */
   b   = amdef("b",&b_a,ccf->numeq_total,1,"DV");
         amzero(&b_a);
   x   = amdef("x",&x_a,ccf->numeq_total,1,"DV");
         amzero(&x_a);
#ifdef PARALLEL
   tmp = amdef("tmp",&tmp_a,ccf->numeq_total,1,"DV");
         amzero(&tmp_a);
   /*--------------------------------------------------------- fill rhs */
   for (i=0; i<rhs->numeq; i++)
   {
      dof      = ccf->update.a.iv[i];
      tmp[dof] = rhs->vec.a.dv[i];
   }
   /*--------------------------------------------- allreduce rhs vector */
   MPI_Allreduce(tmp,b,ccf->numeq_total,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
/*   amdel(&tmp_a);*/
#else
   /*--------------------------------------------------------- fill rhs */
   for (i=0; i<rhs->numeq; i++)
   {
      dof      = ccf->update.a.iv[i];
      b[dof]   = rhs->vec.a.dv[i];
   }
#endif

   n         =  ccf->numeq_total;
   nz        =  ccf->nnz_total;
   nz_loc    =  ccf->nnz;
   Ai        =  ccf->Ai.a.iv;
   Ap        =  ccf->Ap.a.iv;
   Ax        =  ccf->Ax.a.dv;
 
   /*------------------------------ solution is only done on imyrank==0 */   
   if (imyrank==0)
   {
      (void) umfpack_di_symbolic (n, n, Ap, Ai, &symbolic, control, info);
      (void) umfpack_di_numeric (Ap, Ai, Ax, symbolic, &numeric, control, info);
      umfpack_di_free_symbolic (&symbolic);
      (void) umfpack_di_solve (UMFPACK_A, Ap, Ai, Ax, x, b, numeric, control, info);
      umfpack_di_free_numeric (&numeric);
   }/* end of (imyrank==0) */
   /*------------------------------------------------- broadcast result */
#ifdef PARALLEL
   /*MPI_Bcast(b,n,MPI_DOUBLE,0,actintra->MPI_INTRA_COMM);*/
   MPI_Bcast(x,n,MPI_DOUBLE,0,actintra->MPI_INTRA_COMM);
#endif
   /*------------------------------------------------------- get result */
   for (i=0; i<sol->numeq; i++)
   {
      dof              = ccf->update.a.iv[i];
      sol->vec.a.dv[i] = x[dof];
   }
   ccf->ncall++;
   ccf->is_factored=1;

   amdel(&x_a);
   amdel(&b_a);
#ifdef PARALLEL 
   amdel(&tmp_a);
#endif
break;
/*----------------------------------------------------------------------*/
/*                                             end of calculation phase */
/*----------------------------------------------------------------------*/
default:
   dserror("Unknown option for solver call to umfpack");
break;
}
  
#ifdef DEBUG 
dstrc_exit();
#endif
#endif /* end of ifdef UMFPACK */
return;
} /* end of solver_umfpack */
