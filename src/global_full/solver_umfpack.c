#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
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
|  control solver Umpfack                                      mn 05/03 |
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

DOUBLE         t;
INT            status;

ARRAY          b_a;
double        *b;
ARRAY          x_a;
double        *x;
ARRAY          tmp_a;
double        *tmp;

/* some umfpack variables see umfpack manual */
/*double        *control = (double *) NULL;
double        *info    = (double *) NULL;*/
void          *symbolic, *numeric;
double        info [UMFPACK_INFO], control [UMFPACK_CONTROL];

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

   /* get the default control parameters */
   umfpack_di_defaults (control) ;

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
       t = umfpack_timer ( ) ;
       /* get the default control parameters */
       umfpack_di_defaults (control) ;
       control [UMFPACK_PIVOT_TOLERANCE] = 0.2;
       control [UMFPACK_IRSTEP]          =   5;

#ifdef DEBUG 
       control [UMFPACK_PRL] = 5 ;
#endif

       /* symbolic factorization */
       status = umfpack_di_symbolic (n, n, Ap, Ai, &symbolic, control, info);
#ifdef DEBUG 
       if (status < 0)
       {
           umfpack_di_report_info (control, info) ;
           umfpack_di_report_status (control, status) ;
           dserror ("umfpack_di_symbolic failed") ;
       }
#endif

       /* numeric factorization */
       status =  umfpack_di_numeric (Ap, Ai, Ax, symbolic, &numeric, control, info);
#ifdef DEBUG 
       if (status < 0)
       {
           umfpack_di_report_info (control, info) ;
           umfpack_di_report_status (control, status) ;
           dserror ("umfpack_di_numeric failed") ;
       }
#endif

       umfpack_di_free_symbolic (&symbolic);

       /* solve Ax=b */
       status = umfpack_di_solve (UMFPACK_A, Ap, Ai, Ax, x, b, numeric, control, info);
#ifdef DEBUG 
       if (status < 0)
       {
           umfpack_di_report_info (control, info) ;
           umfpack_di_report_status (control, status) ;
           dserror ("umfpack_di_solve failed") ;
       }
#endif

       umfpack_di_free_numeric (&numeric);
       t = umfpack_timer ( ) - t ;
       fprintf(allfiles.out_err,"umfpack solve complete.  Total time: %5.2f (seconds)\n", t);

/*#ifdef DEBUG 
   control [UMFPACK_PRL] = 5 ;
   umfpack_di_report_info (control, info);
#endif*/

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
