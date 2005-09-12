/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

*----------------------------------------------------------------------*/

#ifdef UMFPACK


#include "../headers/standardtypes.h"
#include "../solver/solver.h"
/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;


/*!----------------------------------------------------------------------
\brief controls the solver umfpack

<pre>                                                              mn 05/03
This functions controls the solution with the solver umfpack.
</pre>
\param *actsolv    SOLVAR       (i)   the active solver
\param *actintra   INTRA        (i)   the intra communicator
\param *ccf        CCF          (i)   the system matrix in ccf format
\param *sol        DIST_VECTOR  (o)   the distributed solution vector
\param *rhs        DIST_VECTOR  (i)   the distributed rhs vector
\param options     INT          (i)   flag for init or not

\warning There is nothing special to this routine
\return void
\sa

*----------------------------------------------------------------------*/

void solver_umfpack(struct _SOLVAR         *actsolv,
                    struct _INTRA          *actintra,
                    struct _CCF            *ccf,
                    struct _DIST_VECTOR    *sol,
                    struct _DIST_VECTOR    *rhs,
                    INT                     option)
{
INT            i;
INT            dof;
INT            imyrank;
INT            inprocs;

INT            n;
INT            nz;
INT            nz_loc;
INT           *Ap;
INT           *Ai;
DOUBLE        *Ax;

DOUBLE         t;
INT            status;

ARRAY          b_a;
DOUBLE        *b;
ARRAY          x_a;
DOUBLE        *x;

#ifdef PARALLEL
ARRAY          tmp_a;
DOUBLE        *tmp;
#endif

/* some umfpack variables see umfpack manual */
/*DOUBLE        *control = (DOUBLE *) NULL;
DOUBLE        *info    = (DOUBLE *) NULL;*/
static void          *symbolic, *numeric;
DOUBLE        info [UMFPACK_INFO], control [UMFPACK_CONTROL];

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
   ccf->reuse       = 0;

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
       control [UMFPACK_PIVOT_TOLERANCE] = 1.0;
       control [UMFPACK_IRSTEP]          =   5;

#ifdef DEBUG
       control [UMFPACK_PRL] = 5 ;
#endif
       if (ccf->reuse==0)/* last LU-decomposition is not to be used */
       {
          if (ccf->is_factored==1)/* free pointer *numeric of last decomposition*/
          {
             umfpack_di_free_numeric (&numeric);
          }
          /* symbolic -> (with respect to the mask of matrix A) factorization */

#if defined(LINUX_MUENCH) || defined(HPUX_MUENCH) || defined(WIN_MUENCH) || defined(TX7) || defined(SX8)
       status = umfpack_di_symbolic (n, n, Ap, Ai, Ax, &symbolic, control, info);
#else
       status = umfpack_di_symbolic (n, n, Ap, Ai, &symbolic, control, info);
#endif
#ifdef DEBUG 
          if (status < 0)
          {
             umfpack_di_report_info (control, info) ;
             umfpack_di_report_status (control, status) ;
             dserror("umfpack_di_symbolic failed") ;
          }
#endif
          /* numeric-> (with respect to the values of matrix A) factorization */
          status =  umfpack_di_numeric (Ap, Ai, Ax, symbolic, &numeric, control, info);
#ifdef DEBUG 
          if (status < 0)
          {
             umfpack_di_report_info (control, info) ;
             umfpack_di_report_status (control, status) ;
             dserror("umfpack_di_numeric failed") ;
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
             dserror("umfpack_di_solve failed") ;
          }
#endif
          t = umfpack_timer ( ) - t ;
          fprintf(allfiles.out_err,"umfpack solve complete.  Total time: %5.2f (seconds)\n", t);
          ccf->is_factored = 1;
       }/* if (ccf->reuse==0)*/

       else if (ccf->reuse==1)/* use last LU-decomposition (if only RHS has changed) */
       {
          /* solve Ax=b */
          status = umfpack_di_solve (UMFPACK_A, Ap, Ai, Ax, x, b, numeric, control, info);
#ifdef DEBUG 
          if (status < 0)
          {
             umfpack_di_report_info (control, info) ;
             umfpack_di_report_status (control, status) ;
             dserror("umfpack_di_solve failed") ;
          }
#endif

          t = umfpack_timer ( ) - t ;
          fprintf(allfiles.out_err,"umfpack solve complete.  Total time: %5.2f (seconds)\n", t);
       }/* end else if (ccf->reuse==1)*/

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
   ccf->is_factored = 1;
   /*ccf->reuse       = 1;*/

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
return;
} /* end of solver_umfpack */


#endif /* ifdef UMFPACK */

