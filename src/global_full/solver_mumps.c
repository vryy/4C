#include "../headers/standardtypes.h"
#include "../headers/solution.h"
/*----------------------------------------------------------------------*
 |  control solver lib Lapack                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void solver_mumps( 
                              struct _SOLVAR         *actsolv,
                              struct _INTRA          *actintra,
                              struct _RC_PTR         *rc_ptr,
                              struct _DIST_VECTOR    *sol,
                              struct _DIST_VECTOR    *rhs,
                              int                     option
                             )
{
#ifdef MUMPS_PACKAGE
int            i;
int            dof;
int            info;
int            imyrank;
int            inprocs;

int            job;
int            sym;
int            comm;
int            parproc;
int            n;
int            nz;
int            nz_loc;
int           *irn_loc;
int           *jcn_loc;
double        *A_loc;
int            icntl[20];

MUMPSVARS     *mumpsvars;

ARRAY          b_a;
double        *b;
ARRAY          tmp_a;
double        *tmp;

#ifdef DEBUG 
dstrc_enter("solver_mumps");
#endif
/*----------------------------------------------------------------------*/
imyrank      = actintra->intra_rank;
inprocs      = actintra->intra_nprocs;
mumpsvars    = actsolv->mumpsvars;
/*----------------------------------------------------------------------*/
switch(option)
{
/*----------------------------------------------------------------------*/
/*                                                           init phase */
/*----------------------------------------------------------------------*/
case 1:
   switch(actsolv->solvertyp)
   {
   case mumps_sym:    sym=2; break;
   case mumps_nonsym: sym=0; break;
   default: dserror("Unknown typ of solver"); break;
   }
   /*------------------------------------------- set some other values */
   job     = -1;
   parproc =  1;
#ifdef PARALLEL 
   comm    =  MPI_Comm_c2f(actintra->MPI_INTRA_COMM);
#endif
   n       =  rc_ptr->numeq_total;
   nz      =  rc_ptr->nnz_total;
   nz_loc  =  rc_ptr->nnz;
   irn_loc =  rc_ptr->irn_loc.a.iv;
   jcn_loc =  rc_ptr->jcn_loc.a.iv;
   A_loc   =  rc_ptr->A_loc.a.dv;
   mumps_interface(&job,
                   &parproc,
                   &comm,
                   &sym,
                   &n,
                   &nz,
                   &nz_loc,
                   irn_loc,
                   jcn_loc,
                   A_loc);
   /* set flag, that this matrix has been initialized and is ready for solve */   
   rc_ptr->is_init    = 1;
   rc_ptr->ncall      = 0;
   rc_ptr->is_factored = 0;
break;
/*----------------------------------------------------------------------*/
/*                                                    end of init phase */
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*                                                    calculation phase */
/*----------------------------------------------------------------------*/
case 0:
   /*--------------------------------- allocate rhs and solution vector */
   b   = amdef("b",&b_a,rc_ptr->numeq_total,1,"DV");
         amzero(&b_a);
#ifdef PARALLEL 
   tmp = amdef("tmp",&tmp_a,rc_ptr->numeq_total,1,"DV");
         amzero(&tmp_a);
   /*--------------------------------------------------------- fill rhs */
   for (i=0; i<rhs->numeq; i++)
   {
      dof      = rc_ptr->update.a.iv[i];
      tmp[dof] = rhs->vec.a.dv[i]; 
   }
   /*--------------------------------------------- allreduce rhs vector */
   MPI_Allreduce(tmp,b,rc_ptr->numeq_total,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
   amdel(&tmp_a);
#else
   /*--------------------------------------------------------- fill rhs */
   for (i=0; i<rhs->numeq; i++)
   {
      dof      = rc_ptr->update.a.iv[i];
      b[dof]   = rhs->vec.a.dv[i]; 
   }
#endif
   switch(actsolv->solvertyp)
   {
   case lapack_sym:
   break;
   case lapack_nonsym:/*-------------- nonsymmetric lapack dense solver */
   break;
   default:
      dserror("Unknown typ of solver");
   break;
   }
   /*---------- set flag, that this matrix is factored and count solves */
   rc_ptr->ncall++;
   rc_ptr->is_factored=1;
   /*------------------------------------------------------- get result */
   for (i=0; i<sol->numeq; i++)
   {
      dof              = rc_ptr->update.a.iv[i];
      sol->vec.a.dv[i] = b[dof];
   }
break;
/*----------------------------------------------------------------------*/
/*                                             end of calculation phase */
/*----------------------------------------------------------------------*/
default:
   dserror("Unknown option for solver call to Mumps");
break;   
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
#endif
return;
} /* end of solver_mumps */




