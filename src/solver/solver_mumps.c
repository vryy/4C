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

#ifdef MUMPS_PACKAGE


#include "../headers/standardtypes.h"
#include "../solver/solver.h"

/*----------------------------------------------------------------------*
 |  control solver lib MUMPS                             m.gee 11/01    |
 *----------------------------------------------------------------------*/
void solver_mumps(struct _SOLVAR         *actsolv,
                  struct _INTRA          *actintra,
                  struct _RC_PTR         *rc_ptr,
                  struct _DIST_VECTOR    *sol,
                  struct _DIST_VECTOR    *rhs,
                  INT                     option)
{
INT            i;
INT            dof;
INT            info;
INT            imyrank;
INT            inprocs;

INT            job;
INT            sym;
INT            parproc;
INT            n;
INT            nz;
INT            nz_loc;
INT           *irn_loc;
INT           *jcn_loc;
INT           *irn;
INT           *jcn;
DOUBLE        *A_loc;
INT           *icntl;

MUMPSVARS     *mumpsvars;

ARRAY          b_a;
DOUBLE        *b;
ARRAY          tmp_a;
DOUBLE        *tmp;

#ifdef DEBUG 
dstrc_enter("solver_mumps");
#endif
/*----------------------------------------------------------------------*/
imyrank      = actintra->intra_rank;
inprocs      = actintra->intra_nprocs;
mumpsvars    = actsolv->mumpsvars;
switch(actsolv->solvertyp)
{
case mumps_sym:    sym=2; break;
case mumps_nonsym: sym=0; break;
default: dserror("Unknown typ of solver"); break;
}
/*----------------------------------------------------------------------*/
switch(option)
{
/*----------------------------------------------------------------------*/
/*                                                           init phase */
/*----------------------------------------------------------------------*/
case 1:
   /*--------------------------------------- This will only do on HPUX! */
   /*---------------------------------------- This will only do on SUN! */
#ifdef PARALLEL 
   rc_ptr->comm  =  MPI_Comm_c2f(actintra->MPI_INTRA_COMM);
#endif
   /*------------------- copy the pointer vectors to fortran numbering */
   am_alloc_copy(&(rc_ptr->irn_loc),&(rc_ptr->irn_locf));
   am_alloc_copy(&(rc_ptr->jcn_loc),&(rc_ptr->jcn_locf));
   for (i=0; i<rc_ptr->irn_locf.fdim; i++)
   {
      rc_ptr->irn_locf.a.iv[i]++;
      rc_ptr->jcn_locf.a.iv[i]++;
   }
   if (imyrank==0)
   {
      for (i=0; i<rc_ptr->irn_glob.fdim; i++)
      {
         rc_ptr->irn_glob.a.iv[i]++;
         rc_ptr->jcn_glob.a.iv[i]++;
      }
   }
   /*------------------------------------------ call the solver to init */
   /*--------------------------------- input is dist. assembled matrix */
   /*
   the matrix is rowsum distributed in a_loc, irn_loc and jcn_loc. 
   for the mumps analysis phase, the complete matrix sparsity pattern is given 
   on imyrank==0 and is then deallocated
   */
   job       = -1;                       /* analysis phase */
   parproc   =  1;                       /* imyrank=0 takes part in factorization */
   icntl     =  rc_ptr->icntl;           /* icntl[0..19] are MUMPS options */
   icntl[17] =  3;                       /* see MUMPS manual */
   icntl[2]  =  0;                       /* no output to stdout by solver */
   n         =  rc_ptr->numeq_total;
   nz        =  rc_ptr->nnz_total;
   nz_loc    =  rc_ptr->nnz;
   irn_loc   =  rc_ptr->irn_locf.a.iv;
   jcn_loc   =  rc_ptr->jcn_locf.a.iv;
   A_loc     =  rc_ptr->A_loc.a.dv;
   if (imyrank==0)
   {
   irn       =  rc_ptr->irn_glob.a.iv;
   jcn       =  rc_ptr->jcn_glob.a.iv;
   }
   else
   {
   irn       =  NULL;
   jcn       =  NULL;
   }
   mumps_interface(&job,
                   &parproc,
                   &(rc_ptr->comm),
                   &sym,
                   icntl,
                   &n,
                   &nz,
                   &nz_loc,
                   irn_loc,
                   jcn_loc,
                   irn,
                   jcn,
                   A_loc,
                   NULL);
   /* the global sparsity pattern is no longer needed after analysis phase */
   if (imyrank==0)
   {
   amdel(&(rc_ptr->irn_glob));
   amdel(&(rc_ptr->jcn_glob));
   }
   /* set flag, that this matrix has been initialized and is ready for solve */   
   rc_ptr->is_init     = 1;
   rc_ptr->ncall       = 0;
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
   /*------------------------------------------- set some other values */
   if (rc_ptr->is_factored) 
   {
      job = 3;
   }
   else
   {
      job = 6;
      amcopy(&(rc_ptr->irn_loc),&(rc_ptr->irn_locf));
      amcopy(&(rc_ptr->jcn_loc),&(rc_ptr->jcn_locf));
      for (i=0; i<rc_ptr->irn_locf.fdim; i++)
      {
         rc_ptr->irn_locf.a.iv[i]++;
         rc_ptr->jcn_locf.a.iv[i]++;
      }
   }
   parproc   =  1;
   icntl     =  rc_ptr->icntl;
   icntl[17] =  2;
   icntl[2]  =  0;
   n         =  rc_ptr->numeq_total;
   nz        =  rc_ptr->nnz_total;
   nz_loc    =  rc_ptr->nnz;
   irn_loc   =  rc_ptr->irn_locf.a.iv;
   jcn_loc   =  rc_ptr->jcn_locf.a.iv;
   A_loc     =  rc_ptr->A_loc.a.dv;
   mumps_interface(&job,
                   &parproc,
                   &(rc_ptr->comm),
                   &sym,
                   icntl,
                   &n,
                   &nz,
                   &nz_loc,
                   irn_loc,
                   jcn_loc,
                   NULL,
                   NULL,
                   A_loc,
                   b);
   /*---------- set flag, that this matrix is factored and count solves */
   rc_ptr->ncall++;
   rc_ptr->is_factored=1;
   /*------------------------------------------------- broadcast result */
#ifdef PARALLEL 
   MPI_Bcast(b,n,MPI_DOUBLE,0,actintra->MPI_INTRA_COMM);
#endif
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
return;
} /* end of solver_mumps */



#endif /* ifdef MUMPS_PACKAGE */

