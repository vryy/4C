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
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
/*----------------------------------------------------------------------*
 |  control solver lib Lapack                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void solver_lapack_dense(
                              struct _SOLVAR         *actsolv,
                              struct _INTRA          *actintra,
                              struct _DENSE          *dense,
                              struct _DIST_VECTOR    *sol,
                              struct _DIST_VECTOR    *rhs,
                              INT                     option
                             )
{
INT            i;
INT            dof;
INT            info;
INT            imyrank;
INT            inprocs;
INT            ione=1;
INT            nb=2;

char           trans[1];
char           uplo[1];

LAPACKVARS    *lapackvar;

ARRAY          b_a;
DOUBLE        *b;
#ifdef PARALLEL
ARRAY          tmp_a;
DOUBLE        *tmp;
#endif

#ifdef DEBUG
dstrc_enter("solver_lapack_dense");
#endif
/*----------------------------------------------------------------------*/
imyrank      = actintra->intra_rank;
inprocs      = actintra->intra_nprocs;
lapackvar    = actsolv->lapackvars;
/*----------------------------------------------------------------------*/
switch(option)
{
/*----------------------------------------------------------------------*/
/*                                                           init phase */
/*----------------------------------------------------------------------*/
case 1:
   switch(actsolv->solvertyp)
   {
   case lapack_sym:
      amdef("ipiv",&(dense->ipiv),dense->numeq_total,1,"IV");
      amzero(&(dense->ipiv));
      dense->lwork=nb*(dense->numeq_total);
      amdef("work",&(dense->work),dense->lwork,1,"DV");
      amzero(&(dense->work));
   break;
   case lapack_nonsym:
      amdef("ipiv",&(dense->ipiv),dense->numeq_total,1,"IV");
      amzero(&(dense->ipiv));
   break;
   default:
      dserror("Unknown typ of solver");
   break;
   }
   /* set flag, that this matrix has been initialized and is ready for solve */
   dense->is_init    = 1;
   dense->ncall      = 0;
   dense->is_factored = 0;
break;
/*----------------------------------------------------------------------*/
/*                                                    end of init phase */
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*                                                    calculation phase */
/*----------------------------------------------------------------------*/
case 0:
   /*--------------------------------- allocate rhs and solution vector */
   b   = amdef("b",&b_a,dense->numeq_total,1,"DV");
         amzero(&b_a);
#ifdef PARALLEL
   tmp = amdef("tmp",&tmp_a,dense->numeq_total,1,"DV");
         amzero(&tmp_a);
   /*--------------------------------------------------------- fill rhs */
   for (i=0; i<rhs->numeq; i++)
   {
      dof      = dense->update.a.iv[i];
      tmp[dof] = rhs->vec.a.dv[i];
   }
   /*--------------------------------------------- allreduce rhs vector */
   MPI_Allreduce(tmp,b,dense->numeq_total,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
   amdel(&tmp_a);
#else
   /*--------------------------------------------------------- fill rhs */
   for (i=0; i<rhs->numeq; i++)
   {
      dof      = dense->update.a.iv[i];
      b[dof]   = rhs->vec.a.dv[i];
   }
   /*--------------------------------------------- allreduce rhs vector */
#endif
   switch(actsolv->solvertyp)
   {
   case lapack_sym:
      info=1;
      uplo[0]='L';
      if (dense->is_factored==0)
      {
         dsytrf(
                  uplo,
                  &(dense->numeq_total),
                  dense->A.a.da[0],
                  &(dense->numeq_total),
                  dense->ipiv.a.iv,
                  dense->work.a.dv,
                  &(dense->lwork),
                  &info
               );
         if (info!=0) dserror("Lapack factorization failed");
      }
      info=1;
      dsytrs(
               uplo,
               &(dense->numeq_total),
               &ione,
               dense->A.a.da[0],
               &(dense->numeq_total),
               dense->ipiv.a.iv,
               b,
               &(dense->numeq_total),
               &info
            );
      if (info!=0) dserror("Lapack solve failed");
   break;
   case lapack_nonsym:/*-------------- nonsymmetric lapack dense solver */
      if (dense->is_factored==0)
      {
         info=1;
         dgetrf(
                 &(dense->numeq_total),
                 &(dense->numeq_total),
                 dense->A.a.da[0],
                 &(dense->numeq_total),
                 dense->ipiv.a.iv,
                 &info
               );
         if (info!=0) dserror("Lapack factorization failed");
      }
      info=1;
      trans[0]='N';
#ifndef AZTEC_PACKAGE
      dgetrs(
              trans,
              &(dense->numeq_total),
              &ione,
              dense->A.a.da[0],
              &(dense->numeq_total),
              dense->ipiv.a.iv,
              b,
              &(dense->numeq_total),
              &info
            );
#else
   dserror("solver Lapack conflicts with compilation with -DAZTEC_PACKAGE");
#endif
      if (info!=0) dserror("Lapack solve failed");
   break;
   default:
      dserror("Unknown typ of solver");
   break;
   }
   /*---------- set flag, that this matrix is factored and count solves */
   dense->ncall++;
   dense->is_factored=1;
   /*------------------------------------------------------- get result */
   for (i=0; i<sol->numeq; i++)
   {
      dof      = dense->update.a.iv[i];
      sol->vec.a.dv[i] = b[dof];
   }
break;
/*----------------------------------------------------------------------*/
/*                                             end of calculation phase */
/*----------------------------------------------------------------------*/
default:
   dserror("Unknown option for solver call to Lapack");
break;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of solver_lapack_dense */




