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
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
/* 
the 64bit linker expects an underslash appended to calls for
fortran90 subroutines:
*/
/*----------------------------------------------------------------------*
 |  control solver colsol                                m.gee 01/02    |
 *----------------------------------------------------------------------*/
void solver_colsol(struct _SOLVAR         *actsolv,
                   struct _INTRA          *actintra,
                   struct _SKYMATRIX      *sky,
                   struct _DIST_VECTOR    *sol,
                   struct _DIST_VECTOR    *rhs,
                   INT                     option)
{
INT            i;
INT            dof;
INT            imyrank;
INT            inprocs;

ARRAY          b_a;
DOUBLE        *b;
ARRAY          tmp_a;
DOUBLE        *tmp;

INT            ione=1;
INT            izero=0;
INT            info=6;
INT            kkk;
INT            isc;
INT            nsch;

#ifdef DEBUG 
dstrc_enter("solver_colsol");
#endif
/*----------------------------------------------------------------------*/
imyrank      = actintra->intra_rank;
inprocs      = actintra->intra_nprocs;
/*----------------------------------------------------------------------*/
switch(option)
{
/*----------------------------------------------------------------------*/
/*                                                           init phase */
/*----------------------------------------------------------------------*/
case 1:
   /*---------------------------------- copy maxa to fortran convention */
   am_alloc_copy(&(sky->maxa),&(sky->maxaf));
   for (i=0; i<sky->maxaf.fdim; i++) sky->maxaf.a.iv[i]++;
   /* set flag, that this matrix has been initialized and is ready for solve */   
   sky->is_init     = 1;
   sky->ncall       = 0;
   sky->is_factored = 0;
break;
/*----------------------------------------------------------------------*/
/*                                                    end of init phase */
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*                                                    calculation phase */
/*----------------------------------------------------------------------*/
case 0:
   /*--------------------------------- allocate rhs and solution vector */
   b   = amdef("b",&b_a,sky->numeq_total,1,"DV");
         amzero(&b_a);
#ifdef PARALLEL 
   tmp = amdef("tmp",&tmp_a,sky->numeq_total,1,"DV");
         amzero(&tmp_a);
   /*--------------------------------------------------------- fill rhs */
   for (i=0; i<rhs->numeq; i++)
   {
      dof      = sky->update.a.iv[i];
      tmp[dof] = rhs->vec.a.dv[i]; 
   }
   /*--------------------------------------------- allreduce rhs vector */
   MPI_Allreduce(tmp,b,sky->numeq_total,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
/*   amdel(&tmp_a);*/
#else
   /*--------------------------------------------------------- fill rhs */
   for (i=0; i<rhs->numeq; i++)
   {
      dof      = sky->update.a.iv[i];
      b[dof]   = rhs->vec.a.dv[i]; 
   }
   /*------------------------------------------------------------------ */
#endif
   /*------------------------------ solution is only done on imyrank==0 */   
   if (imyrank==0)
   {
      if (sky->is_factored==0) kkk=3;
      else                     kkk=2;
      colsol(
             sky->A.a.dv,
             b,
             sky->maxaf.a.iv,
             &(sky->numeq_total),
             &(sky->numeq_total),
             &ione,
             &(sky->A.fdim),
             &(sky->maxaf.fdim),
             &ione, 
             &ione,
             &kkk,
             &(sky->det),
             &isc,
             &nsch,
             &izero,
             &info
            );
   }/* end of (imyrank==0) */
   /*------------------------------------------------ distribute result */
#ifdef PARALLEL 
      MPI_Bcast(b,sky->numeq_total,MPI_DOUBLE,0,actintra->MPI_INTRA_COMM);
#endif
   /*------------------------------------------------------- get result */
   for (i=0; i<sol->numeq; i++)
   {
      dof      = sky->update.a.iv[i];
      sol->vec.a.dv[i] = b[dof];
   }
   /*---------- set flag, that this matrix is factored and count solves */
   sky->ncall++;
   sky->is_factored=1;
   /*------------------------------------- deallocate temporary vectors */
   amdel(&b_a);
#ifdef PARALLEL 
   amdel(&tmp_a);
#endif
break;
/*----------------------------------------------------------------------*/
/*                                             end of calculation phase */
/*----------------------------------------------------------------------*/
default:
   dserror("Unknown option for solver call to Colsol");
break;   
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of solver_colsol */




