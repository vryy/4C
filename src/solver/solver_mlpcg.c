#ifdef MLPCG

/*!---------------------------------------------------------------------
\file
\brief contains the multilevel cg solver for shells

<pre>
Maintainer: Michael Gee
            gee@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/gee/
            0771 - 685-6572
</pre>

---------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"

/*! 
\addtogroup MLPCG 
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief the multilevel preconditioner main structure

<pre>                                                         m.gee 09/02    
defined in solver_mlpcg.c
</pre>

*----------------------------------------------------------------------*/
struct _MLPRECOND mlprecond;
/*!----------------------------------------------------------------------
\brief the multilevel preconditioned solver main structure

<pre>                                                         m.gee 09/02    
defined in solver_mlpcg.c
</pre>

*----------------------------------------------------------------------*/
struct _MLSOLVER mlsolver;

/*!---------------------------------------------------------------------
\brief multilevel preconditioned iterative solver                                              

<pre>                                                        m.gee 9/02 
</pre>
\param actsolv    SOLVAR*      (i)   general structure of solver informations                   
\param actintra   INTRA*       (i)   the intra-communicator of this field                  
\param bdcsr      DBCSR_ROOT*  (i)   the dbcsr matrix                 
\param sol        DIST_VECTOR* (o)   the distributed solution vector
\param rhs        DIST_VECTOR* (i)   the distributed right hand side vector
\param option     INT          (i)   option=1 init phase option=0 solve
\return void                                               

------------------------------------------------------------------------*/
void solver_mlpcg(struct _SOLVAR         *actsolv,
                  struct _INTRA          *actintra,
                  struct _DBCSR          *bdcsr,
                  struct _DIST_VECTOR    *sol,
                  struct _DIST_VECTOR    *rhs,
                  INT                     option)
{
MLPCGVARS           *mlpcgvars;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("solver_mlpcg");
#endif
/*----------------------------------------------------------------------*/
mlpcgvars = actsolv->mlpcgvars;
/*----------------------------------------------------------------------*/
switch(option)
{
/*----------------------------------------------------------------------*/
/*                                                           init phase */
/*----------------------------------------------------------------------*/
case 1:
   mlpcg_precond_create(bdcsr,mlpcgvars,actintra);
   mlpcg_solver_create(bdcsr,sol,rhs,mlpcgvars);
   /*------------------------------------- set initial values to matrix */
   mlprecond.ncall    = 0;
   mlprecond.mod      = 0;
   bdcsr->ncall       = 0;
   bdcsr->is_factored = 0;
   bdcsr->is_init     = 1;
break;
/*----------------------------------------------------------------------*/
/*                                                    end of init phase */
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*                                                    calculation phase */
/*----------------------------------------------------------------------*/
case 0:
   if (mlprecond.reuse)
   mlprecond.mod = mlprecond.ncall % mlprecond.reuse;
   else
   mlprecond.mod = 0;
   /*----------------------------- create the multilevel preconditioner */
   mlpcg_precond_init(bdcsr,mlpcgvars,actintra);
   /*-------------------------------------------------- init the solver */
   mlpcg_solver_init(bdcsr,sol,rhs,actintra);
   /*-------------------------------------------------- call the solver */
   mlpcg_pcg(bdcsr,sol,rhs,actintra);
   /*--------------------------------- increment number of solver calls */
   mlprecond.ncall++;
   bdcsr->is_factored = 1;
break;
/*----------------------------------------------------------------------*/
/*                                             end of calculation phase */
/*----------------------------------------------------------------------*/
default:
   dserror("Unknown option for solver call to MLPCG");
break;   
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of solver_mlpcg */


/*! @} (documentation module close)*/
#endif

