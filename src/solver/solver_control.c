#include "../headers/standardtypes.h"
#include "../headers/solution.h"
/*----------------------------------------------------------------------*
 |  routine to control all solver calls                  m.gee 9/01     |
 *----------------------------------------------------------------------*/
void solver_control(
                       struct _SOLVAR         *actsolv,
                       struct _INTRA          *actintra,
                       enum   _SPARSE_TYP     *sysarray_typ,
                       union  _SPARSE_ARRAY   *sysarray,
                       struct _DIST_VECTOR    *sol,
                       struct _DIST_VECTOR    *rhs,
                       int                     option
                      )
{
#ifdef DEBUG 
dstrc_enter("solver_control");
#endif

/*----------------------------------------------------------------------*/
switch(*sysarray_typ)
{
case msr:/*-------------------------------- system matrix is msr matrix */
   solver_az_msr(actsolv,actintra,sysarray->msr,sol,rhs,option);
break;
case parcsr:/*-------------------------- system matrix is parcsr matrix */
   solver_hypre_parcsr(actsolv,actintra,sysarray->parcsr,sol,rhs,option);
break;
case ucchb:/*---------------------------- system matrix is ucchb matrix */
   solver_psuperlu_ucchb(actsolv,actintra,sysarray->ucchb,sol,rhs,option);
break;
case dense:/*---------------------------- system matrix is dense matrix */
   solver_lapack_dense(actsolv,actintra,sysarray->dense,sol,rhs,option);
break;
case rc_ptr:/*---------------------- system matrix is row/column matrix */
   solver_mumps(actsolv,actintra,sysarray->rc_ptr,sol,rhs,option);
break;
default:
   dserror("Unknown format typ of system matrix");
break;   
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of solver_control */




