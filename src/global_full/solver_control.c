#include "../headers/standardtypes.h"
#include "../headers/solution.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure allfiles, which holds all file pointers                    |
 | is defined in input_control_global.c
 *----------------------------------------------------------------------*/
extern struct _FILES  allfiles;
/*----------------------------------------------------------------------*
 | define the global structure solv                                     |
 |                                                                      |
 | global variable *solv, vector of lenght numfld of structures SOLVAR  |
 |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
 struct _SOLVAR  *solv; 
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
double t0,t1;
#ifdef DEBUG 
dstrc_enter("solver_control");
#endif

/*----------------------------------------------------------------------*/
t0=ds_cputime();
/*----------------------------------------------------------------------*/
switch(*sysarray_typ)
{
case mds:/*-------------------------------- system matrix is msr matrix */
   solver_mlib(actsolv,actintra,sysarray->mds,sol,rhs,option);
break;
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
case skymatrix:/*---------------------- system matrix is skyline matrix */
   solver_colsol(actsolv,actintra,sysarray->sky,sol,rhs,option);
break;
case spoolmatrix:/*-------------------- system matrix is skyline matrix */
   solver_spooles(actsolv,actintra,sysarray->spo,sol,rhs,option);
break;
default:
   dserror("Unknown format typ of system matrix");
break;   
}
/*----------------------------------------------------------------------*/
t1=ds_cputime();
t1 -= t0;
fprintf(allfiles.out_err,"Time for this solver call: %f\n",t1);
fflush(allfiles.out_err);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of solver_control */




