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
case ccf:/*------------------ system matrix is compressed column matrix */
   solver_umfpack(actsolv,actintra,sysarray->ccf,sol,rhs,option);
break;
case skymatrix:/*---------------------- system matrix is skyline matrix */
   solver_colsol(actsolv,actintra,sysarray->sky,sol,rhs,option);
break;
case spoolmatrix:/*-------------------- system matrix is spooles matrix */
   solver_spooles(actsolv,actintra,sysarray->spo,sol,rhs,option);
break;
case bdcsr:/*------------------------------ system matrix is csr matrix */
   solver_mlpcg(actsolv,actintra,sysarray->bdcsr,sol,rhs,option);
break;
case oll:/*------------------------------ system matrix is csr matrix */
   solver_oll(actsolv,actintra,sysarray->oll,sol,rhs,option);
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




