#include "../headers/standardtypes.h"
#include "../headers/solution.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure allfiles, which holds all file pointers                    |
 | is defined in input_control_global.c                                 |
 | struct _FILES  allfiles;                                             |
 *----------------------------------------------------------------------*/
extern struct _FILES  allfiles; 
/*--------------------------------------------------------------------------*
 |  control solver lib HYPRE                            m.gee 10/01         |
 |  struct _SOLVAR         *actsolv    ptr to actual solver structure       |
 |  struct _INTRA          *actintra   ptr to field's intra-communicator    |
 |  struct _H_PARCSR       *parcsr     ptr to parcsr system matrix to solve | 
 |  struct _DIST_VECTOR    *sol        ptr to distrib. solution vector      |
 |  struct _DIST_VECTOR    *rhs        ptr to distrib. rhs vector           |
 |  int                     option     =1 init phase                        |
 |                                     =0 calculation phase                 |
 | this routine uses HYPRE V1.6.0 - high performance preconditioner library |
 | the HYPRE header files are included in solution.h                        |
 | for further information see                                              |
 | http://www.llnl.gov/CASC/hypre/                                          |
 | HYPRE User's Manual                                                      |
 | HYPRE Reference Manual                                                   |
 | NOTE: This package is still under extreme construction                   |
 | NOTE2: the solver and preconditioner objects are never destroyed         |
 | This could be reason for growth of memory
 *--------------------------------------------------------------------------*/
void  solver_hypre_parcsr( 
                            struct _SOLVAR         *actsolv,
                            struct _INTRA          *actintra,
                            struct _H_PARCSR       *parcsr,
                            struct _DIST_VECTOR    *sol,
                            struct _DIST_VECTOR    *rhs,
                            int                     option
                           )
{
#ifdef HYPRE_PACKAGE
int                 i;                     /* counter variable */
int                 err;                   /* error flag */

int                 myrank;                /* this processors intra-rank */
int                 nprocs;                /* number of procs in the intra-communicator */

HYPREVARS          *hyprevars;             /* variables for hypre solver read from input file */
SOLVER_TYP         *solvertyp;             /* enum, type of hypre solver */

HYPRE_IJVector      ijvector_rhs;          /* right-hand-side vector structure */
HYPRE_IJVector      ijvector_sol;          /* solution        vector structure */

HYPRE_ParCSRMatrix  parcsr_matrix;         /* parallel compressed sparse row matrix structure */
HYPRE_ParVector     parcsr_vector_rhs;     /* vector that suits the HYPRE_ParCSRMatrix matrix */
HYPRE_ParVector     parcsr_vector_sol;     /* vector that suits the HYPRE_ParCSRMatrix matrix */

#ifdef DEBUG 
dstrc_enter("solver_hypre_parcsr");
#endif
/*----------------------------------------------------------------------*/
myrank = actintra->intra_rank;
nprocs = actintra->intra_nprocs;
/*------------------------ There are lots of solvers in HYPRE, get data */
solvertyp = &(actsolv->solvertyp);
hyprevars = actsolv->hyprevar;
/*----------------------------------------------------------------------*/
switch(option)
{
/*======================================================================*/
/*                                                           init phase */
/*======================================================================*/
case 1:
   err=0;
   /*----------------------------------- create the ij_matrix structure */
   hypre_matrix_create(myrank,parcsr,actintra,&(parcsr->ij_matrix));
   /*---------------------------- set init flag to the parcsr structure */
   parcsr->is_factored=0;
   parcsr->is_init=1;
   parcsr->ncall=0;
   /*-------------------------------------------------------------------*/
   /*---------------------------------------- switch for type of solver */
   switch (*solvertyp)
   {

   /*===================================================================*/
   /*                     parallel algebraic multigrid solver BoomerAMG */
   /*===================================================================*/
   case hypre_amg:
      err=0;
      /*----------------------------------------------- create a solver */
      err=HYPRE_BoomerAMGCreate(&(parcsr->solver));
      if (err) dserror("error in calling solver package HYPRE");
      /*-------------------------------- set additional parameters here */
      hypre_set_params_boomeramg_solver(&(parcsr->solver),hyprevars);
   break;
   /*================================================= end of BoomerAMG */

   /*===================================================================*/
   /*                   parallel GMRES with different preconditioners   */
   /*===================================================================*/
   case hypre_gmres:
      /*----------------------------------------------- create a solver */
      err=HYPRE_ParCSRGMRESCreate(actintra->MPI_INTRA_COMM,&(parcsr->solver));
      if (err) dserror("error in calling solver package HYPRE");
      /*-------------------------------- set additional parameters here */
      hypre_set_params_gmres_solver(&(parcsr->solver),hyprevars);
      /*----------------------------------------- set up preconditioner */
      switch(hyprevars->hypre_prectyp)
      {
      /*--------------------------------------------- no preconditioner */
      case hypreprec_none:
      break;
      /*------------------- use incomplete factorisation package Euclid */
      case hypreprec_euclid:
         /*-------------------------------------- create preconditioner */
         hypre_create_precond_euclid(&(parcsr->precond),actintra,hyprevars);
         /*--------------------------- pass function pointers to solver */     
         err=HYPRE_GMRESSetPrecond(parcsr->solver,
                                   (HYPRE_PtrToSolverFcn) HYPRE_EuclidSolve,
                                   (HYPRE_PtrToSolverFcn) HYPRE_EuclidSetup,
                                   parcsr->precond);
         if (err) dserror("error in calling solver package HYPRE");
         /*-------------------------------------------------------------*/
      break;
      /*------ use approximate inverse preconditioner package parasails */
      case hypreprec_parasails:
         /*-------------------------------------- create preconditioner */
         hypre_create_precond_parasails(&(parcsr->precond),actintra,hyprevars);
         /*--------------------------- pass function pointers to solver */ 
         err+=HYPRE_GMRESSetPrecond(parcsr->solver,
                                    (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSolve,
                                    (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSetup,
                                    parcsr->precond);
         /*-------------------------------------------------------------*/
         if (err) dserror("error in calling solver package HYPRE");
         /*-------------------------------------------------------------*/
      break;
      /*-------------- use algebraic multigrid preconditioner BoomerAMG */
      case hypreprec_amg:
         /*-------------------------------------- create preconditioner */
         hypre_create_precond_boomeramg(&(parcsr->precond),actintra,hyprevars);
         /*--------------------------- pass function pointers to solver */ 
         err=HYPRE_GMRESSetPrecond(parcsr->solver,
                                   (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
                                   (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup,
                                   parcsr->precond);
         /*-------------------------------------------------------------*/
         if (err) dserror("error in calling solver package HYPRE");
         /*-------------------------------------------------------------*/
      break;
      default:
         dserror("Unknown type of preconditioner for HYPRE solver");
      }
      /*----------------------------------------------------------------*/
   break;
   /*===================================================== end of GMRES */

   /*===================================================================*/
   /*                parallel BiCGStab with different preconditioners   */
   /*===================================================================*/
   case hypre_bicgstab:
      /*----------------------------------------------- create a solver */
      err=HYPRE_ParCSRBiCGSTABCreate(actintra->MPI_INTRA_COMM,&(parcsr->solver));
      if (err) dserror("error in calling solver package HYPRE");
      /*------------------------------------ set additional parameters  */
      hypre_set_params_bicgstab_solver(&(parcsr->solver),hyprevars);
      /*----------------------------------------- set up preconditioner */
      switch(hyprevars->hypre_prectyp)
      {
      /*--------------------------------------------- no preconditioner */
      case hypreprec_none:
      break;
      /*------------------- use incomplete factorisation package Euclid */
      case hypreprec_euclid:
         /*-------------------------------------- create preconditioner */
         hypre_create_precond_euclid(&(parcsr->precond),actintra,hyprevars);
         /*--------------------------- pass function pointers to solver */     
         err=HYPRE_BiCGSTABSetPrecond(parcsr->solver,
                                      (HYPRE_PtrToSolverFcn) HYPRE_EuclidSolve,
                                      (HYPRE_PtrToSolverFcn) HYPRE_EuclidSetup,
                                      parcsr->precond);
         /*-------------------------------------------------------------*/
         if (err) dserror("error in calling solver package HYPRE");
         /*-------------------------------------------------------------*/
      break;
      /*------ use approximate inverse preconditioner package parasails */
      case hypreprec_parasails:
         /*-------------------------------------- create preconditioner */
         hypre_create_precond_parasails(&(parcsr->precond),actintra,hyprevars);
         /*--------------------------- pass function pointers to solver */ 
         err=HYPRE_BiCGSTABSetPrecond(parcsr->solver,
                                     (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSolve,
                                     (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSetup,
                                     parcsr->precond);
         /*-------------------------------------------------------------*/
         if (err) dserror("error in calling solver package HYPRE");
         /*-------------------------------------------------------------*/
      break;
      /*-------------- use algebraic multigrid preconditioner BoomerAMG */
      case hypreprec_amg:
         /*-------------------------------------- create preconditioner */
         hypre_create_precond_boomeramg(&(parcsr->precond),actintra,hyprevars);
         /*--------------------------- pass function pointers to solver */ 
         err=HYPRE_BiCGSTABSetPrecond(parcsr->solver,
                                      (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
                                      (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup,
                                      parcsr->precond);
         /*-------------------------------------------------------------*/
         if (err) dserror("error in calling solver package HYPRE");
         /*-------------------------------------------------------------*/
      break;
      default:
         dserror("Unknown type of preconditioner for HYPRE-CG solver");
      }
   break;
   /*================================================== end of BiCGStab */

   /*===================================================================*/
   /*                      parallel CG with different preconditioners   */
   /*===================================================================*/
   case hypre_pcg:
      err=0;
      /*----------------------------------------------- create a solver */
      err=HYPRE_ParCSRPCGCreate(actintra->MPI_INTRA_COMM,&(parcsr->solver));
      if (err) dserror("error in calling solver package HYPRE");
      /*-------------------------------- set additional parameters here */
      hypre_set_params_pcg_solver(&(parcsr->solver),hyprevars);
      /*----------------------------------------- set up preconditioner */
      switch(hyprevars->hypre_prectyp)
      {
      /*--------------------------------------------- no preconditioner */
      case hypreprec_none:
      break;
      /*------------------- use incomplete factorisation package Euclid */
      case hypreprec_euclid:
         /*-------------------------------------- create preconditioner */
         hypre_create_precond_euclid(&(parcsr->precond),actintra,hyprevars);
         /*--------------------------- pass function pointers to solver */     
         err=HYPRE_PCGSetPrecond(parcsr->solver,
                                 (HYPRE_PtrToSolverFcn) HYPRE_EuclidSolve,
                                 (HYPRE_PtrToSolverFcn) HYPRE_EuclidSetup,
                                 parcsr->precond);
         if (err) dserror("error in calling solver package HYPRE");
         /*-------------------------------------------------------------*/
      break;
      /*------ use approximate inverse preconditioner package parasails */
      case hypreprec_parasails:
         /*-------------------------------------- create preconditioner */
         hypre_create_precond_parasails(&(parcsr->precond),actintra,hyprevars);
         /*--------------------------- pass function pointers to solver */ 
         err=HYPRE_PCGSetPrecond(parcsr->solver,
                                 (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSolve,
                                 (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSetup,
                                 parcsr->precond);
         if (err) dserror("error in calling solver package HYPRE");
         /*-------------------------------------------------------------*/
      break;
      /*-------------- use algebraic multigrid preconditioner BoomerAMG */
      case hypreprec_amg:
         /*-------------------------------------- create preconditioner */
         hypre_create_precond_boomeramg(&(parcsr->precond),actintra,hyprevars);
         /*--------------------------- pass function pointers to solver */ 
         err=HYPRE_PCGSetPrecond(parcsr->solver,
                                 (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
                                 (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup,
                                 parcsr->precond);
         if (err) dserror("error in calling solver package HYPRE");
         /*-------------------------------------------------------------*/
      break;
      default:
         dserror("Unknown type of preconditioner for HYPRE-CG solver");
      }
   break;
   /*======================================================== end of CG */

   /*===================================================================*/
   /*                                                         default   */
   /*===================================================================*/
   default:
      dserror("Unknown typ of HYPRE solver");
   break;
   }/*================================ end of switch over HYPRE solvers */
break;
/*======================================================================*/
/*                                                    end of init phase */
/*======================================================================*/










/*======================================================================*/
/*                                                    calculation phase */
/*======================================================================*/
case 0:
   /*--------------------------------------------------- check for init */
   if (parcsr->is_init!=1)
      dserror("solver package HYPRE has not been initialized");
   /*-------------------------------- create a HYPRE vector for the rhs */
   hypre_vector_create(myrank,parcsr,actintra,&ijvector_rhs);
   /*-------------------------------- create a HYPRE vector for the sol */
   hypre_vector_create(myrank,parcsr,actintra,&ijvector_sol);
   /*--------------------------------------put values to the rhs vector */
   err=HYPRE_IJVectorSetValues(ijvector_rhs,rhs->numeq,parcsr->perm.a.ia[myrank],rhs->vec.a.dv);
   if (err) dserror("error in setting up solver package HYPRE");                           
   /*-------------------------------- finalize rhs and sol for solution */
   hypre_vector_assemble(&ijvector_rhs,&parcsr_vector_rhs);
   hypre_vector_assemble(&ijvector_sol,&parcsr_vector_sol);
   /*-------------------------------------- make matrix ready for solve */
   hypre_matrix_assemble(&(parcsr->ij_matrix),&(parcsr_matrix));
   /*-------------------------------------------------------------------*/
   /*---------------------------------------- switch for type of solver */
   switch (*solvertyp)
   {


   /*===================================================================*/
   /*                     parallel algebraic multigrid solver BoomerAMG */
   /*===================================================================*/
   case hypre_amg:
      err=0;
      /*-------------------------------------------------- setup solver */
      err+=HYPRE_BoomerAMGSetup(parcsr->solver,parcsr_matrix,parcsr_vector_rhs,parcsr_vector_sol);
      if (err) dserror("error in calling solver package HYPRE");     
      /*-------------------------------------------- Yes, finally solve */                      
      err=HYPRE_BoomerAMGSolve(parcsr->solver,parcsr_matrix,parcsr_vector_rhs,parcsr_vector_sol);
      if (err && myrank==0) 
      {
         printf("RANK %d: error or noconvergence in calling solver BoomerAMG\n",myrank);     
         fprintf(allfiles.out_err,"RANK %d: error or noconvergence in calling solver BoomerAMG\n",myrank); 
         fflush(stdout);
      }
      /*---------- get number of iterations and final relative residual */
      err=0;
      if (hyprevars->io)
      {
         err+=HYPRE_BoomerAMGGetNumIterations(parcsr->solver,&(hyprevars->numiter));
         err+=HYPRE_BoomerAMGGetFinalRelativeResidualNorm(parcsr->solver,&(hyprevars->resnorm));
         if (myrank==0)
         printf("BoomerAMG: Numiterations=%d ResidualNorm=%E\n",hyprevars->numiter,hyprevars->resnorm);
         fprintf(allfiles.out_err,"ParCG: Numiterations=%d ResidualNorm=%E\n",hyprevars->numiter,hyprevars->resnorm);
      }
      /*----------------------------------------------------------------*/
   break;
   /*================================================= end of BoomerAMG */





   /*===================================================================*/
   /*                   parallel GMRES with different preconditioners   */
   /*===================================================================*/
   case hypre_gmres:
      err=0;
      /*-------------------------------------------------- set up reuse */
      switch(hyprevars->hypre_prectyp)
      {
      case hypreprec_none:
         err=HYPRE_ParCSRGMRESSetup(parcsr->solver,parcsr_matrix,parcsr_vector_rhs,
                                    parcsr_vector_sol);
      break;
      case hypreprec_euclid:
         if (parcsr->is_factored==0)
         err=HYPRE_ParCSRGMRESSetup(parcsr->solver,parcsr_matrix,parcsr_vector_rhs,
                                    parcsr_vector_sol);
      break;
      case hypreprec_parasails:
         if (parcsr->is_factored) HYPRE_ParaSailsSetReuse(parcsr->precond,1);
         else                     HYPRE_ParaSailsSetReuse(parcsr->precond,0);
         err=HYPRE_ParCSRGMRESSetup(parcsr->solver,parcsr_matrix,parcsr_vector_rhs,
                                    parcsr_vector_sol);
      break;
      case hypreprec_amg:
         if (parcsr->is_factored==0)
         err=HYPRE_ParCSRGMRESSetup(parcsr->solver,parcsr_matrix,parcsr_vector_rhs,
                                    parcsr_vector_sol);
      break;
      default:
         dserror("Unknown type of preconditioner for HYPRE solver");
      }
      /*----------------------------------------------------------------*/
      if (err) dserror("error in calling solver package HYPRE");     
      /*-------------------------------------------- Yes, finally solve */                      
      err=HYPRE_ParCSRGMRESSolve(parcsr->solver,parcsr_matrix,parcsr_vector_rhs,
                                 parcsr_vector_sol);
      if (err) 
      {
         printf("RANK %d: error or noconvergence in calling solver GMRES\n",myrank);     
         fprintf(allfiles.out_err,"RANK %d: error or noconvergence in calling solver EuclidCG\n",myrank); 
         fflush(stdout);
      }
      /*---------- get number of iterations and final relative residual */
      err=0;
      if (hyprevars->io)
      {
         err+=HYPRE_ParCSRGMRESGetNumIterations(parcsr->solver,&(hyprevars->numiter));
         err+=HYPRE_ParCSRGMRESGetFinalRelativeResidualNorm(parcsr->solver,&(hyprevars->resnorm));
         if (myrank==0)
         printf("GMRES: Numiterations=%d ResidualNorm=%E\n",hyprevars->numiter,hyprevars->resnorm);
         fprintf(allfiles.out_err,"ParCG: Numiterations=%d ResidualNorm=%E\n",hyprevars->numiter,hyprevars->resnorm);
      }
      /*----------------------------------------------------------------*/
      if (err) dserror("error in calling solver package HYPRE");     
   break;
   /*===================================================== end of GMRES */








   /*===================================================================*/
   /*                parallel BiCGStab with different preconditioners   */
   /*===================================================================*/
   case hypre_bicgstab:
      err=0;
      /*-------------------------------------------------- set up reuse */
      switch(hyprevars->hypre_prectyp)
      {
      case hypreprec_none:
      err=HYPRE_BiCGSTABSetup(parcsr->solver,
                              (HYPRE_Matrix)parcsr_matrix,
                              (HYPRE_Vector)parcsr_vector_rhs,
                              (HYPRE_Vector)parcsr_vector_sol);
      break;
      case hypreprec_euclid:
         if (parcsr->is_factored==0)
      err=HYPRE_BiCGSTABSetup(parcsr->solver,
                              (HYPRE_Matrix)parcsr_matrix,
                              (HYPRE_Vector)parcsr_vector_rhs,
                              (HYPRE_Vector)parcsr_vector_sol);
      break;
      case hypreprec_parasails:
         if (parcsr->is_factored) HYPRE_ParaSailsSetReuse(parcsr->precond,1);
         else                     HYPRE_ParaSailsSetReuse(parcsr->precond,0);
      err=HYPRE_BiCGSTABSetup(parcsr->solver,
                              (HYPRE_Matrix)parcsr_matrix,
                              (HYPRE_Vector)parcsr_vector_rhs,
                              (HYPRE_Vector)parcsr_vector_sol);
      break;
      case hypreprec_amg:
         if (parcsr->is_factored==0)
      err=HYPRE_BiCGSTABSetup(parcsr->solver,
                              (HYPRE_Matrix)parcsr_matrix,
                              (HYPRE_Vector)parcsr_vector_rhs,
                              (HYPRE_Vector)parcsr_vector_sol);
      break;
      default:
         dserror("Unknown type of preconditioner for HYPRE solver");
      }
      /*----------------------------------------------------------------*/
      if (err) dserror("error in calling solver package HYPRE");     
      /*-------------------------------------------- Yes, finally solve */                      
      err=HYPRE_BiCGSTABSolve(parcsr->solver,
                              (HYPRE_Matrix)parcsr_matrix,
                              (HYPRE_Vector)parcsr_vector_rhs,
                              (HYPRE_Vector)parcsr_vector_sol);
      if (err) 
      {
         printf("RANK %d: error or noconvergence in calling solver BiCGStab\n",myrank);     
         fprintf(allfiles.out_err,"RANK %d: error or noconvergence in calling solver EuclidCG\n",myrank); 
         fflush(stdout);
      }
      /*---------- get number of iterations and final relative residual */
      err=0;
      if (hyprevars->io)
      {
         err+=HYPRE_BiCGSTABGetNumIterations(parcsr->solver,&(hyprevars->numiter));
         err+=HYPRE_BiCGSTABGetFinalRelativeResidualNorm(parcsr->solver,&(hyprevars->resnorm));
         if (myrank==0)
         printf("BiCGSTAB: Numiterations=%d ResidualNorm=%E\n",hyprevars->numiter,hyprevars->resnorm);
         fprintf(allfiles.out_err,"ParCG: Numiterations=%d ResidualNorm=%E\n",hyprevars->numiter,hyprevars->resnorm);
      }
      /*----------------------------------------------------------------*/
      if (err) dserror("error in calling solver package HYPRE");     
   break;
   /*================================================== end of BiCGStab */







   /*===================================================================*/
   /*                      parallel CG with different preconditioners   */
   /*===================================================================*/
   case hypre_pcg:
      err=0;
      /*-------------------------------------------------- set up reuse */
      switch(hyprevars->hypre_prectyp)
      {
      case hypreprec_none:
         err=HYPRE_ParCSRPCGSetup(parcsr->solver,parcsr_matrix,parcsr_vector_rhs,
                                  parcsr_vector_sol);
      break;
      case hypreprec_euclid:
         if (parcsr->is_factored==0)
         err=HYPRE_ParCSRPCGSetup(parcsr->solver,parcsr_matrix,parcsr_vector_rhs,
                                  parcsr_vector_sol);
      break;
      case hypreprec_parasails:
         if (parcsr->is_factored) HYPRE_ParaSailsSetReuse(parcsr->precond,1);
         else                     HYPRE_ParaSailsSetReuse(parcsr->precond,0);
         err=HYPRE_ParCSRPCGSetup(parcsr->solver,parcsr_matrix,parcsr_vector_rhs,
                                  parcsr_vector_sol);
      break;
      case hypreprec_amg:
         if (parcsr->is_factored==0)
         err=HYPRE_ParCSRPCGSetup(parcsr->solver,parcsr_matrix,parcsr_vector_rhs,
                                  parcsr_vector_sol);
      break;
      default:
         dserror("Unknown type of preconditioner for HYPRE solver");
      }
      /*----------------------------------------------------------------*/
      if (err) dserror("error in calling solver package HYPRE");     
      /*-------------------------------------------- Yes, finally solve */                      
      err=HYPRE_ParCSRPCGSolve(parcsr->solver,parcsr_matrix,parcsr_vector_rhs,
                               parcsr_vector_sol);
      if (err) 
      {
         printf("RANK %d: error or noconvergence in calling solver EuclidCG\n",myrank);     
         fprintf(allfiles.out_err,"RANK %d: error or noconvergence in calling solver EuclidCG\n",myrank); 
         fflush(stdout);
      }
      /*---------- get number of iterations and final relative residual */
      err=0;
      if (hyprevars->io)
      {
         err+=HYPRE_ParCSRPCGGetNumIterations(parcsr->solver,&(hyprevars->numiter));
         err+=HYPRE_ParCSRPCGGetFinalRelativeResidualNorm(parcsr->solver,&(hyprevars->resnorm));
         if (myrank==0)
         printf("ParCG: Numiterations=%d ResidualNorm=%E\n",hyprevars->numiter,hyprevars->resnorm);
         fprintf(allfiles.out_err,"ParCG: Numiterations=%d ResidualNorm=%E\n",hyprevars->numiter,hyprevars->resnorm);
      }
      /*----------------------------------------------------------------*/
      if (err) dserror("error in calling solver package HYPRE");     
   break;
   /*======================================================== end of CG */





   /*===================================================================*/
   /*                                                         default   */
   /*===================================================================*/
   default:
      dserror("Unknown typ of HYPRE solver");
   break;
   }/*================================ end of switch over HYPRE solvers */








   /*===================================================================*/
   /*                                          tidy up and get solution */
   /*===================================================================*/
   /*------------------------------------------- destroy the rhs object */
   err=HYPRE_IJVectorDestroy(ijvector_rhs);
   /*------------------------------------------------- get the solution */
   err=HYPRE_IJVectorGetValues(
                               ijvector_sol,
                               sol->numeq,
                               parcsr->perm.a.ia[myrank],
                               sol->vec.a.dv 
                              );
    /*------------------------------------- destroy the solution object */
   err=HYPRE_IJVectorDestroy(ijvector_sol);
   /*---------------------------------------- do not destroy the matrix */
   /* NOTE: before adding to the matrix parcsr->ij_matrix again,
            there has to be a call to HYPRE_IJMatrixInitialize.
            This call is done by the routine solserv_zero_mat
            in solver_service.c
   /*-------------------------------------------------------- set flags */
   parcsr->is_factored=1;
   parcsr->ncall++;
   /*-------------------------------------------------------------------*/
break;
/*======================================================================*/
/*                                             end of calculation phase */
/*======================================================================*/
default:
   dserror("Unknown option for solver call to HYPRE");
break;   
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
#endif /* end of ifdef HYPRE_PACKAGE */
return;
} /* end of solver_hypre_parcsr */







#ifdef HYPRE_PACKAGE
/*----------------------------------------------------------------------*
 |  create a vector for HYPRE Package                        m.gee 10/01|
 *----------------------------------------------------------------------*/
void hypre_vector_create(int              myrank,
                           H_PARCSR        *parcsr,
                           INTRA           *actintra,
                           HYPRE_IJVector  *ijvector)
{
int err=0;
int ilower,iupper;
int jlower,jupper;
#ifdef DEBUG 
dstrc_enter("hypre_vector_create");
#endif
/*----------------------------------------------------------------------*/
/*----------------------------- check distribution of permutated matrix */
ilower = jlower = parcsr->perm.a.ia[myrank][0];
iupper = jupper = parcsr->perm.a.ia[myrank][parcsr->perm_sizes.a.iv[myrank]-1];
/*--------------------------------------------------- create the vector */
err+=HYPRE_IJVectorCreate(actintra->MPI_INTRA_COMM,ilower,iupper,ijvector);
err+=HYPRE_IJVectorSetObjectType(*ijvector,HYPRE_PARCSR);
err+=HYPRE_IJVectorInitialize(*ijvector);
if (err) dserror("Creation of HYPRE vector failed");
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of hypre_vector_create */
#endif /* end of ifdef HYPRE_PACKAGE */

#ifdef HYPRE_PACKAGE
/*----------------------------------------------------------------------*
 |  assemble a HYPRE vector and make it ready for solution   m.gee 10/01|
 *----------------------------------------------------------------------*/
void hypre_vector_assemble(HYPRE_IJVector  *ijvector,
                             HYPRE_ParVector *parcsr_vector)
{
int err=0;
#ifdef DEBUG 
dstrc_enter("hypre_vector_assemble");
#endif
/*----------------------------------------------------------------------*/
err+=HYPRE_IJVectorAssemble(*ijvector);
err+=HYPRE_IJVectorGetObject(*ijvector,(void**)parcsr_vector);
if (err) dserror("Assembly of HYPRE vector failed");
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of hypre_vector_assemble */
#endif /* end of ifdef HYPRE_PACKAGE */

#ifdef HYPRE_PACKAGE
/*----------------------------------------------------------------------*
 |  create a vector for HYPRE Package                        m.gee 10/01|
 *----------------------------------------------------------------------*/
void hypre_matrix_create(int              myrank,
                           H_PARCSR        *parcsr,
                           INTRA           *actintra,
                           HYPRE_IJMatrix  *ijmatrix)
{
int err=0;
int ilower,iupper;
int jlower,jupper;
#ifdef DEBUG 
dstrc_enter("hypre_matrix_create");
#endif
/*----------------------------------------------------------------------*/
/*----------------------------- check distribution of permutated matrix */
ilower = jlower = parcsr->perm.a.ia[myrank][0];
iupper = jupper = parcsr->perm.a.ia[myrank][parcsr->perm_sizes.a.iv[myrank]-1];
/*--------------------------------------------------- create the matrix */
err+=HYPRE_IJMatrixCreate(actintra->MPI_INTRA_COMM,ilower,iupper,jlower,jupper,ijmatrix);
/*-------------------------------------------------- set the object typ */
err+=HYPRE_IJMatrixSetObjectType(*ijmatrix,HYPRE_PARCSR);
/*----------------------------------------------------- check for error */
if (err) dserror("Creation of HYPRE matrix failed");
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of hypre_matrix_create */
#endif /* end of ifdef HYPRE_PACKAGE */

#ifdef HYPRE_PACKAGE
/*----------------------------------------------------------------------*
 |  assemble a HYPRE matrix and make it ready for solution   m.gee 10/01|
 *----------------------------------------------------------------------*/
void hypre_matrix_assemble(HYPRE_IJMatrix     *ij_matrix,
                             HYPRE_ParCSRMatrix *parcsr_matrix)
{
int err=0;
#ifdef DEBUG 
dstrc_enter("hypre_matrix_assemble");
#endif
/*----------------------------------------------------------------------*/
err+=HYPRE_IJMatrixAssemble(*ij_matrix);
err+=HYPRE_IJMatrixGetObject(*ij_matrix,(void**)parcsr_matrix);
if (err) dserror("Assembly of HYPRE matrix failed");
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of hypre_matrix_assemble */
#endif /* end of ifdef HYPRE_PACKAGE */


#ifdef HYPRE_PACKAGE
/*----------------------------------------------------------------------*
 |  set parameters for solver boomeramg                      m.gee 10/01|
 *----------------------------------------------------------------------*/
void hypre_set_params_boomeramg_solver(HYPRE_Solver *solver,HYPREVARS *hyprevars)
{
int i;
int err=0;
int gridrelaxtype[4];
int gridnsweep[4];
#ifdef DEBUG 
dstrc_enter("hypre_solver_set_params_boomeramg");
#endif
/*----------------------------------------------------------------------*/
/*---------------------------------------------------------- set output */
if (hyprevars->io)
err+=HYPRE_BoomerAMGSetIOutDat(*solver,3);
else
err+=HYPRE_BoomerAMGSetIOutDat(*solver,0);
/*---------------------------------- set maximum number of grids to use */
err+=HYPRE_BoomerAMGSetMaxLevels(*solver,500);
/*-------------------------------------------------------- set max iter */
err+=HYPRE_BoomerAMGSetMaxIter(*solver,hyprevars->maxiter);
/*------------------------------------------------- set type of measure */
err+=HYPRE_BoomerAMGSetMeasureType(*solver,1);
/*------------------------------------------------------------- set tol */
err+=HYPRE_BoomerAMGSetTol(*solver,hyprevars->tol);
/*-------------------- set strong connection (the smaller the stronger) */
err+=HYPRE_BoomerAMGSetStrongThreshold(*solver,hyprevars->threshold);
/*------------------------------------------------- set max row sum off */
err+=HYPRE_BoomerAMGSetMaxRowSum(*solver,1.0);
/*------------------------------------------------- set good relaxation */
gridrelaxtype[0]=6;
gridrelaxtype[1]=6;
gridrelaxtype[2]=6;
gridrelaxtype[3]=9;
err+=HYPRE_BoomerAMGSetGridRelaxType(*solver,gridrelaxtype);
/*--------------------------------------------------------- set w cycle */
err+=HYPRE_BoomerAMGSetCycleType(*solver,2);
/*--------------------------------------- set number of sweeps per grid */
for (i=0; i<4; i++) gridnsweep[i] = hyprevars->sweep[i];
err+=HYPRE_BoomerAMGSetNumGridSweeps(*solver,gridnsweep);
/*------------------------------------ set type of coarsening algorithm */
err+=HYPRE_BoomerAMGSetCoarsenType(*solver,0);
/*----------------------------------------------------------------------*/
if (err) dserror("error in calling solver package HYPRE");   
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of hypre_solver_set_params_boomeramg */
#endif /* end of ifdef HYPRE_PACKAGE */

#ifdef HYPRE_PACKAGE
/*----------------------------------------------------------------------*
 |  set parameters for solver gmres                          m.gee 10/01|
 *----------------------------------------------------------------------*/
void hypre_set_params_gmres_solver(HYPRE_Solver *solver,HYPREVARS *hyprevars)
{
int i;
int err=0;
#ifdef DEBUG 
dstrc_enter("hypre_set_params_gmres_solver");
#endif
/*----------------------------------------------------------------------*/
/*------------------------------------------------------- set tolerance */
err+=HYPRE_ParCSRGMRESSetTol(*solver,hyprevars->tol);
/*------------------------------------ set maximum number of iterations */
err+=HYPRE_ParCSRGMRESSetMaxIter(*solver,hyprevars->maxiter);
/*------------------------------------------------------------- set i/o */
if (hyprevars->io)
   err+=HYPRE_ParCSRGMRESSetLogging(*solver,3);
else
   err+=HYPRE_ParCSRGMRESSetLogging(*solver,0);
/*--------------------------------------------- set size of krylovspace */
err+=HYPRE_ParCSRGMRESSetKDim(*solver,hyprevars->kryldim);
/*----------------------------------------------------------------------*/
if (err) dserror("error in calling solver package HYPRE");   
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of hypre_set_params_gmres_solver */
#endif /* end of ifdef HYPRE_PACKAGE */

#ifdef HYPRE_PACKAGE
/*----------------------------------------------------------------------*
 |  set parameters for solver bicgstab                       m.gee 10/01|
 *----------------------------------------------------------------------*/
void hypre_set_params_bicgstab_solver(HYPRE_Solver *solver,HYPREVARS *hyprevars)
{
int i;
int err=0;
#ifdef DEBUG 
dstrc_enter("hypre_set_params_bicgstab_solver");
#endif
/*----------------------------------------------------------------------*/
/*------------------------------------------------------- set tolerance */
err+=HYPRE_BiCGSTABSetTol(*solver,hyprevars->tol);
/*------------------------------------ set maximum number of iterations */
err+=HYPRE_BiCGSTABSetMaxIter(*solver,hyprevars->maxiter);
/*------------------------------------------------------------- set i/o */
if (hyprevars->io)
   err+=HYPRE_BiCGSTABSetLogging(*solver,3);
else
   err+=HYPRE_BiCGSTABSetLogging(*solver,0);
/*----------------------------------------------------------------------*/
if (err) dserror("error in calling solver package HYPRE");   
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of hypre_set_params_bicgstab_solver */
#endif /* end of ifdef HYPRE_PACKAGE */

#ifdef HYPRE_PACKAGE
/*----------------------------------------------------------------------*
 |  set parameters for solver pcg                          m.gee 10/01|
 *----------------------------------------------------------------------*/
void hypre_set_params_pcg_solver(HYPRE_Solver *solver,HYPREVARS *hyprevars)
{
int i;
int err=0;
#ifdef DEBUG 
dstrc_enter("hypre_set_params_pcg_solver");
#endif
/*----------------------------------------------------------------------*/
/*------------------------------------------------------- set tolerance */
err+=HYPRE_ParCSRPCGSetTol(*solver,hyprevars->tol);
/*------------------------------------ set maximum number of iterations */
err+=HYPRE_ParCSRPCGSetMaxIter(*solver,hyprevars->maxiter);
/*------------------------------------- use 2-norm in stopping criteria */
err+=HYPRE_ParCSRPCGSetTwoNorm(*solver,1);
/*------------------------------------------------------------- set i/o */
if (hyprevars->io)
   err+=HYPRE_ParCSRPCGSetLogging(*solver,3);
else
   err+=HYPRE_ParCSRPCGSetLogging(*solver,0);
/*----------------------------------------------------------------------*/
if (err) dserror("error in calling solver package HYPRE");   
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of hypre_set_params_pcg_solver */
#endif /* end of ifdef HYPRE_PACKAGE */

#ifdef HYPRE_PACKAGE
/*----------------------------------------------------------------------*
 |  create preconditioner euclid and set parameters          m.gee 10/01|
 *----------------------------------------------------------------------*/
void hypre_create_precond_euclid(HYPRE_Solver *precond, INTRA *actintra, HYPREVARS *hyprevars)
{
int    i;
int    err=0;
int    numarg;
char **args;
#ifdef DEBUG 
dstrc_enter("hypre_create_precond_euclid");
#endif
/*----------------------------------------------------------------------*/
err+=HYPRE_EuclidCreate(actintra->MPI_INTRA_COMM,precond);
/*----------------------------------- allocate char space for arguments */
numarg=4;
args = (char**)calloc(numarg,sizeof(char*));
if (!args) dserror("Allocation for Euclid failed");
for (i=0; i<numarg; i++) 
{
   args[i] = (char*)calloc(15,sizeof(char));
   if (!args[i]) dserror("Allocation for Euclid failed");
}
/*------------------------------------------- set fill level for ilu(k) */
sprintf(args[0],"-level %d",hyprevars->ifill);
/*------------------------------------------------------------- set I/O */
if (hyprevars->io)
{
sprintf(args[1],"-eustats 1");
sprintf(args[2],"-eu_mem 1");
}
else
{
sprintf(args[1],"-eustats 0");
sprintf(args[2],"-eu_mem 0");
}
/*--------------------------------------- set block ILU or parallel ILU */
sprintf(args[3],"-bj 0");
/*-------------------------------------------- pass arguments to Euclid */
err+=HYPRE_EuclidSetParams(*precond,numarg,args);
/*-------------------------------------------- free space for arguments */
for (i=0; i<numarg; i++) free(args[i]);
free(args);
/*----------------------------------------------------------------------*/
if (err) dserror("error in calling solver package HYPRE");   
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of hypre_create_precond_euclid */
#endif /* end of ifdef HYPRE_PACKAGE */

#ifdef HYPRE_PACKAGE
/*----------------------------------------------------------------------*
 |  create preconditioner parasails and set parameters       m.gee 10/01|
 *----------------------------------------------------------------------*/
void hypre_create_precond_parasails(HYPRE_Solver *precond, INTRA *actintra, 
                                      HYPREVARS *hyprevars)
{
int    i;
int    err=0;
#ifdef DEBUG 
dstrc_enter("hypre_create_precond_parasails");
#endif
/*----------------------------------------------------------------------*/
/*----------------------------------------------- create preconditioner */
err+=HYPRE_ParaSailsCreate(actintra->MPI_INTRA_COMM,precond);
/*------------------------------------------------------ set parameters */
err+=HYPRE_ParaSailsSetSym(*precond,hyprevars->parasymm);
err+=HYPRE_ParaSailsSetFilter(*precond,hyprevars->parafilter);
err+=HYPRE_ParaSailsSetParams(*precond,hyprevars->parathresh,
                              hyprevars->paralevel);
/*--------------------------------------------------------------set I/O */
if (hyprevars->io) err+=HYPRE_ParaSailsSetLogging(*precond,1);
/*----------------------------------------------------------------------*/
if (err) dserror("error in calling solver package HYPRE");   
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of hypre_create_precond_parasails */
#endif /* end of ifdef HYPRE_PACKAGE */

#ifdef HYPRE_PACKAGE
/*----------------------------------------------------------------------*
 |  create preconditioner boomeramg and set parameters       m.gee 10/01|
 *----------------------------------------------------------------------*/
void hypre_create_precond_boomeramg(HYPRE_Solver *precond, INTRA *actintra, 
                                      HYPREVARS *hyprevars)
{
int    i;
int    err=0;
int    gridrelaxtype[4];
int    gridnsweep[4];
#ifdef DEBUG  
dstrc_enter("hypre_create_precond_boomeramg");
#endif 
/*----------------------------------------------------------------------*/
/*----------------------------------------------- create preconditioner */
err+=HYPRE_BoomerAMGCreate(precond);
/*------------------------------------------------------ set parameters */
if (hyprevars->io)
err+=HYPRE_BoomerAMGSetIOutDat(*precond,3);
else
err+=HYPRE_BoomerAMGSetIOutDat(*precond,0);
err+=HYPRE_BoomerAMGSetMaxIter(*precond,5);
err+=HYPRE_BoomerAMGSetMeasureType(*precond,1);
err+=HYPRE_BoomerAMGSetStrongThreshold(*precond,hyprevars->threshold);
gridrelaxtype[0]=6;
gridrelaxtype[1]=6;
gridrelaxtype[2]=6;
gridrelaxtype[3]=9;
err+=HYPRE_BoomerAMGSetGridRelaxType(*precond,gridrelaxtype);
/*--------------------------------------------------------- set v cycle */
err+=HYPRE_BoomerAMGSetCycleType(*precond,1);
/*--------------------------------------- set number of sweeps per grid */
for (i=0; i<4; i++) gridnsweep[i] = hyprevars->sweep[i];
err+=HYPRE_BoomerAMGSetNumGridSweeps(*precond,gridnsweep);
/*------------------------------------ set type of coarsening algorithm */
err+=HYPRE_BoomerAMGSetCoarsenType(*precond,0);
/*------------------------------------------------------------- set tol */
err+=HYPRE_BoomerAMGSetTol(*precond,0.0);
/*----------------------------------------------------------------------*/
if (err) dserror("error in calling solver package HYPRE");   
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of hypre_create_precond_boomeramg */
#endif /* end of ifdef HYPRE_PACKAGE */
