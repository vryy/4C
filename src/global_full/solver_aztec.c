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
 |  control solver lib AZTEC                             m.gee 9/01     |
 *----------------------------------------------------------------------*/
void solver_az_msr( 
                      struct _SOLVAR         *actsolv,
                      struct _INTRA          *actintra,
                      struct _AZ_ARRAY_MSR   *msr_array,
                      struct _DIST_VECTOR    *sol,
                      struct _DIST_VECTOR    *rhs,
                      INT                     option
                     )
{
#ifdef AZTEC_PACKAGE
INT         i;
INT         dim;
/* INT         reuse; */
INT         azname;
AZVAR      *azvar;

DOUBLE     *dfrom,*dto;

DOUBLE     *tmpsol;
ARRAY       tmpsol_a;

DOUBLE     *tmprhs;
ARRAY       tmprhs_a;

/* DOUBLE      l2norm; */
#ifdef DEBUG 
dstrc_enter("solver_az_msr");
#endif
/*----------------------------------------------------------------------*/
azvar = actsolv->azvar;
/*----------------------------------------------------------------------*/
switch(option)
{
/*----------------------------------------------------------------------*/
/*                                                           init phase */
/*----------------------------------------------------------------------*/
case 1:
   /*----------- make processor configuration, dependend on parallelism */
   #ifdef PARALLEL
   AZ_set_proc_config(msr_array->proc_config, (MPI_AZComm)(actintra->MPI_INTRA_COMM));
   #else
   AZ_set_proc_config(msr_array->proc_config, AZ_NOT_MPI);
   #endif
   /*------------------------------------- set default value to options */
   AZ_defaults(msr_array->options,msr_array->params);
   /*-------------------------------------- perform check of msr matrix */
   #ifdef DEBUG 
   AZ_check_msr(
                &(msr_array->bindx.a.iv[0]),
                  msr_array->numeq, 
                  msr_array->N_external,
                  AZ_GLOBAL, 
                  msr_array->proc_config
               );
   #endif
   /*----------------------- set options and params from the input file */
   switch(azvar->azsolvertyp)/*--------------------------- set solver */
   {
   case azsolv_CG:
      msr_array->options[AZ_solver] = AZ_cg;
   break;
   case azsolv_GMRES:
      msr_array->options[AZ_solver] = AZ_gmres;
      msr_array->options[AZ_kspace] = azvar->azsub;
   break;
   case azsolv_CGS:
      msr_array->options[AZ_solver] = AZ_cgs;
   break;
   case azsolv_BiCGSTAB:
      msr_array->options[AZ_solver] = AZ_bicgstab;
   break;
   case azsolv_LU:
      msr_array->options[AZ_solver] = AZ_lu;
   break;
   case azsolv_TFQMR:
      msr_array->options[AZ_solver] = AZ_tfqmr;
   break;
   default:
      dserror("No correct solver for Aztec");
   }
   switch(azvar->azprectyp)/*--------------------- set preconditioner */
   {
   case azprec_none:
      msr_array->options[AZ_precond]         = AZ_none;
   break;
   case azprec_ILUT:
      msr_array->options[AZ_precond]         = AZ_dom_decomp;
      msr_array->options[AZ_subdomain_solve] = AZ_ilut;
      msr_array->params[AZ_ilut_fill]        = azvar->azfill;
   break;
   case azprec_ILU:
      msr_array->options[AZ_precond]         = AZ_dom_decomp;
      msr_array->options[AZ_subdomain_solve] = AZ_ilu;
      msr_array->options[AZ_graph_fill]      = azvar->azgfill;
   break;
   case azprec_Jacobi:
      msr_array->options[AZ_precond]         = AZ_Jacobi;
      msr_array->options[AZ_poly_ord]        = azvar->azpoly;
   break;
   case azprec_Neumann:
      msr_array->options[AZ_precond]         = AZ_Neumann;
      msr_array->options[AZ_poly_ord]        = azvar->azpoly;
   break;
   case azprec_Least_Squares:
      msr_array->options[AZ_precond]         = AZ_ls;
      msr_array->options[AZ_poly_ord]        = azvar->azpoly;
   break;
   case azprec_SymmGaussSeidel:
      msr_array->options[AZ_precond]         = AZ_sym_GS;
      msr_array->options[AZ_poly_ord]        = azvar->azpoly;
   break;
   case azprec_LU:
      msr_array->options[AZ_precond]         = AZ_dom_decomp;
      msr_array->options[AZ_subdomain_solve] = AZ_lu;
      msr_array->params[AZ_drop]             = azvar->azdrop;
   break;
   case azprec_RILU:
      msr_array->options[AZ_precond]         = AZ_dom_decomp;
      msr_array->options[AZ_subdomain_solve] = AZ_rilu;
      msr_array->options[AZ_graph_fill]      = azvar->azgfill;
   break;
   case azprec_BILU:
      dserror("Block Preconditioning Bilu cannot be used in MSR format"); 
      msr_array->options[AZ_precond]         = AZ_dom_decomp;
      msr_array->options[AZ_subdomain_solve] = AZ_bilu;
      msr_array->options[AZ_graph_fill]      = azvar->azgfill;
   break;
   case azprec_ICC:
      msr_array->options[AZ_precond]         = AZ_dom_decomp;
      msr_array->options[AZ_subdomain_solve] = AZ_icc;
      msr_array->options[AZ_graph_fill]      = azvar->azgfill;
   break;
   default:
      dserror("No correct preconditioner for Aztec");
   }
   /*---------------------------------------------- set rest of options */
   msr_array->options[AZ_max_iter] = azvar->aziter;
   msr_array->options[AZ_overlap]  = 0;
   msr_array->options[AZ_poly_ord] = azvar->azpoly;
   msr_array->options[AZ_output]   = AZ_none;/*AZ_all;AZ_warnings;AZ_last;10; */
   msr_array->options[AZ_conv]     = AZ_noscaled;
   msr_array->params[AZ_tol]       = azvar->aztol;
   msr_array->params[AZ_drop]      = azvar->azdrop;
   msr_array->options[AZ_scaling]  = AZ_sym_diag; /*AZ_none  */
   /*--------- make backup copy of bindx, as it is permuted in solution */
   /*am_alloc_copy(&(msr_array->bindx),&(msr_array->bindx_backup));*/
   /*-------------------------------------- allocate backup copy of val */
   /*amdef("val_back",&(msr_array->val_backup),msr_array->val.fdim,1,"DV");*/
   /*----------------------------- set NULL-pointers for Amat and Aprec */
   msr_array->Amat  = NULL;
   msr_array->Aprec = NULL;
   msr_array->ncall=0;
   /* set flag, that this matrix has been initialized and is ready for solve */   
   msr_array->is_init=1;
break;
/*----------------------------------------------------------------------*/
/*                                                    end of init phase */
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*                                                    calculation phase */
/*----------------------------------------------------------------------*/
case 0:
#ifdef PERF
  perf_begin(31);
#endif
/*--------------------------------------------- check the reuse feature */
/* NOTE: This is not multifield yet
 * So right now we can use aztec for one field only in a multi field
 * problem (?)
 * There are problems with more than one discretisation (?)
 *
 * Lets assume we need different reuse information for every
 * matrix. Then we'd have to calculate a different number for every
 * matrix in every field. We know the number of matrices in the
 * current field is actsolv->nsysarray, but we have no idea about the
 * other fields. That's why MAXNUMMATRICES is needed.
 *
 * Of course there might be matrices that are never used with this
 * routine. But we cannot know for sure...
 */
    dsassert(MAXNUMMATRICES > actsolv->nsysarray, "Maximum number of matrices exceeded");
    switch (actsolv->fieldtyp)
    {
    case structure: azname=0*MAXNUMMATRICES + actsolv->nsysarray; break;
    case fluid:     azname=1*MAXNUMMATRICES + actsolv->nsysarray; break;
    case ale:       azname=2*MAXNUMMATRICES + actsolv->nsysarray; break;
    case levelset:  azname=3*MAXNUMMATRICES + actsolv->nsysarray; break;
    default:
      dserror("Unknown type of field");
      break;
    }
/*--------------------------------------------- check the reuse feature */
/* if (msr_array->is_factored==0) amcopy(&(msr_array->val),&(msr_array->val_backup)); */
/* else                           amcopy(&(msr_array->val_backup),&(msr_array->val)); */
/*----------------------- transform matrix to processor local numbering */
    /*
     * Transformation and factorization are different. The
     * transformation might be done by the matrix vector product. The
     * factorization is only ever done by the solver.
     */
    if (msr_array->is_transformed==0) {
      msr_array->is_transformed = 1;
      
      if (msr_array->Amat != NULL) {
        AZ_matrix_destroy(&(msr_array->Amat)); msr_array->Amat        =NULL;
      }
      free(msr_array->external);             msr_array->external      =NULL;
      free(msr_array->update_index);         msr_array->update_index  =NULL;
      free(msr_array->extern_index);         msr_array->extern_index  =NULL;
      free(msr_array->data_org);             msr_array->data_org      =NULL;
      
      AZ_transform(msr_array->proc_config,
                   &(msr_array->external),
                   msr_array->bindx.a.iv,
                   msr_array->val.a.dv,
                   msr_array->update.a.iv,
                   &(msr_array->update_index),
                   &(msr_array->extern_index),
                   &(msr_array->data_org),
                   msr_array->numeq,
                   NULL,
                   NULL,
                   NULL,
                   NULL,
                   AZ_MSR_MATRIX);

      /* create Aztec structure AZ_MATRIX */
      msr_array->Amat = AZ_matrix_create(msr_array->data_org[AZ_N_internal]+
                                         msr_array->data_org[AZ_N_border]);

      /* attach dmsr-matrix to this structure */
      AZ_set_MSR(msr_array->Amat, 
                 msr_array->bindx.a.iv, 
                 msr_array->val.a.dv, 
                 msr_array->data_org, 
                 0, 
                 NULL, 
                 AZ_LOCAL);
    }
/*--------------------- save number of external components on this proc */            
msr_array->N_external = msr_array->data_org[AZ_N_external];
/*-------------------------------------------------- reorder rhs-vector */            
tmprhs = am_alloc_copy(&(rhs->vec),&(tmprhs_a));
AZ_reorder_vec(
               tmprhs,
               msr_array->data_org,
               msr_array->update_index,
               NULL
              );
/*--------------------------- reorder initial guess and solution-vector */            
AZ_reorder_vec(
               sol->vec.a.dv,
               msr_array->data_org,
               msr_array->update_index,
               NULL
              );
/*----- allocate temporary solution vector large enough for N_externals */
tmpsol = amdef("tmpsol",&tmpsol_a,(msr_array->numeq+msr_array->N_external),1,"DV");
/*--------------------- copy initial guess to temporary solution vector */
/*             (this looks a bit strange, but it's supposed to be fast) */
dfrom = sol->vec.a.dv;
dto   = tmpsol_a.a.dv;
dim   = sol->vec.fdim;
for (i=0; i<dim; i++) *(dto++) = *(dfrom++);
dto   = &(tmpsol_a.a.dv[dim]);
dim   = tmpsol_a.fdim - sol->vec.fdim;
for (i=0; i<dim; i++) *(dto++) = 0.0;
/*--------------------------------------------- check the reuse feature */
msr_array->data_org[AZ_name]=azname;
/*---------------------------------------------------------- first call */
    if (msr_array->ncall==0) {
      msr_array->options[AZ_pre_calc]  = AZ_calc;
      msr_array->options[AZ_keep_info] = 1;
    }
/*------------------------------------------------------ not first call */
    else {
      if (msr_array->is_factored==0) {
        msr_array->options[AZ_pre_calc] = AZ_recalc;
        msr_array->options[AZ_keep_info] = 1;
      }
      else {
        msr_array->options[AZ_pre_calc] = AZ_reuse;
        msr_array->options[AZ_keep_info] = 1;
      }
    }
/*--------------------------------------------------------- call solver */
#ifdef PERF
  perf_end(31);
#endif
#ifdef PERF
  perf_begin(32);
#endif
/* Let's try several times. Normally the first try should succeed. But
 * there are issues with BiCGSTAB and fluid fields. In case of a
 * breakdown we simply start again and use the current solution as
 * initial guess. See the aztec manual.
 */
for (i=0; i<10; ++i) {
AZ_iterate(
           tmpsol,
           tmprhs,
           msr_array->options,
           msr_array->params,
           msr_array->status,
           msr_array->proc_config,
           msr_array->Amat,
           NULL,
           NULL
          );
if ( (DOUBLE)(msr_array->status[AZ_why]) == AZ_breakdown ) {
  if (actintra->intra_rank==0) {
    printf("Numerical breakdown occured in solver Aztec -> restart aztec with new initial guess\n");
  }

  /* Is this needed? There ought to be some scaling involved... */
  msr_array->options[AZ_pre_calc] = AZ_reuse;
  msr_array->options[AZ_keep_info] = 1;
}
else {
  break;
}
}
#ifdef PERF
  perf_end(32);
#endif
#ifdef PERF
  perf_begin(33);
#endif
/*------------------------------------------------ delete temporary rhs */          
amdel(&tmprhs_a);
/*-------------------------------------------- recover unpermuted bindx */
/* amcopy(&(msr_array->bindx_backup),&(msr_array->bindx)); */
/*------------------------------------------------ invorder solv vector */
AZ_invorder_vec(
                tmpsol,
                msr_array->data_org,
                msr_array->update_index,
                NULL,
                sol->vec.a.dv
               );
/*------------------------------------ delete temporary solution vector */
amdel(&tmpsol_a);
/*---------------------------------------- destroy the Aztec structures */
/* AZ_matrix_destroy(&(msr_array->Amat)); msr_array->Amat          =NULL; */
/* free(msr_array->external);             msr_array->external      =NULL; */
/* free(msr_array->update_index);         msr_array->update_index  =NULL; */
/* free(msr_array->extern_index);         msr_array->extern_index  =NULL; */
/* free(msr_array->data_org);             msr_array->data_org      =NULL; */
/*----------------------------------------- check for success of solver */
if ( (DOUBLE)(msr_array->status[AZ_why]) != AZ_normal )
{
   if (actintra->intra_rank==0)
   {
      /*----------------------------------------------------- breakdown */
      if ( (DOUBLE)(msr_array->status[AZ_why]) == AZ_breakdown )
      dserror("Numerical breakdown occured in solver Aztec -> Abort");
      /*------------------------------------numerical loss of precision */
      if ( (DOUBLE)(msr_array->status[AZ_why]) == AZ_loss )
      {
          printf("RANK 0: AZTEC: Numerical loss of precision occured! continue...\n");
          fprintf(allfiles.out_err,"RANK 0: AZTEC: Numerical loss of precision occured, continue...\n");
      }
      /*------------------------------------------------------ ill cond */
      if ( (DOUBLE)(msr_array->status[AZ_why]) == AZ_ill_cond )
      {
         printf("RANK 0: AZTEC: Preconditioning ill-conditioned or singular,\n");
         printf("               solution is least square ! continue...\n");
         fprintf(allfiles.out_err,"RANK 0: AZTEC: Preconditioning ill-conditioned or singular,\n");
         fprintf(allfiles.out_err,"               solution is least square ! continue...\n");
      }
      /*-------------------------- maximum number of iterations reached */
      if ( (DOUBLE)(msr_array->status[AZ_why]) == AZ_maxits )
      {
         printf("RANK 0: AZTEC: Maximum number of iterations %d reached \n",msr_array->options[AZ_max_iter]);
         fprintf(allfiles.out_err,"RANK 0: AZTEC: Maximum number of iterations %d reached \n",msr_array->options[AZ_max_iter]);
         fflush(allfiles.out_err);
         fflush(stdout);
      }
   }
}
/*------------------------------------ print solver iterations and time */             
if (actintra->intra_rank==0)
{
   if (actsolv->fieldtyp==structure) fprintf(allfiles.out_err,"Structure:\n");
   if (actsolv->fieldtyp==fluid)     fprintf(allfiles.out_err,"Fluid:\n");
   if (actsolv->fieldtyp==ale)       fprintf(allfiles.out_err,"Ale:\n");
   fprintf(allfiles.out_err,"AZTEC: %d unknowns %d iterations %f solving time\n",
   sol->numeq_total,
   (INT)(msr_array->status[AZ_its]),
   msr_array->status[AZ_solve_time]);
}
/*----------------------------------------------------------- set flags */
msr_array->ncall++;
msr_array->is_factored=1;
#ifdef PERF
  perf_end(33);
#endif
break;
/*----------------------------------------------------------------------*/
/*                                             end of calculation phase */
/*----------------------------------------------------------------------*/
default:
   dserror("Unknown option for solver call to Aztec");
break;   
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
#endif /* end of ifdef AZTEC_PACKAGE */
return;
} /* end of solver_az_msr */




