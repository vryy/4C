#include "../headers/standardtypes.h"
#include "../headers/solution.h"
/*----------------------------------------------------------------------*
 |  control solver lib AZTEC                             m.gee 9/01     |
 *----------------------------------------------------------------------*/
void solver_az_msr( 
                      struct _SOLVAR         *actsolv,
                      struct _INTRA          *actintra,
                      struct _AZ_ARRAY_MSR   *msr_array,
                      struct _DIST_VECTOR    *sol,
                      struct _DIST_VECTOR    *rhs,
                      int                     option
                     )
{
#ifdef AZTEC_PACKAGE
int        i;
int        numeq;
int        N_update;
int        nnz;
int        *external;
double    *val;
int       *bindx;
int       *update;
double    *b;
double    *x;
if (par.myrank==0) 
{
N_update  = 3;
nnz       = 13;
val       = (double*)calloc(nnz+1,sizeof(double));
bindx     = (int*)   calloc(nnz+1,sizeof(int));
update    = (int*)   calloc(N_update,sizeof(int));
b         = (double*)calloc(N_update,sizeof(double));
x         = (double*)calloc(N_update,sizeof(double));
update[0] = 0;
update[1] = 1;
update[2] = 3;
bindx[0]  = 4;
bindx[1]  = 7;
bindx[2]  = 9;
bindx[3]  = 14;
bindx[4]  = 1;
bindx[5]  = 3;
bindx[6]  = 4;
bindx[7]  = 0;
bindx[8]  = 3;
bindx[9]  = 0;
bindx[10] = 1;
bindx[11] = 2;
bindx[12] = 4;
bindx[13] = 5;
val[0]  = 6000.0;
val[1]  = 18000.0;
val[2]  = 20.000;
val[3]  = 0.0;
val[4]  = 2000.0;
val[5]  = 3000.0;
val[6]  = 4000.0;
val[7]  = 2000.0;
val[8]  = 5000.0;
val[9]  = 3000.0;
val[10] = 5000.0;
val[11] = 1000.0;
val[12] = 2500.0;
val[13] = 800.0;
b[0]    = 400.0;
b[1]    = 600.0;
b[2]    = 1000.0;
}
if (par.myrank==1) 
{
N_update  = 1;
nnz       = 4;
val       = (double*)calloc(nnz+1,sizeof(double));
bindx     = (int*)   calloc(nnz+1,sizeof(int));
update    = (int*)   calloc(N_update,sizeof(int));
b         = (double*)calloc(N_update,sizeof(double));
x         = (double*)calloc(N_update,sizeof(double));
update[0] = 4;
bindx[0]  = 2;
bindx[1]  = 5;
bindx[2]  = 0;
bindx[3]  = 3;
bindx[4]  = 2;
val[0]  = 12000.0;
val[1]  = 0.0;
val[2]  = 4000.0;
val[3]  = 1500.0;
val[4]  = 2500.0;
b[0]    = 1200.0;
}
if (par.myrank==2) 
{
N_update  = 2;
nnz       = 7;
val       = (double*)calloc(nnz+1,sizeof(double));
bindx     = (int*)   calloc(nnz+1,sizeof(int));
update    = (int*)   calloc(N_update,sizeof(int));
b         = (double*)calloc(N_update,sizeof(double));
x         = (double*)calloc(N_update+200,sizeof(double));
update[0] = 2;
update[1] = 5;
bindx[0]  = 3;
bindx[1]  = 6;
bindx[2]  = 8;
bindx[3]  = 3;
bindx[4]  = 4;
bindx[5]  = 5;
bindx[6]  = 2;
bindx[7]  = 3;
val[0]  = 15000.0;
val[1]  = 25000.0;
val[2]  = 0.0;
val[3]  = 1000.0;
val[4]  = 1500.0;
val[5]  = 1700.0;
val[6]  = 1700.0;
val[7]  = 800.0;
b[0]    = 800.0;
b[1]    = 1400.0;
}



/*AZ_set_proc_config(msr_array->proc_config,(MPI_AZComm)(actintra->MPI_INTRA_COMM));*/
   AZ_set_proc_config(msr_array->proc_config,(MPI_AZComm)MPI_COMM_WORLD);
   
   AZ_defaults(msr_array->options,msr_array->params);

   AZ_check_msr(
                  bindx,
                  N_update, 
                  0,
                  AZ_GLOBAL, 
                  msr_array->proc_config
               );

   msr_array->options[AZ_solver]          = AZ_cg;
   msr_array->options[AZ_precond]         = AZ_dom_decomp;/*AZ_none;*/
   msr_array->options[AZ_subdomain_solve] = AZ_ilu;
   msr_array->options[AZ_graph_fill]      = 1;

   msr_array->options[AZ_max_iter] = 1000;
   msr_array->options[AZ_overlap]  = 0;
   msr_array->options[AZ_poly_ord] = 5;
   msr_array->options[AZ_output]   = AZ_all;/*AZ_warnings;*//*AZ_last;*//*50;*/
   msr_array->options[AZ_conv]     = AZ_r0;
   msr_array->params[AZ_tol]       = 1.0E-6;
   msr_array->params[AZ_drop]      = 0.0;


   AZ_transform(
                msr_array->proc_config,
                &external,
                bindx,
                val,
                update,
                &(msr_array->update_index),
                &(msr_array->extern_index),
                &(msr_array->data_org),
                N_update,
                NULL,
                NULL,
                NULL,
                NULL,
                AZ_MSR_MATRIX
               );


   msr_array->Amat = AZ_matrix_create(
                                      msr_array->data_org[AZ_N_internal]+
                                      msr_array->data_org[AZ_N_border]
                                     );


   AZ_set_MSR(
              msr_array->Amat,
              bindx,
              val,
              msr_array->data_org,
              0,
              NULL,
              AZ_LOCAL
             );


   msr_array->N_external = msr_array->data_org[AZ_N_external];


   AZ_reorder_vec(
                  b,
                  msr_array->data_org,
                  msr_array->update_index,
                  NULL
                 );

   AZ_reorder_vec(
                  x,
                  msr_array->data_org,
                  msr_array->update_index,
                  NULL
                 );


   AZ_iterate(
              x,
              b,
              msr_array->options,
              msr_array->params,
              msr_array->status,
              msr_array->proc_config,
              msr_array->Amat,
              NULL,
              NULL
             );



   AZ_invorder_vec(
                   x,
                   msr_array->data_org,
                   msr_array->update_index,
                   NULL,
                   b
                  );
 
 







#endif /* end of ifdef AZTEC_PACKAGE */
return;
} /* end of solver_az_msr */




