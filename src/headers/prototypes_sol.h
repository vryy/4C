/*----------------------------------------------------------------------*
 | cal_nlnstatic_control.c                               m.gee 11/01    |
 *----------------------------------------------------------------------*/
void stanln();
void conpre(
            FIELD         *actfield,
            SOLVAR        *actsolv,
            PARTITION     *actpart,
            INTRA         *actintra,
            int            kstep,
            int            actsysarray,
            DIST_VECTOR   *rsd,
            DIST_VECTOR   *dispi,
            int            cdof,
            STANLN        *nln_data,
            NR_CONTROLTYP  controltyp
          );
void conequ(
            FIELD         *actfield,
            SOLVAR        *actsolv,
            PARTITION     *actpart,
            INTRA         *actintra,
            int            kstep,
            int           *itnum,
            int            actsysarray,
            DIST_VECTOR   *rsd,
            DIST_VECTOR   *dispi,
            DIST_VECTOR   *re,
            int            cdof,
            STANLN        *nln_data,
            NR_CONTROLTYP  controltyp
          );
void conequ_printhead(int kstep, NR_CONTROLTYP  controltyp, int cdof);
void conequ_printiter(int itnum, double disval, double rlnew, double dinorm,
                     double renorm, double energy, double dnorm, double rrnorm);
/*----------------------------------------------------------------------*
 | global_calelm.c                                       m.gee 11/01    |
 *----------------------------------------------------------------------*/
void calelm(FIELD      *actfield, 
               SOLVAR     *actsolv, 
               PARTITION  *actpart, 
               INTRA      *actintra,
               int         sysarray1,
               int         sysarray2,
               double     *dvec,
               int         global_numeq,
               int         kstep,
               int         calc_option);
void calinit(FIELD      *actfield, 
                PARTITION  *actpart);
/*----------------------------------------------------------------------*
 | global_calrhs.c                                       m.gee 11/01    |
 *----------------------------------------------------------------------*/
void calrhs(FIELD        *actfield, 
               SOLVAR       *actsolv, 
               PARTITION    *actpart, 
               INTRA        *actintra,
               int           actsysarray,
               DIST_VECTOR  *rhs1,
               DIST_VECTOR  *rhs2,
               int           kstep,
               int           calc_option);
void calrhs_nodal_neumann(
                             PARTITION    *actpart,
                             INTRA        *actintra,
                             SPARSE_TYP   *sysarraytyp,
                             SPARSE_ARRAY *sysarray,
                             DIST_VECTOR  *rhs
                            );
void calrhs_ele_neumann(FIELD        *actfield, 
                           SOLVAR       *actsolv,
                           PARTITION    *actpart,
                           INTRA        *actintra,
                           SPARSE_TYP   *sysarraytyp,
                           SPARSE_ARRAY *sysarray,
                           int           actsysarray,
                           DIST_VECTOR  *rhs,
                           int           kstep,
                           int           calc_option);
/*----------------------------------------------------------------------*
 | global_calrhs_nodal.c                                 m.gee 11/01    |
 *----------------------------------------------------------------------*/
void assemble_nn(
                    PARTITION    *actpart,
                    DIST_VECTOR  *rhs,
                    double       *drhs
                   );
/*----------------------------------------------------------------------*
 | global_mask_dense.c                                   m.gee 11/01    |
 *----------------------------------------------------------------------*/
void mask_dense(FIELD         *actfield, 
                  PARTITION     *actpart, 
                  SOLVAR        *actsolv,
                  INTRA         *actintra, 
                  DENSE         *dense);
void  dense_numeq(FIELD         *actfield, 
                    PARTITION    *actpart, 
                    SOLVAR       *actsolv,
                    INTRA        *actintra,
                    int          *numeq);
void  dense_update(FIELD         *actfield, 
                     PARTITION     *actpart, 
                     SOLVAR        *actsolv,
                     INTRA         *actintra,
                     DENSE         *dense);
/*----------------------------------------------------------------------*
 | global_mask_parcsr.c                                   m.gee 11/01    |
 *----------------------------------------------------------------------*/
void  mask_parcsr(FIELD         *actfield, 
                    PARTITION     *actpart, 
                    SOLVAR        *actsolv,
                    INTRA         *actintra, 
                    H_PARCSR  *parcsr);
void parcsr_update(FIELD         *actfield, 
                     PARTITION     *actpart, 
                     SOLVAR        *actsolv,
                     INTRA         *actintra,
                     H_PARCSR      *parcsr);
void parcsr_update_perm(INTRA         *actintra,
                          H_PARCSR      *parcsr);
void  parcsr_nnz_topology(FIELD         *actfield, 
                            PARTITION    *actpart, 
                            SOLVAR       *actsolv,
                            INTRA        *actintra,
                            H_PARCSR     *parcsr,
                            int         **dof_connect);
void parcsr_make_bindx(FIELD         *actfield, 
                          PARTITION     *actpart, 
                          SOLVAR        *actsolv,
                          H_PARCSR      *parcsr,
                          int          **dof_connect);
/*----------------------------------------------------------------------*
 | global_mask_msr.c                                   m.gee 11/01    |
 *----------------------------------------------------------------------*/
void mask_msr(FIELD         *actfield, 
                 PARTITION     *actpart, 
                 SOLVAR        *actsolv,
                 INTRA         *actintra, 
                 AZ_ARRAY_MSR  *msr);
void msr_numeq(FIELD         *actfield, 
                  PARTITION    *actpart, 
                  SOLVAR       *actsolv,
                  INTRA        *actintra,
                  int          *numeq);
void msr_update(FIELD         *actfield, 
                   PARTITION     *actpart, 
                   SOLVAR        *actsolv,
                   INTRA         *actintra,
                   AZ_ARRAY_MSR  *msr);
void msr_nnz_topology(FIELD         *actfield, 
                         PARTITION    *actpart, 
                         SOLVAR       *actsolv,
                         INTRA        *actintra,
                         AZ_ARRAY_MSR *msr,
                         int         **dof_connect);
void dof_in_coupledofs(int dof, PARTITION *actpart, int *iscoupled);
void dof_find_centernode(int dof, PARTITION *actpart, NODE **centernode);
void msr_make_bindx(FIELD         *actfield, 
                       PARTITION     *actpart, 
                       SOLVAR        *actsolv,
                       AZ_ARRAY_MSR  *msr,
                       int          **dof_connect);
/*----------------------------------------------------------------------*
 |  global_mask_ucchb.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void  mask_ucchb(FIELD     *actfield, 
                   PARTITION *actpart, 
                   SOLVAR    *actsolv,
                   INTRA     *actintra, 
                   UCCHB     *ucchb);
void ucchb_numeq(FIELD      *actfield, 
                   PARTITION    *actpart, 
                   SOLVAR       *actsolv,
                   INTRA        *actintra,
                   int          *numeq);
void  ucchb_update(FIELD     *actfield, 
                     PARTITION *actpart, 
                     SOLVAR    *actsolv,
                     INTRA     *actintra,
                     UCCHB     *ucchb);
void  ucchb_nnz_topology(FIELD      *actfield, 
                           PARTITION  *actpart, 
                           SOLVAR     *actsolv,
                           INTRA      *actintra,
                           UCCHB      *ucchb,
                           int       **dof_connect);
void  ucchb_make_a(FIELD         *actfield, 
                     PARTITION     *actpart, 
                     SOLVAR        *actsolv,
                     INTRA         *actintra,
                     UCCHB         *ucchb,
                     int          **dof_connect);
/*----------------------------------------------------------------------*
 |  global_mask_rcptr.c                                  m.gee 01/02    |
 *----------------------------------------------------------------------*/
void mask_rc_ptr(FIELD         *actfield, 
                 PARTITION     *actpart, 
                 SOLVAR        *actsolv,
                 INTRA         *actintra, 
                 RC_PTR        *rc_ptr);
void  rc_ptr_update(FIELD         *actfield, 
                   PARTITION     *actpart, 
                   SOLVAR        *actsolv,
                   INTRA         *actintra,
                   RC_PTR        *rc_ptr);
void  rc_ptr_nnz_topology(FIELD         *actfield, 
                         PARTITION    *actpart, 
                         SOLVAR       *actsolv,
                         INTRA        *actintra,
                         RC_PTR       *rc_ptr,
                         int         **dof_connect);
void  rc_ptr_make_bindx(FIELD         *actfield, 
                       PARTITION     *actpart, 
                       SOLVAR        *actsolv,
                       RC_PTR        *rc_ptr,
                       int          **dof_connect);
void  rc_ptr_make_sparsity(RC_PTR        *rc_ptr);
/*----------------------------------------------------------------------*
 |  solver_add_data.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void assemble(
                 int                sysarray1,
                 struct _ARRAY     *elearray1,
                 int                sysarray2,
                 struct _ARRAY     *elearray2,
                 struct _PARTITION *actpart,
                 struct _SOLVAR    *actsolv,
                 struct _INTRA     *actintra,
                 struct _ELEMENT   *actele,
                 int                option
                );
void init_assembly(
                       struct _PARTITION      *actpart,
                       struct _SOLVAR         *actsolv,
                       struct _INTRA          *actintra,
                       struct _FIELD          *actfield,
                       int                     actsysarray
                     );
void assemble_vec(INTRA        *actintra,
                    SPARSE_TYP   *sysarraytyp,
                    SPARSE_ARRAY *sysarray,
                    DIST_VECTOR  *rhs,
                    double       *drhs,
                    double        factor);
/*----------------------------------------------------------------------*
 |  solver_add_dense.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void  add_dense(struct _PARTITION     *actpart,
                  struct _SOLVAR        *actsolv,
                  struct _INTRA         *actintra,
                  struct _ELEMENT       *actele,
                  struct _DENSE         *dense);
void redundant_dense(
                        PARTITION     *actpart,
                        SOLVAR        *actsolv,
                        INTRA         *actintra,
                        DENSE         *dense
                        );
/*----------------------------------------------------------------------*
 |  solver_add_msr.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void  add_msr(struct _PARTITION     *actpart,
                struct _SOLVAR        *actsolv,
                struct _INTRA         *actintra,
                struct _ELEMENT       *actele,
                struct _AZ_ARRAY_MSR  *msr);
void add_msr_checkcouple(int ii,int **cdofs,int ncdofs,int *iscouple,
                           int *isowner, int nprocs);
void add_msr_sendbuff(int ii,int jj,int i,int j,int ii_owner,int **isend,
                    double **dsend,double **estif, int numsend);
void exchange_coup_msr(
                         PARTITION     *actpart,
                         SOLVAR        *actsolv,
                         INTRA         *actintra,
                         AZ_ARRAY_MSR  *msr
                        );
/*----------------------------------------------------------------------*
 | solver_add_parcsr.c                                   m.gee 11/01    |
 *----------------------------------------------------------------------*/
void  add_parcsr(struct _PARTITION     *actpart,
                   struct _SOLVAR        *actsolv,
                   struct _INTRA         *actintra,
                   struct _ELEMENT       *actele,
                   struct _H_PARCSR      *parcsr);
void add_parcsr_sendbuff(int ii,int jj,int i,int j,int ii_owner,int **isend,
                    double **dsend,double **estif, int numsend);
void add_parcsr_checkcouple(int ii,int **cdofs,int ncdofs,int *iscouple,int *isowner, int nprocs);
void exchange_coup_parcsr(
                             PARTITION     *actpart,
                             SOLVAR        *actsolv,
                             INTRA         *actintra,
                             H_PARCSR      *parcsr
                            );
/*----------------------------------------------------------------------*
 |  solver_add_ucchb.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void  add_ucchb(struct _PARTITION     *actpart,
                  struct _SOLVAR        *actsolv,
                  struct _INTRA         *actintra,
                  struct _ELEMENT       *actele,
                  struct _UCCHB         *ucchb);
void redundant_ucchb(
                        PARTITION     *actpart,
                        SOLVAR        *actsolv,
                        INTRA         *actintra,
                        UCCHB         *ucchb
                        );
/*----------------------------------------------------------------------*
 |  solver_aztec.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void solver_az_msr( 
                      struct _SOLVAR         *actsolv,
                      struct _INTRA          *actintra,
                      struct _AZ_ARRAY_MSR   *msr_array,
                      struct _DIST_VECTOR    *sol,
                      struct _DIST_VECTOR    *rhs,
                      int                     option
                     );
/*----------------------------------------------------------------------*
 |  solver_control.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void solver_control(
                       struct _SOLVAR         *actsolv,
                       struct _INTRA          *actintra,
                       enum   _SPARSE_TYP     *sysarray_typ,
                       union  _SPARSE_ARRAY   *sysarray,
                       struct _DIST_VECTOR    *sol,
                       struct _DIST_VECTOR    *rhs,
                       int                     option
                      );
/*----------------------------------------------------------------------*
 |  solver_hypre.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
#ifdef HYPRE_PACKAGE
void  solver_hypre_parcsr( 
                            struct _SOLVAR         *actsolv,
                            struct _INTRA          *actintra,
                            struct _H_PARCSR       *parcsr,
                            struct _DIST_VECTOR    *sol,
                            struct _DIST_VECTOR    *rhs,
                            int                     option
                           );
void hypre_vector_create(int              myrank,
                           H_PARCSR        *parcsr,
                           INTRA           *actintra,
                           HYPRE_IJVector  *ijvector);
void hypre_vector_assemble(HYPRE_IJVector  *ijvector,
                             HYPRE_ParVector *parcsr_vector);
void hypre_matrix_create(int              myrank,
                           H_PARCSR        *parcsr,
                           INTRA           *actintra,
                           HYPRE_IJMatrix  *ijmatrix);
void hypre_matrix_assemble(HYPRE_IJMatrix     *ij_matrix,
                             HYPRE_ParCSRMatrix *parcsr_matrix);
void hypre_set_params_boomeramg_solver(HYPRE_Solver *solver,HYPREVARS *hyprevars);
void hypre_set_params_gmres_solver(HYPRE_Solver *solver,HYPREVARS *hyprevars);
void hypre_set_params_bicgstab_solver(HYPRE_Solver *solver,HYPREVARS *hyprevars);
void hypre_set_params_pcg_solver(HYPRE_Solver *solver,HYPREVARS *hyprevars);
void hypre_create_precond_euclid(HYPRE_Solver *precond, INTRA *actintra, HYPREVARS *hyprevars);
void hypre_create_precond_parasails(HYPRE_Solver *precond, INTRA *actintra, 
                                      HYPREVARS *hyprevars);
void hypre_create_precond_boomeramg(HYPRE_Solver *precond, INTRA *actintra, 
                                      HYPREVARS *hyprevars);
#endif
/*----------------------------------------------------------------------*
 |  solver_lapack.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void solver_lapack_dense( 
                              struct _SOLVAR         *actsolv,
                              struct _INTRA          *actintra,
                              struct _DENSE          *dense,
                              struct _DIST_VECTOR    *sol,
                              struct _DIST_VECTOR    *rhs,
                              int                     option
                             );
/*----------------------------------------------------------------------*
 |  solver_service.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void solserv_create_vec(
                           DIST_VECTOR         **vector,
                           int                   numvectors,
                           int                   numeq_total,
                           int                   numeq,
                           char                  typstr[]);
void solserv_zero_vec(DIST_VECTOR *disvector);
void solserv_zero_mat(INTRA *actintra, SPARSE_ARRAY *mat,SPARSE_TYP *mattyp);
void solserv_add_vec(DIST_VECTOR *vec_from, 
                        DIST_VECTOR *vec_to);
void solserv_copy_vec(DIST_VECTOR *vec_from, 
                        DIST_VECTOR *vec_to);
void solserv_vecnorm_euclid(INTRA       *actintra,
                            DIST_VECTOR *dist_vec, 
                            double      *result);
void solserv_getele_vec(INTRA       *actintra,
                        SPARSE_TYP   *sysarray_typ,
                        SPARSE_ARRAY *sysarray,
                        DIST_VECTOR *dist_vec,
                        int          indiz, 
                        double      *result);
void solserv_dot_vec(INTRA       *actintra,
                     DIST_VECTOR *dist_vec1,
                     DIST_VECTOR *dist_vec2,
                     double      *dot);
void solserv_scalarprod_vec(DIST_VECTOR *dist_vec,
                            double       scalar);
void solserv_reddistvec(DIST_VECTOR  *distvec,
                          SPARSE_ARRAY *sysarray,
                          SPARSE_TYP   *sysarray_typ,
                          double       *fullvec,
                          int           dim,
                          INTRA        *actintra);
void solserv_result_total(
                          FIELD          *actfield,
                          INTRA          *actintra,
                          DIST_VECTOR    *sol,
                          int             place,
                          SPARSE_ARRAY   *sysarray,
                          SPARSE_TYP     *sysarray_typ
                         );
void solserv_result_incre(
                          FIELD          *actfield,
                          INTRA          *actintra,
                          DIST_VECTOR    *sol,
                          int             place,
                          SPARSE_ARRAY   *sysarray,
                          SPARSE_TYP     *sysarray_typ
                         );
void solserv_result_resid(
                          FIELD          *actfield,
                          INTRA          *actintra,
                          DIST_VECTOR    *sol,
                          int             place,
                          SPARSE_ARRAY   *sysarray,
                          SPARSE_TYP     *sysarray_typ
                         );
/*----------------------------------------------------------------------*
 |  solver_superlu.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void solver_psuperlu_ucchb( 
                      struct _SOLVAR         *actsolv,
                      struct _INTRA          *actintra,
                      struct _UCCHB          *ucchb,
                      struct _DIST_VECTOR    *sol,
                      struct _DIST_VECTOR    *rhs,
                      int                     option
                     );
/*----------------------------------------------------------------------*
 |  input_sol.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void inpctrsol(SOLVAR *solv);
