/*----------------------------------------------------------------------*
 | cal_nlnstatic_control.c                               m.gee 11/01    |
 *----------------------------------------------------------------------*/
void stanln();
void conpre(
            FIELD         *actfield,     /* the actual physical field */
            SOLVAR        *actsolv,      /* the field-corresponding solver */
            PARTITION     *actpart,      /* the partition of the proc */
            INTRA         *actintra,     /* the intra-communicator of this field */
            CALC_ACTION   *action,       /* calculation flag */
            int            kstep,        /* the load or time step we are in */
            int            actsysarray,  /* number of the system matrix in actsolv->sysarray[actsysarray] to be used */
            DIST_VECTOR   *rsd,          /* dist. vector of incremental residual forces used for iteration in conequ */
            DIST_VECTOR   *dispi,        /* dist. vector of incremental displacements */
            int            cdof,         /* number of the dof to be controlled */
            STANLN        *nln_data,     /* data of the Newton-Raphson method */
            NR_CONTROLTYP  controltyp    /* type of control algorithm */
          );
void conequ(
            FIELD         *actfield,      /* the actual physical field */
            SOLVAR        *actsolv,       /* the field-corresponding solver */
            PARTITION     *actpart,       /* the partition of the proc */
            INTRA         *actintra,      /* the intra-communicator of this field */
            CALC_ACTION   *action,        /* calculation flag */
            int            kstep,         /* the load or time step we are in */
            int           *itnum,         /* number of corrector steps taken by this routine */
            int            actsysarray,   /* number of the system matrix in actsolv->sysarray[actsysarray] to be used */
            DIST_VECTOR   *rsd,           /* dist. vector of incremental residual forces used for iteration in conequ */
            DIST_VECTOR   *dispi,         /* dist. vector of incremental displacements */
            DIST_VECTOR   *re,            /* re[0..2] 3 vectors for residual displacements */
            int            cdof,          /* number of dof to be controlled */
            STANLN        *nln_data,      /* data of the Newton-Raphson method */
            NR_CONTROLTYP  controltyp     /* type of control algorithm */
          );
void conequ_printhead(int kstep, NR_CONTROLTYP  controltyp, int cdof, double csp);
void conequ_printiter(int itnum, double disval, double rlnew, double dinorm,
                     double renorm, double energy, double dnorm, double rrnorm);
/*----------------------------------------------------------------------*
 | global_calelm.c                                       m.gee 11/01    |
 *----------------------------------------------------------------------*/
void calelm(FIELD        *actfield,     /* active field */        
            SOLVAR       *actsolv,      /* active SOLVAR */
            PARTITION    *actpart,      /* my partition of this field */
            INTRA        *actintra,     /* my intra-communicator */
            int           sysarray1,    /* number of first sparse system matrix */
            int           sysarray2,    /* number of secnd system matrix, if present, else -1 */
            double       *dvec,         /* global redundant vector passed to elements */
            int           global_numeq, /* size of dvec */
            int           kstep,        /* time in increment step we are in */
            CALC_ACTION  *action);       /* calculation option passed to element routines */

void calinit(FIELD       *actfield,   /* the active physical field */ 
             PARTITION   *actpart,    /* my partition of this field */
             CALC_ACTION *action);
/*----------------------------------------------------------------------*
 | global_calrhs.c                                       m.gee 11/01    |
 *----------------------------------------------------------------------*/
void calrhs(FIELD        *actfield,     /* the active field */
            SOLVAR       *actsolv,      /* the active SOLVAR */
            PARTITION    *actpart,      /* my partition of this field */
            INTRA        *actintra,     /* the field's intra-communicator */
            int           actsysarray,  /* the active sparse array */
            DIST_VECTOR  *rhs1,         /* 2 dist. vectors for rhs */
            DIST_VECTOR  *rhs2,
            int           kstep,        /* actual time or load incremental step */
            CALC_ACTION  *action);       /* action to be passed to element routines */
void calrhs_nodal_neumann(
                             PARTITION    *actpart,
                             INTRA        *actintra,
                             SPARSE_TYP   *sysarraytyp,
                             SPARSE_ARRAY *sysarray,
                             DIST_VECTOR  *rhs
                            );
void calrhs_ele_neumann(FIELD        *actfield,    /* the active field */
                        SOLVAR       *actsolv,     
                        PARTITION    *actpart,
                        INTRA        *actintra,
                        SPARSE_TYP   *sysarraytyp,
                        SPARSE_ARRAY *sysarray,
                        int           actsysarray,
                        DIST_VECTOR  *rhs,
                        int           kstep,
                        CALC_ACTION  *action); /* action to be passsed to element routines */
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
void rc_ptr_red_dof_connect(FIELD        *actfield, 
                            PARTITION    *actpart, 
                            SOLVAR       *actsolv,
                            INTRA        *actintra,
                            RC_PTR       *rc_ptr,
                            int         **dof_connect);
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
                       int          **dof_connect,
                       int           *bindx);
void  rc_ptr_make_sparsity(RC_PTR        *rc_ptr,
                           int           *bindx);
/*----------------------------------------------------------------------*
 |  solver_add_data.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void assemble(
                 int                    sysarray1, /* number of first sparse system matrix */
                 struct _ARRAY         *elearray1, /* pointer to first dense element matrix */
                 int                    sysarray2, /* number of first sparse system matrix or -1 if not given */
                 struct _ARRAY         *elearray2, /* pointer to second dense element matrix or NULL is not present*/
                 struct _PARTITION     *actpart,   /* my partition of theactive field */
                 struct _SOLVAR        *actsolv,   /* the active SOLVAR */
                 struct _INTRA         *actintra,  /* the active intracommunicator */
                 struct _ELEMENT       *actele,    /* the element to assemble */
                 enum _ASSEMBLE_ACTION  assemble_action  /* the assembly option */
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
