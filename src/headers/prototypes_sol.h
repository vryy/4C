/*----------------------------------------------------------------------*
 | cal_nlnstatic_control.c                               m.gee 11/01    |
 *----------------------------------------------------------------------*/
void stanln(void);
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
void increment_controlarc(INTRA         *actintra,
                          SOLVAR        *actsolv,
                          int            actsysarray,
                          DIST_VECTOR   *rsd,
                          DIST_VECTOR   *dispi,
                          DIST_VECTOR   *disp,
                          double         rlnew,
                          double         rlold, 
                          double         stepsize,  
                          double        *rli);
void increment_controldisp(INTRA *actintra,
                          SOLVAR        *actsolv,
                          int            actsysarray,
                          int            cdof,
                          DIST_VECTOR   *rsd,
                          DIST_VECTOR   *dispi,   
                          double        *rli);
/*----------------------------------------------------------------------*
 | cal_nlndyn_struct.c                                   m.gee 02/02    |
 *----------------------------------------------------------------------*/
void dyn_nln_structural(void); 
/*----------------------------------------------------------------------*
 | dyn_service.c                                         m.gee 02/02    |
 *----------------------------------------------------------------------*/
void kefnln_struct(STRUCT_DYN_CALC *dynvar, 
                  STRUCT_DYNAMIC  *sdyn,
                  FIELD           *actfield,
                  SOLVAR          *actsolv,
                  INTRA           *actintra,
                  DIST_VECTOR     *work,
                  int              stiff_array,
                  int              mass_array,
                  int              damp_array);
void pefnln_struct(STRUCT_DYN_CALC *dynvar, 
                  STRUCT_DYNAMIC  *sdyn,
                  FIELD           *actfield,
                  SOLVAR          *actsolv,
                  INTRA           *actintra,
                  DIST_VECTOR     *dispi,
                  DIST_VECTOR     *vel,
                  DIST_VECTOR     *acc,
                  DIST_VECTOR     *work,
                  int              mass_array,
                  int              damp_array);
void dynnle(STRUCT_DYN_CALC *dynvar, STRUCT_DYNAMIC *sdyn, INTRA *actintra, SOLVAR *actsolv, 
           DIST_VECTOR *dispi, /* converged incremental displacements */
           DIST_VECTOR *fie1,  /* internal forces at time t-dt */
           DIST_VECTOR *fie2,  /* internal forces at time t */
           DIST_VECTOR *rhs1,  /* load at time t                      */ 
           DIST_VECTOR *rhs2,  /* load at time t-dt                   */ 
           DIST_VECTOR *work0);
void dyne(STRUCT_DYN_CALC *dynvar,
         INTRA           *actintra,
         SOLVAR          *actsolv,
         int              mass_array,
         DIST_VECTOR     *vel,
         DIST_VECTOR     *work);
void dyn_setconstants(STRUCT_DYN_CALC *dynvar, STRUCT_DYNAMIC *sdyn, double dt);
void dyn_nlnstructupd(STRUCT_DYN_CALC *dynvar, STRUCT_DYNAMIC *sdyn, SOLVAR *actsolv,
                     DIST_VECTOR *sol_old, DIST_VECTOR *sol_new,
                     DIST_VECTOR *rhs_new, DIST_VECTOR *rhs_old,
                     DIST_VECTOR *vel,     DIST_VECTOR *acc,
                     DIST_VECTOR *work0,   DIST_VECTOR *work1,
                     DIST_VECTOR *work2);
void dyn_nlnstruct_outhead(STRUCT_DYN_CALC *dynvar, STRUCT_DYNAMIC *sdyn);
void dyn_nlnstruct_outstep(STRUCT_DYN_CALC *dynvar, STRUCT_DYNAMIC *sdyn, int numiter);
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
void calreduce(FIELD       *actfield, /* the active field */
               PARTITION   *actpart,  /* my partition of this field */
               INTRA       *actintra, /* the field's intra-communicator */
               CALC_ACTION *action,   /* action for element routines */
               int          kstep);    /* the actual time or incremental step */
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
 |  global_mask_skyline.c                                m.gee 02/02    |
 *----------------------------------------------------------------------*/
void mask_skyline(FIELD         *actfield, 
                  PARTITION     *actpart, 
                  SOLVAR        *actsolv,
                  INTRA         *actintra, 
                  SKYMATRIX     *sky);
void  skyline_update(FIELD         *actfield, 
                    PARTITION     *actpart, 
                    SOLVAR        *actsolv,
                    INTRA         *actintra,
                    SKYMATRIX     *sky);
void  skyline_nnz_topology(FIELD      *actfield, 
                         PARTITION    *actpart, 
                         SOLVAR       *actsolv,
                         INTRA        *actintra,
                         SKYMATRIX    *sky,
                         int         **dof_connect);
void   skyline_make_red_dof_connect(FIELD         *actfield, 
                                   PARTITION     *actpart, 
                                   SOLVAR        *actsolv,
                                   INTRA         *actintra,
                                   SKYMATRIX     *sky,
                                   int          **dof_connect,
                                   ARRAY         *red_dof_connect);
void  skyline_make_sparsity(SKYMATRIX  *sky, ARRAY *red_dof_connect);
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
                  struct _DENSE         *dense1,
                  struct _DENSE         *dense2);
void redundant_dense(
                        PARTITION     *actpart,
                        SOLVAR        *actsolv,
                        INTRA         *actintra,
                        DENSE         *dense1,
                        DENSE         *dense2
                        );
/*----------------------------------------------------------------------*
 |  solver_add_msr.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void  add_msr(struct _PARTITION     *actpart,
                struct _SOLVAR        *actsolv,
                struct _INTRA         *actintra,
                struct _ELEMENT       *actele,
                struct _AZ_ARRAY_MSR  *msr1,
                struct _AZ_ARRAY_MSR  *msr2);
void add_msr_checkcouple(int ii,int **cdofs,int ncdofs,int *iscouple,
                           int *isowner, int nprocs);
void add_msr_sendbuff(int ii,int jj,int i,int j,int ii_owner,int **isend,
                    double **dsend,double **estif, int numsend);
void exchange_coup_msr(
                         PARTITION     *actpart,
                         SOLVAR        *actsolv,
                         INTRA         *actintra,
                         AZ_ARRAY_MSR  *msr);
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
 |  solver_add_skyline.c                                 m.gee 09/01    |
 *----------------------------------------------------------------------*/
void  add_skyline(struct _PARTITION     *actpart,
                  struct _SOLVAR        *actsolv,
                  struct _INTRA         *actintra,
                  struct _ELEMENT       *actele,
                  struct _SKYMATRIX     *sky1,
                  struct _SKYMATRIX     *sky2);
void redundant_skyline(
                        PARTITION     *actpart,
                        SOLVAR        *actsolv,
                        INTRA         *actintra,
                        SKYMATRIX     *sky1,
                        SKYMATRIX     *sky2
                        );
/*----------------------------------------------------------------------*
 |  solver_add_rc_ptr.c                                  m.gee 02/02    |
 *----------------------------------------------------------------------*/
void  add_rc_ptr(struct _PARTITION     *actpart,
                struct _SOLVAR        *actsolv,
                struct _INTRA         *actintra,
                struct _ELEMENT       *actele,
                struct _RC_PTR        *rc_ptr);
void add_rcptr_sendbuff(int ii,int jj,int i,int j,int ii_owner,int **isend,
                    double **dsend,double **estif, int numsend);
void exchange_coup_rc_ptr(
                         PARTITION     *actpart,
                         SOLVAR        *actsolv,
                         INTRA         *actintra,
                         RC_PTR        *rc_ptr
                        );
/*----------------------------------------------------------------------*
 |  solver_colsol.c                                       m.gee 02/02    |
 *----------------------------------------------------------------------*/
void solver_colsol(struct _SOLVAR         *actsolv,
                   struct _INTRA          *actintra,
                   struct _SKYMATRIX      *sky,
                   struct _DIST_VECTOR    *sol,
                   struct _DIST_VECTOR    *rhs,
                   int                     option);
/*----------------------------------------------------------------------*
 |  solver_hypre.c                                       m.gee 02/02    |
 *----------------------------------------------------------------------*/
void  solver_hypre_parcsr(struct _SOLVAR         *actsolv,
                          struct _INTRA          *actintra,
                          struct _H_PARCSR       *parcsr,
                          struct _DIST_VECTOR    *sol,
                          struct _DIST_VECTOR    *rhs,
                          int                     option);
/*----------------------------------------------------------------------*
 |  solver_mumps.c                                       m.gee 02/02    |
 *----------------------------------------------------------------------*/
void solver_mumps(struct _SOLVAR         *actsolv,
                  struct _INTRA          *actintra,
                  struct _RC_PTR         *rc_ptr,
                  struct _DIST_VECTOR    *sol,
                  struct _DIST_VECTOR    *rhs,
                  int                     option);
/*----------------------------------------------------------------------*
 |  solver_aztec.c                                       m.gee 11/01    |
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
 |  solver_control.c                                     m.gee 11/01    |
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
 |  solver_hypre.c                                       m.gee 11/01    |
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
 |  solver_lapack.c                                      m.gee 11/01    |
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
 |  solver_service.c                                     m.gee 11/01    |
 *----------------------------------------------------------------------*/
/* A = A * factor */
void solserv_scal_mat(SPARSE_TYP *Atyp,SPARSE_ARRAY *A,double factor);
/* A = A + B * factor */
void solserv_add_mat(INTRA *actintra,SPARSE_TYP *Atyp,SPARSE_ARRAY *A,SPARSE_TYP *Btyp,SPARSE_ARRAY *B,double factor);
/* extracts numeq and numeq_total from ditributed sparse matrix */
void solserv_getmatdims(SPARSE_ARRAY mat,SPARSE_TYP mattyp,int *numeq, int *numeq_total);
/* initializes matrix by zero */
void solserv_zero_mat(INTRA *actintra,SPARSE_ARRAY *mat,SPARSE_TYP *mattyp);
/* copies the mask of a sparse matrix, space is allocated */
void solserv_alloc_cp_sparsemask(INTRA *actintra,SPARSE_TYP *typfrom,SPARSE_ARRAY *matfrom,SPARSE_TYP *typto,SPARSE_ARRAY *matto);
/*internal routine called by solserv_alloc_cp_sparsemask */
void solserv_cp_rc_ptrmask(INTRA *actintra, RC_PTR *from, RC_PTR *to);
/*internal routine called by solserv_alloc_cp_sparsemask */
void solserv_cp_ucchbmask(UCCHB *from, UCCHB *to);
/*internal routine called by solserv_alloc_cp_sparsemask */
void solserv_cp_skymask(SKYMATRIX *from, SKYMATRIX *to);
/*internal routine called by solserv_alloc_cp_sparsemask */
void solserv_cp_msrmask(AZ_ARRAY_MSR *from, AZ_ARRAY_MSR *to);
/*internal routine called by solserv_alloc_cp_sparsemask */
void solserv_cp_densemask(DENSE *from, DENSE *to);
/* performs matrix vector product */
void solserv_sparsematvec(INTRA *actintra,DIST_VECTOR *result,SPARSE_ARRAY *mat,SPARSE_TYP *mattyp,DIST_VECTOR *vec);
/*internal routine called by solserv_sparsematvec */
void solserv_matvec_msr(INTRA *actintra,AZ_ARRAY_MSR *msr,double *work1,double *work2);
/*internal routine called by solserv_sparsematvec */
void solserv_matvec_sky(INTRA *actintra,SKYMATRIX *sky,double *work1,double *work2);
/*internal routine called by solserv_sparsematvec */
void solserv_matvec_dense(INTRA *actintra,DENSE *dense,double *work1,double *work2);

/*----------------------------------------------------------------------*
 |  solver_service2.c                                    m.gee 02/02    |
 *----------------------------------------------------------------------*/
/* create and allocate a vector of DIST_VECTORS */
void solserv_create_vec(DIST_VECTOR **vector,int numvectors,int numeq_total,
                        int numeq,char typstr[]);
/* delete a vector of DIST_VECTORS */
void solserv_del_vec(DIST_VECTOR **vector,int numvectors);
/* init a DIST_VECTOR by zero */
void solserv_zero_vec(DIST_VECTOR *disvector);
/* perform a = a *  b * factor */
void solserv_add_vec(DIST_VECTOR *vec_from,DIST_VECTOR *vec_to,double factor);
/* copy a to b */
void solserv_copy_vec(DIST_VECTOR *vec_from,DIST_VECTOR *vec_to);
/* make euclidian vector norm of a DIST_VECTOR */
void solserv_vecnorm_euclid(INTRA *actintra,DIST_VECTOR *dist_vec,double *result);
/* make Linf norm of a vector (absolute maximum value) */
void solserv_vecnorm_Linf(INTRA *actintra,DIST_VECTOR *dist_vec,double *result);
/* extract redundantly a certain entry a[i] with a given i */
void solserv_getele_vec(INTRA*actintra,SPARSE_TYP *sysarray_typ,
                        SPARSE_ARRAY *sysarray,DIST_VECTOR *dist_vec,
                        int indiz,double *result);
/* perform scalar = a * b */
void solserv_dot_vec(INTRA *actintra,DIST_VECTOR *dist_vec1,
                     DIST_VECTOR *dist_vec2,double *dot);
/* perform a = a * scalar */
void solserv_scalarprod_vec(DIST_VECTOR *dist_vec,double scalar);
/* extract a full-sizez redundant vector from a distributed vector */
void solserv_reddistvec(DIST_VECTOR *distvec,SPARSE_ARRAY *sysarray,
                        SPARSE_TYP *sysarray_typ,double *fullvec,
                        int dim,INTRA *actintra);
/* make a disrtibuted vector (by a format given through a sparse matrix) from a redundant full vector */
void solserv_distribdistvec(DIST_VECTOR  *distvec,SPARSE_ARRAY *sysarray,
                            SPARSE_TYP *sysarray_typ,double *fullvec,
                            int dim,INTRA *actintra);
/* returns values a[i] to the structures NODE.sol in a certain row of sol */
void solserv_result_total(FIELD *actfield,INTRA *actintra,DIST_VECTOR *sol,
                          int place,SPARSE_ARRAY *sysarray,SPARSE_TYP *sysarray_typ);

/* returns values a[i] to the structures NODE.sol_increment in a certain row of sol_increment */
void solserv_result_incre(FIELD *actfield,INTRA *actintra,DIST_VECTOR *sol,
                          int place,SPARSE_ARRAY *sysarray,SPARSE_TYP *sysarray_typ);

/* returns values a[i] to the structures NODE.sol_residual in a certain row of sol_residual */
void solserv_result_resid(FIELD *actfield,INTRA *actintra,DIST_VECTOR *sol,
                          int place,SPARSE_ARRAY *sysarray,SPARSE_TYP *sysarray_typ);
/*----------------------------------------------------------------------*
 |  solver_superlu.c                                     m.gee 11/01    |
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
