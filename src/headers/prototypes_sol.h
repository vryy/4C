/*!---------------------------------------------------------------------
\file
\brief 

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

---------------------------------------------------------------------*/
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
            INT            kstep,        /* the load or time step we are in */
            INT            actsysarray,  /* number of the system matrix in actsolv->sysarray[actsysarray] to be used */
            DIST_VECTOR   *rsd,          /* dist. vector of incremental residual forces used for iteration in conequ */
            DIST_VECTOR   *dispi,        /* dist. vector of incremental displacements */
            INT            cdof,         /* number of the dof to be controlled */
            STANLN        *nln_data,     /* data of the Newton-Raphson method */
            NR_CONTROLTYP  controltyp,    /* type of control algorithm */
            CONTAINER     *container     /*!< contains variables defined in container.h */
          );

void conequ(
            FIELD         *actfield,      /* the actual physical field */
            SOLVAR        *actsolv,       /* the field-corresponding solver */
            PARTITION     *actpart,       /* the partition of the proc */
            INTRA         *actintra,      /* the intra-communicator of this field */
            CALC_ACTION   *action,        /* calculation flag */
            INT            kstep,         /* the load or time step we are in */
            INT           *itnum,         /* number of corrector steps taken by this routine */
            INT            actsysarray,   /* number of the system matrix in actsolv->sysarray[actsysarray] to be used */
            DIST_VECTOR   *rsd,           /* dist. vector of incremental residual forces used for iteration in conequ */
            DIST_VECTOR   *dispi,         /* dist. vector of incremental displacements */
            DIST_VECTOR   *re,            /* re[0..2] 3 vectors for residual displacements */
            INT            cdof,          /* number of dof to be controlled */
            INT           *reldof,        /* numbers of dofs for output */
            STANLN        *nln_data,      /* data of the Newton-Raphson method */
            NR_CONTROLTYP  controltyp,     /* type of control algorithm */
            CONTAINER     *container      /*!< contains variables defined in container.h */
          );
          
void conequ_printhead(INT kstep, NR_CONTROLTYP  controltyp, INT cdof, DOUBLE csp);
void conequ_printiter(INT itnum, DOUBLE disval, DOUBLE rlnew, DOUBLE dinorm,
                     DOUBLE renorm, DOUBLE energy, DOUBLE dnorm, DOUBLE rrnorm);
void increment_controlarc(INTRA         *actintra,
                          SOLVAR        *actsolv,
                          INT            actsysarray,
                          DIST_VECTOR   *rsd,
                          DIST_VECTOR   *dispi,
                          DIST_VECTOR   *disp,
                          DOUBLE         rlnew,
                          DOUBLE         rlold, 
                          DOUBLE         stepsize,  
                          DOUBLE        *rli);
void increment_controldisp(INTRA *actintra,
                          SOLVAR        *actsolv,
                          INT            actsysarray,
                          INT            cdof,
                          DIST_VECTOR   *rsd,
                          DIST_VECTOR   *dispi,   
                          DOUBLE        *rli);
/*----------------------------------------------------------------------*
 | cal_nlndyn_struct.c                                   m.gee 02/02    |
 *----------------------------------------------------------------------*/
void dyn_nln_structural(void); 
/*----------------------------------------------------------------------*
 | cal_nlndyn_stru_expl.c                                m.gee 05/02    |
 *----------------------------------------------------------------------*/
void dyn_nln_stru_expl(void); 
/*----------------------------------------------------------------------*
 | dyn_service.c                                         m.gee 02/02    |
 *----------------------------------------------------------------------*/
void dyn_eigen(FIELD *actfield, PARTITION *actpart, SOLVAR *actsolv,
               INTRA *actintra, INT stiff, INT mass);
void kefnln_struct(STRUCT_DYN_CALC *dynvar, 
                  STRUCT_DYNAMIC  *sdyn,
                  FIELD           *actfield,
                  SOLVAR          *actsolv,
                  INTRA           *actintra,
                  DIST_VECTOR     *work,
                  INT              stiff_array,
                  INT              mass_array,
                  INT              damp_array);
void dyn_keff_expl(INTRA *actintra,
                   SPARSE_TYP *sysarray_typ, SPARSE_ARRAY *sysarray,
                   INT stiff_array, INT mass_array, INT damp_array,
                   STRUCT_DYN_CALC *dynvar, STRUCT_DYNAMIC *sdyn);
void pefnln_struct(STRUCT_DYN_CALC *dynvar, 
                  STRUCT_DYNAMIC  *sdyn,
                  FIELD           *actfield,
                  SOLVAR          *actsolv,
                  INTRA           *actintra,
                  DIST_VECTOR     *dispi,
                  DIST_VECTOR     *vel,
                  DIST_VECTOR     *acc,
                  DIST_VECTOR     *work,
                  INT              mass_array,
                  INT              damp_array);
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
         INT              mass_array,
         DIST_VECTOR     *vel,
         DIST_VECTOR     *work);
void dyn_setconstants(STRUCT_DYN_CALC *dynvar, STRUCT_DYNAMIC *sdyn, DOUBLE dt);
void dyn_setconstants_expl(STRUCT_DYN_CALC *dynvar, STRUCT_DYNAMIC *sdyn, DOUBLE dt);
void dyn_nlnstructupd(FIELD *actfield,      STRUCT_DYN_CALC *dynvar, 
                      STRUCT_DYNAMIC *sdyn, SOLVAR *actsolv,
                      DIST_VECTOR *sol_old, DIST_VECTOR *sol_new,
                      DIST_VECTOR *rhs_new, DIST_VECTOR *rhs_old,
                      DIST_VECTOR *vel,     DIST_VECTOR *acc,
                      DIST_VECTOR *work0,   DIST_VECTOR *work1,
                      DIST_VECTOR *work2);
void dyn_nlnstruct_outhead(STRUCT_DYN_CALC *dynvar, STRUCT_DYNAMIC *sdyn);
void dyn_nlnstru_outhead_expl(void);
void dyn_nlnstruct_outstep(STRUCT_DYN_CALC *dynvar, STRUCT_DYNAMIC *sdyn, INT numiter, DOUBLE dt);
void assemble_dirich_dyn(ELEMENT *actele, ARRAY *estif_global, 
                         ARRAY *emass_global, CONTAINER *container);
void dyn_epot(FIELD *actfield, INT disnum, INTRA *actintra, STRUCT_DYN_CALC *dynvar, DOUBLE *deltaepot);
void dyn_ekin(FIELD *actfield, SOLVAR *actsolv, PARTITION *actpart, INTRA *actintra, CALC_ACTION *action,
             CONTAINER *container, INT stiff_array, INT mass_array);
void dyn_eout(STRUCT_DYN_CALC *dynvar, STRUCT_DYNAMIC *sdyn, INTRA *actintra, SOLVAR *actsolv, 
              DIST_VECTOR *dispi, /* converged incremental displacements */
              DIST_VECTOR *rhs1,  /* load at time t                      */ 
              DIST_VECTOR *rhs2,  /* load at time t-dt                   */ 
              DIST_VECTOR *work0);
void dyn_ekin_local(ELEMENT *actele,ARRAY *emass, CONTAINER  *container);

/*----------------------------------------------------------------------*
 | global_calelm.c                                       m.gee 11/01    |
 *----------------------------------------------------------------------*/
void calelm(FIELD        *actfield,     /* active field */        
            SOLVAR       *actsolv,      /* active SOLVAR */
            PARTITION    *actpart,      /* my partition of this field */
            INTRA        *actintra,     /* my intra-communicator */
            INT           sysarray1,    /* number of first sparse system matrix */
            INT           sysarray2,    /* number of secnd system matrix, if present, else -1 */
            CONTAINER    *container,    /*!< contains variables defined in container.h */
            CALC_ACTION  *action);       /* calculation option passed to element routines */            
            
void calinit(FIELD       *actfield,   /* the active physical field */ 
             PARTITION   *actpart,    /* my partition of this field */
             CALC_ACTION *action,
             CONTAINER   *container); /*!< contains variables defined in container.h */

void calreduce(FIELD       *actfield, /* the active field */
               PARTITION   *actpart,  /* my partition of this field */
               INTRA       *actintra, /* the field's intra-communicator */
               CALC_ACTION *action,   /* action for element routines */
               CONTAINER   *container); /*!< contains variables defined in container.h */               
/*----------------------------------------------------------------------*
 | global_calrhs.c                                       m.gee 11/01    |
 *----------------------------------------------------------------------*/
void calrhs(FIELD        *actfield,     /* the active field */
            SOLVAR       *actsolv,      /* the active SOLVAR */
            PARTITION    *actpart,      /* my partition of this field */
            INTRA        *actintra,     /* the field's intra-communicator */
            INT           actsysarray,  /* the active sparse array */
            DIST_VECTOR  *rhs1,         /* 2 dist. vectors for rhs */
            CALC_ACTION  *action,       /* action to be passed to element routines */
            CONTAINER    *container);    /*!< contains variables defined in container.h */
void rhs_point_neum(DOUBLE *rhs, INT dimrhs, PARTITION *actpart);     
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
                    INT          *numeq);
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
                            INT         **dof_connect);
void parcsr_make_bindx(FIELD         *actfield, 
                          PARTITION     *actpart, 
                          SOLVAR        *actsolv,
                          H_PARCSR      *parcsr,
                          INT          **dof_connect);
/*----------------------------------------------------------------------*
 | global_mask_mds.c                                        al 02/03    |
 *----------------------------------------------------------------------*/
void mask_mds(FIELD        *actfield, 
              PARTITION    *actpart, 
              SOLVAR       *actsolv,
              INTRA        *actintra, 
              ML_ARRAY_MDS *mds);
/*----------------------------------------------------------------------*
 | global_mask_msr.c                                   m.gee 11/01    |
 *----------------------------------------------------------------------*/
void mask_msr(FIELD         *actfield, 
              PARTITION     *actpart, 
              SOLVAR        *actsolv,
              INTRA         *actintra, 
              AZ_ARRAY_MSR  *msr,
	      INT            actdis);
void msr_numeq(FIELD         *actfield, 
                  PARTITION    *actpart, 
                  SOLVAR       *actsolv,
                  INTRA        *actintra,
                  INT          *numeq);
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
                         INT         **dof_connect);
void dof_in_coupledofs(INT dof, PARTITION *actpart, INT *iscoupled);
void dof_find_centernode(INT dof, PARTITION *actpart, NODE **centernode);
void msr_make_bindx(FIELD         *actfield, 
                       PARTITION     *actpart, 
                       SOLVAR        *actsolv,
                       AZ_ARRAY_MSR  *msr,
                       INT          **dof_connect);
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
                   INT          *numeq);
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
                           INT       **dof_connect);
void  ucchb_make_a(FIELD         *actfield, 
                     PARTITION     *actpart, 
                     SOLVAR        *actsolv,
                     INTRA         *actintra,
                     UCCHB         *ucchb,
                     INT          **dof_connect);
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
                            INT         **dof_connect);
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
                         INT         **dof_connect);
void  rc_ptr_make_bindx(FIELD         *actfield, 
                       PARTITION     *actpart, 
                       SOLVAR        *actsolv,
                       RC_PTR        *rc_ptr,
                       INT          **dof_connect,
                       INT           *bindx);
void  rc_ptr_make_sparsity(RC_PTR        *rc_ptr,
                           INT           *bindx);
/*----------------------------------------------------------------------*
 |  global_mask_ccf.c                              s.offermanns 04/02    |
 *----------------------------------------------------------------------*/
void mask_ccf(FIELD         *actfield, 
              PARTITION     *actpart, 
              SOLVAR        *actsolv,
              INTRA         *actintra, 
              CCF        *ccf);
void  ccf_red_dof_connect(FIELD        *actfield, 
                          PARTITION    *actpart,
                          SOLVAR       *actsolv,
                          INTRA        *actintra,
                          CCF          *ccf,
                          INT         **dof_connect,
                          ARRAY        *red_dof_connect);
void  ccf_update(FIELD         *actfield, 
                 PARTITION     *actpart, 
                 SOLVAR        *actsolv,
                 INTRA         *actintra,
                 CCF           *ccf);
void  ccf_nnz_topology(FIELD        *actfield, 
                       PARTITION    *actpart, 
                       SOLVAR       *actsolv,
                       INTRA        *actintra,
                       CCF          *ccf,
                       INT         **dof_connect);
void  ccf_make_bindx(FIELD         *actfield, 
                     PARTITION     *actpart, 
                     SOLVAR        *actsolv,
                     CCF           *ccf,
                     INT           *bindx,
                     ARRAY         *red_dof_connect);
void  ccf_make_sparsity(CCF        *ccf,
                        INT        *bindx);
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
                         INT         **dof_connect);
void   skyline_make_red_dof_connect(FIELD         *actfield, 
                                   PARTITION     *actpart, 
                                   SOLVAR        *actsolv,
                                   INTRA         *actintra,
                                   SKYMATRIX     *sky,
                                   INT          **dof_connect,
                                   ARRAY         *red_dof_connect);
void  skyline_make_sparsity(SKYMATRIX  *sky, ARRAY *red_dof_connect);
/*----------------------------------------------------------------------*
 |  global_mask_spooles.c                                m.gee 05/02    |
 *----------------------------------------------------------------------*/
void mask_spooles(FIELD         *actfield, 
                  PARTITION     *actpart, 
                  SOLVAR        *actsolv,
                  INTRA         *actintra, 
                  SPOOLMAT      *spo);
void     spo_make_sparsity(SPOOLMAT        *spo,
                           INT           *bindx);
void    spo_make_bindx(FIELD         *actfield, 
                       PARTITION     *actpart, 
                       SOLVAR        *actsolv,
                       SPOOLMAT      *spo,
                       INT          **dof_connect,
                       INT           *bindx);
void  spo_nnz_topology(FIELD         *actfield, 
                       PARTITION    *actpart, 
                       SOLVAR       *actsolv,
                       INTRA        *actintra,
                       SPOOLMAT     *spo,
                       INT         **dof_connect);
void spo_update(FIELD         *actfield, 
                PARTITION     *actpart, 
                SOLVAR        *actsolv,
                INTRA         *actintra,
                SPOOLMAT      *spo);
/*----------------------------------------------------------------------*
 |  solver_add_data.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void assemble(
                 INT                    sysarray1, /* number of first sparse system matrix */
                 struct _ARRAY         *elearray1, /* pointer to first dense element matrix */
                 INT                    sysarray2, /* number of first sparse system matrix or -1 if not given */
                 struct _ARRAY         *elearray2, /* pointer to second dense element matrix or NULL is not present*/
                 struct _PARTITION     *actpart,   /* my partition of theactive field */
                 struct _SOLVAR        *actsolv,   /* the active SOLVAR */
                 struct _INTRA         *actintra,  /* the active intracommunicator */
                 struct _ELEMENT       *actele,    /* the element to assemble */
                 enum _ASSEMBLE_ACTION  assemble_action, /* the assembly option */
                 CONTAINER             *container);  /*!< contains variables defined in container.h */                 
void init_assembly(
                       struct _PARTITION      *actpart,
                       struct _SOLVAR         *actsolv,
                       struct _INTRA          *actintra,
                       struct _FIELD          *actfield,
                       INT                     actsysarray,
		       INT                     actndis
                     );
void assemble_vec(INTRA        *actintra,
                    SPARSE_TYP   *sysarraytyp,
                    SPARSE_ARRAY *sysarray,
                    DIST_VECTOR  *rhs,
                    DOUBLE       *drhs,
                    DOUBLE        factor);
void assemble_intforce(ELEMENT *actele, ARRAY *elevec_a, CONTAINER *container, INTRA *actintra);
void assemble_dirich(ELEMENT *actele,ARRAY *estif_global,CONTAINER *container);
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
void add_msr_checkcouple(INT ii,INT **cdofs,INT ncdofs,INT *iscouple,
                           INT *isowner, INT nprocs);
void add_msr_sendbuff(INT ii,INT jj,INT i,INT j,INT ii_owner,INT **isend,
                    DOUBLE **dsend,DOUBLE **estif, INT numsend);
void exchange_coup_msr(
                         PARTITION     *actpart,
                         SOLVAR        *actsolv,
                         INTRA         *actintra,
                         AZ_ARRAY_MSR  *msr);
/*----------------------------------------------------------------------*
 |  solver_add_mlib.c                                    m.gee 11/01    |
 *----------------------------------------------------------------------*/
INT  add_mds(
    struct _PARTITION     *actpart,
    struct _SOLVAR        *actsolv,
    struct _ELEMENT       *actele,
    struct _ML_ARRAY_MDS  *mds);

/*----------------------------------------------------------------------*
 | solver_add_parcsr.c                                   m.gee 11/01    |
 *----------------------------------------------------------------------*/
void  add_parcsr(struct _PARTITION     *actpart,
                   struct _SOLVAR        *actsolv,
                   struct _INTRA         *actintra,
                   struct _ELEMENT       *actele,
                   struct _H_PARCSR      *parcsr);
void add_parcsr_sendbuff(INT ii,INT jj,INT i,INT j,INT ii_owner,INT **isend,
                    DOUBLE **dsend,DOUBLE **estif, INT numsend);
void add_parcsr_checkcouple(INT ii,INT **cdofs,INT ncdofs,INT *iscouple,INT *isowner, INT nprocs);
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
                struct _RC_PTR        *rc_ptr1,
                struct _RC_PTR        *rc_ptr2);
void add_rcptr_sendbuff(INT ii,INT jj,INT i,INT j,INT ii_owner,INT **isend,
                    DOUBLE **dsend,DOUBLE **estif, INT numsend);
void exchange_coup_rc_ptr(
                         PARTITION     *actpart,
                         SOLVAR        *actsolv,
                         INTRA         *actintra,
                         RC_PTR        *rc_ptr
                        );
/*----------------------------------------------------------------------*
 |  solver_add_ccf.c                              s.offermanns 02/02    |
 *----------------------------------------------------------------------*/
void  add_ccf(struct _PARTITION       *actpart,
                struct _SOLVAR        *actsolv,
                struct _INTRA         *actintra,
                struct _ELEMENT       *actele,
                struct _CCF           *ccf,
                struct _CCF           *ccf2);
void redundant_ccf(struct _PARTITION *actpart,
                   struct _SOLVAR    *actsolv,
                   struct _INTRA     *actintra,
                   struct _CCF       *ccf1,
                   struct _CCF       *ccf2);
/*----------------------------------------------------------------------*
 |  solver_add_spooles.c                                 m.gee 05/02    |
 *----------------------------------------------------------------------*/
void  add_spo(struct _PARTITION     *actpart,
              struct _SOLVAR        *actsolv,
              struct _INTRA         *actintra,
              struct _ELEMENT       *actele,
              struct _SPOOLMAT      *spo1,
              struct _SPOOLMAT      *spo2);
void add_spo_sendbuff(INT ii,INT jj,INT i,INT j,INT ii_owner,INT **isend,
                      DOUBLE **dsend,DOUBLE **estif, INT numsend);
void exchange_coup_spo(
                         PARTITION     *actpart,
                         SOLVAR        *actsolv,
                         INTRA         *actintra,
                         SPOOLMAT      *spo
                        );
void add_val_spo(INT ii,INT index, INT jj, struct _SPOOLMAT *spo, DOUBLE val, INTRA *actintra);
void set_val_spo(INT ii,INT index, INT jj, struct _SPOOLMAT *spo, DOUBLE val, INTRA *actintra);
void close_spooles_matrix(struct _SPOOLMAT *spo, INTRA *actintra);
void add_spooles_matrix(struct _SPOOLMAT *to, struct _SPOOLMAT *from,
                       DOUBLE factor, INT init, INTRA *actintra);
/*----------------------------------------------------------------------*
 |  solver_colsol.c                                       m.gee 02/02    |
 *----------------------------------------------------------------------*/
void solver_colsol(struct _SOLVAR         *actsolv,
                   struct _INTRA          *actintra,
                   struct _SKYMATRIX      *sky,
                   struct _DIST_VECTOR    *sol,
                   struct _DIST_VECTOR    *rhs,
                   INT                     option);
/*----------------------------------------------------------------------*
 |  solver_mlib.c                                         m.gee 02/02    |
 *----------------------------------------------------------------------*/
void solver_mlib( 
    struct _SOLVAR         *actsolv,
    struct _INTRA          *actintra,
    struct _ML_ARRAY_MDS   *mds,
    struct _DIST_VECTOR    *sol,
    struct _DIST_VECTOR    *rhs,
    INT                     option
    );

/*----------------------------------------------------------------------*
 |  solver_hypre.c                                       m.gee 02/02    |
 *----------------------------------------------------------------------*/
void  solver_hypre_parcsr(struct _SOLVAR         *actsolv,
                          struct _INTRA          *actintra,
                          struct _H_PARCSR       *parcsr,
                          struct _DIST_VECTOR    *sol,
                          struct _DIST_VECTOR    *rhs,
                          INT                     option);
/*----------------------------------------------------------------------*
 |  solver_mumps.c                                       m.gee 02/02    |
 *----------------------------------------------------------------------*/
void solver_mumps(struct _SOLVAR         *actsolv,
                  struct _INTRA          *actintra,
                  struct _RC_PTR         *rc_ptr,
                  struct _DIST_VECTOR    *sol,
                  struct _DIST_VECTOR    *rhs,
                  INT                     option);
/*----------------------------------------------------------------------*
 |  solver_umfpack.c                              s.offermanns 02/02    |
 *----------------------------------------------------------------------*/
void solver_umfpack(struct _SOLVAR         *actsolv,
                    struct _INTRA          *actintra,
                    struct _CCF            *ccf,
                    struct _DIST_VECTOR    *sol,
                    struct _DIST_VECTOR    *rhs,
                    INT                     option);
/*----------------------------------------------------------------------*
 |  solver_aztec.c                                       m.gee 11/01    |
 *----------------------------------------------------------------------*/
void solver_az_msr( 
                      struct _SOLVAR         *actsolv,
                      struct _INTRA          *actintra,
                      struct _AZ_ARRAY_MSR   *msr_array,
                      struct _DIST_VECTOR    *sol,
                      struct _DIST_VECTOR    *rhs,
                      INT                     option
                     );
/*----------------------------------------------------------------------*
 |  solver_spooles.c                                     m.gee 05/02    |
 *----------------------------------------------------------------------*/
void solver_spooles( 
                      struct _SOLVAR         *actsolv,
                      struct _INTRA          *actintra,
                      struct _SPOOLMAT       *spo,
                      struct _DIST_VECTOR    *sol,
                      struct _DIST_VECTOR    *rhs,
                      INT                     option
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
                       INT                     option
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
                            INT                     option
                           );
void hypre_vector_create(INT              myrank,
                           H_PARCSR        *parcsr,
                           INTRA           *actintra,
                           HYPRE_IJVector  *ijvector);
void hypre_vector_assemble(HYPRE_IJVector  *ijvector,
                             HYPRE_ParVector *parcsr_vector);
void hypre_matrix_create(INT              myrank,
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
                              INT                     option
                             );
/*----------------------------------------------------------------------*
 |  solver_service.c                                     m.gee 11/01    |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 |  make A = A * factor                                      m.gee 02/02|
 | SPARSE_TYP *Atyp (i)        type of sparse matrix                    |
 | SPARSE_ARRAY *A  (i/o)      sparse matrix                            |
 | DOUBLE factor    (i)        factor                                   |
 *----------------------------------------------------------------------*/
void solserv_scal_mat(SPARSE_TYP *Atyp,SPARSE_ARRAY *A,DOUBLE factor);
/*----------------------------------------------------------------------*
 |  make A = A + B * factor                                  m.gee 02/02|
 | INTRA *actintra  (i)     intra-communicator the matrices live on     |
 | SPARSE_TYP *Atyp (i)     type of sparse matrix                       |
 | SPARSE_ARRAY *A  (i/o)   sparse matrix                               |
 | SPARSE_TYP *Btyp (i)     type of sparse matrix                       |
 | SPARSE_ARRAY *B  (i)     sparse matrix                               |
 | DOUBLE factor    (i)     factor                                      |
 *----------------------------------------------------------------------*/
void solserv_add_mat(INTRA *actintra,
                     SPARSE_TYP *Atyp,SPARSE_ARRAY *A,
                     SPARSE_TYP *Btyp, SPARSE_ARRAY *B,
                     DOUBLE factor);
/*----------------------------------------------------------------------*
 |  get dimensions of sparse matrix                          m.gee 02/02|
 | SPARSE_ARRAY mat  (i)     sparse matrix                              |
 | SPARSE_TYP mattyp (i)     type of sparse matrix                      |
 | INT *numeq        (o)     proc-local dimension of sparse matrix      |
 | INT *numeq_total  (o)     global dimension of sparse matrix          |
 *----------------------------------------------------------------------*/
void solserv_getmatdims(SPARSE_ARRAY mat,SPARSE_TYP mattyp,
                        INT *numeq, INT *numeq_total);
/*----------------------------------------------------------------------*
 |  init a distributed matrix to zero - collective call !    m.gee 10/01|
 | INTRA *actintra  (i)     intra-communicator the matrices live on     |
 | SPARSE_ARRAY *mat  (i/o)    sparse matrix                            |
 | SPARSE_TYP *mattyp (i)    type of sparse matrix                      |
 *----------------------------------------------------------------------*/
void solserv_zero_mat(INTRA *actintra,SPARSE_ARRAY *mat,
                                      SPARSE_TYP *mattyp);
/*----------------------------------------------------------------------*
 |  copies the sparsity mask of a system matrix              m.gee 02/02|
 |  for the new sparsity mask, a suitable structure is allocated        |
 | INTRA        *actintra (i)  intra-communicator the matrices live on  |
 | SPARSE_TYP   *typfrom  (i)  type of sparsity mask to be copied       |
 | SPARSE_ARRAY *matfrom  (i)  sparsity mask to be copied               |
 | SPARSE_TYP   *typto    (o)  type of sparsity mask to be copied to    |
 | SPARSE_ARRAY *matto    (o)  sparsity mask to be allocated and copied |
 *----------------------------------------------------------------------*/
void solserv_alloc_cp_sparsemask(INTRA *actintra,
                                 SPARSE_TYP *typfrom,SPARSE_ARRAY *matfrom,
                                 SPARSE_TYP *typto,  SPARSE_ARRAY *matto);
/*----------------------------------------------------------------------*
 |  make matrix vector multiplication                        m.gee 02/02|
 | INTRA        *actintra (i)  intra-communicator the matrices live on  |
 | DIST_VECTOR  *result   (o)  result = mat * vec                       |
 | SPARSE_ARRAY *mat      (i)  sparse matrix                            |
 | SPARSE_TYP   *mattyp   (i)  type of sparse matrix                    |
 | DIST_VECTOR  *vec      (i)  vector to be multiplied with             |
 *----------------------------------------------------------------------*/
void solserv_sparsematvec(INTRA *actintra,DIST_VECTOR *result,
                          SPARSE_ARRAY *mat,SPARSE_TYP *mattyp,
                          DIST_VECTOR *vec);

/*----------------------------------------------------------------------*
 |  solver_service2.c                                    m.gee 02/02    |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 |  create number of distributed vectors - collective call ! m.gee 10/01|
 |  DIST_VECTOR **vector (i/o) adress of pointer a vector of            |
 |                             DIST_VECTORs will be allocated to        |
 |  INT numvectors       (i)   number of DIST_VECTORs to allocate       |
 |  INT numeq_total      (i)   proc-global dimension of the DIST_VECTORs|
 |  INT numeq            (i)   proc_local  dimension of the DIST_VECTORs|
 |  char typstr[]        (i)   ="DV" for DOUBLE-DIST_VECTORs            |
 |  the values in the DIST_VECTORs is NOT initialized                   |
 *----------------------------------------------------------------------*/
void solserv_create_vec(DIST_VECTOR **vector,INT numvectors,
                        INT numeq_total,
                        INT numeq,        char typstr[]);
/*----------------------------------------------------------------------*
 |   delete number of distributed vectors - collective call ! m.gee 2/02|
 |  DIST_VECTOR **vector (i/o) adress of pointer a vector of            |
 |                             DIST_VECTORs is allocated to             |
 |  INT numvectors       (i)   number of DIST_VECTORs to free           |
 |  the routine frees all DIST_VECTORs in vector and sets vector=NULL   |
 *----------------------------------------------------------------------*/
void solserv_del_vec(DIST_VECTOR **vector,INT numvectors);
/*----------------------------------------------------------------------*
 |  init a distributed vector to zero - collective call !    m.gee 10/01|
 |  DIST_VECTOR *disvector (i/o) adress of a DIST_VECTOR to be set to 0.0|
 *----------------------------------------------------------------------*/
void solserv_zero_vec(DIST_VECTOR *disvector);
/*----------------------------------------------------------------------*
 |  add contents of the vector vec_from to vec_to            m.gee 10/01|
 |  vec_to->vec.a.dv[i] += vec_from->vec.a.dv[i]*factor                 |
 |  DIST_VECTOR *vec_from (i)   vector to be added to another vector    |
 |  DIST_VECTOR *vec_to   (i/o) vector to be added to                   |
 |  DOUBLE factor         (i)   scaling factor                          |
 *----------------------------------------------------------------------*/
void solserv_add_vec(DIST_VECTOR *vec_from,DIST_VECTOR *vec_to,
                     DOUBLE factor);
/*----------------------------------------------------------------------*
 |  copy contents of the vector vec_from to vec_to           m.gee 11/01|
 |  vec_to->vec.a.dv[i] = vec_from->vec.a.dv[i]                         |
 |  DIST_VECTOR *vec_from (i)   vector to be copied to another vector   |
 |  DIST_VECTOR *vec_to   (i/o) vector to be copied to                  |
 |  user must assure matching dimensions and types                      |
 *----------------------------------------------------------------------*/
void solserv_copy_vec(DIST_VECTOR *vec_from,DIST_VECTOR *vec_to);
/*----------------------------------------------------------------------*
 |  make euclidian norm of a distributed vector              m.gee 11/01|
 |  result = sqrt( sumof(vec[i]*vec[i]) )                               |
 |  INTRA *actintra (i) intra-communicator the DIST_VECTOR lives on     |
 |  DIST_VECTOR *dist_vec (i) vector to make norm of                    |
 |  DOUBLE *result        (o) norm of the vector                        |
 *----------------------------------------------------------------------*/
void solserv_vecnorm_euclid(INTRA *actintra,DIST_VECTOR *dist_vec,
                            DOUBLE *result);
/*----------------------------------------------------------------------*
 |  find  absolute maximum value in a vector (Linf-Norm)     m.gee 02/02|
 |  *result = MAX( ABS(vec[i]) )                                        |
 |  INTRA *actintra (i) intra-communicator the DIST_VECTOR lives on     |
 |  DIST_VECTOR *dist_vec (i) vector to make norm of                    |
 |  DOUBLE *result        (o) norm of the vector                        |
 *----------------------------------------------------------------------*/
void solserv_vecnorm_Linf(INTRA *actintra,DIST_VECTOR *dist_vec,
                          DOUBLE *result);
/*----------------------------------------------------------------------*
 |  get a certain entry from a distr. vector to all procs    m.gee 11/01|
 |  returns the value of dof indiz in the vector dist_vec on all procs  |
 |  INTRA *actintra (i) intra-communicator the DIST_VECTOR lives on     |
 |  SPARSE_TYP *sysarray_typ (i) sparsity typ of vector-matching matrix |
 |  SPARSE_ARRAY *sysarray   (i) sparse matrix the vector matches in    |
 |                               distribution                           | 
 |  DIST_VECTOR  *dist_vec   (i) vector the value shall be taken from   |
 |  INT           indiz      (i) field-local (unsupported) dof number   |
 |  DOUBLE       *result     (o) value in vector at the given dof       |
 |                               returned redundant on all procs        |
 *----------------------------------------------------------------------*/
void solserv_getele_vec(INTRA*actintra,SPARSE_TYP *sysarray_typ,
                        SPARSE_ARRAY *sysarray,DIST_VECTOR *dist_vec,
                        INT indiz,DOUBLE *result);
/*----------------------------------------------------------------------*
 |  make dot product between 2 distr. vectors                m.gee 11/01|
 |  *dot = sumover_i( vec1[i]*vec2[i] )                                 |
 |  INTRA *actintra (i) intra-communicator the DIST_VECTORs live on     |
 |  DIST_VECTOR *dist_vec1 (i) first vector to be multiplied            |
 |  DIST_VECTOR *dist_vec2 (i) scnd  vector to be multiplied            |
 |  DOUBLE      *dot       (o) result of vector-vector multiplication   |
 |                             returned redundant on all procs          |
 *----------------------------------------------------------------------*/
void solserv_dot_vec(INTRA *actintra,DIST_VECTOR *dist_vec1,
                     DIST_VECTOR *dist_vec2,DOUBLE *dot);
/*----------------------------------------------------------------------*
 |  make product between scalar and distr. vector            m.gee 11/01|
 |  vec[i] = vec[i] * scalar                                            |
 |  DIST_VECTOR *dist_vec (i/o) vector to be multiplied by scalar       |
 |  DOUBLE       scalar   (o)   scalar value                            |
 *----------------------------------------------------------------------*/
void solserv_scalarprod_vec(DIST_VECTOR *dist_vec,DOUBLE scalar);
/*----------------------------------------------------------------------*
 |  Allreduce a distributed vector in an INTRACOMM           m.gee 10/01|
 |  This is a collective call!                                          |
 |  distributed vector to full redundant vector                         |
 |                                                                      |
 |  note that the disributed vectors match a certain type of sparse     |
 |  matrix in the layout of distribution and values. This means, that   |
 |  the value of a certain dof are NOT in distvec->vec.a.dv[dof]!!!!!   |
 |                                                                      |
 |  the redundant vector fullvec holds values of a certain dof in       |
 |  fullvec[dof]                                                        |
 |                                                                      |
 |  the values in the given DIST_VECTOR are copied to a vector of       |
 |  size numeq_total, which is redundant on all procs                   |
 |  DIST_VECTOR *distvec (i) DIST_VECTORto be 'allreduced'              |
 |  SPARSE_ARRAY *sysarray (i) sparse matrix matching the distribution  |
 |                             of distvec                               |
 |  SPARSE_TYP *sysarray_typ (i) type of sparse matrix                  |
 |  DOUBLE *fullvec (o) vector of lenght numeq_total will be holding    |
 |                      the values from distvec in correct dof-ordering:|
 |                      fullvec[dof] = value of a certain dof           |
 |  INTRA *actintra (i) intra-communicator the DIST_VECTOR lives on     |
 *----------------------------------------------------------------------*/
void solserv_reddistvec(DIST_VECTOR *distvec,SPARSE_ARRAY *sysarray,
                        SPARSE_TYP *sysarray_typ,DOUBLE *fullvec,
                        INT dim,INTRA *actintra);
/*----------------------------------------------------------------------*
 |  distribute a full redundant vector                       m.gee 02/02|
 |  This is a collective call!                                          |
 |  full redundant vector to distributed vector                         |
 |  this routine is the inverse of solserv_reddistvec                   |
 |  It copies the values in a vector fullvec of size numeq_total, that  |
 |  is ordered such that fullvec[dof] = value at dof                    |
 |  to a distributed vector matching a certain sparse matrix in         |
 |  distribution of values. Note that in the distvec the values         |
 |  value_at_dof are NOT in distvec->vec.a.dv[dof] !!!!                 |
 |                                                                      |
 |  DIST_VECTOR *distvec (o) DIST_VECTOR to be copied to                |
 |  SPARSE_ARRAY *sysarray (i) sparse matrix matching the distribution  |
 |                             of distvec                               |
 |  SPARSE_TYP *sysarray_typ (i) type of sparse matrix                  |
 |  DOUBLE *fullvec (o) vector of lenght numeq_total  holding           |
 |                      the values  correct dof-ordering:               |
 |                      fullvec[dof] = value of a certain dof           |
 |  INTRA *actintra (i) intra-communicator the DIST_VECTOR lives on     |
 *----------------------------------------------------------------------*/
void solserv_distribdistvec(DIST_VECTOR  *distvec,SPARSE_ARRAY *sysarray,
                            SPARSE_TYP *sysarray_typ,DOUBLE *fullvec,
                            INT dim,INTRA *actintra);
/*----------------------------------------------------------------------*
 |  Put the results of a DIST_VECTOR to the nodes in a       m.gee 10/01|
 |  certain place  in ARRAY sol                                         |
 |  Result has to be allreduced and are put to the whole                |
 |  field on each procs                                                 |
 |  FIELD *actfield (i) the active field                                |
 |  INTRA *actintra (i) intra-communicator the DIST_VECTOR lives on     |
 |  DIST_VECTOR *sol (i) vector of values to be put to the nodes        |
 |  INT place        (i) place in the ARRAY node->sol where to put the  |
 |                       values. Every structure NODE has an ARRAY sol  |
 |                       of type sol.a.da[place][0..numdf-1]            |
 |                       if place >= actual dimensions of the ARRAY sol |
 |                       sol is enlarged                                |
 |  SPARSE_ARRAY *sysarray (i) sparse matrix matching the distribution  |
 |                             of DIST_VECTOR *sol                      |
 |  SPARSE_TYP *sysarray_typ (i) type of sparse matrix                  |
 |                                                                      |
 *----------------------------------------------------------------------*/
void solserv_result_total(FIELD *actfield,INTRA *actintra,DIST_VECTOR *sol,
                          INT place,SPARSE_ARRAY *sysarray,
                          SPARSE_TYP *sysarray_typ);

/*----------------------------------------------------------------------*
 |  Put the results of a DIST_VECTOR to the nodes in a       m.gee 11/01|
 |  certain place in ARRAY sol_increment                                |
 |  Result have to bee allreduced and are put to the whole              |
 |  field on each proc                                                  |
 |  Functionality is the same as in solserv_result_total                |
 *----------------------------------------------------------------------*/
void solserv_result_incre(FIELD *actfield,INTRA *actintra,DIST_VECTOR *sol,
                          int place,SPARSE_ARRAY *sysarray,
                          SPARSE_TYP *sysarray_typ,int ndis);
/*----------------------------------------------------------------------*
 |  Put the results of a DIST_VECTOR to the nodes in a       m.gee 11/01|
 |  certain place in ARRAY sol_residual                                 |
 |  Result have to bee allreduced and are put to the whole              |
 |  field on each proc                                                  |
 |  Functionality is the same as in solserv_result_total                |
 *----------------------------------------------------------------------*/
void solserv_result_resid(FIELD *actfield,INTRA *actintra,DIST_VECTOR *sol,
                          INT place,SPARSE_ARRAY *sysarray,SPARSE_TYP *sysarray_typ);
/*----------------------------------------------------------------------*
 |  Put the results of a DIST_VECTOR to the nodes in a       genk 01/03 |
 |  certain place in ARRAY sol_mf                                       |
 |  Result have to bee allreduced and are put to the whole              |
 |  field on each proc                                                  |
 |  Functionality is the same as in solserv_result_total                |
 *----------------------------------------------------------------------*/
void solserv_result_mf(FIELD *actfield,INTRA *actintra,DIST_VECTOR *sol,
                          INT place,SPARSE_ARRAY *sysarray,
                          SPARSE_TYP *sysarray_typ);
/*----------------------------------------------------------------------*
 |  solver_service3.c                                    m.gee 04/03    |
 *----------------------------------------------------------------------*/
void solserv_sol_zero(FIELD *actfield, INT disnum, INT arraynum, INT place);
void solserv_sol_copy(FIELD *actfield, INT disnum, INT arrayfrom, INT arrayto, 
                      INT from, INT to);
void solserv_sol_add(FIELD *actfield, INT disnum, INT arrayfrom, INT arrayto, 
                      INT from, INT to, DOUBLE fac);
void solserv_sol_localassemble(INTRA *actintra, ELEMENT *actele, DOUBLE *localvec, INT arraynum,
                              INT place);
void solserv_putdirich_to_dof(FIELD *actfield, INT disnum, INT arraynum, DOUBLE scale, 
                              INT place);
void solserv_adddirich(FIELD *actfield, INT disnum, INT arraynum,
                              INT from1,INT from2,INT to,
                              DOUBLE facfrom1, DOUBLE facfrom2);
void solserv_assdirich_fac(FIELD *actfield, INT disnum, INT arraynum,
                           INT from1,INT from2,INT to, 
                           DOUBLE facfrom1, DOUBLE facfrom2);
void solserv_cpdirich(FIELD *actfield, INT disnum, INT arraynum,
                      INT from,INT to);
void solserv_zerodirich(FIELD *actfield, INT disnum, INT arraynum, INT place);
/*----------------------------------------------------------------------*
 |  solver_superlu.c                                     m.gee 11/01    |
 *----------------------------------------------------------------------*/
void solver_psuperlu_ucchb( 
                      struct _SOLVAR         *actsolv,
                      struct _INTRA          *actintra,
                      struct _UCCHB          *ucchb,
                      struct _DIST_VECTOR    *sol,
                      struct _DIST_VECTOR    *rhs,
                      INT                     option
                     );

		     
/*----------------------------------------------------------------------*
 |  input_sol.c                                          m.gee 11/01    |
 *----------------------------------------------------------------------*/
void inpctrsol(SOLVAR *solv);

/*----------------------------------------------------------------------*
 |  restart_control.c                                    m.gee 02/02    |
 *----------------------------------------------------------------------*/
void restart_write_nlnstructdyn(STRUCT_DYNAMIC  *sdyn,                  
                                STRUCT_DYN_CALC *dynvar,
                                FIELD           *actfield,
                                PARTITION       *actpart,
                                INTRA           *actintra,
                                CALC_ACTION     *action,
                                INT nrhs,  DIST_VECTOR *rhs,
                                INT nsol,  DIST_VECTOR *sol,
                                INT ndis,  DIST_VECTOR *dispi,
                                INT nvel,  DIST_VECTOR *vel,
                                INT nacc,  DIST_VECTOR *acc,
                                INT nfie,  DIST_VECTOR *fie,
                                INT nwork, DIST_VECTOR *work,
                                ARRAY *intforce_a,
                                ARRAY *dirich_a,
                                CONTAINER    *container);
void restart_read_nlnstructdyn(INT restart,
                               STRUCT_DYNAMIC  *sdyn,
                               STRUCT_DYN_CALC *dynvar,
                               FIELD           *actfield,
                               PARTITION       *actpart,
                               INTRA           *actintra,
                               CALC_ACTION     *action,
                               INT nrhs,  DIST_VECTOR *rhs,
                               INT nsol,  DIST_VECTOR *sol,
                               INT ndis,  DIST_VECTOR *dispi,
                               INT nvel,  DIST_VECTOR *vel,
                               INT nacc,  DIST_VECTOR *acc,
                               INT nfie,  DIST_VECTOR *fie,
                               INT nwork, DIST_VECTOR *work,
                               ARRAY *intforce_a,
                               ARRAY *dirich_a,
                               CONTAINER    *container);
void restart_write_nlnstructstat(STATIC_VAR     *statvar,               
                STANLN          *nln_data,  
                FIELD           *actfield,  
                PARTITION       *actpart,   
                INTRA           *actintra,  
                CALC_ACTION     *action,    
                INT kstep,                  
                INT nrhs,  DIST_VECTOR *rhs,
                INT nsol,  DIST_VECTOR *sol,
                INT ndis,  DIST_VECTOR *dispi,
                CONTAINER    *container) ;	
void restart_read_nlnstructstat(INT restart,   
                STATIC_VAR     *statvar,                         
                STANLN          *nln_data,     
                FIELD           *actfield,     
                PARTITION       *actpart,      
                INTRA           *actintra,     
                CALC_ACTION     *action,       
                INT nrhs,  DIST_VECTOR *rhs,   
                INT nsol,  DIST_VECTOR *sol,   
                INT ndis,  DIST_VECTOR *dispi, 
                CONTAINER    *container);
void restart_write_fluiddyn(FLUID_DYNAMIC   *fdyn,                  
                            FIELD	    *actfield,
                            PARTITION	    *actpart,
                            INTRA	    *actintra,
			    CALC_ACTION     *action,
			    CONTAINER       *container);
void restart_read_fluiddyn(INT restart,
                           FLUID_DYNAMIC   *fdyn,
                           FIELD	   *actfield,
                           PARTITION	   *actpart,
                           INTRA	   *actintra,
			   CALC_ACTION     *action,
			   CONTAINER       *container);
void restart_write_aledyn(ALE_DYNAMIC       *adyn,                  
                            FIELD	    *actfield,
                            PARTITION	    *actpart,
                            INTRA	    *actintra);
void restart_read_aledyn(INT restart,
                           ALE_DYNAMIC   *adyn,
                           FIELD	   *actfield,
                           PARTITION	   *actpart,
                           INTRA	   *actintra);
void restart_write_fsidyn(FSI_DYNAMIC       *fsidyn);
void restart_read_fsidyn(INT restart,
                         FSI_DYNAMIC   *fsidyn);
/*---------------------------------------------------------------------*
 | routine to find the maximum value of a distributed vector           |
 | ab =  0 absolut maximum value                                       |
 | ab =  1 maxium value                                                |
 | ab = -1 minimum value                                               |
 |                                                         genk 03/02  |
 *---------------------------------------------------------------------*/ 
void solserv_dmax_distvec(
			  DIST_VECTOR  *distvec,
			  DOUBLE *res,   /* result */
			  INT ab        /* flag */
			  );
/*----------------------------------------------------------------------*
 | ale_dyn_control.c                                          mn 06/02  |
 *----------------------------------------------------------------------*/
void dyn_ale(void);
/*----------------------------------------------------------------------*
 | ale_rhs.c                                                  mn 06/02  |
 *----------------------------------------------------------------------*/
void ale_rhs(FIELD *actfield, SOLVAR *actsolv, PARTITION *actpart,
             INTRA *actintra, INT sysarray1, INT sysarray2, DOUBLE *dirich,
             INT global_numeq, INT kstep, CONTAINER *container, CALC_ACTION *action);
/*----------------------------------------------------------------------*
 | mlpcg prototypes                                          m.gee9/02  |
 *----------------------------------------------------------------------*/
#include "prototypes_mlpcg.h"
void solver_mlpcg(struct _SOLVAR         *actsolv,
                  struct _INTRA          *actintra,
                  struct _DBCSR          *bdcsr,
                  struct _DIST_VECTOR    *sol,
                  struct _DIST_VECTOR    *rhs,
                  INT                     option);
void mlpcg_solver_create(DBCSR       *bdcsr, 
                         DIST_VECTOR *sol,
                         DIST_VECTOR *rhs,
                         MLPCGVARS   *mlpcgvars);
void mlpcg_pcg(DBCSR       *bdcsr, 
               DIST_VECTOR *sol,
               DIST_VECTOR *rhs,
               INTRA       *actintra);
void mlpcg_solver_init(DBCSR       *bdcsr, 
                       DIST_VECTOR *sol,
                       DIST_VECTOR *rhs,
                       INTRA       *actintra);
void mask_bdcsr(FIELD         *actfield, 
                PARTITION     *actpart, 
                SOLVAR        *actsolv,
                INTRA         *actintra, 
                DBCSR         *bdcsr);
void bdcsr_make_csr(FIELD         *actfield, 
                    PARTITION     *actpart, 
                    SOLVAR        *actsolv,
                    DBCSR         *bdcsr,
                    INT          **dof_connect);
void  bdcsr_nnz_topology(FIELD         *actfield, 
                         PARTITION     *actpart, 
                         SOLVAR        *actsolv,
                         INTRA         *actintra,
                         DBCSR         *bdcsr,
                         INT          **dof_connect);
void bdcsr_update(FIELD          *actfield, 
                  PARTITION      *actpart, 
                  SOLVAR         *actsolv,
                  INTRA          *actintra,
                  DBCSR          *bdcsr);
void bdcsr_numeq(FIELD         *actfield, 
                 PARTITION     *actpart, 
                 SOLVAR        *actsolv,
                 INTRA         *actintra,
                 INT           *numeq);
void  add_bdcsr(struct _PARTITION     *actpart,
                struct _SOLVAR        *actsolv,
                struct _INTRA         *actintra,
                struct _ELEMENT       *actele,
                struct _DBCSR         *bdcsr1,
                struct _DBCSR         *bdcsr2);
/*----------------------------------------------------------------------*
 | solver_copy_matrices                                      genk 11/02 |
 *----------------------------------------------------------------------*/
void solver_copy_csr(
                       SPARSE_ARRAY  *amatrix,
                       SPARSE_TYP    *amatrix_typ,
                       DBCSR         *amatrix_csr,
		       INT            numeq_total
                     );

/*----------------------------------------------------------------------*
 |  s8_contact1.c                                         m.gee 3/03    |
 *----------------------------------------------------------------------*/
#ifdef S8CONTACT
void s8_contact_detection(FIELD        *actfield, 
                          INTRA        *actintra, 
                          SPARSE_ARRAY *matrix, 
                          SPARSE_TYP   *matrix_type, 
                          DOUBLE       *cforce,
                          INT          *iscontact,
                          DOUBLE       *maxdt);
#endif

/************************************************************************
 | fluid_service.c                                                      |
 ************************************************************************/
void fluid_result_incre(  
                          FIELD             *actfield,    
                          INTRA             *actintra,   
			  DIST_VECTOR       *sol,        
                          INT                place,      
			  SPARSE_ARRAY      *sysarray,      
			  SPARSE_TYP        *sysarray_typ,
			  DOUBLE            *vrat,        
			  DOUBLE            *prat,
                          DOUBLE            *grat          
		       );		        

/************************************************************************
 | fluid_liftdrag.c                                                     |
 ************************************************************************/
void fluid_liftdrag(
    INT            init,
    CALC_ACTION   *action,
    CONTAINER     *container,
    FIELD         *actfield,
    SOLVAR        *actsolv,
    PARTITION     *actpart,
    INTRA         *actintra);

/************************************************************************
 | fluid_service_tu.c                                                   |
 ************************************************************************/
void fluid_result_incre_tu(FIELD         *actfield,    
                           INTRA         *actintra,   
			         DIST_VECTOR   *sol,        
                           INT            place,      
			         SPARSE_ARRAY  *sysarray,      
			         SPARSE_TYP    *sysarray_typ,
			         DOUBLE        *kapepsrat,        
		               FLUID_DYNAMIC *fdyn,
                           DOUBLE         lower_limit_kappa,
                           DOUBLE         lower_limit_eps          
		              );

void fluid_eddy_update(FIELD         *actfield, 
                       DIST_VECTOR   *sol   
                       );

void fluid_lenght_update(FIELD         *actfield, 
                         DIST_VECTOR   *sol,   
		             DOUBLE        *lenghtrat
                        );
/************************************************************************
 | fluid_service_tu_1.c                                                 |
 ************************************************************************/
void fluid_result_incre_tu_1( FIELD       *actfield,    
                              INTRA         *actintra,   
                              DIST_VECTOR   *sol,        
                              INT            place,      
                              SPARSE_ARRAY  *sysarray,      
                              SPARSE_TYP    *sysarray_typ,
                              DOUBLE        *kapomegarat,
                              DOUBLE         lower_limit_kappa,
                              DOUBLE         lower_limit_omega         
		              );

void fluid_eddy_update_1(FIELD         *actfield, 
                         DIST_VECTOR   *sol   
                        );

void fluid_lenght_update_1(FIELD         *actfield, 
                          DIST_VECTOR   *sol,   
		              DOUBLE        *lenghtrat
                         );
/* -------------------------------------------------------------------- *
 *   global_oll_add.c                                          mn 05/03 *
 * -------------------------------------------------------------------- */
void add_oll_sendbuff(
    INT                   ii,
    INT                   jj,
    INT                   i,
    INT                   j,
    INT                   ii_owner,
    INT                 **isend,
    DOUBLE              **dsend,
    DOUBLE              **estif,
    INT                   numsend);
void add_oll_checkcouple(
    INT                   ii,
    INT                 **cdofs,
    INT                   ncdofs,
    INT                  *iscouple,
    INT                  *isowner,
    INT                   nprocs);
void  add_oll(
    struct _PARTITION     *actpart,
    struct _INTRA         *actintra,
    struct _ELEMENT       *actele,
    struct _OLL           *oll1,
    struct _OLL           *oll2);
void exchange_coup_oll(
    PARTITION             *actpart,
    INTRA                 *actintra,
    struct _OLL           *oll);
/* -------------------------------------------------------------------- *
 *   solver_oll.c                                              mn 05/03 *
 * -------------------------------------------------------------------- */
void solver_oll(
    struct _SOLVAR         *actsolv,
    struct _INTRA          *actintra,
    struct _OLL            *oll,
    struct _DIST_VECTOR    *sol,
    struct _DIST_VECTOR    *rhs,
    INT                     option);
/* -------------------------------------------------------------------- *
 *   solver_oll_aztec.c                                        mn 05/03 *
 * -------------------------------------------------------------------- */
void solver_az_oll( 
    struct _SOLVAR         *actsolv,
    struct _INTRA          *actintra,
    struct _AZ_ARRAY_MSR   *msr_array,
    struct _DIST_VECTOR    *sol,
    struct _DIST_VECTOR    *rhs,
    INT                     option);
/* -------------------------------------------------------------------- *
 *   solver_oll_spooles.c                                      mn 05/03 *
 * -------------------------------------------------------------------- */
void solver_spo_oll( 
    struct _SOLVAR         *actsolv,
    struct _INTRA          *actintra,
    struct _OLL            *oll,
    struct _SPOOLMAT       *spo,
    struct _DIST_VECTOR    *sol,
    struct _DIST_VECTOR    *rhs,
    INT                     option);
/* -------------------------------------------------------------------- *
 *   global_oll_service.c                                      mn 05/03 *
 * -------------------------------------------------------------------- */
void oll_numeq(
    struct _FIELD          *actfield, 
    struct _PARTITION      *actpart, 
    struct _INTRA          *actintra,
    INT                     dis,
    INT                    *numeq);
void oll_dof_in_coupledofs(
    INT dof,
    PARTITION *actpart,
    INT *iscoupled,
    INT  dis);
void oll_dof_find_centernode(
    INT dof,
    PARTITION *actpart,
    NODE **centernode,
    INT dis);
void oll_nnz_topology(
    struct _FIELD         *actfield, 
    struct _PARTITION     *actpart, 
    struct _INTRA         *actintra,
    struct _OLL           *oll,
    INT            dis);
void oll_update(
    struct _FIELD          *actfield, 
    struct _PARTITION      *actpart, 
    struct _INTRA          *actintra,
    INT                     dis,
    struct _OLL            *oll);
INT oll_getindex(
    INT                     dof,
    INT                    *update,
    INT                     length);
void oll_print(
    struct _OLL            *oll,
    INT                     n_max);
void oll_pattern(
    struct _OLL            *oll,
    INT                     n_max);
void oll_gnupattern(
    struct _OLL            *oll
    );

/* -------------------------------------------------------------------- *
 *   global_oll_service2.c                                     mn 05/03 *
 * -------------------------------------------------------------------- */
void oll_open(struct _OLL       *matrix,
              INT        numeq,
              INT        numeq_total,
              struct _FIELD     *actfield,
              struct _PARTITION *actpart,
              struct _INTRA     *actintra,
              INT        dis);
/*void oll_getentry(
    struct _OLL            *matrix,
    struct _MATENTRY      **ret);*/
void oll_zero(
    struct _OLL            *oll);
void oll_add(
    struct _OLL            *oll1,
    struct _OLL            *oll2,
    DOUBLE                  factor);
void oll_scal(
    struct _OLL            *oll,
    DOUBLE                  factor);
void oll_cp_mask(
    struct _OLL            *from,
    struct _OLL            *to);
void oll_delete(
    struct _OLL            *matrix);
void oll_setval(
    struct _OLL            *matrix,
    INT                     actrow,
    INT                     actcol,
    DOUBLE val);
void oll_addval(
    struct _OLL            *matrix,
    INT                     actrow,
    INT                     actcol,
    DOUBLE                  val);
void oll_addrow(
    struct _OLL            *matrix, 
    INT                     actrow,
    INT                     lm[],
    DOUBLE                  val[],
    INT                     nd);
void oll_to_sky(
    struct _OLL            *oll,
    union  _SPARSE_ARRAY   *sysarray);
void oll_to_spo(
    struct _OLL            *oll,
    union  _SPARSE_ARRAY   *sysarray);
void oll_to_msr(
    struct _OLL            *oll,
    union  _SPARSE_ARRAY   *sysarray);
void oll_to_ccf(
                struct _OLL     *oll,
                union  _SPARSE_ARRAY    *sysarray);
void oll_copy(
    struct _OLL            *from,
    struct _OLL            *to);

