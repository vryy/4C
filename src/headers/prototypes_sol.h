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
            NR_CONTROLTYP  controltyp,    /* type of control algorithm */
            CONTAINER     *container     /*!< contains variables defined in container.h */
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
            NR_CONTROLTYP  controltyp,     /* type of control algorithm */
            CONTAINER     *container      /*!< contains variables defined in container.h */
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
 | cal_nlndyn_stru_expl.c                                m.gee 05/02    |
 *----------------------------------------------------------------------*/
void dyn_nln_stru_expl(void); 
/*----------------------------------------------------------------------*
 | dyn_service.c                                         m.gee 02/02    |
 *----------------------------------------------------------------------*/
void dyn_eigen(FIELD *actfield, PARTITION *actpart, SOLVAR *actsolv,
               INTRA *actintra, int stiff, int mass);
void kefnln_struct(STRUCT_DYN_CALC *dynvar, 
                  STRUCT_DYNAMIC  *sdyn,
                  FIELD           *actfield,
                  SOLVAR          *actsolv,
                  INTRA           *actintra,
                  DIST_VECTOR     *work,
                  int              stiff_array,
                  int              mass_array,
                  int              damp_array);
void dyn_keff_expl(INTRA *actintra,
                   SPARSE_TYP *sysarray_typ, SPARSE_ARRAY *sysarray,
                   int stiff_array, int mass_array, int damp_array,
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
void dyn_setconstants_expl(STRUCT_DYN_CALC *dynvar, STRUCT_DYNAMIC *sdyn, double dt);
void dyn_nlnstructupd(FIELD *actfield,      STRUCT_DYN_CALC *dynvar, 
                      STRUCT_DYNAMIC *sdyn, SOLVAR *actsolv,
                      DIST_VECTOR *sol_old, DIST_VECTOR *sol_new,
                      DIST_VECTOR *rhs_new, DIST_VECTOR *rhs_old,
                      DIST_VECTOR *vel,     DIST_VECTOR *acc,
                      DIST_VECTOR *work0,   DIST_VECTOR *work1,
                      DIST_VECTOR *work2);
void dyn_nlnstruct_outhead(STRUCT_DYN_CALC *dynvar, STRUCT_DYNAMIC *sdyn);
void dyn_nlnstru_outhead_expl(void);
void dyn_nlnstruct_outstep(STRUCT_DYN_CALC *dynvar, STRUCT_DYNAMIC *sdyn, int numiter, double dt);
void assemble_dirich_dyn(ELEMENT *actele, ARRAY *estif_global, 
                         ARRAY *emass_global, CONTAINER *container);
void dyn_epot(FIELD *actfield, int disnum, INTRA *actintra, STRUCT_DYN_CALC *dynvar, double *deltaepot);
void dyn_ekin(FIELD *actfield, SOLVAR *actsolv, PARTITION *actpart, INTRA *actintra, CALC_ACTION *action,
             CONTAINER *container, int stiff_array, int mass_array);
void dyn_eout(STRUCT_DYN_CALC *dynvar, STRUCT_DYNAMIC *sdyn, INTRA *actintra, SOLVAR *actsolv, 
              DIST_VECTOR *dispi, /* converged incremental displacements */
              DIST_VECTOR *rhs1,  /* load at time t                      */ 
              DIST_VECTOR *rhs2,  /* load at time t-dt                   */ 
              DIST_VECTOR *work0);
void dyn_ekin_local(ELEMENT *actele,ARRAY *emass, CONTAINER  *container);

/*----------------------------------------------------------------------*
 | fluid_service.c                                         genk 04/02   |
 *----------------------------------------------------------------------*/
/*!---------------------------------------------------------------------                                         
\brief storing results in solution history

<pre>                                                         genk 05/02

in this routine the results in the DIST_VECTOR are put to the nodes in
a certain place in ARRAY sol_increment.
Result has to be allreduced and is put to the whole field on each proc.
If necassary the norms for the iteration check of the nonlinear iteration
scheme are calculated	   
			     
</pre>   
\param **actfield      FIELD	      (i)    actual field       
\param  *actintra      INTRA	      (i)    actual intra comm.
\param	*sol 	       DIST_VECTOR    (i)    solution vector
\param	 place         int	      (i)    place in sol_incr.
\param	*sysarray      SPARSE_ARRAY   (i)
\param	*sysarray_typ  SPARSE_TYP     (i)
\param	*vrat          double	      (o)    vel.  conv. ratio
\param	*prat          double	      (o)    pre.  conv. ratio
\param	*fdyn	       FLUID_DYNAMIC	     	
\return void 

------------------------------------------------------------------------*/
void fluid_result_incre(FIELD         *actfield,    
                        INTRA         *actintra,   
			DIST_VECTOR   *sol,        
                        int            place,      
			SPARSE_ARRAY  *sysarray,      
			SPARSE_TYP    *sysarray_typ,
			double        *vrat,        
			double        *prat,       
			FLUID_DYNAMIC *fdyn           
		       ); 
/*----------------------------------------------------------------------*
 | global_calelm.c                                       m.gee 11/01    |
 *----------------------------------------------------------------------*/
void calelm(FIELD        *actfield,     /* active field */        
            SOLVAR       *actsolv,      /* active SOLVAR */
            PARTITION    *actpart,      /* my partition of this field */
            INTRA        *actintra,     /* my intra-communicator */
            int           sysarray1,    /* number of first sparse system matrix */
            int           sysarray2,    /* number of secnd system matrix, if present, else -1 */
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
            int           actsysarray,  /* the active sparse array */
            DIST_VECTOR  *rhs1,         /* 2 dist. vectors for rhs */
            CALC_ACTION  *action,       /* action to be passed to element routines */
            CONTAINER    *container);    /*!< contains variables defined in container.h */
void rhs_point_neum(double *rhs, int dimrhs, PARTITION *actpart);     
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
	      int            actdis);
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
                          int         **dof_connect,
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
                       int         **dof_connect);
void  ccf_make_bindx(FIELD         *actfield, 
                     PARTITION     *actpart, 
                     SOLVAR        *actsolv,
                     CCF           *ccf,
                     int           *bindx,
                     ARRAY         *red_dof_connect);
void  ccf_make_sparsity(CCF        *ccf,
                        int        *bindx);
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
 |  global_mask_spooles.c                                m.gee 05/02    |
 *----------------------------------------------------------------------*/
void mask_spooles(FIELD         *actfield, 
                  PARTITION     *actpart, 
                  SOLVAR        *actsolv,
                  INTRA         *actintra, 
                  SPOOLMAT      *spo);
void     spo_make_sparsity(SPOOLMAT        *spo,
                           int           *bindx);
void    spo_make_bindx(FIELD         *actfield, 
                       PARTITION     *actpart, 
                       SOLVAR        *actsolv,
                       SPOOLMAT      *spo,
                       int          **dof_connect,
                       int           *bindx);
void  spo_nnz_topology(FIELD         *actfield, 
                       PARTITION    *actpart, 
                       SOLVAR       *actsolv,
                       INTRA        *actintra,
                       SPOOLMAT     *spo,
                       int         **dof_connect);
void spo_update(FIELD         *actfield, 
                PARTITION     *actpart, 
                SOLVAR        *actsolv,
                INTRA         *actintra,
                SPOOLMAT      *spo);
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
                 enum _ASSEMBLE_ACTION  assemble_action, /* the assembly option */
                 CONTAINER             *container);  /*!< contains variables defined in container.h */                 
void init_assembly(
                       struct _PARTITION      *actpart,
                       struct _SOLVAR         *actsolv,
                       struct _INTRA          *actintra,
                       struct _FIELD          *actfield,
                       int                     actsysarray,
		       int                     actndis
                     );
void assemble_vec(INTRA        *actintra,
                    SPARSE_TYP   *sysarraytyp,
                    SPARSE_ARRAY *sysarray,
                    DIST_VECTOR  *rhs,
                    double       *drhs,
                    double        factor);
void assemble_intforce(ELEMENT *actele, ARRAY *elevec_a, CONTAINER *container);
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
                struct _RC_PTR        *rc_ptr1,
                struct _RC_PTR        *rc_ptr2);
void add_rcptr_sendbuff(int ii,int jj,int i,int j,int ii_owner,int **isend,
                    double **dsend,double **estif, int numsend);
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
void add_spo_sendbuff(int ii,int jj,int i,int j,int ii_owner,int **isend,
                      double **dsend,double **estif, int numsend);
void exchange_coup_spo(
                         PARTITION     *actpart,
                         SOLVAR        *actsolv,
                         INTRA         *actintra,
                         SPOOLMAT      *spo
                        );
void add_val_spo(int ii,int jj, struct _SPOOLMAT *spo, double val, INTRA *actintra);
void set_val_spo(int ii,int jj, struct _SPOOLMAT *spo, double val, INTRA *actintra);
void close_spooles_matrix(struct _SPOOLMAT *spo, INTRA *actintra);
void add_spooles_matrix(struct _SPOOLMAT *to, struct _SPOOLMAT *from,
                       double factor, int init, INTRA *actintra);
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
 |  solver_umfpack.c                              s.offermanns 02/02    |
 *----------------------------------------------------------------------*/
void solver_umfpack(struct _SOLVAR         *actsolv,
                    struct _INTRA          *actintra,
                    struct _CCF            *ccf,
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
 |  solver_spooles.c                                     m.gee 05/02    |
 *----------------------------------------------------------------------*/
void solver_spooles( 
                      struct _SOLVAR         *actsolv,
                      struct _INTRA          *actintra,
                      struct _SPOOLMAT       *spo,
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
/*----------------------------------------------------------------------*
 |  make A = A * factor                                      m.gee 02/02|
 | SPARSE_TYP *Atyp (i)        type of sparse matrix                    |
 | SPARSE_ARRAY *A  (i/o)      sparse matrix                            |
 | double factor    (i)        factor                                   |
 *----------------------------------------------------------------------*/
void solserv_scal_mat(SPARSE_TYP *Atyp,SPARSE_ARRAY *A,double factor);
/*----------------------------------------------------------------------*
 |  make A = A + B * factor                                  m.gee 02/02|
 | INTRA *actintra  (i)     intra-communicator the matrices live on     |
 | SPARSE_TYP *Atyp (i)     type of sparse matrix                       |
 | SPARSE_ARRAY *A  (i/o)   sparse matrix                               |
 | SPARSE_TYP *Btyp (i)     type of sparse matrix                       |
 | SPARSE_ARRAY *B  (i)     sparse matrix                               |
 | double factor    (i)     factor                                      |
 *----------------------------------------------------------------------*/
void solserv_add_mat(INTRA *actintra,
                     SPARSE_TYP *Atyp,SPARSE_ARRAY *A,
                     SPARSE_TYP *Btyp, SPARSE_ARRAY *B,
                     double factor);
/*----------------------------------------------------------------------*
 |  get dimensions of sparse matrix                          m.gee 02/02|
 | SPARSE_ARRAY mat  (i)     sparse matrix                              |
 | SPARSE_TYP mattyp (i)     type of sparse matrix                      |
 | int *numeq        (o)     proc-local dimension of sparse matrix      |
 | int *numeq_total  (o)     global dimension of sparse matrix          |
 *----------------------------------------------------------------------*/
void solserv_getmatdims(SPARSE_ARRAY mat,SPARSE_TYP mattyp,
                        int *numeq, int *numeq_total);
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
 |  int numvectors       (i)   number of DIST_VECTORs to allocate       |
 |  int numeq_total      (i)   proc-global dimension of the DIST_VECTORs|
 |  int numeq            (i)   proc_local  dimension of the DIST_VECTORs|
 |  char typstr[]        (i)   ="DV" for double-DIST_VECTORs            |
 |  the values in the DIST_VECTORs is NOT initialized                   |
 *----------------------------------------------------------------------*/
void solserv_create_vec(DIST_VECTOR **vector,int numvectors,
                        int numeq_total,
                        int numeq,        char typstr[]);
/*----------------------------------------------------------------------*
 |   delete number of distributed vectors - collective call ! m.gee 2/02|
 |  DIST_VECTOR **vector (i/o) adress of pointer a vector of            |
 |                             DIST_VECTORs is allocated to             |
 |  int numvectors       (i)   number of DIST_VECTORs to free           |
 |  the routine frees all DIST_VECTORs in vector and sets vector=NULL   |
 *----------------------------------------------------------------------*/
void solserv_del_vec(DIST_VECTOR **vector,int numvectors);
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
 |  double factor         (i)   scaling factor                          |
 *----------------------------------------------------------------------*/
void solserv_add_vec(DIST_VECTOR *vec_from,DIST_VECTOR *vec_to,
                     double factor);
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
 |  double *result        (o) norm of the vector                        |
 *----------------------------------------------------------------------*/
void solserv_vecnorm_euclid(INTRA *actintra,DIST_VECTOR *dist_vec,
                            double *result);
/*----------------------------------------------------------------------*
 |  find  absolute maximum value in a vector (Linf-Norm)     m.gee 02/02|
 |  *result = MAX( ABS(vec[i]) )                                        |
 |  INTRA *actintra (i) intra-communicator the DIST_VECTOR lives on     |
 |  DIST_VECTOR *dist_vec (i) vector to make norm of                    |
 |  double *result        (o) norm of the vector                        |
 *----------------------------------------------------------------------*/
void solserv_vecnorm_Linf(INTRA *actintra,DIST_VECTOR *dist_vec,
                          double *result);
/*----------------------------------------------------------------------*
 |  get a certain entry from a distr. vector to all procs    m.gee 11/01|
 |  returns the value of dof indiz in the vector dist_vec on all procs  |
 |  INTRA *actintra (i) intra-communicator the DIST_VECTOR lives on     |
 |  SPARSE_TYP *sysarray_typ (i) sparsity typ of vector-matching matrix |
 |  SPARSE_ARRAY *sysarray   (i) sparse matrix the vector matches in    |
 |                               distribution                           | 
 |  DIST_VECTOR  *dist_vec   (i) vector the value shall be taken from   |
 |  int           indiz      (i) field-local (unsupported) dof number   |
 |  double       *result     (o) value in vector at the given dof       |
 |                               returned redundant on all procs        |
 *----------------------------------------------------------------------*/
void solserv_getele_vec(INTRA*actintra,SPARSE_TYP *sysarray_typ,
                        SPARSE_ARRAY *sysarray,DIST_VECTOR *dist_vec,
                        int indiz,double *result);
/*----------------------------------------------------------------------*
 |  make dot product between 2 distr. vectors                m.gee 11/01|
 |  *dot = sumover_i( vec1[i]*vec2[i] )                                 |
 |  INTRA *actintra (i) intra-communicator the DIST_VECTORs live on     |
 |  DIST_VECTOR *dist_vec1 (i) first vector to be multiplied            |
 |  DIST_VECTOR *dist_vec2 (i) scnd  vector to be multiplied            |
 |  double      *dot       (o) result of vector-vector multiplication   |
 |                             returned redundant on all procs          |
 *----------------------------------------------------------------------*/
void solserv_dot_vec(INTRA *actintra,DIST_VECTOR *dist_vec1,
                     DIST_VECTOR *dist_vec2,double *dot);
/*----------------------------------------------------------------------*
 |  make product between scalar and distr. vector            m.gee 11/01|
 |  vec[i] = vec[i] * scalar                                            |
 |  DIST_VECTOR *dist_vec (i/o) vector to be multiplied by scalar       |
 |  double       scalar   (o)   scalar value                            |
 *----------------------------------------------------------------------*/
void solserv_scalarprod_vec(DIST_VECTOR *dist_vec,double scalar);
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
 |  double *fullvec (o) vector of lenght numeq_total will be holding    |
 |                      the values from distvec in correct dof-ordering:|
 |                      fullvec[dof] = value of a certain dof           |
 |  INTRA *actintra (i) intra-communicator the DIST_VECTOR lives on     |
 *----------------------------------------------------------------------*/
void solserv_reddistvec(DIST_VECTOR *distvec,SPARSE_ARRAY *sysarray,
                        SPARSE_TYP *sysarray_typ,double *fullvec,
                        int dim,INTRA *actintra);
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
 |  double *fullvec (o) vector of lenght numeq_total  holding           |
 |                      the values  correct dof-ordering:               |
 |                      fullvec[dof] = value of a certain dof           |
 |  INTRA *actintra (i) intra-communicator the DIST_VECTOR lives on     |
 *----------------------------------------------------------------------*/
void solserv_distribdistvec(DIST_VECTOR  *distvec,SPARSE_ARRAY *sysarray,
                            SPARSE_TYP *sysarray_typ,double *fullvec,
                            int dim,INTRA *actintra);
/*----------------------------------------------------------------------*
 |  Put the results of a DIST_VECTOR to the nodes in a       m.gee 10/01|
 |  certain place  in ARRAY sol                                         |
 |  Result has to be allreduced and are put to the whole                |
 |  field on each procs                                                 |
 |  FIELD *actfield (i) the active field                                |
 |  INTRA *actintra (i) intra-communicator the DIST_VECTOR lives on     |
 |  DIST_VECTOR *sol (i) vector of values to be put to the nodes        |
 |  int place        (i) place in the ARRAY node->sol where to put the  |
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
                          int place,SPARSE_ARRAY *sysarray,
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
                          SPARSE_TYP *sysarray_typ);
/*----------------------------------------------------------------------*
 |  Put the results of a DIST_VECTOR to the nodes in a       m.gee 11/01|
 |  certain place in ARRAY sol_residual                                 |
 |  Result have to bee allreduced and are put to the whole              |
 |  field on each proc                                                  |
 |  Functionality is the same as in solserv_result_total                |
 *----------------------------------------------------------------------*/
void solserv_result_resid(FIELD *actfield,INTRA *actintra,DIST_VECTOR *sol,
                          int place,SPARSE_ARRAY *sysarray,SPARSE_TYP *sysarray_typ);
/*----------------------------------------------------------------------*
 |  solver_service3.c                                    m.gee 04/03    |
 *----------------------------------------------------------------------*/
void solserv_sol_zero(FIELD *actfield, int disnum, int arraynum, int place);
void solserv_sol_copy(FIELD *actfield, int disnum, int arrayfrom, int arrayto, 
                      int from, int to);
void solserv_sol_add(FIELD *actfield, int disnum, int arrayfrom, int arrayto, 
                      int from, int to, double fac);
void solserv_sol_localassemble(INTRA *actintra, ELEMENT *actele, double *localvec, int arraynum,
                              int place);
void solserv_putdirich_to_dof(FIELD *actfield, int disnum, int arraynum, double scale, 
                              int place);
void solserv_adddirich(FIELD *actfield, int disnum, int arraynum,
                              int from1,int from2,int to,
                              double facfrom1, double facfrom2);
void solserv_assdirich_fac(FIELD *actfield, int disnum, int arraynum,
                           int from1,int from2,int to, 
                           double facfrom1, double facfrom2);
void solserv_cpdirich(FIELD *actfield, int disnum, int arraynum,
                      int from,int to);
void solserv_zerodirich(FIELD *actfield, int disnum, int arraynum, int place);
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
 |  input_sol.c                                          m.gee 11/01    |
 *----------------------------------------------------------------------*/
void inpctrsol(SOLVAR *solv);

/*----------------------------------------------------------------------*
 |  restart_control.c                                    m.gee 05/02    |
 *----------------------------------------------------------------------*/
void restart_write_nlnstructdyn(STRUCT_DYNAMIC  *sdyn,
                                STRUCT_DYN_CALC *dynvar,
                                FIELD           *actfield,
                                PARTITION       *actpart,
                                INTRA           *actintra,
                                CALC_ACTION     *action,
                                int nrhs,  DIST_VECTOR *rhs,
                                int nsol,  DIST_VECTOR *sol,
                                int ndis,  DIST_VECTOR *dispi,
                                int nvel,  DIST_VECTOR *vel,
                                int nacc,  DIST_VECTOR *acc,
                                int nfie,  DIST_VECTOR *fie,
                                int nwork, DIST_VECTOR *work,
                                ARRAY *intforce_a,
                                ARRAY *dirich_a,
                                CONTAINER    *container);     /*!< contains variables defined in container.h */
void restart_read_nlnstructdyn(int restart,
                               STRUCT_DYNAMIC  *sdyn,
                               STRUCT_DYN_CALC *dynvar,
                               FIELD           *actfield,
                               PARTITION       *actpart,
                               INTRA           *actintra,
                               CALC_ACTION     *action,
                               int nrhs,  DIST_VECTOR *rhs,
                               int nsol,  DIST_VECTOR *sol,
                               int ndis,  DIST_VECTOR *dispi,
                               int nvel,  DIST_VECTOR *vel,
                               int nacc,  DIST_VECTOR *acc,
                               int nfie,  DIST_VECTOR *fie,
                               int nwork, DIST_VECTOR *work,
                               ARRAY *intforce_a,
                               ARRAY *dirich_a,
                               CONTAINER    *container);     /*!< contains variables defined in container.h */ 
void restart_write_nlnstructstat(STATIC_VAR     *statvar,/*------------ static input --*/                  
             STANLN          *nln_data,  /*-- control variables for global NR-Iterat --*/
             FIELD           *actfield,  /*---------------------------- actual field --*/
             PARTITION       *actpart,   /*------------------------ actual partition --*/
             INTRA           *actintra,  /*---------------- actual intra comunicator --*/
             CALC_ACTION     *action,    /*---------- element action = write-restart --*/
             int kstep,                  /*------------------------ actual load step --*/
             int nrhs,  DIST_VECTOR *rhs,/*-- Fext processorpart of actual load step --*/
             int nsol,  DIST_VECTOR *sol,/* solution processorpart of actual load step */
             int ndis,  DIST_VECTOR *dispi,/*- displacement processorpart  --"--     --*/
             CONTAINER    *container);     /*!< contains variables defined in container.h */
void restart_read_nlnstructstat(int restart,   /*------------------------ restart step ??? --*/
                STATIC_VAR     *statvar,       /*---------------------------- static input --*/                  
                STANLN          *nln_data,     /*-- control variables for global NR-Iterat --*/
                FIELD           *actfield,     /*---------------------------- actual field --*/
                PARTITION       *actpart,      /*------------------------ actual partition --*/
                INTRA           *actintra,     /*---------------- actual intra comunicator --*/
                CALC_ACTION     *action,       /*---------- element action = write-restart --*/
                int nrhs,  DIST_VECTOR *rhs,   /*-- Fext processorpart of actual load step --*/
                int nsol,  DIST_VECTOR *sol,   /*-- solution processorpart     --"--       --*/
                int ndis,  DIST_VECTOR *dispi, /*-- displacement processorpart  --"--     --*/
                 CONTAINER    *container);       /*!< contains variables defined in container.h */
/*---------------------------------------------------------------------*
 | routine to find the maximum value of a distributed vector           |
 | ab =  0 absolut maximum value                                       |
 | ab =  1 maxium value                                                |
 | ab = -1 minimum value                                               |
 |                                                         genk 03/02  |
 *---------------------------------------------------------------------*/ 
void solserv_dmax_distvec(
			  DIST_VECTOR  *distvec,
			  double *res,   /* result */
			  int ab        /* flag */
			  );
/*----------------------------------------------------------------------*
 | ale_calelm.c                                               mn 06/02  |
 *----------------------------------------------------------------------*/
void ale_calelm(FIELD *actfield, SOLVAR *actsolv, PARTITION *actpart, 
                INTRA *actintra, int sysarray1, int sysarray2, CALC_ACTION  *action);
/*----------------------------------------------------------------------*
 | ale_rhs.c                                                  mn 06/02  |
 *----------------------------------------------------------------------*/
void ale_rhs(FIELD *actfield, SOLVAR *actsolv, PARTITION *actpart,
             INTRA *actintra, int sysarray1, int sysarray2, double *dirich,
             int global_numeq, int kstep, CONTAINER *container, CALC_ACTION *action);
/*----------------------------------------------------------------------*
 | mlpcg prototypes                                          m.gee9/02  |
 *----------------------------------------------------------------------*/
#include "prototypes_mlpcg.h"
void solver_mlpcg(struct _SOLVAR         *actsolv,
                  struct _INTRA          *actintra,
                  struct _DBCSR          *bdcsr,
                  struct _DIST_VECTOR    *sol,
                  struct _DIST_VECTOR    *rhs,
                  int                     option);
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
                    int          **dof_connect);
void  bdcsr_nnz_topology(FIELD         *actfield, 
                         PARTITION     *actpart, 
                         SOLVAR        *actsolv,
                         INTRA         *actintra,
                         DBCSR         *bdcsr,
                         int          **dof_connect);
void bdcsr_update(FIELD          *actfield, 
                  PARTITION      *actpart, 
                  SOLVAR         *actsolv,
                  INTRA          *actintra,
                  DBCSR          *bdcsr);
void bdcsr_numeq(FIELD         *actfield, 
                 PARTITION     *actpart, 
                 SOLVAR        *actsolv,
                 INTRA         *actintra,
                 int           *numeq);
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
		       int            numeq_total
                     );
/*! 
\addtogroup CONTACT 
*//*! @{ (documentation module open)*/
/*----------------------------------------------------------------------*
 |  s8_contact1.c                                         m.gee 3/03    |
 *----------------------------------------------------------------------*/
#ifdef S8CONTACT
void s8_contact_detection(FIELD        *actfield, 
                          INTRA        *actintra, 
                          SPARSE_ARRAY *matrix, 
                          SPARSE_TYP   *matrix_type, 
                          double       *cforce,
                          int          *iscontact,
                          double       *maxdt);
#endif
/*! @} (documentation module close)*/
