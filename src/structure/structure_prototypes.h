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


#ifndef STRUCTURE_PROTOTYPES_H
#define STRUCTURE_PROTOTYPES_H

/*----------------------------------------------------------------------*
  | stru_static_lin.c                                     m.gee 11/01    |
 *----------------------------------------------------------------------*/
void calsta(void);
void stalin(void);

/*----------------------------------------------------------------------*
  | stru_static_service.c                                 m.gee 11/01    |
 *----------------------------------------------------------------------*/
void calstatserv_findcontroldof(
    FIELD     *actfield,
    INT        control_node_global,
    INT        control_dof,
    NODE     **node,
    INT       *cdof); 

void calstatserv_findreldofs(
    FIELD     *actfield,
    INT       *reldisnode_ID,
    INT       *reldis_dof,
    INT        num_reldis,
    INT       *reldof); 

void get_stepsize(
    INT         kstep,
    STATIC_VAR *statvar); 


/*----------------------------------------------------------------------*
  | stru_static_nln.c                                     m.gee 11/01    |
 *----------------------------------------------------------------------*/
void stanln(void);

void conpre(
    FIELD         *actfield,     /* the actual physical field */
    SOLVAR        *actsolv,      /* the field-corresponding solver */
    PARTITION     *actpart,      /* the partition of the proc */
    INTRA         *actintra,     /* the intra-communicator of this field */
    CALC_ACTION   *action,       /* calculation flag */
    INT            kstep,        /* the load or time step we are in */
    INT            actsysarray,  /* number of the system matrix in 
                                    actsolv->sysarray[actsysarray] to be used */
    DIST_VECTOR   *rsd,          /* dist. vector of incremental residual forces
                                    used for iteration in conequ */
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
    INT            actsysarray,   /* number of the system matrix in 
                                     actsolv->sysarray[actsysarray] to be used */
    DIST_VECTOR   *rsd,           /* dist. vector of incremental residual forces used
                                     for iteration in conequ */
    DIST_VECTOR   *dispi,         /* dist. vector of incremental displacements */
    DIST_VECTOR   *re,            /* re[0..2] 3 vectors for residual displacements */
    INT            cdof,          /* number of dof to be controlled */
    INT           *reldof,        /* numbers of dofs for output */
    STANLN        *nln_data,      /* data of the Newton-Raphson method */
    NR_CONTROLTYP  controltyp,    /* type of control algorithm */
    CONTAINER     *container      /*!< contains variables defined in container.h */
    );

void conequ_printhead(
    INT            kstep,
    NR_CONTROLTYP  controltyp,
    INT            cdof,
    DOUBLE         csp);

void conequ_printiter(
    INT            itnum,
    DOUBLE         disval,
    DOUBLE         rlnew,
    DOUBLE         dinorm,
    DOUBLE         renorm,
    DOUBLE         energy,
    DOUBLE         dnorm,
    DOUBLE         rrnorm);

void increment_controlarc(
    INTRA         *actintra,
    SOLVAR        *actsolv,
    INT            actsysarray,
    DIST_VECTOR   *rsd,
    DIST_VECTOR   *dispi,
    DIST_VECTOR   *disp,
    DOUBLE         rlnew,
    DOUBLE         rlold, 
    DOUBLE         stepsize,  
    DOUBLE        *rli);

void increment_controldisp(
    INTRA         *actintra,
    SOLVAR        *actsolv,
    INT            actsysarray,
    INT            cdof,
    DIST_VECTOR   *rsd,
    DIST_VECTOR   *dispi,   
    DOUBLE        *rli);


/*----------------------------------------------------------------------*
  | stru_dyn_nln.c                                        m.gee 02/02    |
 *----------------------------------------------------------------------*/
void dyn_nln_structural(void); 


/*----------------------------------------------------------------------*
  | stru_dyn_nln_expl.c                                   m.gee 05/02    |
 *----------------------------------------------------------------------*/
void dyn_nln_stru_expl(void); 


/*----------------------------------------------------------------------*
  | stru_dyn_service.c                                    m.gee 02/02    |
 *----------------------------------------------------------------------*/
void dyn_eigen(
    FIELD           *actfield,
    PARTITION       *actpart,
    SOLVAR          *actsolv,
    INTRA           *actintra,
    INT              stiff,
    INT              mass);

void kefnln_struct(
    STRUCT_DYN_CALC *dynvar, 
    STRUCT_DYNAMIC  *sdyn,
    FIELD           *actfield,
    SOLVAR          *actsolv,
    INTRA           *actintra,
    DIST_VECTOR     *work,
    INT              stiff_array,
    INT              mass_array,
    INT              damp_array);

void dyn_keff_expl(
    INTRA           *actintra,
    SPARSE_TYP      *sysarray_typ,
    SPARSE_ARRAY    *sysarray,
    INT              stiff_array,
    INT              mass_array,
    INT              damp_array,
    STRUCT_DYN_CALC *dynvar,
    STRUCT_DYNAMIC  *sdyn);

void pefnln_struct(
    STRUCT_DYN_CALC *dynvar, 
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

void dynnle(
    STRUCT_DYN_CALC *dynvar,
    STRUCT_DYNAMIC  *sdyn,
    INTRA           *actintra,
    SOLVAR          *actsolv, 
    DIST_VECTOR     *dispi, /* converged incremental displacements */
    DIST_VECTOR     *fie1,  /* internal forces at time t-dt */
    DIST_VECTOR     *fie2,  /* internal forces at time t */
    DIST_VECTOR     *rhs1,  /* load at time t                      */ 
    DIST_VECTOR     *rhs2,  /* load at time t-dt                   */ 
    DIST_VECTOR     *work0);

void dyne(
    STRUCT_DYN_CALC *dynvar,
    INTRA           *actintra,
    SOLVAR          *actsolv,
    INT              mass_array,
    DIST_VECTOR     *vel,
    DIST_VECTOR     *work);

void dyn_setconstants(
    STRUCT_DYN_CALC *dynvar,
    STRUCT_DYNAMIC  *sdyn,
    DOUBLE           dt);

void dyn_setconstants_expl(
    STRUCT_DYN_CALC *dynvar,
    STRUCT_DYNAMIC  *sdyn,
    DOUBLE           dt);

void dyn_nlnstructupd(
    FIELD           *actfield,
    STRUCT_DYN_CALC *dynvar, 
    STRUCT_DYNAMIC  *sdyn,
    SOLVAR          *actsolv,
    DIST_VECTOR     *sol_old,
    DIST_VECTOR     *sol_new,
    DIST_VECTOR     *rhs_new,
    DIST_VECTOR     *rhs_old,
    DIST_VECTOR     *vel,
    DIST_VECTOR     *acc,
    DIST_VECTOR     *work0,
    DIST_VECTOR     *work1,
    DIST_VECTOR     *work2);

void dyn_nlnstruct_outhead(
    STRUCT_DYN_CALC *dynvar,
    STRUCT_DYNAMIC  *sdyn);

void dyn_nlnstru_outhead_expl(void);

void dyn_nlnstruct_outstep(
    STRUCT_DYN_CALC *dynvar, 
    STRUCT_DYNAMIC  *sdyn,
    INT              numiter,
    DOUBLE           dt);

void assemble_dirich_dyn(
    ELEMENT         *actele,
    ARRAY           *estif_global, 
    ARRAY           *emass_global,
    CONTAINER       *container);

void dyn_epot(
    FIELD           *actfield,
    INT              disnum,
    INTRA           *actintra,
    STRUCT_DYN_CALC *dynvar,
    DOUBLE          *deltaepot);

void dyn_ekin(
    FIELD           *actfield,
    SOLVAR          *actsolv,
    PARTITION       *actpart,
    INTRA           *actintra,
    CALC_ACTION     *action,
    CONTAINER       *container,
    INT              stiff_array,
    INT              mass_array);

void dyn_eout(
    STRUCT_DYN_CALC *dynvar,
    STRUCT_DYNAMIC  *sdyn,
    INTRA           *actintra,
    SOLVAR          *actsolv, 
    DIST_VECTOR     *dispi, /* converged incremental displacements */
    DIST_VECTOR     *rhs1,  /* load at time t                      */ 
    DIST_VECTOR     *rhs2,  /* load at time t-dt                   */ 
    DIST_VECTOR     *work0);

void dyn_ekin_local(
    ELEMENT         *actele,
    ARRAY           *emass,
    CONTAINER       *container);



#endif
