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

#include "../solver/solver.h"

/* include all the other prototypes */
#include "../pss_full/pss_prototypes.h"
#include "../output/output_prototypes.h"
#include "../math/math_prototypes.h"
#include "../parallel/parallel_prototypes.h"
#include "../visual/visual_prototypes.h"


/*----------------------------------------------------------------------*
 |  main_ccarat.c                                        m.gee 11/01    |
 *----------------------------------------------------------------------*/
INT main(INT argc, char *argv[]);

/*----------------------------------------------------------------------*
 |  global_ass_dof.c                                     m.gee 11/01    |
 *----------------------------------------------------------------------*/
void assign_dof(FIELD *actfield);

/*----------------------------------------------------------------------*
 |  global_ass_dof_ndis.c                                 genk 08/02    |
 *----------------------------------------------------------------------*/
void assign_dof_ndis(FIELD *actfield);

/*----------------------------------------------------------------------*
 |  global_cal_control.c                                 m.gee 11/01    |
 *----------------------------------------------------------------------*/
void ntacal(void);

/*----------------------------------------------------------------------*
 |  global_cal_control.c                                 genk 10/03     |
 *----------------------------------------------------------------------*/
void global_result_test(void); 

/*----------------------------------------------------------------------*
 |  cal_dyn_control.c                                    m.gee 11/01    |
 *----------------------------------------------------------------------*/
void caldyn(void);

/*----------------------------------------------------------------------*
 | cal_static_control.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void calsta(void);
void stalin(void);

/*----------------------------------------------------------------------*
 | cal_static_service.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void calstatserv_findcontroldof(
    FIELD     *actfield,
                                INT        control_node_global,
                                INT        control_dof,
                                NODE     **node,
                                INT       *cdof); 

/*----------------------------------------------------------------------*
 | cal_static_service.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void calstatserv_findreldofs(
    FIELD     *actfield,
                             INT       *reldisnode_ID,
                             INT       *reldis_dof,
                             INT        num_reldis,
                             INT       *reldof); 

/*----------------------------------------------------------------------*
 | cal_static_service.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void get_stepsize(
    INT         kstep,
                  STATIC_VAR *statvar); 

/*----------------------------------------------------------------------*
 | global_control.c                                      m.gee 11/01    |
 *----------------------------------------------------------------------*/
void ntam(INT argc, char *argv[]);


/*----------------------------------------------------------------------*
 | global_init_control.c                                 m.gee 11/01    |
 *----------------------------------------------------------------------*/
void ntaini(INT argc, char *argv[]);

/*----------------------------------------------------------------------*
 | global_inp_control.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void ntainp(void);

/*----------------------------------------------------------------------*
 | global_monitoring.c                                    genk 01/03    |
 *----------------------------------------------------------------------*/
void monitoring(FIELD *actfield,INT numf, INT actpos, INT actstep, DOUBLE time); 

/*----------------------------------------------------------------------*
 |  machine_hpux.c                                       m.gee 11/01    |
 *----------------------------------------------------------------------*/
void ntadev(INT argc, char *argv[]);

/*----------------------------------------------------------------------*
 |  map_node_find.c                                      m.gee 11/01    |
 *----------------------------------------------------------------------*/
void iscouple_find_node_comp(
    NODE     *actnode, 
    FIELD    *searchfield, 
    NODE    **partnernode,
    INT       coupleID,
    INT       dof);

void cheque_distance(
    DOUBLE   *x1,
    DOUBLE   *x2,
    DOUBLE    tol,
    INT      *ierr);

void find_assign_coupset(
    FIELD    *actfield, 
    INT       coupleID, 
    INT      *counter);
/*----------------------------------------------------------------------*
 |  dyn_timecurve.c                                      m.gee 02/02    |
 *----------------------------------------------------------------------*/
void dyn_init_curve(INT actcurve,
                   INT    nstep,
                   DOUBLE dt,
                   DOUBLE maxtime);
void dyn_facfromcurve(INT actcurve,
                   DOUBLE T,
                   DOUBLE *fac);
DOUBLE dyn_facexplcurve(INT actcurve,
                      DOUBLE T);		   
/*----------------------------------------------------------------------*
 |  inherit_insidedesign.c                                  m.gee 3/02  |
 *----------------------------------------------------------------------*/
void inherit_dirich_coup_indesign(void);
/*----------------------------------------------------------------------*
 |  inherit_design_dis.c                                    m.gee 3/02  |
 *----------------------------------------------------------------------*/
void inherit_design_dis_dirichlet(DISCRET *actdis);
void inherit_design_dis_couple(DISCRET *actdis);
void inherit_design_dis_fsicouple(DISCRET *actdis);
void inherit_design_dis_freesurf(DISCRET *actdis);
void inherit_design_dis_neum(DISCRET *actdis);
/*----------------------------------------------------------------------*
 |  inherit_design_ele.c                                   chfoe 01/04  |
 *----------------------------------------------------------------------*/
void inherit_design_ele(DISCRET *actdis);
/*----------------------------------------------------------------------*
 |  input_conditions.c                                  m.gee 11/01     |
 *----------------------------------------------------------------------*/
void inp_conditions(void);
/*----------------------------------------------------------------------*
 |  input_control_global.c                                  m.gee 11/01 |
 *----------------------------------------------------------------------*/
void inpctr(void);
void inpctrprob(void);
void inpctrdyn(void);
void inpctrstat(void);
void inpctreig(void);
void inpctr_dyn_struct(STRUCT_DYNAMIC *sdyn);
void inpctr_dyn_ale(ALE_DYNAMIC *adyn);
void inpctr_eig_struct(ALLEIG *alleig);
/*!---------------------------------------------------------------------
\brief input of the FLUID DYNAMIC block in the input-file

<pre>                                                         genk 03/02

In this routine the data in the FLUID DYNAMIC block of the input file
are read and stored in fdyn	       

</pre>
\param  *fdyn 	  FLUID_DATA       (o)	   
\return void                                                                       

------------------------------------------------------------------------*/
void inpctr_dyn_fluid(FLUID_DYNAMIC *fdyn);
/*!---------------------------------------------------------------------
\brief input of the FSI DYNAMIC block in the input-file

<pre>                                                         genk 09/02

In this routine the data in the FSI DYNAMIC block of the input file
are read and stored in fsidyn	       

</pre>
\param  *fsidyn 	  FSI_DATA       (o)	   
\return void                                                                       

------------------------------------------------------------------------*/
void inpctr_dyn_fsi(FSI_DYNAMIC *fsidyn);
/*----------------------------------------------------------------------*
 |  input_ctr_head.c                                  m.gee 11/01       |
 *----------------------------------------------------------------------*/
void inpctrhed(void);
void inptrace(void);
/*----------------------------------------------------------------------*
 |  input_curves.c                                  m.gee 11/01         |
 *----------------------------------------------------------------------*/
void inp_cond_curve(void);
void inp_read_curve(char *string);

/*----------------------------------------------------------------------*
    input_funct.c
 *----------------------------------------------------------------------*/
void inp_cond_funct(void);


/*----------------------------------------------------------------------*
 |  input_design.c                                  m.gee 11/01         |
 *----------------------------------------------------------------------*/
void inpdesign(void);
void inp_dnode(void);
void read_1_dnode(DNODE *dnode, INT readId);
void inp_dline(void);
void read_1_dline(DLINE *dline, INT readId);
void inp_dsurface(void);
void read_1_dsurf(DSURF *dsurf, INT readId);
void inp_dvolume(void);
void read_1_dvol(DVOL *dvol, INT readId);
void inp_designsize(void);
/*----------------------------------------------------------------------*
 |  input_design_top.c                                  m.gee 11/01     |
 *---------------------------------------------------------------------*/
void inpdesign_topology_design(void);
void inpdesign_topology_fe(void);
/*----------------------------------------------------------------------*
 |  input_material.c                                  m.gee 11/01       |
 *----------------------------------------------------------------------*/
void inp_material(void);
void inp_multimat(void);
/*----------------------------------------------------------------------*
 |  input_mesh.c                                  m.gee 11/01           |
 *----------------------------------------------------------------------*/
void inpfield(void);
void inp_assign_nodes(DISCRET *actdis);
void inpdis(FIELD *actfield);
void inpnodes(void);
void inp_struct_field(FIELD *structfield);
void inp_fluid_field(FIELD *fluidfield);
void inp_ale_field(FIELD *alefield);
/*----------------------------------------------------------------------*
 |  input_monitor                                         genk 01/03    |
 *----------------------------------------------------------------------*/
void inp_monitor(void);
/*----------------------------------------------------------------------*
 |  input_resultdescr                                       uk 05/04    |
 *----------------------------------------------------------------------*/
void inp_resultdescr(void);
/*----------------------------------------------------------------------*
 |  input_topology.c                                  m.gee 11/01       |
 *----------------------------------------------------------------------*/
void inp_topology(DISCRET *actdis);
void inp_detailed_topology(DISCRET   *actdis);


#ifdef CHECK_MAX
/* ====================================================================
 * file: check_max_sizes.c
 * ==================================================================== */
/*!----------------------------------------------------------------------
\brief check the values of the max sizes

<pre>                                                              mn 04/04
This routine determines the optimal values for maxele, maxnod, maxdofpernode
and maxgauss and compares those to the given values of the respective defines.
</pre>

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: --- 

*----------------------------------------------------------------------*/
void check_max_sizes(
    );


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
/*!---------------------------------------------------------------------                                         
\brief input of optimization data 

<pre>                                                          al  05/01      
</pre>  		 
\return void                                                                       

------------------------------------------------------------------------*/
void inpctropt(void);



/*!---------------------------------------------------------------------
\brief control execution of optimization

<pre>                                                          al  05/01      
</pre>  		 
\return void                                                                       

------------------------------------------------------------------------*/
void caloptmain(void);







#endif

