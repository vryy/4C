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

#ifndef PROTOTYPES_H
#define PROTOTYPES_H


#include "../solver/solver.h"

/* include all the other prototypes */
#include "../pss_full/pss_prototypes.h"
#include "../input_full/input_prototypes.h"
#include "../output/output_prototypes.h"
#include "../math/math_prototypes.h"
#include "../parallel/parallel_prototypes.h"
#include "../visual/visual_prototypes.h"
#include "../structure/structure_prototypes.h"


/*----------------------------------------------------------------------*
  |  main_ccarat.c                                      m.gee 11/01    |
 *----------------------------------------------------------------------*/
INT main(
    INT            argc,
    char          *argv[]);


/*----------------------------------------------------------------------*
  |  global_ass_dof.c                                   m.gee 11/01    |
 *----------------------------------------------------------------------*/
void assign_dof(
    FIELD         *actfield);


/*----------------------------------------------------------------------*
  |  global_ass_dof_ndis.c                               genk 08/02    |
 *----------------------------------------------------------------------*/
void assign_dof_ndis(
    FIELD         *actfield);


/*----------------------------------------------------------------------*
  |  global_cal_control.c                               m.gee 11/01    |
 *----------------------------------------------------------------------*/
void ntacal(void);


/*----------------------------------------------------------------------*
  |  global_cal_control.c                               genk 10/03     |
 *----------------------------------------------------------------------*/
void global_result_test(void); 


/*----------------------------------------------------------------------*
  |  global_dyn_control.c                               m.gee 11/01    |
 *----------------------------------------------------------------------*/
void caldyn(void);


/*----------------------------------------------------------------------*
  | global_control.c                                    m.gee 11/01    |
 *----------------------------------------------------------------------*/
void ntam(
    INT            argc,
    char          *argv[]);


/*----------------------------------------------------------------------*
  | global_init_control.c                               m.gee 11/01    |
 *----------------------------------------------------------------------*/
void ntaini(
    INT            argc,
    char          *argv[]);


/*----------------------------------------------------------------------*
  | global_inp_control.c                                m.gee 11/01    |
 *----------------------------------------------------------------------*/
void ntainp(void);


/*----------------------------------------------------------------------*
  | global_monitoring.c                                  genk 01/03    |
 *----------------------------------------------------------------------*/
void monitoring(
    FIELD         *actfield,
    INT            numf,
    INT            actpos,
    INT            actstep,
    DOUBLE         time); 


/*----------------------------------------------------------------------*
  | global_open_files.c                                 m.gee 11/01    |
 *----------------------------------------------------------------------*/
void ntadev(
    INT            argc,
    char          *argv[]);


/*----------------------------------------------------------------------*
  |  global_node_find.c                                 m.gee 11/01    |
 *----------------------------------------------------------------------*/
void iscouple_find_node_comp(
    NODE          *actnode, 
    FIELD         *searchfield, 
    NODE         **partnernode,
    INT            coupleID,
    INT            dof);

void cheque_distance(
    DOUBLE        *x1,
    DOUBLE        *x2,
    DOUBLE         tol,
    INT           *ierr);

void find_assign_coupset(
    FIELD         *actfield, 
    INT            coupleID, 
    INT           *counter);


/*----------------------------------------------------------------------*
  | global_timecurve.c                                  m.gee 02/02    |
 *----------------------------------------------------------------------*/
void dyn_init_curve(
    INT            actcurve,
    INT            nstep,
    DOUBLE         dt,
    DOUBLE         maxtime);

void dyn_facfromcurve(
    INT            actcurve,
    DOUBLE         T,
    DOUBLE        *fac);

DOUBLE dyn_facexplcurve(
    INT            actcurve,
    DOUBLE         T);		   


/*----------------------------------------------------------------------*
  | global_calelm.c                                     m.gee 11/01    |
 *----------------------------------------------------------------------*/
void calelm(
    FIELD         *actfield,     /* active field */        
    SOLVAR        *actsolv,      /* active SOLVAR */
    PARTITION     *actpart,      /* my partition of this field */
    INTRA         *actintra,     /* my intra-communicator */
    INT            sysarray1,    /* number of first sparse system matrix */
    INT            sysarray2,    /* number of secnd system matrix, if present, else -1 */
    CONTAINER     *container,    /*!< contains variables defined in container.h */
    CALC_ACTION   *action);       /* calculation option passed to element routines */            

void calinit(
    FIELD         *actfield,   /* the active physical field */ 
    PARTITION     *actpart,    /* my partition of this field */
    CALC_ACTION   *action,
    CONTAINER     *container); /*!< contains variables defined in container.h */

void calreduce(
    FIELD         *actfield,    /* the active field */
    PARTITION     *actpart,     /* my partition of this field */
    INTRA         *actintra,    /* the field's intra-communicator */
    CALC_ACTION   *action,      /* action for element routines */
    CONTAINER     *container);  /*!< contains variables defined in container.h */               


/*----------------------------------------------------------------------*
  | global_calrhs.c                                     m.gee 11/01    |
 *----------------------------------------------------------------------*/
void calrhs(
    FIELD         *actfield,     /* the active field */
    SOLVAR        *actsolv,      /* the active SOLVAR */
    PARTITION     *actpart,      /* my partition of this field */
    INTRA         *actintra,     /* the field's intra-communicator */
    INT            actsysarray,  /* the active sparse array */
    DIST_VECTOR   *rhs1,         /* 2 dist. vectors for rhs */
    CALC_ACTION   *action,       /* action to be passed to element routines */
    CONTAINER     *container);   /*!< contains variables defined in container.h */

void rhs_point_neum(
    DOUBLE        *rhs,
    INT            dimrhs,
    PARTITION     *actpart);     


/*----------------------------------------------------------------------*
  | global_check_max.c                                        mn 04/04 |
 *----------------------------------------------------------------------*/
#ifdef CHECK_MAX
void check_max_sizes(void);
#endif





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

