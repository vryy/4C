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

#ifndef INPUT_PROTOTYPES_H
#define INPUT_PROTOTYPES_H


/*----------------------------------------------------------------------*
  |  inherit_insidedesign.c                                m.gee 3/02  |
 *----------------------------------------------------------------------*/
void inherit_dirich_coup_indesign(void);


/*----------------------------------------------------------------------*
  |  inherit_design_dis.c                                  m.gee 3/02  |
 *----------------------------------------------------------------------*/
void inherit_design_dis_dirichlet(
    FIELD         *actfield,
    DISCRET       *actdis);

void inherit_design_dis_couple(
    DISCRET       *actdis);

void inherit_design_dis_fsicouple(
    DISCRET       *actdis);

void inherit_design_dis_ssicouple(
    DISCRET       *actdis);

void inherit_design_dis_freesurf(
    DISCRET       *actdis);

void inherit_design_dis_neum(
    DISCRET       *actdis);

void inherit_design_dis_slipdirich(DISCRET *actdis);

/*----------------------------------------------------------------------*
  |  inherit_design_ele.c                                   chfoe 01/04  |
 *----------------------------------------------------------------------*/
void inherit_design_ele(
    DISCRET       *actdis);


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

void inpctr_dyn_struct(
    STRUCT_DYNAMIC *sdyn);

void inpctr_dyn_ale(
    ALE_DYNAMIC    *adyn);

void inpctr_eig_struct(
    ALLEIG         *alleig);

void inpctr_dyn_fluid(
    FLUID_DYNAMIC  *fdyn);

void inpctr_dyn_fsi(
    FSI_DYNAMIC    *fsidyn);

void inpctr_dyn_tsi(
    TSI_DYNAMIC    *tsidyn);


/*----------------------------------------------------------------------*
  |  input_ctr_head.c                                  m.gee 11/01       |
 *----------------------------------------------------------------------*/
void inpctrhed(void);


/*----------------------------------------------------------------------*
  |  input_curves.c                                  m.gee 11/01         |
 *----------------------------------------------------------------------*/
void inp_cond_curve(void);

void inp_read_curve(
    char          *string);


/*----------------------------------------------------------------------*
  input_funct.c
 *----------------------------------------------------------------------*/
void inp_cond_funct(void);


/*----------------------------------------------------------------------*
  |  input_design.c                                  m.gee 11/01         |
 *----------------------------------------------------------------------*/
void inpdesign(void);

void inp_dnode(void);

void read_1_dnode(
    DNODE         *dnode,
    INT            readId);

void inp_dline(void);

void read_1_dline(
    DLINE         *dline,
    INT            readId);

void inp_dsurface(void);

void read_1_dsurf(
    DSURF         *dsurf,
    INT            readId);

void inp_dvolume(void);

void read_1_dvol(
    DVOL          *dvol,
    INT            readId);

void inp_designsize(void);


/*----------------------------------------------------------------------*
  |  input_design_top.c                                  m.gee 11/01     |
 *---------------------------------------------------------------------*/
void inpdesign_topology_design(void);

void inpdesign_topology_fe(
    DISCRET        *actdis
    );


/*----------------------------------------------------------------------*
  |  input_material.c                                  m.gee 11/01       |
 *----------------------------------------------------------------------*/
void inp_material(void);

void inp_mat_nonconstparam(CHAR* param_w,  /* parameter key word */
                           MAT_PARAM_INTPOL* param_ipl,  /* parameter interpolation */
                           INT* param_n,   /* number of parameter components */
                           DOUBLE** param);  /* parameter vector */

void inp_multimat(void);

/*----------------------------------------------------------------------*
 |  input_surface_energy.c                               lw 11/06       |
 *----------------------------------------------------------------------*/
void input_surface_energy(void);

/*----------------------------------------------------------------------*
  |  input_mesh.c                                  m.gee 11/01           |
 *----------------------------------------------------------------------*/
void inpfield(void);

void inp_assign_nodes(
    DISCRET       *actdis);

void inpdis(
    FIELD         *actfield);

void inpnodes(void);

void inp_struct_field(
    FIELD         *structfield);

void inp_fluid_field(
    FIELD         *fluidfield);

void inp_ale_field(
    FIELD         *alefield);

void inp_therm_field(
    FIELD         *thermfield);


/*----------------------------------------------------------------------*
 | input of submesh: input_submesh.c                         ah 4/04    |
 *----------------------------------------------------------------------*/
void inp_submesh(void);

/*----------------------------------------------------------------------*
  |  input_monitor                                         genk 01/03    |
 *----------------------------------------------------------------------*/
void inp_monitor(void);

/*----------------------------------------------------------------------*
  |  input_locsys                                         genk 04/04    |
 *----------------------------------------------------------------------*/
void inp_cond_locsys(void);


/*----------------------------------------------------------------------*
  |  input_resultdescr                                       uk 05/04    |
 *----------------------------------------------------------------------*/
void inp_resultdescr(void);


/*----------------------------------------------------------------------*
  |  input_topology.c                                  m.gee 11/01       |
 *----------------------------------------------------------------------*/
void inp_topology(
    DISCRET       *actdis);

void inp_detailed_topology(
    DISCRET       *actdis);


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


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void input_ReadGlobalParameterList();

#endif

