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
    DISCRET       *actdis);

void inherit_design_dis_couple(
    DISCRET       *actdis);

void inherit_design_dis_fsicouple(
    DISCRET       *actdis);

void inherit_design_dis_freesurf(
    DISCRET       *actdis);

void inherit_design_dis_neum(
    DISCRET       *actdis);


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


/*!---------------------------------------------------------------------
  \brief input of the FLUID DYNAMIC block in the input-file

  <pre>                                                         genk 03/02

  In this routine the data in the FLUID DYNAMIC block of the input file
  are read and stored in fdyn	       

  </pre>
  \param  *fdyn 	  FLUID_DATA       (o)	   
  \return void                                                                       

  ------------------------------------------------------------------------*/
void inpctr_dyn_fluid(
    FLUID_DYNAMIC  *fdyn);


/*!---------------------------------------------------------------------
  \brief input of the FSI DYNAMIC block in the input-file

  <pre>                                                         genk 09/02

  In this routine the data in the FSI DYNAMIC block of the input file
  are read and stored in fsidyn	       

  </pre>
  \param  *fsidyn 	  FSI_DATA       (o)	   
  \return void                                                                       

  ------------------------------------------------------------------------*/
void inpctr_dyn_fsi(
    FSI_DYNAMIC    *fsidyn);


/*----------------------------------------------------------------------*
  |  input_ctr_head.c                                  m.gee 11/01       |
 *----------------------------------------------------------------------*/
void inpctrhed(void);

void inptrace(void);


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


#endif

