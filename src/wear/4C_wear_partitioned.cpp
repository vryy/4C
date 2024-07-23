/*----------------------------------------------------------------------*/
/*! \file

\brief  Basis of all structure approaches with ale
        (Lagrangian step followed by Shape Evolution step )

\level 2


*/
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                                  farah 11/13 |
 *----------------------------------------------------------------------*/
#include "4C_wear_partitioned.hpp"

#include "4C_adapter_ale_wear.hpp"
#include "4C_ale_utils_mapextractor.hpp"
#include "4C_contact_abstract_strategy.hpp"
#include "4C_contact_defines.hpp"
#include "4C_contact_element.hpp"
#include "4C_contact_friction_node.hpp"
#include "4C_contact_integrator.hpp"
#include "4C_contact_interface.hpp"
#include "4C_contact_lagrange_strategy_wear.hpp"
#include "4C_contact_node.hpp"
#include "4C_contact_wear_interface.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_coupling_adapter_volmortar.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_fs3i_biofilm_fsi_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_ale.hpp"
#include "4C_inpar_wear.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_utils_densematrix_communication.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_linear_solver_method.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_mortar_manager_base.hpp"
#include "4C_so3_hex20.hpp"
#include "4C_so3_hex27.hpp"
#include "4C_so3_hex8.hpp"
#include "4C_so3_tet10.hpp"
#include "4C_so3_tet4.hpp"
#include "4C_structure_aux.hpp"
#include "4C_utils_parameter_list.hpp"
#include "4C_w1.hpp"
#include "4C_wear_utils.hpp"

#include <Epetra_SerialComm.h>
#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | constructor (public)                                     farah 05/13 |
 *----------------------------------------------------------------------*/
Wear::Partitioned::Partitioned(const Epetra_Comm& comm) : Algorithm(comm)
{
  const int ndim = Global::Problem::instance()->n_dim();

  // create ale-struct coupling
  const Epetra_Map* structdofmap = structure_field()->discretization()->node_row_map();
  const Epetra_Map* aledofmap = ale_field().discretization()->node_row_map();

  if (Core::UTILS::IntegralValue<bool>(Global::Problem::instance()->wear_params(), "MATCHINGGRID"))
  {
    // if there are two identical nodes (i.e. for initial contact) the nodes matching creates an
    // error !!!
    coupalestru_ = Teuchos::rcp(new Core::Adapter::Coupling());
    Teuchos::rcp_dynamic_cast<Core::Adapter::Coupling>(coupalestru_)
        ->setup_coupling(*ale_field().discretization(), *structure_field()->discretization(),
            *aledofmap, *structdofmap, ndim);
  }
  else
  {
    // Scheme: non matching meshes --> volumetric mortar coupling...
    coupalestru_ = Teuchos::rcp(new Core::Adapter::MortarVolCoupl());

    // projection ale -> structure : all ndim dofs (displacements)
    std::vector<int> coupleddof12 = std::vector<int>(ndim, 1);

    // projection structure -> ale : all ndim dofs (displacements)
    std::vector<int> coupleddof21 = std::vector<int>(ndim, 1);

    std::pair<int, int> dofset12(0, 0);
    std::pair<int, int> dofset21(0, 0);

    // init coupling
    Teuchos::rcp_dynamic_cast<Core::Adapter::MortarVolCoupl>(coupalestru_)
        ->init(ndim, Global::Problem::instance()->get_dis("ale"),
            Global::Problem::instance()->get_dis("structure"), &coupleddof12, &coupleddof21,
            &dofset12, &dofset21, Teuchos::null, false);

    // redistribute discretizations to meet needs of volmortar coupling
    //    Teuchos::rcp_dynamic_cast<Adapter::MortarVolCoupl>(coupalestru_)->redistribute();

    // setup projection matrices
    Teuchos::rcp_dynamic_cast<Core::Adapter::MortarVolCoupl>(coupalestru_)
        ->setup(Global::Problem::instance()->volmortar_params(),
            Global::Problem::instance()->cut_general_params());
  }

  // create interface coupling
  coupstrualei_ = Teuchos::rcp(new Core::Adapter::Coupling());
  coupstrualei_->setup_condition_coupling(*structure_field()->discretization(),
      structure_field()->interface()->ale_wear_cond_map(), *ale_field().discretization(),
      ale_field().interface()->Map(ale_field().interface()->cond_ale_wear), "AleWear", ndim);

  // initialize intern variables for wear
  wearnp_i_ = Teuchos::rcp(
      new Epetra_Vector(*ale_field().interface()->Map(ale_field().interface()->cond_ale_wear)),
      true);
  wearnp_ip_ = Teuchos::rcp(
      new Epetra_Vector(*ale_field().interface()->Map(ale_field().interface()->cond_ale_wear)),
      true);
  wearincr_ = Teuchos::rcp(
      new Epetra_Vector(*ale_field().interface()->Map(ale_field().interface()->cond_ale_wear)),
      true);
  delta_ale_ = Teuchos::rcp(new Epetra_Vector(ale_field().dispnp()->Map(), true));
  ale_i_ = Teuchos::rcp(new Epetra_Vector(ale_field().dispnp()->Map(), true));

  alepara_ = Global::Problem::instance()->ale_dynamic_params();
}


/*----------------------------------------------------------------------*
 | general time loop                                        farah 10/13 |
 *----------------------------------------------------------------------*/
void Wear::Partitioned::time_loop()
{
  // get wear paramter list
  const Teuchos::ParameterList& wearpara = Global::Problem::instance()->wear_params();
  double timeratio = wearpara.get<double>("WEAR_TIMERATIO");

  int counter = -1;
  bool alestep = false;

  // time loop
  while (not_finished())
  {
    if ((int)(step() / timeratio) > counter)
    {
      counter++;
      alestep = true;
    }

    if (Core::UTILS::IntegralValue<Inpar::Wear::WearCoupAlgo>(wearpara, "WEAR_COUPALGO") ==
        Inpar::Wear::wear_stagg)
      time_loop_stagg(alestep);
    else if (Core::UTILS::IntegralValue<Inpar::Wear::WearCoupAlgo>(wearpara, "WEAR_COUPALGO") ==
             Inpar::Wear::wear_iterstagg)
      time_loop_iter_stagg();
    else
      FOUR_C_THROW("Wear::TimeLoop: Algorithm not provided!");

    alestep = false;
  }  // time loop
}


/*----------------------------------------------------------------------*
 | time loop for staggered coupling                         farah 11/13 |
 *----------------------------------------------------------------------*/
void Wear::Partitioned::time_loop_iter_stagg()
{
  // counter and print header
  increment_time_and_step();
  print_header();

  // prepare time step for both fields
  prepare_time_step();

  bool converged = false;  // converged state?
  int iter = 0;            // iteration counter

  // stactic cast of mortar strategy to contact strategy
  Mortar::StrategyBase& strategy = cmtman_->get_strategy();
  Wear::LagrangeStrategyWear& cstrategy = static_cast<Wear::LagrangeStrategyWear&>(strategy);

  // reset waccu, wold and wcurr...
  cstrategy.update_wear_discret_iterate(false);

  /*************************************************************
   * Nonlinear iterations between Structure and ALE:           *
   * 1. Solve structure + contact to get wear                  *
   * 2. Apply wear increment (i+1 - i) onto ALE (add function) *
   * 3. Employ ALE disp incr (i+1 - i) and spat disp i to get  *
   *    abs mat disp for timestep n+1                          *
   * 4. Upadate spat disp from i to i+1                        *
   * 5. Check for convergence                                  *
   * 6. store ALE disp i = i+1                                 *
   *************************************************************/
  while (converged == false)
  {
    // 1. solution
    structure_field()->solve();

    // 2. wear as interface displacements in ale dofs
    Teuchos::RCP<Epetra_Vector> idisale_s, idisale_m;
    interface_disp(idisale_s, idisale_m);

    // merge the both wear vectors for master and slave side to one global vector
    merge_wear(idisale_s, idisale_m, wearincr_);

    // coupling of struct/mortar and ale dofs
    disp_coupling(wearincr_);

    if (get_comm().MyPID() == 0)
      std::cout << "========================= ALE STEP =========================" << std::endl;

    // do ale step
    ale_step(wearincr_);

    // 3. application of mesh displacements to structural field,
    // update material displacements
    update_mat_conf();

    // 4. update dispnp
    update_spat_conf();

    // 5. convergence check fot current iteration
    converged = convergence_check(iter);

    // store old wear
    cstrategy.update_wear_discret_iterate(true);

    ++iter;
  }  // end nonlin loop

  // update for structure and ale
  update();

  // output for structure and ale
  output();

  return;
}


/*----------------------------------------------------------------------*
 | time loop for oneway coupling                            farah 11/13 |
 *----------------------------------------------------------------------*/
void Wear::Partitioned::time_loop_stagg(bool alestep)
{
  // stactic cast of mortar strategy to contact strategy
  Mortar::StrategyBase& strategy = cmtman_->get_strategy();
  Wear::LagrangeStrategyWear& cstrategy = static_cast<Wear::LagrangeStrategyWear&>(strategy);

  // counter and print header
  increment_time_and_step();
  print_header();

  // prepare time step for both fields
  prepare_time_step();

  /********************************************************************/
  /* START LAGRANGE STEP                                              */
  /* structural lagrange step with contact                            */
  /********************************************************************/

  // solution
  structure_field()->solve();

  if (alestep)
  {
    if (get_comm().MyPID() == 0)
      std::cout << "========================= ALE STEP =========================" << std::endl;

    /********************************************************************/
    /* COUPLING                                                         */
    /* Wear from structure solve as dirichlet for ALE                   */
    /********************************************************************/

    // wear as interface displacements in interface dofs
    Teuchos::RCP<Epetra_Vector> idisale_s, idisale_m, idisale_global;
    interface_disp(idisale_s, idisale_m);

    // merge the both wear vectors for master and slave side to one global vector
    merge_wear(idisale_s, idisale_m, idisale_global);

    // coupling of struct/mortar and ale dofs
    disp_coupling(idisale_global);

    /********************************************************************/
    /* Shape Evolution STEP                                             */
    /* 1. mesh displacements due to wear from ALE system                */
    /* 2. mapping of results from "old" to "new" mesh                   */
    /********************************************************************/

    // do ale step
    ale_step(idisale_global);

    // update material displacements
    update_mat_conf();

    // update spatial displacements
    update_spat_conf();

    // reset wear
    cstrategy.update_wear_discret_iterate(false);
  }
  else
  {
    cstrategy.update_wear_discret_accumulation();
  }

  /********************************************************************/
  /* FINISH STEP:                                                     */
  /* Update and Write Output                                          */
  /********************************************************************/

  // update for structure and ale
  update();

  // output for structure and ale
  output();

  return;
}


/*----------------------------------------------------------------------*
 | prepare time step for ale and structure                  farah 11/13 |
 *----------------------------------------------------------------------*/
bool Wear::Partitioned::convergence_check(int iter)
{
  double Wincr = 0.0;
  double ALEincr = 0.0;
  wearincr_->Norm2(&Wincr);
  delta_ale_->Norm2(&ALEincr);

  if (get_comm().MyPID() == 0)
  {
    std::cout << "-----------------"
              << " Step " << iter + 1 << " --------------------" << std::endl;
    std::cout << "Wear incr.= " << Wincr << "         ALE incr.= " << ALEincr << std::endl;
    std::cout << "---------------------------------------------" << std::endl;
  }

  // TODO tolerance from input!!!
  // check loads
  if (abs(Wincr) < 1e-8 and abs(ALEincr) < 1e-8)
  {
    // reset vectors
    ale_i_->PutScalar(0.0);
    delta_ale_->PutScalar(0.0);

    return true;
  }

  if (iter > 50)
    FOUR_C_THROW(
        "Staggered solution scheme for ale-wear problem unconverged within 50 nonlinear iteration "
        "steps!");

  return false;
}


/*----------------------------------------------------------------------*
 | prepare time step for ale and structure                  farah 11/13 |
 *----------------------------------------------------------------------*/
void Wear::Partitioned::prepare_time_step()
{
  // predict and solve structural system
  structure_field()->prepare_time_step();

  // prepare ale output: increase time step
  ale_field().prepare_time_step();

  return;
}


/*----------------------------------------------------------------------*
 | update ale and structure                                 farah 11/13 |
 *---------------------------------------------------------------- ------*/
void Wear::Partitioned::update()
{
  // update at time step
  structure_field()->update();

  // update
  ale_field().update();

  return;
}


/*----------------------------------------------------------------------*
 | update spatial displacements                             farah 05/14 |
 *----------------------------------------------------------------------*/
void Wear::Partitioned::update_spat_conf()
{
  // mesh displacement from solution of ALE field in structural dofs
  // first perform transformation from ale to structure dofs
  Teuchos::RCP<Epetra_Vector> disalenp = ale_to_structure(ale_field().dispnp());
  Teuchos::RCP<Epetra_Vector> disalen = ale_to_structure(ale_field().dispn());

  // get structure dispnp vector
  Teuchos::RCP<Epetra_Vector> dispnp =
      structure_field()->write_access_dispnp();  // change to ExtractDispn() for overlap

  // get info about wear conf
  Inpar::Wear::WearShapeEvo wconf = Core::UTILS::IntegralValue<Inpar::Wear::WearShapeEvo>(
      Global::Problem::instance()->wear_params(), "WEAR_SHAPE_EVO");

  // for shape evol in spat conf
  if (wconf == Inpar::Wear::wear_se_sp)
  {
    int err = 0;
    // update per absolute vector
    err = dispnp->Update(1.0, *disalenp, 0.0);
    if (err != 0) FOUR_C_THROW("update wrong!");
  }
  // for shape evol in mat conf
  else if (wconf == Inpar::Wear::wear_se_mat)
  {
    // set state
    (structure_field()->discretization())->set_state(0, "displacement", dispnp);

    // set state
    (structure_field()->discretization())
        ->set_state(0, "material_displacement", structure_field()->disp_mat());

    // loop over all row nodes to fill graph
    for (int k = 0; k < structure_field()->discretization()->num_my_row_nodes(); ++k)
    {
      int gid = structure_field()->discretization()->node_row_map()->GID(k);
      Core::Nodes::Node* node = structure_field()->discretization()->g_node(gid);
      Core::Elements::Element** ElementPtr = node->elements();
      int numelement = node->num_element();

      const int numdof = structure_field()->discretization()->num_dof(node);

      // create Xmat for 3D problems
      std::vector<double> Xspatial(numdof);
      std::vector<double> Xmat(numdof);

      for (int dof = 0; dof < numdof; ++dof)
      {
        int dofgid = structure_field()->discretization()->dof(node, dof);
        int doflid = (dispnp->Map()).LID(dofgid);
        Xmat[dof] = node->x()[dof] + (*structure_field()->disp_mat())[doflid];
      }

      // create updated  Xspatial --> via nonlinear interpolation between nodes (like gp projection)
      advection_map(Xspatial.data(), Xmat.data(), ElementPtr, numelement, false);

      // store in dispmat
      for (int dof = 0; dof < numdof; ++dof)
      {
        int dofgid = structure_field()->discretization()->dof(node, dof);
        int doflid = (dispnp->Map()).LID(dofgid);
        (*dispnp)[doflid] = Xspatial[dof] - node->x()[dof];
      }
    }  // end row node loop
  }
  else
  {
    FOUR_C_THROW("Unknown wear configuration!");
  }

  return;
}


/*----------------------------------------------------------------------*
 | output ale and structure                                 farah 11/13 |
 *----------------------------------------------------------------------*/
void Wear::Partitioned::output()
{
  // calculate stresses, strains, energies
  constexpr bool force_prepare = false;
  structure_field()->prepare_output(force_prepare);

  // write strcture output to screen and files
  structure_field()->output();

  // output ale
  ale_field().output();

  return;
}


/*----------------------------------------------------------------------*
 | Perform Coupling from struct/mortar to ale dofs          farah 05/13 |
 | This is necessary due to the parallel redistribution                 |
 | of the contact interface                                             |
 *----------------------------------------------------------------------*/
void Wear::Partitioned::disp_coupling(Teuchos::RCP<Epetra_Vector>& disinterface)
{
  // Teuchos::RCP<Epetra_Vector> aledofs = Teuchos::rcp(new
  // Epetra_Vector(*ale_field().Interface()->Map(ale_field().Interface()->cond_ale_wear)),true);
  Teuchos::RCP<Epetra_Vector> strudofs =
      Teuchos::rcp(new Epetra_Vector(*structure_field()->interface()->ale_wear_cond_map()), true);

  // change the parallel distribution from mortar interface to structure
  Core::LinAlg::export_to(*disinterface, *strudofs);

  // perform coupling to ale dofs
  disinterface.reset();
  disinterface = coupstrualei_->master_to_slave(strudofs);

  return;
}


/*----------------------------------------------------------------------*
 | Merge wear from slave and master surface to one          farah 06/13 |
 | wear vector                                                          |
 *----------------------------------------------------------------------*/
void Wear::Partitioned::merge_wear(Teuchos::RCP<Epetra_Vector>& disinterface_s,
    Teuchos::RCP<Epetra_Vector>& disinterface_m, Teuchos::RCP<Epetra_Vector>& disinterface_g)
{
  Mortar::StrategyBase& strategy = cmtman_->get_strategy();
  CONTACT::AbstractStrategy& cstrategy = static_cast<CONTACT::AbstractStrategy&>(strategy);
  std::vector<Teuchos::RCP<CONTACT::Interface>> interface = cstrategy.contact_interfaces();
  Teuchos::RCP<Wear::WearInterface> winterface =
      Teuchos::rcp_dynamic_cast<Wear::WearInterface>(interface[0]);
  if (winterface == Teuchos::null) FOUR_C_THROW("Casting to WearInterface returned null!");

  disinterface_g = Teuchos::rcp(new Epetra_Vector(*winterface->discret().dof_row_map()), true);
  Teuchos::RCP<Epetra_Vector> auxvector =
      Teuchos::rcp(new Epetra_Vector(*winterface->discret().dof_row_map()), true);

  Core::LinAlg::export_to(*disinterface_s, *disinterface_g);
  Core::LinAlg::export_to(*disinterface_m, *auxvector);

  int err = 0;
  err = disinterface_g->Update(1.0, *auxvector, true);
  if (err != 0) FOUR_C_THROW("update wrong!");

  return;
}


/*----------------------------------------------------------------------*
 | Vector of interface displacements in ALE dofs            farah 05/13 |
 | Currently just for 1 interface                                       |
 *----------------------------------------------------------------------*/
void Wear::Partitioned::interface_disp(
    Teuchos::RCP<Epetra_Vector>& disinterface_s, Teuchos::RCP<Epetra_Vector>& disinterface_m)
{
  // get info about wear side
  Inpar::Wear::WearSide wside = Core::UTILS::IntegralValue<Inpar::Wear::WearSide>(
      Global::Problem::instance()->wear_params(), "WEAR_SIDE");

  // get info about wear type
  Inpar::Wear::WearType wtype = Core::UTILS::IntegralValue<Inpar::Wear::WearType>(
      Global::Problem::instance()->wear_params(), "WEARTYPE");

  // get info about wear coeff conf
  Inpar::Wear::WearCoeffConf wcoeffconf = Core::UTILS::IntegralValue<Inpar::Wear::WearCoeffConf>(
      Global::Problem::instance()->wear_params(), "WEARCOEFF_CONF");

  if (interfaces_.size() > 1)
    FOUR_C_THROW("Wear algorithm not able to handle more than 1 interface yet!");

  //------------------------------------------------
  // Wear coefficient constant in material config: -
  //------------------------------------------------
  if (wcoeffconf == Inpar::Wear::wear_coeff_mat)
  {
    // redistribute int. according to spatial interfaces!
    redistribute_mat_interfaces();

    // 1. pull back slave wear to material conf.
    wear_pull_back_slave(disinterface_s);

    // 2. pull back master wear to material conf.
    if (wside == Inpar::Wear::wear_both)
    {
      wear_pull_back_master(disinterface_m);
    }
    else
    {
      // zeroes
      Teuchos::RCP<Epetra_Map> masterdofs = interfaces_[0]->master_row_dofs();
      disinterface_m = Teuchos::rcp(new Epetra_Vector(*masterdofs, true));
    }
  }
  //------------------------------------------------
  // Wear coefficient constant in spatial config:  -
  //------------------------------------------------
  else if (wcoeffconf == Inpar::Wear::wear_coeff_sp)
  {
    // postproc wear for spatial conf.
    wear_spatial_slave(disinterface_s);

    if (wside == Inpar::Wear::wear_both and wtype == Inpar::Wear::wear_primvar)
    {
      wear_spatial_master(disinterface_m);
    }
    else if (wside == Inpar::Wear::wear_both and wtype == Inpar::Wear::wear_intstate)
    {
      // redistribute int. according to spatial interfaces!
      redistribute_mat_interfaces();
      wear_spatial_master_map(disinterface_s, disinterface_m);
    }
    else
    {
      // zeroes
      Teuchos::RCP<Epetra_Map> masterdofs = interfaces_[0]->master_row_dofs();
      disinterface_m = Teuchos::rcp(new Epetra_Vector(*masterdofs, true));
    }
  }
  //------------------------------------------------
  // ERROR                                         -
  //------------------------------------------------
  else
  {
    FOUR_C_THROW("Chosen wear configuration not supported!");
  }

  return;
}


/*----------------------------------------------------------------------*
 | Wear in spatial conf.                                    farah 09/14 |
 *----------------------------------------------------------------------*/
void Wear::Partitioned::wear_spatial_master_map(
    Teuchos::RCP<Epetra_Vector>& disinterface_s, Teuchos::RCP<Epetra_Vector>& disinterface_m)
{
  if (disinterface_s == Teuchos::null) FOUR_C_THROW("no slave wear for mapping!");

  // stactic cast of mortar strategy to contact strategy
  Mortar::StrategyBase& strategy = cmtman_->get_strategy();
  Wear::LagrangeStrategyWear& cstrategy = static_cast<Wear::LagrangeStrategyWear&>(strategy);

  for (int i = 0; i < (int)interfaces_.size(); ++i)
  {
    Teuchos::RCP<Wear::WearInterface> winterface =
        Teuchos::rcp_dynamic_cast<Wear::WearInterface>(interfacesMat_[i]);
    if (winterface == Teuchos::null) FOUR_C_THROW("Casting to WearInterface returned null!");

    Teuchos::RCP<Epetra_Map> masterdofs = interfaces_[i]->master_row_dofs();
    Teuchos::RCP<Epetra_Map> slavedofs = interfaces_[i]->slave_row_dofs();
    Teuchos::RCP<Epetra_Map> activedofs = interfaces_[i]->active_dofs();

    disinterface_m = Teuchos::rcp(new Epetra_Vector(*masterdofs, true));

    // different wear coefficients on both sides...
    double wearcoeff_s = interfaces_[i]->interface_params().get<double>("WEARCOEFF", 0.0);
    double wearcoeff_m = interfaces_[i]->interface_params().get<double>("WEARCOEFF_MASTER", 0.0);
    if (wearcoeff_s < 1e-12) FOUR_C_THROW("wcoeff negative!!!");

    double fac = wearcoeff_m / (wearcoeff_s);

    Teuchos::RCP<Epetra_Vector> wear_master = Teuchos::rcp(new Epetra_Vector(*masterdofs, true));

    cstrategy.m_matrix()->multiply(true, *disinterface_s, *wear_master);

    // 1. set state to material displacement state
    winterface->set_state(
        Mortar::state_new_displacement, *structure_field()->write_access_dispnp());

    // 2. initialize
    winterface->initialize();

    // 3. calc N and areas
    winterface->set_element_areas();
    winterface->evaluate_nodal_normals();

    // 6. init data container for d2 mat
    const Teuchos::RCP<Epetra_Map> masternodesmat =
        Core::LinAlg::AllreduceEMap(*(winterface->master_row_nodes()));

    for (int i = 0; i < masternodesmat->NumMyElements();
         ++i)  // for (int i=0;i<MasterRowNodes()->NumMyElements();++i)
    {
      int gid = masternodesmat->GID(i);
      Core::Nodes::Node* node = winterface->discret().g_node(gid);
      if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
      CONTACT::FriNode* cnode = dynamic_cast<CONTACT::FriNode*>(node);

      if (cnode->is_slave() == false)
      {
        // reset nodal Mortar maps
        for (int j = 0; j < (int)((cnode->wear_data().get_d2()).size()); ++j)
          (cnode->wear_data().get_d2())[j].clear();

        (cnode->wear_data().get_d2()).resize(0);
      }
    }

    // 8. evaluate dmat
    Teuchos::RCP<Core::LinAlg::SparseMatrix> dmat = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
        *masterdofs, 100, true, false, Core::LinAlg::SparseMatrix::FE_MATRIX));

    for (int j = 0; j < winterface->master_col_elements()->NumMyElements(); ++j)
    {
      int gid = winterface->master_col_elements()->GID(j);
      Core::Elements::Element* ele = winterface->discret().g_element(gid);
      if (!ele) FOUR_C_THROW("Cannot find ele with gid %", gid);
      CONTACT::Element* cele = dynamic_cast<CONTACT::Element*>(ele);

      Teuchos::RCP<CONTACT::Integrator> integrator = Teuchos::rcp(
          new CONTACT::Integrator(winterface->interface_params(), cele->shape(), get_comm()));

      integrator->integrate_d(*cele, get_comm());
    }

    // 10. assemble dmat
    winterface->assemble_d2(*dmat);

    // 12. complete dmat
    dmat->complete();

    Teuchos::ParameterList solvparams;
    Core::UTILS::AddEnumClassToParameterList<Core::LinearSolver::SolverType>(
        "SOLVER", Core::LinearSolver::SolverType::umfpack, solvparams);
    Core::LinAlg::Solver solver(solvparams, get_comm(),
        Global::Problem::instance()->solver_params_callback(),
        Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
            Global::Problem::instance()->io_params(), "VERBOSITY"));

    Core::LinAlg::SolverParams solver_params;
    solver_params.refactor = true;
    solver.solve(dmat->epetra_operator(), disinterface_m, wear_master, solver_params);

    disinterface_m->Scale(-fac);
  }

  // bye
  return;
}


/*----------------------------------------------------------------------*
 | Wear in spatial conf.                                    farah 09/14 |
 *----------------------------------------------------------------------*/
void Wear::Partitioned::wear_spatial_master(Teuchos::RCP<Epetra_Vector>& disinterface_m)
{
  // get info about wear conf
  Inpar::Wear::WearTimeScale wtime = Core::UTILS::IntegralValue<Inpar::Wear::WearTimeScale>(
      Global::Problem::instance()->wear_params(), "WEAR_TIMESCALE");

  for (int i = 0; i < (int)interfaces_.size(); ++i)
  {
    Teuchos::RCP<Epetra_Map> masterdofs = interfaces_[i]->master_row_dofs();
    disinterface_m = Teuchos::rcp(new Epetra_Vector(*masterdofs, true));

    // FIRST: get the wear values and the normal directions for the interface
    // loop over all slave row nodes on the current interface
    for (int j = 0; j < interfaces_[i]->master_row_nodes()->NumMyElements(); ++j)
    {
      int gid = interfaces_[i]->master_row_nodes()->GID(j);
      Core::Nodes::Node* node = interfaces_[i]->discret().g_node(gid);
      if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
      CONTACT::FriNode* frinode = dynamic_cast<CONTACT::FriNode*>(node);

      // be aware of problem dimension
      int numdof = frinode->num_dof();
      if (dim_ != numdof) FOUR_C_THROW("Inconsistency Dim <-> NumDof");

      // nodal normal vector and wear
      double nn[3];
      double wear = 0.0;

      for (int j = 0; j < 3; ++j) nn[j] = frinode->mo_data().n()[j];

      if (wtime == Inpar::Wear::wear_time_different)
      {
        if (abs(frinode->wear_data().wcurr()[0] + frinode->wear_data().waccu()[0]) > 1e-12)
          wear = frinode->wear_data().wcurr()[0] + frinode->wear_data().waccu()[0];
        else
          wear = 0.0;
      }
      else
      {
        if (abs(frinode->wear_data().wcurr()[0]) > 1e-12)
          wear = frinode->wear_data().wcurr()[0];
        else
          wear = 0.0;
      }


      // find indices for DOFs of current node in Epetra_Vector
      // and put node values (normal and tangential stress components) at these DOFs
      std::vector<int> locindex(dim_);

      for (int dof = 0; dof < dim_; ++dof)
      {
        locindex[dof] = (disinterface_m->Map()).LID(frinode->dofs()[dof]);
        (*disinterface_m)[locindex[dof]] = -wear * nn[dof];
      }
    }
  }

  // bye
  return;
}


/*----------------------------------------------------------------------*
 | Wear in spatial conf.                                    farah 09/14 |
 *----------------------------------------------------------------------*/
void Wear::Partitioned::wear_spatial_slave(Teuchos::RCP<Epetra_Vector>& disinterface_s)
{
  // stactic cast of mortar strategy to contact strategy
  Mortar::StrategyBase& strategy = cmtman_->get_strategy();
  Wear::LagrangeStrategyWear& cstrategy = static_cast<Wear::LagrangeStrategyWear&>(strategy);

  Inpar::Wear::WearType wtype = Core::UTILS::IntegralValue<Inpar::Wear::WearType>(
      Global::Problem::instance()->wear_params(), "WEARTYPE");

  Inpar::Wear::WearTimInt wtimint = Core::UTILS::IntegralValue<Inpar::Wear::WearTimInt>(
      Global::Problem::instance()->wear_params(), "WEARTIMINT");

  Inpar::Wear::WearTimeScale wtime = Core::UTILS::IntegralValue<Inpar::Wear::WearTimeScale>(
      Global::Problem::instance()->wear_params(), "WEAR_TIMESCALE");

  if (!(wtype == Inpar::Wear::wear_intstate and wtimint == Inpar::Wear::wear_impl))
    cstrategy.store_nodal_quantities(Mortar::StrategyBase::weightedwear);

  for (int i = 0; i < (int)interfaces_.size(); ++i)
  {
    Teuchos::RCP<Epetra_Map> slavedofs = interfaces_[i]->slave_row_dofs();
    Teuchos::RCP<Epetra_Map> activedofs = interfaces_[i]->active_dofs();

    // additional spatial displacements
    disinterface_s = Teuchos::rcp(new Epetra_Vector(*slavedofs, true));

    // FIRST: get the wear values and the normal directions for the interface
    // loop over all slave row nodes on the current interface
    for (int j = 0; j < interfaces_[i]->slave_row_nodes()->NumMyElements(); ++j)
    {
      int gid = interfaces_[i]->slave_row_nodes()->GID(j);
      Core::Nodes::Node* node = interfaces_[i]->discret().g_node(gid);
      if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
      CONTACT::FriNode* frinode = dynamic_cast<CONTACT::FriNode*>(node);

      // be aware of problem dimension
      int numdof = frinode->num_dof();
      if (dim_ != numdof) FOUR_C_THROW("Inconsistency Dim <-> NumDof");

      // nodal normal vector and wear
      double nn[3];
      double wear = 0.0;

      for (int j = 0; j < 3; ++j) nn[j] = frinode->mo_data().n()[j];

      if (wtype == Inpar::Wear::wear_primvar)
      {
        if (wtime == Inpar::Wear::wear_time_different)
        {
          if (abs(frinode->wear_data().wcurr()[0] + frinode->wear_data().waccu()[0]) > 1e-12)
            wear = frinode->wear_data().wcurr()[0] + frinode->wear_data().waccu()[0];
          else
            wear = 0.0;
        }
        else
        {
          if (abs(frinode->wear_data().wcurr()[0]) > 1e-12)
            wear = frinode->wear_data().wcurr()[0];
          else
            wear = 0.0;
        }
      }
      else if (wtype == Inpar::Wear::wear_intstate)
      {
        wear = frinode->wear_data().weighted_wear();
      }

      // find indices for DOFs of current node in Epetra_Vector
      // and put node values (normal and tangential stress components) at these DOFs
      std::vector<int> locindex(dim_);

      for (int dof = 0; dof < dim_; ++dof)
      {
        locindex[dof] = (disinterface_s->Map()).LID(frinode->dofs()[dof]);
        (*disinterface_s)[locindex[dof]] = -wear * nn[dof];
      }
    }

    // un-weight for internal state approach
    if (wtype == Inpar::Wear::wear_intstate)
    {
      Teuchos::RCP<Core::LinAlg::SparseMatrix> daa, dai, dia, dii;
      Teuchos::RCP<Epetra_Map> gidofs;
      Core::LinAlg::SplitMatrix2x2(
          cstrategy.d_matrix(), activedofs, gidofs, activedofs, gidofs, daa, dai, dia, dii);

      Teuchos::RCP<Epetra_Vector> wear_vectora = Teuchos::rcp(new Epetra_Vector(*activedofs, true));
      Teuchos::RCP<Epetra_Vector> wear_vectori = Teuchos::rcp(new Epetra_Vector(*gidofs));
      Core::LinAlg::split_vector(
          *slavedofs, *disinterface_s, activedofs, wear_vectora, gidofs, wear_vectori);

      Teuchos::RCP<Epetra_Vector> zref = Teuchos::rcp(new Epetra_Vector(*activedofs));

      // solve with default solver
      Teuchos::ParameterList solvparams;
      Core::UTILS::AddEnumClassToParameterList<Core::LinearSolver::SolverType>(
          "SOLVER", Core::LinearSolver::SolverType::umfpack, solvparams);
      Core::LinAlg::Solver solver(solvparams, get_comm(),
          Global::Problem::instance()->solver_params_callback(),
          Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
              Global::Problem::instance()->io_params(), "VERBOSITY"));

      if (activedofs->NumMyElements())
      {
        Core::LinAlg::SolverParams solver_params;
        solver_params.refactor = true;

        solver.solve(daa->epetra_operator(), zref, wear_vectora, solver_params);
      }

      // different wear coefficients on both sides...
      double wearcoeff_s = interfaces_[0]->interface_params().get<double>("WEARCOEFF", 0.0);
      double wearcoeff_m = interfaces_[0]->interface_params().get<double>("WEARCOEFF_MASTER", 0.0);
      if (wearcoeff_s < 1e-12) FOUR_C_THROW("wcoeff negative!!!");
      double fac = wearcoeff_s / (wearcoeff_s + wearcoeff_m);
      zref->Scale(fac);

      disinterface_s = Teuchos::rcp(new Epetra_Vector(*slavedofs));
      Core::LinAlg::export_to(*zref, *disinterface_s);
    }
  }

  // bye
  return;
}


/*----------------------------------------------------------------------*
 | Redistribute material interfaces acc. to cur interf.     farah 09/14 |
 *----------------------------------------------------------------------*/
void Wear::Partitioned::redistribute_mat_interfaces()
{
  // barrier
  get_comm().Barrier();

  // loop over all interfaces
  for (int m = 0; m < (int)interfaces_.size(); ++m)
  {
    int redistglobal = 0;
    int redistlocal = 0;
    if (interfaces_[m]->is_redistributed()) redistlocal++;

    get_comm().SumAll(&redistlocal, &redistglobal, 1);
    get_comm().Barrier();


    if (redistglobal > 0)
    {
      if (get_comm().MyPID() == 0)
      {
        std::cout << "===========================================" << std::endl;
        std::cout << "=======    Redistribute Mat. Int.   =======" << std::endl;
        std::cout << "===========================================" << std::endl;
      }
      Teuchos::RCP<Wear::WearInterface> winterface =
          Teuchos::rcp_dynamic_cast<Wear::WearInterface>(interfacesMat_[m]);

      // export nodes and elements to the row map
      winterface->discret().export_row_nodes(*interfaces_[m]->discret().node_row_map());
      winterface->discret().export_row_elements(*interfaces_[m]->discret().element_row_map());

      // export nodes and elements to the column map (create ghosting)
      winterface->discret().export_column_nodes(*interfaces_[m]->discret().node_col_map());
      winterface->discret().export_column_elements(*interfaces_[m]->discret().element_col_map());

      winterface->fill_complete(true);
      winterface->print_parallel_distribution();

      if (get_comm().MyPID() == 0)
      {
        std::cout << "===========================================" << std::endl;
        std::cout << "==============     Done!     ==============" << std::endl;
        std::cout << "===========================================" << std::endl;
      }
    }
  }

  // barrier
  get_comm().Barrier();
  return;
}

/*----------------------------------------------------------------------*
 | Pull-Back wear: W = w * ds/dS * N                        farah 09/14 |
 *----------------------------------------------------------------------*/
void Wear::Partitioned::wear_pull_back_slave(Teuchos::RCP<Epetra_Vector>& disinterface_s)
{
  // stactic cast of mortar strategy to contact strategy
  Mortar::StrategyBase& strategy = cmtman_->get_strategy();
  Wear::LagrangeStrategyWear& cstrategy = dynamic_cast<Wear::LagrangeStrategyWear&>(strategy);

  Inpar::Wear::WearType wtype = Core::UTILS::IntegralValue<Inpar::Wear::WearType>(
      Global::Problem::instance()->wear_params(), "WEARTYPE");

  Inpar::Wear::WearTimInt wtimint = Core::UTILS::IntegralValue<Inpar::Wear::WearTimInt>(
      Global::Problem::instance()->wear_params(), "WEARTIMINT");

  Inpar::Wear::WearTimeScale wtime = Core::UTILS::IntegralValue<Inpar::Wear::WearTimeScale>(
      Global::Problem::instance()->wear_params(), "WEAR_TIMESCALE");

  if (!(wtype == Inpar::Wear::wear_intstate and wtimint == Inpar::Wear::wear_impl))
    cstrategy.store_nodal_quantities(Mortar::StrategyBase::weightedwear);

  // loop over all interfaces
  for (int m = 0; m < (int)interfaces_.size(); ++m)
  {
    Teuchos::RCP<Wear::WearInterface> winterface =
        Teuchos::rcp_dynamic_cast<Wear::WearInterface>(interfaces_[m]);
    if (winterface == Teuchos::null) FOUR_C_THROW("Casting to WearInterface returned null!");

    // get slave row dofs as map
    Teuchos::RCP<Epetra_Map> slavedofs = winterface->slave_row_dofs();
    // additional spatial displacements
    disinterface_s = Teuchos::rcp(new Epetra_Vector(*slavedofs, true));

    // call material interfaces and evaluate!
    // 1. set state to material displacement state
    interfacesMat_[m]->set_state(Mortar::state_new_displacement, *structure_field()->disp_mat());

    // 2. initialize
    interfacesMat_[m]->initialize();

    // 3. calc N and areas
    interfacesMat_[m]->set_element_areas();
    interfacesMat_[m]->evaluate_nodal_normals();

    // 4. calc -w*N
    for (int j = 0; j < winterface->slave_row_nodes()->NumMyElements(); ++j)
    {
      int gid = winterface->slave_row_nodes()->GID(j);
      Core::Nodes::Node* node = winterface->discret().g_node(gid);
      if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
      CONTACT::FriNode* frinode = dynamic_cast<CONTACT::FriNode*>(node);

      int gidm = interfacesMat_[m]->slave_row_nodes()->GID(j);
      Core::Nodes::Node* nodem = interfacesMat_[m]->discret().g_node(gidm);
      if (!nodem) FOUR_C_THROW("Cannot find node with gid %", gidm);
      CONTACT::FriNode* frinodem = dynamic_cast<CONTACT::FriNode*>(nodem);

      // be aware of problem dimension
      int numdof = frinode->num_dof();
      if (dim_ != numdof) FOUR_C_THROW("Inconsistency Dim <-> NumDof");

      // nodal normal vector and wear
      double nn[3];
      double wear = 0.0;

      // get material normal
      for (int j = 0; j < 3; ++j) nn[j] = frinodem->mo_data().n()[j];

      // primary variable approach:
      if (wtype == Inpar::Wear::wear_primvar)
      {
        if (wtime == Inpar::Wear::wear_time_different)
        {
          if (abs(frinode->wear_data().wcurr()[0] + frinode->wear_data().waccu()[0]) > 1e-12)
            wear = frinode->wear_data().wcurr()[0] + frinode->wear_data().waccu()[0];
          else
            wear = 0.0;
        }
        else
        {
          if (abs(frinode->wear_data().wcurr()[0]) > 1e-12)
            wear = frinode->wear_data().wcurr()[0];
          else
            wear = 0.0;
        }
      }
      // internal state variable approach:
      else if (wtype == Inpar::Wear::wear_intstate)
      {
        wear = frinode->wear_data().weighted_wear();
      }

      // find indices for DOFs of current node in Epetra_Vector
      // and put node values (normal and tangential stress components) at these DOFs
      std::vector<int> locindex(dim_);

      for (int dof = 0; dof < dim_; ++dof)
      {
        locindex[dof] = (disinterface_s->Map()).LID(frinode->dofs()[dof]);
        (*disinterface_s)[locindex[dof]] = -wear * nn[dof];
      }
    }

    // 5. evaluate dmat
    Teuchos::RCP<Core::LinAlg::SparseMatrix> dmat =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*slavedofs, 10));

    for (int j = 0; j < interfacesMat_[m]->slave_col_elements()->NumMyElements(); ++j)
    {
      int gid = interfacesMat_[m]->slave_col_elements()->GID(j);
      Core::Elements::Element* ele = interfacesMat_[m]->discret().g_element(gid);
      if (!ele) FOUR_C_THROW("Cannot find ele with gid %", gid);
      CONTACT::Element* cele = dynamic_cast<CONTACT::Element*>(ele);

      Teuchos::RCP<CONTACT::Integrator> integrator = Teuchos::rcp(new CONTACT::Integrator(
          interfacesMat_[m]->interface_params(), cele->shape(), get_comm()));

      integrator->integrate_d(*cele, get_comm());
    }

    // 6. assemble dmat
    interfacesMat_[m]->assemble_d(*dmat);

    // 7. complete dmat
    dmat->complete();

    // 8. area trafo:
    if (wtype == Inpar::Wear::wear_primvar)
    {
      // multiply current D matrix with current wear
      Teuchos::RCP<Epetra_Vector> forcecurr = Teuchos::rcp(new Epetra_Vector(*slavedofs));
      cstrategy.d_matrix()->multiply(false, *disinterface_s, *forcecurr);

      // LM in reference / current configuration
      Teuchos::RCP<Epetra_Vector> zref = Teuchos::rcp(new Epetra_Vector(*slavedofs));

      // solve with default solver
      Teuchos::ParameterList solvparams;
      Core::UTILS::AddEnumClassToParameterList<Core::LinearSolver::SolverType>(
          "SOLVER", Core::LinearSolver::SolverType::umfpack, solvparams);
      Core::LinAlg::Solver solver(solvparams, get_comm(),
          Global::Problem::instance()->solver_params_callback(),
          Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
              Global::Problem::instance()->io_params(), "VERBOSITY"));

      Core::LinAlg::SolverParams solver_params;
      solver_params.refactor = true;
      solver.solve(dmat->epetra_operator(), zref, forcecurr, solver_params);


      // store reference LM into global vector and nodes
      disinterface_s = zref;
    }
    else if (wtype == Inpar::Wear::wear_intstate)
    {
      Teuchos::RCP<Epetra_Vector> zref = Teuchos::rcp(new Epetra_Vector(*slavedofs));

      // solve with default solver
      Teuchos::ParameterList solvparams;
      Core::UTILS::AddEnumClassToParameterList<Core::LinearSolver::SolverType>(
          "SOLVER", Core::LinearSolver::SolverType::umfpack, solvparams);
      Core::LinAlg::Solver solver(solvparams, get_comm(),
          Global::Problem::instance()->solver_params_callback(),
          Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
              Global::Problem::instance()->io_params(), "VERBOSITY"));

      Core::LinAlg::SolverParams solver_params;
      solver_params.refactor = true;
      solver.solve(dmat->epetra_operator(), zref, disinterface_s, solver_params);


      // store reference LM into global vector and nodes
      disinterface_s = zref;

      // different wear coefficients on both sides...
      double wearcoeff_s = interfaces_[0]->interface_params().get<double>("WEARCOEFF", 0.0);
      double wearcoeff_m = interfaces_[0]->interface_params().get<double>("WEARCOEFF_MASTER", 0.0);
      if (wearcoeff_s < 1e-12) FOUR_C_THROW("wcoeff negative!!!");

      double fac = wearcoeff_s / (wearcoeff_s + wearcoeff_m);
      disinterface_s->Scale(fac);
    }
    else
    {
      FOUR_C_THROW("wrong wear type!");
    }
  }

  // bye
  return;
}


/*----------------------------------------------------------------------*
 | Pull-Back wear: W = w * ds/dS * N                        farah 09/14 |
 *----------------------------------------------------------------------*/
void Wear::Partitioned::wear_pull_back_master(Teuchos::RCP<Epetra_Vector>& disinterface_m)
{
  Inpar::Wear::WearType wtype = Core::UTILS::IntegralValue<Inpar::Wear::WearType>(
      Global::Problem::instance()->wear_params(), "WEARTYPE");

  Inpar::Wear::WearTimeScale wtime = Core::UTILS::IntegralValue<Inpar::Wear::WearTimeScale>(
      Global::Problem::instance()->wear_params(), "WEAR_TIMESCALE");

  // loop over all interfaces
  for (int m = 0; m < (int)interfaces_.size(); ++m)
  {
    Teuchos::RCP<Wear::WearInterface> winterface =
        Teuchos::rcp_dynamic_cast<Wear::WearInterface>(interfaces_[m]);
    if (winterface == Teuchos::null) FOUR_C_THROW("Casting to WearInterface returned null!");

    Teuchos::RCP<Wear::WearInterface> winterfaceMat =
        Teuchos::rcp_dynamic_cast<Wear::WearInterface>(interfacesMat_[m]);
    if (winterfaceMat == Teuchos::null) FOUR_C_THROW("Casting to WearInterface returned null!");

    // get slave row dofs as map
    Teuchos::RCP<Epetra_Map> masterdofs = winterface->master_row_dofs();
    // additional spatial displacements
    disinterface_m = Teuchos::rcp(new Epetra_Vector(*masterdofs, true));

    // call material interfaces and evaluate!
    // 1. set state to material displacement state
    winterfaceMat->set_state(Mortar::state_new_displacement, *structure_field()->disp_mat());

    // 2. initialize
    winterfaceMat->initialize();

    // 3. calc N and areas
    winterfaceMat->set_element_areas();
    winterfaceMat->evaluate_nodal_normals();

    // 4. calc -w*N
    for (int j = 0; j < winterface->master_row_nodes()->NumMyElements(); ++j)
    {
      int gid = winterface->master_row_nodes()->GID(j);
      Core::Nodes::Node* node = winterface->discret().g_node(gid);
      if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
      CONTACT::FriNode* frinode = dynamic_cast<CONTACT::FriNode*>(node);

      int gidm = interfacesMat_[m]->master_row_nodes()->GID(j);
      Core::Nodes::Node* nodem = interfacesMat_[m]->discret().g_node(gidm);
      if (!nodem) FOUR_C_THROW("Cannot find node with gid %", gidm);
      CONTACT::FriNode* frinodem = dynamic_cast<CONTACT::FriNode*>(nodem);

      // be aware of problem dimension
      int numdof = frinode->num_dof();
      if (dim_ != numdof) FOUR_C_THROW("Inconsistency Dim <-> NumDof");

      // nodal normal vector and wear
      double nn[3];
      double wear = 0.0;

      // get material normal
      for (int j = 0; j < 3; ++j) nn[j] = frinodem->mo_data().n()[j];

      if (wtype == Inpar::Wear::wear_primvar)
      {
        if (wtime == Inpar::Wear::wear_time_different)
        {
          if (abs(frinode->wear_data().wcurr()[0] + frinode->wear_data().waccu()[0]) > 1e-12)
            wear = frinode->wear_data().wcurr()[0] + frinode->wear_data().waccu()[0];
          else
            wear = 0.0;
        }
        else
        {
          if (abs(frinode->wear_data().wcurr()[0]) > 1e-12)
            wear = frinode->wear_data().wcurr()[0];
          else
            wear = 0.0;
        }
      }
      else if (wtype == Inpar::Wear::wear_intstate)
      {
        wear = frinode->wear_data().weighted_wear();
      }

      // find indices for DOFs of current node in Epetra_Vector
      // and put node values (normal and tangential stress components) at these DOFs
      std::vector<int> locindex(dim_);

      for (int dof = 0; dof < dim_; ++dof)
      {
        locindex[dof] = (disinterface_m->Map()).LID(frinode->dofs()[dof]);
        (*disinterface_m)[locindex[dof]] = -wear * nn[dof];
      }
    }

    // 5. init data container for d2 curr
    const Teuchos::RCP<Epetra_Map> masternodes =
        Core::LinAlg::AllreduceEMap(*(winterface->master_row_nodes()));

    for (int i = 0; i < masternodes->NumMyElements();
         ++i)  // for (int i=0;i<MasterRowNodes()->NumMyElements();++i)
    {
      int gid = masternodes->GID(i);
      Core::Nodes::Node* node = winterface->discret().g_node(gid);
      if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
      CONTACT::FriNode* cnode = dynamic_cast<CONTACT::FriNode*>(node);

      if (cnode->is_slave() == false)
      {
        // reset nodal Mortar maps
        for (int j = 0; j < (int)((cnode->wear_data().get_d2()).size()); ++j)
          (cnode->wear_data().get_d2())[j].clear();

        (cnode->wear_data().get_d2()).resize(0);
      }
    }

    // 6. init data container for d2 mat
    const Teuchos::RCP<Epetra_Map> masternodesmat =
        Core::LinAlg::AllreduceEMap(*(winterfaceMat->master_row_nodes()));

    for (int i = 0; i < masternodesmat->NumMyElements();
         ++i)  // for (int i=0;i<MasterRowNodes()->NumMyElements();++i)
    {
      int gid = masternodesmat->GID(i);
      Core::Nodes::Node* node = winterfaceMat->discret().g_node(gid);
      if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
      CONTACT::FriNode* cnode = dynamic_cast<CONTACT::FriNode*>(node);

      if (cnode->is_slave() == false)
      {
        // reset nodal Mortar maps
        for (int j = 0; j < (int)((cnode->wear_data().get_d2()).size()); ++j)
          (cnode->wear_data().get_d2())[j].clear();

        (cnode->wear_data().get_d2()).resize(0);
      }
    }

    // 7. evaluate dcur
    Teuchos::RCP<Core::LinAlg::SparseMatrix> dcur = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
        *masterdofs, 100, true, false, Core::LinAlg::SparseMatrix::FE_MATRIX));
    for (int j = 0; j < winterface->master_col_elements()->NumMyElements(); ++j)
    {
      int gid = winterface->master_col_elements()->GID(j);
      Core::Elements::Element* ele = winterface->discret().g_element(gid);
      if (!ele) FOUR_C_THROW("Cannot find ele with gid %", gid);
      CONTACT::Element* cele = dynamic_cast<CONTACT::Element*>(ele);

      Teuchos::RCP<CONTACT::Integrator> integrator = Teuchos::rcp(
          new CONTACT::Integrator(winterface->interface_params(), cele->shape(), get_comm()));

      integrator->integrate_d(*cele, get_comm());
    }

    // 8. evaluate dmat
    Teuchos::RCP<Core::LinAlg::SparseMatrix> dmat = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
        *masterdofs, 100, true, false, Core::LinAlg::SparseMatrix::FE_MATRIX));

    for (int j = 0; j < winterfaceMat->master_col_elements()->NumMyElements(); ++j)
    {
      int gid = winterfaceMat->master_col_elements()->GID(j);
      Core::Elements::Element* ele = winterfaceMat->discret().g_element(gid);
      if (!ele) FOUR_C_THROW("Cannot find ele with gid %", gid);
      CONTACT::Element* cele = dynamic_cast<CONTACT::Element*>(ele);

      Teuchos::RCP<CONTACT::Integrator> integrator = Teuchos::rcp(
          new CONTACT::Integrator(winterfaceMat->interface_params(), cele->shape(), get_comm()));

      integrator->integrate_d(*cele, get_comm());
    }

    // 9. assemble dcur
    winterface->assemble_d2(*dcur);

    // 10. assemble dmat
    winterfaceMat->assemble_d2(*dmat);

    // 11. complete dcur
    dcur->complete();

    // 12. complete dmat
    dmat->complete();

    // 13. area trafo:
    if (wtype == Inpar::Wear::wear_primvar)
    {
      // multiply current D matrix with current wear
      Teuchos::RCP<Epetra_Vector> forcecurr = Teuchos::rcp(new Epetra_Vector(*masterdofs));
      dcur->multiply(false, *disinterface_m, *forcecurr);

      // LM in reference / current configuration
      Teuchos::RCP<Epetra_Vector> zref = Teuchos::rcp(new Epetra_Vector(*masterdofs));

      // solve with default solver
      Teuchos::ParameterList solvparams;
      Core::UTILS::AddEnumClassToParameterList<Core::LinearSolver::SolverType>(
          "SOLVER", Core::LinearSolver::SolverType::umfpack, solvparams);
      Core::LinAlg::Solver solver(solvparams, get_comm(),
          Global::Problem::instance()->solver_params_callback(),
          Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
              Global::Problem::instance()->io_params(), "VERBOSITY"));

      Core::LinAlg::SolverParams solver_params;
      solver_params.refactor = true;
      solver.solve(dmat->epetra_operator(), zref, forcecurr, solver_params);


      // store reference LM into global vector and nodes
      disinterface_m = zref;
    }
    else if (wtype == Inpar::Wear::wear_intstate)
    {
      FOUR_C_THROW("not working yet!");
      Teuchos::RCP<Epetra_Vector> zref = Teuchos::rcp(new Epetra_Vector(*masterdofs));

      // solve with default solver
      Teuchos::ParameterList solvparams;
      Core::UTILS::AddEnumClassToParameterList<Core::LinearSolver::SolverType>(
          "SOLVER", Core::LinearSolver::SolverType::umfpack, solvparams);
      Core::LinAlg::Solver solver(solvparams, get_comm(),
          Global::Problem::instance()->solver_params_callback(),
          Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
              Global::Problem::instance()->io_params(), "VERBOSITY"));

      Core::LinAlg::SolverParams solver_params;
      solver_params.refactor = true;
      solver.solve(dmat->epetra_operator(), zref, disinterface_m, solver_params);

      // store reference LM into global vector and nodes
      disinterface_m = zref;
    }
    else
    {
      FOUR_C_THROW("wrong wear type!");
    }
  }

  // bye
  return;
}


/*----------------------------------------------------------------------*
 | Application of mesh displacement to material conf         farah 04/15|
 *----------------------------------------------------------------------*/
void Wear::Partitioned::update_mat_conf()
{
  // mesh displacement from solution of ALE field in structural dofs
  // first perform transformation from ale to structure dofs
  Teuchos::RCP<Epetra_Vector> disalenp = ale_to_structure(ale_field().dispnp());

  // vector of current spatial displacements
  Teuchos::RCP<const Epetra_Vector> dispnp =
      structure_field()->dispnp();  // change to ExtractDispn() for overlap

  // material displacements
  Teuchos::RCP<Epetra_Vector> dismat = Teuchos::rcp(new Epetra_Vector(dispnp->Map()), true);

  // set state
  (structure_field()->discretization())->set_state(0, "displacement", dispnp);

  // set state
  (structure_field()->discretization())
      ->set_state(0, "material_displacement", structure_field()->disp_mat());

  // get info about wear conf
  Inpar::Wear::WearShapeEvo wconf = Core::UTILS::IntegralValue<Inpar::Wear::WearShapeEvo>(
      Global::Problem::instance()->wear_params(), "WEAR_SHAPE_EVO");

  // if shape evol. in mat conf: ale dispnp = material displ.
  if (wconf == Inpar::Wear::wear_se_mat)
  {
    // just information for user
    int err = 0;
    err = delta_ale_->Update(-1.0, *ale_i_, 0.0);
    if (err != 0) FOUR_C_THROW("update wrong!");
    err = delta_ale_->Update(1.0, *ale_field().dispnp(), 1.0);
    if (err != 0) FOUR_C_THROW("update wrong!");
    err = ale_i_->Update(1.0, *ale_field().dispnp(), 0.0);
    if (err != 0) FOUR_C_THROW("update wrong!");

    // important vector to update mat conf
    Teuchos::RCP<Epetra_Vector> dismat_struct =
        Teuchos::rcp(new Epetra_Vector(dispnp->Map()), true);

    Core::LinAlg::export_to(*disalenp, *dismat_struct);

    err = dismat->Update(1.0, *dismat_struct, 0.0);
    if (err != 0) FOUR_C_THROW("update wrong!");
  }
  // if shape evol. in spat conf: advection map!
  else if (wconf == Inpar::Wear::wear_se_sp)
  {
    int err = 0;
    err = disalenp->Update(-1.0, *dispnp, 1.0);
    if (err != 0) FOUR_C_THROW("update wrong!");
    err = delta_ale_->Update(1.0, *structure_to_ale(disalenp), 0.0);
    if (err != 0) FOUR_C_THROW("update wrong!");

    // loop over all row nodes to fill graph
    for (int k = 0; k < structure_field()->discretization()->num_my_row_nodes(); ++k)
    {
      int gid = structure_field()->discretization()->node_row_map()->GID(k);
      Core::Nodes::Node* node = structure_field()->discretization()->g_node(gid);
      Core::Elements::Element** ElementPtr = node->elements();
      int numelement = node->num_element();

      const int numdof = structure_field()->discretization()->num_dof(node);

      // create Xmat for 3D problems
      std::vector<double> XMat(numdof);
      std::vector<double> XMesh(numdof);

      for (int dof = 0; dof < numdof; ++dof)
      {
        int dofgid = structure_field()->discretization()->dof(node, dof);
        int doflid = (dispnp->Map()).LID(dofgid);
        XMesh[dof] = node->x()[dof] + (*dispnp)[doflid] + (*disalenp)[doflid];
      }

      // create updated  XMat --> via nonlinear interpolation between nodes (like gp projection)
      advection_map(XMat.data(), XMesh.data(), ElementPtr, numelement, true);

      // store in dispmat
      for (int dof = 0; dof < numdof; ++dof)
      {
        int dofgid = structure_field()->discretization()->dof(node, dof);
        int doflid = (dispnp->Map()).LID(dofgid);
        (*dismat)[doflid] = XMat[dof] - node->x()[dof];
      }
    }  // end row node loop
  }

  // apply material displacements to structural field
  // if advection map is not succesful --> use old xmat
  structure_field()->apply_dis_mat(dismat);

  // bye
  return;
}


/*----------------------------------------------------------------------*
 | material coordinates evaluated from spatial ones         farah 12/13 |
 *----------------------------------------------------------------------*/
void Wear::Partitioned::advection_map(double* Xtarget,  // out
    double* Xsource,                                    // in
    Core::Elements::Element** ElementPtr,               // in
    int numelements,                                    // in
    bool spatialtomaterial)                             // in
{
  // get problem dimension
  const int ndim = Global::Problem::instance()->n_dim();

  // define source and target configuration
  std::string sourceconf;
  std::string targetconf;

  if (spatialtomaterial)
  {
    sourceconf = "displacement";
    targetconf = "material_displacement";
  }
  else
  {
    sourceconf = "material_displacement";
    targetconf = "displacement";
  }

  // found element the spatial coordinate lies in
  bool found = false;

  // parameter space coordinates
  double e[3];
  double ge1 = 1e12;
  double ge2 = 1e12;
  double ge3 = 1e12;
  int gele = 0;

  // get state
  Teuchos::RCP<const Epetra_Vector> dispsource =
      (structure_field()->discretization())->get_state(sourceconf);
  Teuchos::RCP<const Epetra_Vector> disptarget =
      (structure_field()->discretization())->get_state(targetconf);

  // loop over adjacent elements
  for (int jele = 0; jele < numelements; jele++)
  {
    // get element
    Core::Elements::Element* actele = ElementPtr[jele];

    // get element location vector, dirichlet flags and ownerships
    Core::Elements::Element::LocationArray la(1);
    actele->location_vector(*(structure_field()->discretization()), la, false);

    if (ndim == 2)
    {
      if (actele->shape() == Core::FE::CellType::quad4)
        Wear::UTILS::av<Core::FE::CellType::quad4>(
            actele, Xtarget, Xsource, dispsource, disptarget, la[0].lm_, found, e);
      else if (actele->shape() == Core::FE::CellType::quad8)
        Wear::UTILS::av<Core::FE::CellType::quad8>(
            actele, Xtarget, Xsource, dispsource, disptarget, la[0].lm_, found, e);
      else if (actele->shape() == Core::FE::CellType::quad9)
        Wear::UTILS::av<Core::FE::CellType::quad9>(
            actele, Xtarget, Xsource, dispsource, disptarget, la[0].lm_, found, e);
      else if (actele->shape() == Core::FE::CellType::tri3)
        Wear::UTILS::av<Core::FE::CellType::tri3>(
            actele, Xtarget, Xsource, dispsource, disptarget, la[0].lm_, found, e);
      else if (actele->shape() == Core::FE::CellType::tri6)
        Wear::UTILS::av<Core::FE::CellType::tri6>(
            actele, Xtarget, Xsource, dispsource, disptarget, la[0].lm_, found, e);
      else
        FOUR_C_THROW("shape function not supported!");

      // checks if the spatial coordinate lies within this element
      // if yes, returns the material displacements
      // w1ele->AdvectionMapElement(XMat1,XMat2,XMesh1,XMesh2,disp,dispmat, la,found,e1,e2);

      // if parameter space coord. 'e' does not lie within any element (i.e. found = false),
      // then jele is the element lying closest near the considered spatial point.
      if (found == false)
      {
        if (abs(ge1) > 1.0 and abs(e[0]) < abs(ge1))
        {
          ge1 = e[0];
          gele = jele;
        }
        if (abs(ge2) > 1.0 and abs(e[1]) < abs(ge2))
        {
          ge2 = e[1];
          gele = jele;
        }
      }
    }
    else
    {
      if (actele->element_type() == Discret::ELEMENTS::SoHex8Type::instance())
        Wear::UTILS::av<Core::FE::CellType::hex8>(
            actele, Xtarget, Xsource, dispsource, disptarget, la[0].lm_, found, e);
      else if (actele->element_type() == Discret::ELEMENTS::SoHex20Type::instance())
        Wear::UTILS::av<Core::FE::CellType::hex20>(
            actele, Xtarget, Xsource, dispsource, disptarget, la[0].lm_, found, e);
      else if (actele->element_type() == Discret::ELEMENTS::SoHex27Type::instance())
        Wear::UTILS::av<Core::FE::CellType::hex27>(
            actele, Xtarget, Xsource, dispsource, disptarget, la[0].lm_, found, e);
      else if (actele->element_type() == Discret::ELEMENTS::SoTet4Type::instance())
        Wear::UTILS::av<Core::FE::CellType::tet4>(
            actele, Xtarget, Xsource, dispsource, disptarget, la[0].lm_, found, e);
      else if (actele->element_type() == Discret::ELEMENTS::SoTet10Type::instance())
        Wear::UTILS::av<Core::FE::CellType::tet10>(
            actele, Xtarget, Xsource, dispsource, disptarget, la[0].lm_, found, e);
      else
        FOUR_C_THROW("element type not supported!");

      // if parameter space coord. 'e' does not lie within any element (i.e. found = false),
      // then 'gele = jele' is the element lying closest near the considered spatial point.
      if (found == false)
      {
        if (abs(ge1) > 1.0 and abs(e[0]) < abs(ge1))
        {
          ge1 = e[0];
          gele = jele;
        }
        if (abs(ge2) > 1.0 and abs(e[1]) < abs(ge2))
        {
          ge2 = e[1];
          gele = jele;
        }
        if (abs(ge3) > 1.0 and abs(e[2]) < abs(ge3))
        {
          ge3 = e[2];
          gele = jele;
        }
      }
    }

    // leave when element is found
    if (found == true) return;
  }  // end loop over adj elements

  // ****************************************
  //  if displ not into elements: get
  //  Xtarget from closest element 'gele'
  // ****************************************
  Core::Elements::Element* actele = ElementPtr[gele];

  // get element location vector, dirichlet flags and ownerships
  Core::Elements::Element::LocationArray la(1);
  actele->location_vector(*(structure_field()->discretization()), la, false);

  if (ndim == 2)
  {
    if (actele->shape() == Core::FE::CellType::quad4)
      Wear::UTILS::av<Core::FE::CellType::quad4>(
          actele, Xtarget, Xsource, dispsource, disptarget, la[0].lm_, found, e);
    else if (actele->shape() == Core::FE::CellType::quad8)
      Wear::UTILS::av<Core::FE::CellType::quad8>(
          actele, Xtarget, Xsource, dispsource, disptarget, la[0].lm_, found, e);
    else if (actele->shape() == Core::FE::CellType::quad9)
      Wear::UTILS::av<Core::FE::CellType::quad9>(
          actele, Xtarget, Xsource, dispsource, disptarget, la[0].lm_, found, e);
    else if (actele->shape() == Core::FE::CellType::tri3)
      Wear::UTILS::av<Core::FE::CellType::tri3>(
          actele, Xtarget, Xsource, dispsource, disptarget, la[0].lm_, found, e);
    else if (actele->shape() == Core::FE::CellType::tri6)
      Wear::UTILS::av<Core::FE::CellType::tri6>(
          actele, Xtarget, Xsource, dispsource, disptarget, la[0].lm_, found, e);
    else
      FOUR_C_THROW("shape function not supported!");
  }
  else
  {
    if (actele->element_type() == Discret::ELEMENTS::SoHex8Type::instance())
      Wear::UTILS::av<Core::FE::CellType::hex8>(
          actele, Xtarget, Xsource, dispsource, disptarget, la[0].lm_, found, e);
    else if (actele->element_type() == Discret::ELEMENTS::SoHex20Type::instance())
      Wear::UTILS::av<Core::FE::CellType::hex20>(
          actele, Xtarget, Xsource, dispsource, disptarget, la[0].lm_, found, e);
    else if (actele->element_type() == Discret::ELEMENTS::SoHex27Type::instance())
      Wear::UTILS::av<Core::FE::CellType::hex27>(
          actele, Xtarget, Xsource, dispsource, disptarget, la[0].lm_, found, e);
    else if (actele->element_type() == Discret::ELEMENTS::SoTet4Type::instance())
      Wear::UTILS::av<Core::FE::CellType::tet4>(
          actele, Xtarget, Xsource, dispsource, disptarget, la[0].lm_, found, e);
    else if (actele->element_type() == Discret::ELEMENTS::SoTet10Type::instance())
      Wear::UTILS::av<Core::FE::CellType::tet10>(
          actele, Xtarget, Xsource, dispsource, disptarget, la[0].lm_, found, e);
    else
      FOUR_C_THROW("element type not supported!");
  }

  // bye
  return;
}


/*----------------------------------------------------------------------*
 | Perform ALE step                                         farah 11/13 |
 *----------------------------------------------------------------------*/
void Wear::Partitioned::ale_step(Teuchos::RCP<Epetra_Vector> idisale_global)
{
  // get info about ale dynamic
  // Inpar::ALE::AleDynamic aletype =
  Core::UTILS::IntegralValue<Inpar::ALE::AleDynamic>(params_ale(), "ALE_TYPE");

  // get info about wear conf
  Inpar::Wear::WearShapeEvo wconf = Core::UTILS::IntegralValue<Inpar::Wear::WearShapeEvo>(
      Global::Problem::instance()->wear_params(), "WEAR_SHAPE_EVO");

  //  if(aletype != Inpar::ALE::solid)
  //    FOUR_C_THROW("ERORR: Chosen ALE type not supported!");

  if (wconf == Inpar::Wear::wear_se_sp)
  {
    //    Teuchos::RCP<Epetra_Vector> dispnpstru = structure_to_ale(
    //        structure_field()->Dispnp());
    //
    //    FS3I::Biofilm::UTILS::updateMaterialConfigWithALE_Disp(
    //        ale_field().write_access_discretization(),
    //        dispnpstru );
    //
    //    ale_field().WriteAccessDispnp()->Update(0.0, *(dispnpstru), 0.0);
    //
    //    // application of interface displacements as dirichlet conditions
    //    //ale_field().apply_interface_displacements(idisale_global);
    //
    //    // solve time step
    //    ale_field().TimeStep(ALE::UTILS::MapExtractor::dbc_set_wear);
    //
    //    ale_field().WriteAccessDispnp()->Update(1.0, *(dispnpstru), 1.0);


    Teuchos::RCP<Epetra_Vector> dispnpstru = structure_to_ale(structure_field()->dispnp());

    ale_field().write_access_dispnp()->Update(1.0, *(dispnpstru), 0.0);

    // application of interface displacements as dirichlet conditions
    ale_field().apply_interface_displacements(idisale_global);

    // solve time step
    ale_field().time_step(ALE::UTILS::MapExtractor::dbc_set_wear);
  }
  // classical lin in mat. conf --> not correct at all
  else if (wconf == Inpar::Wear::wear_se_mat)
  {
    // application of interface displacements as dirichlet conditions
    ale_field().apply_interface_displacements(idisale_global);

    // solve time step
    ale_field().time_step(ALE::UTILS::MapExtractor::dbc_set_wear);
  }
  else
    FOUR_C_THROW("Chosen wear configuration not supported!");

  return;
}


/*----------------------------------------------------------------------*
 | transform from ale to structure map                      farah 11/13 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Wear::Partitioned::ale_to_structure(
    Teuchos::RCP<Epetra_Vector> vec) const
{
  return coupalestru_->master_to_slave(vec);
}


/*----------------------------------------------------------------------*
 | transform from ale to structure map                      farah 11/13 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Wear::Partitioned::ale_to_structure(
    Teuchos::RCP<const Epetra_Vector> vec) const
{
  return coupalestru_->master_to_slave(vec);
}


/*----------------------------------------------------------------------*
 | transform from ale to structure map                      farah 11/13 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Wear::Partitioned::structure_to_ale(
    Teuchos::RCP<Epetra_Vector> vec) const
{
  return coupalestru_->slave_to_master(vec);
}


/*----------------------------------------------------------------------*
 | transform from ale to structure map                      farah 11/13 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Wear::Partitioned::structure_to_ale(
    Teuchos::RCP<const Epetra_Vector> vec) const
{
  return coupalestru_->slave_to_master(vec);
}


/*----------------------------------------------------------------------*
 | read restart information for given time step (public)    farah 10/13 |
 *----------------------------------------------------------------------*/
void Wear::Partitioned::read_restart(int step)
{
  structure_field()->read_restart(step);
  ale_field().read_restart(step);
  set_time_step(structure_field()->time_old(), step);

  return;
}
/*----------------------------------------------------------------------*/

FOUR_C_NAMESPACE_CLOSE
