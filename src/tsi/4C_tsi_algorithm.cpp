/*----------------------------------------------------------------------*/
/*! \file

\brief  Basis of all TSI algorithms that perform a coupling between the linear
        momentum equation and the heat conduction equation
\level 2
*/

/*----------------------------------------------------------------------*
 | definitions                                               dano 12/09 |
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                                   dano 12/09 |
 *----------------------------------------------------------------------*/
#include "4C_tsi_algorithm.hpp"

#include "4C_adapter_str_factory.hpp"
#include "4C_adapter_str_structure_new.hpp"
#include "4C_adapter_str_wrapper.hpp"
#include "4C_adapter_thermo.hpp"
#include "4C_contact_lagrange_strategy.hpp"
#include "4C_contact_lagrange_strategy_tsi.hpp"
#include "4C_contact_meshtying_contact_bridge.hpp"
#include "4C_contact_nitsche_strategy_tsi.hpp"
#include "4C_contact_strategy_factory.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_coupling_adapter_mortar.hpp"
#include "4C_coupling_adapter_volmortar.hpp"
#include "4C_coupling_volmortar_utils.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_tsi.hpp"
#include "4C_io.hpp"
#include "4C_mortar_multifield_coupling.hpp"
#include "4C_so3_base.hpp"
#include "4C_solid_3D_ele.hpp"
#include "4C_structure_new_model_evaluator_contact.hpp"
#include "4C_structure_new_model_evaluator_structure.hpp"
#include "4C_thermo_element.hpp"
#include "4C_tsi_defines.hpp"
#include "4C_tsi_utils.hpp"

FOUR_C_NAMESPACE_OPEN

//! Note: The order of calling the two BaseAlgorithm-constructors is
//! important here! In here control file entries are written. And these entries
//! define the order in which the filters handle the Discretizations, which in
//! turn defines the dof number ordering of the Discretizations.
/*----------------------------------------------------------------------*
 | constructor (public)                                      dano 12/09 |
 *----------------------------------------------------------------------*/
TSI::Algorithm::Algorithm(const Epetra_Comm& comm)
    : AlgorithmBase(comm, Global::Problem::instance()->tsi_dynamic_params()),
      dispnp_(Teuchos::null),
      tempnp_(Teuchos::null),
      matchinggrid_(Core::UTILS::IntegralValue<bool>(
          Global::Problem::instance()->tsi_dynamic_params(), "MATCHINGGRID")),
      volcoupl_(Teuchos::null)
{
  // access the structural discretization
  Teuchos::RCP<Core::FE::Discretization> structdis =
      Global::Problem::instance()->get_dis("structure");
  // access the thermo discretization
  Teuchos::RCP<Core::FE::Discretization> thermodis = Global::Problem::instance()->get_dis("thermo");

  // get the problem instance
  Global::Problem* problem = Global::Problem::instance();
  // get the restart step
  const int restart = problem->restart();

  if (!matchinggrid_)
  {
    // Scheme: non matching meshes --> volumetric mortar coupling...
    volcoupl_ = Teuchos::rcp(new Core::Adapter::MortarVolCoupl());

    Teuchos::RCP<Core::VolMortar::UTILS::DefaultMaterialStrategy> materialstrategy =
        Teuchos::rcp(new TSI::UTILS::TSIMaterialStrategy());
    // init coupling adapter projection matrices
    volcoupl_->init(Global::Problem::instance()->n_dim(), structdis, thermodis, nullptr, nullptr,
        nullptr, nullptr, materialstrategy);
    // redistribute discretizations to meet needs of volmortar coupling
    Teuchos::ParameterList binning_params = Global::Problem::instance()->binning_strategy_params();
    Core::UTILS::AddEnumClassToParameterList<Core::FE::ShapeFunctionType>(
        "spatial_approximation_type", Global::Problem::instance()->spatial_approximation_type(),
        binning_params);
    auto element_filter = [](const Core::Elements::Element* element)
    { return Core::Binstrategy::Utils::SpecialElement::none; };
    auto rigid_sphere_radius = [](const Core::Elements::Element* element) { return 0.0; };
    auto correct_beam_center_node = [](const Core::Nodes::Node* node) { return node; };
    volcoupl_->redistribute(binning_params, Global::Problem::instance()->output_control_file(),
        element_filter, rigid_sphere_radius, correct_beam_center_node);
    // setup projection matrices
    volcoupl_->setup(Global::Problem::instance()->volmortar_params(),
        Global::Problem::instance()->cut_general_params());
  }

  if (Core::UTILS::IntegralValue<Inpar::Solid::IntegrationStrategy>(
          Global::Problem::instance()->structural_dynamic_params(), "INT_STRATEGY") ==
      Inpar::Solid::int_old)
    FOUR_C_THROW("old structural time integration no longer supported in tsi");
  else
  {
    Teuchos::RCP<Adapter::ThermoBaseAlgorithm> thermo =
        Teuchos::rcp(new Adapter::ThermoBaseAlgorithm(
            Global::Problem::instance()->tsi_dynamic_params(), thermodis));
    thermo_ = thermo->thermo_fieldrcp();

    //  // access structural dynamic params list which will be possibly modified while creating the
    //  time integrator
    const Teuchos::ParameterList& sdyn = Global::Problem::instance()->structural_dynamic_params();
    Teuchos::RCP<Adapter::StructureBaseAlgorithmNew> adapterbase_ptr =
        Adapter::build_structure_algorithm(sdyn);
    adapterbase_ptr->init(Global::Problem::instance()->tsi_dynamic_params(),
        const_cast<Teuchos::ParameterList&>(sdyn), structdis);

    // set the temperature; Monolithic does this in it's own constructor with potentially
    // redistributed discretizations
    if (Core::UTILS::IntegralValue<Inpar::TSI::SolutionSchemeOverFields>(
            Global::Problem::instance()->tsi_dynamic_params(), "COUPALGO") !=
        Inpar::TSI::Monolithic)
    {
      if (matchinggrid_)
        structdis->set_state(1, "temperature", thermo_field()->tempnp());
      else
        structdis->set_state(
            1, "temperature", volcoupl_->apply_vector_mapping12(thermo_field()->tempnp()));
    }

    adapterbase_ptr->setup();
    structure_ =
        Teuchos::rcp_dynamic_cast<Adapter::StructureWrapper>(adapterbase_ptr->structure_field());

    if (restart && Core::UTILS::IntegralValue<Inpar::TSI::SolutionSchemeOverFields>(
                       Global::Problem::instance()->tsi_dynamic_params(), "COUPALGO") ==
                       Inpar::TSI::Monolithic)
      structure_->setup();

    structure_field()->discretization()->clear_state(true);
  }

  // initialise displacement field needed for output()
  // (get noderowmap of discretisation for creating this multivector)
  // TODO: why nds 0 and not 1????
  dispnp_ = Teuchos::rcp(
      new Epetra_MultiVector(*(thermo_field()->discretization()->node_row_map()), 3, true));
  tempnp_ = Teuchos::rcp(
      new Epetra_MultiVector(*(structure_field()->discretization()->node_row_map()), 1, true));

  // setup coupling object for matching discretization
  if (matchinggrid_)
  {
    coupST_ = Teuchos::rcp(new Core::Adapter::Coupling());
    coupST_->setup_coupling(*structure_field()->discretization(), *thermo_field()->discretization(),
        *structure_field()->discretization()->node_row_map(),
        *thermo_field()->discretization()->node_row_map(), 1, true);
  }

  // setup mortar coupling
  if (Global::Problem::instance()->get_problem_type() == Core::ProblemType::tsi)
  {
    Core::Conditions::Condition* mrtrcond =
        structure_field()->discretization()->get_condition("MortarMulti");
    if (mrtrcond != nullptr)
    {
      mortar_coupling_ = Teuchos::rcp(new Mortar::MultiFieldCoupling());
      mortar_coupling_->push_back_coupling(
          structure_field()->discretization(), 0, std::vector<int>(3, 1));
      mortar_coupling_->push_back_coupling(
          thermo_field()->discretization(), 0, std::vector<int>(1, 1));
    }
  }

  // reset states
  structure_field()->discretization()->clear_state(true);
  thermo_field()->discretization()->clear_state(true);

  return;
}



/*----------------------------------------------------------------------*
 | update (protected)                                        dano 12/09 |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::update()
{
  apply_thermo_coupling_state(thermo_field()->tempnp());
  structure_field()->update();
  thermo_field()->update();
  if (contact_strategy_lagrange_ != Teuchos::null)
    contact_strategy_lagrange_->update((structure_field()->dispnp()));
  return;
}


/*----------------------------------------------------------------------*
 | output (protected)                                        dano 12/09 |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::output(bool forced_writerestart)
{
  // Note: The order in the output is important here!

  // In here control file entries are written. And these entries define the
  // order in which the filters handle the Discretizations, which in turn
  // defines the dof number ordering of the Discretizations.

  // call the TSI parameter list
  const Teuchos::ParameterList& tsidyn = Global::Problem::instance()->tsi_dynamic_params();
  // Get the parameters for the Newton iteration
  int upres = tsidyn.get<int>("RESULTSEVRY");
  int uprestart = tsidyn.get<int>("RESTARTEVRY");

  //========================
  // output for thermofield:
  //========================
  apply_struct_coupling_state(structure_field()->dispnp(), structure_field()->velnp());
  thermo_field()->output(forced_writerestart);

  // communicate the deformation to the thermal field,
  // current displacements are contained in Dispn()
  if (forced_writerestart == true and
      ((upres != 0 and (step() % upres == 0)) or ((uprestart != 0) and (step() % uprestart == 0))))
  {
    // displacement has already been written into thermo field for this step
  }
  else if ((upres != 0 and (step() % upres == 0)) or
           ((uprestart != 0) and (step() % uprestart == 0)) or forced_writerestart == true)
  {
    if (matchinggrid_)
    {
      output_deformation_in_thr(structure_field()->dispn(), structure_field()->discretization());

      thermo_field()->disc_writer()->write_vector("displacement", dispnp_, Core::IO::nodevector);
    }
    else
    {
      Teuchos::RCP<const Epetra_Vector> dummy =
          volcoupl_->apply_vector_mapping21(structure_field()->dispnp());

      // determine number of space dimensions
      const int numdim = Global::Problem::instance()->n_dim();

      int err(0);

      // loop over all local nodes of thermal discretisation
      for (int lnodeid = 0; lnodeid < (thermo_field()->discretization()->num_my_row_nodes());
           lnodeid++)
      {
        Core::Nodes::Node* thermnode = thermo_field()->discretization()->l_row_node(lnodeid);
        std::vector<int> thermnodedofs_1 = thermo_field()->discretization()->dof(1, thermnode);

        // now we transfer displacment dofs only
        for (int index = 0; index < numdim; ++index)
        {
          // global and processor's local fluid dof ID
          const int sgid = thermnodedofs_1[index];
          const int slid = thermo_field()->discretization()->dof_row_map(1)->LID(sgid);


          // get value of corresponding displacement component
          double disp = (*dummy)[slid];
          // insert velocity value into node-based vector
          err = dispnp_->ReplaceMyValue(lnodeid, index, disp);
          if (err != 0) FOUR_C_THROW("error while inserting a value into dispnp_");
        }

        // for security reasons in 1D or 2D problems:
        // set zeros for all unused velocity components
        for (int index = numdim; index < 3; ++index)
        {
          err = dispnp_->ReplaceMyValue(lnodeid, index, 0.0);
          if (err != 0) FOUR_C_THROW("error while inserting a value into dispnp_");
        }
      }  // for lnodid

      thermo_field()->disc_writer()->write_vector("displacement", dispnp_, Core::IO::nodevector);
    }
  }


  //===========================
  // output for structurefield:
  //===========================
  apply_thermo_coupling_state(thermo_field()->tempnp());
  structure_field()->output(forced_writerestart);

  // mapped temperatures for structure field
  if ((upres != 0 and (step() % upres == 0)) or ((uprestart != 0) and (step() % uprestart == 0)) or
      forced_writerestart == true)
    if (not matchinggrid_)
    {
      //************************************************************************************
      Teuchos::RCP<const Epetra_Vector> dummy1 =
          volcoupl_->apply_vector_mapping12(thermo_field()->tempnp());

      // loop over all local nodes of thermal discretisation
      for (int lnodeid = 0; lnodeid < (structure_field()->discretization()->num_my_row_nodes());
           lnodeid++)
      {
        Core::Nodes::Node* structnode = structure_field()->discretization()->l_row_node(lnodeid);
        std::vector<int> structdofs = structure_field()->discretization()->dof(1, structnode);

        // global and processor's local structure dof ID
        const int sgid = structdofs[0];
        const int slid = structure_field()->discretization()->dof_row_map(1)->LID(sgid);

        // get value of corresponding displacement component
        double temp = (*dummy1)[slid];
        // insert velocity value into node-based vector
        int err = tempnp_->ReplaceMyValue(lnodeid, 0, temp);
        if (err != 0) FOUR_C_THROW("error while inserting a value into tempnp_");
      }  // for lnodid

      structure_field()->discretization()->writer()->write_vector(
          "struct_temperature", tempnp_, Core::IO::nodevector);
    }


  // reset states
  structure_field()->discretization()->clear_state(true);
  thermo_field()->discretization()->clear_state(true);
}  // output()


/*----------------------------------------------------------------------*
 | communicate the displacement vector to THR field          dano 12/11 |
 | enable visualisation of thermal variables on deformed body           |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::output_deformation_in_thr(
    Teuchos::RCP<const Epetra_Vector> dispnp, Teuchos::RCP<Core::FE::Discretization> structdis)
{
  if (dispnp == Teuchos::null) FOUR_C_THROW("Got null pointer for displacements");

  int err(0);

  // get dofrowmap of structural discretisation
  const Epetra_Map* structdofrowmap = structdis->dof_row_map(0);

  // loop over all local nodes of thermal discretisation
  for (int lnodeid = 0; lnodeid < (thermo_field()->discretization()->num_my_row_nodes()); lnodeid++)
  {
    // Here we rely on the fact that the thermal discretisation is a clone of
    // the structural mesh.
    // => a thermal node has the same local (and global) ID as its corresponding
    // structural node!

    // get the processor's local structural node with the same lnodeid
    Core::Nodes::Node* structlnode = structdis->l_row_node(lnodeid);
    // get the degrees of freedom associated with this structural node
    std::vector<int> structnodedofs = structdis->dof(0, structlnode);
    // determine number of space dimensions
    const int numdim = Global::Problem::instance()->n_dim();

    // now we transfer displacment dofs only
    for (int index = 0; index < numdim; ++index)
    {
      // global and processor's local fluid dof ID
      const int sgid = structnodedofs[index];
      const int slid = structdofrowmap->LID(sgid);

      // get value of corresponding displacement component
      double disp = (*dispnp)[slid];
      // insert velocity value into node-based vector
      err = dispnp_->ReplaceMyValue(lnodeid, index, disp);
      if (err != 0) FOUR_C_THROW("error while inserting a value into dispnp_");
    }

    // for security reasons in 1D or 2D problems:
    // set zeros for all unused velocity components
    for (int index = numdim; index < 3; ++index)
    {
      err = dispnp_->ReplaceMyValue(lnodeid, index, 0.0);
      if (err != 0) FOUR_C_THROW("error while inserting a value into dispnp_");
    }

  }  // for lnodid

  return;

}  // output_deformation_in_thr()


/*----------------------------------------------------------------------*
 | calculate velocities                                      dano 12/10 |
 | like interface_velocity(disp) in FSI::DirichletNeumann                |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> TSI::Algorithm::calc_velocity(
    Teuchos::RCP<const Epetra_Vector> dispnp)
{
  Teuchos::RCP<Epetra_Vector> vel = Teuchos::null;
  // copy D_n onto V_n+1
  vel = Teuchos::rcp(new Epetra_Vector(*(structure_field()->dispn())));
  // calculate velocity with timestep Dt()
  //  V_n+1^k = (D_n+1^k - D_n) / Dt
  vel->Update(1. / dt(), *dispnp, -1. / dt());

  return vel;
}  // calc_velocity()


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::apply_thermo_coupling_state(
    Teuchos::RCP<const Epetra_Vector> temp, Teuchos::RCP<const Epetra_Vector> temp_res)
{
  if (matchinggrid_)
  {
    if (temp != Teuchos::null)
      structure_field()->discretization()->set_state(1, "temperature", temp);
    if (temp_res != Teuchos::null)
      structure_field()->discretization()->set_state(1, "residual temperature", temp_res);
  }
  else
  {
    if (temp != Teuchos::null)
      structure_field()->discretization()->set_state(
          1, "temperature", volcoupl_->apply_vector_mapping12(temp));
  }

  // set new temperatures to contact
  {
    if (contact_strategy_lagrange_ != Teuchos::null)
      contact_strategy_lagrange_->set_state(
          Mortar::state_temperature, *coupST_()->slave_to_master(thermo_field()->tempnp()));
    if (contact_strategy_nitsche_ != Teuchos::null)
      contact_strategy_nitsche_->set_state(Mortar::state_temperature, *thermo_field()->tempnp());
  }
}  // apply_thermo_coupling_state()


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::apply_struct_coupling_state(
    Teuchos::RCP<const Epetra_Vector> disp, Teuchos::RCP<const Epetra_Vector> vel)
{
  if (matchinggrid_)
  {
    if (disp != Teuchos::null) thermo_field()->discretization()->set_state(1, "displacement", disp);
    if (vel != Teuchos::null) thermo_field()->discretization()->set_state(1, "velocity", vel);
  }
  else
  {
    if (disp != Teuchos::null)
      thermo_field()->discretization()->set_state(
          1, "displacement", volcoupl_->apply_vector_mapping21(disp));
    if (vel != Teuchos::null)
      thermo_field()->discretization()->set_state(
          1, "velocity", volcoupl_->apply_vector_mapping21(vel));
  }
}  // apply_struct_coupling_state()


/*----------------------------------------------------------------------*/
void TSI::Algorithm::prepare_contact_strategy()
{
  Inpar::CONTACT::SolvingStrategy stype =
      Core::UTILS::IntegralValue<Inpar::CONTACT::SolvingStrategy>(
          Global::Problem::instance()->contact_dynamic_params(), "STRATEGY");

  if (stype == Inpar::CONTACT::solution_nitsche)
  {
    if (Core::UTILS::IntegralValue<Inpar::Solid::IntegrationStrategy>(
            Global::Problem::instance()->structural_dynamic_params(), "INT_STRATEGY") !=
        Inpar::Solid::int_standard)
      FOUR_C_THROW("thermo-mechanical contact only with new structural time integration");

    if (coupST_ == Teuchos::null) FOUR_C_THROW("coupST_ not yet here");

    Solid::MODELEVALUATOR::Contact& a = static_cast<Solid::MODELEVALUATOR::Contact&>(
        structure_field()->model_evaluator(Inpar::Solid::model_contact));
    contact_strategy_nitsche_ =
        Teuchos::rcp_dynamic_cast<CONTACT::NitscheStrategyTsi>(a.strategy_ptr(), false);
    contact_strategy_nitsche_->enable_redistribution();

    thermo_->set_nitsche_contact_strategy(contact_strategy_nitsche_);

    return;
  }

  else if (stype == Inpar::CONTACT::solution_lagmult)
  {
    if (structure_field()->have_model(Inpar::Solid::model_contact))
      FOUR_C_THROW(
          "structure should not have a Lagrange strategy ... as long as condensed"
          "contact formulations are not moved to the new structural time integration");

    std::vector<Core::Conditions::Condition*> ccond(0);
    structure_field()->discretization()->get_condition("Contact", ccond);
    if (ccond.size() == 0) return;

    // ---------------------------------------------------------------------
    // create the contact factory
    // ---------------------------------------------------------------------
    CONTACT::STRATEGY::Factory factory;
    factory.init(structure_field()->discretization());
    factory.setup();

    // check the problem dimension
    factory.check_dimension();

    // create some local variables (later to be stored in strategy)
    std::vector<Teuchos::RCP<CONTACT::Interface>> interfaces;
    Teuchos::ParameterList cparams;

    // read and check contact input parameters
    factory.read_and_check_input(cparams);

    // ---------------------------------------------------------------------
    // build the contact interfaces
    // ---------------------------------------------------------------------
    // FixMe Would be great, if we get rid of these poro parameters...
    bool poroslave = false;
    bool poromaster = false;
    factory.build_interfaces(cparams, interfaces, poroslave, poromaster);

    // ---------------------------------------------------------------------
    // build the solver strategy object
    // ---------------------------------------------------------------------
    contact_strategy_lagrange_ = Teuchos::rcp_dynamic_cast<CONTACT::LagrangeStrategyTsi>(
        factory.build_strategy(cparams, poroslave, poromaster, 1e8, interfaces), true);

    // build the search tree
    factory.build_search_tree(interfaces);

    // print final screen output
    factory.print(interfaces, contact_strategy_lagrange_, cparams);

    // ---------------------------------------------------------------------
    // final touches to the contact strategy
    // ---------------------------------------------------------------------

    contact_strategy_lagrange_->store_dirichlet_status(structure_field()->get_dbc_map_extractor());

    Teuchos::RCP<Epetra_Vector> zero_disp =
        Teuchos::rcp(new Epetra_Vector(*structure_field()->dof_row_map(), true));
    contact_strategy_lagrange_->set_state(Mortar::state_new_displacement, *zero_disp);
    contact_strategy_lagrange_->save_reference_state(zero_disp);
    contact_strategy_lagrange_->evaluate_reference_state();
    contact_strategy_lagrange_->inttime_init();
    contact_strategy_lagrange_->set_time_integration_info(structure_field()->tim_int_param(),
        Core::UTILS::IntegralValue<Inpar::Solid::DynamicType>(
            Global::Problem::instance()->structural_dynamic_params(), "DYNAMICTYP"));
    contact_strategy_lagrange_->redistribute_contact(
        structure_field()->dispn(), structure_field()->veln());

    if (contact_strategy_lagrange_ != Teuchos::null)
    {
      contact_strategy_lagrange_->set_alphaf_thermo(
          Global::Problem::instance()->thermal_dynamic_params());
      contact_strategy_lagrange_->set_coupling(coupST_);
    }
  }
}

FOUR_C_NAMESPACE_CLOSE
