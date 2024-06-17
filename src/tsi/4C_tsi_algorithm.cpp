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
#include "4C_structure_new_model_evaluator_contact.hpp"
#include "4C_structure_new_model_evaluator_structure.hpp"
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
    : AlgorithmBase(comm, Global::Problem::Instance()->TSIDynamicParams()),
      dispnp_(Teuchos::null),
      tempnp_(Teuchos::null),
      matchinggrid_(Core::UTILS::IntegralValue<bool>(
          Global::Problem::Instance()->TSIDynamicParams(), "MATCHINGGRID")),
      volcoupl_(Teuchos::null)
{
  // access the structural discretization
  Teuchos::RCP<Core::FE::Discretization> structdis =
      Global::Problem::Instance()->GetDis("structure");
  // access the thermo discretization
  Teuchos::RCP<Core::FE::Discretization> thermodis = Global::Problem::Instance()->GetDis("thermo");

  // get the problem instance
  Global::Problem* problem = Global::Problem::Instance();
  // get the restart step
  const int restart = problem->Restart();

  if (!matchinggrid_)
  {
    // Scheme: non matching meshes --> volumetric mortar coupling...
    volcoupl_ = Teuchos::rcp(new Core::Adapter::MortarVolCoupl());

    Teuchos::RCP<Core::VolMortar::UTILS::DefaultMaterialStrategy> materialstrategy =
        Teuchos::rcp(new TSI::UTILS::TSIMaterialStrategy());
    // init coupling adapter projection matrices
    volcoupl_->Init(Global::Problem::Instance()->NDim(), structdis, thermodis, nullptr, nullptr,
        nullptr, nullptr, materialstrategy);
    // redistribute discretizations to meet needs of volmortar coupling
    Teuchos::ParameterList binning_params = Global::Problem::Instance()->binning_strategy_params();
    Core::UTILS::AddEnumClassToParameterList<Core::FE::ShapeFunctionType>(
        "spatial_approximation_type", Global::Problem::Instance()->spatial_approximation_type(),
        binning_params);
    volcoupl_->Redistribute(binning_params, Global::Problem::Instance()->OutputControlFile());
    // setup projection matrices
    volcoupl_->Setup(Global::Problem::Instance()->VolmortarParams());
  }

  if (Core::UTILS::IntegralValue<Inpar::STR::IntegrationStrategy>(
          Global::Problem::Instance()->structural_dynamic_params(), "INT_STRATEGY") ==
      Inpar::STR::int_old)
    FOUR_C_THROW("old structural time integration no longer supported in tsi");
  else
  {
    Teuchos::RCP<Adapter::ThermoBaseAlgorithm> thermo =
        Teuchos::rcp(new Adapter::ThermoBaseAlgorithm(
            Global::Problem::Instance()->TSIDynamicParams(), thermodis));
    thermo_ = thermo->ThermoFieldrcp();

    //  // access structural dynamic params list which will be possibly modified while creating the
    //  time integrator
    const Teuchos::ParameterList& sdyn = Global::Problem::Instance()->structural_dynamic_params();
    Teuchos::RCP<Adapter::StructureBaseAlgorithmNew> adapterbase_ptr =
        Adapter::build_structure_algorithm(sdyn);
    adapterbase_ptr->Init(Global::Problem::Instance()->TSIDynamicParams(),
        const_cast<Teuchos::ParameterList&>(sdyn), structdis);

    // set the temperature; Monolithic does this in it's own constructor with potentially
    // redistributed discretizations
    if (Core::UTILS::IntegralValue<Inpar::TSI::SolutionSchemeOverFields>(
            Global::Problem::Instance()->TSIDynamicParams(), "COUPALGO") != Inpar::TSI::Monolithic)
    {
      if (matchinggrid_)
        structdis->set_state(1, "temperature", ThermoField()->Tempnp());
      else
        structdis->set_state(
            1, "temperature", volcoupl_->apply_vector_mapping12(ThermoField()->Tempnp()));
    }

    adapterbase_ptr->Setup();
    structure_ =
        Teuchos::rcp_dynamic_cast<Adapter::StructureWrapper>(adapterbase_ptr->structure_field());

    if (restart &&
        Core::UTILS::IntegralValue<Inpar::TSI::SolutionSchemeOverFields>(
            Global::Problem::Instance()->TSIDynamicParams(), "COUPALGO") == Inpar::TSI::Monolithic)
      structure_->Setup();

    structure_field()->discretization()->ClearState(true);
  }

  // initialise displacement field needed for Output()
  // (get noderowmap of discretisation for creating this multivector)
  // TODO: why nds 0 and not 1????
  dispnp_ = Teuchos::rcp(
      new Epetra_MultiVector(*(ThermoField()->discretization()->NodeRowMap()), 3, true));
  tempnp_ = Teuchos::rcp(
      new Epetra_MultiVector(*(structure_field()->discretization()->NodeRowMap()), 1, true));

  // setup coupling object for matching discretization
  if (matchinggrid_)
  {
    coupST_ = Teuchos::rcp(new Core::Adapter::Coupling());
    coupST_->setup_coupling(*structure_field()->discretization(), *ThermoField()->discretization(),
        *structure_field()->discretization()->NodeRowMap(),
        *ThermoField()->discretization()->NodeRowMap(), 1, true);
  }

  // setup mortar coupling
  if (Global::Problem::Instance()->GetProblemType() == Core::ProblemType::tsi)
  {
    Core::Conditions::Condition* mrtrcond =
        structure_field()->discretization()->GetCondition("MortarMulti");
    if (mrtrcond != nullptr)
    {
      mortar_coupling_ = Teuchos::rcp(new Mortar::MultiFieldCoupling());
      mortar_coupling_->PushBackCoupling(
          structure_field()->discretization(), 0, std::vector<int>(3, 1));
      mortar_coupling_->PushBackCoupling(
          ThermoField()->discretization(), 0, std::vector<int>(1, 1));
    }
  }

  // reset states
  structure_field()->discretization()->ClearState(true);
  ThermoField()->discretization()->ClearState(true);

  return;
}



/*----------------------------------------------------------------------*
 | update (protected)                                        dano 12/09 |
 *----------------------------------------------------------------------*/
void TSI::Algorithm::update()
{
  apply_thermo_coupling_state(ThermoField()->Tempnp());
  structure_field()->Update();
  ThermoField()->Update();
  if (contact_strategy_lagrange_ != Teuchos::null)
    contact_strategy_lagrange_->Update((structure_field()->Dispnp()));
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
  const Teuchos::ParameterList& tsidyn = Global::Problem::Instance()->TSIDynamicParams();
  // Get the parameters for the Newton iteration
  int upres = tsidyn.get<int>("RESULTSEVRY");
  int uprestart = tsidyn.get<int>("RESTARTEVRY");

  //========================
  // output for thermofield:
  //========================
  apply_struct_coupling_state(structure_field()->Dispnp(), structure_field()->Velnp());
  ThermoField()->Output(forced_writerestart);

  // communicate the deformation to the thermal field,
  // current displacements are contained in Dispn()
  if (forced_writerestart == true and
      ((upres != 0 and (Step() % upres == 0)) or ((uprestart != 0) and (Step() % uprestart == 0))))
  {
    // displacement has already been written into thermo field for this step
  }
  else if ((upres != 0 and (Step() % upres == 0)) or
           ((uprestart != 0) and (Step() % uprestart == 0)) or forced_writerestart == true)
  {
    if (matchinggrid_)
    {
      output_deformation_in_thr(structure_field()->Dispn(), structure_field()->discretization());

      ThermoField()->DiscWriter()->write_vector("displacement", dispnp_, Core::IO::nodevector);
    }
    else
    {
      Teuchos::RCP<const Epetra_Vector> dummy =
          volcoupl_->apply_vector_mapping21(structure_field()->Dispnp());

      // determine number of space dimensions
      const int numdim = Global::Problem::Instance()->NDim();

      int err(0);

      // loop over all local nodes of thermal discretisation
      for (int lnodeid = 0; lnodeid < (ThermoField()->discretization()->NumMyRowNodes()); lnodeid++)
      {
        Core::Nodes::Node* thermnode = ThermoField()->discretization()->lRowNode(lnodeid);
        std::vector<int> thermnodedofs_1 = ThermoField()->discretization()->Dof(1, thermnode);

        // now we transfer displacment dofs only
        for (int index = 0; index < numdim; ++index)
        {
          // global and processor's local fluid dof ID
          const int sgid = thermnodedofs_1[index];
          const int slid = ThermoField()->discretization()->dof_row_map(1)->LID(sgid);


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

      ThermoField()->DiscWriter()->write_vector("displacement", dispnp_, Core::IO::nodevector);
    }
  }


  //===========================
  // output for structurefield:
  //===========================
  apply_thermo_coupling_state(ThermoField()->Tempnp());
  structure_field()->Output(forced_writerestart);

  // mapped temperatures for structure field
  if ((upres != 0 and (Step() % upres == 0)) or ((uprestart != 0) and (Step() % uprestart == 0)) or
      forced_writerestart == true)
    if (not matchinggrid_)
    {
      //************************************************************************************
      Teuchos::RCP<const Epetra_Vector> dummy1 =
          volcoupl_->apply_vector_mapping12(ThermoField()->Tempnp());

      // loop over all local nodes of thermal discretisation
      for (int lnodeid = 0; lnodeid < (structure_field()->discretization()->NumMyRowNodes());
           lnodeid++)
      {
        Core::Nodes::Node* structnode = structure_field()->discretization()->lRowNode(lnodeid);
        std::vector<int> structdofs = structure_field()->discretization()->Dof(1, structnode);

        // global and processor's local structure dof ID
        const int sgid = structdofs[0];
        const int slid = structure_field()->discretization()->dof_row_map(1)->LID(sgid);

        // get value of corresponding displacement component
        double temp = (*dummy1)[slid];
        // insert velocity value into node-based vector
        int err = tempnp_->ReplaceMyValue(lnodeid, 0, temp);
        if (err != 0) FOUR_C_THROW("error while inserting a value into tempnp_");
      }  // for lnodid

      structure_field()->discretization()->Writer()->write_vector(
          "struct_temperature", tempnp_, Core::IO::nodevector);
    }


  // reset states
  structure_field()->discretization()->ClearState(true);
  ThermoField()->discretization()->ClearState(true);
}  // Output()


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
  for (int lnodeid = 0; lnodeid < (ThermoField()->discretization()->NumMyRowNodes()); lnodeid++)
  {
    // Here we rely on the fact that the thermal discretisation is a clone of
    // the structural mesh.
    // => a thermal node has the same local (and global) ID as its corresponding
    // structural node!

    // get the processor's local structural node with the same lnodeid
    Core::Nodes::Node* structlnode = structdis->lRowNode(lnodeid);
    // get the degrees of freedom associated with this structural node
    std::vector<int> structnodedofs = structdis->Dof(0, structlnode);
    // determine number of space dimensions
    const int numdim = Global::Problem::Instance()->NDim();

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
  vel = Teuchos::rcp(new Epetra_Vector(*(structure_field()->Dispn())));
  // calculate velocity with timestep Dt()
  //  V_n+1^k = (D_n+1^k - D_n) / Dt
  vel->Update(1. / Dt(), *dispnp, -1. / Dt());

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
          Mortar::state_temperature, *coupST_()->SlaveToMaster(ThermoField()->Tempnp()));
    if (contact_strategy_nitsche_ != Teuchos::null)
      contact_strategy_nitsche_->set_state(Mortar::state_temperature, *ThermoField()->Tempnp());
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
    if (disp != Teuchos::null) ThermoField()->discretization()->set_state(1, "displacement", disp);
    if (vel != Teuchos::null) ThermoField()->discretization()->set_state(1, "velocity", vel);
  }
  else
  {
    if (disp != Teuchos::null)
      ThermoField()->discretization()->set_state(
          1, "displacement", volcoupl_->apply_vector_mapping21(disp));
    if (vel != Teuchos::null)
      ThermoField()->discretization()->set_state(
          1, "velocity", volcoupl_->apply_vector_mapping21(vel));
  }
}  // apply_struct_coupling_state()


/*----------------------------------------------------------------------*/
void TSI::Algorithm::prepare_contact_strategy()
{
  Inpar::CONTACT::SolvingStrategy stype =
      Core::UTILS::IntegralValue<Inpar::CONTACT::SolvingStrategy>(
          Global::Problem::Instance()->contact_dynamic_params(), "STRATEGY");

  if (stype == Inpar::CONTACT::solution_nitsche)
  {
    if (Core::UTILS::IntegralValue<Inpar::STR::IntegrationStrategy>(
            Global::Problem::Instance()->structural_dynamic_params(), "INT_STRATEGY") !=
        Inpar::STR::int_standard)
      FOUR_C_THROW("thermo-mechanical contact only with new structural time integration");

    if (coupST_ == Teuchos::null) FOUR_C_THROW("coupST_ not yet here");

    STR::MODELEVALUATOR::Contact& a = static_cast<STR::MODELEVALUATOR::Contact&>(
        structure_field()->ModelEvaluator(Inpar::STR::model_contact));
    contact_strategy_nitsche_ =
        Teuchos::rcp_dynamic_cast<CONTACT::NitscheStrategyTsi>(a.strategy_ptr(), false);
    contact_strategy_nitsche_->enable_redistribution();

    thermo_->set_nitsche_contact_strategy(contact_strategy_nitsche_);

    return;
  }

  else if (stype == Inpar::CONTACT::solution_lagmult)
  {
    if (structure_field()->HaveModel(Inpar::STR::model_contact))
      FOUR_C_THROW(
          "structure should not have a Lagrange strategy ... as long as condensed"
          "contact formulations are not moved to the new structural time integration");

    std::vector<Core::Conditions::Condition*> ccond(0);
    structure_field()->discretization()->GetCondition("Contact", ccond);
    if (ccond.size() == 0) return;

    // ---------------------------------------------------------------------
    // create the contact factory
    // ---------------------------------------------------------------------
    CONTACT::STRATEGY::Factory factory;
    factory.Init(structure_field()->discretization());
    factory.Setup();

    // check the problem dimension
    factory.CheckDimension();

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
    factory.BuildInterfaces(cparams, interfaces, poroslave, poromaster);

    // ---------------------------------------------------------------------
    // build the solver strategy object
    // ---------------------------------------------------------------------
    contact_strategy_lagrange_ = Teuchos::rcp_dynamic_cast<CONTACT::LagrangeStrategyTsi>(
        factory.BuildStrategy(cparams, poroslave, poromaster, 1e8, interfaces), true);

    // build the search tree
    factory.BuildSearchTree(interfaces);

    // print final screen output
    factory.Print(interfaces, contact_strategy_lagrange_, cparams);

    // ---------------------------------------------------------------------
    // final touches to the contact strategy
    // ---------------------------------------------------------------------

    contact_strategy_lagrange_->store_dirichlet_status(structure_field()->GetDBCMapExtractor());

    Teuchos::RCP<Epetra_Vector> zero_disp =
        Teuchos::rcp(new Epetra_Vector(*structure_field()->dof_row_map(), true));
    contact_strategy_lagrange_->set_state(Mortar::state_new_displacement, *zero_disp);
    contact_strategy_lagrange_->SaveReferenceState(zero_disp);
    contact_strategy_lagrange_->evaluate_reference_state();
    contact_strategy_lagrange_->Inttime_init();
    contact_strategy_lagrange_->set_time_integration_info(structure_field()->TimIntParam(),
        Core::UTILS::IntegralValue<Inpar::STR::DynamicType>(
            Global::Problem::Instance()->structural_dynamic_params(), "DYNAMICTYP"));
    contact_strategy_lagrange_->RedistributeContact(
        structure_field()->Dispn(), structure_field()->Veln());

    if (contact_strategy_lagrange_ != Teuchos::null)
    {
      contact_strategy_lagrange_->SetAlphafThermo(
          Global::Problem::Instance()->thermal_dynamic_params());
      contact_strategy_lagrange_->SetCoupling(coupST_);
    }
  }
}

FOUR_C_NAMESPACE_CLOSE
