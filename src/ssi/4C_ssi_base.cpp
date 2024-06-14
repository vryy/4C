/*----------------------------------------------------------------------*/
/*! \file
 \brief base class for all scalar structure algorithms

 \level 1


 *------------------------------------------------------------------------------------------------*/

#include "4C_ssi_base.hpp"

#include "4C_adapter_scatra_base_algorithm.hpp"
#include "4C_adapter_str_factory.hpp"
#include "4C_adapter_str_ssiwrapper.hpp"
#include "4C_adapter_str_structure_new.hpp"
#include "4C_contact_nitsche_strategy_ssi.hpp"
#include "4C_coupling_volmortar.hpp"
#include "4C_fem_general_utils_createdis.hpp"
#include "4C_global_data.hpp"
#include "4C_global_data_read.hpp"
#include "4C_inpar_ssi.hpp"
#include "4C_io_control.hpp"
#include "4C_io_inputreader.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_membrane.hpp"
#include "4C_rebalance_binning_based.hpp"
#include "4C_scatra_timint_implicit.hpp"
#include "4C_scatra_timint_meshtying_strategy_s2i.hpp"
#include "4C_scatra_utils.hpp"
#include "4C_so3_base.hpp"
#include "4C_solid_3D_ele.hpp"
#include "4C_ssi_clonestrategy.hpp"
#include "4C_ssi_coupling.hpp"
#include "4C_ssi_partitioned.hpp"
#include "4C_ssi_resulttest.hpp"
#include "4C_ssi_str_model_evaluator_partitioned.hpp"
#include "4C_ssi_utils.hpp"
#include "4C_structure_new_model_evaluator_contact.hpp"
#include "4C_utils_function_of_time.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
SSI::SSIBase::SSIBase(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams)
    : AlgorithmBase(comm, globaltimeparams),
      diff_time_step_size_(
          static_cast<int>(Core::UTILS::IntegralValue<int>(globaltimeparams, "DIFFTIMESTEPSIZE"))),
      fieldcoupling_(Teuchos::getIntegralValue<Inpar::SSI::FieldCoupling>(
          Global::Problem::Instance()->SSIControlParams(), "FIELDCOUPLING")),
      is_scatra_manifold_(
          Core::UTILS::IntegralValue<bool>(globaltimeparams.sublist("MANIFOLD"), "ADD_MANIFOLD")),
      is_manifold_meshtying_(Core::UTILS::IntegralValue<bool>(
          globaltimeparams.sublist("MANIFOLD"), "MESHTYING_MANIFOLD")),
      is_s2i_kinetic_with_pseudo_contact_(
          check_s2_i_kinetics_condition_for_pseudo_contact("structure")),
      macro_scale_(Global::Problem::Instance()->Materials()->FirstIdByType(
                       Core::Materials::m_scatra_multiscale) != -1 or
                   Global::Problem::Instance()->Materials()->FirstIdByType(
                       Core::Materials::m_newman_multiscale) != -1),
      ssiinterfacecontact_(
          Global::Problem::Instance()->GetDis("structure")->GetCondition("SSIInterfaceContact") !=
          nullptr),
      ssiinterfacemeshtying_(Global::Problem::Instance()
                                 ->GetDis("structure")
                                 ->GetCondition("ssi_interface_meshtying") != nullptr),
      temperature_funct_num_(
          Global::Problem::Instance()->ELCHControlParams().get<int>("TEMPERATURE_FROM_FUNCT")),
      use_old_structure_(Global::Problem::Instance()->structural_dynamic_params().get<std::string>(
                             "INT_STRATEGY") == "Old")
{
  // Keep this constructor empty!
  // First do everything on the more basic objects like the discretizations, like e.g.
  // redistribution of elements. Only then call the setup to this class. This will call the setup to
  // all classes in the inheritance hierarchy. This way, this class may also override a method that
  // is called during setup() in a base class.
}

/*----------------------------------------------------------------------*
 | Init this class                                          rauch 08/16 |
 *----------------------------------------------------------------------*/
void SSI::SSIBase::init(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& scatraparams, const Teuchos::ParameterList& structparams,
    const std::string& struct_disname, const std::string& scatra_disname, bool isAle)
{
  // reset the setup flag
  set_is_setup(false);

  // do discretization specific setup (e.g. clone discr. scatra from structure)
  InitDiscretizations(comm, struct_disname, scatra_disname,
      Core::UTILS::IntegralValue<bool>(globaltimeparams, "REDISTRIBUTE_SOLID"));

  init_time_integrators(
      globaltimeparams, scatraparams, structparams, struct_disname, scatra_disname, isAle);

  const RedistributionType redistribution_type = InitFieldCoupling(struct_disname);

  if (redistribution_type != SSI::RedistributionType::none) Redistribute(redistribution_type);

  check_ssi_flags();

  check_ssi_interface_conditions(struct_disname);

  // set isinit_ flag true
  set_is_init(true);
}

/*----------------------------------------------------------------------*
 | Setup this class                                         rauch 08/16 |
 *----------------------------------------------------------------------*/
void SSI::SSIBase::setup()
{
  // check initialization
  check_is_init();

  // set up helper class for field coupling
  ssicoupling_->setup();

  // in case of an ssi  multi scale formulation we need to set the displacement here
  auto dummy_vec = Teuchos::rcp(
      new Epetra_Vector(*Global::Problem::Instance()->GetDis("structure")->dof_row_map(), true));
  ssicoupling_->set_mesh_disp(ScaTraBaseAlgorithm(), dummy_vec);

  // set up scalar transport field
  ScaTraField()->setup();
  if (is_sca_tra_manifold()) ScaTraManifold()->setup();

  // only relevant for new structural time integration
  // only if adapter base has not already been set up outside
  if (not use_old_structure_ and not struct_adapterbase_ptr_->is_setup())
  {
    // set up structural model evaluator
    setup_model_evaluator();

    // pass initial scalar field to structural discretization to correctly compute initial
    // accelerations
    if (Teuchos::getIntegralValue<Inpar::SSI::SolutionSchemeOverFields>(
            Global::Problem::Instance()->SSIControlParams(), "COUPALGO") !=
        Inpar::SSI::SolutionSchemeOverFields::ssi_OneWay_SolidToScatra)
      ssicoupling_->SetScalarField(
          *Global::Problem::Instance()->GetDis("structure"), ScaTraField()->Phinp(), 1);

    if (macro_scale_)
    {
      ScaTraField()->calc_mean_micro_concentration();
      ssicoupling_->SetScalarFieldMicro(
          *Global::Problem::Instance()->GetDis("structure"), ScaTraField()->PhinpMicro(), 2);
    }

    //   temperature is non primary variable. Only set, if function for temperature is given
    if (temperature_funct_num_ != -1)
    {
      temperature_vector_ = Teuchos::rcp(new Epetra_Vector(
          *Global::Problem::Instance()->GetDis("structure")->dof_row_map(2), true));

      temperature_vector_->PutScalar(
          Global::Problem::Instance()
              ->FunctionById<Core::UTILS::FunctionOfTime>(temperature_funct_num_ - 1)
              .evaluate(Time()));

      ssicoupling_->SetTemperatureField(
          *Global::Problem::Instance()->GetDis("structure"), temperature_vector_);
    }

    // set up structural base algorithm
    struct_adapterbase_ptr_->setup();

    // get wrapper and cast it to specific type
    // do not do so, in case the wrapper has already been set from outside
    if (structure_ == Teuchos::null)
      structure_ = Teuchos::rcp_dynamic_cast<Adapter::SSIStructureWrapper>(
          struct_adapterbase_ptr_->structure_field());

    if (structure_ == Teuchos::null)
    {
      FOUR_C_THROW(
          "No valid pointer to Adapter::SSIStructureWrapper !\n"
          "Either cast failed, or no valid wrapper was set using set_structure_wrapper(...) !");
    }
  }

  // for old structural time integration
  else if (use_old_structure_)
    structure_->setup();

  if (is_s2_i_kinetics_with_pseudo_contact())
  {
    auto dummy_stress_state =
        Teuchos::rcp(new Epetra_Vector(*structure_field()->discretization()->dof_row_map(2), true));
    ssicoupling_->set_mechanical_stress_state(*ScaTraField()->discretization(), dummy_stress_state,
        ScaTraField()->nds_two_tensor_quantity());
  }

  // check maps from scalar transport and structure discretizations
  if (ScaTraField()->dof_row_map()->NumGlobalElements() == 0)
    FOUR_C_THROW("Scalar transport discretization does not have any degrees of freedom!");
  if (structure_->dof_row_map()->NumGlobalElements() == 0)
    FOUR_C_THROW("Structure discretization does not have any degrees of freedom!");

  // set up materials
  ssicoupling_->assign_material_pointers(
      structure_->discretization(), ScaTraField()->discretization());

  // set up scatra-scatra interface coupling
  if (ssi_interface_meshtying())
  {
    ssi_structure_meshtying_ = Teuchos::rcp(new SSI::UTILS::SSIMeshTying(
        "ssi_interface_meshtying", structure_->discretization(), true, true));

    // extract meshtying strategy for scatra-scatra interface coupling on scatra discretization
    meshtying_strategy_s2i_ =
        Teuchos::rcp_dynamic_cast<const ScaTra::MeshtyingStrategyS2I>(ScaTraField()->Strategy());

    // safety checks
    if (meshtying_strategy_s2i_ == Teuchos::null)
      FOUR_C_THROW("Invalid scatra-scatra interface coupling strategy!");
  }

  // construct vector of zeroes
  zeros_structure_ = Core::LinAlg::CreateVector(*structure_->dof_row_map());

  // set flag
  set_is_setup(true);
}

/*----------------------------------------------------------------------*
 | Setup the discretizations                                rauch 08/16 |
 *----------------------------------------------------------------------*/
void SSI::SSIBase::InitDiscretizations(const Epetra_Comm& comm, const std::string& struct_disname,
    const std::string& scatra_disname, bool redistribute_struct_dis)
{
  Global::Problem* problem = Global::Problem::Instance();

  auto structdis = problem->GetDis(struct_disname);
  auto scatradis = problem->GetDis(scatra_disname);

  if (redistribute_struct_dis)
  {
    Teuchos::ParameterList binning_params = Global::Problem::Instance()->binning_strategy_params();
    Core::UTILS::AddEnumClassToParameterList<Core::FE::ShapeFunctionType>(
        "spatial_approximation_type", Global::Problem::Instance()->spatial_approximation_type(),
        binning_params);
    auto element_filter = [](const Core::Elements::Element* element)
    { return BINSTRATEGY::UTILS::SpecialElement::none; };
    auto rigid_sphere_radius = [](const Core::Elements::Element* element) { return 0.0; };
    auto correct_beam_center_node = [](const Core::Nodes::Node* node) { return node; };
    Core::Rebalance::RebalanceDiscretizationsByBinning(binning_params,
        Global::Problem::Instance()->OutputControlFile(), {structdis}, element_filter,
        rigid_sphere_radius, correct_beam_center_node, false);
  }

  if (scatradis->NumGlobalNodes() == 0)
  {
    if (fieldcoupling_ != Inpar::SSI::FieldCoupling::volume_match and
        fieldcoupling_ != Inpar::SSI::FieldCoupling::volumeboundary_match)
    {
      FOUR_C_THROW(
          "If 'FIELDCOUPLING' is NOT 'volume_matching' or 'volumeboundary_matching' in the SSI "
          "CONTROL section cloning of the scatra discretization from the structure discretization "
          "is not supported!");
    }

    // fill scatra discretization by cloning structure discretization
    Core::FE::CloneDiscretization<SSI::ScatraStructureCloneStrategy>(
        structdis, scatradis, Global::Problem::Instance()->CloningMaterialMap());
    scatradis->fill_complete();

    // create discretization for scatra manifold based on SSISurfaceManifold condition
    if (is_sca_tra_manifold())
    {
      auto scatra_manifold_dis = problem->GetDis("scatra_manifold");
      Core::FE::CloneDiscretizationFromCondition<SSI::ScatraStructureCloneStrategyManifold>(
          *structdis, *scatra_manifold_dis, "SSISurfaceManifold",
          Global::Problem::Instance()->CloningMaterialMap());

      // clone conditions. Needed this way, as many conditions are cloned from SSISurfaceManifold.
      std::vector<std::map<std::string, std::string>> conditions_to_copy = {
          {std::make_pair("SSISurfaceManifold", "SSISurfaceManifold")},
          {std::make_pair("ScaTraManifoldInitfield", "Initfield")},
          {std::make_pair("ManifoldDirichlet", "Dirichlet")}};

      //! in case of no mesh tying between manifolds: partition manifold domains
      if (!is_manifold_meshtying_)
      {
        std::map<std::string, std::string> temp_map = {
            std::make_pair("SSISurfaceManifold", "ScatraPartitioning")};
        conditions_to_copy.emplace_back(temp_map);
      }

      const auto output_scalar_type = Core::UTILS::IntegralValue<Inpar::ScaTra::OutputScalarType>(
          problem->scalar_transport_dynamic_params(), "OUTPUTSCALARS");
      if (output_scalar_type == Inpar::ScaTra::outputscalars_condition or
          output_scalar_type == Inpar::ScaTra::outputscalars_entiredomain_condition)
      {
        std::map<std::string, std::string> tempmap = {
            std::make_pair("SSISurfaceManifold", "TotalAndMeanScalar")};

        conditions_to_copy.emplace_back(tempmap);
      }

      Core::FE::DiscretizationCreatorBase creator;
      for (const auto& condition_to_copy : conditions_to_copy)
        creator.CopyConditions(*structdis, *scatra_manifold_dis, condition_to_copy);

      scatra_manifold_dis->fill_complete();

      //! in case of mesh tying between manifolds: unite manifold domains -> create new
      //! ScatraPartitioning condition
      if (is_manifold_meshtying_)
      {
        // create vector of all node GIDs (all procs) of manifold dis
        int num_my_nodes = scatra_manifold_dis->NodeRowMap()->NumMyElements();
        std::vector<int> my_node_ids(num_my_nodes);
        for (int lid = 0; lid < num_my_nodes; ++lid)
          my_node_ids[lid] = scatra_manifold_dis->NodeRowMap()->GID(lid);

        int max_num_nodes = 0;
        Comm().MaxAll(&num_my_nodes, &max_num_nodes, 1);

        // resize vector and fill with place holders (-1)
        my_node_ids.resize(max_num_nodes, -1);

        std::vector<int> glob_node_ids(max_num_nodes * Comm().NumProc(), -1);
        Comm().GatherAll(
            my_node_ids.data(), glob_node_ids.data(), static_cast<int>(my_node_ids.size()));

        // remove place holders (-1)
        glob_node_ids.erase(
            std::remove(glob_node_ids.begin(), glob_node_ids.end(), -1), glob_node_ids.end());

        // create new condition
        const int num_conditions = static_cast<int>(scatra_manifold_dis->GetAllConditions().size());
        auto cond = Teuchos::rcp(new Core::Conditions::Condition(num_conditions + 1,
            Core::Conditions::ScatraPartitioning, true, Core::Conditions::geometry_type_surface));
        cond->parameters().Add("ConditionID", 0);
        cond->SetNodes(glob_node_ids);

        scatra_manifold_dis->SetCondition("ScatraPartitioning", cond);

        scatra_manifold_dis->fill_complete();
      }
    }
  }
  else
  {
    if (fieldcoupling_ == Inpar::SSI::FieldCoupling::volume_match)
    {
      FOUR_C_THROW(
          "Reading a TRANSPORT discretization from the .dat file for the input parameter "
          "'FIELDCOUPLING volume_matching' in the SSI CONTROL section is not supported! As this "
          "coupling relies on matching node (and sometimes element) IDs, the ScaTra discretization "
          "is cloned from the structure discretization. Delete the ScaTra discretization from your "
          "input file.");
    }

    // copy conditions
    // this is actually only needed for copying TRANSPORT DIRICHLET/NEUMANN CONDITIONS
    // as standard DIRICHLET/NEUMANN CONDITIONS
    SSI::ScatraStructureCloneStrategy clonestrategy;
    const auto conditions_to_copy = clonestrategy.conditions_to_copy();
    Core::FE::DiscretizationCreatorBase creator;
    creator.CopyConditions(*scatradis, *scatradis, conditions_to_copy);

    // safety check, since it is not reasonable to have SOLIDSCATRA or SOLIDPOROP1 Elements with a
    // ScaTra::ImplType != 'impltype_undefined' if they are not cloned! Therefore loop over all
    // structure elements and check the impltype
    for (int i = 0; i < structdis->NumMyColElements(); ++i)
    {
      if (clonestrategy.GetImplType(structdis->lColElement(i)) != Inpar::ScaTra::impltype_undefined)
      {
        FOUR_C_THROW(
            "A TRANSPORT discretization is read from the .dat file, which is fine since the scatra "
            "discretization is not cloned from the structure discretization. But in the STRUCTURE "
            "ELEMENTS section of the .dat file an ImplType that is NOT 'Undefined' is prescribed "
            "which does not make sense if you don't want to clone the structure discretization. "
            "Change the ImplType to 'Undefined' or decide to clone the scatra discretization from "
            "the structure discretization.");
      }
    }
  }
  // read in the micro field, has to be done after cloning of the scatra discretization
  auto input_file_name = problem->OutputControlFile()->input_file_name();
  Core::IO::DatFileReader local_reader(input_file_name);
  Global::ReadMicroFields(*problem, local_reader);
}

/*----------------------------------------------------------------------*
 | Setup ssi coupling object                                rauch 08/16 |
 *----------------------------------------------------------------------*/
SSI::RedistributionType SSI::SSIBase::InitFieldCoupling(const std::string& struct_disname)
{
  // initialize return variable
  RedistributionType redistribution_required = SSI::RedistributionType::none;

  auto scatra_integrator = ScaTraBaseAlgorithm()->ScaTraField();

  // safety check
  {
    // check for ssi coupling condition
    std::vector<Core::Conditions::Condition*> ssicoupling;
    scatra_integrator->discretization()->GetCondition("SSICoupling", ssicoupling);
    const bool havessicoupling = (ssicoupling.size() > 0);

    if (havessicoupling and (fieldcoupling_ != Inpar::SSI::FieldCoupling::boundary_nonmatch and
                                fieldcoupling_ != Inpar::SSI::FieldCoupling::volumeboundary_match))
    {
      FOUR_C_THROW(
          "SSICoupling condition only valid in combination with FIELDCOUPLING set to "
          "'boundary_nonmatching' or 'volumeboundary_matching' in SSI DYNAMIC section. ");
    }

    if (fieldcoupling_ == Inpar::SSI::FieldCoupling::volume_nonmatch)
    {
      const Teuchos::ParameterList& volmortarparams =
          Global::Problem::Instance()->VolmortarParams();
      if (Core::UTILS::IntegralValue<Core::VolMortar::CouplingType>(
              volmortarparams, "COUPLINGTYPE") != Core::VolMortar::couplingtype_coninter)
      {
        FOUR_C_THROW(
            "Volmortar coupling only tested for consistent interpolation, "
            "i.e. 'COUPLINGTYPE consint' in VOLMORTAR COUPLING section. Try other couplings "
            "at "
            "own "
            "risk.");
      }
    }
    if (is_sca_tra_manifold() and fieldcoupling_ != Inpar::SSI::FieldCoupling::volumeboundary_match)
      FOUR_C_THROW("Solving manifolds only in combination with matching volumes and boundaries");
  }

  // build SSI coupling class
  switch (fieldcoupling_)
  {
    case Inpar::SSI::FieldCoupling::volume_match:
      ssicoupling_ = Teuchos::rcp(new SSICouplingMatchingVolume());
      break;
    case Inpar::SSI::FieldCoupling::volume_nonmatch:
      ssicoupling_ = Teuchos::rcp(new SSICouplingNonMatchingVolume());
      // redistribution is still performed inside
      redistribution_required = SSI::RedistributionType::binning;
      break;
    case Inpar::SSI::FieldCoupling::boundary_nonmatch:
      ssicoupling_ = Teuchos::rcp(new SSICouplingNonMatchingBoundary());
      break;
    case Inpar::SSI::FieldCoupling::volumeboundary_match:
      ssicoupling_ = Teuchos::rcp(new SSICouplingMatchingVolumeAndBoundary());
      redistribution_required = SSI::RedistributionType::match;
      break;
    default:
      FOUR_C_THROW("unknown type of field coupling for SSI!");
  }

  // initialize coupling objects including dof sets
  Global::Problem* problem = Global::Problem::Instance();
  ssicoupling_->init(problem->NDim(), problem->GetDis(struct_disname), Teuchos::rcp(this, false));

  return redistribution_required;
}

/*----------------------------------------------------------------------*
 | read restart information for given time step (public)   vuong 01/12  |
 *----------------------------------------------------------------------*/
void SSI::SSIBase::read_restart(int restart)
{
  if (restart)
  {
    structure_->read_restart(restart);

    const Teuchos::ParameterList& ssidyn = Global::Problem::Instance()->SSIControlParams();
    const bool restartfromstructure =
        Core::UTILS::IntegralValue<int>(ssidyn, "RESTART_FROM_STRUCTURE");

    if (not restartfromstructure)  // standard restart
    {
      ScaTraField()->read_restart(restart);
      if (is_sca_tra_manifold()) ScaTraManifold()->read_restart(restart);
    }
    else  // restart from structure simulation
    {
      // Since there is no restart output for the scatra field available, we only have to fix the
      // time and step counter
      ScaTraField()->SetTimeStep(structure_->TimeOld(), restart);
      if (is_sca_tra_manifold()) ScaTraManifold()->SetTimeStep(structure_->TimeOld(), restart);
    }

    SetTimeStep(structure_->TimeOld(), restart);
  }

  // Material pointers to other field were deleted during read_restart().
  // They need to be reset.
  ssicoupling_->assign_material_pointers(
      structure_->discretization(), ScaTraField()->discretization());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIBase::TestResults(const Epetra_Comm& comm) const
{
  Global::Problem* problem = Global::Problem::Instance();

  problem->AddFieldTest(structure_->CreateFieldTest());
  problem->AddFieldTest(ScaTraBaseAlgorithm()->create_sca_tra_field_test());
  if (is_sca_tra_manifold())
    problem->AddFieldTest(sca_tra_manifold_base_algorithm()->create_sca_tra_field_test());
  problem->AddFieldTest(Teuchos::rcp(new SSI::SSIResultTest(Teuchos::rcp(this, false))));
  problem->TestAll(comm);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIBase::set_struct_solution(Teuchos::RCP<const Epetra_Vector> disp,
    Teuchos::RCP<const Epetra_Vector> vel, const bool set_mechanical_stress)
{
  // safety checks
  check_is_init();
  check_is_setup();

  set_mesh_disp(disp);
  set_velocity_fields(vel);

  if (set_mechanical_stress)
    set_mechanical_stress_state(modelevaluator_ssi_base_->get_mechanical_stress_state());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIBase::SetScatraSolution(Teuchos::RCP<const Epetra_Vector> phi) const
{
  // safety checks
  check_is_init();
  check_is_setup();

  ssicoupling_->SetScalarField(*structure_field()->discretization(), phi, 1);

  // set state for contact evaluation
  if (contact_strategy_nitsche_ != Teuchos::null) set_ssi_contact_states(phi);
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void SSI::SSIBase::set_ssi_contact_states(Teuchos::RCP<const Epetra_Vector> phi) const
{
  contact_strategy_nitsche_->set_state(Mortar::state_scalar, *phi);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIBase::set_micro_scatra_solution(Teuchos::RCP<const Epetra_Vector> phi) const
{
  // safety checks
  check_is_init();
  check_is_setup();

  ssicoupling_->SetScalarFieldMicro(*structure_field()->discretization(), phi, 2);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIBase::evaluate_and_set_temperature_field()
{
  // temperature is non primary variable. Only set, if function for temperature is given
  if (temperature_funct_num_ != -1)
  {
    // evaluate temperature at current time and put to scalar
    const double temperature =
        Global::Problem::Instance()
            ->FunctionById<Core::UTILS::FunctionOfTime>(temperature_funct_num_ - 1)
            .evaluate(Time());
    temperature_vector_->PutScalar(temperature);

    // set temperature vector to structure discretization
    ssicoupling_->SetTemperatureField(*structure_->discretization(), temperature_vector_);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIBase::set_velocity_fields(Teuchos::RCP<const Epetra_Vector> vel)
{
  // safety checks
  check_is_init();
  check_is_setup();

  ssicoupling_->set_velocity_fields(ScaTraBaseAlgorithm(), zeros_structure_, vel);
  if (is_sca_tra_manifold())
    ssicoupling_->set_velocity_fields(sca_tra_manifold_base_algorithm(), zeros_structure_, vel);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIBase::set_mechanical_stress_state(
    Teuchos::RCP<const Epetra_Vector> mechanical_stress_state) const
{
  check_is_init();
  check_is_setup();

  ssicoupling_->set_mechanical_stress_state(*ScaTraField()->discretization(),
      mechanical_stress_state, ScaTraField()->nds_two_tensor_quantity());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIBase::set_mesh_disp(Teuchos::RCP<const Epetra_Vector> disp)
{
  // safety checks
  check_is_init();
  check_is_setup();

  ssicoupling_->set_mesh_disp(ScaTraBaseAlgorithm(), disp);
  if (is_sca_tra_manifold()) ssicoupling_->set_mesh_disp(sca_tra_manifold_base_algorithm(), disp);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIBase::check_ssi_flags() const
{
  if (ScaTraField()->S2IKinetics())
  {
    if (!(SSIInterfaceContact() or ssi_interface_meshtying()))
    {
      FOUR_C_THROW(
          "You defined an 'S2IKinetics' condition in the input-file. However, neither an "
          "'SSIInterfaceContact' condition nor an 'ssi_interface_meshtying' condition defined. "
          "This "
          "is not reasonable!");
    }
  }

  const bool is_nitsche_penalty_adaptive(Core::UTILS::IntegralValue<int>(
      Global::Problem::Instance()->contact_dynamic_params(), "NITSCHE_PENALTY_ADAPTIVE"));

  if (SSIInterfaceContact() and is_nitsche_penalty_adaptive)
    FOUR_C_THROW("Adaptive nitsche penalty parameter currently not supported!");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIBase::set_dt_from_sca_tra_to_structure()
{
  structure_field()->set_dt(ScaTraField()->Dt());
  structure_field()->SetTimen(ScaTraField()->Time());
  structure_field()->post_update();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIBase::set_dt_from_sca_tra_to_manifold()
{
  ScaTraManifold()->set_dt(ScaTraField()->Dt());
  ScaTraManifold()->SetTimeStep(ScaTraField()->Time(), ScaTraField()->Step());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIBase::set_dt_from_sca_tra_to_ssi()
{
  // set values for this SSI algorithm
  SetTimeStep(ScaTraField()->Time(), Step());
  set_dt(ScaTraField()->Dt());

  // set values for other fields
  set_dt_from_sca_tra_to_structure();
  if (is_sca_tra_manifold()) set_dt_from_sca_tra_to_manifold();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIBase::Redistribute(const RedistributionType redistribution_type)
{
  Global::Problem* problem = Global::Problem::Instance();

  auto structdis = problem->GetDis("structure");
  auto scatradis = problem->GetDis("scatra");
  if (redistribution_type == SSI::RedistributionType::match and !is_sca_tra_manifold())
  {
    // first we bin the scatra discretization
    std::vector<Teuchos::RCP<Core::FE::Discretization>> dis;
    dis.push_back(scatradis);
    Teuchos::ParameterList binning_params = Global::Problem::Instance()->binning_strategy_params();
    Core::UTILS::AddEnumClassToParameterList<Core::FE::ShapeFunctionType>(
        "spatial_approximation_type", Global::Problem::Instance()->spatial_approximation_type(),
        binning_params);
    auto element_filter = [](const Core::Elements::Element* element)
    { return BINSTRATEGY::UTILS::SpecialElement::none; };
    auto rigid_sphere_radius = [](const Core::Elements::Element* element) { return 0.0; };
    auto correct_beam_center_node = [](const Core::Nodes::Node* node) { return node; };
    Core::Rebalance::RebalanceDiscretizationsByBinning(binning_params,
        Global::Problem::Instance()->OutputControlFile(), dis, element_filter, rigid_sphere_radius,
        correct_beam_center_node, false);

    Core::Rebalance::MatchElementDistributionOfMatchingConditionedElements(
        *scatradis, *scatradis, "ScatraHeteroReactionMaster", "ScatraHeteroReactionSlave");

    // now we redistribute the structure dis to match the scatra dis
    Core::Rebalance::MatchElementDistributionOfMatchingDiscretizations(*scatradis, *structdis);
  }
  else if (redistribution_type == SSI::RedistributionType::binning)
  {
    // create vector of discr.
    std::vector<Teuchos::RCP<Core::FE::Discretization>> dis;
    dis.push_back(structdis);
    dis.push_back(scatradis);

    Teuchos::ParameterList binning_params = Global::Problem::Instance()->binning_strategy_params();
    Core::UTILS::AddEnumClassToParameterList<Core::FE::ShapeFunctionType>(
        "spatial_approximation_type", Global::Problem::Instance()->spatial_approximation_type(),
        binning_params);
    auto element_filter = [](const Core::Elements::Element* element)
    { return BINSTRATEGY::UTILS::SpecialElement::none; };
    auto rigid_sphere_radius = [](const Core::Elements::Element* element) { return 0.0; };
    auto correct_beam_center_node = [](const Core::Nodes::Node* node) { return node; };
    Core::Rebalance::RebalanceDiscretizationsByBinning(binning_params,
        Global::Problem::Instance()->OutputControlFile(), dis, element_filter, rigid_sphere_radius,
        correct_beam_center_node, false);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<ScaTra::ScaTraTimIntImpl> SSI::SSIBase::ScaTraField() const
{
  return scatra_base_algorithm_->ScaTraField();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<ScaTra::ScaTraTimIntImpl> SSI::SSIBase::ScaTraManifold() const
{
  return scatra_manifold_base_algorithm_->ScaTraField();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIBase::init_time_integrators(const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& scatraparams, const Teuchos::ParameterList& structparams,
    const std::string& struct_disname, const std::string& scatra_disname, const bool isAle)
{
  // get the global problem
  auto* problem = Global::Problem::Instance();

  // time parameter handling
  // In case of different time stepping, time params have to be read from single field sections.
  // In case of equal time step size for all fields the time params are controlled solely by the
  // problem section (e.g. ---SSI DYNAMIC or ---CELL DYNAMIC).
  const auto* structtimeparams = &globaltimeparams;
  const auto* scatratimeparams = &globaltimeparams;
  if (diff_time_step_size_)
  {
    structtimeparams = &structparams;
    scatratimeparams = &scatraparams;
  }

  // we do not construct a structure, in case it was built externally and handed into this object
  if (struct_adapterbase_ptr_ == Teuchos::null)
  {
    // access the structural discretization
    auto structdis = problem->GetDis(struct_disname);

    // build structure based on new structural time integration
    if (structparams.get<std::string>("INT_STRATEGY") == "Standard")
    {
      struct_adapterbase_ptr_ = Adapter::build_structure_algorithm(structparams);

      // initialize structure base algorithm
      struct_adapterbase_ptr_->init(
          *structtimeparams, const_cast<Teuchos::ParameterList&>(structparams), structdis);
    }
    // build structure based on old structural time integration
    else if (structparams.get<std::string>("INT_STRATEGY") == "Old")
    {
      auto structure = Teuchos::rcp(new Adapter::StructureBaseAlgorithm(
          *structtimeparams, const_cast<Teuchos::ParameterList&>(structparams), structdis));
      structure_ = Teuchos::rcp_dynamic_cast<Adapter::SSIStructureWrapper>(
          structure->structure_field(), true);
      if (structure_ == Teuchos::null)
        FOUR_C_THROW("cast from Adapter::Structure to Adapter::SSIStructureWrapper failed");
    }
    else
    {
      FOUR_C_THROW(
          "Unknown time integration requested!\n"
          "Set parameter INT_STRATEGY to Standard in ---STRUCTURAL DYNAMIC section!\n"
          "If you want to use yet unsupported elements or algorithms,\n"
          "set INT_STRATEGY to Old in ---STRUCUTRAL DYNAMIC section!");
    }
  }

  // create and initialize scatra base algorithm.
  // scatra time integrator constructed and initialized inside.
  // mesh is written inside. cloning must happen before!
  scatra_base_algorithm_ = Teuchos::rcp(new Adapter::ScaTraBaseAlgorithm(*scatratimeparams,
      SSI::UTILS::ModifyScaTraParams(scatraparams),
      problem->SolverParams(scatraparams.get<int>("LINEAR_SOLVER")), scatra_disname, isAle));

  ScaTraBaseAlgorithm()->init();

  // create and initialize scatra base algorithm for manifolds
  if (is_sca_tra_manifold())
  {
    scatra_manifold_base_algorithm_ = Teuchos::rcp(new Adapter::ScaTraBaseAlgorithm(
        *scatratimeparams,
        SSI::UTILS::CloneScaTraManifoldParams(scatraparams, globaltimeparams.sublist("MANIFOLD")),
        problem->SolverParams(globaltimeparams.sublist("MANIFOLD").get<int>("LINEAR_SOLVER")),
        "scatra_manifold", isAle));

    sca_tra_manifold_base_algorithm()->init();
  }

  // do checks if adaptive time stepping is activated
  if (Core::UTILS::IntegralValue<bool>(globaltimeparams, "ADAPTIVE_TIMESTEPPING"))
    check_adaptive_time_stepping(scatraparams, structparams);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool SSI::SSIBase::do_calculate_initial_potential_field() const
{
  const auto ssi_params = Global::Problem::Instance()->SSIControlParams();
  const bool init_pot_calc =
      Core::UTILS::IntegralValue<bool>(ssi_params.sublist("ELCH"), "INITPOTCALC");

  return init_pot_calc and is_elch_sca_tra_tim_int_type();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool SSI::SSIBase::is_elch_sca_tra_tim_int_type() const
{
  const auto ssi_params = Global::Problem::Instance()->SSIControlParams();
  const auto scatra_type =
      Teuchos::getIntegralValue<Inpar::SSI::ScaTraTimIntType>(ssi_params, "SCATRATIMINTTYPE");

  return scatra_type == Inpar::SSI::ScaTraTimIntType::elch;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool SSI::SSIBase::IsRestart() const
{
  // get the global problem
  const auto* problem = Global::Problem::Instance();

  const int restartstep = problem->restart();

  return (restartstep > 0);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIBase::check_adaptive_time_stepping(
    const Teuchos::ParameterList& scatraparams, const Teuchos::ParameterList& structparams)
{
  // safety check: adaptive time stepping in one of the sub problems
  if (!Core::UTILS::IntegralValue<bool>(scatraparams, "ADAPTIVE_TIMESTEPPING"))
  {
    FOUR_C_THROW(
        "Must provide adaptive time stepping algorithm in one of the sub problems. (Currently "
        "just ScaTra)");
  }
  if (Core::UTILS::IntegralValue<int>(structparams.sublist("TIMEADAPTIVITY"), "KIND") !=
      Inpar::STR::timada_kind_none)
    FOUR_C_THROW("Adaptive time stepping in SSI currently just from ScaTra");
  if (Core::UTILS::IntegralValue<int>(structparams, "DYNAMICTYP") == Inpar::STR::dyna_ab2)
    FOUR_C_THROW("Currently, only one step methods are allowed for adaptive time stepping");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool SSI::SSIBase::check_s2_i_kinetics_condition_for_pseudo_contact(
    const std::string& struct_disname) const
{
  bool is_s2i_kinetic_with_pseudo_contact = false;

  auto structdis = Global::Problem::Instance()->GetDis(struct_disname);
  // get all s2i kinetics conditions
  std::vector<Core::Conditions::Condition*> s2ikinetics_conditons(0, nullptr);
  structdis->GetCondition("S2IKinetics", s2ikinetics_conditons);
  // get all ssi contact conditions
  std::vector<Core::Conditions::Condition*> ssi_contact_conditions;
  structdis->GetCondition("SSIInterfaceContact", ssi_contact_conditions);
  for (auto* s2ikinetics_cond : s2ikinetics_conditons)
  {
    if ((s2ikinetics_cond->parameters().get<int>("interface side") == Inpar::S2I::side_slave) and
        (s2ikinetics_cond->parameters().get<int>("kinetic model") !=
            Inpar::S2I::kinetics_nointerfaceflux) and
        (s2ikinetics_cond->parameters().get<int>("is_pseudo_contact") == 1))
    {
      is_s2i_kinetic_with_pseudo_contact = true;
      const int s2i_kinetics_condition_id = s2ikinetics_cond->parameters().get<int>("ConditionID");

      for (auto* contact_condition : ssi_contact_conditions)
      {
        if (contact_condition->parameters().get<int>("ConditionID") == s2i_kinetics_condition_id)
        {
          FOUR_C_THROW(
              "Pseudo contact formulation of s2i kinetics conditions does not make sense in "
              "combination with resolved contact formulation. Set the respective is_pseudo_contact "
              "flag to '0'");
        }
      }
    }
  }

  const bool do_output_cauchy_stress =
      Core::UTILS::IntegralValue<Inpar::STR::StressType>(
          Global::Problem::Instance()->IOParams(), "STRUCT_STRESS") == Inpar::STR::stress_cauchy;

  if (is_s2i_kinetic_with_pseudo_contact and !do_output_cauchy_stress)
  {
    FOUR_C_THROW(
        "Consideration fo pseudo contact with 'S2IKinetics' condition only possible when Cauchy "
        "stress output is written.");
  }

  return is_s2i_kinetic_with_pseudo_contact;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIBase::check_ssi_interface_conditions(const std::string& struct_disname) const
{
  // access the structural discretization
  auto structdis = Global::Problem::Instance()->GetDis(struct_disname);

  if (ssi_interface_meshtying())
    ScaTra::ScaTraUtils::CheckConsistencyWithS2IKineticsCondition(
        "ssi_interface_meshtying", structdis);

  // check scatra-structure-interaction contact condition
  if (SSIInterfaceContact())
  {
    // get ssi condition to be tested
    std::vector<Core::Conditions::Condition*> ssiconditions;
    structdis->GetCondition("SSIInterfaceContact", ssiconditions);
    SSI::UTILS::CheckConsistencyOfSSIInterfaceContactCondition(ssiconditions, structdis);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIBase::SetupSystem()
{
  if (ssiinterfacemeshtying_)
    ssi_structure_mesh_tying()->check_slave_side_has_dirichlet_conditions(
        structure_field()->GetDBCMapExtractor()->CondMap());
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void SSI::SSIBase::setup_model_evaluator()
{
  // register the model evaluator if s2i condition with pseudo contact is available
  if (is_s2_i_kinetics_with_pseudo_contact())
  {
    modelevaluator_ssi_base_ = Teuchos::rcp(new STR::MODELEVALUATOR::BaseSSI());
    structure_base_algorithm()->register_model_evaluator(
        "Basic Coupling Model", modelevaluator_ssi_base_);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::SSIBase::setup_contact_strategy()
{
  // get the contact solution strategy
  auto contact_solution_type = Core::UTILS::IntegralValue<Inpar::CONTACT::SolvingStrategy>(
      Global::Problem::Instance()->contact_dynamic_params(), "STRATEGY");

  if (contact_solution_type == Inpar::CONTACT::solution_nitsche)
  {
    if (Core::UTILS::IntegralValue<Inpar::STR::IntegrationStrategy>(
            Global::Problem::Instance()->structural_dynamic_params(), "INT_STRATEGY") !=
        Inpar::STR::int_standard)
      FOUR_C_THROW("ssi contact only with new structural time integration");

    // get the contact model evaluator and store a pointer to the strategy
    auto& model_evaluator_contact = dynamic_cast<STR::MODELEVALUATOR::Contact&>(
        structure_field()->ModelEvaluator(Inpar::STR::model_contact));
    contact_strategy_nitsche_ = Teuchos::rcp_dynamic_cast<CONTACT::NitscheStrategySsi>(
        model_evaluator_contact.strategy_ptr(), true);
  }
  else
    FOUR_C_THROW("Only Nitsche contact implemented for SSI problems at the moment!");
}

FOUR_C_NAMESPACE_CLOSE
