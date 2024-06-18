/*----------------------------------------------------------------------*/
/*! \file
 \brief base class for all scalar structure algorithms

 \level 2

 *------------------------------------------------------------------------------------------------*/

#include "4C_ssti_algorithm.hpp"

#include "4C_adapter_scatra_base_algorithm.hpp"
#include "4C_adapter_str_factory.hpp"
#include "4C_adapter_str_ssiwrapper.hpp"
#include "4C_adapter_str_structure_new.hpp"
#include "4C_fem_general_utils_createdis.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_ssti.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_scatra_timint_implicit.hpp"
#include "4C_scatra_timint_meshtying_strategy_s2i.hpp"
#include "4C_scatra_utils.hpp"
#include "4C_ssi_utils.hpp"
#include "4C_ssti_monolithic.hpp"
#include "4C_ssti_resulttest.hpp"
#include "4C_ssti_utils.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
SSTI::SSTIAlgorithm::SSTIAlgorithm(
    const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams)
    : AlgorithmBase(comm, globaltimeparams),
      iter_(0),
      scatra_(Teuchos::null),
      structure_(Teuchos::null),
      struct_adapterbase_ptr_(Teuchos::null),
      thermo_(Teuchos::null),
      meshtying_strategy_scatra_(Teuchos::null),
      meshtying_strategy_thermo_(Teuchos::null),
      ssti_structure_meshtying_(Teuchos::null),
      interfacemeshtying_(Global::Problem::Instance()
                              ->GetDis("structure")
                              ->GetCondition("SSTIInterfaceMeshtying") != nullptr),
      isinit_(false),
      issetup_(false)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::SSTIAlgorithm::init(const Epetra_Comm& comm,
    const Teuchos::ParameterList& sstitimeparams, const Teuchos::ParameterList& scatraparams,
    const Teuchos::ParameterList& thermoparams, const Teuchos::ParameterList& structparams)
{
  // reset the setup flag
  issetup_ = false;

  // get the global problem
  Global::Problem* problem = Global::Problem::Instance();

  problem->GetDis("structure")->fill_complete(true, true, true);
  problem->GetDis("scatra")->fill_complete(true, true, true);
  problem->GetDis("thermo")->fill_complete(true, true, true);

  // clone scatra discretization from structure discretization first. Afterwards, clone thermo
  // discretization from scatra discretization
  clone_discretizations(comm);

  Teuchos::RCP<Core::FE::Discretization> structuredis = problem->GetDis("structure");
  Teuchos::RCP<Core::FE::Discretization> scatradis = problem->GetDis("scatra");
  Teuchos::RCP<Core::FE::Discretization> thermodis = problem->GetDis("thermo");

  // safety check
  if (structparams.get<std::string>("INT_STRATEGY") == "Old")
    FOUR_C_THROW("Old structural time integration is not supported");

  struct_adapterbase_ptr_ = Adapter::build_structure_algorithm(structparams);

  // initialize structure base algorithm
  struct_adapterbase_ptr_->init(
      sstitimeparams, const_cast<Teuchos::ParameterList&>(structparams), structuredis);

  // create and initialize scatra problem and thermo problem
  scatra_ = Teuchos::rcp(
      new Adapter::ScaTraBaseAlgorithm(sstitimeparams, SSI::UTILS::ModifyScaTraParams(scatraparams),
          problem->SolverParams(scatraparams.get<int>("LINEAR_SOLVER")), "scatra", true));
  scatra_->init();
  scatra_->ScaTraField()->set_number_of_dof_set_displacement(1);
  scatra_->ScaTraField()->set_number_of_dof_set_velocity(1);
  scatra_->ScaTraField()->set_number_of_dof_set_thermo(2);
  thermo_ = Teuchos::rcp(new Adapter::ScaTraBaseAlgorithm(sstitimeparams,
      clone_thermo_params(scatraparams, thermoparams),
      problem->SolverParams(thermoparams.get<int>("LINEAR_SOLVER")), "thermo", true));
  thermo_->init();
  thermo_->ScaTraField()->set_number_of_dof_set_displacement(1);
  thermo_->ScaTraField()->set_number_of_dof_set_velocity(1);
  thermo_->ScaTraField()->set_number_of_dof_set_sca_tra(2);
  thermo_->ScaTraField()->set_number_of_dof_set_thermo(3);

  // distribute dofsets among subproblems
  Teuchos::RCP<Core::DOFSets::DofSetInterface> scatradofset = scatradis->GetDofSetProxy();
  Teuchos::RCP<Core::DOFSets::DofSetInterface> structdofset = structuredis->GetDofSetProxy();
  Teuchos::RCP<Core::DOFSets::DofSetInterface> thermodofset = thermodis->GetDofSetProxy();
  if (scatradis->AddDofSet(structdofset) != 1) FOUR_C_THROW("unexpected dof sets in scatra field");
  if (scatradis->AddDofSet(thermodofset) != 2) FOUR_C_THROW("unexpected dof sets in scatra field");
  if (structuredis->AddDofSet(scatradofset) != 1)
    FOUR_C_THROW("unexpected dof sets in structure field");
  if (structuredis->AddDofSet(thermodofset) != 2)
    FOUR_C_THROW("unexpected dof sets in structure field");
  if (thermodis->AddDofSet(structdofset) != 1) FOUR_C_THROW("unexpected dof sets in thermo field");
  if (thermodis->AddDofSet(scatradofset) != 2) FOUR_C_THROW("unexpected dof sets in thermo field");
  if (thermodis->AddDofSet(thermodofset) != 3) FOUR_C_THROW("unexpected dof sets in thermo field");

  // is adaptive time stepping activated?
  if (Core::UTILS::IntegralValue<bool>(sstitimeparams, "ADAPTIVE_TIMESTEPPING"))
  {
    // safety check: adaptive time stepping in one of the subproblems?
    if (!Core::UTILS::IntegralValue<bool>(scatraparams, "ADAPTIVE_TIMESTEPPING"))
      FOUR_C_THROW(
          "Must provide adaptive time stepping in one of the subproblems. (Currently just ScaTra)");
    if (Core::UTILS::IntegralValue<int>(structparams.sublist("TIMEADAPTIVITY"), "KIND") !=
        Inpar::STR::timada_kind_none)
      FOUR_C_THROW("Adaptive time stepping in SSI currently just from ScaTra");
    if (Core::UTILS::IntegralValue<int>(structparams, "DYNAMICTYP") == Inpar::STR::dyna_ab2)
      FOUR_C_THROW("Currently, only one step methods are allowed for adaptive time stepping");
  }

  // now we can finally fill our discretizations
  // reinitialization of the structural elements is
  // vital for parallelization here!
  problem->GetDis("structure")->fill_complete(true, true, true);
  problem->GetDis("scatra")->fill_complete(true, false, true);
  problem->GetDis("thermo")->fill_complete(true, false, true);

  isinit_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::SSTIAlgorithm::setup()
{
  // get the global problem
  Global::Problem* problem = Global::Problem::Instance();

  // check initialization
  check_is_init();

  // set up scatra and thermo problem
  scatra_->ScaTraField()->setup();
  thermo_->ScaTraField()->setup();

  // pass initial scalar field to structural discretization to correctly compute initial
  // accelerations
  problem->GetDis("structure")->set_state(1, "scalarfield", scatra_->ScaTraField()->Phinp());
  problem->GetDis("structure")->set_state(2, "tempfield", thermo_->ScaTraField()->Phinp());

  // set up structural base algorithm
  struct_adapterbase_ptr_->setup();

  // get wrapper and cast it to specific type
  if (structure_ == Teuchos::null)
    structure_ = Teuchos::rcp_dynamic_cast<Adapter::SSIStructureWrapper>(
        struct_adapterbase_ptr_->structure_field());
  if (structure_ == Teuchos::null)
    FOUR_C_THROW("No valid pointer to Adapter::SSIStructureWrapper !");

  // check maps from subproblems
  if (scatra_->ScaTraField()->dof_row_map()->NumGlobalElements() == 0)
    FOUR_C_THROW("Scalar transport discretization does not have any degrees of freedom!");
  if (thermo_->ScaTraField()->dof_row_map()->NumGlobalElements() == 0)
    FOUR_C_THROW("Scalar transport discretization does not have any degrees of freedom!");
  if (structure_->dof_row_map()->NumGlobalElements() == 0)
    FOUR_C_THROW("Structure discretization does not have any degrees of freedom!");

  // set up materials
  assign_material_pointers();

  // set up scatra-scatra interface coupling
  if (interface_meshtying())
  {
    // check for consistent parameterization of these conditions
    ScaTra::ScaTraUtils::CheckConsistencyWithS2IKineticsCondition(
        "SSTIInterfaceMeshtying", structure_field()->discretization());

    // extract meshtying strategy for scatra-scatra interface coupling on scatra discretization
    meshtying_strategy_scatra_ = Teuchos::rcp_dynamic_cast<const ScaTra::MeshtyingStrategyS2I>(
        scatra_->ScaTraField()->Strategy());

    // safety checks
    if (meshtying_strategy_scatra_ == Teuchos::null)
      FOUR_C_THROW("Invalid scatra-scatra interface coupling strategy!");
    if (meshtying_strategy_scatra_->CouplingType() != Inpar::S2I::coupling_matching_nodes)
      FOUR_C_THROW("SSTI only implemented for interface coupling with matching interface nodes!");

    // extract meshtying strategy for scatra-scatra interface coupling on thermo discretization
    meshtying_strategy_thermo_ = Teuchos::rcp_dynamic_cast<const ScaTra::MeshtyingStrategyS2I>(
        thermo_->ScaTraField()->Strategy());
    if (meshtying_strategy_thermo_ == Teuchos::null)
      FOUR_C_THROW("Invalid scatra-scatra interface coupling strategy!");
    if (meshtying_strategy_thermo_->CouplingType() != Inpar::S2I::coupling_matching_nodes)
      FOUR_C_THROW("SSTI only implemented for interface coupling with matching interface nodes!");

    // setup everything for SSTI structure meshtying
    ssti_structure_meshtying_ = Teuchos::rcp(new SSI::UTILS::SSIMeshTying(
        "SSTIInterfaceMeshtying", structure_->discretization(), true, true));
  }

  issetup_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::SSTIAlgorithm::clone_discretizations(const Epetra_Comm& comm)
{
  // The structure discretization is received from the input.
  // Then, the scatra discretization is cloned.
  // Then, the thermo discretization is cloned.

  Global::Problem* problem = Global::Problem::Instance();

  Teuchos::RCP<Core::FE::Discretization> structdis = problem->GetDis("structure");
  Teuchos::RCP<Core::FE::Discretization> scatradis = problem->GetDis("scatra");
  Teuchos::RCP<Core::FE::Discretization> thermodis = problem->GetDis("thermo");

  if (scatradis->NumGlobalNodes() == 0)
  {
    Core::FE::CloneDiscretization<SSTI::SSTIScatraStructureCloneStrategy>(
        structdis, scatradis, Global::Problem::Instance()->CloningMaterialMap());
    scatradis->fill_complete();
    Core::FE::CloneDiscretization<SSTI::SSTIScatraThermoCloneStrategy>(
        scatradis, thermodis, Global::Problem::Instance()->CloningMaterialMap());
    thermodis->fill_complete();
  }
  else
    FOUR_C_THROW("Only matching nodes in SSTI");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::SSTIAlgorithm::read_restart(int restart)
{
  structure_field()->read_restart(restart);
  ScaTraField()->read_restart(restart);
  ThermoField()->read_restart(restart);

  SetTimeStep(structure_->TimeOld(), restart);

  // Material pointers to other field were deleted during read_restart(). They need to be reset.
  assign_material_pointers();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSTI::SSTIAlgorithm::TestResults(const Epetra_Comm& comm) const
{
  Global::Problem* problem = Global::Problem::Instance();

  problem->AddFieldTest(structure_->CreateFieldTest());
  problem->AddFieldTest(scatra_->create_sca_tra_field_test());
  problem->AddFieldTest(thermo_->create_sca_tra_field_test());
  problem->AddFieldTest(Teuchos::rcp(new SSTI::SSTIResultTest(*this)));
  problem->TestAll(comm);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSTI::SSTIAlgorithm::distribute_structure_solution()
{
  ScaTraField()->ApplyMeshMovement(structure_->Dispnp());
  ThermoField()->ApplyMeshMovement(structure_->Dispnp());

  // convective velocity is set to zero
  const auto convective_velocity = Core::LinAlg::CreateVector(*structure_->dof_row_map());

  ScaTraField()->set_velocity_field(
      convective_velocity, Teuchos::null, structure_->Velnp(), Teuchos::null);
  ThermoField()->set_velocity_field(
      convective_velocity, Teuchos::null, structure_->Velnp(), Teuchos::null);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSTI::SSTIAlgorithm::distribute_scatra_solution()
{
  structure_field()->discretization()->set_state(1, "scalarfield", scatra_->ScaTraField()->Phinp());
  ThermoField()->discretization()->set_state(2, "scatra", scatra_->ScaTraField()->Phinp());

  if (interfacemeshtying_)
  {
    // pass master-side scatra degrees of freedom to thermo discretization
    const Teuchos::RCP<Epetra_Vector> imasterphinp =
        Core::LinAlg::CreateVector(*ScaTraField()->discretization()->dof_row_map(), true);
    meshtying_strategy_scatra_->InterfaceMaps()->InsertVector(
        meshtying_strategy_scatra_->CouplingAdapter()->MasterToSlave(
            meshtying_strategy_scatra_->InterfaceMaps()->ExtractVector(*ScaTraField()->Phinp(), 2)),
        1, imasterphinp);
    ThermoField()->discretization()->set_state(2, "imasterscatra", imasterphinp);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSTI::SSTIAlgorithm::distribute_thermo_solution()
{
  structure_field()->discretization()->set_state(2, "tempfield", thermo_->ScaTraField()->Phinp());

  ScaTraField()->discretization()->set_state(2, "thermo", thermo_->ScaTraField()->Phinp());

  if (interfacemeshtying_)
  {
    // extract master side temperatures and copy to slave side dof map
    const Teuchos::RCP<Epetra_Vector> imastertempnp =
        Core::LinAlg::CreateVector(*ThermoField()->discretization()->dof_row_map(), true);
    meshtying_strategy_thermo_->InterfaceMaps()->InsertVector(
        meshtying_strategy_thermo_->CouplingAdapter()->MasterToSlave(
            meshtying_strategy_thermo_->InterfaceMaps()->ExtractVector(*ThermoField()->Phinp(), 2)),
        1, imastertempnp);

    // extract slave side temperatures
    const Teuchos::RCP<Epetra_Vector> islavetempnp =
        Core::LinAlg::CreateVector(*ThermoField()->discretization()->dof_row_map(), true);
    meshtying_strategy_thermo_->InterfaceMaps()->InsertVector(
        meshtying_strategy_thermo_->InterfaceMaps()->ExtractVector(*ThermoField()->Phinp(), 1), 1,
        islavetempnp);

    // set master side temperature to thermo discretization
    ThermoField()->discretization()->set_state(3, "imastertemp", imastertempnp);

    // set master and slave side temperature to scatra discretization
    ScaTraField()->discretization()->set_state(2, "islavetemp", islavetempnp);
    ScaTraField()->discretization()->set_state(2, "imastertemp", imastertempnp);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSTI::SSTIAlgorithm::distribute_solution_all_fields()
{
  distribute_scatra_solution();
  distribute_structure_solution();
  distribute_thermo_solution();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSTI::SSTIAlgorithm::distribute_dt_from_sca_tra()
{
  // get adapted time, timestep, and step (incremented)
  const double newtime = ScaTraField()->Time();
  const double newtimestep = ScaTraField()->Dt();
  const int newstep = Step();

  // change time step size of thermo according to ScaTra
  ThermoField()->set_dt(newtimestep);
  // change current time and time step number of thermo field according to ScaTra
  // The incremental reduction is is needed, because time and step are incremented in
  // prepare_time_step() of thermo field
  ThermoField()->SetTimeStep(newtime - newtimestep, newstep - 1);

  // change current time and time step number of structure according to ScaTra
  structure_field()->set_dt(newtimestep);
  structure_field()->SetTimen(newtime);
  structure_field()->post_update();

  // change current time and time step of this algorithm according to ScaTra
  SetTimeStep(newtime, newstep);
  set_dt(newtimestep);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSTI::SSTIAlgorithm::assign_material_pointers()
{
  // scatra - structure
  const int numscatraelements = ScaTraField()->discretization()->NumMyColElements();
  for (int i = 0; i < numscatraelements; ++i)
  {
    Core::Elements::Element* scatratele = ScaTraField()->discretization()->lColElement(i);
    const int gid = scatratele->Id();

    Core::Elements::Element* structele = structure_field()->discretization()->gElement(gid);

    // for coupling we add the source material to the target element and vice versa
    scatratele->AddMaterial(structele->Material());
    structele->AddMaterial(scatratele->Material());
  }

  // thermo - scatra
  const int numthermoelements = ThermoField()->discretization()->NumMyColElements();
  for (int i = 0; i < numthermoelements; ++i)
  {
    Core::Elements::Element* thermotele = ThermoField()->discretization()->lColElement(i);
    const int gid = thermotele->Id();

    Core::Elements::Element* scatraele = ScaTraField()->discretization()->gElement(gid);

    // for coupling we add the source material to the target element and vice versa
    thermotele->AddMaterial(scatraele->Material());
    scatraele->AddMaterial(thermotele->Material());
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<ScaTra::ScaTraTimIntImpl> SSTI::SSTIAlgorithm::ScaTraField() const
{
  return scatra_->ScaTraField();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<ScaTra::ScaTraTimIntImpl> SSTI::SSTIAlgorithm::ThermoField() const
{
  return thermo_->ScaTraField();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSTI::SSTIAlgorithm::check_is_init()
{
  if (not isinit_) FOUR_C_THROW("init(...) was not called.");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::ParameterList SSTI::SSTIAlgorithm::clone_thermo_params(
    const Teuchos::ParameterList& scatraparams, const Teuchos::ParameterList& thermoparams)
{
  auto thermoparams_copy = Teuchos::ParameterList(scatraparams);

  switch (Teuchos::getIntegralValue<Inpar::ScaTra::InitialField>(thermoparams, "INITIALFIELD"))
  {
    case Inpar::ScaTra::initfield_field_by_function:
    {
      thermoparams_copy.set<std::string>("INITIALFIELD", "field_by_function");
      break;
    }
    case Inpar::ScaTra::initfield_field_by_condition:
    {
      thermoparams_copy.set<std::string>("INITIALFIELD", "field_by_condition");
      break;
    }
    default:
      FOUR_C_THROW("Initial field type for thermo not supported");
  }

  thermoparams_copy.set<int>("INITFUNCNO", thermoparams.get<int>("INITTHERMOFUNCT"));
  thermoparams_copy.sublist("S2I COUPLING").set<std::string>("SLAVEONLY", "No");

  if (Core::UTILS::IntegralValue<Inpar::ScaTra::OutputScalarType>(scatraparams, "OUTPUTSCALARS") !=
      Inpar::ScaTra::outputscalars_none)
    thermoparams_copy.set<bool>("output_file_name_discretization", true);

  // adaptive time stepping only from scatra
  thermoparams_copy.set<std::string>("ADAPTIVE_TIMESTEPPING", "No");

  return thermoparams_copy;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<SSTI::SSTIAlgorithm> SSTI::BuildSSTI(Inpar::SSTI::SolutionScheme coupling,
    const Epetra_Comm& comm, const Teuchos::ParameterList& sstiparams)
{
  Teuchos::RCP<SSTI::SSTIAlgorithm> ssti = Teuchos::null;
  switch (coupling)
  {
    case Inpar::SSTI::SolutionScheme::monolithic:
    {
      ssti = Teuchos::rcp(new SSTI::SSTIMono(comm, sstiparams));
      break;
    }
    default:
      FOUR_C_THROW("unknown coupling algorithm for SSTI!");
  }
  return ssti;
}

FOUR_C_NAMESPACE_CLOSE
