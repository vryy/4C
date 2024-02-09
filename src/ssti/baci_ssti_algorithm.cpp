/*----------------------------------------------------------------------*/
/*! \file
 \brief base class for all scalar structure algorithms

 \level 2

 *------------------------------------------------------------------------------------------------*/

#include "baci_ssti_algorithm.hpp"

#include "baci_adapter_scatra_base_algorithm.hpp"
#include "baci_adapter_str_factory.hpp"
#include "baci_adapter_str_ssiwrapper.hpp"
#include "baci_adapter_str_structure_new.hpp"
#include "baci_global_data.hpp"
#include "baci_inpar_ssti.hpp"
#include "baci_lib_utils_createdis.hpp"
#include "baci_linalg_utils_sparse_algebra_create.hpp"
#include "baci_linalg_utils_sparse_algebra_manipulation.hpp"
#include "baci_scatra_timint_implicit.hpp"
#include "baci_scatra_timint_meshtying_strategy_s2i.hpp"
#include "baci_scatra_utils.hpp"
#include "baci_ssi_utils.hpp"
#include "baci_ssti_monolithic.hpp"
#include "baci_ssti_resulttest.hpp"
#include "baci_ssti_utils.hpp"

BACI_NAMESPACE_OPEN

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
      interfacemeshtying_(GLOBAL::Problem::Instance()
                              ->GetDis("structure")
                              ->GetCondition("SSTIInterfaceMeshtying") != nullptr),
      isinit_(false),
      issetup_(false)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::SSTIAlgorithm::Init(const Epetra_Comm& comm,
    const Teuchos::ParameterList& sstitimeparams, const Teuchos::ParameterList& scatraparams,
    const Teuchos::ParameterList& thermoparams, const Teuchos::ParameterList& structparams)
{
  // reset the setup flag
  issetup_ = false;

  // get the global problem
  GLOBAL::Problem* problem = GLOBAL::Problem::Instance();

  problem->GetDis("structure")->FillComplete(true, true, true);
  problem->GetDis("scatra")->FillComplete(true, true, true);
  problem->GetDis("thermo")->FillComplete(true, true, true);

  // clone scatra discretization from structure discretization first. Afterwards, clone thermo
  // discretization from scatra discretization
  CloneDiscretizations(comm);

  Teuchos::RCP<DRT::Discretization> structuredis = problem->GetDis("structure");
  Teuchos::RCP<DRT::Discretization> scatradis = problem->GetDis("scatra");
  Teuchos::RCP<DRT::Discretization> thermodis = problem->GetDis("thermo");

  // safety check
  if (structparams.get<std::string>("INT_STRATEGY") == "Old")
    dserror("Old structural time integration is not supported");

  struct_adapterbase_ptr_ = ADAPTER::BuildStructureAlgorithm(structparams);

  // initialize structure base algorithm
  struct_adapterbase_ptr_->Init(
      sstitimeparams, const_cast<Teuchos::ParameterList&>(structparams), structuredis);

  // create and initialize scatra problem and thermo problem
  scatra_ = Teuchos::rcp(
      new ADAPTER::ScaTraBaseAlgorithm(sstitimeparams, SSI::UTILS::ModifyScaTraParams(scatraparams),
          problem->SolverParams(scatraparams.get<int>("LINEAR_SOLVER")), "scatra", true));
  scatra_->Init();
  scatra_->ScaTraField()->SetNumberOfDofSetDisplacement(1);
  scatra_->ScaTraField()->SetNumberOfDofSetVelocity(1);
  scatra_->ScaTraField()->SetNumberOfDofSetThermo(2);
  thermo_ = Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm(sstitimeparams,
      CloneThermoParams(scatraparams, thermoparams),
      problem->SolverParams(thermoparams.get<int>("LINEAR_SOLVER")), "thermo", true));
  thermo_->Init();
  thermo_->ScaTraField()->SetNumberOfDofSetDisplacement(1);
  thermo_->ScaTraField()->SetNumberOfDofSetVelocity(1);
  thermo_->ScaTraField()->SetNumberOfDofSetScaTra(2);
  thermo_->ScaTraField()->SetNumberOfDofSetThermo(3);

  // distribute dofsets among subproblems
  Teuchos::RCP<DRT::DofSetInterface> scatradofset = scatradis->GetDofSetProxy();
  Teuchos::RCP<DRT::DofSetInterface> structdofset = structuredis->GetDofSetProxy();
  Teuchos::RCP<DRT::DofSetInterface> thermodofset = thermodis->GetDofSetProxy();
  if (scatradis->AddDofSet(structdofset) != 1) dserror("unexpected dof sets in scatra field");
  if (scatradis->AddDofSet(thermodofset) != 2) dserror("unexpected dof sets in scatra field");
  if (structuredis->AddDofSet(scatradofset) != 1) dserror("unexpected dof sets in structure field");
  if (structuredis->AddDofSet(thermodofset) != 2) dserror("unexpected dof sets in structure field");
  if (thermodis->AddDofSet(structdofset) != 1) dserror("unexpected dof sets in thermo field");
  if (thermodis->AddDofSet(scatradofset) != 2) dserror("unexpected dof sets in thermo field");
  if (thermodis->AddDofSet(thermodofset) != 3) dserror("unexpected dof sets in thermo field");

  // is adaptive time stepping activated?
  if (INPUT::IntegralValue<bool>(sstitimeparams, "ADAPTIVE_TIMESTEPPING"))
  {
    // safety check: adaptive time stepping in one of the subproblems?
    if (!INPUT::IntegralValue<bool>(scatraparams, "ADAPTIVE_TIMESTEPPING"))
      dserror(
          "Must provide adaptive time stepping in one of the subproblems. (Currently just ScaTra)");
    if (INPUT::IntegralValue<int>(structparams.sublist("TIMEADAPTIVITY"), "KIND") !=
        INPAR::STR::timada_kind_none)
      dserror("Adaptive time stepping in SSI currently just from ScaTra");
    if (INPUT::IntegralValue<int>(structparams, "DYNAMICTYP") == INPAR::STR::dyna_ab2)
      dserror("Currently, only one step methods are allowed for adaptive time stepping");
  }

  // now we can finally fill our discretizations
  // reinitialization of the structural elements is
  // vital for parallelization here!
  problem->GetDis("structure")->FillComplete(true, true, true);
  problem->GetDis("scatra")->FillComplete(true, false, true);
  problem->GetDis("thermo")->FillComplete(true, false, true);

  isinit_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::SSTIAlgorithm::Setup()
{
  // get the global problem
  GLOBAL::Problem* problem = GLOBAL::Problem::Instance();

  // check initialization
  CheckIsInit();

  // set up scatra and thermo problem
  scatra_->ScaTraField()->Setup();
  thermo_->ScaTraField()->Setup();

  // pass initial scalar field to structural discretization to correctly compute initial
  // accelerations
  problem->GetDis("structure")->SetState(1, "scalarfield", scatra_->ScaTraField()->Phinp());
  problem->GetDis("structure")->SetState(2, "tempfield", thermo_->ScaTraField()->Phinp());

  // set up structural base algorithm
  struct_adapterbase_ptr_->Setup();

  // get wrapper and cast it to specific type
  if (structure_ == Teuchos::null)
    structure_ = Teuchos::rcp_dynamic_cast<ADAPTER::SSIStructureWrapper>(
        struct_adapterbase_ptr_->StructureField());
  if (structure_ == Teuchos::null) dserror("No valid pointer to ADAPTER::SSIStructureWrapper !");

  // check maps from subproblems
  if (scatra_->ScaTraField()->DofRowMap()->NumGlobalElements() == 0)
    dserror("Scalar transport discretization does not have any degrees of freedom!");
  if (thermo_->ScaTraField()->DofRowMap()->NumGlobalElements() == 0)
    dserror("Scalar transport discretization does not have any degrees of freedom!");
  if (structure_->DofRowMap()->NumGlobalElements() == 0)
    dserror("Structure discretization does not have any degrees of freedom!");

  // set up materials
  AssignMaterialPointers();

  // set up scatra-scatra interface coupling
  if (InterfaceMeshtying())
  {
    // check for consistent parameterization of these conditions
    SCATRA::SCATRAUTILS::CheckConsistencyWithS2IKineticsCondition(
        "SSTIInterfaceMeshtying", StructureField()->Discretization());

    // extract meshtying strategy for scatra-scatra interface coupling on scatra discretization
    meshtying_strategy_scatra_ = Teuchos::rcp_dynamic_cast<const SCATRA::MeshtyingStrategyS2I>(
        scatra_->ScaTraField()->Strategy());

    // safety checks
    if (meshtying_strategy_scatra_ == Teuchos::null)
      dserror("Invalid scatra-scatra interface coupling strategy!");
    if (meshtying_strategy_scatra_->CouplingType() != INPAR::S2I::coupling_matching_nodes)
      dserror("SSTI only implemented for interface coupling with matching interface nodes!");

    // extract meshtying strategy for scatra-scatra interface coupling on thermo discretization
    meshtying_strategy_thermo_ = Teuchos::rcp_dynamic_cast<const SCATRA::MeshtyingStrategyS2I>(
        thermo_->ScaTraField()->Strategy());
    if (meshtying_strategy_thermo_ == Teuchos::null)
      dserror("Invalid scatra-scatra interface coupling strategy!");
    if (meshtying_strategy_thermo_->CouplingType() != INPAR::S2I::coupling_matching_nodes)
      dserror("SSTI only implemented for interface coupling with matching interface nodes!");

    // setup everything for SSTI structure meshtying
    ssti_structure_meshtying_ = Teuchos::rcp(new SSI::UTILS::SSIMeshTying(
        "SSTIInterfaceMeshtying", structure_->Discretization(), true, true));
  }

  issetup_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::SSTIAlgorithm::CloneDiscretizations(const Epetra_Comm& comm)
{
  // The structure discretization is received from the input.
  // Then, the scatra discretization is cloned.
  // Then, the thermo discretization is cloned.

  GLOBAL::Problem* problem = GLOBAL::Problem::Instance();

  Teuchos::RCP<DRT::Discretization> structdis = problem->GetDis("structure");
  Teuchos::RCP<DRT::Discretization> scatradis = problem->GetDis("scatra");
  Teuchos::RCP<DRT::Discretization> thermodis = problem->GetDis("thermo");

  if (scatradis->NumGlobalNodes() == 0)
  {
    DRT::UTILS::CloneDiscretization<SSTI::SSTIScatraStructureCloneStrategy>(structdis, scatradis);
    scatradis->FillComplete();
    DRT::UTILS::CloneDiscretization<SSTI::SSTIScatraThermoCloneStrategy>(scatradis, thermodis);
    thermodis->FillComplete();
  }
  else
    dserror("Only matching nodes in SSTI");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::SSTIAlgorithm::ReadRestart(int restart)
{
  StructureField()->ReadRestart(restart);
  ScaTraField()->ReadRestart(restart);
  ThermoField()->ReadRestart(restart);

  SetTimeStep(structure_->TimeOld(), restart);

  // Material pointers to other field were deleted during ReadRestart(). They need to be reset.
  AssignMaterialPointers();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSTI::SSTIAlgorithm::TestResults(const Epetra_Comm& comm) const
{
  GLOBAL::Problem* problem = GLOBAL::Problem::Instance();

  problem->AddFieldTest(structure_->CreateFieldTest());
  problem->AddFieldTest(scatra_->CreateScaTraFieldTest());
  problem->AddFieldTest(thermo_->CreateScaTraFieldTest());
  problem->AddFieldTest(Teuchos::rcp(new SSTI::SSTIResultTest(*this)));
  problem->TestAll(comm);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSTI::SSTIAlgorithm::DistributeStructureSolution()
{
  ScaTraField()->ApplyMeshMovement(structure_->Dispnp());
  ThermoField()->ApplyMeshMovement(structure_->Dispnp());

  // convective velocity is set to zero
  const auto convective_velocity = CORE::LINALG::CreateVector(*structure_->DofRowMap());

  ScaTraField()->SetVelocityField(
      convective_velocity, Teuchos::null, structure_->Velnp(), Teuchos::null);
  ThermoField()->SetVelocityField(
      convective_velocity, Teuchos::null, structure_->Velnp(), Teuchos::null);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSTI::SSTIAlgorithm::DistributeScatraSolution()
{
  StructureField()->Discretization()->SetState(1, "scalarfield", scatra_->ScaTraField()->Phinp());
  ThermoField()->Discretization()->SetState(2, "scatra", scatra_->ScaTraField()->Phinp());

  if (interfacemeshtying_)
  {
    // pass master-side scatra degrees of freedom to thermo discretization
    const Teuchos::RCP<Epetra_Vector> imasterphinp =
        CORE::LINALG::CreateVector(*ScaTraField()->Discretization()->DofRowMap(), true);
    meshtying_strategy_scatra_->InterfaceMaps()->InsertVector(
        meshtying_strategy_scatra_->CouplingAdapter()->MasterToSlave(
            meshtying_strategy_scatra_->InterfaceMaps()->ExtractVector(*ScaTraField()->Phinp(), 2)),
        1, imasterphinp);
    ThermoField()->Discretization()->SetState(2, "imasterscatra", imasterphinp);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSTI::SSTIAlgorithm::DistributeThermoSolution()
{
  StructureField()->Discretization()->SetState(2, "tempfield", thermo_->ScaTraField()->Phinp());

  ScaTraField()->Discretization()->SetState(2, "thermo", thermo_->ScaTraField()->Phinp());

  if (interfacemeshtying_)
  {
    // extract master side temperatures and copy to slave side dof map
    const Teuchos::RCP<Epetra_Vector> imastertempnp =
        CORE::LINALG::CreateVector(*ThermoField()->Discretization()->DofRowMap(), true);
    meshtying_strategy_thermo_->InterfaceMaps()->InsertVector(
        meshtying_strategy_thermo_->CouplingAdapter()->MasterToSlave(
            meshtying_strategy_thermo_->InterfaceMaps()->ExtractVector(*ThermoField()->Phinp(), 2)),
        1, imastertempnp);

    // extract slave side temperatures
    const Teuchos::RCP<Epetra_Vector> islavetempnp =
        CORE::LINALG::CreateVector(*ThermoField()->Discretization()->DofRowMap(), true);
    meshtying_strategy_thermo_->InterfaceMaps()->InsertVector(
        meshtying_strategy_thermo_->InterfaceMaps()->ExtractVector(*ThermoField()->Phinp(), 1), 1,
        islavetempnp);

    // set master side temperature to thermo discretization
    ThermoField()->Discretization()->SetState(3, "imastertemp", imastertempnp);

    // set master and slave side temperature to scatra discretization
    ScaTraField()->Discretization()->SetState(2, "islavetemp", islavetempnp);
    ScaTraField()->Discretization()->SetState(2, "imastertemp", imastertempnp);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSTI::SSTIAlgorithm::DistributeSolutionAllFields()
{
  DistributeScatraSolution();
  DistributeStructureSolution();
  DistributeThermoSolution();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSTI::SSTIAlgorithm::DistributeDtFromScaTra()
{
  // get adapted time, timestep, and step (incremented)
  const double newtime = ScaTraField()->Time();
  const double newtimestep = ScaTraField()->Dt();
  const int newstep = Step();

  // change time step size of thermo according to ScaTra
  ThermoField()->SetDt(newtimestep);
  // change current time and time step number of thermo field according to ScaTra
  // The incremental reduction is is needed, because time and step are incremented in
  // PrepareTimeStep() of thermo field
  ThermoField()->SetTimeStep(newtime - newtimestep, newstep - 1);

  // change current time and time step number of structure according to ScaTra
  StructureField()->SetDt(newtimestep);
  StructureField()->SetTimen(newtime);
  StructureField()->PostUpdate();

  // change current time and time step of this algorithm according to ScaTra
  SetTimeStep(newtime, newstep);
  SetDt(newtimestep);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSTI::SSTIAlgorithm::AssignMaterialPointers()
{
  // scatra - structure
  const int numscatraelements = ScaTraField()->Discretization()->NumMyColElements();
  for (int i = 0; i < numscatraelements; ++i)
  {
    DRT::Element* scatratele = ScaTraField()->Discretization()->lColElement(i);
    const int gid = scatratele->Id();

    DRT::Element* structele = StructureField()->Discretization()->gElement(gid);

    // for coupling we add the source material to the target element and vice versa
    scatratele->AddMaterial(structele->Material());
    structele->AddMaterial(scatratele->Material());
  }

  // thermo - scatra
  const int numthermoelements = ThermoField()->Discretization()->NumMyColElements();
  for (int i = 0; i < numthermoelements; ++i)
  {
    DRT::Element* thermotele = ThermoField()->Discretization()->lColElement(i);
    const int gid = thermotele->Id();

    DRT::Element* scatraele = ScaTraField()->Discretization()->gElement(gid);

    // for coupling we add the source material to the target element and vice versa
    thermotele->AddMaterial(scatraele->Material());
    scatraele->AddMaterial(thermotele->Material());
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<SCATRA::ScaTraTimIntImpl> SSTI::SSTIAlgorithm::ScaTraField() const
{
  return scatra_->ScaTraField();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<SCATRA::ScaTraTimIntImpl> SSTI::SSTIAlgorithm::ThermoField() const
{
  return thermo_->ScaTraField();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSTI::SSTIAlgorithm::CheckIsInit()
{
  if (not isinit_) dserror("Init(...) was not called.");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::ParameterList SSTI::SSTIAlgorithm::CloneThermoParams(
    const Teuchos::ParameterList& scatraparams, const Teuchos::ParameterList& thermoparams)
{
  auto thermoparams_copy = Teuchos::ParameterList(scatraparams);

  switch (Teuchos::getIntegralValue<INPAR::SCATRA::InitialField>(thermoparams, "INITIALFIELD"))
  {
    case INPAR::SCATRA::initfield_field_by_function:
    {
      thermoparams_copy.set<std::string>("INITIALFIELD", "field_by_function");
      break;
    }
    case INPAR::SCATRA::initfield_field_by_condition:
    {
      thermoparams_copy.set<std::string>("INITIALFIELD", "field_by_condition");
      break;
    }
    default:
      dserror("Initial field type for thermo not supported");
  }

  thermoparams_copy.set<int>("INITFUNCNO", thermoparams.get<int>("INITTHERMOFUNCT"));
  thermoparams_copy.sublist("S2I COUPLING").set<std::string>("SLAVEONLY", "No");

  if (INPUT::IntegralValue<INPAR::SCATRA::OutputScalarType>(scatraparams, "OUTPUTSCALARS") !=
      INPAR::SCATRA::outputscalars_none)
    thermoparams_copy.set<bool>("output_file_name_discretization", true);

  // adaptive time stepping only from scatra
  thermoparams_copy.set<std::string>("ADAPTIVE_TIMESTEPPING", "No");

  return thermoparams_copy;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<SSTI::SSTIAlgorithm> SSTI::BuildSSTI(INPAR::SSTI::SolutionScheme coupling,
    const Epetra_Comm& comm, const Teuchos::ParameterList& sstiparams)
{
  Teuchos::RCP<SSTI::SSTIAlgorithm> ssti = Teuchos::null;
  switch (coupling)
  {
    case INPAR::SSTI::SolutionScheme::monolithic:
    {
      ssti = Teuchos::rcp(new SSTI::SSTIMono(comm, sstiparams));
      break;
    }
    default:
      dserror("unknown coupling algorithm for SSTI!");
  }
  return ssti;
}

BACI_NAMESPACE_CLOSE
