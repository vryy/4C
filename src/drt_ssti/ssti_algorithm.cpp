/*----------------------------------------------------------------------*/
/*! \file
 \brief base class for all scalar structure algorithms

 \level 1


 *------------------------------------------------------------------------------------------------*/

#include "ssti_algorithm.H"

#include "ssti_monolithic.H"
#include "ssti_utils.H"

#include "../drt_adapter/ad_str_factory.H"
#include "../drt_adapter/ad_str_ssiwrapper.H"
#include "../drt_adapter/ad_str_structure_new.H"

#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"

#include "../drt_inpar/inpar_ssti.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils_createdis.H"

#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_scatra/scatra_timint_meshtying_strategy_s2i.H"

#include "../drt_ssi/ssi_clonestrategy.H"

#include "../drt_sti/sti_clonestrategy.H"

#include "../linalg/linalg_utils_sparse_algebra_create.H"
#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"

#include <Teuchos_TimeMonitor.hpp>
#include <Epetra_Time.h>

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
      icoup_structure_(Teuchos::null),
      meshtying_strategy_scatra_(Teuchos::null),
      meshtying_strategy_thermo_(Teuchos::null),
      interfacemeshtying_(
          DRT::Problem::Instance()->GetDis("structure")->GetCondition("SSIInterfaceMeshtying") !=
          nullptr),
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
  DRT::Problem* problem = DRT::Problem::Instance();

  problem->GetDis("structure")->FillComplete(true, true, true);
  problem->GetDis("scatra")->FillComplete(true, true, true);
  problem->GetDis("thermo")->FillComplete(true, true, true);

  // clone scatra discretization from structure discretization first. Afterwards, clone thermo
  // discretization from scatra discretization
  CloneDiscretizations(comm);

  Teuchos::RCP<DRT::Discretization> structuredis = problem->GetDis("structure");
  Teuchos::RCP<DRT::Discretization> scatradis = problem->GetDis("scatra");
  Teuchos::RCP<DRT::Discretization> thermodis = problem->GetDis("thermo");

  if (structparams.get<std::string>("INT_STRATEGY") == "Standard")
  {
    struct_adapterbase_ptr_ = ADAPTER::STR::BuildStructureAlgorithm(structparams);

    // initialize structure base algorithm
    struct_adapterbase_ptr_->Init(
        sstitimeparams, const_cast<Teuchos::ParameterList&>(structparams), structuredis);
  }
  else if (structparams.get<std::string>("INT_STRATEGY") == "Old")
    dserror("Old structural time integration is not supported");
  else
    dserror("Unknown time integration requested!");

  // create and initialize scatra problem and thermo problem
  scatra_ = Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm());
  scatra_->Init(sstitimeparams, scatraparams,
      problem->SolverParams(scatraparams.get<int>("LINEAR_SOLVER")), "scatra", true);
  thermo_ = Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm());
  thermo_->Init(sstitimeparams, CloneThermoParams(scatraparams, thermoparams),
      problem->SolverParams(thermoparams.get<int>("LINEAR_SOLVER")), "thermo", true);

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
  if (DRT::INPUT::IntegralValue<bool>(sstitimeparams, "ADAPTIVE_TIMESTEPPING"))
  {
    // safety check: adaptive time stepping in one of the subproblems?
    if (!DRT::INPUT::IntegralValue<bool>(scatraparams, "ADAPTIVE_TIMESTEPPING"))
      dserror(
          "Must provide adaptive time stepping in one of the subproblems. (Currently just ScaTra)");
    if (DRT::INPUT::IntegralValue<int>(structparams.sublist("TIMEADAPTIVITY"), "KIND") !=
        INPAR::STR::timada_kind_none)
      dserror("Adaptive time stepping in SSI currently just from ScaTra");
    if (DRT::INPUT::IntegralValue<int>(structparams, "DYNAMICTYP") == INPAR::STR::dyna_ab2)
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
  DRT::Problem* problem = DRT::Problem::Instance();

  // check initialization
  CheckIsInit();

  // set up scatra and thermo problem
  scatra_->ScaTraField()->Setup();
  thermo_->ScaTraField()->Setup();

  // only if adapter base has not already been set up outside
  if (!struct_adapterbase_ptr_->IsSetup())
  {
    // pass initial scalar field to structural discretization to correctly compute initial
    // accelerations
    problem->GetDis("structure")->SetState(1, "scalarfield", scatra_->ScaTraField()->Phinp());
    problem->GetDis("structure")->SetState(2, "tempfield", thermo_->ScaTraField()->Phinp());

    // set up structural base algorithm
    struct_adapterbase_ptr_->Setup();

    // get wrapper and cast it to specific type
    // do not do so, in case the wrapper has already been set from outside
    if (structure_ == Teuchos::null)
      structure_ = Teuchos::rcp_dynamic_cast<::ADAPTER::SSIStructureWrapper>(
          struct_adapterbase_ptr_->StructureField());

    if (structure_ == Teuchos::null) dserror("No valid pointer to ADAPTER::SSIStructureWrapper !");
  }

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
    CheckConsistencyInterfaceMeshtyingCondition();

    // set up scatra-scatra interface coupling adapter for structure field
    SetupCouplingAdapterStructure();

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
  }

  issetup_ = true;
}

/*----------------------------------------------------------------------------------*
 *----------------------------------------------------------------------------------*/
void SSTI::SSTIAlgorithm::SetupCouplingAdapterStructure()
{
  // initialize integer vectors for global IDs of master-side and slave-side interface nodes on
  // structure discretization
  std::vector<int> inodegidvec_master;
  std::vector<int> inodegidvec_slave;

  // extract scatra-scatra interface coupling conditions from structure discretization
  std::vector<DRT::Condition*> conditions(0, nullptr);
  structure_->Discretization()->GetCondition("S2ICoupling", conditions);

  // loop over all conditions
  for (const auto& condition : conditions)
  {
    // extract interface side associated with current condition
    const int side = condition->GetInt("interface side");

    // extract nodes associated with current condition
    const std::vector<int>* const inodegids = condition->Nodes();

    for (const int inodegid : *inodegids)
    {
      // insert global ID of current node into associated vector only if node is owned by current
      // processor need to make sure that node is stored on current processor, otherwise cannot
      // resolve "->Owner()"
      if (structure_->Discretization()->HaveGlobalNode(inodegid) and
          structure_->Discretization()->gNode(inodegid)->Owner() ==
              structure_->Discretization()->Comm().MyPID())
        side == INPAR::S2I::side_master ? inodegidvec_master.push_back(inodegid)
                                        : inodegidvec_slave.push_back(inodegid);
    }
  }

  // remove potential duplicates from vectors
  std::sort(inodegidvec_master.begin(), inodegidvec_master.end());
  inodegidvec_master.erase(
      unique(inodegidvec_master.begin(), inodegidvec_master.end()), inodegidvec_master.end());
  std::sort(inodegidvec_slave.begin(), inodegidvec_slave.end());
  inodegidvec_slave.erase(
      unique(inodegidvec_slave.begin(), inodegidvec_slave.end()), inodegidvec_slave.end());

  // setup scatra-scatra interface coupling adapter for structure field
  icoup_structure_ = Teuchos::rcp(new ADAPTER::Coupling());
  icoup_structure_->SetupCoupling(*structure_->Discretization(), *structure_->Discretization(),
      inodegidvec_master, inodegidvec_slave, DRT::Problem::Instance()->NDim(), true, 1.0e-8);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::SSTIAlgorithm::CloneDiscretizations(const Epetra_Comm& comm)
{
  // The structure discretization is received from the input.
  // Then, the scatra discretization is cloned.
  // Then, the thermo discretization is cloned.

  DRT::Problem* problem = DRT::Problem::Instance();

  Teuchos::RCP<DRT::Discretization> structdis = problem->GetDis("structure");
  Teuchos::RCP<DRT::Discretization> scatradis = problem->GetDis("scatra");
  Teuchos::RCP<DRT::Discretization> thermodis = problem->GetDis("thermo");

  if (scatradis->NumGlobalNodes() == 0)
  {
    DRT::UTILS::CloneDiscretization<SSTI::SSTIScatraStructureCloneStrategy>(structdis, scatradis);
    scatradis->FillComplete();
    DRT::UTILS::CloneDiscretization<STI::ScatraThermoCloneStrategy>(scatradis, thermodis);
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
  DRT::Problem* problem = DRT::Problem::Instance();

  problem->AddFieldTest(structure_->CreateFieldTest());
  problem->AddFieldTest(scatra_->CreateScaTraFieldTest());
  problem->AddFieldTest(thermo_->CreateScaTraFieldTest());
  problem->TestAll(comm);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSTI::SSTIAlgorithm::SetStructureSolution()
{
  ScaTraField()->ApplyMeshMovement(structure_->Dispnp(), 1);
  ThermoField()->ApplyMeshMovement(structure_->Dispnp(), 1);

  const auto structure_dofs_zero = LINALG::CreateVector(*structure_->DofRowMap());

  ScaTraField()->SetVelocityField(
      structure_dofs_zero, Teuchos::null, structure_->Velnp(), Teuchos::null, 1);
  ThermoField()->SetVelocityField(
      structure_dofs_zero, Teuchos::null, structure_->Velnp(), Teuchos::null, 1);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSTI::SSTIAlgorithm::SetScatraSolution()
{
  StructureField()->Discretization()->SetState(1, "scalarfield", scatra_->ScaTraField()->Phinp());
  ThermoField()->Discretization()->SetState(2, "scatra", scatra_->ScaTraField()->Phinp());

  if (interfacemeshtying_)
  {
    // pass master-side scatra degrees of freedom to thermo discretization
    const Teuchos::RCP<Epetra_Vector> imasterphinp =
        LINALG::CreateVector(*ScaTraField()->Discretization()->DofRowMap(), true);
    meshtying_strategy_scatra_->InterfaceMaps()->InsertVector(
        meshtying_strategy_scatra_->CouplingAdapter()->MasterToSlave(
            meshtying_strategy_scatra_->InterfaceMaps()->ExtractVector(*ScaTraField()->Phinp(), 2)),
        1, imasterphinp);
    ThermoField()->Discretization()->SetState(2, "imasterscatra", imasterphinp);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSTI::SSTIAlgorithm::SetThermoSolution()
{
  StructureField()->Discretization()->SetState(2, "tempfield", thermo_->ScaTraField()->Phinp());

  if (interfacemeshtying_)
  {
    // pass master-side thermo degrees of freedom to thermo discretization
    const Teuchos::RCP<Epetra_Vector> imastertempnp =
        LINALG::CreateVector(*ThermoField()->Discretization()->DofRowMap(), true);
    meshtying_strategy_thermo_->InterfaceMaps()->InsertVector(
        meshtying_strategy_thermo_->CouplingAdapter()->MasterToSlave(
            meshtying_strategy_thermo_->InterfaceMaps()->ExtractVector(*ThermoField()->Phinp(), 2)),
        1, imastertempnp);
    ThermoField()->Discretization()->SetState(3, "imastertemp", imastertempnp);

    // averaged temperature at interface (from slave and master side)
    auto itempnp = Teuchos::rcp(new Epetra_Vector(*imastertempnp));
    meshtying_strategy_thermo_->InterfaceMaps()->AddVector(
        meshtying_strategy_thermo_->InterfaceMaps()->ExtractVector(
            *thermo_->ScaTraField()->Phinp(), 1),
        1, itempnp, 1.0);
    itempnp->Scale(0.5);

    ScaTraField()->Discretization()->SetState(2, "itemp", itempnp);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSTI::SSTIAlgorithm::SetSolutionAllFields()
{
  SetScatraSolution();
  SetStructureSolution();
  SetThermoSolution();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSTI::SSTIAlgorithm::CheckConsistencyInterfaceMeshtyingCondition()
{
  // access the structural discretization
  Teuchos::RCP<DRT::Discretization> structdis = structure_->Discretization();

  // get ssi and s2i conditions to be tested
  std::vector<DRT::Condition*> ssiconditions;
  structdis->GetCondition("SSIInterfaceMeshtying", ssiconditions);
  std::vector<DRT::Condition*> s2iconditions;
  structdis->GetCondition("S2ICoupling", s2iconditions);

  // loop over all ssi conditions and check for a consistent initialization of the s2i conditions
  for (const auto* ssicondition : ssiconditions)
  {
    bool matchingconditions(false);
    bool isslave(true);
    const int s2icouplingid = ssicondition->GetInt("S2ICouplingID");
    const auto* side = ssicondition->Get<std::string>("Side");
    // check interface side
    if (*side == "Slave")
      isslave = true;
    else if (*side == "Master")
      isslave = false;
    else
    {
      dserror(
          "Interface side of SSIInterfaceMeshtying condition not recognized, has to be either "
          "'Slave' or 'Master'");
    }

    // loop over all s2i conditions to find the one that is matching the current ssi condition
    for (const auto* s2icondition : s2iconditions)
    {
      const int s2iconditionid = s2icondition->GetInt("ConditionID");
      // only do further checks if Ids match
      if (s2icouplingid != s2iconditionid) continue;

      // check the interface side
      switch (s2icondition->GetInt("interface side"))
      {
        case INPAR::S2I::side_slave:
        {
          if (isslave) matchingconditions = DRT::UTILS::HaveSameNodes(ssicondition, s2icondition);
          break;
        }
        case INPAR::S2I::side_master:
        {
          if (!isslave) matchingconditions = DRT::UTILS::HaveSameNodes(ssicondition, s2icondition);
          break;
        }
        default:
        {
          dserror("interface side of 'S2iCondition' has to be either 'Slave' or 'Master'");
          break;
        }
      }
    }

    if (not matchingconditions)
    {
      dserror(
          "Did not find 'S2ICoupling' condition with ID: %i and interface side: %s as defined in "
          "the 'SSIInterfaceMeshtying' condition",
          s2icouplingid, side->c_str());
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSTI::SSTIAlgorithm::DistributeDtFromScaTra()
{
  // get adapted time, timestep, and step (incremented)
  const double newtime = ScaTraField()->Time();
  const double newtimestep = ScaTraField()->Dt();
  const int newstep = Step();

  // change current time and time step of thermo according to ScaTra
  ThermoField()->SetDt(newtimestep);
  ThermoField()->SetTimeStep(newtime, newstep);

  // change current time and time step of structure according to ScaTra
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
const Teuchos::RCP<SCATRA::ScaTraTimIntImpl> SSTI::SSTIAlgorithm::ScaTraField() const
{
  return scatra_->ScaTraField();
};

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
const Teuchos::RCP<SCATRA::ScaTraTimIntImpl> SSTI::SSTIAlgorithm::ThermoField() const
{
  return thermo_->ScaTraField();
};

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSTI::SSTIAlgorithm::CheckIsInit()
{
  if (not isinit_) dserror("Init(...) was not called.");
};

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::ParameterList SSTI::SSTIAlgorithm::CloneThermoParams(
    const Teuchos::ParameterList& scatraparams, const Teuchos::ParameterList& thermoparams)
{
  auto* thermoparams_copy = new Teuchos::ParameterList(scatraparams);

  switch (Teuchos::getIntegralValue<INPAR::SCATRA::InitialField>(thermoparams, "INITIALFIELD"))
  {
    case INPAR::SCATRA::initfield_field_by_function:
    {
      thermoparams_copy->set<std::string>("INITIALFIELD", "field_by_function");
      break;
    }
    case INPAR::SCATRA::initfield_field_by_condition:
    {
      thermoparams_copy->set<std::string>("INITIALFIELD", "field_by_condition");
      break;
    }
    default:
      dserror("Initial field type for thermo not supported");
      break;
  }

  thermoparams_copy->set<int>("INITFUNCNO", thermoparams.get<int>("INITTHERMOFUNCT"));
  thermoparams_copy->sublist("S2I COUPLING").set<std::string>("SLAVEONLY", "No");
  thermoparams_copy->set<std::string>("OUTPUTSCALARS", "none");

  // adaptive time stepping only from scatra
  thermoparams_copy->set<std::string>("ADAPTIVE_TIMESTEPPING", "No");

  return *thermoparams_copy;
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
      break;
  }
  return ssti;
}
