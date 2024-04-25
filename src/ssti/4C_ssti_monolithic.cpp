/*--------------------------------------------------------------------------*/
/*! \file
\brief monolithic scalar-structure interaction

\level 2

*--------------------------------------------------------------------------*/

#include "4C_ssti_monolithic.hpp"

#include "4C_adapter_scatra_base_algorithm.hpp"
#include "4C_adapter_str_ssiwrapper.hpp"
#include "4C_adapter_str_structure_new.hpp"
#include "4C_global_data.hpp"
#include "4C_io_control.hpp"
#include "4C_lib_discret.hpp"
#include "4C_linalg_equilibrate.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_scatra_timint_implicit.hpp"
#include "4C_scatra_timint_meshtying_strategy_s2i.hpp"
#include "4C_ssi_monolithic_evaluate_OffDiag.hpp"
#include "4C_ssti_algorithm.hpp"
#include "4C_ssti_monolithic_assemble_strategy.hpp"
#include "4C_ssti_monolithic_evaluate_OffDiag.hpp"
#include "4C_ssti_utils.hpp"
#include "4C_sti_monolithic_evaluate_OffDiag.hpp"

FOUR_C_NAMESPACE_OPEN

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
SSTI::SSTIMono::SSTIMono(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams)
    : SSTIAlgorithm(comm, globaltimeparams),
      increment_(Teuchos::null),
      residual_(Teuchos::null),
      solver_(Teuchos::rcp(new CORE::LINALG::Solver(
          GLOBAL::Problem::Instance()->SolverParams(
              globaltimeparams.sublist("MONOLITHIC").get<int>("LINEAR_SOLVER")),
          comm))),
      scatrastructureoffdiagcoupling_(Teuchos::null),
      scatrathermooffdiagcoupling_(Teuchos::null),
      thermostructureoffdiagcoupling_(Teuchos::null),
      dtassemble_(0.0),
      dtevaluate_(0.0),
      dtnewton_(0.0),
      dtsolve_(0.0),
      timer_(Teuchos::rcp(new Teuchos::Time("SSTI_Monolithic", true))),
      equilibration_method_{Teuchos::getIntegralValue<CORE::LINALG::EquilibrationMethod>(
                                globaltimeparams.sublist("MONOLITHIC"), "EQUILIBRATION"),
          Teuchos::getIntegralValue<CORE::LINALG::EquilibrationMethod>(
              globaltimeparams.sublist("MONOLITHIC"), "EQUILIBRATION_SCATRA"),
          Teuchos::getIntegralValue<CORE::LINALG::EquilibrationMethod>(
              globaltimeparams.sublist("MONOLITHIC"), "EQUILIBRATION_STRUCTURE"),
          Teuchos::getIntegralValue<CORE::LINALG::EquilibrationMethod>(
              globaltimeparams.sublist("MONOLITHIC"), "EQUILIBRATION_THERMO")},
      matrixtype_(Teuchos::getIntegralValue<CORE::LINALG::MatrixType>(
          globaltimeparams.sublist("MONOLITHIC"), "MATRIXTYPE")),
      convcheck_(Teuchos::rcp(new SSTI::ConvCheckMono(globaltimeparams))),
      ssti_maps_mono_(Teuchos::null),
      ssti_matrices_(Teuchos::null),
      strategy_assemble_(Teuchos::null),
      strategy_equilibration_(Teuchos::null)
{
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSTI::SSTIMono::AssembleMatAndRHS()
{
  double starttime = timer_->wallTime();

  ssti_matrices_->SystemMatrix()->Zero();

  // assemble blocks of subproblems into system matrix
  strategy_assemble_->AssembleScatra(
      ssti_matrices_->SystemMatrix(), ScaTraField()->SystemMatrixOperator());
  strategy_assemble_->AssembleStructure(
      ssti_matrices_->SystemMatrix(), StructureField()->SystemMatrix());
  strategy_assemble_->AssembleThermo(
      ssti_matrices_->SystemMatrix(), ThermoField()->SystemMatrixOperator());

  // assemble domain contributions from coupling into system matrix
  strategy_assemble_->AssembleScatraStructure(ssti_matrices_->SystemMatrix(),
      ssti_matrices_->ScaTraStructureDomain(), ssti_matrices_->ScaTraStructureInterface());
  strategy_assemble_->AssembleStructureScatra(
      ssti_matrices_->SystemMatrix(), ssti_matrices_->StructureScaTraDomain());
  strategy_assemble_->AssembleThermoStructure(ssti_matrices_->SystemMatrix(),
      ssti_matrices_->ThermoStructureDomain(), ssti_matrices_->ThermoStructureInterface());
  strategy_assemble_->AssembleStructureThermo(
      ssti_matrices_->SystemMatrix(), ssti_matrices_->StructureThermoDomain());
  strategy_assemble_->AssembleThermoScatra(ssti_matrices_->SystemMatrix(),
      ssti_matrices_->ThermoScaTraDomain(), ssti_matrices_->ThermoScaTraInterface());
  strategy_assemble_->AssembleScatraThermoDomain(
      ssti_matrices_->SystemMatrix(), ssti_matrices_->ScaTraThermoDomain());

  // assemble interface contributions from coupling into system matrix
  if (InterfaceMeshtying())
  {
    strategy_assemble_->AssembleScatraThermoInterface(
        ssti_matrices_->SystemMatrix(), ssti_matrices_->ScaTraThermoInterface());
  }

  // apply meshtying on structural linearizations
  strategy_assemble_->ApplyMeshtyingSystemMatrix(ssti_matrices_->SystemMatrix());

  // finalize global system matrix
  ssti_matrices_->SystemMatrix()->Complete();

  // apply Dirichlet conditions
  ssti_matrices_->SystemMatrix()->ApplyDirichlet(*ScaTraField()->DirichMaps()->CondMap(), true);
  ssti_matrices_->SystemMatrix()->ApplyDirichlet(*ThermoField()->DirichMaps()->CondMap(), true);
  strategy_assemble_->ApplyStructuralDBCSystemMatrix(ssti_matrices_->SystemMatrix());

  // assemble RHS
  strategy_assemble_->AssembleRHS(
      residual_, ScaTraField()->Residual(), StructureField()->RHS(), ThermoField()->Residual());

  double mydt = timer_->wallTime() - starttime;
  Comm().MaxAll(&mydt, &dtassemble_, 1);
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void SSTI::SSTIMono::BuildNullSpaces()
{
  // build null spaces for scatra and thermo
  switch (ScaTraField()->MatrixType())
  {
    case CORE::LINALG::MatrixType::block_condition:
    case CORE::LINALG::MatrixType::block_condition_dof:
    {
      ScaTraField()->BuildBlockNullSpaces(
          solver_, GetBlockPositions(Subproblem::scalar_transport).at(0));
      ThermoField()->BuildBlockNullSpaces(solver_, GetBlockPositions(Subproblem::thermo).at(0));
      break;
    }
    case CORE::LINALG::MatrixType::sparse:
    {
      // equip smoother for scatra matrix block with empty parameter sub lists to trigger null space
      // computation
      std::ostringstream scatrablockstr;
      scatrablockstr << GetBlockPositions(Subproblem::scalar_transport).at(0) + 1;
      Teuchos::ParameterList& blocksmootherparamsscatra =
          solver_->Params().sublist("Inverse" + scatrablockstr.str());

      blocksmootherparamsscatra.sublist("Belos Parameters");
      blocksmootherparamsscatra.sublist("MueLu Parameters");

      // equip smoother for scatra matrix block with null space associated with all degrees of
      // freedom on scatra discretization
      ScaTraField()->Discretization()->ComputeNullSpaceIfNecessary(blocksmootherparamsscatra);

      std::ostringstream thermoblockstr;
      thermoblockstr << GetBlockPositions(Subproblem::thermo).at(0) + 1;
      Teuchos::ParameterList& blocksmootherparamsthermo =
          solver_->Params().sublist("Inverse" + thermoblockstr.str());
      blocksmootherparamsthermo.sublist("Belos Parameters");
      blocksmootherparamsthermo.sublist("MueLu Parameters");

      // equip smoother for scatra matrix block with null space associated with all degrees of
      // freedom on scatra discretization
      ThermoField()->Discretization()->ComputeNullSpaceIfNecessary(blocksmootherparamsthermo);
      break;
    }
    default:
    {
      FOUR_C_THROW("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }
  // build null spaces for structure
  {
    // store number of matrix block associated with structural field as string
    std::stringstream iblockstr;
    iblockstr << GetBlockPositions(Subproblem::structure).at(0) + 1;

    // equip smoother for structural matrix block with empty parameter sub lists to trigger null
    // space computation
    Teuchos::ParameterList& blocksmootherparams =
        solver_->Params().sublist("Inverse" + iblockstr.str());
    blocksmootherparams.sublist("Belos Parameters");
    blocksmootherparams.sublist("MueLu Parameters");

    // equip smoother for structural matrix block with null space associated with all degrees of
    // freedom on structural discretization
    StructureField()->Discretization()->ComputeNullSpaceIfNecessary(blocksmootherparams);
  }
}  // SSTI::SSTI_Mono::BuildNullSpaces

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSTI::SSTIMono::Init(const Epetra_Comm& comm, const Teuchos::ParameterList& sstitimeparams,
    const Teuchos::ParameterList& scatraparams, const Teuchos::ParameterList& thermoparams,
    const Teuchos::ParameterList& structparams)
{
  // check input parameters for scalar transport field
  if (CORE::UTILS::IntegralValue<INPAR::SCATRA::VelocityField>(scatraparams, "VELOCITYFIELD") !=
      INPAR::SCATRA::velocity_Navier_Stokes)
    FOUR_C_THROW("Invalid type of velocity field for scalar-structure interaction!");

  // call base class routine
  SSTIAlgorithm::Init(comm, sstitimeparams, scatraparams, thermoparams, structparams);
}


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSTI::SSTIMono::Output()
{
  // print finish line of convergence table to screen
  if (Comm().MyPID() == 0)
  {
    std::cout << "+------------+-------------------+--------------+--------------+--------------+--"
                 "------------+--------------+--------------+--------------+--------------+--------"
                 "------+"
              << std::endl;
    std::cout << "| Computation time for this timestep: " << std::setw(10) << TimeStatistics()[2]
              << "                                                                                 "
                 "                                       |"
              << std::endl;
    std::cout << "+--------------------------------------------------------------------------------"
                 "---------------------------------------------------------------------------------"
                 "------+"
              << std::endl;
  }

  ScaTraField()->CheckAndWriteOutputAndRestart();
  ThermoField()->CheckAndWriteOutputAndRestart();
  StructureField()->Output();
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSTI::SSTIMono::PrepareTimeStep()
{
  // update time and time step
  IncrementTimeAndStep();

  DistributeSolutionAllFields();

  // in first time step: solve to get initital derivatives
  ScaTraField()->PrepareTimeStep();

  // if adaptive time stepping and different time step size: calculate time step in scatra
  // (PrepareTimeStep() of Scatra) and pass to structure and thermo
  if (ScaTraField()->TimeStepAdapted()) DistributeDtFromScaTra();

  // in first time step: solve to get initital derivatives
  ThermoField()->PrepareTimeStep();

  // pass scalar transport degrees of freedom to structural discretization
  // has to be called AFTER ScaTraField()->PrepareTimeStep() to ensure
  // consistent scalar transport state vector with valid Dirichlet conditions
  StructureField()->PrepareTimeStep();

  ScaTraField()->PrintTimeStepInfo();
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSTI::SSTIMono::Setup()
{
  // call base class routine
  SSTIAlgorithm::Setup();

  // safety checks
  if (ScaTraField()->NumScal() != 1)
  {
    FOUR_C_THROW(
        "Since the ssti_monolithic framework is only implemented for usage in combination with "
        "volume change laws 'MAT_InelasticDefgradLinScalarIso' or  "
        "'MAT_InelasticDefgradLinScalarAniso' so far and these laws are implemented for only one "
        "transported scalar at the moment it is not reasonable to use them with more than one "
        "transported scalar. So you need to cope with it or change implementation! ;-)");
  }
  if (equilibration_method_.global != CORE::LINALG::EquilibrationMethod::local and
      (equilibration_method_.structure != CORE::LINALG::EquilibrationMethod::none or
          equilibration_method_.scatra != CORE::LINALG::EquilibrationMethod::none or
          equilibration_method_.thermo != CORE::LINALG::EquilibrationMethod::none))
    FOUR_C_THROW("Either global equilibration or local equilibration");

  if (matrixtype_ == CORE::LINALG::MatrixType::sparse and
      (equilibration_method_.structure != CORE::LINALG::EquilibrationMethod::none or
          equilibration_method_.scatra != CORE::LINALG::EquilibrationMethod::none or
          equilibration_method_.thermo != CORE::LINALG::EquilibrationMethod::none))
    FOUR_C_THROW("Block based equilibration only for block matrices");

  const bool equilibration_scatra_initial = CORE::UTILS::IntegralValue<bool>(
      GLOBAL::Problem::Instance()->SSTIControlParams().sublist("MONOLITHIC"),
      "EQUILIBRATION_INIT_SCATRA");
  const bool calc_initial_pot = CORE::UTILS::IntegralValue<bool>(
      GLOBAL::Problem::Instance()->ELCHControlParams(), "INITPOTCALC");

  if (!equilibration_scatra_initial and
      ScaTraField()->EquilibrationMethod() != CORE::LINALG::EquilibrationMethod::none)
  {
    FOUR_C_THROW(
        "You are within the monolithic SSTI framework but activated a pure scatra equilibration "
        "method. Delete this from 'SCALAR TRANSPORT DYNAMIC' section and set it in 'SSTI "
        "CONTROL/MONOLITHIC' instead.");
  }
  if (equilibration_scatra_initial and
      ScaTraField()->EquilibrationMethod() == CORE::LINALG::EquilibrationMethod::none)
  {
    FOUR_C_THROW(
        "You selected to equilibrate equations of initial potential but did not specify any "
        "equilibration method in ScaTra.");
  }
  if (equilibration_scatra_initial and !calc_initial_pot)
  {
    FOUR_C_THROW(
        "You selected to equilibrate equations of initial potential but did not activate "
        "INITPOTCALC in ELCH CONTROL");
  }

  if (!ScaTraField()->IsIncremental())
    FOUR_C_THROW("Must have incremental solution approach for monolithic SSTI!");
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSTI::SSTIMono::SetupSystem()
{
  if (InterfaceMeshtying())
    SSTIStructureMeshTying()->CheckSlaveSideHasDirichletConditions(
        StructureField()->GetDBCMapExtractor()->CondMap());

  // Setup all kind of maps
  ssti_maps_mono_ = Teuchos::rcp(new SSTI::SSTIMapsMono(*this));

  // initialize global increment vector for Newton-Raphson iteration
  increment_ = CORE::LINALG::CreateVector(*ssti_maps_mono_->MapsSubProblems()->FullMap(), true);

  // initialize global residual vector
  residual_ = CORE::LINALG::CreateVector(*ssti_maps_mono_->MapsSubProblems()->FullMap(), true);

  if (matrixtype_ == CORE::LINALG::MatrixType::block_field)
  {
    if (!solver_->Params().isSublist("AMGnxn Parameters"))
      FOUR_C_THROW(
          "Global system matrix with block structure requires AMGnxn block preconditioner!");

    // feed AMGnxn block preconditioner with null space information for each block of global
    // block system matrix
    BuildNullSpaces();
  }

  // initialize submatrices and system matrix
  ssti_matrices_ = Teuchos::rcp(new SSTI::SSTIMatrices(
      ssti_maps_mono_, matrixtype_, ScaTraField()->MatrixType(), InterfaceMeshtying()));

  // initialize strategy for assembly
  strategy_assemble_ = SSTI::BuildAssembleStrategy(
      Teuchos::rcp(this, false), matrixtype_, ScaTraField()->MatrixType());

  // initialize evaluation objects for coupling between subproblems
  scatrastructureoffdiagcoupling_ =
      Teuchos::rcp(new SSI::ScatraStructureOffDiagCouplingSSTI(ssti_maps_mono_->BlockMapStructure(),
          ssti_maps_mono_->MapsSubProblems()->Map(GetProblemPosition(Subproblem::scalar_transport)),
          ssti_maps_mono_->MapsSubProblems()->Map(GetProblemPosition(Subproblem::structure)),
          SSTIStructureMeshTying(), MeshtyingScatra(), ScaTraField(), StructureField()));

  thermostructureoffdiagcoupling_ = Teuchos::rcp(new SSTI::ThermoStructureOffDiagCoupling(
      ssti_maps_mono_->BlockMapStructure(), ssti_maps_mono_->BlockMapThermo(),
      ssti_maps_mono_->MapsSubProblems()->Map(GetProblemPosition(Subproblem::structure)),
      ssti_maps_mono_->MapsSubProblems()->Map(GetProblemPosition(Subproblem::thermo)),
      SSTIStructureMeshTying(), MeshtyingThermo(), StructureField(), ThermoFieldBase()));

  // Note: STI evaluation of off diagonal coupling is designed to use interface maps for the
  // interface coupling matrices. In SSTI we always use the full maps and thus hand in the same map
  // multiple times for both domain and interface contributions.
  scatrathermooffdiagcoupling_ = Teuchos::rcp(
      new STI::ScatraThermoOffDiagCouplingMatchingNodes(ssti_maps_mono_->BlockMapThermo(),
          ssti_maps_mono_->BlockMapThermo(), ssti_maps_mono_->BlockMapThermo(),
          ssti_maps_mono_->MapsSubProblems()->Map(GetProblemPosition(Subproblem::scalar_transport)),
          ssti_maps_mono_->MapsSubProblems()->Map(GetProblemPosition(Subproblem::thermo)),
          ssti_maps_mono_->MapsSubProblems()->Map(GetProblemPosition(Subproblem::scalar_transport)),
          ssti_maps_mono_->MapsSubProblems()->Map(GetProblemPosition(Subproblem::thermo)), true,
          MeshtyingScatra(), MeshtyingThermo(), ScaTraFieldBase(), ThermoFieldBase()));

  // initialize equilibration class
  strategy_equilibration_ = CORE::LINALG::BuildEquilibration(
      matrixtype_, GetBlockEquilibration(), AllMaps()->MapsSubProblems()->FullMap());
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSTI::SSTIMono::NewtonLoop()
{
  double starttime = timer_->wallTime();

  // initialize counter for Newton-Raphson iteration
  ResetIter();

  // start Newton-Raphson iteration
  while (true)
  {
    PrepareNewtonStep();

    ssti_matrices_->UnCompleteCouplingMatrices();

    EvaluateSubproblems();

    ssti_matrices_->CompleteCouplingMatrices();

    AssembleMatAndRHS();

    if (convcheck_->Converged(*this)) break;

    LinearSolve();

    UpdateIterStates();
  }

  double mydt = timer_->wallTime() - starttime;
  Comm().MaxAll(&mydt, &dtnewton_, 1);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSTI::SSTIMono::Timeloop()
{
  // output initial scalar transport solution to screen and files
  if (Step() == 0)
  {
    DistributeSolutionAllFields();

    ScaTraField()->PrepareTimeLoop();
    ThermoField()->PrepareTimeLoop();
  }
  // time loop
  while (NotFinished() and ScaTraField()->NotFinished())
  {
    PrepareTimeStep();

    NewtonLoop();

    constexpr bool force_prepare = false;
    StructureField()->PrepareOutput(force_prepare);

    Update();

    Output();
  }
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void SSTI::SSTIMono::Update()
{
  ScaTraField()->Update();
  ThermoField()->Update();
  StructureField()->Update();
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> SSTI::SSTIMono::ExtractSubIncrement(Subproblem sub)
{
  Teuchos::RCP<Epetra_Vector> subincrement(Teuchos::null);
  switch (sub)
  {
    case Subproblem::structure:
    {
      // First, extract increment from domain and master side
      subincrement = ssti_maps_mono_->MapsSubProblems()->ExtractVector(
          increment_, GetProblemPosition(Subproblem::structure));

      // Second, copy master side displacements and increments to slave side for meshtying
      if (InterfaceMeshtying())
      {
        for (const auto& meshtying : SSTIStructureMeshTying()->MeshTyingHandlers())
        {
          auto coupling_adapter = meshtying->SlaveMasterCoupling();
          auto coupling_map_extractor = meshtying->SlaveMasterExtractor();

          // displacements
          coupling_map_extractor->InsertVector(
              coupling_adapter->MasterToSlave(
                  coupling_map_extractor->ExtractVector(StructureField()->Dispnp(), 2)),
              1, StructureField()->WriteAccessDispnp());
          StructureField()->SetState(StructureField()->WriteAccessDispnp());
          // increments
          coupling_map_extractor->InsertVector(
              coupling_adapter->MasterToSlave(
                  coupling_map_extractor->ExtractVector(subincrement, 2)),
              1, subincrement);
        }
      }
      break;
    }
    case Subproblem::scalar_transport:
    {
      subincrement = ssti_maps_mono_->MapsSubProblems()->ExtractVector(
          increment_, GetProblemPosition(Subproblem::scalar_transport));
      break;
    }
    case Subproblem::thermo:
    {
      subincrement = ssti_maps_mono_->MapsSubProblems()->ExtractVector(
          increment_, GetProblemPosition(Subproblem::thermo));
      break;
    }
    default:
    {
      FOUR_C_THROW("Unknown type of subproblem in SSTI");
    }
  }
  return subincrement;
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void SSTI::SSTIMono::EvaluateSubproblems()
{
  double starttime = timer_->wallTime();

  // clear all matrices from previous Newton iteration
  ssti_matrices_->ClearMatrices();

  // needed to communicate to NOX state
  StructureField()->SetState(StructureField()->WriteAccessDispnp());

  // distribute solution from all fields to each other
  DistributeSolutionAllFields();

  // evaluate all subproblems
  StructureField()->Evaluate();
  ScaTraField()->PrepareLinearSolve();
  ThermoField()->PrepareLinearSolve();

  // evaluate domain contributions from coupling
  scatrastructureoffdiagcoupling_->EvaluateOffDiagBlockScatraStructureDomain(
      ssti_matrices_->ScaTraStructureDomain());
  scatrastructureoffdiagcoupling_->EvaluateOffDiagBlockStructureScatraDomain(
      ssti_matrices_->StructureScaTraDomain());
  thermostructureoffdiagcoupling_->EvaluateOffDiagBlockThermoStructureDomain(
      ssti_matrices_->ThermoStructureDomain());
  thermostructureoffdiagcoupling_->EvaluateOffDiagBlockStructureThermoDomain(
      ssti_matrices_->StructureThermoDomain());
  scatrathermooffdiagcoupling_->EvaluateOffDiagBlockThermoScatraDomain(
      ssti_matrices_->ThermoScaTraDomain());
  scatrathermooffdiagcoupling_->EvaluateOffDiagBlockScatraThermoDomain(
      ssti_matrices_->ScaTraThermoDomain());

  // evaluate interface contributions from coupling
  if (InterfaceMeshtying())
  {
    scatrastructureoffdiagcoupling_->EvaluateOffDiagBlockScatraStructureInterface(
        ssti_matrices_->ScaTraStructureInterface());
    thermostructureoffdiagcoupling_->EvaluateOffDiagBlockThermoStructureInterface(
        ssti_matrices_->ThermoStructureInterface());
    scatrathermooffdiagcoupling_->EvaluateOffDiagBlockThermoScatraInterface(
        ssti_matrices_->ThermoScaTraInterface());
    scatrathermooffdiagcoupling_->EvaluateOffDiagBlockScatraThermoInterface(
        ssti_matrices_->ScaTraThermoInterface());
  }

  double mydt = timer_->wallTime() - starttime;
  Comm().MaxAll(&mydt, &dtevaluate_, 1);
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void SSTI::SSTIMono::LinearSolve()
{
  double starttime = timer_->wallTime();

  increment_->PutScalar(0.0);

  if (!ssti_matrices_->SystemMatrix()->Filled())
    FOUR_C_THROW("Complete() has not been called on global system matrix yet!");

  strategy_equilibration_->EquilibrateSystem(
      ssti_matrices_->SystemMatrix(), residual_, AllMaps()->BlockMapSystemMatrix());

  CORE::LINALG::SolverParams solver_params;
  solver_params.refactor = true;
  solver_params.reset = Iter() == 1;

  solver_->Solve(
      ssti_matrices_->SystemMatrix()->EpetraOperator(), increment_, residual_, solver_params);

  strategy_equilibration_->UnequilibrateIncrement(increment_);

  double mydt = timer_->wallTime() - starttime;
  Comm().MaxAll(&mydt, &dtsolve_, 1);
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void SSTI::SSTIMono::UpdateIterStates()
{
  ScaTraField()->UpdateIter(ExtractSubIncrement(Subproblem::scalar_transport));
  ScaTraField()->ComputeIntermediateValues();

  ThermoField()->UpdateIter(ExtractSubIncrement(Subproblem::thermo));
  ThermoField()->ComputeIntermediateValues();

  StructureField()->UpdateStateIncrementally(ExtractSubIncrement(Subproblem::structure));
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void SSTI::SSTIMono::PrepareNewtonStep()
{
  // update iteration counter
  IncrementIter();

  // reset timer
  timer_->reset();

  ssti_matrices_->SystemMatrix()->Zero();
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
std::vector<int> SSTI::SSTIMono::GetBlockPositions(Subproblem subproblem) const
{
  if (matrixtype_ == CORE::LINALG::MatrixType::sparse)
    FOUR_C_THROW("Sparse matrices have just one block");

  auto block_position = std::vector<int>(0);

  switch (subproblem)
  {
    case Subproblem::structure:
    {
      if (ScaTraField()->MatrixType() == CORE::LINALG::MatrixType::sparse)
        block_position.emplace_back(1);
      else
        block_position.emplace_back(ScaTraField()->BlockMaps()->NumMaps());
      break;
    }
    case Subproblem::scalar_transport:
    {
      if (ScaTraField()->MatrixType() == CORE::LINALG::MatrixType::sparse)
        block_position.emplace_back(0);
      else

      {
        for (int i = 0; i < static_cast<int>(ScaTraField()->BlockMaps()->NumMaps()); ++i)
          block_position.emplace_back(i);
      }
      break;
    }
    case Subproblem::thermo:
    {
      if (ThermoField()->MatrixType() == CORE::LINALG::MatrixType::sparse)
        block_position.emplace_back(2);
      else
      {
        for (int i = 0; i < static_cast<int>(ThermoField()->BlockMaps()->NumMaps()); ++i)
          block_position.emplace_back(ScaTraField()->BlockMaps()->NumMaps() + 1 + i);
      }
      break;
    }
    default:
    {
      FOUR_C_THROW("Unknown type of subproblem");
      break;
    }
  }

  return block_position;
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
int SSTI::SSTIMono::GetProblemPosition(Subproblem subproblem) const
{
  int position = -1;

  switch (subproblem)
  {
    case Subproblem::structure:
    {
      position = 1;
      break;
    }
    case Subproblem::scalar_transport:
    {
      position = 0;
      break;
    }
    case Subproblem::thermo:
    {
      position = 2;
      break;
    }
    default:
    {
      FOUR_C_THROW("Unknown type of subproblem");
      break;
    }
  }

  return position;
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
std::vector<CORE::LINALG::EquilibrationMethod> SSTI::SSTIMono::GetBlockEquilibration()
{
  std::vector<CORE::LINALG::EquilibrationMethod> equilibration_method_vector;
  switch (matrixtype_)
  {
    case CORE::LINALG::MatrixType::sparse:
    {
      equilibration_method_vector =
          std::vector<CORE::LINALG::EquilibrationMethod>(1, equilibration_method_.global);
      break;
    }
    case CORE::LINALG::MatrixType::block_field:
    {
      if (equilibration_method_.global != CORE::LINALG::EquilibrationMethod::local)
      {
        equilibration_method_vector =
            std::vector<CORE::LINALG::EquilibrationMethod>(1, equilibration_method_.global);
      }
      else if (equilibration_method_.structure == CORE::LINALG::EquilibrationMethod::none and
               equilibration_method_.scatra == CORE::LINALG::EquilibrationMethod::none and
               equilibration_method_.thermo == CORE::LINALG::EquilibrationMethod::none)
      {
        equilibration_method_vector = std::vector<CORE::LINALG::EquilibrationMethod>(
            1, CORE::LINALG::EquilibrationMethod::none);
      }
      else
      {
        auto block_positions_scatra = GetBlockPositions(Subproblem::scalar_transport);
        auto block_position_structure = GetBlockPositions(Subproblem::structure);
        auto block_positions_thermo = GetBlockPositions(Subproblem::thermo);

        equilibration_method_vector = std::vector<CORE::LINALG::EquilibrationMethod>(
            block_positions_scatra.size() + block_position_structure.size() +
                block_positions_thermo.size(),
            CORE::LINALG::EquilibrationMethod::none);

        for (const int block_position_scatra : block_positions_scatra)
          equilibration_method_vector.at(block_position_scatra) = equilibration_method_.scatra;

        equilibration_method_vector.at(block_position_structure.at(0)) =
            equilibration_method_.structure;

        for (const int block_position_thermo : block_positions_thermo)
          equilibration_method_vector.at(block_position_thermo) = equilibration_method_.thermo;
      }

      break;
    }
    default:
    {
      FOUR_C_THROW("Invalid matrix type associated with system matrix field!");
      break;
    }
  }

  return equilibration_method_vector;
}

FOUR_C_NAMESPACE_CLOSE
