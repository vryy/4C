/*--------------------------------------------------------------------------*/
/*! \file
\brief monolithic scalar-structure interaction

\level 2

*--------------------------------------------------------------------------*/

#include "ssti_monolithic.H"

#include "ssti_algorithm.H"
#include "ssti_monolithic_assemble_strategy.H"
#include "ssti_monolithic_evaluate_OffDiag.H"
#include "ssti_utils.H"

#include "../drt_adapter/ad_str_ssiwrapper.H"
#include "../drt_adapter/ad_str_structure_new.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"

#include "../drt_io/io_control.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_scatra/scatra_timint_meshtying_strategy_s2i.H"

#include "../drt_ssi/ssi_monolithic_evaluate_OffDiag.H"

#include "../drt_sti/sti_monolithic_evaluate_OffDiag.H"

#include "../linalg/linalg_equilibrate.H"
#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_utils_sparse_algebra_create.H"
#include "../linalg/linalg_solver.H"

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
SSTI::SSTIMono::SSTIMono(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams)
    : SSTIAlgorithm(comm, globaltimeparams),
      increment_(Teuchos::null),
      itermax_(globaltimeparams.get<int>("ITEMAX")),
      itertol_(globaltimeparams.sublist("MONOLITHIC").get<double>("CONVTOL")),
      residual_(Teuchos::null),
      restol_(globaltimeparams.sublist("MONOLITHIC").get<double>("ABSTOLRES")),
      solver_(Teuchos::rcp(
          new LINALG::Solver(DRT::Problem::Instance()->SolverParams(
                                 globaltimeparams.sublist("MONOLITHIC").get<int>("LINEAR_SOLVER")),
              comm, DRT::Problem::Instance()->ErrorFile()->Handle()))),
      scatrastructureoffdiagcoupling_(Teuchos::null),
      scatrathermooffdiagcoupling_(Teuchos::null),
      thermostructureoffdiagcoupling_(Teuchos::null),
      dtassemble_(0.0),
      dtevaluate_(0.0),
      dtnewton_(0.0),
      dtsolve_(0.0),
      timer_(Teuchos::rcp(new Epetra_Time(comm))),
      equilibration_method_(Teuchos::getIntegralValue<LINALG::EquilibrationMethod>(
          globaltimeparams.sublist("MONOLITHIC"), "EQUILIBRATION")),
      matrixtype_(Teuchos::getIntegralValue<LINALG::MatrixType>(
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
  double starttime = timer_->WallTime();

  // assemble blocks of subproblems into system matrix
  strategy_assemble_->AssembleScatraDomain(
      ssti_matrices_->SystemMatrix(), ScaTraField()->SystemMatrixOperator());
  strategy_assemble_->AssembleStructureDomain(
      ssti_matrices_->SystemMatrix(), StructureField()->SystemMatrix());
  strategy_assemble_->AssembleThermoDomain(
      ssti_matrices_->SystemMatrix(), ThermoField()->SystemMatrixOperator());

  // assemble domain contributions from coupling into system matrix
  strategy_assemble_->AssembleScatraStructureDomain(
      ssti_matrices_->SystemMatrix(), ssti_matrices_->ScaTraStructureDomain());
  strategy_assemble_->AssembleStructureScatraDomain(
      ssti_matrices_->SystemMatrix(), ssti_matrices_->StructureScaTraDomain());
  strategy_assemble_->AssembleThermoStructureDomain(
      ssti_matrices_->SystemMatrix(), ssti_matrices_->ThermoStructureDomain());
  strategy_assemble_->AssembleStructureThermoDomain(
      ssti_matrices_->SystemMatrix(), ssti_matrices_->StructureThermoDomain());
  strategy_assemble_->AssembleThermoScatraDomain(
      ssti_matrices_->SystemMatrix(), ssti_matrices_->ThermoScaTraDomain());

  // assemble interface contributions from coupling into system matrix
  if (InterfaceMeshtying())
  {
    strategy_assemble_->AssembleScatraStructureInterface(
        ssti_matrices_->SystemMatrix(), ssti_matrices_->ScaTraStructureInterface());
    strategy_assemble_->AssembleThermoStructureInterface(
        ssti_matrices_->SystemMatrix(), ssti_matrices_->ThermoStructureInterface());
    strategy_assemble_->AssembleScatraThermoInterface(
        ssti_matrices_->SystemMatrix(), ssti_matrices_->ScaTraThermoInterface());
    strategy_assemble_->AssembleThermoScatraInterface(
        ssti_matrices_->SystemMatrix(), ssti_matrices_->ThermoScaTraIntefrace());
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

  double mydt = timer_->WallTime() - starttime;
  Comm().MaxAll(&mydt, &dtassemble_, 1);
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void SSTI::SSTIMono::BuildNullSpaces()
{
  // build null spaces for scatra and thermo
  switch (ScaTraField()->MatrixType())
  {
    case LINALG::MatrixType::block_condition:
    case LINALG::MatrixType::block_condition_dof:
    {
      ScaTraField()->BuildBlockNullSpaces(
          solver_, GetBlockPositions(Subproblem::scalar_transport)->at(0));
      ThermoField()->BuildBlockNullSpaces(solver_, GetBlockPositions(Subproblem::thermo)->at(0));
      break;
    }
    case LINALG::MatrixType::sparse:
    {
      // equip smoother for scatra matrix block with empty parameter sub lists to trigger null space
      // computation
      std::ostringstream scatrablockstr;
      scatrablockstr << GetBlockPositions(Subproblem::scalar_transport)->at(0) + 1;
      Teuchos::ParameterList& blocksmootherparamsscatra =
          solver_->Params().sublist("Inverse" + scatrablockstr.str());

      blocksmootherparamsscatra.sublist("Aztec Parameters");
      blocksmootherparamsscatra.sublist("MueLu Parameters");

      // equip smoother for scatra matrix block with null space associated with all degrees of
      // freedom on scatra discretization
      ScaTraField()->Discretization()->ComputeNullSpaceIfNecessary(blocksmootherparamsscatra);

      std::ostringstream thermoblockstr;
      thermoblockstr << GetBlockPositions(Subproblem::thermo)->at(0) + 1;
      Teuchos::ParameterList& blocksmootherparamsthermo =
          solver_->Params().sublist("Inverse" + thermoblockstr.str());
      blocksmootherparamsthermo.sublist("Aztec Parameters");
      blocksmootherparamsthermo.sublist("MueLu Parameters");

      // equip smoother for scatra matrix block with null space associated with all degrees of
      // freedom on scatra discretization
      ThermoField()->Discretization()->ComputeNullSpaceIfNecessary(blocksmootherparamsthermo);
      break;
    }
    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }
  // build null spaces for structure
  {
    // store number of matrix block associated with structural field as string
    std::stringstream iblockstr;
    iblockstr << GetBlockPositions(Subproblem::structure)->at(0) + 1;

    // equip smoother for structural matrix block with empty parameter sub lists to trigger null
    // space computation
    Teuchos::ParameterList& blocksmootherparams =
        solver_->Params().sublist("Inverse" + iblockstr.str());
    blocksmootherparams.sublist("Aztec Parameters");
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
  if (DRT::INPUT::IntegralValue<INPAR::SCATRA::VelocityField>(scatraparams, "VELOCITYFIELD") !=
      INPAR::SCATRA::velocity_Navier_Stokes)
    dserror("Invalid type of velocity field for scalar-structure interaction!");

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

  ScaTraField()->Output();
  ThermoField()->Output();
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
    dserror(
        "Since the ssti_monolithic framework is only implemented for usage in combination with "
        "volume change laws 'MAT_InelasticDefgradLinScalarIso' or  "
        "'MAT_InelasticDefgradLinScalarAniso' so far and these laws are implemented for only one "
        "transported scalar at the moment it is not reasonable to use them with more than one "
        "transported scalar. So you need to cope with it or change implementation! ;-)");
  }
  if (ScaTraField()->EquilibrationMethod() != LINALG::EquilibrationMethod::none)
  {
    dserror(
        "You are within the monolithic SSTI framework but activated a pure scatra equilibration "
        "method. Delete this from 'SCALAR TRANSPORT DYNAMIC' section and set it in 'SSTI "
        "CONTROL/MONOLITHIC' instead.");
  }
  if (!ScaTraField()->IsIncremental())
    dserror("Must have incremental solution approach for monolithic scalar-structure interaction!");
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSTI::SSTIMono::SetupSystem()
{
  // Setup all kind of maps
  ssti_maps_mono_ = Teuchos::rcp(new SSTI::SSTIMapsMono(*this));

  // initialize global increment vector for Newton-Raphson iteration
  increment_ = LINALG::CreateVector(*ssti_maps_mono_->MapsSubproblems()->FullMap(), true);

  // initialize global residual vector
  residual_ = LINALG::CreateVector(*ssti_maps_mono_->MapsSubproblems()->FullMap(), true);

  if (matrixtype_ == LINALG::MatrixType::block_field)
  {
    if (!solver_->Params().isSublist("AMGnxn Parameters"))
      dserror("Global system matrix with block structure requires AMGnxn block preconditioner!");

    // feed AMGnxn block preconditioner with null space information for each block of global
    // block system matrix
    BuildNullSpaces();
  }
  // setup interface maps with master and slace side dofs (fill only, if interface condition
  // available)
  Teuchos::RCP<Epetra_Map> interface_map_scatra(Teuchos::null);
  Teuchos::RCP<Epetra_Map> interface_map_thermo(Teuchos::null);
  Teuchos::RCP<LINALG::MultiMapExtractor> blockmapthermointerface(Teuchos::null);
  Teuchos::RCP<LINALG::MultiMapExtractor> blockmapthermointerfaceslave(Teuchos::null);
  Teuchos::RCP<LINALG::MultiMapExtractor> blockmapscatrainterface(Teuchos::null);

  if (InterfaceMeshtying())
  {
    interface_map_scatra = ssti_maps_mono_->MapInterface(MeshtyingScatra());
    interface_map_thermo = ssti_maps_mono_->MapInterface(MeshtyingThermo());

    switch (ScaTraField()->MatrixType())
    {
      case LINALG::MatrixType::block_condition:
      {
        blockmapscatrainterface = ssti_maps_mono_->MapsInterfaceBlocks(MeshtyingScatra(),
            LINALG::MatrixType::block_condition, ssti_maps_mono_->MapsScatra()->NumMaps());

        blockmapthermointerface = ssti_maps_mono_->MapsInterfaceBlocks(MeshtyingThermo(),
            LINALG::MatrixType::block_condition, ssti_maps_mono_->MapsThermo()->NumMaps());
        blockmapthermointerfaceslave = ssti_maps_mono_->MapsInterfaceBlocksSlave(MeshtyingThermo(),
            LINALG::MatrixType::block_condition, ssti_maps_mono_->MapsThermo()->NumMaps());
        break;
      }
      case LINALG::MatrixType::sparse:
      {
        blockmapthermointerface = ssti_maps_mono_->MapsInterfaceBlocks(MeshtyingThermo(),
            LINALG::MatrixType::sparse, ssti_maps_mono_->MapsThermo()->NumMaps());
        blockmapthermointerfaceslave = ssti_maps_mono_->MapsInterfaceBlocksSlave(MeshtyingThermo(),
            LINALG::MatrixType::sparse, ssti_maps_mono_->MapsThermo()->NumMaps());
        break;
      }
      default:
      {
        dserror("Invalid matrix type associated with scalar transport field!");
        break;
      }
    }
  }

  // initialize submatrices and system matrix
  ssti_matrices_ = Teuchos::rcp(new SSTI::SSTIMatrices(ssti_maps_mono_, matrixtype_,
      ScaTraField()->MatrixType(), interface_map_scatra, interface_map_thermo,
      blockmapscatrainterface, blockmapthermointerface, InterfaceMeshtying()));

  // initialize strategy for assembly
  ADAPTER::CouplingSlaveConverter converter(*CouplingAdapterStructure());
  strategy_assemble_ = SSTI::BuildAssembleStrategy(
      Teuchos::rcp(this, false), converter, matrixtype_, ScaTraField()->MatrixType());

  // initialize evaluation objects for coupling betwee subproblems
  scatrastructureoffdiagcoupling_ = Teuchos::rcp(
      new SSI::ScatraStructureOffDiagCoupling(ssti_maps_mono_->MapsInterfaceStructure(),
          ssti_maps_mono_->MapsSubproblems()->Map(GetProblemPosition(Subproblem::scalar_transport)),
          ssti_maps_mono_->MapsSubproblems()->Map(GetProblemPosition(Subproblem::structure)),
          CouplingAdapterStructure(), ssti_maps_mono_->MapInterface(MeshtyingScatra()),
          MeshtyingScatra(), ScaTraFieldBase(), StructureField()));

  thermostructureoffdiagcoupling_ = Teuchos::rcp(new SSTI::ThermoStructureOffDiagCoupling(
      ssti_maps_mono_->MapsInterfaceStructure(), ssti_maps_mono_->MapsThermo(),
      ssti_maps_mono_->MapsSubproblems()->Map(GetProblemPosition(Subproblem::structure)),
      ssti_maps_mono_->MapsSubproblems()->Map(GetProblemPosition(Subproblem::thermo)),
      CouplingAdapterStructure(), interface_map_thermo, MeshtyingThermo(), StructureField(),
      ThermoFieldBase()));

  scatrathermooffdiagcoupling_ = Teuchos::rcp(new STI::ScatraThermoOffDiagCouplingMatchingNodes(
      ssti_maps_mono_->MapsThermo(), blockmapthermointerface, blockmapthermointerfaceslave,
      ssti_maps_mono_->MapsSubproblems()->Map(GetProblemPosition(Subproblem::scalar_transport)),
      ssti_maps_mono_->MapsSubproblems()->Map(GetProblemPosition(Subproblem::thermo)),
      ssti_maps_mono_->MapInterface(MeshtyingScatra()), interface_map_thermo, true,
      MeshtyingScatra(), MeshtyingThermo(), ScaTraFieldBase(), ThermoFieldBase()));

  // initialize equilibration class
  strategy_equilibration_ = LINALG::BuildEquilibration(
      matrixtype_, equilibration_method_, AllMaps()->MapsSubproblems()->FullMap());
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSTI::SSTIMono::NewtonLoop()
{
  double starttime = timer_->WallTime();

  // initialize counter for Newton-Raphson iteration
  iter_ = 0;

  // start Newton-Raphson iteration
  while (true)
  {
    PrepareNewtonStep();

    EvaluateSubproblems();

    AssembleMatAndRHS();

    if (convcheck_->Converged(*this)) break;

    LinearSolve();

    UpdateIterStates();
  }

  double mydt = timer_->WallTime() - starttime;
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

    ScaTraField()->Output();
    ThermoField()->Output();
  }
  // time loop
  while (NotFinished() and ScaTraField()->NotFinished())
  {
    PrepareTimeStep();

    NewtonLoop();

    StructureField()->PrepareOutput();

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
      subincrement = ssti_maps_mono_->MapsSubproblems()->ExtractVector(
          increment_, GetProblemPosition(Subproblem::structure));

      // Second, copy master side displacements and increments to slave side for meshtying
      if (InterfaceMeshtying())
      {
        // displacements
        ssti_maps_mono_->MapsInterfaceStructure()->InsertVector(
            CouplingAdapterStructure()->MasterToSlave(
                ssti_maps_mono_->MapsInterfaceStructure()->ExtractVector(
                    StructureField()->Dispnp(), 1)),
            0, StructureField()->WriteAccessDispnp());

        // increments
        StructureField()->SetState(StructureField()->WriteAccessDispnp());
        ssti_maps_mono_->MapsInterfaceStructure()->InsertVector(
            CouplingAdapterStructure()->MasterToSlave(
                ssti_maps_mono_->MapsInterfaceStructure()->ExtractVector(subincrement, 1)),
            0, subincrement);
      }
      break;
    }
    case Subproblem::scalar_transport:
    {
      subincrement = ssti_maps_mono_->MapsSubproblems()->ExtractVector(
          increment_, GetProblemPosition(Subproblem::scalar_transport));
      break;
    }
    case Subproblem::thermo:
    {
      subincrement = ssti_maps_mono_->MapsSubproblems()->ExtractVector(
          increment_, GetProblemPosition(Subproblem::thermo));
      break;
    }
    default:
    {
      dserror("Unknown type of subproblem in SSTI");
      break;
    }
  }
  return subincrement;
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void SSTI::SSTIMono::EvaluateSubproblems()
{
  double starttime = timer_->WallTime();

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

  // evaluate interface contributions from coupling
  if (InterfaceMeshtying())
  {
    scatrastructureoffdiagcoupling_->EvaluateOffDiagBlockScatraStructureInterface(
        ssti_matrices_->ScaTraStructureInterface());
    thermostructureoffdiagcoupling_->EvaluateOffDiagBlockThermoStructureInterface(
        ssti_matrices_->ThermoStructureInterface());
    scatrathermooffdiagcoupling_->EvaluateOffDiagBlockThermoScatraInterface(
        ssti_matrices_->ThermoScaTraIntefrace());
    scatrathermooffdiagcoupling_->EvaluateOffDiagBlockScatraThermoInterface(
        ssti_matrices_->ScaTraThermoInterface());
  }

  double mydt = timer_->WallTime() - starttime;
  Comm().MaxAll(&mydt, &dtevaluate_, 1);
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void SSTI::SSTIMono::LinearSolve()
{
  double starttime = timer_->WallTime();

  increment_->PutScalar(0.0);

  if (!ssti_matrices_->SystemMatrix()->Filled())
    dserror("Complete() has not been called on global system matrix yet!");

  strategy_equilibration_->EquilibrateSystem(
      ssti_matrices_->SystemMatrix(), residual_, *AllMaps()->MapsSystemMatrixSubblocks());

  solver_->Solve(
      ssti_matrices_->SystemMatrix()->EpetraOperator(), increment_, residual_, true, iter_ == 1);

  strategy_equilibration_->UnequilibrateIncrement(increment_);

  double mydt = timer_->WallTime() - starttime;
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
  ++iter_;

  // reset timer
  timer_->ResetStartTime();

  ssti_matrices_->SystemMatrix()->Zero();
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
Teuchos::RCP<std::vector<int>> SSTI::SSTIMono::GetBlockPositions(Subproblem subproblem) const
{
  if (matrixtype_ == LINALG::MatrixType::sparse) dserror("Sparse matrices have just one block");

  Teuchos::RCP<std::vector<int>> block_position = Teuchos::rcp(new std::vector<int>(0));

  switch (subproblem)
  {
    case Subproblem::structure:
    {
      if (ScaTraField()->MatrixType() == LINALG::MatrixType::sparse)
        block_position->emplace_back(1);
      else
        block_position->emplace_back(ScaTraField()->BlockMaps().NumMaps());
      break;
    }
    case Subproblem::scalar_transport:
    {
      if (ScaTraField()->MatrixType() == LINALG::MatrixType::sparse)
        block_position->emplace_back(0);
      else

      {
        for (int i = 0; i < static_cast<int>(ScaTraField()->BlockMaps().NumMaps()); ++i)
          block_position->emplace_back(i);
      }
      break;
    }
    case Subproblem::thermo:
    {
      if (ThermoField()->MatrixType() == LINALG::MatrixType::sparse)
        block_position->emplace_back(2);
      else
      {
        for (int i = 0; i < static_cast<int>(ThermoField()->BlockMaps().NumMaps()); ++i)
          block_position->emplace_back(ScaTraField()->BlockMaps().NumMaps() + 1 + i);
      }
      break;
    }
    default:
    {
      dserror("Unknown type of subproblem");
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
      dserror("Unknown type of subproblem");
      break;
    }
  }

  return position;
}