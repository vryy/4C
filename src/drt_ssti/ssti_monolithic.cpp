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
      systemmatrix_(Teuchos::null),
      scatrastructuredomain_(Teuchos::null),
      scatrastructureinterface_(Teuchos::null),
      scatrathermodomain_(Teuchos::null),
      scatrathermointerface_(Teuchos::null),
      structurescatradomain_(Teuchos::null),
      structurethermodomain_(Teuchos::null),
      thermoscatradomain_(Teuchos::null),
      thermoscatrainterface_(Teuchos::null),
      thermostructuredomain_(Teuchos::null),
      thermostructureinterface_(Teuchos::null),
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
      ssti_maps_mono_(Teuchos::null),
      strategy_assemble_(Teuchos::null),
      convcheck_(Teuchos::rcp(new SSTI::ConvCheckMono(globaltimeparams)))
{
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSTI::SSTIMono::AssembleMatAndRHS()
{
  double starttime = timer_->WallTime();

  // assemble blocks of subproblems into system matrix
  strategy_assemble_->AssembleScatraDomain(systemmatrix_, ScaTraField()->SystemMatrixOperator());
  strategy_assemble_->AssembleStructureDomain(systemmatrix_, StructureField()->SystemMatrix());
  strategy_assemble_->AssembleThermoDomain(systemmatrix_, ThermoField()->SystemMatrixOperator());

  // assemble domain contributions from coupling into system matrix
  strategy_assemble_->AssembleScatraStructureDomain(systemmatrix_, scatrastructuredomain_);
  strategy_assemble_->AssembleStructureScatraDomain(systemmatrix_, structurescatradomain_);
  strategy_assemble_->AssembleThermoStructureDomain(systemmatrix_, thermostructuredomain_);
  strategy_assemble_->AssembleStructureThermoDomain(systemmatrix_, structurethermodomain_);
  strategy_assemble_->AssembleThermoScatraDomain(systemmatrix_, thermoscatradomain_);

  // assemble interface contributions from coupling into system matrix
  if (InterfaceMeshtying())
  {
    strategy_assemble_->AssembleScatraStructureInterface(systemmatrix_, scatrastructureinterface_);
    strategy_assemble_->AssembleThermoStructureInterface(systemmatrix_, thermostructureinterface_);
    strategy_assemble_->AssembleScatraThermoInterface(systemmatrix_, scatrathermointerface_);
    strategy_assemble_->AssembleThermoScatraInterface(systemmatrix_, thermoscatrainterface_);
  }

  // apply meshtying on structural linearizations
  strategy_assemble_->ApplyMeshtyingSystemMatrix(systemmatrix_);

  // finalize global system matrix
  systemmatrix_->Complete();

  // apply Dirichlet conditions
  systemmatrix_->ApplyDirichlet(*ScaTraField()->DirichMaps()->CondMap(), true);
  systemmatrix_->ApplyDirichlet(*ThermoField()->DirichMaps()->CondMap(), true);
  strategy_assemble_->ApplyStructuralDBCSystemMatrix(systemmatrix_);

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
      ScaTraField()->BuildBlockNullSpaces(solver_, 0);
      ThermoField()->BuildBlockNullSpaces(solver_, ssti_maps_mono_->MapsScatra()->NumMaps() + 1);
      break;
    }
    case LINALG::MatrixType::sparse:
    {
      // equip smoother for scatra matrix block with empty parameter sub lists to trigger null space
      // computation
      Teuchos::ParameterList& blocksmootherparams = solver_->Params().sublist("Inverse1");
      blocksmootherparams.sublist("Aztec Parameters");
      blocksmootherparams.sublist("MueLu Parameters");

      // equip smoother for scatra matrix block with null space associated with all degrees of
      // freedom on scatra discretization
      ScaTraField()->Discretization()->ComputeNullSpaceIfNecessary(blocksmootherparams);

      Teuchos::ParameterList& blocksmootherparamsthermo = solver_->Params().sublist("Inverse3");
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
    iblockstr << ssti_maps_mono_->MapsScatra()->NumMaps() + 1;

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
  ssti_maps_mono_ = Teuchos::rcp(new SSTI::SSTIMapsMono(StructureField(), ScaTraField(),
      ThermoField(), CouplingAdapterStructure(), InterfaceMeshtying()));

  // initialize global increment vector for Newton-Raphson iteration
  increment_ = LINALG::CreateVector(*ssti_maps_mono_->MapsSubproblems()->FullMap(), true);

  // initialize global residual vector
  residual_ = LINALG::CreateVector(*ssti_maps_mono_->MapsSubproblems()->FullMap(), true);

  // perform initializations associated with global system matrix
  switch (matrixtype_)
  {
    case LINALG::MatrixType::block_field:
    {
      if (!solver_->Params().isSublist("AMGnxn Parameters"))
        dserror("Global system matrix with block structure requires AMGnxn block preconditioner!");

      // initialize global system matrix
      systemmatrix_ =
          Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
              *ssti_maps_mono_->MapsSystemMatrixSubblocks(),
              *ssti_maps_mono_->MapsSystemMatrixSubblocks(), 81, false, true));

      // feed AMGnxn block preconditioner with null space information for each block of global
      // block system matrix
      BuildNullSpaces();

      break;
    }

    case LINALG::MatrixType::sparse:
    {
      if (ScaTraField()->SystemMatrix() == Teuchos::null)
        dserror("Incompatible matrix type associated with scalar transport field!");

      // initialize global system matrix
      systemmatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(
          *ssti_maps_mono_->MapsSubproblems()->FullMap(), 27, false, true));

      break;
    }

    default:
    {
      dserror("Type of global system matrix for scalar-structure interaction not recognized!");
      break;
    }
  }

  // setup interface maps with master and slace side dofs
  Teuchos::RCP<Epetra_Map> interface_map_scatra = ssti_maps_mono_->MapInterface(MeshtyingScatra());
  Teuchos::RCP<Epetra_Map> interface_map_thermo = ssti_maps_mono_->MapInterface(MeshtyingThermo());
  Teuchos::RCP<LINALG::MultiMapExtractor> blockmapthermointerface(Teuchos::null);
  Teuchos::RCP<LINALG::MultiMapExtractor> blockmapthermointerfaceslave(Teuchos::null);

  switch (ScaTraField()->MatrixType())
  {
    case LINALG::MatrixType::block_condition:
    {
      blockmapthermointerface = ssti_maps_mono_->MapsInterfaceBlocks(MeshtyingThermo(),
          LINALG::MatrixType::block_condition, ssti_maps_mono_->MapsThermo()->NumMaps());
      blockmapthermointerfaceslave = ssti_maps_mono_->MapsInterfaceBlocksSlave(MeshtyingThermo(),
          LINALG::MatrixType::block_condition, ssti_maps_mono_->MapsThermo()->NumMaps());
      break;
    }
    case LINALG::MatrixType::sparse:
    {
      blockmapthermointerface = ssti_maps_mono_->MapsInterfaceBlocks(
          MeshtyingThermo(), LINALG::MatrixType::sparse, ssti_maps_mono_->MapsThermo()->NumMaps());
      blockmapthermointerfaceslave = ssti_maps_mono_->MapsInterfaceBlocksSlave(
          MeshtyingThermo(), LINALG::MatrixType::sparse, ssti_maps_mono_->MapsThermo()->NumMaps());
      break;
    }
    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // setup blocks for coupling matrices
  switch (ScaTraField()->MatrixType())
  {
    case LINALG::MatrixType::block_condition:
    {
      Teuchos::RCP<LINALG::MultiMapExtractor> blockmapscatrainterface =
          ssti_maps_mono_->MapsInterfaceBlocks(MeshtyingScatra(),
              LINALG::MatrixType::block_condition, ssti_maps_mono_->MapsScatra()->NumMaps());

      Teuchos::RCP<LINALG::MultiMapExtractor> blockmapthermointerface =
          ssti_maps_mono_->MapsInterfaceBlocks(MeshtyingThermo(),
              LINALG::MatrixType::block_condition, ssti_maps_mono_->MapsThermo()->NumMaps());

      scatrastructuredomain_ =
          Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
              *ssti_maps_mono_->MapsStructure(), ScaTraField()->BlockMaps(), 81, false, true));

      structurescatradomain_ =
          Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
              ScaTraField()->BlockMaps(), *ssti_maps_mono_->MapsStructure(), 81, false, true));

      structurethermodomain_ =
          Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
              ThermoField()->BlockMaps(), *ssti_maps_mono_->MapsStructure(), 81, false, true));

      thermostructuredomain_ =
          Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
              *ssti_maps_mono_->MapsStructure(), ThermoField()->BlockMaps(), 81, false, true));

      scatrathermodomain_ =
          Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
              ThermoField()->BlockMaps(), ScaTraField()->BlockMaps(), 81, false, true));

      thermoscatradomain_ =
          Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
              ScaTraField()->BlockMaps(), ThermoField()->BlockMaps(), 81, false, true));

      if (InterfaceMeshtying())
      {
        scatrastructureinterface_ =
            Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
                *ssti_maps_mono_->MapsStructure(), *blockmapscatrainterface, 81, false, true));
        thermostructureinterface_ =
            Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
                *ssti_maps_mono_->MapsStructure(), *blockmapthermointerface, 81, false, true));

        scatrathermointerface_ =
            Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
                *ssti_maps_mono_->MapsThermo(), *blockmapscatrainterface, 81, false, true));

        thermoscatrainterface_ =
            Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
                *ssti_maps_mono_->MapsScatra(), *blockmapthermointerface, 81, false, true));
      }
      break;
    }
    case LINALG::MatrixType::sparse:
    {
      scatrastructuredomain_ =
          Teuchos::rcp(new LINALG::SparseMatrix(*ScaTraField()->DofRowMap(), 27, false, true));

      structurescatradomain_ =
          Teuchos::rcp(new LINALG::SparseMatrix(*StructureField()->DofRowMap(), 27, false, true));

      structurethermodomain_ =
          Teuchos::rcp(new LINALG::SparseMatrix(*StructureField()->DofRowMap(), 27, false, true));

      thermostructuredomain_ =
          Teuchos::rcp(new LINALG::SparseMatrix(*ThermoField()->DofRowMap(), 27, false, true));

      scatrathermodomain_ =
          Teuchos::rcp(new LINALG::SparseMatrix(*ScaTraField()->DofRowMap(), 27, false, true));

      thermoscatradomain_ =
          Teuchos::rcp(new LINALG::SparseMatrix(*ThermoField()->DofRowMap(), 27, false, true));

      if (InterfaceMeshtying())
      {
        scatrastructureinterface_ =
            Teuchos::rcp(new LINALG::SparseMatrix(*interface_map_scatra, 27, false, true));

        thermostructureinterface_ =
            Teuchos::rcp(new LINALG::SparseMatrix(*interface_map_thermo, 27, false, true));

        scatrathermointerface_ =
            Teuchos::rcp(new LINALG::SparseMatrix(*interface_map_scatra, 27, false, true));

        thermoscatrainterface_ =
            Teuchos::rcp(new LINALG::SparseMatrix(*interface_map_thermo, 27, false, true));
      }
      break;
    }
    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // initialize strategy for assembly
  ADAPTER::CouplingSlaveConverter converter(*CouplingAdapterStructure());
  strategy_assemble_ = SSTI::BuildAssembleStrategy(
      Teuchos::rcp(this, false), converter, matrixtype_, ScaTraField()->MatrixType());

  // initialize evaluation objects for coupling betwee subproblems
  scatrastructureoffdiagcoupling_ = Teuchos::rcp(
      new SSI::ScatraStructureOffDiagCoupling(ssti_maps_mono_->MapsInterfaceStructure(),
          ssti_maps_mono_->MapsSubproblems()->Map(0), ssti_maps_mono_->MapsSubproblems()->Map(1),
          CouplingAdapterStructure(), ssti_maps_mono_->MapInterface(MeshtyingScatra()),
          MeshtyingScatra(), ScaTraFieldBase(), StructureField()));

  thermostructureoffdiagcoupling_ = Teuchos::rcp(
      new SSTI::ThermoStructureOffDiagCoupling(ssti_maps_mono_->MapsInterfaceStructure(),
          ssti_maps_mono_->MapsThermo(), ssti_maps_mono_->MapsSubproblems()->Map(1),
          ssti_maps_mono_->MapsSubproblems()->Map(2), CouplingAdapterStructure(),
          interface_map_thermo, MeshtyingThermo(), StructureField(), ThermoFieldBase()));

  scatrathermooffdiagcoupling_ = Teuchos::rcp(new STI::ScatraThermoOffDiagCouplingMatchingNodes(
      ssti_maps_mono_->MapsThermo(), blockmapthermointerface, blockmapthermointerfaceslave,
      ssti_maps_mono_->MapsSubproblems()->Map(0), ssti_maps_mono_->MapsSubproblems()->Map(2),
      ssti_maps_mono_->MapInterface(MeshtyingScatra()), interface_map_thermo, true,
      MeshtyingScatra(), MeshtyingThermo(), ScaTraFieldBase(), ThermoFieldBase()));

  // initialize equilibration class
  equilibration_ = LINALG::BuildEquilibration(
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
      subincrement = ssti_maps_mono_->MapsSubproblems()->ExtractVector(increment_, 1);

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
      subincrement = ssti_maps_mono_->MapsSubproblems()->ExtractVector(increment_, 0);
      break;
    }
    case Subproblem::thermo:
    {
      subincrement = ssti_maps_mono_->MapsSubproblems()->ExtractVector(increment_, 2);
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

  // consider order of 'set state'-calls. Velocity is updated inside Structure

  DistributeScatraSolution();
  DistributeThermoSolution();

  StructureField()->Evaluate();

  DistributeStructureSolution();
  ScaTraField()->PrepareLinearSolve();
  ThermoField()->PrepareLinearSolve();

  // evaluate domain contributions from coupling
  scatrastructureoffdiagcoupling_->EvaluateOffDiagBlockScatraStructureDomain(
      scatrastructuredomain_);
  scatrastructureoffdiagcoupling_->EvaluateOffDiagBlockStructureScatraDomain(
      structurescatradomain_);
  thermostructureoffdiagcoupling_->EvaluateOffDiagBlockThermoStructureDomain(
      thermostructuredomain_);
  thermostructureoffdiagcoupling_->EvaluateOffDiagBlockStructureThermoDomain(
      structurethermodomain_);
  scatrathermooffdiagcoupling_->EvaluateOffDiagBlockThermoScatraDomain(thermoscatradomain_);

  // evaluate interface contributions from coupling
  if (InterfaceMeshtying())
  {
    scatrastructureoffdiagcoupling_->EvaluateOffDiagBlockScatraStructureInterface(
        scatrastructureinterface_);
    thermostructureoffdiagcoupling_->EvaluateOffDiagBlockThermoStructureInterface(
        thermostructureinterface_);
    scatrathermooffdiagcoupling_->EvaluateOffDiagBlockThermoScatraInterface(thermoscatrainterface_);
    scatrathermooffdiagcoupling_->EvaluateOffDiagBlockScatraThermoInterface(scatrathermointerface_);
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

  if (!systemmatrix_->Filled())
    dserror("Complete() has not been called on global system matrix yet!");

  equilibration_->EquilibrateSystem(
      systemmatrix_, residual_, *AllMaps()->MapsSystemMatrixSubblocks());

  solver_->Solve(systemmatrix_->EpetraOperator(), increment_, residual_, true, iter_ == 1);

  equilibration_->UnequilibrateIncrement(increment_);

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

  systemmatrix_->Zero();
}
