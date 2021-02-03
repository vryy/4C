/*--------------------------------------------------------------------------*/
/*! \file
\brief monolithic scalar-structure interaction

\level 2

*/
/*--------------------------------------------------------------------------*/
#include "ssi_monolithic.H"

#include "ssi_coupling.H"
#include "ssi_monolithic_assemble_strategy.H"
#include "ssi_monolithic_convcheck_strategies.H"
#include "ssi_monolithic_evaluate_OffDiag.H"
#include "ssi_resulttest.H"
#include "ssi_str_model_evaluator_monolithic.H"
#include "ssi_utils.H"

#include "../drt_adapter/ad_str_ssiwrapper.H"
#include "../drt_adapter/ad_str_structure_new.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"

#include "../drt_inpar/inpar_scatra.H"

#include "../drt_io/io_control.H"

#include "../drt_lib/drt_assemblestrategy.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_locsys.H"

#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_scatra/scatra_timint_meshtying_strategy_s2i.H"

#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_matrixtransform.H"
#include "../linalg/linalg_equilibrate.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils_sparse_algebra_assemble.H"
#include "../linalg/linalg_utils_sparse_algebra_create.H"
#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"

#include <Epetra_Time.h>


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
SSI::SSIMono::SSIMono(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams)
    : SSIBase(comm, globaltimeparams),
      dtele_(0.0),
      dtsolve_(0.0),
      equilibration_method_{Teuchos::getIntegralValue<LINALG::EquilibrationMethod>(
                                globaltimeparams.sublist("MONOLITHIC"), "EQUILIBRATION"),
          Teuchos::getIntegralValue<LINALG::EquilibrationMethod>(
              globaltimeparams.sublist("MONOLITHIC"), "EQUILIBRATION_SCATRA"),
          Teuchos::getIntegralValue<LINALG::EquilibrationMethod>(
              globaltimeparams.sublist("MONOLITHIC"), "EQUILIBRATION_STRUCTURE")},
      increment_(Teuchos::null),
      map_structure_(Teuchos::null),
      maps_scatra_(Teuchos::null),
      maps_sub_problems_(Teuchos::null),
      maps_systemmatrix_(Teuchos::null),
      matrixtype_(Teuchos::getIntegralValue<LINALG::MatrixType>(
          globaltimeparams.sublist("MONOLITHIC"), "MATRIXTYPE")),
      meshtying_strategy_s2i_(Teuchos::null),
      residual_(Teuchos::null),
      scatrastructureOffDiagcoupling_(Teuchos::null),
      solver_(Teuchos::rcp(
          new LINALG::Solver(DRT::Problem::Instance()->SolverParams(
                                 globaltimeparams.sublist("MONOLITHIC").get<int>("LINEAR_SOLVER")),
              comm, DRT::Problem::Instance()->ErrorFile()->Handle()))),
      ssi_matrices_(Teuchos::null),
      strategy_assemble_(Teuchos::null),
      strategy_convcheck_(Teuchos::null),
      strategy_equilibration_(Teuchos::null),
      timer_(Teuchos::rcp(new Epetra_Time(comm)))

{
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SSIMono::AssembleMatAndRHS()
{
  // assemble scatra block into system matrix
  strategy_assemble_->AssembleScatraDomain(
      ssi_matrices_->SystemMatrix(), ScaTraField()->SystemMatrixOperator());

  // assemble scatra-strucutre block (domain contributions) into system matrix
  strategy_assemble_->AssembleScatraStructureDomain(
      ssi_matrices_->SystemMatrix(), ssi_matrices_->ScaTraStructureDomain());

  // assemble scatra-strucutre block (interface contributions) into system matrix
  if (SSIInterfaceMeshtying())
    strategy_assemble_->AssembleScatraStructureInterface(
        ssi_matrices_->SystemMatrix(), ssi_matrices_->ScaTraStructureInterface());

  // assemble structure-scatra block (domain contributions) into system matrix
  strategy_assemble_->AssembleStructureScatraDomain(
      ssi_matrices_->SystemMatrix(), ssi_matrices_->StructureScaTraDomain());

  // assemble structure block into system matrix
  strategy_assemble_->AssembleStructureDomain(
      ssi_matrices_->SystemMatrix(), StructureField()->SystemMatrix());

  if (IsScaTraManifold())
  {
    // assemble manifold block into system matrix
    strategy_assemble_->AssembleScaTraManifoldDomain(
        ssi_matrices_->SystemMatrix(), ScaTraManifold()->SystemMatrixOperator());

    // assemble manifold-structure block into system matrix
    strategy_assemble_->AssembleScaTraManifoldStructuredDomain(
        ssi_matrices_->SystemMatrix(), ssi_matrices_->ScaTraManifoldStructureDomain());
  }

  // apply meshtying
  strategy_assemble_->ApplyMeshtyingSystemMatrix(ssi_matrices_->SystemMatrix());

  // finalize global system matrix
  ssi_matrices_->SystemMatrix()->Complete();

  // apply scatra Dirichlet
  ssi_matrices_->SystemMatrix()->ApplyDirichlet(*ScaTraField()->DirichMaps()->CondMap(), true);

  // apply structural Dirichlet conditions
  strategy_assemble_->ApplyStructuralDBCSystemMatrix(ssi_matrices_->SystemMatrix());

  // assemble monolithic RHS
  strategy_assemble_->AssembleRHS(residual_, ScaTraField()->Residual(), StructureField()->RHS(),
      IsScaTraManifold() ? ScaTraManifold()->Residual() : Teuchos::null);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SSIMono::EvaluateSubproblems()
{
  // needed to communicate to NOX state
  StructureField()->SetState(StructureField()->WriteAccessDispnp());

  // pass structural degrees of freedom to scalar transport discretization
  SetStructSolution(StructureField()->Dispnp(), StructureField()->Velnp());

  // pass scalar transport degrees of freedom to structural discretization
  SetScatraSolution(ScaTraField()->Phinp());

  // evaluate temperature from function and set to structural discretization
  EvaluateAndSetTemperatureField();

  // build system matrix and residual for structure field
  StructureField()->Evaluate();

  // build system matrix and residual for scalar transport field
  ScaTraField()->PrepareLinearSolve();

  // build system matrix and residual for scalar transport field on manifold
  if (IsScaTraManifold()) ScaTraManifold()->PrepareLinearSolve();

  // evaluate off-diagonal scatra-structure block (domain contributions) of global system matrix
  scatrastructureOffDiagcoupling_->EvaluateOffDiagBlockScatraStructureDomain(
      ssi_matrices_->ScaTraStructureDomain());

  // evaluate off-diagonal scatra-structure block (interface contributions) of global system matrix
  if (SSIInterfaceMeshtying())
    scatrastructureOffDiagcoupling_->EvaluateOffDiagBlockScatraStructureInterface(
        ssi_matrices_->ScaTraStructureInterface());

  // evaluate off-diagonal structure-scatra block (we only have domain contributions so far) of
  // global system matrix
  scatrastructureOffDiagcoupling_->EvaluateOffDiagBlockStructureScatraDomain(
      ssi_matrices_->StructureScaTraDomain());

  if (IsScaTraManifold())
  {
    // evaluate off-diagonal manifold-structure block of global system matrix
    scatrastructureOffDiagcoupling_->EvaluateOffDiagBlockScatraManifoldStructureDomain(
        ssi_matrices_->ScaTraManifoldStructureDomain());
  }
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void SSI::SSIMono::BuildNullSpaces() const
{
  switch (ScaTraField()->MatrixType())
  {
    case LINALG::MatrixType::block_condition:
    case LINALG::MatrixType::block_condition_dof:
    {
      // equip smoother for scatra matrix blocks with null space
      ScaTraField()->BuildBlockNullSpaces(
          solver_, GetBlockPositions(Subproblem::scalar_transport)->at(0));
      if (IsScaTraManifold())
      {
        ScaTraManifold()->BuildBlockNullSpaces(
            solver_, GetBlockPositions(Subproblem::manifold)->at(0));
      }
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

      if (IsScaTraManifold())
      {
        std::ostringstream scatramanifoldblockstr;
        scatramanifoldblockstr << GetBlockPositions(Subproblem::manifold)->at(0) + 1;
        Teuchos::ParameterList& blocksmootherparamsscatramanifold =
            solver_->Params().sublist("Inverse" + scatramanifoldblockstr.str());
        blocksmootherparamsscatramanifold.sublist("Aztec Parameters");
        blocksmootherparamsscatramanifold.sublist("MueLu Parameters");

        // equip smoother for scatra matrix block with null space associated with all degrees of
        // freedom on scatra discretization
        ScaTraManifold()->Discretization()->ComputeNullSpaceIfNecessary(
            blocksmootherparamsscatramanifold);
      }

      break;
    }

    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // store number of matrix block associated with structural field as string
  std::stringstream iblockstr;
  iblockstr << GetBlockPositions(Subproblem::structure)->at(0) + 1;

  // equip smoother for structural matrix block with empty parameter sub lists to trigger null space
  // computation
  Teuchos::ParameterList& blocksmootherparams =
      solver_->Params().sublist("Inverse" + iblockstr.str());
  blocksmootherparams.sublist("Aztec Parameters");
  blocksmootherparams.sublist("MueLu Parameters");

  // equip smoother for structural matrix block with null space associated with all degrees of
  // freedom on structural discretization
  StructureField()->Discretization()->ComputeNullSpaceIfNecessary(blocksmootherparams);
}  // SSI::SSIMono::BuildNullSpaces


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
const Teuchos::RCP<const Epetra_Map>& SSI::SSIMono::DofRowMap() const
{
  return MapsSubProblems()->FullMap();
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SSIMono::Init(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& scatraparams, const Teuchos::ParameterList& structparams,
    const std::string& struct_disname, const std::string& scatra_disname, bool isAle)
{
  // check input parameters for scalar transport field
  if (DRT::INPUT::IntegralValue<INPAR::SCATRA::VelocityField>(scatraparams, "VELOCITYFIELD") !=
      INPAR::SCATRA::velocity_Navier_Stokes)
    dserror("Invalid type of velocity field for scalar-structure interaction!");

  // initialize strategy for Newton-Raphson convergence check
  switch (
      Teuchos::getIntegralValue<INPAR::SSI::ScaTraTimIntType>(globaltimeparams, "SCATRATIMINTTYPE"))
  {
    case INPAR::SSI::ScaTraTimIntType::elch:
    {
      if (IsScaTraManifold())
      {
        strategy_convcheck_ =
            Teuchos::rcp(new SSI::SSIMono::ConvCheckStrategyElchScaTraManifold(globaltimeparams));
      }
      else
        strategy_convcheck_ =
            Teuchos::rcp(new SSI::SSIMono::ConvCheckStrategyElch(globaltimeparams));
      break;
    }

    case INPAR::SSI::ScaTraTimIntType::standard:
    {
      strategy_convcheck_ = Teuchos::rcp(new SSI::SSIMono::ConvCheckStrategyStd(globaltimeparams));
      break;
    }

    default:
    {
      dserror("Type of scalar transport time integrator currently not supported!");
      break;
    }
  }

  // call base class routine
  SSIBase::Init(
      comm, globaltimeparams, scatraparams, structparams, struct_disname, scatra_disname, isAle);
}


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SSIMono::Output()
{
  // output scalar transport field
  ScaTraField()->Output();
  if (IsScaTraManifold()) ScaTraManifold()->Output();

  // output structure field
  StructureField()->Output();
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SSIMono::PrepareTimeStep()
{
  // update time and time step
  IncrementTimeAndStep();

  // pass structural degrees of freedom to scalar transport discretization
  SetStructSolution(StructureField()->Dispnp(), StructureField()->Velnp());

  // prepare time step for scalar transport field
  ScaTraField()->PrepareTimeStep();
  if (IsScaTraManifold()) ScaTraManifold()->PrepareTimeStep();

  // if adaptive time stepping and different time step size: calculate time step in scatra
  // (PrepareTimeStep() of Scatra) and pass to structure
  if (ScaTraField()->TimeStepAdapted()) SetDtFromScaTraToStructure();

  // pass scalar transport degrees of freedom to structural discretization
  // has to be called AFTER ScaTraField()->PrepareTimeStep() to ensure
  // consistent scalar transport state vector with valid Dirichlet conditions
  SetScatraSolution(ScaTraField()->Phinp());

  // evaluate temperature from function and set to structural discretization
  EvaluateAndSetTemperatureField();

  // prepare time step for structural field
  StructureField()->PrepareTimeStep();

  // print time step information to screen
  ScaTraField()->PrintTimeStepInfo();
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SSIMono::Setup()
{
  // call base class routine
  SSIBase::Setup();

  // safety checks
  if (ScaTraField()->NumScal() != 1)
  {
    dserror(
        "Since the ssi_monolithic framework is only implemented for usage in combination with "
        "volume change laws 'MAT_InelasticDefgradLinScalarIso' or "
        "'MAT_InelasticDefgradLinScalarAniso' so far and these laws are implemented for only "
        "one transported scalar at the moment it is not reasonable to use them with more than one "
        "transported scalar. So you need to cope with it or change implementation! ;-)");
  }

  const bool equilibration_scatra_initial = DRT::INPUT::IntegralValue<bool>(
      DRT::Problem::Instance()->SSIControlParams().sublist("MONOLITHIC"),
      "EQUILIBRATION_INIT_SCATRA");
  const bool calc_initial_pot =
      DRT::INPUT::IntegralValue<bool>(DRT::Problem::Instance()->ELCHControlParams(), "INITPOTCALC");

  if (!equilibration_scatra_initial and
      ScaTraField()->EquilibrationMethod() != LINALG::EquilibrationMethod::none)
  {
    dserror(
        "You are within the monolithic solid scatra interaction framework but activated a pure "
        "scatra equilibration method. Delete this from 'SCALAR TRANSPORT DYNAMIC' section and set "
        "it in 'SSI CONTROL/MONOLITHIC' instead.");
  }
  if (equilibration_scatra_initial and
      ScaTraField()->EquilibrationMethod() == LINALG::EquilibrationMethod::none)
  {
    dserror(
        "You selected to equilibrate equations of initial potential but did not specify any "
        "equilibration method in ScaTra.");
  }
  if (equilibration_scatra_initial and !calc_initial_pot)
  {
    dserror(
        "You selected to equilibrate equations of initial potential but did not activate "
        "INITPOTCALC in ELCH CONTROL");
  }

  if (equilibration_method_.global != LINALG::EquilibrationMethod::local and
      (equilibration_method_.structure != LINALG::EquilibrationMethod::none or
          equilibration_method_.scatra != LINALG::EquilibrationMethod::none))
    dserror("Either global equilibration or local equilibration");

  if (matrixtype_ == LINALG::MatrixType::sparse and
      (equilibration_method_.structure != LINALG::EquilibrationMethod::none or
          equilibration_method_.scatra != LINALG::EquilibrationMethod::none))
    dserror("Block based equilibration only for block matrices");

  if (!ScaTraField()->IsIncremental())
    dserror("Must have incremental solution approach for monolithic scalar-structure interaction!");

  // set up scatra-scatra interface meshtying if set in input file
  if (SSIInterfaceMeshtying())
  {
    // extract meshtying strategy for scatra-scatra interface coupling on scatra discretization
    meshtying_strategy_s2i_ =
        Teuchos::rcp_dynamic_cast<const SCATRA::MeshtyingStrategyS2I>(ScaTraField()->Strategy());

    // safety checks
    if (meshtying_strategy_s2i_ == Teuchos::null)
      dserror("Invalid scatra-scatra interface coupling strategy!");
    if (meshtying_strategy_s2i_->CouplingType() != INPAR::S2I::coupling_matching_nodes)
    {
      dserror(
          "Monolithic scalar-structure interaction only implemented for scatra-scatra "
          "interface coupling with matching interface nodes!");
    }
  }
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SSIMono::SetupSystem()
{
  // merge slave and master side block maps for interface matrix for scatra
  Teuchos::RCP<Epetra_Map> interface_map_scatra(Teuchos::null);

  if (SSIInterfaceMeshtying())
  {
    // check whether slave-side degrees of freedom are Dirichlet-free
    std::vector<Teuchos::RCP<const Epetra_Map>> maps(2, Teuchos::null);
    maps[0] = InterfaceCouplingAdapterStructure()->SlaveDofMap();
    maps[1] = StructureField()->GetDBCMapExtractor()->CondMap();
    if (LINALG::MultiMapExtractor::IntersectMaps(maps)->NumGlobalElements() > 0)
      dserror("Must not apply Dirichlet conditions to slave-side structural displacements!");

    interface_map_scatra = LINALG::MultiMapExtractor::MergeMaps(
        {meshtying_strategy_s2i_->CouplingAdapter()->MasterDofMap(),
            meshtying_strategy_s2i_->CouplingAdapter()->SlaveDofMap()});
  }

  // initialize global map extractor
  std::vector<Teuchos::RCP<const Epetra_Map>> partial_maps;
  Teuchos::RCP<const Epetra_Map> merged_map;
  partial_maps.emplace_back(Teuchos::rcp(new Epetra_Map(*ScaTraField()->DofRowMap())));
  partial_maps.emplace_back(Teuchos::rcp(new Epetra_Map(*StructureField()->DofRowMap())));
  if (IsScaTraManifold())
  {
    partial_maps.emplace_back(Teuchos::rcp(new Epetra_Map(*ScaTraManifold()->DofRowMap())));
    auto temp_map = LINALG::MergeMap(partial_maps[0], partial_maps[1], false);
    merged_map = LINALG::MergeMap(temp_map, partial_maps[2], false);
  }
  else
    merged_map = LINALG::MergeMap(partial_maps[0], partial_maps[1], false);

  maps_sub_problems_ = Teuchos::rcp(new LINALG::MultiMapExtractor(*merged_map, partial_maps));
  // check global map extractor
  maps_sub_problems_->CheckForValidMapExtractor();

  // initialize global increment vector for Newton-Raphson iteration
  increment_ = LINALG::CreateVector(*DofRowMap(), true);

  // initialize global residual vector
  residual_ = LINALG::CreateVector(*DofRowMap(), true);

  // initialize map extractors associated with blocks of global system matrix
  switch (ScaTraField()->MatrixType())
  {
    // one single main-diagonal matrix block associated with scalar transport field
    case LINALG::MatrixType::sparse:
    {
      maps_systemmatrix_ = MapsSubProblems();
      break;
    }

    // several main-diagonal matrix blocks associated with scalar transport field
    case LINALG::MatrixType::block_condition:
    case LINALG::MatrixType::block_condition_dof:
    {
      std::vector<Teuchos::RCP<const Epetra_Map>> maps_systemmatrix;

      // store an RCP to the block maps of the scatra field
      maps_scatra_ = Teuchos::rcpFromRef(ScaTraField()->BlockMaps());
      maps_scatra_->CheckForValidMapExtractor();

      if (IsScaTraManifold())
      {
        auto maps_scatra_manifold = Teuchos::rcpFromRef(ScaTraManifold()->BlockMaps());
        maps_scatra_manifold->CheckForValidMapExtractor();
        maps_systemmatrix.resize(GetBlockPositions(Subproblem::scalar_transport)->size() +
                                 GetBlockPositions(Subproblem::structure)->size() +
                                 GetBlockPositions(Subproblem::manifold)->size());

        for (int imap = 0; imap < static_cast<int>(GetBlockPositions(Subproblem::manifold)->size());
             ++imap)
        {
          maps_systemmatrix[GetBlockPositions(Subproblem::manifold)->at(imap)] =
              maps_scatra_manifold->Map(imap);
        }
      }
      else
      {
        // extract maps underlying main-diagonal matrix blocks associated with scalar transport
        // field
        maps_systemmatrix.resize(GetBlockPositions(Subproblem::scalar_transport)->size() +
                                 GetBlockPositions(Subproblem::structure)->size());
      }

      for (int imap = 0;
           imap < static_cast<int>(GetBlockPositions(Subproblem::scalar_transport)->size()); ++imap)
      {
        maps_systemmatrix[GetBlockPositions(Subproblem::scalar_transport)->at(imap)] =
            maps_scatra_->Map(imap);
      }

      // extract map underlying single main-diagonal matrix block associated with structural field
      maps_systemmatrix[GetBlockPositions(Subproblem::structure)->at(0)] =
          StructureField()->DofRowMap();

      // initialize map extractor associated with blocks of global system matrix
      maps_systemmatrix_ =
          Teuchos::rcp(new LINALG::MultiMapExtractor(*DofRowMap(), maps_systemmatrix));

      // initialize map extractor associated with all degrees of freedom inside structural field
      map_structure_ = Teuchos::rcp(
          new LINALG::MultiMapExtractor(*StructureField()->Discretization()->DofRowMap(),
              std::vector<Teuchos::RCP<const Epetra_Map>>(1, StructureField()->DofRowMap())));

      // safety check
      map_structure_->CheckForValidMapExtractor();

      break;
    }

    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // safety check
  maps_systemmatrix_->CheckForValidMapExtractor();

  // perform initializations associated with global system matrix
  switch (matrixtype_)
  {
    case LINALG::MatrixType::block_field:
    {
      // safety check
      if (!solver_->Params().isSublist("AMGnxn Parameters"))
        dserror("Global system matrix with block structure requires AMGnxn block preconditioner!");

      // feed AMGnxn block preconditioner with null space information for each block of global block
      // system matrix
      BuildNullSpaces();

      break;
    }

    case LINALG::MatrixType::sparse:
    {
      // safety check
      if (ScaTraField()->SystemMatrix() == Teuchos::null)
        dserror("Incompatible matrix type associated with scalar transport field!");
      break;
    }

    default:
    {
      dserror("Type of global system matrix for scalar-structure interaction not recognized!");
      break;
    }
  }

  // initialize subblocks and system matrix
  ssi_matrices_ = Teuchos::rcp(new SSI::UTILS::SSIMatrices(*this, interface_map_scatra));

  // initialize strategy for assembly
  strategy_assemble_ = SSI::BuildAssembleStrategy(
      Teuchos::rcp(this, false), matrixtype_, ScaTraField()->MatrixType());

  // initialize object, that performs evaluations of OD coupling
  if (IsScaTraManifold())
  {
    scatrastructureOffDiagcoupling_ = Teuchos::rcp(new SSI::ScatraManifoldStructureOffDiagCoupling(
        MapStructure(), MapsSubProblems()->Map(GetProblemPosition(Subproblem::scalar_transport)),
        MapsSubProblems()->Map(GetProblemPosition(Subproblem::structure)),
        MapsSubProblems()->Map(GetProblemPosition(Subproblem::manifold)),
        MapStructureManifold()->Map(0), InterfaceCouplingAdapterStructure(),
        InterfaceCouplingAdapterStructure3DomainIntersection(), interface_map_scatra,
        meshtying_strategy_s2i_, ScaTraBaseAlgorithm(), ScaTraManifoldBaseAlgorithm(),
        StructureField(), Meshtying3DomainIntersection()));
  }
  else
  {
    scatrastructureOffDiagcoupling_ = Teuchos::rcp(new SSI::ScatraStructureOffDiagCoupling(
        MapStructure(), MapsSubProblems()->Map(GetProblemPosition(Subproblem::scalar_transport)),
        MapsSubProblems()->Map(GetProblemPosition(Subproblem::structure)),
        InterfaceCouplingAdapterStructure(), InterfaceCouplingAdapterStructure3DomainIntersection(),
        interface_map_scatra, meshtying_strategy_s2i_, ScaTraBaseAlgorithm(), StructureField(),
        Meshtying3DomainIntersection()));
  }
  // instantiate appropriate equilibration class
  strategy_equilibration_ = LINALG::BuildEquilibration(
      matrixtype_, GetBlockEquilibration(), MapsSubProblems()->FullMap());
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void SSI::SSIMono::SetupModelEvaluator() const
{
  // construct and register structural model evaluator if necessary

  const bool do_output_stress =
      DRT::INPUT::IntegralValue<INPAR::STR::StressType>(
          DRT::Problem::Instance()->IOParams(), "STRUCT_STRESS") != INPAR::STR::stress_none;
  const bool smooth_output_interface_stress = DRT::INPUT::IntegralValue<bool>(
      DRT::Problem::Instance()->SSIControlParams().sublist("MONOLITHIC"),
      "SMOOTH_OUTPUT_INTERFACE_STRESS");

  if (Meshtying3DomainIntersection() and smooth_output_interface_stress)
    dserror("Smoothing of interface stresses not implemented for triple meshtying.");

  if (smooth_output_interface_stress and !do_output_stress)
    dserror("Smoothing of interface stresses only when stress output is written.");

  if (do_output_stress and SSIInterfaceMeshtying())
  {
    StructureBaseAlgorithm()->RegisterModelEvaluator("Monolithic Coupling Model",
        Teuchos::rcp(new STR::MODELEVALUATOR::MonolithicSSI(
            Teuchos::rcp(this, false), smooth_output_interface_stress)));
  }
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void SSI::SSIMono::SolveLinearSystem()
{
  strategy_equilibration_->EquilibrateSystem(
      ssi_matrices_->SystemMatrix(), residual_, *MapsSystemMatrix());

  // solve global system of equations
  // Dirichlet boundary conditions have already been applied to global system of equations
  solver_->Solve(ssi_matrices_->SystemMatrix()->EpetraOperator(), increment_, residual_, true,
      IterationCount() == 1);

  strategy_equilibration_->UnequilibrateIncrement(increment_);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SSIMono::NewtonLoop()
{
  // reset counter for Newton-Raphson iteration
  ResetIterationCount();

  // start Newton-Raphson iteration
  while (true)
  {
    // update iteration counter
    IncrementIterationCount();

    // reset timer
    timer_->ResetStartTime();

    // store time before evaluating elements and assembling global system of equations
    double time = timer_->WallTime();

    // evaluate subproblems and get all matrices and right-hand-sides
    EvaluateSubproblems();

    // assemble global system of equations
    AssembleMatAndRHS();

    // determine time needed for evaluating elements and assembling global system of
    // equations, and take maximum over all processors via communication
    double mydtele = timer_->WallTime() - time;
    Comm().MaxAll(&mydtele, &dtele_, 1);

    // safety check
    if (!ssi_matrices_->SystemMatrix()->Filled())
      dserror("Complete() has not been called on global system matrix yet!");

    // check termination criterion for Newton-Raphson iteration
    if (strategy_convcheck_->ExitNewtonRaphson(*this)) break;

    // initialize global increment vector
    increment_->PutScalar(0.0);

    // store time before solving global system of equations
    time = timer_->WallTime();

    SolveLinearSystem();

    // determine time needed for solving global system of equations,
    // and take maximum over all processors via communication
    double mydtsolve = timer_->WallTime() - time;
    Comm().MaxAll(&mydtsolve, &dtsolve_, 1);

    // output performance statistics associated with linear solver into text file if
    // applicable
    if (DRT::INPUT::IntegralValue<bool>(
            *ScaTraField()->ScatraParameterList(), "OUTPUTLINSOLVERSTATS"))
      ScaTraField()->OutputLinSolverStats(
          *solver_, dtsolve_, Step(), IterationCount(), residual_->Map().NumGlobalElements());

    // update states for next Newton iteration
    UpdateIterScaTra();
    UpdateIterStructure();

  }  // Newton-Raphson iteration
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SSIMono::Timeloop()
{
  // output initial scalar transport solution to screen and files
  if (Step() == 0)
  {
    SetStructSolution(StructureField()->Dispnp(), StructureField()->Velnp());
    ScaTraField()->Output();
    if (IsScaTraManifold()) ScaTraManifold()->Output();
  }

  // time loop
  while (NotFinished() and ScaTraField()->NotFinished())
  {
    // prepare time step
    PrepareTimeStep();

    // store time before calling nonlinear solver
    const double time = timer_->WallTime();

    // evaluate time step
    NewtonLoop();

    // determine time spent by nonlinear solver and take maximum over all processors via
    // communication
    double mydtnonlinsolve(timer_->WallTime() - time), dtnonlinsolve(0.);
    Comm().MaxAll(&mydtnonlinsolve, &dtnonlinsolve, 1);

    // output performance statistics associated with nonlinear solver into *.csv file if
    // applicable
    if (DRT::INPUT::IntegralValue<int>(
            *ScaTraField()->ScatraParameterList(), "OUTPUTNONLINSOLVERSTATS"))
      ScaTraField()->OutputNonlinSolverStats(IterationCount(), dtnonlinsolve, Step(), Comm());

    // prepare structure output
    StructureField()->PrepareOutput();

    // update scalar transport and structure fields
    Update();

    // output solution to screen and files
    Output();
  }
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void SSI::SSIMono::Update()
{
  // update scalar transport field
  ScaTraField()->Update();
  if (IsScaTraManifold()) ScaTraManifold()->Update();

  // update structure field
  StructureField()->Update();
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void SSI::SSIMono::UpdateIterScaTra()
{
  // update scalar transport field
  ScaTraField()->UpdateIter(MapsSubProblems()->ExtractVector(
      increment_, GetProblemPosition(Subproblem::scalar_transport)));
  ScaTraField()->ComputeIntermediateValues();

  if (IsScaTraManifold())
  {
    ScaTraManifold()->UpdateIter(
        MapsSubProblems()->ExtractVector(increment_, GetProblemPosition(Subproblem::manifold)));
    ScaTraManifold()->ComputeIntermediateValues();
  }
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void SSI::SSIMono::UpdateIterStructure()
{
  // set up structural increment vector
  const Teuchos::RCP<Epetra_Vector> increment_structure =
      MapsSubProblems()->ExtractVector(increment_, GetProblemPosition(Subproblem::structure));

  // consider structural meshtying. Copy master increments and displacements to slave side.
  if (SSIInterfaceMeshtying())
  {
    MapsCoupStruct()->InsertVector(
        InterfaceCouplingAdapterStructure()->MasterToSlave(
            MapsCoupStruct()->ExtractVector(StructureField()->Dispnp(), 2)),
        1, StructureField()->WriteAccessDispnp());
    StructureField()->SetState(StructureField()->WriteAccessDispnp());
    MapsCoupStruct()->InsertVector(InterfaceCouplingAdapterStructure()->MasterToSlave(
                                       MapsCoupStruct()->ExtractVector(increment_structure, 2)),
        1, increment_structure);

    if (Meshtying3DomainIntersection())
    {
      MapsCoupStruct3DomainIntersection()->InsertVector(
          InterfaceCouplingAdapterStructure3DomainIntersection()->MasterToSlave(
              MapsCoupStruct3DomainIntersection()->ExtractVector(StructureField()->Dispnp(), 2)),
          1, StructureField()->WriteAccessDispnp());
      StructureField()->SetState(StructureField()->WriteAccessDispnp());
      MapsCoupStruct3DomainIntersection()->InsertVector(
          InterfaceCouplingAdapterStructure3DomainIntersection()->MasterToSlave(
              MapsCoupStruct3DomainIntersection()->ExtractVector(increment_structure, 2)),
          1, increment_structure);
    }
  }

  // update displacement of structure field
  StructureField()->UpdateStateIncrementally(increment_structure);
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
Teuchos::RCP<std::vector<int>> SSI::SSIMono::GetBlockPositions(Subproblem subproblem) const
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
        for (int i = 0; i < ScaTraField()->BlockMaps().NumMaps(); ++i)
          block_position->emplace_back(i);
      }
      break;
    }
    case Subproblem::manifold:
    {
      if (ScaTraManifold()->MatrixType() == LINALG::MatrixType::sparse)
        block_position->emplace_back(2);
      else
      {
        for (int i = 0; i < static_cast<int>(ScaTraManifold()->BlockMaps().NumMaps()); ++i)
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
int SSI::SSIMono::GetProblemPosition(Subproblem subproblem) const
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
    case Subproblem::manifold:
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

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
Teuchos::RCP<std::vector<LINALG::EquilibrationMethod>> SSI::SSIMono::GetBlockEquilibration()
{
  Teuchos::RCP<std::vector<LINALG::EquilibrationMethod>> equilibration_method_vector;
  switch (matrixtype_)
  {
    case LINALG::MatrixType::sparse:
    {
      equilibration_method_vector = Teuchos::rcp(
          new std::vector<LINALG::EquilibrationMethod>(1, equilibration_method_.global));
      break;
    }
    case LINALG::MatrixType::block_field:
    {
      if (equilibration_method_.global != LINALG::EquilibrationMethod::local)
      {
        equilibration_method_vector = Teuchos::rcp(
            new std::vector<LINALG::EquilibrationMethod>(1, equilibration_method_.global));
      }
      else if (equilibration_method_.structure == LINALG::EquilibrationMethod::none and
               equilibration_method_.scatra == LINALG::EquilibrationMethod::none)
      {
        equilibration_method_vector = Teuchos::rcp(
            new std::vector<LINALG::EquilibrationMethod>(1, LINALG::EquilibrationMethod::none));
      }
      else
      {
        auto block_positions_scatra = GetBlockPositions(Subproblem::scalar_transport);
        auto block_position_structure = GetBlockPositions(Subproblem::structure);
        auto block_positions_scatra_manifold =
            IsScaTraManifold() ? GetBlockPositions(Subproblem::manifold) : Teuchos::null;

        equilibration_method_vector = Teuchos::rcp(new std::vector<LINALG::EquilibrationMethod>(
            block_positions_scatra->size() + block_position_structure->size() +
                (IsScaTraManifold() ? block_positions_scatra_manifold->size() : 0),
            LINALG::EquilibrationMethod::none));

        for (const int block_position_scatra : *block_positions_scatra)
          equilibration_method_vector->at(block_position_scatra) = equilibration_method_.scatra;

        equilibration_method_vector->at(block_position_structure->at(0)) =
            equilibration_method_.structure;

        if (IsScaTraManifold())
        {
          for (const int block_position_scatra_manifold : *block_positions_scatra_manifold)
          {
            equilibration_method_vector->at(block_position_scatra_manifold) =
                equilibration_method_.scatra;
          }
        }
      }

      break;
    }
    default:
    {
      dserror("Invalid matrix type associated with system matrix field!");
      break;
    }
  }
  return equilibration_method_vector;
}