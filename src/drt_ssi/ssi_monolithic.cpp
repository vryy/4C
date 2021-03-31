/*--------------------------------------------------------------------------*/
/*! \file
\brief monolithic scalar-structure interaction

\level 2

*/
/*--------------------------------------------------------------------------*/
#include "ssi_monolithic.H"

#include "ssi_coupling.H"
#include "ssi_manifold_flux_evaluator.H"
#include "ssi_monolithic_assemble_strategy.H"
#include "ssi_monolithic_contact_strategy.H"
#include "ssi_monolithic_convcheck_strategies.H"
#include "ssi_monolithic_dbc_handler.H"
#include "ssi_monolithic_evaluate_OffDiag.H"
#include "ssi_monolithic_meshtying_strategy.H"
#include "ssi_resulttest.H"
#include "ssi_str_model_evaluator_monolithic.H"
#include "ssi_utils.H"

#include "../drt_adapter/ad_str_ssiwrapper.H"
#include "../drt_adapter/ad_str_structure_new.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"

#include "../drt_contact/contact_nitsche_strategy_ssi.H"

#include "../drt_inpar/inpar_scatra.H"
#include "../drt_inpar/inpar_ssi.H"

#include "../drt_io/io_control.H"

#include "../drt_lib/drt_assemblestrategy.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_locsys.H"

#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_scatra/scatra_timint_meshtying_strategy_s2i.H"

#include "../drt_structure_new/str_model_evaluator_contact.H"

#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_matrixtransform.H"
#include "../linalg/linalg_equilibrate.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils_sparse_algebra_assemble.H"
#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"

#include <Epetra_Time.h>


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
SSI::SSIMono::SSIMono(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams)
    : SSIBase(comm, globaltimeparams),
      contact_strategy_nitsche_(Teuchos::null),
      dbc_handler_(Teuchos::null),
      dtele_(0.0),
      dtsolve_(0.0),
      equilibration_method_{Teuchos::getIntegralValue<LINALG::EquilibrationMethod>(
                                globaltimeparams.sublist("MONOLITHIC"), "EQUILIBRATION"),
          Teuchos::getIntegralValue<LINALG::EquilibrationMethod>(
              globaltimeparams.sublist("MONOLITHIC"), "EQUILIBRATION_SCATRA"),
          Teuchos::getIntegralValue<LINALG::EquilibrationMethod>(
              globaltimeparams.sublist("MONOLITHIC"), "EQUILIBRATION_STRUCTURE")},
      manifoldscatraflux_(Teuchos::null),
      map_structure_(Teuchos::null),
      maps_scatra_(Teuchos::null),
      maps_sub_problems_(Teuchos::null),
      maps_systemmatrix_(Teuchos::null),
      matrixtype_(Teuchos::getIntegralValue<LINALG::MatrixType>(
          globaltimeparams.sublist("MONOLITHIC"), "MATRIXTYPE")),
      scatrastructureOffDiagcoupling_(Teuchos::null),
      solver_(Teuchos::rcp(
          new LINALG::Solver(DRT::Problem::Instance()->SolverParams(
                                 globaltimeparams.sublist("MONOLITHIC").get<int>("LINEAR_SOLVER")),
              comm, DRT::Problem::Instance()->ErrorFile()->Handle()))),
      ssi_matrices_(Teuchos::null),
      ssi_vectors_(Teuchos::null),
      strategy_assemble_(Teuchos::null),
      strategy_contact_(Teuchos::null),
      strategy_convcheck_(Teuchos::null),
      strategy_equilibration_(Teuchos::null),
      strategy_meshtying_(Teuchos::null),
      timer_(Teuchos::rcp(new Epetra_Time(comm)))
{
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void SSI::SSIMono::ApplyContactToSubProblems()
{
  strategy_contact_->ApplyContactToScatraResidual(ssi_vectors_->ScatraResidual());

  strategy_contact_->ApplyContactToScatraScatra(ssi_matrices_->ScaTraMatrix());

  strategy_contact_->ApplyContactToScatraStructure(ssi_matrices_->ScaTraStructureMatrix());

  strategy_contact_->ApplyContactToStructureScatra(ssi_matrices_->StructureScaTraMatrix());
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void SSI::SSIMono::ApplyDBCToSystem()
{
  // apply Dirichlet boundary conditions to global system matrix
  dbc_handler_->ApplyDBCToSystemMatrix(ssi_matrices_->SystemMatrix());

  // apply Dirichlet boundary conditions to global RHS
  dbc_handler_->ApplyDBCToRHS(ssi_vectors_->Residual());
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void SSI::SSIMono::ApplyMeshtyingToSubProblems()
{
  if (SSIInterfaceMeshtying())
  {
    if (IsScaTraManifold())
    {
      strategy_meshtying_->ApplyMeshtyingToScatraManifoldStructure(
          ssi_matrices_->ScaTraManifoldStructureMatrix());

      strategy_meshtying_->ApplyMeshtyingToScatraManifoldStructure(
          manifoldscatraflux_->MatrixManifoldStructure());

      strategy_meshtying_->ApplyMeshtyingToScatraStructure(
          manifoldscatraflux_->MatrixScaTraStructure());
    }

    strategy_meshtying_->ApplyMeshtyingToScatraStructure(ssi_matrices_->ScaTraStructureMatrix());

    strategy_meshtying_->ApplyMeshtyingToStructureMatrix(
        *ssi_matrices_->StructureMatrix(), StructureField()->SystemMatrix());

    strategy_meshtying_->ApplyMeshtyingToStructureScatra(ssi_matrices_->StructureScaTraMatrix());

    ssi_vectors_->StructureResidual()->Update(
        1.0, strategy_meshtying_->ApplyMeshtyingToStructureRHS(StructureField()->RHS()), 1.0);
  }
  // copy the structure residual and matrix if we do not have a mesh tying problem
  else
  {
    ssi_vectors_->StructureResidual()->Update(1.0, *(StructureField()->RHS()), 1.0);
    ssi_matrices_->StructureMatrix()->Add(*StructureField()->SystemMatrix(), false, 1.0, 1.0);
  }
  // call fill complete on ssi structure matrix
  ssi_matrices_->StructureMatrix()->Complete();
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SSIMono::AssembleMatAndRHS()
{
  AssembleMatScaTra();

  AssembleMatStructure();

  if (IsScaTraManifold()) AssembleMatScaTraManifold();

  // finalize global system matrix
  ssi_matrices_->SystemMatrix()->Complete();

  // assemble monolithic RHS
  strategy_assemble_->AssembleRHS(ssi_vectors_->Residual(), ssi_vectors_->ScatraResidual(),
      ssi_vectors_->StructureResidual(),
      IsScaTraManifold() ? ScaTraManifold()->Residual() : Teuchos::null,
      IsScaTraManifold() ? manifoldscatraflux_->RHSManifold() : Teuchos::null,
      IsScaTraManifold() ? manifoldscatraflux_->RHSScaTra() : Teuchos::null);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SSIMono::AssembleMatScaTra()
{
  // assemble scatra-scatra block into system matrix
  strategy_assemble_->AssembleScatraScatra(
      ssi_matrices_->SystemMatrix(), ssi_matrices_->ScaTraMatrix());

  // assemble scatra-structure block into system matrix
  strategy_assemble_->AssembleScatraStructure(
      ssi_matrices_->SystemMatrix(), ssi_matrices_->ScaTraStructureMatrix());
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SSIMono::AssembleMatScaTraManifold()
{
  // assemble scatra manifold - scatra manifold block into system matrix
  strategy_assemble_->AssembleScatramanifoldScatramanifold(
      ssi_matrices_->SystemMatrix(), ScaTraManifold()->SystemMatrixOperator());

  // assemble scatra manifold-structure block into system matrix
  strategy_assemble_->AssembleScatramanifoldStructure(
      ssi_matrices_->SystemMatrix(), ssi_matrices_->ScaTraManifoldStructureMatrix());

  // assemble contributions from scatra - scatra manifold coupling: derivs. of manifold side w.r.t.
  // manifold side
  strategy_assemble_->AssembleScatramanifoldScatramanifold(
      ssi_matrices_->SystemMatrix(), manifoldscatraflux_->SystemMatrixManifold());

  // assemble contributions from scatra - scatra manifold coupling: derivs. of scatra side w.r.t.
  // scatra side
  strategy_assemble_->AssembleScatraScatra(
      ssi_matrices_->SystemMatrix(), manifoldscatraflux_->SystemMatrixScaTra());

  // assemble contributions from scatra - scatra manifold coupling: derivs. of manifold side w.r.t.
  // scatra side
  strategy_assemble_->AssembleScatraScatramanifold(
      ssi_matrices_->SystemMatrix(), manifoldscatraflux_->MatrixScaTraManifold());

  // assemble contributions from scatra - scatra manifold coupling: derivs. of scatra side w.r.t.
  // manifold side
  strategy_assemble_->AssembleScatramanifoldScatra(
      ssi_matrices_->SystemMatrix(), manifoldscatraflux_->MatrixManifoldScatra());

  strategy_assemble_->AssembleScatramanifoldStructure(
      ssi_matrices_->SystemMatrix(), manifoldscatraflux_->MatrixManifoldStructure());

  strategy_assemble_->AssembleScatraStructure(
      ssi_matrices_->SystemMatrix(), manifoldscatraflux_->MatrixScaTraStructure());
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SSIMono::AssembleMatStructure()
{  // assemble structure-scatra block into system matrix
  strategy_assemble_->AssembleStructureScatra(
      ssi_matrices_->SystemMatrix(), ssi_matrices_->StructureScaTraMatrix());

  // assemble structure-structure block into system matrix
  strategy_assemble_->AssembleStructureStructure(
      ssi_matrices_->SystemMatrix(), ssi_matrices_->StructureMatrix());
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SSIMono::EvaluateSubproblems()
{
  // clear all matrices and residuals from previous Newton iteration
  ssi_matrices_->ClearMatrices();
  ssi_vectors_->ClearResiduals();

  // needed to communicate to NOX state
  StructureField()->SetState(StructureField()->WriteAccessDispnp());

  // distribute states to other fields
  SetStructSolution(StructureField()->Dispnp(), StructureField()->Velnp());
  SetScatraSolution(ScaTraField()->Phinp());
  if (IsScaTraManifold()) SetScatraManifoldSolution(ScaTraManifold()->Phinp());

  // evaluate temperature from function and set to structural discretization
  EvaluateAndSetTemperatureField();

  // build system matrix and residual for structure field
  StructureField()->Evaluate();

  // build system matrix and residual for scalar transport field
  EvaluateScaTra();

  // build system matrix and residual for scalar transport field on manifold
  if (IsScaTraManifold()) EvaluateScaTraManifold();

  // build all off diagonal matrices
  EvaluateOffDiagContributions();

  // apply mesh tying to sub problems
  ApplyMeshtyingToSubProblems();

  // apply contact contributions to sub problems
  if (SSIInterfaceContact()) ApplyContactToSubProblems();
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void SSI::SSIMono::EvaluateOffDiagContributions()
{
  // evaluate off-diagonal scatra-structure block (domain contributions) of global system matrix
  scatrastructureOffDiagcoupling_->EvaluateOffDiagBlockScatraStructureDomain(
      ssi_matrices_->ScaTraStructureMatrix());

  // evaluate off-diagonal scatra-structure block (interface contributions) of global system matrix
  if (SSIInterfaceMeshtying())
    scatrastructureOffDiagcoupling_->EvaluateOffDiagBlockScatraStructureInterface(
        ssi_matrices_->ScaTraStructureMatrix());

  // evaluate off-diagonal structure-scatra block (we only have domain contributions so far) of
  // global system matrix
  scatrastructureOffDiagcoupling_->EvaluateOffDiagBlockStructureScatraDomain(
      ssi_matrices_->StructureScaTraMatrix());

  if (IsScaTraManifold())
  {
    // evaluate off-diagonal manifold-structure block of global system matrix
    scatrastructureOffDiagcoupling_->EvaluateOffDiagBlockScatraManifoldStructureDomain(
        ssi_matrices_->ScaTraManifoldStructureMatrix());
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

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::SSIMono::SetupContactStrategy()
{
  // get the contact solution strategy
  auto contact_solution_type = DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(
      DRT::Problem::Instance()->ContactDynamicParams(), "STRATEGY");

  if (contact_solution_type == INPAR::CONTACT::solution_nitsche)
  {
    if (DRT::INPUT::IntegralValue<INPAR::STR::IntegrationStrategy>(
            DRT::Problem::Instance()->StructuralDynamicParams(), "INT_STRATEGY") !=
        INPAR::STR::int_standard)
      dserror("ssi contact only with new structural time integration");

    // get the contact model evaluator and store a pointer to the strategy
    auto& model_evaluator_contact = dynamic_cast<STR::MODELEVALUATOR::Contact&>(
        StructureField()->ModelEvaluator(INPAR::STR::model_contact));
    contact_strategy_nitsche_ = Teuchos::rcp_dynamic_cast<CONTACT::CoNitscheStrategySsi>(
        model_evaluator_contact.StrategyPtr(), true);
  }
  else
    dserror("Only Nitsche contact implemented for SSI problems at the moment!");
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
  if (IsScaTraManifold())
  {
    // domain output
    ScaTraManifold()->Output();
    // coupling output
    if (manifoldscatraflux_->DoOutput()) manifoldscatraflux_->Output();
  }

  // output structure field
  StructureField()->Output();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::SSIMono::ReadRestart(int restart)
{
  // call base class
  SSIBase::ReadRestart(restart);

  // do ssi contact specific tasks
  if (SSIInterfaceContact())
  {
    SetupContactStrategy();
    SetSSIContactStates(ScaTraField()->Phinp());
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::SSIMono::ReadRestartfromTime(double restarttime)
{
  // call base class
  SSIBase::ReadRestartfromTime(restarttime);

  // do ssi contact specific tasks
  if (SSIInterfaceContact())
  {
    SetupContactStrategy();
    SetSSIContactStates(ScaTraField()->Phinp());
  }
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
  if (IsScaTraManifold()) SetScatraManifoldSolution(ScaTraManifold()->Phinp());

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

  if (SSIInterfaceMeshtying() and
      MeshtyingStrategyS2I()->CouplingType() != INPAR::S2I::coupling_matching_nodes)
  {
    dserror(
        "Monolithic scalar-structure interaction only implemented for scatra-scatra "
        "interface coupling with matching interface nodes!");
  }

  if (SSIInterfaceContact() and !IsRestart()) SetupContactStrategy();
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
        {MeshtyingStrategyS2I()->CouplingAdapter()->MasterDofMap(),
            MeshtyingStrategyS2I()->CouplingAdapter()->SlaveDofMap()});
  }

  // initialize global map extractor
  std::vector<Teuchos::RCP<const Epetra_Map>> partial_maps(
      IsScaTraManifold() ? 3 : 2, Teuchos::null);
  Teuchos::RCP<const Epetra_Map> merged_map;

  partial_maps[GetProblemPosition(Subproblem::scalar_transport)] =
      Teuchos::rcp(new Epetra_Map(*ScaTraField()->DofRowMap()));
  partial_maps[GetProblemPosition(Subproblem::structure)] =
      Teuchos::rcp(new Epetra_Map(*StructureField()->DofRowMap()));
  if (IsScaTraManifold())
  {
    partial_maps[GetProblemPosition(Subproblem::manifold)] =
        Teuchos::rcp(new Epetra_Map(*ScaTraManifold()->DofRowMap()));
    auto temp_map = LINALG::MergeMap(partial_maps[0], partial_maps[1], false);
    merged_map = LINALG::MergeMap(temp_map, partial_maps[2], false);
  }
  else
    merged_map = LINALG::MergeMap(partial_maps[0], partial_maps[1], false);

  maps_sub_problems_ = Teuchos::rcp(new LINALG::MultiMapExtractor(*merged_map, partial_maps));
  // check global map extractor
  maps_sub_problems_->CheckForValidMapExtractor();

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

  // initialize sub blocks and system matrix
  ssi_matrices_ = Teuchos::rcp(new SSI::UTILS::SSIMatrices(*this));

  // initialize residual and increment vectors
  ssi_vectors_ = Teuchos::rcp(new SSI::UTILS::SSIVectors(*this));

  // initialize strategy for assembly
  strategy_assemble_ = SSI::BuildAssembleStrategy(*this, matrixtype_, ScaTraField()->MatrixType());

  if (IsScaTraManifold())
  {
    // initialize object, that performs evaluations of OD coupling
    scatrastructureOffDiagcoupling_ = Teuchos::rcp(new SSI::ScatraManifoldStructureOffDiagCoupling(
        MapStructure(), MapsSubProblems()->Map(GetProblemPosition(Subproblem::scalar_transport)),
        MapsSubProblems()->Map(GetProblemPosition(Subproblem::structure)),
        MapsSubProblems()->Map(GetProblemPosition(Subproblem::manifold)),
        MapStructureOnScaTraManifold()->Map(0), InterfaceCouplingAdapterStructure(),
        InterfaceCouplingAdapterStructure3DomainIntersection(), interface_map_scatra,
        MeshtyingStrategyS2I(), ScaTraBaseAlgorithm(), ScaTraManifoldBaseAlgorithm(),
        StructureField(), Meshtying3DomainIntersection()));

    // initialize object, that performs evaluations of scatra - scatra on manifold coupling
    manifoldscatraflux_ = Teuchos::rcp(new SSI::ScaTraManifoldScaTraFluxEvaluator(*this));
  }
  else
  {
    scatrastructureOffDiagcoupling_ = Teuchos::rcp(new SSI::ScatraStructureOffDiagCoupling(
        MapStructure(), MapsSubProblems()->Map(GetProblemPosition(Subproblem::scalar_transport)),
        MapsSubProblems()->Map(GetProblemPosition(Subproblem::structure)),
        InterfaceCouplingAdapterStructure(), InterfaceCouplingAdapterStructure3DomainIntersection(),
        interface_map_scatra, MeshtyingStrategyS2I(), ScaTraBaseAlgorithm(), StructureField(),
        Meshtying3DomainIntersection()));
  }
  // instantiate appropriate equilibration class
  strategy_equilibration_ = LINALG::BuildEquilibration(
      matrixtype_, GetBlockEquilibration(), MapsSubProblems()->FullMap());

  // instantiate appropriate contact class
  strategy_contact_ = SSI::BuildContactStrategy(*this, ScaTraField()->MatrixType());

  // instantiate appropriate mesh tying class
  strategy_meshtying_ =
      SSI::BuildMeshtyingStrategy(*this, matrixtype_, ScaTraField()->MatrixType());

  // instantiate Dirichlet boundary condition handler class
  dbc_handler_ = SSI::BuildDBCHandler(Teuchos::rcp(this, false), matrixtype_);
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
void SSI::SSIMono::SetScatraSolution(Teuchos::RCP<const Epetra_Vector> phi) const
{
  // call base class
  SSIBase::SetScatraSolution(phi);

  // set state for contact evaluation
  SetSSIContactStates(phi);
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void SSI::SSIMono::SetSSIContactStates(Teuchos::RCP<const Epetra_Vector> phi) const
{
  if (contact_strategy_nitsche_ != Teuchos::null)
    contact_strategy_nitsche_->SetState(MORTAR::state_scalar, *phi);
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void SSI::SSIMono::SolveLinearSystem()
{
  strategy_equilibration_->EquilibrateSystem(
      ssi_matrices_->SystemMatrix(), ssi_vectors_->Residual(), *MapsSystemMatrix());

  // solve global system of equations
  // Dirichlet boundary conditions have already been applied to global system of equations
  solver_->Solve(ssi_matrices_->SystemMatrix()->EpetraOperator(), ssi_vectors_->Increment(),
      ssi_vectors_->Residual(), true, IterationCount() == 1);

  strategy_equilibration_->UnequilibrateIncrement(ssi_vectors_->Increment());
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

    // evaluate sub problems and get all matrices and right-hand-sides
    EvaluateSubproblems();

    // assemble global system of equations
    AssembleMatAndRHS();

    // apply the Dirichlet boundary conditions to global system
    ApplyDBCToSystem();

    // determine time needed for evaluating elements and assembling global system of
    // equations, and take maximum over all processors via communication
    double mydtele = timer_->WallTime() - time;
    Comm().MaxAll(&mydtele, &dtele_, 1);

    // safety check
    if (!ssi_matrices_->SystemMatrix()->Filled())
      dserror("Complete() has not been called on global system matrix yet!");

    // check termination criterion for Newton-Raphson iteration
    if (strategy_convcheck_->ExitNewtonRaphson(*this)) break;

    // clear the global increment vector
    ssi_vectors_->ClearIncrement();

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
      ScaTraField()->OutputLinSolverStats(*solver_, dtsolve_, Step(), IterationCount(),
          ssi_vectors_->Residual()->Map().NumGlobalElements());

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

    PrepareOutput();

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
      ssi_vectors_->Increment(), GetProblemPosition(Subproblem::scalar_transport)));
  ScaTraField()->ComputeIntermediateValues();

  if (IsScaTraManifold())
  {
    ScaTraManifold()->UpdateIter(MapsSubProblems()->ExtractVector(
        ssi_vectors_->Increment(), GetProblemPosition(Subproblem::manifold)));
    ScaTraManifold()->ComputeIntermediateValues();
  }
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void SSI::SSIMono::UpdateIterStructure()
{
  // set up structural increment vector
  const Teuchos::RCP<Epetra_Vector> increment_structure = MapsSubProblems()->ExtractVector(
      ssi_vectors_->Increment(), GetProblemPosition(Subproblem::structure));

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

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void SSI::SSIMono::EvaluateScaTra()
{
  // evaluate the scatra field
  ScaTraField()->PrepareLinearSolve();

  // copy the matrix to the corresponding ssi matrix and complete it such that additional
  // contributions like contact contributions can be added before assembly
  ssi_matrices_->ScaTraMatrix()->Add(*ScaTraField()->SystemMatrixOperator(), false, 1.0, 1.0);
  ssi_matrices_->ScaTraMatrix()->Complete();

  // copy the residual to the corresponding ssi vector to enable application of contact
  // contributions before assembly
  ssi_vectors_->ScatraResidual()->Update(1.0, *ScaTraField()->Residual(), 1.0);
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void SSI::SSIMono::EvaluateScaTraManifold()
{
  // evaluate single problem
  ScaTraManifold()->PrepareLinearSolve();

  // evaluate coupling fluxes
  manifoldscatraflux_->Evaluate();
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void SSI::SSIMono::PrepareOutput()
{
  StructureField()->PrepareOutput();

  // prepare output of coupling sctra manifold - scatra
  if (IsScaTraManifold() and manifoldscatraflux_->DoOutput())
    manifoldscatraflux_->EvaluateScaTraManifoldInflow();
}