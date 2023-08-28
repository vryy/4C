/*--------------------------------------------------------------------------*/
/*! \file
\brief monolithic scalar-structure interaction

\level 2

*/
/*--------------------------------------------------------------------------*/
#include "baci_ssi_monolithic.H"

#include "baci_adapter_scatra_base_algorithm.H"
#include "baci_adapter_str_ssiwrapper.H"
#include "baci_adapter_str_structure_new.H"
#include "baci_contact_nitsche_strategy_ssi.H"
#include "baci_inpar_ssi.H"
#include "baci_io_control.H"
#include "baci_lib_assemblestrategy.H"
#include "baci_lib_globalproblem.H"
#include "baci_lib_locsys.H"
#include "baci_linalg_equilibrate.H"
#include "baci_linalg_mapextractor.H"
#include "baci_linalg_matrixtransform.H"
#include "baci_linalg_utils_sparse_algebra_assemble.H"
#include "baci_linalg_utils_sparse_algebra_create.H"
#include "baci_linalg_utils_sparse_algebra_manipulation.H"
#include "baci_linalg_utils_sparse_algebra_print.H"
#include "baci_linear_solver_method_linalg.H"
#include "baci_scatra_timint_elch.H"
#include "baci_scatra_timint_meshtying_strategy_s2i.H"
#include "baci_ssi_coupling.H"
#include "baci_ssi_manifold_utils.H"
#include "baci_ssi_monolithic_assemble_strategy.H"
#include "baci_ssi_monolithic_contact_strategy.H"
#include "baci_ssi_monolithic_convcheck_strategies.H"
#include "baci_ssi_monolithic_dbc_handler.H"
#include "baci_ssi_monolithic_evaluate_OffDiag.H"
#include "baci_ssi_monolithic_meshtying_strategy.H"
#include "baci_ssi_utils.H"
#include "baci_structure_new_model_evaluator_contact.H"


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
SSI::SSIMono::SSIMono(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams)
    : SSIBase(comm, globaltimeparams),
      contact_strategy_nitsche_(Teuchos::null),
      dbc_handler_(Teuchos::null),
      dt_eval_(0.0),
      dt_solve_(0.0),
      equilibration_method_{Teuchos::getIntegralValue<CORE::LINALG::EquilibrationMethod>(
                                globaltimeparams.sublist("MONOLITHIC"), "EQUILIBRATION"),
          Teuchos::getIntegralValue<CORE::LINALG::EquilibrationMethod>(
              globaltimeparams.sublist("MONOLITHIC"), "EQUILIBRATION_SCATRA"),
          Teuchos::getIntegralValue<CORE::LINALG::EquilibrationMethod>(
              globaltimeparams.sublist("MONOLITHIC"), "EQUILIBRATION_STRUCTURE")},
      manifoldscatraflux_(Teuchos::null),
      matrixtype_(Teuchos::getIntegralValue<CORE::LINALG::MatrixType>(
          globaltimeparams.sublist("MONOLITHIC"), "MATRIXTYPE")),
      print_matlab_(DRT::INPUT::IntegralValue<bool>(
          globaltimeparams.sublist("MONOLITHIC"), "PRINT_MAT_RHS_MAP_MATLAB")),
      scatrastructureOffDiagcoupling_(Teuchos::null),
      solver_(Teuchos::rcp(new CORE::LINALG::Solver(
          DRT::Problem::Instance()->SolverParams(
              globaltimeparams.sublist("MONOLITHIC").get<int>("LINEAR_SOLVER")),
          comm, DRT::Problem::Instance()->ErrorFile()->Handle()))),
      ssi_maps_(Teuchos::null),
      ssi_matrices_(Teuchos::null),
      ssi_vectors_(Teuchos::null),
      strategy_assemble_(Teuchos::null),
      strategy_contact_(Teuchos::null),
      strategy_convcheck_(Teuchos::null),
      strategy_equilibration_(Teuchos::null),
      strategy_manifold_meshtying_(Teuchos::null),
      strategy_meshtying_(Teuchos::null),
      timer_(Teuchos::rcp(new Teuchos::Time("SSI_Mono", true)))
{
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void SSI::SSIMono::ApplyContactToSubProblems()
{
  // uncomplete matrices; we need to do this here since in contact simulations the dofs that
  // interact with each other can change and thus the graph of the matrix can also change.
  ssi_matrices_->ScaTraMatrix()->UnComplete();
  ssi_matrices_->ScaTraStructureMatrix()->UnComplete();
  ssi_matrices_->StructureScaTraMatrix()->UnComplete();

  // add contributions
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
bool SSI::SSIMono::IsUncompleteOfMatricesNecessaryForMeshTying() const
{
  // check for first iteration in calculation of initial time derivative
  if (IterationCount() == 0 and Step() == 0 and !DoCalculateInitialPotentialField()) return true;

  if (IterationCount() <= 2)
  {
    // check for first iteration in standard Newton loop
    if (Step() == 1 and !DoCalculateInitialPotentialField()) return true;

    // check for first iterations in calculation of initial potential field
    if (Step() == 0 and DoCalculateInitialPotentialField()) return true;

    // check for first iteration in restart simulations
    if (IsRestart())
    {
      auto* problem = DRT::Problem::Instance();
      // restart based on time step
      if (Step() == problem->Restart() + 1) return true;
      // restart based on time
      if (Time() == problem->RestartTime() + Dt()) return true;
    }
  }

  // if we have at least one contact interface the dofs that are in contact can change and therefore
  // also the matrices have to be uncompleted
  return SSIInterfaceContact();
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void SSI::SSIMono::ApplyMeshtyingToSubProblems()
{
  TEUCHOS_FUNC_TIME_MONITOR("SSI mono: apply mesh tying");
  if (SSIInterfaceMeshtying())
  {
    // check if matrices are filled because they have to be for the below methods
    if (!ssi_matrices_->StructureScaTraMatrix()->Filled())
      ssi_matrices_->CompleteStructureScaTraMatrix();
    if (!ssi_matrices_->ScaTraStructureMatrix()->Filled())
      ssi_matrices_->CompleteScaTraStructureMatrix();

    strategy_meshtying_->ApplyMeshtyingToScatraStructure(
        ssi_matrices_->ScaTraStructureMatrix(), IsUncompleteOfMatricesNecessaryForMeshTying());

    strategy_meshtying_->ApplyMeshtyingToStructureMatrix(*ssi_matrices_->StructureMatrix(),
        StructureField()->SystemMatrix(), IsUncompleteOfMatricesNecessaryForMeshTying());

    strategy_meshtying_->ApplyMeshtyingToStructureScatra(
        ssi_matrices_->StructureScaTraMatrix(), IsUncompleteOfMatricesNecessaryForMeshTying());

    ssi_vectors_->StructureResidual()->Update(
        1.0, strategy_meshtying_->ApplyMeshtyingToStructureRHS(StructureField()->RHS()), 1.0);

    if (IsScaTraManifold())
    {
      if (!ssi_matrices_->ScaTraManifoldStructureMatrix()->Filled())
        ssi_matrices_->CompleteScaTraManifoldStructureMatrix();
      if (!manifoldscatraflux_->MatrixManifoldStructure()->Filled())
        manifoldscatraflux_->CompleteMatrixManifoldStructure();
      if (!manifoldscatraflux_->MatrixScaTraStructure()->Filled())
        manifoldscatraflux_->CompleteMatrixScaTraStructure();

      strategy_meshtying_->ApplyMeshtyingToScatraManifoldStructure(
          ssi_matrices_->ScaTraManifoldStructureMatrix(),
          IsUncompleteOfMatricesNecessaryForMeshTying());

      strategy_meshtying_->ApplyMeshtyingToScatraManifoldStructure(
          manifoldscatraflux_->MatrixManifoldStructure(),
          IsUncompleteOfMatricesNecessaryForMeshTying());

      strategy_meshtying_->ApplyMeshtyingToScatraStructure(
          manifoldscatraflux_->MatrixScaTraStructure(),
          IsUncompleteOfMatricesNecessaryForMeshTying());
    }
  }
  // copy the structure residual and matrix if we do not have a mesh tying problem
  else
  {
    ssi_vectors_->StructureResidual()->Update(1.0, *(StructureField()->RHS()), 1.0);
    ssi_matrices_->StructureMatrix()->Add(*StructureField()->SystemMatrix(), false, 1.0, 1.0);
  }
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SSIMono::ApplyManifoldMeshtying()
{
  if (!manifoldscatraflux_->SystemMatrixManifold()->Filled())
    manifoldscatraflux_->SystemMatrixManifold()->Complete();

  if (!ssi_matrices_->ScaTraManifoldStructureMatrix()->Filled())
    ssi_matrices_->CompleteScaTraManifoldStructureMatrix();

  if (!manifoldscatraflux_->MatrixManifoldStructure()->Filled())
    manifoldscatraflux_->CompleteMatrixManifoldStructure();

  if (!manifoldscatraflux_->MatrixScaTraManifold()->Filled())
    manifoldscatraflux_->CompleteMatrixScaTraManifold();

  if (!manifoldscatraflux_->MatrixManifoldScatra()->Filled())
    manifoldscatraflux_->CompleteMatrixManifoldScaTra();

  // apply mesh tying to...
  // manifold - manifold
  strategy_manifold_meshtying_->ApplyMeshtyingToManifoldMatrix(
      ssi_matrices_->ManifoldMatrix(), ScaTraManifold()->SystemMatrixOperator());
  strategy_manifold_meshtying_->ApplyMeshtyingToManifoldMatrix(
      ssi_matrices_->ManifoldMatrix(), manifoldscatraflux_->SystemMatrixManifold());

  // manifold - structure
  strategy_manifold_meshtying_->ApplyMeshtyingToManifoldStructureMatrix(
      ssi_matrices_->ScaTraManifoldStructureMatrix(),
      manifoldscatraflux_->MatrixManifoldStructure(),
      IsUncompleteOfMatricesNecessaryForMeshTying());

  // scatra - manifold
  strategy_manifold_meshtying_->ApplyMeshtyingToScatraManifoldMatrix(
      ssi_matrices_->ScaTraScaTraManifoldMatrix(), manifoldscatraflux_->MatrixScaTraManifold(),
      IsUncompleteOfMatricesNecessaryForMeshTying());

  // manifold - scatra
  strategy_manifold_meshtying_->ApplyMeshtyingToManifoldScatraMatrix(
      ssi_matrices_->ScaTraManifoldScaTraMatrix(), manifoldscatraflux_->MatrixManifoldScatra());

  // RHS
  strategy_manifold_meshtying_->ApplyMeshTyingToManifoldRHS(ssi_vectors_->ManifoldResidual());
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SSIMono::AssembleMatAndRHS()
{
  TEUCHOS_FUNC_TIME_MONITOR("SSI mono: assemble global system");

  AssembleMatScaTra();

  AssembleMatStructure();

  if (IsScaTraManifold()) AssembleMatScaTraManifold();

  // finalize global system matrix
  ssi_matrices_->SystemMatrix()->Complete();

  // assemble monolithic RHS
  strategy_assemble_->AssembleRHS(ssi_vectors_->Residual(), ssi_vectors_->ScatraResidual(),
      ssi_vectors_->StructureResidual(),
      IsScaTraManifold() ? ssi_vectors_->ManifoldResidual() : Teuchos::null);
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
      ssi_matrices_->SystemMatrix(), ssi_matrices_->ManifoldMatrix());

  // assemble scatra manifold-structure block into system matrix
  strategy_assemble_->AssembleScatramanifoldStructure(
      ssi_matrices_->SystemMatrix(), ssi_matrices_->ScaTraManifoldStructureMatrix());

  // assemble contributions from scatra - scatra manifold coupling: derivs. of scatra side w.r.t.
  // scatra side
  strategy_assemble_->AssembleScatraScatra(
      ssi_matrices_->SystemMatrix(), manifoldscatraflux_->SystemMatrixScaTra());

  // assemble contributions from scatra - scatra manifold coupling: derivs. of manifold side w.r.t.
  // scatra side
  strategy_assemble_->AssembleScatraScatramanifold(
      ssi_matrices_->SystemMatrix(), ssi_matrices_->ScaTraScaTraManifoldMatrix());

  // assemble contributions from scatra - scatra manifold coupling: derivs. of scatra side w.r.t.
  // manifold side
  strategy_assemble_->AssembleScatramanifoldScatra(
      ssi_matrices_->SystemMatrix(), ssi_matrices_->ScaTraManifoldScaTraMatrix());

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
  TEUCHOS_FUNC_TIME_MONITOR("SSI mono: evaluate sub problems");

  // clear all matrices and residuals from previous Newton iteration
  ssi_matrices_->ClearMatrices();
  ssi_vectors_->ClearResiduals();

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

  // apply mesh tying to manifold domains
  if (IsScaTraManifold()) ApplyManifoldMeshtying();

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
    case CORE::LINALG::MatrixType::block_condition:
    case CORE::LINALG::MatrixType::block_condition_dof:
    {
      // equip smoother for scatra matrix blocks with null space
      ScaTraField()->BuildBlockNullSpaces(
          solver_, ssi_maps_->GetBlockPositions(Subproblem::scalar_transport)->at(0));
      if (IsScaTraManifold())
      {
        ScaTraManifold()->BuildBlockNullSpaces(
            solver_, ssi_maps_->GetBlockPositions(Subproblem::manifold)->at(0));
      }
      break;
    }

    case CORE::LINALG::MatrixType::sparse:
    {
      // equip smoother for scatra matrix block with empty parameter sub lists to trigger null space
      // computation
      std::ostringstream scatrablockstr;
      scatrablockstr << ssi_maps_->GetBlockPositions(Subproblem::scalar_transport)->at(0) + 1;
      Teuchos::ParameterList& blocksmootherparamsscatra =
          solver_->Params().sublist("Inverse" + scatrablockstr.str());
      blocksmootherparamsscatra.sublist("Belos Parameters");
      blocksmootherparamsscatra.sublist("MueLu Parameters");

      // equip smoother for scatra matrix block with null space associated with all degrees of
      // freedom on scatra discretization
      ScaTraField()->Discretization()->ComputeNullSpaceIfNecessary(blocksmootherparamsscatra);

      if (IsScaTraManifold())
      {
        std::ostringstream scatramanifoldblockstr;
        scatramanifoldblockstr << ssi_maps_->GetBlockPositions(Subproblem::manifold)->at(0) + 1;
        Teuchos::ParameterList& blocksmootherparamsscatramanifold =
            solver_->Params().sublist("Inverse" + scatramanifoldblockstr.str());
        blocksmootherparamsscatramanifold.sublist("Belos Parameters");
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
  iblockstr << ssi_maps_->GetBlockPositions(Subproblem::structure)->at(0) + 1;

  // equip smoother for structural matrix block with empty parameter sub lists to trigger null space
  // computation
  Teuchos::ParameterList& blocksmootherparams =
      solver_->Params().sublist("Inverse" + iblockstr.str());
  blocksmootherparams.sublist("Belos Parameters");
  blocksmootherparams.sublist("MueLu Parameters");

  // equip smoother for structural matrix block with null space associated with all degrees of
  // freedom on structural discretization
  StructureField()->Discretization()->ComputeNullSpaceIfNecessary(blocksmootherparams);
}  // SSI::SSIMono::BuildNullSpaces

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SSIMono::CompleteSubproblemMatrices()
{
  ssi_matrices_->ScaTraMatrix()->Complete();
  ssi_matrices_->CompleteScaTraStructureMatrix();
  ssi_matrices_->CompleteStructureScaTraMatrix();
  ssi_matrices_->StructureMatrix()->Complete();

  if (IsScaTraManifold())
  {
    ssi_matrices_->ManifoldMatrix()->Complete();
    ssi_matrices_->CompleteScaTraManifoldStructureMatrix();
    ssi_matrices_->CompleteScaTraScaTraManifoldMatrix();
    ssi_matrices_->CompleteScaTraManifoldScaTraMatrix();

    manifoldscatraflux_->CompleteMatrixManifoldScaTra();
    manifoldscatraflux_->CompleteMatrixManifoldStructure();
    manifoldscatraflux_->CompleteMatrixScaTraManifold();
    manifoldscatraflux_->CompleteMatrixScaTraStructure();
    manifoldscatraflux_->CompleteSystemMatrixManifold();
    manifoldscatraflux_->CompleteSystemMatrixScaTra();
  }
}

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

  if (print_matlab_) PrintSystemMatrixRHSToMatLabFormat();
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
void SSI::SSIMono::PrepareTimeLoop()
{
  SetStructSolution(
      StructureField()->Dispnp(), StructureField()->Velnp(), IsS2IKineticsWithPseudoContact());

  // calculate initial potential field if needed
  if (DoCalculateInitialPotentialField()) CalcInitialPotentialField();

  // calculate initial time derivatives
  CalcInitialTimeDerivative();

  ScaTraField()->Output();
  if (IsScaTraManifold()) ScaTraManifold()->Output();
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SSIMono::PrepareTimeStep()
{
  // update time and time step
  IncrementTimeAndStep();

  // pass structural degrees of freedom to scalar transport discretization
  SetStructSolution(
      StructureField()->Dispnp(), StructureField()->Velnp(), IsS2IKineticsWithPseudoContact());

  // prepare time step for scalar transport field
  ScaTraField()->PrepareTimeStep();
  if (IsScaTraManifold()) ScaTraManifold()->PrepareTimeStep();

  // if adaptive time stepping and different time step size: calculate time step in scatra
  // (PrepareTimeStep() of Scatra) and pass to other fields
  if (ScaTraField()->TimeStepAdapted()) SetDtFromScaTraToSSI();

  // pass scalar transport degrees of freedom to structural discretization
  // has to be called AFTER ScaTraField()->PrepareTimeStep() to ensure
  // consistent scalar transport state vector with valid Dirichlet conditions
  SetScatraSolution(ScaTraField()->Phinp());
  if (IsScaTraManifold()) SetScatraManifoldSolution(ScaTraManifold()->Phinp());

  // evaluate temperature from function and set to structural discretization
  EvaluateAndSetTemperatureField();

  // prepare time step for structural field
  StructureField()->PrepareTimeStep();

  // StructureField()->PrepareTimeStep() evaluates the DBC displaements on the master side. Now, the
  // master side displacements are copied to slave side to consider non zero DBC values in the first
  // Newton step on the slave side in case of interface mesh tying
  if (SSIInterfaceMeshtying())
  {
    for (const auto& meshtying : SSIStructureMeshTying()->MeshTyingHandlers())
    {
      auto coupling_adapter = meshtying->SlaveMasterCoupling();
      auto coupling_map_extractor = meshtying->SlaveMasterExtractor();

      // displacements
      coupling_map_extractor->InsertVector(
          coupling_adapter->MasterToSlave(
              coupling_map_extractor->ExtractVector(StructureField()->Dispnp(), 2)),
          1, StructureField()->WriteAccessDispnp());
      StructureField()->SetState(StructureField()->WriteAccessDispnp());
    }
  }

  PrintTimeStepInfo();
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
  const auto ssi_params = DRT::Problem::Instance()->SSIControlParams();

  const bool calc_initial_pot_elch =
      DRT::INPUT::IntegralValue<bool>(DRT::Problem::Instance()->ELCHControlParams(), "INITPOTCALC");
  const bool calc_initial_pot_ssi =
      DRT::INPUT::IntegralValue<bool>(ssi_params.sublist("ELCH"), "INITPOTCALC");

  if (ScaTraField()->EquilibrationMethod() != CORE::LINALG::EquilibrationMethod::none)
  {
    dserror(
        "You are within the monolithic solid scatra interaction framework but activated a pure "
        "scatra equilibration method. Delete this from 'SCALAR TRANSPORT DYNAMIC' section and set "
        "it in 'SSI CONTROL/MONOLITHIC' instead.");
  }
  if (equilibration_method_.global != CORE::LINALG::EquilibrationMethod::local and
      (equilibration_method_.structure != CORE::LINALG::EquilibrationMethod::none or
          equilibration_method_.scatra != CORE::LINALG::EquilibrationMethod::none))
    dserror("Either global equilibration or local equilibration");

  if (matrixtype_ == CORE::LINALG::MatrixType::sparse and
      (equilibration_method_.structure != CORE::LINALG::EquilibrationMethod::none or
          equilibration_method_.scatra != CORE::LINALG::EquilibrationMethod::none))
    dserror("Block based equilibration only for block matrices");

  if (!DRT::INPUT::IntegralValue<int>(
          DRT::Problem::Instance()->ScalarTransportDynamicParams(), "SKIPINITDER"))
  {
    dserror(
        "Initial derivatives are already calculated in monolithic SSI. Enable 'SKIPINITDER' in the "
        "input file.");
  }

  if (calc_initial_pot_elch)
    dserror("Initial potential is calculated by SSI. Disable in Elch section.");
  if (calc_initial_pot_ssi and Teuchos::getIntegralValue<INPAR::SSI::ScaTraTimIntType>(ssi_params,
                                   "SCATRATIMINTTYPE") != INPAR::SSI::ScaTraTimIntType::elch)
    dserror("Calculation of initial potential only in case of Elch");

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
  SSI::SSIBase::SetupSystem();

  // setup the ssi maps object
  ssi_maps_ = Teuchos::rcp(new SSI::UTILS::SSIMaps(*this));

  // perform initializations associated with global system matrix
  switch (matrixtype_)
  {
    case CORE::LINALG::MatrixType::block_field:
    {
      // safety check
      if (!solver_->Params().isSublist("AMGnxn Parameters"))
        dserror("Global system matrix with block structure requires AMGnxn block preconditioner!");

      // feed AMGnxn block preconditioner with null space information for each block of global
      // block system matrix
      BuildNullSpaces();

      break;
    }

    case CORE::LINALG::MatrixType::sparse:
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
  ssi_matrices_ = Teuchos::rcp(new SSI::UTILS::SSIMatrices(
      ssi_maps_, matrixtype_, ScaTraField()->MatrixType(), IsScaTraManifold()));

  // initialize residual and increment vectors
  ssi_vectors_ = Teuchos::rcp(new SSI::UTILS::SSIVectors(ssi_maps_, IsScaTraManifold()));

  // initialize strategy for assembly
  strategy_assemble_ = SSI::BuildAssembleStrategy(
      ssi_maps_, IsScaTraManifold(), matrixtype_, ScaTraField()->MatrixType());

  if (IsScaTraManifold())
  {
    // initialize object, that performs evaluations of OD coupling
    scatrastructureOffDiagcoupling_ = Teuchos::rcp(new SSI::ScatraManifoldStructureOffDiagCoupling(
        BlockMapStructure(), SSIMaps()->StructureDofRowMap(), SSIStructureMeshTying(),
        MeshtyingStrategyS2I(), ScaTraField(), ScaTraManifold(), StructureField()));

    // initialize object, that performs evaluations of scatra - scatra on manifold coupling
    manifoldscatraflux_ = Teuchos::rcp(new SSI::ScaTraManifoldScaTraFluxEvaluator(*this));

    // initialize object, that performs meshtying between manifold domains
    strategy_manifold_meshtying_ =
        SSI::BuildManifoldMeshTyingStrategy(ScaTraManifold()->Discretization(), ssi_maps_,
            IsScaTraManifoldMeshtying(), ScaTraManifold()->MatrixType());
  }
  else
  {
    scatrastructureOffDiagcoupling_ = Teuchos::rcp(new SSI::ScatraStructureOffDiagCoupling(
        BlockMapStructure(), SSIMaps()->StructureDofRowMap(), SSIStructureMeshTying(),
        MeshtyingStrategyS2I(), ScaTraField(), StructureField()));
  }
  // instantiate appropriate equilibration class
  strategy_equilibration_ = CORE::LINALG::BuildEquilibration(
      matrixtype_, GetBlockEquilibration(), MapsSubProblems()->FullMap());

  // instantiate appropriate contact class
  strategy_contact_ =
      SSI::BuildContactStrategy(CoNitscheStrategySsi(), ssi_maps_, ScaTraField()->MatrixType());

  // instantiate appropriate mesh tying class
  strategy_meshtying_ = SSI::BuildMeshtyingStrategy(
      IsScaTraManifold(), ScaTraField()->MatrixType(), ssi_maps_, SSIStructureMeshTying());

  // instantiate Dirichlet boundary condition handler class
  dbc_handler_ = SSI::BuildDBCHandler(IsScaTraManifold(), matrixtype_, ScaTraField(),
      IsScaTraManifold() ? ScaTraManifold() : Teuchos::null, ssi_maps_, StructureField());
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void SSI::SSIMono::SetScatraSolution(Teuchos::RCP<const Epetra_Vector> phi) const
{
  // call base class
  SSIBase::SetScatraSolution(phi);

  // set state for contact evaluation
  if (contact_strategy_nitsche_ != Teuchos::null) SetSSIContactStates(phi);
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void SSI::SSIMono::SetSSIContactStates(Teuchos::RCP<const Epetra_Vector> phi) const
{
  contact_strategy_nitsche_->SetState(MORTAR::state_scalar, *phi);
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void SSI::SSIMono::SolveLinearSystem()
{
  TEUCHOS_FUNC_TIME_MONITOR("SSI mono: solve linear system");
  strategy_equilibration_->EquilibrateSystem(
      ssi_matrices_->SystemMatrix(), ssi_vectors_->Residual(), BlockMapSystemMatrix());

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
  TEUCHOS_FUNC_TIME_MONITOR("SSI mono: solve Newton loop");
  // reset counter for Newton-Raphson iteration
  ResetIterationCount();

  // start Newton-Raphson iteration
  while (true)
  {
    // update iteration counter
    IncrementIterationCount();

    timer_->reset();

    // store time before evaluating elements and assembling global system of equations
    const double time_before_evaluate = timer_->wallTime();

    // set solution from last Newton step to all fields
    DistributeSolutionAllFields();

    // evaluate sub problems and get all matrices and right-hand-sides
    EvaluateSubproblems();

    // complete the sub problem matrices
    CompleteSubproblemMatrices();

    // assemble global system of equations
    AssembleMatAndRHS();

    // apply the Dirichlet boundary conditions to global system
    ApplyDBCToSystem();

    // time needed for evaluating elements and assembling global system of equations
    double my_evaluation_time = timer_->wallTime() - time_before_evaluate;
    Comm().MaxAll(&my_evaluation_time, &dt_eval_, 1);

    // safety check
    if (!ssi_matrices_->SystemMatrix()->Filled())
      dserror("Complete() has not been called on global system matrix yet!");

    // check termination criterion for Newton-Raphson iteration
    if (strategy_convcheck_->ExitNewtonRaphson(*this)) break;

    // clear the global increment vector
    ssi_vectors_->ClearIncrement();

    // store time before solving global system of equations
    const double time_before_solving = timer_->wallTime();

    SolveLinearSystem();

    // time needed for solving global system of equations
    double my_solve_time = timer_->wallTime() - time_before_solving;
    Comm().MaxAll(&my_solve_time, &dt_solve_, 1);

    // output performance statistics associated with linear solver into text file if
    // applicable
    if (DRT::INPUT::IntegralValue<bool>(
            *ScaTraField()->ScatraParameterList(), "OUTPUTLINSOLVERSTATS"))
      ScaTraField()->OutputLinSolverStats(*solver_, dt_solve_, Step(), IterationCount(),
          ssi_vectors_->Residual()->Map().NumGlobalElements());

    // update states for next Newton iteration
    UpdateIterScaTra();
    UpdateIterStructure();
  }
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SSIMono::Timeloop()
{
  if (Step() == 0) PrepareTimeLoop();

  // time loop
  while (NotFinished() and ScaTraField()->NotFinished())
  {
    TEUCHOS_FUNC_TIME_MONITOR("SSI mono: solve time step");
    // prepare time step
    PrepareTimeStep();

    // store time before calling nonlinear solver
    const double time = timer_->wallTime();

    // evaluate time step
    NewtonLoop();

    // determine time spent by nonlinear solver and take maximum over all processors via
    // communication
    double mydtnonlinsolve(timer_->wallTime() - time), dtnonlinsolve(0.);
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
  strategy_convcheck_->PrintNonConvergedSteps(Comm().MyPID());
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
      ssi_vectors_->Increment(), UTILS::SSIMaps::GetProblemPosition(Subproblem::scalar_transport)));
  ScaTraField()->ComputeIntermediateValues();

  if (IsScaTraManifold())
  {
    auto increment_manifold = MapsSubProblems()->ExtractVector(
        ssi_vectors_->Increment(), UTILS::SSIMaps::GetProblemPosition(Subproblem::manifold));

    // reconstruct slave side solution from master side
    if (IsScaTraManifoldMeshtying())
    {
      for (const auto& meshtying :
          strategy_manifold_meshtying_->SSIMeshTying()->MeshTyingHandlers())
      {
        auto coupling_adapter = meshtying->SlaveMasterCoupling();
        auto multimap = meshtying->SlaveMasterExtractor();

        auto master_dofs = multimap->ExtractVector(increment_manifold, 2);
        auto master_dofs_to_slave = coupling_adapter->MasterToSlave(master_dofs);
        multimap->InsertVector(master_dofs_to_slave, 1, increment_manifold);
      }
    }

    ScaTraManifold()->UpdateIter(increment_manifold);
    ScaTraManifold()->ComputeIntermediateValues();
  }
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void SSI::SSIMono::UpdateIterStructure()
{
  // set up structural increment vector
  const Teuchos::RCP<Epetra_Vector> increment_structure = MapsSubProblems()->ExtractVector(
      ssi_vectors_->Increment(), UTILS::SSIMaps::GetProblemPosition(Subproblem::structure));

  // consider structural meshtying. Copy master increments and displacements to slave side.
  if (SSIInterfaceMeshtying())
  {
    for (const auto& meshtying : SSIStructureMeshTying()->MeshTyingHandlers())
    {
      auto coupling_adapter = meshtying->SlaveMasterCoupling();
      auto coupling_map_extractor = meshtying->SlaveMasterExtractor();

      // displacements
      coupling_map_extractor->InsertVector(
          coupling_adapter->MasterToSlave(
              coupling_map_extractor->ExtractVector(StructureField()->Dispnp(), 2)),
          1, StructureField()->WriteAccessDispnp());
      StructureField()->SetState(StructureField()->WriteAccessDispnp());

      // increment
      coupling_map_extractor->InsertVector(
          coupling_adapter->MasterToSlave(
              coupling_map_extractor->ExtractVector(increment_structure, 2)),
          1, increment_structure);
    }
  }

  // update displacement of structure field
  StructureField()->UpdateStateIncrementally(increment_structure);
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
std::vector<CORE::LINALG::EquilibrationMethod> SSI::SSIMono::GetBlockEquilibration()
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
               equilibration_method_.scatra == CORE::LINALG::EquilibrationMethod::none)
      {
        equilibration_method_vector = std::vector<CORE::LINALG::EquilibrationMethod>(
            1, CORE::LINALG::EquilibrationMethod::none);
      }
      else
      {
        auto block_positions_scatra = ssi_maps_->GetBlockPositions(Subproblem::scalar_transport);
        auto block_position_structure = ssi_maps_->GetBlockPositions(Subproblem::structure);
        auto block_positions_scatra_manifold =
            IsScaTraManifold() ? ssi_maps_->GetBlockPositions(Subproblem::manifold) : Teuchos::null;

        equilibration_method_vector = std::vector<CORE::LINALG::EquilibrationMethod>(
            block_positions_scatra->size() + block_position_structure->size() +
                (IsScaTraManifold() ? block_positions_scatra_manifold->size() : 0),
            CORE::LINALG::EquilibrationMethod::none);

        for (const int block_position_scatra : *block_positions_scatra)
          equilibration_method_vector.at(block_position_scatra) = equilibration_method_.scatra;

        equilibration_method_vector.at(block_position_structure->at(0)) =
            equilibration_method_.structure;

        if (IsScaTraManifold())
        {
          for (const int block_position_scatra_manifold : *block_positions_scatra_manifold)
          {
            equilibration_method_vector.at(block_position_scatra_manifold) =
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

  ssi_vectors_->ManifoldResidual()->Update(1.0, *ScaTraManifold()->Residual(), 1.0);

  // evaluate coupling fluxes
  manifoldscatraflux_->Evaluate();

  ssi_vectors_->ManifoldResidual()->Update(1.0, *manifoldscatraflux_()->RHSManifold(), 1.0);
  ssi_vectors_->ScatraResidual()->Update(1.0, *manifoldscatraflux_()->RHSScaTra(), 1.0);
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void SSI::SSIMono::PrepareOutput()
{
  constexpr bool force_prepare = false;
  StructureField()->PrepareOutput(force_prepare);

  // prepare output of coupling sctra manifold - scatra
  if (IsScaTraManifold() and manifoldscatraflux_->DoOutput())
  {
    DistributeSolutionAllFields();
    manifoldscatraflux_->EvaluateScaTraManifoldInflow();
  }
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void SSI::SSIMono::DistributeSolutionAllFields(const bool restore_velocity)
{
  // has to be called before the call of 'SetStructSolution()' to have updated stress/strain states
  if (IsS2IKineticsWithPseudoContact()) StructureField()->DetermineStressStrain();

  // clear all states before redistributing the new states
  StructureField()->Discretization()->ClearState(true);
  ScaTraField()->Discretization()->ClearState(true);
  if (IsScaTraManifold()) ScaTraManifold()->Discretization()->ClearState(true);

  // needed to communicate to NOX state
  if (restore_velocity)
  {
    auto vel_temp = *StructureField()->Velnp();
    StructureField()->SetState(StructureField()->WriteAccessDispnp());
    StructureField()->WriteAccessVelnp()->Update(1.0, vel_temp, 0.0);
  }
  else
    StructureField()->SetState(StructureField()->WriteAccessDispnp());

  // distribute states to other fields
  SetStructSolution(
      StructureField()->Dispnp(), StructureField()->Velnp(), IsS2IKineticsWithPseudoContact());
  SetScatraSolution(ScaTraField()->Phinp());
  if (IsScaTraManifold()) SetScatraManifoldSolution(ScaTraManifold()->Phinp());
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void SSI::SSIMono::CalcInitialPotentialField()
{
  const auto equpot = DRT::INPUT::IntegralValue<INPAR::ELCH::EquPot>(
      DRT::Problem::Instance()->ELCHControlParams(), "EQUPOT");
  if (equpot != INPAR::ELCH::equpot_divi and equpot != INPAR::ELCH::equpot_enc_pde and
      equpot != INPAR::ELCH::equpot_enc_pde_elim)
  {
    dserror(
        "Initial potential field cannot be computed for chosen closing equation for electric "
        "potential!");
  }

  // store initial velocity to restore them afterwards
  auto init_velocity = *StructureField()->Velnp();

  // cast scatra time integrators to elch to call elch specific methods
  auto scatra_elch = Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntElch>(ScaTraField());
  auto manifold_elch = IsScaTraManifold()
                           ? Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntElch>(ScaTraManifold())
                           : Teuchos::null;
  if (scatra_elch == Teuchos::null or (IsScaTraManifold() and manifold_elch == Teuchos::null))
    dserror("Cast to Elch time integrator faild. Scatra is not an Elch problem");

  // prepare specific time integrators
  scatra_elch->PreCalcInitialPotentialField();
  if (IsScaTraManifold()) manifold_elch->PreCalcInitialPotentialField();

  auto scatra_elch_splitter = ScaTraField()->Splitter();
  auto manifold_elch_splitter = IsScaTraManifold() ? ScaTraManifold()->Splitter() : Teuchos::null;

  ResetIterationCount();

  while (true)
  {
    IncrementIterationCount();

    timer_->reset();

    // store time before evaluating elements and assembling global system of equations
    const double time_before_evaluate = timer_->wallTime();

    // prepare full SSI system
    DistributeSolutionAllFields(true);
    EvaluateSubproblems();

    // complete the sub problem matrices
    CompleteSubproblemMatrices();

    AssembleMatAndRHS();
    ApplyDBCToSystem();

    // apply artificial Dirichlet boundary conditions to system of equations (on concentration
    // dofs and on structure dofs)
    Teuchos::RCP<Epetra_Map> pseudo_dbc_map;
    if (IsScaTraManifold())
    {
      auto conc_map = CORE::LINALG::MergeMap(
          scatra_elch_splitter->OtherMap(), manifold_elch_splitter->OtherMap());
      pseudo_dbc_map = CORE::LINALG::MergeMap(conc_map, StructureField()->DofRowMap());
    }
    else
    {
      pseudo_dbc_map =
          CORE::LINALG::MergeMap(scatra_elch_splitter->OtherMap(), StructureField()->DofRowMap());
    }

    auto dbc_zeros = Teuchos::rcp(new Epetra_Vector(*pseudo_dbc_map, true));

    auto rhs = ssi_vectors_->Residual();
    CORE::LINALG::ApplyDirichletToSystem(*ssi_matrices_->SystemMatrix(), *ssi_vectors_->Increment(),
        *rhs, *dbc_zeros, *pseudo_dbc_map);
    ssi_vectors_->Residual()->Update(1.0, *rhs, 0.0);

    // time needed for evaluating elements and assembling global system of equations
    double my_evaluation_time = timer_->wallTime() - time_before_evaluate;
    Comm().MaxAll(&my_evaluation_time, &dt_eval_, 1);

    if (strategy_convcheck_->ExitNewtonRaphsonInitPotCalc(*this)) break;

    // solve for potential increments
    ssi_vectors_->ClearIncrement();

    // store time before solving global system of equations
    const double time_before_solving = timer_->wallTime();

    SolveLinearSystem();

    // time needed for solving global system of equations
    double my_solve_time = timer_->wallTime() - time_before_solving;
    Comm().MaxAll(&my_solve_time, &dt_solve_, 1);

    // update potential dofs in scatra and manifold fields
    UpdateIterScaTra();

    // copy initial state vector
    ScaTraField()->Phin()->Update(1.0, *ScaTraField()->Phinp(), 0.0);
    if (IsScaTraManifold()) ScaTraManifold()->Phin()->Update(1.0, *ScaTraManifold()->Phinp(), 0.0);

    // update state vectors for intermediate time steps (only for generalized alpha)
    ScaTraField()->ComputeIntermediateValues();
    if (IsScaTraManifold()) ScaTraManifold()->ComputeIntermediateValues();
  }

  scatra_elch->PostCalcInitialPotentialField();
  if (IsScaTraManifold()) manifold_elch->PostCalcInitialPotentialField();

  StructureField()->WriteAccessVelnp()->Update(1.0, init_velocity, 0.0);
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void SSI::SSIMono::CalcInitialTimeDerivative()
{
  // store initial velocity to restore them afterwards
  auto init_velocity = *StructureField()->Velnp();

  const bool is_elch = IsElchScaTraTimIntType();

  // prepare specific time integrators
  ScaTraField()->PreCalcInitialTimeDerivative();
  if (IsScaTraManifold()) ScaTraManifold()->PreCalcInitialTimeDerivative();

  auto scatra_elch_splitter = is_elch ? ScaTraField()->Splitter() : Teuchos::null;
  auto manifold_elch_splitter =
      (is_elch and IsScaTraManifold()) ? ScaTraManifold()->Splitter() : Teuchos::null;

  // initial screen output
  if (Comm().MyPID() == 0)
  {
    std::cout << "Calculating initial time derivative of state variables on discretization "
              << ScaTraField()->Discretization()->Name();
    if (IsScaTraManifold())
      std::cout << " and discretization " << ScaTraManifold()->Discretization()->Name();
    std::cout << std::endl;
  }

  // evaluate Dirichlet and Neumann boundary conditions
  ScaTraField()->ApplyBCToSystem();
  if (IsScaTraManifold()) ScaTraManifold()->ApplyBCToSystem();

  // clear history values (this is the first step)
  ScaTraField()->Hist()->PutScalar(0.0);
  if (IsScaTraManifold()) ScaTraManifold()->Hist()->PutScalar(0.0);

  // In a first step, we assemble the standard global system of equations (we need the residual)
  DistributeSolutionAllFields(true);
  EvaluateSubproblems();

  // complete the sub problem matrices
  CompleteSubproblemMatrices();

  AssembleMatAndRHS();
  ApplyDBCToSystem();

  // prepare mass matrices of sub problems and global system
  auto massmatrix_scatra =
      ScaTraField()->MatrixType() == CORE::LINALG::MatrixType::sparse
          ? Teuchos::rcp_dynamic_cast<CORE::LINALG::SparseOperator>(
                UTILS::SSIMatrices::SetupSparseMatrix(ScaTraField()->DofRowMap()))
          : Teuchos::rcp_dynamic_cast<CORE::LINALG::SparseOperator>(
                UTILS::SSIMatrices::SetupBlockMatrix(
                    ScaTraField()->BlockMaps(), ScaTraField()->BlockMaps()));

  auto massmatrix_manifold =
      IsScaTraManifold()
          ? (ScaTraManifold()->MatrixType() == CORE::LINALG::MatrixType::sparse
                    ? Teuchos::rcp_dynamic_cast<CORE::LINALG::SparseOperator>(
                          UTILS::SSIMatrices::SetupSparseMatrix(ScaTraManifold()->DofRowMap()))
                    : Teuchos::rcp_dynamic_cast<CORE::LINALG::SparseOperator>(
                          UTILS::SSIMatrices::SetupBlockMatrix(
                              ScaTraManifold()->BlockMaps(), ScaTraManifold()->BlockMaps())))
          : Teuchos::null;

  auto massmatrix_system = MatrixType() == CORE::LINALG::MatrixType::sparse
                               ? Teuchos::rcp_dynamic_cast<CORE::LINALG::SparseOperator>(
                                     UTILS::SSIMatrices::SetupSparseMatrix(DofRowMap()))
                               : Teuchos::rcp_dynamic_cast<CORE::LINALG::SparseOperator>(
                                     UTILS::SSIMatrices::SetupBlockMatrix(
                                         BlockMapSystemMatrix(), BlockMapSystemMatrix()));

  // fill ones on main diag of structure block (not solved)
  auto ones_struct = Teuchos::rcp(new Epetra_Vector(*StructureField()->DofRowMap(), true));
  ones_struct->PutScalar(1.0);
  MatrixType() == CORE::LINALG::MatrixType::sparse
      ? CORE::LINALG::InsertMyRowDiagonalIntoUnfilledMatrix(
            *CORE::LINALG::CastToSparseMatrixAndCheckSuccess(massmatrix_system), *ones_struct)
      : CORE::LINALG::InsertMyRowDiagonalIntoUnfilledMatrix(
            CORE::LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(massmatrix_system)
                ->Matrix(ssi_maps_->GetBlockPositions(Subproblem::structure)->at(0),
                    ssi_maps_->GetBlockPositions(Subproblem::structure)->at(0)),
            *ones_struct);

  // extract residuals of scatra and manifold from global residual
  auto rhs_scatra = Teuchos::rcp(new Epetra_Vector(*ScaTraField()->DofRowMap(), true));
  auto rhs_manifold = IsScaTraManifold()
                          ? Teuchos::rcp(new Epetra_Vector(*ScaTraManifold()->DofRowMap(), true))
                          : Teuchos::null;

  rhs_scatra->Update(1.0,
      *MapsSubProblems()->ExtractVector(ssi_vectors_->Residual(),
          UTILS::SSIMaps::GetProblemPosition(Subproblem::scalar_transport)),
      0.0);
  if (IsScaTraManifold())
  {
    rhs_manifold->Update(1.0,
        *MapsSubProblems()->ExtractVector(
            ssi_vectors_->Residual(), UTILS::SSIMaps::GetProblemPosition(Subproblem::manifold)),
        0.0);
  }

  // In a second step, we need to modify the assembled system of equations, since we want to solve
  // M phidt^0 = f^n - K\phi^n - C(u_n)\phi^n
  // In particular, we need to replace the global system matrix by a global mass matrix,
  // and we need to remove all transient contributions associated with time discretization from the
  // global residual vector.

  // Evaluate mass matrix and modify residual
  ScaTraField()->EvaluateInitialTimeDerivative(massmatrix_scatra, rhs_scatra);
  if (IsScaTraManifold())
    ScaTraManifold()->EvaluateInitialTimeDerivative(massmatrix_manifold, rhs_manifold);

  // assemble global mass matrix from
  switch (MatrixType())
  {
    case CORE::LINALG::MatrixType::block_field:
    {
      switch (ScaTraField()->MatrixType())
      {
        case CORE::LINALG::MatrixType::block_condition:
        case CORE::LINALG::MatrixType::block_condition_dof:
        {
          auto massmatrix_system_block =
              CORE::LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(massmatrix_system);

          auto massmatrix_scatra_block =
              CORE::LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(massmatrix_scatra);

          auto positions_scatra = ssi_maps_->GetBlockPositions(Subproblem::scalar_transport);

          for (int i = 0; i < static_cast<int>(positions_scatra->size()); ++i)
          {
            const int position_scatra = positions_scatra->at(i);
            massmatrix_system_block->Matrix(position_scatra, position_scatra)
                .Add(massmatrix_scatra_block->Matrix(i, i), false, 1.0, 1.0);
          }
          if (IsScaTraManifold())
          {
            auto positions_manifold = ssi_maps_->GetBlockPositions(Subproblem::manifold);

            auto massmatrix_manifold_block =
                CORE::LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(massmatrix_manifold);

            for (int i = 0; i < static_cast<int>(positions_manifold->size()); ++i)
            {
              const int position_manifold = positions_manifold->at(i);
              massmatrix_system_block->Matrix(position_manifold, position_manifold)
                  .Add(massmatrix_manifold_block->Matrix(i, i), false, 1.0, 1.0);
            }
          }

          break;
        }

        case CORE::LINALG::MatrixType::sparse:
        {
          auto massmatrix_system_block =
              CORE::LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(massmatrix_system);

          const int position_scatra =
              ssi_maps_->GetBlockPositions(Subproblem::scalar_transport)->at(0);

          massmatrix_system_block->Matrix(position_scatra, position_scatra)
              .Add(*CORE::LINALG::CastToSparseMatrixAndCheckSuccess(massmatrix_scatra), false, 1.0,
                  1.0);

          if (IsScaTraManifold())
          {
            const int position_manifold = ssi_maps_->GetBlockPositions(Subproblem::manifold)->at(0);

            massmatrix_system_block->Matrix(position_manifold, position_manifold)
                .Add(*CORE::LINALG::CastToSparseMatrixAndCheckSuccess(massmatrix_manifold), false,
                    1.0, 1.0);
          }
          break;
        }

        default:
        {
          dserror("Invalid matrix type associated with scalar transport field!");
          break;
        }
      }
      massmatrix_system->Complete();
      break;
    }
    case CORE::LINALG::MatrixType::sparse:
    {
      auto massmatrix_system_sparse =
          CORE::LINALG::CastToSparseMatrixAndCheckSuccess(massmatrix_system);
      massmatrix_system_sparse->Add(
          *CORE::LINALG::CastToSparseMatrixAndCheckSuccess(massmatrix_scatra), false, 1.0, 1.0);

      if (IsScaTraManifold())
        massmatrix_system_sparse->Add(
            *CORE::LINALG::CastToSparseMatrixAndCheckSuccess(massmatrix_manifold), false, 1.0, 1.0);

      massmatrix_system->Complete(*DofRowMap(), *DofRowMap());
      break;
    }
    default:
    {
      dserror("Type of global system matrix for scalar-structure interaction not recognized!");
      break;
    }
  }

  // reconstruct global residual from partial residuals
  auto rhs_system = Teuchos::RCP<Epetra_Vector>(new Epetra_Vector(*DofRowMap(), true));
  MapsSubProblems()->InsertVector(
      rhs_scatra, UTILS::SSIMaps::GetProblemPosition(Subproblem::scalar_transport), rhs_system);
  if (IsScaTraManifold())
    MapsSubProblems()->InsertVector(
        rhs_manifold, UTILS::SSIMaps::GetProblemPosition(Subproblem::manifold), rhs_system);

  // apply artificial Dirichlet boundary conditions to system of equations to non-transported
  // scalars and structure
  Teuchos::RCP<Epetra_Map> pseudo_dbc_map;
  if (IsScaTraManifold() and is_elch)
  {
    auto conc_map =
        CORE::LINALG::MergeMap(scatra_elch_splitter->CondMap(), manifold_elch_splitter->CondMap());
    pseudo_dbc_map = CORE::LINALG::MergeMap(conc_map, StructureField()->DofRowMap());
  }
  else if (is_elch)
  {
    pseudo_dbc_map =
        CORE::LINALG::MergeMap(scatra_elch_splitter->CondMap(), StructureField()->DofRowMap());
  }
  else
    pseudo_dbc_map = Teuchos::rcp(new Epetra_Map(*StructureField()->DofRowMap()));

  auto dbc_zeros = Teuchos::rcp(new Epetra_Vector(*pseudo_dbc_map, true));

  // temporal derivative of transported scalars
  auto phidtnp_system = Teuchos::RCP<Epetra_Vector>(new Epetra_Vector(*DofRowMap(), true));
  CORE::LINALG::ApplyDirichletToSystem(
      *massmatrix_system, *phidtnp_system, *rhs_system, *dbc_zeros, *(pseudo_dbc_map));

  // solve global system of equations for initial time derivative of state variables
  solver_->Solve(massmatrix_system->EpetraOperator(), phidtnp_system, rhs_system, true, true);

  // copy solution to sub problmes
  auto phidtnp_scatra = MapsSubProblems()->ExtractVector(
      phidtnp_system, UTILS::SSIMaps::GetProblemPosition(Subproblem::scalar_transport));
  ScaTraField()->Phidtnp()->Update(1.0, *phidtnp_scatra, 0.0);
  ScaTraField()->Phidtn()->Update(1.0, *phidtnp_scatra, 0.0);

  if (IsScaTraManifold())
  {
    auto phidtnp_manifold = MapsSubProblems()->ExtractVector(
        phidtnp_system, UTILS::SSIMaps::GetProblemPosition(Subproblem::manifold));
    ScaTraManifold()->Phidtnp()->Update(1.0, *phidtnp_manifold, 0.0);
    ScaTraManifold()->Phidtn()->Update(1.0, *phidtnp_manifold, 0.0);
  }

  // reset solver
  solver_->Reset();

  ScaTraField()->PostCalcInitialTimeDerivative();
  if (IsScaTraManifold()) ScaTraManifold()->PostCalcInitialTimeDerivative();

  StructureField()->WriteAccessVelnp()->Update(1.0, init_velocity, 0.0);
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> SSI::SSIMono::MapsSubProblems() const
{
  return ssi_maps_->MapsSubProblems();
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> SSI::SSIMono::BlockMapScaTra() const
{
  return ssi_maps_->BlockMapScaTra();
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> SSI::SSIMono::BlockMapScaTraManifold() const
{
  return ssi_maps_->BlockMapScaTraManifold();
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> SSI::SSIMono::BlockMapStructure() const
{
  return ssi_maps_->BlockMapStructure();
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> SSI::SSIMono::BlockMapSystemMatrix() const
{
  return ssi_maps_->BlockMapSystemMatrix();
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void SSI::SSIMono::PrintTimeStepInfo()
{
  if (Comm().MyPID() == 0)
  {
    std::cout << std::endl
              << "TIME: " << std::setw(11) << std::setprecision(4) << std::scientific << Time()
              << "/" << MaxTime() << "  DT = " << Dt() << "  STEP = " << Step() << "/" << NStep()
              << std::endl;
  }
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void SSI::SSIMono::PrintSystemMatrixRHSToMatLabFormat()
{
  // print system matrix
  switch (matrixtype_)
  {
    case CORE::LINALG::MatrixType::block_field:
    {
      auto block_matrix = CORE::LINALG::CastToConstBlockSparseMatrixBaseAndCheckSuccess(
          ssi_matrices_->SystemMatrix());

      for (int row = 0; row < block_matrix->Rows(); ++row)
      {
        for (int col = 0; col < block_matrix->Cols(); ++col)
        {
          std::ostringstream filename;
          filename << DRT::Problem::Instance()->OutputControlFile()->FileName()
                   << "_block_system_matrix_" << row << "_" << col << ".csv";

          CORE::LINALG::PrintMatrixInMatlabFormat(
              filename.str(), *block_matrix->Matrix(row, col).EpetraMatrix(), true);
        }
      }
      break;
    }

    case CORE::LINALG::MatrixType::sparse:
    {
      auto sparse_matrix =
          CORE::LINALG::CastToConstSparseMatrixAndCheckSuccess(ssi_matrices_->SystemMatrix());

      const std::string filename =
          DRT::Problem::Instance()->OutputControlFile()->FileName() + "_sparse_system_matrix.csv";

      CORE::LINALG::PrintMatrixInMatlabFormat(filename, *sparse_matrix->EpetraMatrix(), true);
      break;
    }

    default:
    {
      dserror("Type of global system matrix for scalar-structure interaction not recognized!");
      break;
    }
  }

  // print rhs
  {
    const std::string filename =
        DRT::Problem::Instance()->OutputControlFile()->FileName() + "_system_vector.csv";
    CORE::LINALG::PrintVectorInMatlabFormat(filename, *ssi_vectors_->Residual(), true);
  }

  // print full map
  {
    const std::string filename =
        DRT::Problem::Instance()->OutputControlFile()->FileName() + "_full_map.csv";
    CORE::LINALG::PrintMapInMatlabFormat(filename, *ssi_maps_->MapSystemMatrix(), true);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIMono::SetScatraManifoldSolution(Teuchos::RCP<const Epetra_Vector> phi)
{
  // scatra values on master side copied to manifold
  auto manifold_on_scatra =
      CORE::LINALG::CreateVector(*ScaTraField()->Discretization()->DofRowMap(), true);

  for (const auto& coup : manifoldscatraflux_->ScaTraManifoldCouplings())
  {
    auto manifold_cond = coup->ManifoldMapExtractor()->ExtractCondVector(*phi);
    auto manifold_on_scatra_cond = coup->CouplingAdapter()->SlaveToMaster(manifold_cond);
    coup->ScaTraMapExtractor()->InsertCondVector(manifold_on_scatra_cond, manifold_on_scatra);
  }
  ScaTraField()->Discretization()->SetState(0, "manifold_on_scatra", manifold_on_scatra);
}