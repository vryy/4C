/*--------------------------------------------------------------------------*/
/*! \file
\brief monolithic scalar-structure interaction

\level 2

*/
/*--------------------------------------------------------------------------*/
#include "4C_ssi_monolithic.hpp"

#include "4C_adapter_scatra_base_algorithm.hpp"
#include "4C_adapter_str_ssiwrapper.hpp"
#include "4C_adapter_str_structure_new.hpp"
#include "4C_contact_nitsche_strategy_ssi.hpp"
#include "4C_discretization_condition_locsys.hpp"
#include "4C_discretization_fem_general_assemblestrategy.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_ssi.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_equilibrate.hpp"
#include "4C_linalg_mapextractor.hpp"
#include "4C_linalg_matrixtransform.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_linalg_utils_sparse_algebra_print.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_scatra_timint_elch.hpp"
#include "4C_scatra_timint_meshtying_strategy_s2i.hpp"
#include "4C_ssi_contact_strategy.hpp"
#include "4C_ssi_coupling.hpp"
#include "4C_ssi_manifold_utils.hpp"
#include "4C_ssi_monolithic_assemble_strategy.hpp"
#include "4C_ssi_monolithic_convcheck_strategies.hpp"
#include "4C_ssi_monolithic_dbc_handler.hpp"
#include "4C_ssi_monolithic_evaluate_OffDiag.hpp"
#include "4C_ssi_monolithic_meshtying_strategy.hpp"
#include "4C_ssi_utils.hpp"
#include "4C_structure_new_model_evaluator_contact.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
SSI::SsiMono::SsiMono(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams)
    : SSIBase(comm, globaltimeparams),
      equilibration_method_{Teuchos::getIntegralValue<CORE::LINALG::EquilibrationMethod>(
                                globaltimeparams.sublist("MONOLITHIC"), "EQUILIBRATION"),
          Teuchos::getIntegralValue<CORE::LINALG::EquilibrationMethod>(
              globaltimeparams.sublist("MONOLITHIC"), "EQUILIBRATION_SCATRA"),
          Teuchos::getIntegralValue<CORE::LINALG::EquilibrationMethod>(
              globaltimeparams.sublist("MONOLITHIC"), "EQUILIBRATION_STRUCTURE")},
      matrixtype_(Teuchos::getIntegralValue<CORE::LINALG::MatrixType>(
          globaltimeparams.sublist("MONOLITHIC"), "MATRIXTYPE")),
      print_matlab_(CORE::UTILS::IntegralValue<bool>(
          globaltimeparams.sublist("MONOLITHIC"), "PRINT_MAT_RHS_MAP_MATLAB")),
      relax_lin_solver_tolerance_(
          globaltimeparams.sublist("MONOLITHIC").get<double>("RELAX_LIN_SOLVER_TOLERANCE")),
      relax_lin_solver_iter_step_(
          globaltimeparams.sublist("MONOLITHIC").get<int>("RELAX_LIN_SOLVER_STEP")),
      solver_(Teuchos::rcp(new CORE::LINALG::Solver(
          GLOBAL::Problem::Instance()->SolverParams(
              globaltimeparams.sublist("MONOLITHIC").get<int>("LINEAR_SOLVER")),
          comm))),
      timer_(Teuchos::rcp(new Teuchos::Time("SSI_Mono", true)))
{
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void SSI::SsiMono::apply_contact_to_sub_problems()
{
  // uncomplete matrices; we need to do this here since in contact simulations the dofs that
  // interact with each other can change and thus the graph of the matrix can also change.
  ssi_matrices_->ScaTraMatrix()->UnComplete();
  ssi_matrices_->sca_tra_structure_matrix()->UnComplete();
  ssi_matrices_->structure_sca_tra_matrix()->UnComplete();

  // add contributions
  strategy_contact_->apply_contact_to_scatra_residual(ssi_vectors_->ScatraResidual());
  strategy_contact_->apply_contact_to_scatra_scatra(ssi_matrices_->ScaTraMatrix());
  strategy_contact_->apply_contact_to_scatra_structure(ssi_matrices_->sca_tra_structure_matrix());
  strategy_contact_->apply_contact_to_structure_scatra(ssi_matrices_->structure_sca_tra_matrix());
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void SSI::SsiMono::apply_dbc_to_system()
{
  // apply Dirichlet boundary conditions to global system matrix
  dbc_handler_->apply_dbc_to_system_matrix(ssi_matrices_->SystemMatrix());

  // apply Dirichlet boundary conditions to global RHS
  dbc_handler_->ApplyDBCToRHS(ssi_vectors_->Residual());
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
bool SSI::SsiMono::is_uncomplete_of_matrices_necessary_for_mesh_tying() const
{
  // check for first iteration in calculation of initial time derivative
  if (IterationCount() == 0 and Step() == 0 and !do_calculate_initial_potential_field())
    return true;

  if (IterationCount() <= 2)
  {
    // check for first iteration in standard Newton loop
    if (Step() == 1 and !do_calculate_initial_potential_field()) return true;

    // check for first iterations in calculation of initial potential field
    if (Step() == 0 and do_calculate_initial_potential_field()) return true;

    // check for first iteration in restart simulations
    if (IsRestart())
    {
      auto* problem = GLOBAL::Problem::Instance();
      // restart based on time step
      if (Step() == problem->Restart() + 1) return true;
    }
  }

  // if we have at least one contact interface the dofs that are in contact can change and therefore
  // also the matrices have to be uncompleted
  return SSIInterfaceContact();
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void SSI::SsiMono::apply_meshtying_to_sub_problems()
{
  TEUCHOS_FUNC_TIME_MONITOR("SSI mono: apply mesh tying");
  if (ssi_interface_meshtying())
  {
    // check if matrices are filled because they have to be for the below methods
    if (!ssi_matrices_->structure_sca_tra_matrix()->Filled())
      ssi_matrices_->complete_structure_sca_tra_matrix();
    if (!ssi_matrices_->sca_tra_structure_matrix()->Filled())
      ssi_matrices_->complete_sca_tra_structure_matrix();

    strategy_meshtying_->apply_meshtying_to_scatra_structure(
        ssi_matrices_->sca_tra_structure_matrix(),
        is_uncomplete_of_matrices_necessary_for_mesh_tying());

    strategy_meshtying_->apply_meshtying_to_structure_matrix(*ssi_matrices_->StructureMatrix(),
        structure_field()->SystemMatrix(), is_uncomplete_of_matrices_necessary_for_mesh_tying());

    strategy_meshtying_->apply_meshtying_to_structure_scatra(
        ssi_matrices_->structure_sca_tra_matrix(),
        is_uncomplete_of_matrices_necessary_for_mesh_tying());

    ssi_vectors_->StructureResidual()->Update(
        1.0, strategy_meshtying_->apply_meshtying_to_structure_rhs(structure_field()->RHS()), 1.0);

    if (is_sca_tra_manifold())
    {
      if (!ssi_matrices_->sca_tra_manifold_structure_matrix()->Filled())
        ssi_matrices_->complete_sca_tra_manifold_structure_matrix();
      if (!manifoldscatraflux_->matrix_manifold_structure()->Filled())
        manifoldscatraflux_->complete_matrix_manifold_structure();
      if (!manifoldscatraflux_->matrix_sca_tra_structure()->Filled())
        manifoldscatraflux_->complete_matrix_sca_tra_structure();

      strategy_meshtying_->apply_meshtying_to_scatra_manifold_structure(
          ssi_matrices_->sca_tra_manifold_structure_matrix(),
          is_uncomplete_of_matrices_necessary_for_mesh_tying());

      strategy_meshtying_->apply_meshtying_to_scatra_manifold_structure(
          manifoldscatraflux_->matrix_manifold_structure(),
          is_uncomplete_of_matrices_necessary_for_mesh_tying());

      strategy_meshtying_->apply_meshtying_to_scatra_structure(
          manifoldscatraflux_->matrix_sca_tra_structure(),
          is_uncomplete_of_matrices_necessary_for_mesh_tying());
    }
  }
  // copy the structure residual and matrix if we do not have a mesh tying problem
  else
  {
    ssi_vectors_->StructureResidual()->Update(1.0, *(structure_field()->RHS()), 1.0);
    ssi_matrices_->StructureMatrix()->Add(*structure_field()->SystemMatrix(), false, 1.0, 1.0);
  }
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SsiMono::apply_manifold_meshtying()
{
  if (!manifoldscatraflux_->system_matrix_manifold()->Filled())
    manifoldscatraflux_->system_matrix_manifold()->Complete();

  if (!ssi_matrices_->sca_tra_manifold_structure_matrix()->Filled())
    ssi_matrices_->complete_sca_tra_manifold_structure_matrix();

  if (!manifoldscatraflux_->matrix_manifold_structure()->Filled())
    manifoldscatraflux_->complete_matrix_manifold_structure();

  if (!manifoldscatraflux_->matrix_sca_tra_manifold()->Filled())
    manifoldscatraflux_->complete_matrix_sca_tra_manifold();

  if (!manifoldscatraflux_->matrix_manifold_scatra()->Filled())
    manifoldscatraflux_->complete_matrix_manifold_sca_tra();

  // apply mesh tying to...
  // manifold - manifold
  strategy_manifold_meshtying_->apply_meshtying_to_manifold_matrix(
      ssi_matrices_->ManifoldMatrix(), ScaTraManifold()->system_matrix_operator());
  strategy_manifold_meshtying_->apply_meshtying_to_manifold_matrix(
      ssi_matrices_->ManifoldMatrix(), manifoldscatraflux_->system_matrix_manifold());

  // manifold - structure
  strategy_manifold_meshtying_->apply_meshtying_to_manifold_structure_matrix(
      ssi_matrices_->sca_tra_manifold_structure_matrix(),
      manifoldscatraflux_->matrix_manifold_structure(),
      is_uncomplete_of_matrices_necessary_for_mesh_tying());

  // scatra - manifold
  strategy_manifold_meshtying_->apply_meshtying_to_scatra_manifold_matrix(
      ssi_matrices_->sca_tra_sca_tra_manifold_matrix(),
      manifoldscatraflux_->matrix_sca_tra_manifold(),
      is_uncomplete_of_matrices_necessary_for_mesh_tying());

  // manifold - scatra
  strategy_manifold_meshtying_->apply_meshtying_to_manifold_scatra_matrix(
      ssi_matrices_->sca_tra_manifold_sca_tra_matrix(),
      manifoldscatraflux_->matrix_manifold_scatra());

  // RHS
  strategy_manifold_meshtying_->apply_mesh_tying_to_manifold_rhs(ssi_vectors_->ManifoldResidual());
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SsiMono::assemble_mat_and_rhs()
{
  TEUCHOS_FUNC_TIME_MONITOR("SSI mono: assemble global system");

  assemble_mat_sca_tra();

  assemble_mat_structure();

  if (is_sca_tra_manifold()) assemble_mat_sca_tra_manifold();

  // finalize global system matrix
  ssi_matrices_->SystemMatrix()->Complete();

  // assemble monolithic RHS
  strategy_assemble_->AssembleRHS(ssi_vectors_->Residual(), ssi_vectors_->ScatraResidual(),
      ssi_vectors_->StructureResidual(),
      is_sca_tra_manifold() ? ssi_vectors_->ManifoldResidual() : Teuchos::null);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SsiMono::assemble_mat_sca_tra()
{
  // assemble scatra-scatra block into system matrix
  strategy_assemble_->assemble_scatra_scatra(
      ssi_matrices_->SystemMatrix(), ssi_matrices_->ScaTraMatrix());

  // assemble scatra-structure block into system matrix
  strategy_assemble_->assemble_scatra_structure(
      ssi_matrices_->SystemMatrix(), ssi_matrices_->sca_tra_structure_matrix());
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SsiMono::assemble_mat_sca_tra_manifold()
{
  // assemble scatra manifold - scatra manifold block into system matrix
  strategy_assemble_->assemble_scatramanifold_scatramanifold(
      ssi_matrices_->SystemMatrix(), ssi_matrices_->ManifoldMatrix());

  // assemble scatra manifold-structure block into system matrix
  strategy_assemble_->assemble_scatramanifold_structure(
      ssi_matrices_->SystemMatrix(), ssi_matrices_->sca_tra_manifold_structure_matrix());

  // assemble contributions from scatra - scatra manifold coupling: derivs. of scatra side w.r.t.
  // scatra side
  strategy_assemble_->assemble_scatra_scatra(
      ssi_matrices_->SystemMatrix(), manifoldscatraflux_->SystemMatrixScaTra());

  // assemble contributions from scatra - scatra manifold coupling: derivs. of manifold side w.r.t.
  // scatra side
  strategy_assemble_->assemble_scatra_scatramanifold(
      ssi_matrices_->SystemMatrix(), ssi_matrices_->sca_tra_sca_tra_manifold_matrix());

  // assemble contributions from scatra - scatra manifold coupling: derivs. of scatra side w.r.t.
  // manifold side
  strategy_assemble_->assemble_scatramanifold_scatra(
      ssi_matrices_->SystemMatrix(), ssi_matrices_->sca_tra_manifold_sca_tra_matrix());

  strategy_assemble_->assemble_scatra_structure(
      ssi_matrices_->SystemMatrix(), manifoldscatraflux_->matrix_sca_tra_structure());
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SsiMono::assemble_mat_structure()
{  // assemble structure-scatra block into system matrix
  strategy_assemble_->assemble_structure_scatra(
      ssi_matrices_->SystemMatrix(), ssi_matrices_->structure_sca_tra_matrix());

  // assemble structure-structure block into system matrix
  strategy_assemble_->assemble_structure_structure(
      ssi_matrices_->SystemMatrix(), ssi_matrices_->StructureMatrix());
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SsiMono::evaluate_subproblems()
{
  TEUCHOS_FUNC_TIME_MONITOR("SSI mono: evaluate sub problems");

  // clear all matrices and residuals from previous Newton iteration
  ssi_matrices_->ClearMatrices();
  ssi_vectors_->ClearResiduals();

  // evaluate temperature from function and set to structural discretization
  evaluate_and_set_temperature_field();

  // build system matrix and residual for structure field
  structure_field()->Evaluate();

  // build system matrix and residual for scalar transport field
  evaluate_sca_tra();

  // build system matrix and residual for scalar transport field on manifold
  if (is_sca_tra_manifold()) evaluate_sca_tra_manifold();

  // build all off diagonal matrices
  evaluate_off_diag_contributions();

  // apply mesh tying to sub problems
  apply_meshtying_to_sub_problems();

  // apply mesh tying to manifold domains
  if (is_sca_tra_manifold()) apply_manifold_meshtying();

  // apply contact contributions to sub problems
  if (SSIInterfaceContact()) apply_contact_to_sub_problems();
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void SSI::SsiMono::evaluate_off_diag_contributions()
{
  // evaluate off-diagonal scatra-structure block (domain contributions) of global system matrix
  scatrastructure_off_diagcoupling_->evaluate_off_diag_block_scatra_structure_domain(
      ssi_matrices_->sca_tra_structure_matrix());

  // evaluate off-diagonal scatra-structure block (interface contributions) of global system matrix
  if (ssi_interface_meshtying())
    scatrastructure_off_diagcoupling_->evaluate_off_diag_block_scatra_structure_interface(
        ssi_matrices_->sca_tra_structure_matrix());

  // evaluate off-diagonal structure-scatra block (we only have domain contributions so far) of
  // global system matrix
  scatrastructure_off_diagcoupling_->evaluate_off_diag_block_structure_scatra_domain(
      ssi_matrices_->structure_sca_tra_matrix());

  if (is_sca_tra_manifold())
  {
    // evaluate off-diagonal manifold-structure block of global system matrix
    scatrastructure_off_diagcoupling_->evaluate_off_diag_block_scatra_manifold_structure_domain(
        ssi_matrices_->sca_tra_manifold_structure_matrix());
  }
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void SSI::SsiMono::build_null_spaces() const
{
  switch (ScaTraField()->MatrixType())
  {
    case CORE::LINALG::MatrixType::block_condition:
    case CORE::LINALG::MatrixType::block_condition_dof:
    {
      // equip smoother for scatra matrix blocks with null space
      ScaTraField()->build_block_null_spaces(
          solver_, ssi_maps_->GetBlockPositions(Subproblem::scalar_transport).at(0));
      if (is_sca_tra_manifold())
      {
        ScaTraManifold()->build_block_null_spaces(
            solver_, ssi_maps_->GetBlockPositions(Subproblem::manifold).at(0));
      }
      break;
    }

    case CORE::LINALG::MatrixType::sparse:
    {
      // equip smoother for scatra matrix block with empty parameter sub lists to trigger null space
      // computation
      std::ostringstream scatrablockstr;
      scatrablockstr << ssi_maps_->GetBlockPositions(Subproblem::scalar_transport).at(0) + 1;
      Teuchos::ParameterList& blocksmootherparamsscatra =
          solver_->Params().sublist("Inverse" + scatrablockstr.str());
      blocksmootherparamsscatra.sublist("Belos Parameters");
      blocksmootherparamsscatra.sublist("MueLu Parameters");

      // equip smoother for scatra matrix block with null space associated with all degrees of
      // freedom on scatra discretization
      ScaTraField()->discretization()->compute_null_space_if_necessary(blocksmootherparamsscatra);

      if (is_sca_tra_manifold())
      {
        std::ostringstream scatramanifoldblockstr;
        scatramanifoldblockstr << ssi_maps_->GetBlockPositions(Subproblem::manifold).at(0) + 1;
        Teuchos::ParameterList& blocksmootherparamsscatramanifold =
            solver_->Params().sublist("Inverse" + scatramanifoldblockstr.str());
        blocksmootherparamsscatramanifold.sublist("Belos Parameters");
        blocksmootherparamsscatramanifold.sublist("MueLu Parameters");

        // equip smoother for scatra matrix block with null space associated with all degrees of
        // freedom on scatra discretization
        ScaTraManifold()->discretization()->compute_null_space_if_necessary(
            blocksmootherparamsscatramanifold);
      }

      break;
    }

    default:
    {
      FOUR_C_THROW("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // store number of matrix block associated with structural field as string
  std::stringstream iblockstr;
  iblockstr << ssi_maps_->GetBlockPositions(Subproblem::structure).at(0) + 1;

  // equip smoother for structural matrix block with empty parameter sub lists to trigger null space
  // computation
  Teuchos::ParameterList& blocksmootherparams =
      solver_->Params().sublist("Inverse" + iblockstr.str());
  blocksmootherparams.sublist("Belos Parameters");
  blocksmootherparams.sublist("MueLu Parameters");

  // equip smoother for structural matrix block with null space associated with all degrees of
  // freedom on structural discretization
  structure_field()->discretization()->compute_null_space_if_necessary(blocksmootherparams);
}  // SSI::ssi_mono::build_null_spaces

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SsiMono::complete_subproblem_matrices()
{
  ssi_matrices_->ScaTraMatrix()->Complete();
  ssi_matrices_->complete_sca_tra_structure_matrix();
  ssi_matrices_->complete_structure_sca_tra_matrix();
  ssi_matrices_->StructureMatrix()->Complete();

  if (is_sca_tra_manifold())
  {
    ssi_matrices_->ManifoldMatrix()->Complete();
    ssi_matrices_->complete_sca_tra_manifold_structure_matrix();
    ssi_matrices_->complete_sca_tra_sca_tra_manifold_matrix();
    ssi_matrices_->complete_sca_tra_manifold_sca_tra_matrix();

    manifoldscatraflux_->complete_matrix_manifold_sca_tra();
    manifoldscatraflux_->complete_matrix_manifold_structure();
    manifoldscatraflux_->complete_matrix_sca_tra_manifold();
    manifoldscatraflux_->complete_matrix_sca_tra_structure();
    manifoldscatraflux_->complete_system_matrix_manifold();
    manifoldscatraflux_->complete_system_matrix_sca_tra();
  }
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
const Teuchos::RCP<const Epetra_Map>& SSI::SsiMono::dof_row_map() const
{
  return MapsSubProblems()->FullMap();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::SsiMono::setup_contact_strategy()
{
  // get the contact solution strategy
  auto contact_solution_type = CORE::UTILS::IntegralValue<INPAR::CONTACT::SolvingStrategy>(
      GLOBAL::Problem::Instance()->contact_dynamic_params(), "STRATEGY");

  if (contact_solution_type == INPAR::CONTACT::solution_nitsche)
  {
    if (CORE::UTILS::IntegralValue<INPAR::STR::IntegrationStrategy>(
            GLOBAL::Problem::Instance()->structural_dynamic_params(), "INT_STRATEGY") !=
        INPAR::STR::int_standard)
      FOUR_C_THROW("ssi contact only with new structural time integration");

    // get the contact model evaluator and store a pointer to the strategy
    auto& model_evaluator_contact = dynamic_cast<STR::MODELEVALUATOR::Contact&>(
        structure_field()->ModelEvaluator(INPAR::STR::model_contact));
    contact_strategy_nitsche_ = Teuchos::rcp_dynamic_cast<CONTACT::NitscheStrategySsi>(
        model_evaluator_contact.StrategyPtr(), true);
  }
  else
    FOUR_C_THROW("Only Nitsche contact implemented for SSI problems at the moment!");
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SsiMono::Init(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& scatraparams, const Teuchos::ParameterList& structparams,
    const std::string& struct_disname, const std::string& scatra_disname, bool isAle)
{
  // check input parameters for scalar transport field
  if (CORE::UTILS::IntegralValue<INPAR::SCATRA::VelocityField>(scatraparams, "VELOCITYFIELD") !=
      INPAR::SCATRA::velocity_Navier_Stokes)
    FOUR_C_THROW("Invalid type of velocity field for scalar-structure interaction!");

  if (CORE::UTILS::IntegralValue<INPAR::STR::DynamicType>(structparams, "DYNAMICTYP") ==
      INPAR::STR::DynamicType::dyna_statics)
    FOUR_C_THROW(
        "Mass conservation is not fulfilled if 'Statics' time integration is chosen since the "
        "deformation velocities are incorrectly calculated.\n"
        "Use 'NEGLECTINERTIA Yes' in combination with another time integration scheme instead!");

  // initialize strategy for Newton-Raphson convergence check
  switch (
      Teuchos::getIntegralValue<INPAR::SSI::ScaTraTimIntType>(globaltimeparams, "SCATRATIMINTTYPE"))
  {
    case INPAR::SSI::ScaTraTimIntType::elch:
    {
      if (is_sca_tra_manifold())
      {
        strategy_convcheck_ =
            Teuchos::rcp(new SSI::SsiMono::ConvCheckStrategyElchScaTraManifold(globaltimeparams));
      }
      else
        strategy_convcheck_ =
            Teuchos::rcp(new SSI::SsiMono::ConvCheckStrategyElch(globaltimeparams));
      break;
    }

    case INPAR::SSI::ScaTraTimIntType::standard:
    {
      strategy_convcheck_ = Teuchos::rcp(new SSI::SsiMono::ConvCheckStrategyStd(globaltimeparams));
      break;
    }

    default:
    {
      FOUR_C_THROW("Type of scalar transport time integrator currently not supported!");
      break;
    }
  }

  // call base class routine
  SSIBase::Init(
      comm, globaltimeparams, scatraparams, structparams, struct_disname, scatra_disname, isAle);
}


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SsiMono::output()
{
  // output scalar transport field
  ScaTraField()->check_and_write_output_and_restart();
  if (is_sca_tra_manifold())
  {
    // domain output
    ScaTraManifold()->check_and_write_output_and_restart();
    // coupling output
    if (manifoldscatraflux_->DoOutput()) manifoldscatraflux_->Output();
  }

  // output structure field
  structure_field()->Output();

  if (print_matlab_) print_system_matrix_rhs_to_mat_lab_format();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::SsiMono::read_restart(int restart)
{
  // call base class
  SSIBase::read_restart(restart);

  // do ssi contact specific tasks
  if (SSIInterfaceContact())
  {
    setup_contact_strategy();
    set_ssi_contact_states(ScaTraField()->Phinp());
  }
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SsiMono::prepare_time_loop()
{
  set_struct_solution(structure_field()->Dispnp(), structure_field()->Velnp(),
      is_s2_i_kinetics_with_pseudo_contact());

  // calculate initial potential field if needed
  if (do_calculate_initial_potential_field()) calc_initial_potential_field();

  // calculate initial time derivatives
  calc_initial_time_derivative();

  ScaTraField()->prepare_time_loop();
  if (is_sca_tra_manifold()) ScaTraManifold()->prepare_time_loop();
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SsiMono::prepare_time_step()
{
  // update time and time step
  increment_time_and_step();

  // pass structural degrees of freedom to scalar transport discretization
  set_struct_solution(structure_field()->Dispnp(), structure_field()->Velnp(),
      is_s2_i_kinetics_with_pseudo_contact());

  // prepare time step for scalar transport field
  ScaTraField()->prepare_time_step();
  if (is_sca_tra_manifold()) ScaTraManifold()->prepare_time_step();

  // if adaptive time stepping and different time step size: calculate time step in scatra
  // (prepare_time_step() of Scatra) and pass to other fields
  if (ScaTraField()->TimeStepAdapted()) set_dt_from_sca_tra_to_ssi();

  // pass scalar transport degrees of freedom to structural discretization
  // has to be called AFTER ScaTraField()->prepare_time_step() to ensure
  // consistent scalar transport state vector with valid Dirichlet conditions
  SetScatraSolution(ScaTraField()->Phinp());
  if (is_sca_tra_manifold()) set_scatra_manifold_solution(ScaTraManifold()->Phinp());

  // evaluate temperature from function and set to structural discretization
  evaluate_and_set_temperature_field();

  // prepare time step for structural field
  structure_field()->prepare_time_step();

  // structure_field()->prepare_time_step() evaluates the DBC displaements on the master side. Now,
  // the master side displacements are copied to slave side to consider non zero DBC values in the
  // first Newton step on the slave side in case of interface mesh tying
  if (ssi_interface_meshtying())
  {
    for (const auto& meshtying : ssi_structure_mesh_tying()->MeshTyingHandlers())
    {
      auto coupling_adapter = meshtying->SlaveMasterCoupling();
      auto coupling_map_extractor = meshtying->slave_master_extractor();

      // displacements
      coupling_map_extractor->InsertVector(
          coupling_adapter->MasterToSlave(
              coupling_map_extractor->ExtractVector(structure_field()->Dispnp(), 2)),
          1, structure_field()->WriteAccessDispnp());
      structure_field()->set_state(structure_field()->WriteAccessDispnp());
    }
  }

  print_time_step_info();
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SsiMono::Setup()
{
  // call base class routine
  SSIBase::Setup();

  // safety checks
  if (ScaTraField()->NumScal() != 1)
  {
    FOUR_C_THROW(
        "Since the ssi_monolithic framework is only implemented for usage in combination with "
        "volume change laws 'MAT_InelasticDefgradLinScalarIso' or "
        "'MAT_InelasticDefgradLinScalarAniso' so far and these laws are implemented for only "
        "one transported scalar at the moment it is not reasonable to use them with more than one "
        "transported scalar. So you need to cope with it or change implementation! ;-)");
  }
  const auto ssi_params = GLOBAL::Problem::Instance()->SSIControlParams();

  const bool calc_initial_pot_elch = CORE::UTILS::IntegralValue<bool>(
      GLOBAL::Problem::Instance()->ELCHControlParams(), "INITPOTCALC");
  const bool calc_initial_pot_ssi =
      CORE::UTILS::IntegralValue<bool>(ssi_params.sublist("ELCH"), "INITPOTCALC");

  if (ScaTraField()->EquilibrationMethod() != CORE::LINALG::EquilibrationMethod::none)
  {
    FOUR_C_THROW(
        "You are within the monolithic solid scatra interaction framework but activated a pure "
        "scatra equilibration method. Delete this from 'SCALAR TRANSPORT DYNAMIC' section and set "
        "it in 'SSI CONTROL/MONOLITHIC' instead.");
  }
  if (equilibration_method_.global != CORE::LINALG::EquilibrationMethod::local and
      (equilibration_method_.structure != CORE::LINALG::EquilibrationMethod::none or
          equilibration_method_.scatra != CORE::LINALG::EquilibrationMethod::none))
    FOUR_C_THROW("Either global equilibration or local equilibration");

  if (matrixtype_ == CORE::LINALG::MatrixType::sparse and
      (equilibration_method_.structure != CORE::LINALG::EquilibrationMethod::none or
          equilibration_method_.scatra != CORE::LINALG::EquilibrationMethod::none))
    FOUR_C_THROW("Block based equilibration only for block matrices");

  if (!CORE::UTILS::IntegralValue<int>(
          GLOBAL::Problem::Instance()->scalar_transport_dynamic_params(), "SKIPINITDER"))
  {
    FOUR_C_THROW(
        "Initial derivatives are already calculated in monolithic SSI. Enable 'SKIPINITDER' in the "
        "input file.");
  }

  if (calc_initial_pot_elch)
    FOUR_C_THROW("Initial potential is calculated by SSI. Disable in Elch section.");
  if (calc_initial_pot_ssi and Teuchos::getIntegralValue<INPAR::SSI::ScaTraTimIntType>(ssi_params,
                                   "SCATRATIMINTTYPE") != INPAR::SSI::ScaTraTimIntType::elch)
    FOUR_C_THROW("Calculation of initial potential only in case of Elch");

  if (!ScaTraField()->IsIncremental())
    FOUR_C_THROW(
        "Must have incremental solution approach for monolithic scalar-structure interaction!");

  if (ssi_interface_meshtying() and
      meshtying_strategy_s2_i()->CouplingType() != INPAR::S2I::coupling_matching_nodes)
  {
    FOUR_C_THROW(
        "Monolithic scalar-structure interaction only implemented for scatra-scatra "
        "interface coupling with matching interface nodes!");
  }

  if (SSIInterfaceContact() and !IsRestart()) setup_contact_strategy();
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SsiMono::SetupSystem()
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
        FOUR_C_THROW(
            "Global system matrix with block structure requires AMGnxn block preconditioner!");

      // feed AMGnxn block preconditioner with null space information for each block of global
      // block system matrix
      build_null_spaces();

      break;
    }

    case CORE::LINALG::MatrixType::sparse:
    {
      // safety check
      if (ScaTraField()->SystemMatrix() == Teuchos::null)
        FOUR_C_THROW("Incompatible matrix type associated with scalar transport field!");
      break;
    }

    default:
    {
      FOUR_C_THROW("Type of global system matrix for scalar-structure interaction not recognized!");
      break;
    }
  }

  // initialize sub blocks and system matrix
  ssi_matrices_ = Teuchos::rcp(new SSI::UTILS::SSIMatrices(
      ssi_maps_, matrixtype_, ScaTraField()->MatrixType(), is_sca_tra_manifold()));

  // initialize residual and increment vectors
  ssi_vectors_ = Teuchos::rcp(new SSI::UTILS::SSIVectors(ssi_maps_, is_sca_tra_manifold()));

  // initialize strategy for assembly
  strategy_assemble_ = SSI::BuildAssembleStrategy(
      ssi_maps_, is_sca_tra_manifold(), matrixtype_, ScaTraField()->MatrixType());

  if (is_sca_tra_manifold())
  {
    // initialize object, that performs evaluations of OD coupling
    scatrastructure_off_diagcoupling_ =
        Teuchos::rcp(new SSI::ScatraManifoldStructureOffDiagCoupling(BlockMapStructure(),
            SSIMaps()->StructureDofRowMap(), ssi_structure_mesh_tying(), meshtying_strategy_s2_i(),
            ScaTraField(), ScaTraManifold(), structure_field()));

    // initialize object, that performs evaluations of scatra - scatra on manifold coupling
    manifoldscatraflux_ = Teuchos::rcp(new SSI::ScaTraManifoldScaTraFluxEvaluator(*this));

    // initialize object, that performs meshtying between manifold domains
    strategy_manifold_meshtying_ =
        SSI::BuildManifoldMeshTyingStrategy(ScaTraManifold()->discretization(), ssi_maps_,
            is_sca_tra_manifold_meshtying(), ScaTraManifold()->MatrixType());
  }
  else
  {
    scatrastructure_off_diagcoupling_ = Teuchos::rcp(new SSI::ScatraStructureOffDiagCoupling(
        BlockMapStructure(), SSIMaps()->StructureDofRowMap(), ssi_structure_mesh_tying(),
        meshtying_strategy_s2_i(), ScaTraField(), structure_field()));
  }
  // instantiate appropriate equilibration class
  strategy_equilibration_ = CORE::LINALG::BuildEquilibration(
      matrixtype_, get_block_equilibration(), MapsSubProblems()->FullMap());

  // instantiate appropriate contact class
  strategy_contact_ =
      SSI::BuildContactStrategy(nitsche_strategy_ssi(), ssi_maps_, ScaTraField()->MatrixType());

  // instantiate appropriate mesh tying class
  strategy_meshtying_ = SSI::BuildMeshtyingStrategy(
      is_sca_tra_manifold(), ScaTraField()->MatrixType(), ssi_maps_, ssi_structure_mesh_tying());

  // instantiate Dirichlet boundary condition handler class
  dbc_handler_ = SSI::BuildDBCHandler(is_sca_tra_manifold(), matrixtype_, ScaTraField(),
      is_sca_tra_manifold() ? ScaTraManifold() : Teuchos::null, ssi_maps_, structure_field());
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void SSI::SsiMono::SetScatraSolution(Teuchos::RCP<const Epetra_Vector> phi) const
{
  // call base class
  SSIBase::SetScatraSolution(phi);

  // set state for contact evaluation
  if (contact_strategy_nitsche_ != Teuchos::null) set_ssi_contact_states(phi);
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void SSI::SsiMono::set_ssi_contact_states(Teuchos::RCP<const Epetra_Vector> phi) const
{
  contact_strategy_nitsche_->set_state(MORTAR::state_scalar, *phi);
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void SSI::SsiMono::solve_linear_system()
{
  TEUCHOS_FUNC_TIME_MONITOR("SSI mono: solve linear system");
  strategy_equilibration_->EquilibrateSystem(
      ssi_matrices_->SystemMatrix(), ssi_vectors_->Residual(), block_map_system_matrix());

  // solve global system of equations
  // Dirichlet boundary conditions have already been applied to global system of equations
  CORE::LINALG::SolverParams solver_params;
  solver_params.refactor = true;
  solver_params.reset = IterationCount() == 1;
  if (relax_lin_solver_iter_step_ > 0)
  {
    solver_->ResetTolerance();
    if (IterationCount() <= relax_lin_solver_iter_step_)
    {
      solver_params.tolerance = solver_->GetTolerance() * relax_lin_solver_tolerance_;
    }
  }
  solver_->Solve(ssi_matrices_->SystemMatrix()->EpetraOperator(), ssi_vectors_->Increment(),
      ssi_vectors_->Residual(), solver_params);

  strategy_equilibration_->unequilibrate_increment(ssi_vectors_->Increment());
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SsiMono::newton_loop()
{
  TEUCHOS_FUNC_TIME_MONITOR("SSI mono: solve Newton loop");
  // reset counter for Newton-Raphson iteration
  ResetIterationCount();

  // start Newton-Raphson iteration
  while (true)
  {
    // update iteration counter
    increment_iteration_count();

    timer_->reset();

    // store time before evaluating elements and assembling global system of equations
    const double time_before_evaluate = timer_->wallTime();

    // set solution from last Newton step to all fields
    distribute_solution_all_fields();

    // evaluate sub problems and get all matrices and right-hand-sides
    evaluate_subproblems();

    // complete the sub problem matrices
    complete_subproblem_matrices();

    // assemble global system of equations
    assemble_mat_and_rhs();

    // apply the Dirichlet boundary conditions to global system
    apply_dbc_to_system();

    // time needed for evaluating elements and assembling global system of equations
    double my_evaluation_time = timer_->wallTime() - time_before_evaluate;
    Comm().MaxAll(&my_evaluation_time, &dt_eval_, 1);

    // safety check
    if (!ssi_matrices_->SystemMatrix()->Filled())
      FOUR_C_THROW("Complete() has not been called on global system matrix yet!");

    // check termination criterion for Newton-Raphson iteration
    if (strategy_convcheck_->exit_newton_raphson(*this)) break;

    // clear the global increment vector
    ssi_vectors_->ClearIncrement();

    // store time before solving global system of equations
    const double time_before_solving = timer_->wallTime();

    solve_linear_system();

    // time needed for solving global system of equations
    double my_solve_time = timer_->wallTime() - time_before_solving;
    Comm().MaxAll(&my_solve_time, &dt_solve_, 1);

    // output performance statistics associated with linear solver into text file if
    // applicable
    if (CORE::UTILS::IntegralValue<bool>(
            *ScaTraField()->ScatraParameterList(), "OUTPUTLINSOLVERSTATS"))
      ScaTraField()->output_lin_solver_stats(*solver_, dt_solve_, Step(), IterationCount(),
          ssi_vectors_->Residual()->Map().NumGlobalElements());

    // update states for next Newton iteration
    update_iter_sca_tra();
    update_iter_structure();
  }
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SsiMono::Timeloop()
{
  if (Step() == 0) prepare_time_loop();

  // time loop
  while (NotFinished() and ScaTraField()->NotFinished())
  {
    TEUCHOS_FUNC_TIME_MONITOR("SSI mono: solve time step");
    // prepare time step
    prepare_time_step();

    // store time before calling nonlinear solver
    const double time = timer_->wallTime();

    // evaluate time step
    newton_loop();

    // determine time spent by nonlinear solver and take maximum over all processors via
    // communication
    double mydtnonlinsolve(timer_->wallTime() - time), dtnonlinsolve(0.);
    Comm().MaxAll(&mydtnonlinsolve, &dtnonlinsolve, 1);

    // output performance statistics associated with nonlinear solver into *.csv file if
    // applicable
    if (CORE::UTILS::IntegralValue<int>(
            *ScaTraField()->ScatraParameterList(), "OUTPUTNONLINSOLVERSTATS"))
      ScaTraField()->output_nonlin_solver_stats(IterationCount(), dtnonlinsolve, Step(), Comm());

    prepare_output();

    // update scalar transport and structure fields
    update();

    // output solution to screen and files
    output();
  }
  strategy_convcheck_->print_non_converged_steps(Comm().MyPID());
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void SSI::SsiMono::update()
{
  // update scalar transport field
  ScaTraField()->Update();
  if (is_sca_tra_manifold()) ScaTraManifold()->Update();

  // update structure field
  structure_field()->Update();
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void SSI::SsiMono::update_iter_sca_tra()
{
  // update scalar transport field
  ScaTraField()->UpdateIter(MapsSubProblems()->ExtractVector(
      ssi_vectors_->Increment(), UTILS::SSIMaps::GetProblemPosition(Subproblem::scalar_transport)));
  ScaTraField()->compute_intermediate_values();

  if (is_sca_tra_manifold())
  {
    auto increment_manifold = MapsSubProblems()->ExtractVector(
        ssi_vectors_->Increment(), UTILS::SSIMaps::GetProblemPosition(Subproblem::manifold));

    // reconstruct slave side solution from master side
    if (is_sca_tra_manifold_meshtying())
    {
      for (const auto& meshtying :
          strategy_manifold_meshtying_->SSIMeshTying()->MeshTyingHandlers())
      {
        auto coupling_adapter = meshtying->SlaveMasterCoupling();
        auto multimap = meshtying->slave_master_extractor();

        auto master_dofs = multimap->ExtractVector(increment_manifold, 2);
        auto master_dofs_to_slave = coupling_adapter->MasterToSlave(master_dofs);
        multimap->InsertVector(master_dofs_to_slave, 1, increment_manifold);
      }
    }

    ScaTraManifold()->UpdateIter(increment_manifold);
    ScaTraManifold()->compute_intermediate_values();
  }
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void SSI::SsiMono::update_iter_structure()
{
  // set up structural increment vector
  const Teuchos::RCP<Epetra_Vector> increment_structure = MapsSubProblems()->ExtractVector(
      ssi_vectors_->Increment(), UTILS::SSIMaps::GetProblemPosition(Subproblem::structure));

  // consider structural meshtying. Copy master increments and displacements to slave side.
  if (ssi_interface_meshtying())
  {
    for (const auto& meshtying : ssi_structure_mesh_tying()->MeshTyingHandlers())
    {
      auto coupling_adapter = meshtying->SlaveMasterCoupling();
      auto coupling_map_extractor = meshtying->slave_master_extractor();

      // displacements
      coupling_map_extractor->InsertVector(
          coupling_adapter->MasterToSlave(
              coupling_map_extractor->ExtractVector(structure_field()->Dispnp(), 2)),
          1, structure_field()->WriteAccessDispnp());
      structure_field()->set_state(structure_field()->WriteAccessDispnp());

      // increment
      coupling_map_extractor->InsertVector(
          coupling_adapter->MasterToSlave(
              coupling_map_extractor->ExtractVector(increment_structure, 2)),
          1, increment_structure);
    }
  }

  // update displacement of structure field
  structure_field()->update_state_incrementally(increment_structure);
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
std::vector<CORE::LINALG::EquilibrationMethod> SSI::SsiMono::get_block_equilibration()
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
        if (is_sca_tra_manifold())
        {
          auto block_positions_scatra_manifold = ssi_maps_->GetBlockPositions(Subproblem::manifold);

          equilibration_method_vector = std::vector<CORE::LINALG::EquilibrationMethod>(
              block_positions_scatra.size() + block_position_structure.size() +
                  block_positions_scatra_manifold.size(),
              CORE::LINALG::EquilibrationMethod::none);
        }
        else
        {
          equilibration_method_vector = std::vector<CORE::LINALG::EquilibrationMethod>(
              block_positions_scatra.size() + block_position_structure.size(),
              CORE::LINALG::EquilibrationMethod::none);
        }


        for (const int block_position_scatra : block_positions_scatra)
          equilibration_method_vector.at(block_position_scatra) = equilibration_method_.scatra;

        equilibration_method_vector.at(block_position_structure.at(0)) =
            equilibration_method_.structure;

        if (is_sca_tra_manifold())
        {
          for (const int block_position_scatra_manifold :
              ssi_maps_->GetBlockPositions(Subproblem::manifold))
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
      FOUR_C_THROW("Invalid matrix type associated with system matrix field!");
      break;
    }
  }
  return equilibration_method_vector;
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void SSI::SsiMono::evaluate_sca_tra()
{
  // evaluate the scatra field
  ScaTraField()->PrepareLinearSolve();

  // copy the matrix to the corresponding ssi matrix and complete it such that additional
  // contributions like contact contributions can be added before assembly
  ssi_matrices_->ScaTraMatrix()->Add(*ScaTraField()->system_matrix_operator(), false, 1.0, 1.0);

  // copy the residual to the corresponding ssi vector to enable application of contact
  // contributions before assembly
  ssi_vectors_->ScatraResidual()->Update(1.0, *ScaTraField()->Residual(), 1.0);
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void SSI::SsiMono::evaluate_sca_tra_manifold()
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
void SSI::SsiMono::prepare_output()
{
  constexpr bool force_prepare = false;
  structure_field()->prepare_output(force_prepare);

  // prepare output of coupling sctra manifold - scatra
  if (is_sca_tra_manifold() and manifoldscatraflux_->DoOutput())
  {
    distribute_solution_all_fields();
    manifoldscatraflux_->evaluate_sca_tra_manifold_inflow();
  }
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void SSI::SsiMono::distribute_solution_all_fields(const bool restore_velocity)
{
  // has to be called before the call of 'set_struct_solution()' to have updated stress/strain
  // states
  if (is_s2_i_kinetics_with_pseudo_contact()) structure_field()->determine_stress_strain();

  // clear all states before redistributing the new states
  structure_field()->discretization()->ClearState(true);
  ScaTraField()->discretization()->ClearState(true);
  if (is_sca_tra_manifold()) ScaTraManifold()->discretization()->ClearState(true);

  // needed to communicate to NOX state
  if (restore_velocity)
  {
    auto vel_temp = *structure_field()->Velnp();
    structure_field()->set_state(structure_field()->WriteAccessDispnp());
    structure_field()->WriteAccessVelnp()->Update(1.0, vel_temp, 0.0);
  }
  else
    structure_field()->set_state(structure_field()->WriteAccessDispnp());

  // distribute states to other fields
  set_struct_solution(structure_field()->Dispnp(), structure_field()->Velnp(),
      is_s2_i_kinetics_with_pseudo_contact());
  SetScatraSolution(ScaTraField()->Phinp());
  if (is_sca_tra_manifold()) set_scatra_manifold_solution(ScaTraManifold()->Phinp());
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void SSI::SsiMono::calc_initial_potential_field()
{
  const auto equpot = CORE::UTILS::IntegralValue<INPAR::ELCH::EquPot>(
      GLOBAL::Problem::Instance()->ELCHControlParams(), "EQUPOT");
  if (equpot != INPAR::ELCH::equpot_divi and equpot != INPAR::ELCH::equpot_enc_pde and
      equpot != INPAR::ELCH::equpot_enc_pde_elim)
  {
    FOUR_C_THROW(
        "Initial potential field cannot be computed for chosen closing equation for electric "
        "potential!");
  }

  // store initial velocity to restore them afterwards
  auto init_velocity = *structure_field()->Velnp();

  // cast scatra time integrators to elch to call elch specific methods
  auto scatra_elch = Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntElch>(ScaTraField());
  auto manifold_elch = is_sca_tra_manifold()
                           ? Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntElch>(ScaTraManifold())
                           : Teuchos::null;
  if (scatra_elch == Teuchos::null or (is_sca_tra_manifold() and manifold_elch == Teuchos::null))
    FOUR_C_THROW("Cast to Elch time integrator faild. Scatra is not an Elch problem");

  // prepare specific time integrators
  scatra_elch->pre_calc_initial_potential_field();
  if (is_sca_tra_manifold()) manifold_elch->pre_calc_initial_potential_field();

  auto scatra_elch_splitter = ScaTraField()->Splitter();
  auto manifold_elch_splitter =
      is_sca_tra_manifold() ? ScaTraManifold()->Splitter() : Teuchos::null;

  ResetIterationCount();

  while (true)
  {
    increment_iteration_count();

    timer_->reset();

    // store time before evaluating elements and assembling global system of equations
    const double time_before_evaluate = timer_->wallTime();

    // prepare full SSI system
    distribute_solution_all_fields(true);
    evaluate_subproblems();

    // complete the sub problem matrices
    complete_subproblem_matrices();

    assemble_mat_and_rhs();
    apply_dbc_to_system();

    // apply artificial Dirichlet boundary conditions to system of equations (on concentration
    // dofs and on structure dofs)
    Teuchos::RCP<Epetra_Map> pseudo_dbc_map;
    if (is_sca_tra_manifold())
    {
      auto conc_map = CORE::LINALG::MergeMap(
          scatra_elch_splitter->OtherMap(), manifold_elch_splitter->OtherMap());
      pseudo_dbc_map = CORE::LINALG::MergeMap(conc_map, structure_field()->dof_row_map());
    }
    else
    {
      pseudo_dbc_map = CORE::LINALG::MergeMap(
          scatra_elch_splitter->OtherMap(), structure_field()->dof_row_map());
    }

    auto dbc_zeros = Teuchos::rcp(new Epetra_Vector(*pseudo_dbc_map, true));

    auto rhs = ssi_vectors_->Residual();
    CORE::LINALG::apply_dirichlet_to_system(*ssi_matrices_->SystemMatrix(),
        *ssi_vectors_->Increment(), *rhs, *dbc_zeros, *pseudo_dbc_map);
    ssi_vectors_->Residual()->Update(1.0, *rhs, 0.0);

    // time needed for evaluating elements and assembling global system of equations
    double my_evaluation_time = timer_->wallTime() - time_before_evaluate;
    Comm().MaxAll(&my_evaluation_time, &dt_eval_, 1);

    if (strategy_convcheck_->exit_newton_raphson_init_pot_calc(*this)) break;

    // solve for potential increments
    ssi_vectors_->ClearIncrement();

    // store time before solving global system of equations
    const double time_before_solving = timer_->wallTime();

    solve_linear_system();

    // time needed for solving global system of equations
    double my_solve_time = timer_->wallTime() - time_before_solving;
    Comm().MaxAll(&my_solve_time, &dt_solve_, 1);

    // update potential dofs in scatra and manifold fields
    update_iter_sca_tra();

    // copy initial state vector
    ScaTraField()->Phin()->Update(1.0, *ScaTraField()->Phinp(), 0.0);
    if (is_sca_tra_manifold())
      ScaTraManifold()->Phin()->Update(1.0, *ScaTraManifold()->Phinp(), 0.0);

    // update state vectors for intermediate time steps (only for generalized alpha)
    ScaTraField()->compute_intermediate_values();
    if (is_sca_tra_manifold()) ScaTraManifold()->compute_intermediate_values();
  }

  scatra_elch->post_calc_initial_potential_field();
  if (is_sca_tra_manifold()) manifold_elch->post_calc_initial_potential_field();

  structure_field()->WriteAccessVelnp()->Update(1.0, init_velocity, 0.0);
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void SSI::SsiMono::calc_initial_time_derivative()
{
  // store initial velocity to restore them afterwards
  auto init_velocity = *structure_field()->Velnp();

  const bool is_elch = is_elch_sca_tra_tim_int_type();

  // prepare specific time integrators
  ScaTraField()->pre_calc_initial_time_derivative();
  if (is_sca_tra_manifold()) ScaTraManifold()->pre_calc_initial_time_derivative();

  auto scatra_elch_splitter = is_elch ? ScaTraField()->Splitter() : Teuchos::null;
  auto manifold_elch_splitter =
      (is_elch and is_sca_tra_manifold()) ? ScaTraManifold()->Splitter() : Teuchos::null;

  // initial screen output
  if (Comm().MyPID() == 0)
  {
    std::cout << "Calculating initial time derivative of state variables on discretization "
              << ScaTraField()->discretization()->Name();
    if (is_sca_tra_manifold())
      std::cout << " and discretization " << ScaTraManifold()->discretization()->Name();
    std::cout << std::endl;
  }

  // evaluate Dirichlet and Neumann boundary conditions
  ScaTraField()->ApplyBCToSystem();
  if (is_sca_tra_manifold()) ScaTraManifold()->ApplyBCToSystem();

  // clear history values (this is the first step)
  ScaTraField()->Hist()->PutScalar(0.0);
  if (is_sca_tra_manifold()) ScaTraManifold()->Hist()->PutScalar(0.0);

  // In a first step, we assemble the standard global system of equations (we need the residual)
  distribute_solution_all_fields(true);
  evaluate_subproblems();

  // complete the sub problem matrices
  complete_subproblem_matrices();

  assemble_mat_and_rhs();
  apply_dbc_to_system();

  // prepare mass matrices of sub problems and global system
  auto massmatrix_scatra =
      ScaTraField()->MatrixType() == CORE::LINALG::MatrixType::sparse
          ? Teuchos::rcp_dynamic_cast<CORE::LINALG::SparseOperator>(
                UTILS::SSIMatrices::setup_sparse_matrix(ScaTraField()->dof_row_map()))
          : Teuchos::rcp_dynamic_cast<CORE::LINALG::SparseOperator>(
                UTILS::SSIMatrices::setup_block_matrix(
                    ScaTraField()->BlockMaps(), ScaTraField()->BlockMaps()));

  auto massmatrix_manifold =
      is_sca_tra_manifold()
          ? (ScaTraManifold()->MatrixType() == CORE::LINALG::MatrixType::sparse
                    ? Teuchos::rcp_dynamic_cast<CORE::LINALG::SparseOperator>(
                          UTILS::SSIMatrices::setup_sparse_matrix(ScaTraManifold()->dof_row_map()))
                    : Teuchos::rcp_dynamic_cast<CORE::LINALG::SparseOperator>(
                          UTILS::SSIMatrices::setup_block_matrix(
                              ScaTraManifold()->BlockMaps(), ScaTraManifold()->BlockMaps())))
          : Teuchos::null;

  auto massmatrix_system = MatrixType() == CORE::LINALG::MatrixType::sparse
                               ? Teuchos::rcp_dynamic_cast<CORE::LINALG::SparseOperator>(
                                     UTILS::SSIMatrices::setup_sparse_matrix(dof_row_map()))
                               : Teuchos::rcp_dynamic_cast<CORE::LINALG::SparseOperator>(
                                     UTILS::SSIMatrices::setup_block_matrix(
                                         block_map_system_matrix(), block_map_system_matrix()));

  // fill ones on main diag of structure block (not solved)
  auto ones_struct = Teuchos::rcp(new Epetra_Vector(*structure_field()->dof_row_map(), true));
  ones_struct->PutScalar(1.0);
  MatrixType() == CORE::LINALG::MatrixType::sparse
      ? CORE::LINALG::InsertMyRowDiagonalIntoUnfilledMatrix(
            *CORE::LINALG::CastToSparseMatrixAndCheckSuccess(massmatrix_system), *ones_struct)
      : CORE::LINALG::InsertMyRowDiagonalIntoUnfilledMatrix(
            CORE::LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(massmatrix_system)
                ->Matrix(ssi_maps_->GetBlockPositions(Subproblem::structure).at(0),
                    ssi_maps_->GetBlockPositions(Subproblem::structure).at(0)),
            *ones_struct);

  // extract residuals of scatra and manifold from global residual
  auto rhs_scatra = Teuchos::rcp(new Epetra_Vector(*ScaTraField()->dof_row_map(), true));
  auto rhs_manifold = is_sca_tra_manifold()
                          ? Teuchos::rcp(new Epetra_Vector(*ScaTraManifold()->dof_row_map(), true))
                          : Teuchos::null;

  rhs_scatra->Update(1.0,
      *MapsSubProblems()->ExtractVector(ssi_vectors_->Residual(),
          UTILS::SSIMaps::GetProblemPosition(Subproblem::scalar_transport)),
      0.0);
  if (is_sca_tra_manifold())
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
  ScaTraField()->evaluate_initial_time_derivative(massmatrix_scatra, rhs_scatra);
  if (is_sca_tra_manifold())
    ScaTraManifold()->evaluate_initial_time_derivative(massmatrix_manifold, rhs_manifold);

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

          for (int i = 0; i < static_cast<int>(positions_scatra.size()); ++i)
          {
            const int position_scatra = positions_scatra.at(i);
            massmatrix_system_block->Matrix(position_scatra, position_scatra)
                .Add(massmatrix_scatra_block->Matrix(i, i), false, 1.0, 1.0);
          }
          if (is_sca_tra_manifold())
          {
            auto positions_manifold = ssi_maps_->GetBlockPositions(Subproblem::manifold);

            auto massmatrix_manifold_block =
                CORE::LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(massmatrix_manifold);

            for (int i = 0; i < static_cast<int>(positions_manifold.size()); ++i)
            {
              const int position_manifold = positions_manifold.at(i);
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
              ssi_maps_->GetBlockPositions(Subproblem::scalar_transport).at(0);

          massmatrix_system_block->Matrix(position_scatra, position_scatra)
              .Add(*CORE::LINALG::CastToSparseMatrixAndCheckSuccess(massmatrix_scatra), false, 1.0,
                  1.0);

          if (is_sca_tra_manifold())
          {
            const int position_manifold = ssi_maps_->GetBlockPositions(Subproblem::manifold).at(0);

            massmatrix_system_block->Matrix(position_manifold, position_manifold)
                .Add(*CORE::LINALG::CastToSparseMatrixAndCheckSuccess(massmatrix_manifold), false,
                    1.0, 1.0);
          }
          break;
        }

        default:
        {
          FOUR_C_THROW("Invalid matrix type associated with scalar transport field!");
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

      if (is_sca_tra_manifold())
        massmatrix_system_sparse->Add(
            *CORE::LINALG::CastToSparseMatrixAndCheckSuccess(massmatrix_manifold), false, 1.0, 1.0);

      massmatrix_system->Complete(*dof_row_map(), *dof_row_map());
      break;
    }
    default:
    {
      FOUR_C_THROW("Type of global system matrix for scalar-structure interaction not recognized!");
      break;
    }
  }

  // reconstruct global residual from partial residuals
  auto rhs_system = Teuchos::RCP<Epetra_Vector>(new Epetra_Vector(*dof_row_map(), true));
  MapsSubProblems()->InsertVector(
      rhs_scatra, UTILS::SSIMaps::GetProblemPosition(Subproblem::scalar_transport), rhs_system);
  if (is_sca_tra_manifold())
    MapsSubProblems()->InsertVector(
        rhs_manifold, UTILS::SSIMaps::GetProblemPosition(Subproblem::manifold), rhs_system);

  // apply artificial Dirichlet boundary conditions to system of equations to non-transported
  // scalars and structure
  Teuchos::RCP<Epetra_Map> pseudo_dbc_map;
  if (is_sca_tra_manifold() and is_elch)
  {
    auto conc_map =
        CORE::LINALG::MergeMap(scatra_elch_splitter->CondMap(), manifold_elch_splitter->CondMap());
    pseudo_dbc_map = CORE::LINALG::MergeMap(conc_map, structure_field()->dof_row_map());
  }
  else if (is_elch)
  {
    pseudo_dbc_map =
        CORE::LINALG::MergeMap(scatra_elch_splitter->CondMap(), structure_field()->dof_row_map());
  }
  else
    pseudo_dbc_map = Teuchos::rcp(new Epetra_Map(*structure_field()->dof_row_map()));

  auto dbc_zeros = Teuchos::rcp(new Epetra_Vector(*pseudo_dbc_map, true));

  // temporal derivative of transported scalars
  auto phidtnp_system = Teuchos::RCP<Epetra_Vector>(new Epetra_Vector(*dof_row_map(), true));
  CORE::LINALG::apply_dirichlet_to_system(
      *massmatrix_system, *phidtnp_system, *rhs_system, *dbc_zeros, *(pseudo_dbc_map));

  // solve global system of equations for initial time derivative of state variables
  CORE::LINALG::SolverParams solver_params;
  solver_params.refactor = true;
  solver_params.reset = true;
  solver_->Solve(massmatrix_system->EpetraOperator(), phidtnp_system, rhs_system, solver_params);

  // copy solution to sub problmes
  auto phidtnp_scatra = MapsSubProblems()->ExtractVector(
      phidtnp_system, UTILS::SSIMaps::GetProblemPosition(Subproblem::scalar_transport));
  ScaTraField()->Phidtnp()->Update(1.0, *phidtnp_scatra, 0.0);
  ScaTraField()->Phidtn()->Update(1.0, *phidtnp_scatra, 0.0);

  if (is_sca_tra_manifold())
  {
    auto phidtnp_manifold = MapsSubProblems()->ExtractVector(
        phidtnp_system, UTILS::SSIMaps::GetProblemPosition(Subproblem::manifold));
    ScaTraManifold()->Phidtnp()->Update(1.0, *phidtnp_manifold, 0.0);
    ScaTraManifold()->Phidtn()->Update(1.0, *phidtnp_manifold, 0.0);
  }

  // reset solver
  solver_->Reset();

  ScaTraField()->post_calc_initial_time_derivative();
  if (is_sca_tra_manifold()) ScaTraManifold()->post_calc_initial_time_derivative();

  structure_field()->WriteAccessVelnp()->Update(1.0, init_velocity, 0.0);
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> SSI::SsiMono::MapsSubProblems() const
{
  return ssi_maps_->MapsSubProblems();
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> SSI::SsiMono::BlockMapScaTra() const
{
  return ssi_maps_->BlockMapScaTra();
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> SSI::SsiMono::block_map_sca_tra_manifold() const
{
  return ssi_maps_->block_map_sca_tra_manifold();
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> SSI::SsiMono::BlockMapStructure() const
{
  return ssi_maps_->BlockMapStructure();
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> SSI::SsiMono::block_map_system_matrix() const
{
  return ssi_maps_->block_map_system_matrix();
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void SSI::SsiMono::print_time_step_info()
{
  if (Comm().MyPID() == 0)
  {
    std::cout << std::endl
              << "TIME: " << std::setw(11) << std::setprecision(4) << std::scientific << Time()
              << "/" << max_time() << "  DT = " << Dt() << "  STEP = " << Step() << "/" << n_step()
              << std::endl;
  }
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void SSI::SsiMono::print_system_matrix_rhs_to_mat_lab_format()
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
          filename << GLOBAL::Problem::Instance()->OutputControlFile()->FileName()
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

      const std::string filename = GLOBAL::Problem::Instance()->OutputControlFile()->FileName() +
                                   "_sparse_system_matrix.csv";

      CORE::LINALG::PrintMatrixInMatlabFormat(filename, *sparse_matrix->EpetraMatrix(), true);
      break;
    }

    default:
    {
      FOUR_C_THROW("Type of global system matrix for scalar-structure interaction not recognized!");
      break;
    }
  }

  // print rhs
  {
    const std::string filename =
        GLOBAL::Problem::Instance()->OutputControlFile()->FileName() + "_system_vector.csv";
    CORE::LINALG::PrintVectorInMatlabFormat(filename, *ssi_vectors_->Residual(), true);
  }

  // print full map
  {
    const std::string filename =
        GLOBAL::Problem::Instance()->OutputControlFile()->FileName() + "_full_map.csv";
    CORE::LINALG::PrintMapInMatlabFormat(filename, *ssi_maps_->MapSystemMatrix(), true);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SsiMono::set_scatra_manifold_solution(Teuchos::RCP<const Epetra_Vector> phi)
{
  // scatra values on master side copied to manifold
  auto manifold_on_scatra =
      CORE::LINALG::CreateVector(*ScaTraField()->discretization()->dof_row_map(), true);

  for (const auto& coup : manifoldscatraflux_->sca_tra_manifold_couplings())
  {
    auto manifold_cond = coup->manifold_map_extractor()->ExtractCondVector(*phi);
    auto manifold_on_scatra_cond = coup->CouplingAdapter()->SlaveToMaster(manifold_cond);
    coup->ScaTraMapExtractor()->InsertCondVector(manifold_on_scatra_cond, manifold_on_scatra);
  }
  ScaTraField()->discretization()->set_state(0, "manifold_on_scatra", manifold_on_scatra);
}
FOUR_C_NAMESPACE_CLOSE
