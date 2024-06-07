/*----------------------------------------------------------------------*/
/*! \file

\brief Standard solution strategy for standard scalar transport problems (without meshtying)

\level 2


*----------------------------------------------------------------------*/
#include "4C_scatra_timint_meshtying_strategy_std.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_scatra_timint_implicit.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | constructor                                               fang 12/14 |
 *----------------------------------------------------------------------*/
ScaTra::MeshtyingStrategyStd::MeshtyingStrategyStd(ScaTra::ScaTraTimIntImpl* scatratimint)
    : MeshtyingStrategyBase(scatratimint)
{
  return;
}  // ScaTra::MeshtyingStrategyStd::MeshtyingStrategyStd


/*-----------------------------------------------------------------------*
 | return global map of degrees of freedom                    fang 02/18 |
 *-----------------------------------------------------------------------*/
const Epetra_Map& ScaTra::MeshtyingStrategyStd::dof_row_map() const
{
  return *scatratimint_->dof_row_map();
}


/*----------------------------------------------------------------------*
 | dummy meshtying evaluate for standard scalar transport    fang 12/14 |
 *----------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyStd::EvaluateMeshtying()
{
  return;
}  // ScaTra::MeshtyingStrategyStd::EvaluateMeshtying


/*----------------------------------------------------------------------*
 | setup meshtying objects                                   fang 02/16 |
 *----------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyStd::setup_meshtying() { return; }


/*----------------------------------------------------------------------*
 | init meshtying objects                                   rauch 09/16 |
 *----------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyStd::InitMeshtying()
{
  // instantiate strategy for Newton-Raphson convergence check
  init_conv_check_strategy();
  return;
}

/*-----------------------------------------------------------------------------*
 | solve linear system of equations for standard scalar transport   fang 12/14 |
 *-----------------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyStd::Solve(
    const Teuchos::RCP<Core::LinAlg::Solver>& solver,                //!< solver
    const Teuchos::RCP<Core::LinAlg::SparseOperator>& systemmatrix,  //!< system matrix
    const Teuchos::RCP<Epetra_Vector>& increment,                    //!< increment vector
    const Teuchos::RCP<Epetra_Vector>& residual,                     //!< residual vector
    const Teuchos::RCP<Epetra_Vector>& phinp,                        //!< state vector at time n+1
    const int iteration,  //!< number of current Newton-Raphson iteration
    Core::LinAlg::SolverParams& solver_params) const
{
  solver_params.refactor = true;
  solver_params.reset = iteration == 1;
  solver->Solve(systemmatrix->EpetraOperator(), increment, residual, solver_params);

  return;
}  // ScaTra::MeshtyingStrategyStd::Solve


/*-------------------------------------------------------------------------*
 | return linear solver for global system of linear equations   fang 01/18 |
 *-------------------------------------------------------------------------*/
const Core::LinAlg::Solver& ScaTra::MeshtyingStrategyStd::Solver() const
{
  if (scatratimint_->Solver() == Teuchos::null) FOUR_C_THROW("Invalid linear solver!");
  return *scatratimint_->Solver();
}  // ScaTra::MeshtyingStrategyStd::Solver()


/*------------------------------------------------------------------------*
 | instantiate strategy for Newton-Raphson convergence check   fang 02/16 |
 *------------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyStd::init_conv_check_strategy()
{
  if (scatratimint_->MicroScale())
    convcheckstrategy_ = Teuchos::rcp(new ScaTra::ConvCheckStrategyStdMicroScale(
        scatratimint_->ScatraParameterList()->sublist("NONLINEAR")));
  else if (Global::Problem::Instance()->GetProblemType() == Core::ProblemType::poromultiphasescatra)
    convcheckstrategy_ = Teuchos::rcp(new ScaTra::ConvCheckStrategyPoroMultiphaseScatra(
        scatratimint_->ScatraParameterList()->sublist("NONLINEAR")));
  else
    convcheckstrategy_ = Teuchos::rcp(new ScaTra::ConvCheckStrategyStd(
        scatratimint_->ScatraParameterList()->sublist("NONLINEAR")));

  return;
}  // ScaTra::MeshtyingStrategyStd::init_conv_check_strategy

FOUR_C_NAMESPACE_CLOSE
