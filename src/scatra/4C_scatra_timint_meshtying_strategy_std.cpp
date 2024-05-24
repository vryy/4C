/*----------------------------------------------------------------------*/
/*! \file

\brief Standard solution strategy for standard scalar transport problems (without meshtying)

\level 2


*----------------------------------------------------------------------*/
#include "4C_scatra_timint_meshtying_strategy_std.hpp"

#include "4C_global_data.hpp"
#include "4C_lib_discret.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_scatra_timint_implicit.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | constructor                                               fang 12/14 |
 *----------------------------------------------------------------------*/
SCATRA::MeshtyingStrategyStd::MeshtyingStrategyStd(SCATRA::ScaTraTimIntImpl* scatratimint)
    : MeshtyingStrategyBase(scatratimint)
{
  return;
}  // SCATRA::MeshtyingStrategyStd::MeshtyingStrategyStd


/*-----------------------------------------------------------------------*
 | return global map of degrees of freedom                    fang 02/18 |
 *-----------------------------------------------------------------------*/
const Epetra_Map& SCATRA::MeshtyingStrategyStd::dof_row_map() const
{
  return *scatratimint_->dof_row_map();
}


/*----------------------------------------------------------------------*
 | dummy meshtying evaluate for standard scalar transport    fang 12/14 |
 *----------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyStd::EvaluateMeshtying()
{
  return;
}  // SCATRA::MeshtyingStrategyStd::EvaluateMeshtying


/*----------------------------------------------------------------------*
 | setup meshtying objects                                   fang 02/16 |
 *----------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyStd::SetupMeshtying() { return; }


/*----------------------------------------------------------------------*
 | init meshtying objects                                   rauch 09/16 |
 *----------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyStd::InitMeshtying()
{
  // instantiate strategy for Newton-Raphson convergence check
  init_conv_check_strategy();
  return;
}

/*-----------------------------------------------------------------------------*
 | solve linear system of equations for standard scalar transport   fang 12/14 |
 *-----------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyStd::Solve(
    const Teuchos::RCP<CORE::LINALG::Solver>& solver,                //!< solver
    const Teuchos::RCP<CORE::LINALG::SparseOperator>& systemmatrix,  //!< system matrix
    const Teuchos::RCP<Epetra_Vector>& increment,                    //!< increment vector
    const Teuchos::RCP<Epetra_Vector>& residual,                     //!< residual vector
    const Teuchos::RCP<Epetra_Vector>& phinp,                        //!< state vector at time n+1
    const int iteration,  //!< number of current Newton-Raphson iteration
    CORE::LINALG::SolverParams& solver_params) const
{
  solver_params.refactor = true;
  solver_params.reset = iteration == 1;
  solver->Solve(systemmatrix->EpetraOperator(), increment, residual, solver_params);

  return;
}  // SCATRA::MeshtyingStrategyStd::Solve


/*-------------------------------------------------------------------------*
 | return linear solver for global system of linear equations   fang 01/18 |
 *-------------------------------------------------------------------------*/
const CORE::LINALG::Solver& SCATRA::MeshtyingStrategyStd::Solver() const
{
  if (scatratimint_->Solver() == Teuchos::null) FOUR_C_THROW("Invalid linear solver!");
  return *scatratimint_->Solver();
}  // SCATRA::MeshtyingStrategyStd::Solver()


/*------------------------------------------------------------------------*
 | instantiate strategy for Newton-Raphson convergence check   fang 02/16 |
 *------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyStd::init_conv_check_strategy()
{
  if (scatratimint_->MicroScale())
    convcheckstrategy_ = Teuchos::rcp(new SCATRA::ConvCheckStrategyStdMicroScale(
        scatratimint_->ScatraParameterList()->sublist("NONLINEAR")));
  else if (GLOBAL::Problem::Instance()->GetProblemType() ==
           GLOBAL::ProblemType::poromultiphasescatra)
    convcheckstrategy_ = Teuchos::rcp(new SCATRA::ConvCheckStrategyPoroMultiphaseScatra(
        scatratimint_->ScatraParameterList()->sublist("NONLINEAR")));
  else
    convcheckstrategy_ = Teuchos::rcp(new SCATRA::ConvCheckStrategyStd(
        scatratimint_->ScatraParameterList()->sublist("NONLINEAR")));

  return;
}  // SCATRA::MeshtyingStrategyStd::init_conv_check_strategy

FOUR_C_NAMESPACE_CLOSE
