/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid-fluid meshtying strategy for standard scalar transport problems

\level 2


*----------------------------------------------------------------------*/
#include "4C_scatra_timint_meshtying_strategy_fluid.hpp"

#include "4C_fluid_meshtying.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_sparseoperator.hpp"
#include "4C_scatra_timint_implicit.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | constructor                                               fang 12/14 |
 *----------------------------------------------------------------------*/
ScaTra::MeshtyingStrategyFluid::MeshtyingStrategyFluid(ScaTra::ScaTraTimIntImpl* scatratimint)
    : MeshtyingStrategyBase(scatratimint),
      meshtying_(Teuchos::null),
      type_(Inpar::FLUID::no_meshtying)
{
  return;
}  // ScaTra::MeshtyingStrategyFluid::MeshtyingStrategyFluid


/*-----------------------------------------------------------------------*
 | return global map of degrees of freedom                    fang 02/18 |
 *-----------------------------------------------------------------------*/
const Epetra_Map& ScaTra::MeshtyingStrategyFluid::dof_row_map() const
{
  return *scatratimint_->dof_row_map();
}


/*----------------------------------------------------------------------*
 | evaluate fluid-fluid meshtying                            fang 12/14 |
 *----------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyFluid::evaluate_meshtying()
{
  // need to complete system matrix due to subsequent matrix-matrix multiplications
  scatratimint_->system_matrix_operator()->complete();

  // evaluate fluid-fluid meshtying
  meshtying_->prepare_meshtying_system(
      scatratimint_->system_matrix_operator(), scatratimint_->residual(), scatratimint_->phinp());

  return;
}  // ScaTra::MeshtyingStrategyFluid::evaluate_meshtying


/*----------------------------------------------------------------------*
 | include Dirichlet conditions into condensation            fang 12/14 |
 *----------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyFluid::include_dirichlet_in_condensation() const
{
  meshtying_->include_dirichlet_in_condensation(scatratimint_->phinp(), scatratimint_->phin());

  return;
}  // ScaTra::MeshtyingStrategyFluid::include_dirichlet_in_condensation()


/*----------------------------------------------------------------------*
 | perform setup of fluid-fluid meshtying                    fang 12/14 |
 *----------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyFluid::setup_meshtying()
{
  // safety check
  if (scatratimint_->num_scal() < 1)
    FOUR_C_THROW("Number of transported scalars not correctly set!");

  // define coupling and initialize system matrix
  std::vector<int> coupleddof(scatratimint_->num_scal(), 1);

  meshtying_->setup_meshtying(coupleddof);
}  // ScaTra::MeshtyingStrategyFluid::setup_meshtying


/*----------------------------------------------------------------------*
 | perform init of fluid-fluid meshtying                    rauch 09/16 |
 *----------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyFluid::init_meshtying()
{
  // instantiate strategy for Newton-Raphson convergence check
  init_conv_check_strategy();

  // Important: Meshtying for scalar transport is not well tested!
  // get meshtying type
  type_ = Core::UTILS::IntegralValue<Inpar::FLUID::MeshTying>(
      *(scatratimint_->scatra_parameter_list()), "MESHTYING");

  // safety checks
  if (type_ == Inpar::FLUID::condensed_bmat)
    FOUR_C_THROW(
        "The 2x2 block solver algorithm for a block matrix system has not been activated yet. Just "
        "do it!");

  // setup meshtying
  meshtying_ = Teuchos::rcp(new FLD::Meshtying(scatratimint_->discretization(),
      *(scatratimint_->solver()), type_, Global::Problem::instance()->n_dim()));

  return;
}  // ScaTra::MeshtyingStrategyFluid::InitMeshtying


/*----------------------------------------------------------------------*
 | initialize system matrix for fluid-fluid meshtying        fang 12/14 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::SparseOperator> ScaTra::MeshtyingStrategyFluid::init_system_matrix()
    const
{
  return meshtying_->init_system_matrix();
}  // ScaTra::MeshtyingStrategyFluid::init_system_matrix


/*-------------------------------------------------------------------------*
 | solve linear system of equations for fluid-fluid meshtying   fang 12/14 |
 *-------------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyFluid::solve(
    const Teuchos::RCP<Core::LinAlg::Solver>& solver,                //!< solver
    const Teuchos::RCP<Core::LinAlg::SparseOperator>& systemmatrix,  //!< system matrix
    const Teuchos::RCP<Epetra_Vector>& increment,                    //!< increment vector
    const Teuchos::RCP<Epetra_Vector>& residual,                     //!< residual vector
    const Teuchos::RCP<Epetra_Vector>& phinp,                        //!< state vector at time n+1
    const int iteration,  //!< number of current Newton-Raphson iteration
    Core::LinAlg::SolverParams& solver_params) const
{
  meshtying_->solve_meshtying(
      *solver, systemmatrix, increment, residual, phinp, iteration, solver_params);

  return;
}  // ScaTra::MeshtyingStrategyFluid::Solve


/*-------------------------------------------------------------------------*
 | return linear solver for global system of linear equations   fang 01/18 |
 *-------------------------------------------------------------------------*/
const Core::LinAlg::Solver& ScaTra::MeshtyingStrategyFluid::solver() const
{
  if (scatratimint_->solver() == Teuchos::null) FOUR_C_THROW("Invalid linear solver!");
  return *scatratimint_->solver();
}


/*------------------------------------------------------------------------*
 | instantiate strategy for Newton-Raphson convergence check   fang 02/16 |
 *------------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyFluid::init_conv_check_strategy()
{
  convcheckstrategy_ = Teuchos::rcp(new ScaTra::ConvCheckStrategyStd(
      scatratimint_->scatra_parameter_list()->sublist("NONLINEAR")));

  return;
}  // ScaTra::MeshtyingStrategyFluid::init_conv_check_strategy

FOUR_C_NAMESPACE_CLOSE
