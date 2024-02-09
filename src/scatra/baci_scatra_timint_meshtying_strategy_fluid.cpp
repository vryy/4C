/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid-fluid meshtying strategy for standard scalar transport problems

\level 2


*----------------------------------------------------------------------*/
#include "baci_scatra_timint_meshtying_strategy_fluid.hpp"

#include "baci_fluid_meshtying.hpp"
#include "baci_global_data.hpp"
#include "baci_linalg_sparseoperator.hpp"
#include "baci_scatra_timint_implicit.hpp"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | constructor                                               fang 12/14 |
 *----------------------------------------------------------------------*/
SCATRA::MeshtyingStrategyFluid::MeshtyingStrategyFluid(SCATRA::ScaTraTimIntImpl* scatratimint)
    : MeshtyingStrategyBase(scatratimint),
      meshtying_(Teuchos::null),
      type_(INPAR::FLUID::no_meshtying)
{
  return;
}  // SCATRA::MeshtyingStrategyFluid::MeshtyingStrategyFluid


/*-----------------------------------------------------------------------*
 | return global map of degrees of freedom                    fang 02/18 |
 *-----------------------------------------------------------------------*/
const Epetra_Map& SCATRA::MeshtyingStrategyFluid::DofRowMap() const
{
  return *scatratimint_->DofRowMap();
}


/*----------------------------------------------------------------------*
 | evaluate fluid-fluid meshtying                            fang 12/14 |
 *----------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyFluid::EvaluateMeshtying()
{
  // need to complete system matrix due to subsequent matrix-matrix multiplications
  scatratimint_->SystemMatrixOperator()->Complete();

  // evaluate fluid-fluid meshtying
  meshtying_->PrepareMeshtyingSystem(
      scatratimint_->SystemMatrixOperator(), scatratimint_->Residual(), scatratimint_->Phinp());

  return;
}  // SCATRA::MeshtyingStrategyFluid::EvaluateMeshtying


/*----------------------------------------------------------------------*
 | include Dirichlet conditions into condensation            fang 12/14 |
 *----------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyFluid::IncludeDirichletInCondensation() const
{
  meshtying_->IncludeDirichletInCondensation(scatratimint_->Phinp(), scatratimint_->Phin());

  return;
}  // SCATRA::MeshtyingStrategyFluid::IncludeDirichletInCondensation()


/*----------------------------------------------------------------------*
 | perform setup of fluid-fluid meshtying                    fang 12/14 |
 *----------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyFluid::SetupMeshtying()
{
  // safety check
  if (scatratimint_->NumScal() < 1) dserror("Number of transported scalars not correctly set!");

  // define coupling and initialize system matrix
  std::vector<int> coupleddof(scatratimint_->NumScal(), 1);

  meshtying_->SetupMeshtying(coupleddof);
}  // SCATRA::MeshtyingStrategyFluid::SetupMeshtying


/*----------------------------------------------------------------------*
 | perform init of fluid-fluid meshtying                    rauch 09/16 |
 *----------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyFluid::InitMeshtying()
{
  // instantiate strategy for Newton-Raphson convergence check
  InitConvCheckStrategy();

  // Important: Meshtying for scalar transport is not well tested!
  // get meshtying type
  type_ = INPUT::IntegralValue<INPAR::FLUID::MeshTying>(
      *(scatratimint_->ScatraParameterList()), "MESHTYING");

  // safety checks
  if (type_ == INPAR::FLUID::condensed_bmat)
    dserror(
        "The 2x2 block solver algorithm for a block matrix system has not been activated yet. Just "
        "do it!");

  // setup meshtying
  meshtying_ = Teuchos::rcp(new FLD::Meshtying(scatratimint_->Discretization(),
      *(scatratimint_->Solver()), type_, GLOBAL::Problem::Instance()->NDim()));

  return;
}  // SCATRA::MeshtyingStrategyFluid::InitMeshtying


/*----------------------------------------------------------------------*
 | initialize system matrix for fluid-fluid meshtying        fang 12/14 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<CORE::LINALG::SparseOperator> SCATRA::MeshtyingStrategyFluid::InitSystemMatrix() const
{
  return meshtying_->InitSystemMatrix();
}  // SCATRA::MeshtyingStrategyFluid::InitSystemMatrix


/*-------------------------------------------------------------------------*
 | solve linear system of equations for fluid-fluid meshtying   fang 12/14 |
 *-------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyFluid::Solve(
    const Teuchos::RCP<CORE::LINALG::Solver>& solver,                //!< solver
    const Teuchos::RCP<CORE::LINALG::SparseOperator>& systemmatrix,  //!< system matrix
    const Teuchos::RCP<Epetra_Vector>& increment,                    //!< increment vector
    const Teuchos::RCP<Epetra_Vector>& residual,                     //!< residual vector
    const Teuchos::RCP<Epetra_Vector>& phinp,                        //!< state vector at time n+1
    const int iteration,  //!< number of current Newton-Raphson iteration
    CORE::LINALG::SolverParams& solver_params) const
{
  meshtying_->SolveMeshtying(
      *solver, systemmatrix, increment, residual, phinp, iteration, solver_params);

  return;
}  // SCATRA::MeshtyingStrategyFluid::Solve


/*-------------------------------------------------------------------------*
 | return linear solver for global system of linear equations   fang 01/18 |
 *-------------------------------------------------------------------------*/
const CORE::LINALG::Solver& SCATRA::MeshtyingStrategyFluid::Solver() const
{
  if (scatratimint_->Solver() == Teuchos::null) dserror("Invalid linear solver!");
  return *scatratimint_->Solver();
}


/*------------------------------------------------------------------------*
 | instantiate strategy for Newton-Raphson convergence check   fang 02/16 |
 *------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyFluid::InitConvCheckStrategy()
{
  convcheckstrategy_ = Teuchos::rcp(
      new SCATRA::ConvCheckStrategyStd(scatratimint_->ScatraParameterList()->sublist("NONLINEAR")));

  return;
}  // SCATRA::MeshtyingStrategyFluid::InitConvCheckStrategy

BACI_NAMESPACE_CLOSE
