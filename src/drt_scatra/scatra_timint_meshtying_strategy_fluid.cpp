/*!----------------------------------------------------------------------
\file scatra_timint_meshtying_strategy_fluid.cpp

\brief Fluid-fluid meshtying strategy for standard scalar transport problems

\level 2

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
</pre>

*----------------------------------------------------------------------*/
#include "scatra_timint_meshtying_strategy_fluid.H"

#include "../drt_fluid/fluid_meshtying.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_scatra/scatra_timint_implicit.H"

#include "../linalg/linalg_sparseoperator.H"

/*----------------------------------------------------------------------*
 | constructor                                               fang 12/14 |
 *----------------------------------------------------------------------*/
SCATRA::MeshtyingStrategyFluid::MeshtyingStrategyFluid(
    SCATRA::ScaTraTimIntImpl* scatratimint
    ) :
MeshtyingStrategyBase(scatratimint),
meshtying_(Teuchos::null),
type_(INPAR::FLUID::no_meshtying)
{
  return;
} // SCATRA::MeshtyingStrategyFluid::MeshtyingStrategyFluid


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
  meshtying_->PrepareMeshtyingSystem(scatratimint_->SystemMatrixOperator(),scatratimint_->Residual(),scatratimint_->Phinp());

  return;
} // SCATRA::MeshtyingStrategyFluid::EvaluateMeshtying


/*----------------------------------------------------------------------*
 | include Dirichlet conditions into condensation            fang 12/14 |
 *----------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyFluid::IncludeDirichletInCondensation() const
{
  meshtying_->IncludeDirichletInCondensation(scatratimint_->Phinp(),scatratimint_->Phin());

  return;
} // SCATRA::MeshtyingStrategyFluid::IncludeDirichletInCondensation()


/*----------------------------------------------------------------------*
 | perform setup of fluid-fluid meshtying                    fang 12/14 |
 *----------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyFluid::SetupMeshtying()
{
  return;
} // SCATRA::MeshtyingStrategyFluid::SetupMeshtying


/*----------------------------------------------------------------------*
 | perform init of fluid-fluid meshtying                    rauch 09/16 |
 *----------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyFluid::InitMeshtying()
{
  // instantiate strategy for Newton-Raphson convergence check
  InitConvCheckStrategy();

  // Important: Meshtying for scalar transport is not well tested!
  // get meshtying type
  type_ = DRT::INPUT::IntegralValue<INPAR::FLUID::MeshTying>(*(scatratimint_->ScatraParameterList()),"MESHTYING");

  // safety checks
  if(type_ == INPAR::FLUID::condensed_bmat)
    dserror("The 2x2 block solver algorithm for a block matrix system has not been activated yet. Just do it!");

  // setup meshtying
  meshtying_ = Teuchos::rcp(new FLD::Meshtying(scatratimint_->Discretization(),*(scatratimint_->Solver()),type_,DRT::Problem::Instance()->NDim()));

  return;
} // SCATRA::MeshtyingStrategyFluid::InitMeshtying


/*----------------------------------------------------------------------*
 | initialize system matrix for fluid-fluid meshtying        fang 12/14 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseOperator> SCATRA::MeshtyingStrategyFluid::InitSystemMatrix() const
{
  // safety check
  if(scatratimint_->NumScal() < 1)
    dserror("Number of transported scalars not correctly set!");

  // define coupling and initialize system matrix
  std::vector<int> coupleddof(scatratimint_->NumScal(),1);

  return meshtying_->Setup(coupleddof);
} // SCATRA::MeshtyingStrategyFluid::InitSystemMatrix


/*-------------------------------------------------------------------------*
 | solve linear system of equations for fluid-fluid meshtying   fang 12/14 |
 *-------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyFluid::Solve(
    const Teuchos::RCP<LINALG::Solver>&            solver,         //!< solver
    const Teuchos::RCP<LINALG::SparseOperator>&    systemmatrix,   //!< system matrix
    const Teuchos::RCP<Epetra_Vector>&             increment,      //!< increment vector
    const Teuchos::RCP<Epetra_Vector>&             residual,       //!< residual vector
    const Teuchos::RCP<Epetra_Vector>&             phinp,          //!< state vector at time n+1
    const int&                                     iteration,      //!< number of current Newton-Raphson iteration
    const Teuchos::RCP<LINALG::KrylovProjector>&   projector       //!< Krylov projector
    ) const
{
  meshtying_->SolveMeshtying(*solver,systemmatrix,increment,residual,phinp,iteration,projector);

  return;
} // SCATRA::MeshtyingStrategyFluid::Solve


/*-------------------------------------------------------------------------*
 | return linear solver for global system of linear equations   fang 01/18 |
 *-------------------------------------------------------------------------*/
const LINALG::Solver& SCATRA::MeshtyingStrategyFluid::Solver() const
{
  if(scatratimint_->Solver() == Teuchos::null)
    dserror("Invalid linear solver!");
  return *scatratimint_->Solver();
}


/*------------------------------------------------------------------------*
 | instantiate strategy for Newton-Raphson convergence check   fang 02/16 |
 *------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyFluid::InitConvCheckStrategy()
{
  convcheckstrategy_ = Teuchos::rcp(new SCATRA::ConvCheckStrategyStd(scatratimint_->ScatraParameterList()->sublist("NONLINEAR")));

  return;
} // SCATRA::MeshtyingStrategyFluid::InitConvCheckStrategy
