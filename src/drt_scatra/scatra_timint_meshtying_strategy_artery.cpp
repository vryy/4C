/*----------------------------------------------------------------------*/
/*!
 \file scatra_timint_meshtying_strategy_artery.cpp

 \brief routines for coupling between 1D arterial network and 2D/3D
        scatra-algorithm

   \level 3

   \maintainer  Johannes Kremheller
                kremheller@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/

#include "scatra_timint_meshtying_strategy_artery.H"
#include "scatra_timint_implicit.H"

#include "../drt_adapter/adapter_scatra_base_algorithm.H"
#include "../drt_adapter/ad_art_net.H"
#include "../linalg/linalg_solver.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_poromultiphase_scatra/poromultiphase_scatra_artery_coupling_nodebased.H"
#include "../drt_poromultiphase_scatra/poromultiphase_scatra_utils.H"


/*----------------------------------------------------------------------*
 | constructor                                         kremheller 04/18 |
 *----------------------------------------------------------------------*/
SCATRA::MeshtyingStrategyArtery::MeshtyingStrategyArtery(
    SCATRA::ScaTraTimIntImpl* scatratimint  //!< scalar transport time integrator
    )
    : MeshtyingStrategyBase(scatratimint)
{
}

/*----------------------------------------------------------------------*
 | init                                                kremheller 04/18 |
 *----------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyArtery::InitMeshtying()
{
  // instantiate strategy for Newton-Raphson convergence check
  InitConvCheckStrategy();

  const Teuchos::ParameterList& globaltimeparams =
      DRT::Problem::Instance()->PoroMultiPhaseScatraDynamicParams();
  const Teuchos::ParameterList& myscatraparams =
      DRT::Problem::Instance()->ScalarTransportDynamicParams();
  if (DRT::INPUT::IntegralValue<INPAR::SCATRA::VelocityField>(myscatraparams, "VELOCITYFIELD") !=
      INPAR::SCATRA::velocity_zero)
    dserror("set your velocity field to zero!");

  // artery scatra problem
  Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> art_scatra =
      Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm());

  // initialize the base algo.
  // scatra time integrator is constructed and initialized inside.
  art_scatra->Init(globaltimeparams, myscatraparams,
      DRT::Problem::Instance()->SolverParams(myscatraparams.get<int>("LINEAR_SOLVER")),
      "artery_scatra", false);

  // only now we must call Setup() on the scatra time integrator.
  // all objects relying on the parallel distribution are
  // created and pointers are set.
  // calls Setup() on the scatra time integrator inside.
  art_scatra->ScaTraField()->Setup();
  DRT::Problem::Instance()->AddFieldTest(art_scatra->CreateScaTraFieldTest());

  // set the time integrator
  SetArteryScatraTimeIntegrator(art_scatra->ScaTraField());

  // get the two discretizations
  artscatradis_ = artscatratimint_->Discretization();
  scatradis_ = scatratimint_->Discretization();

  if (scatratimint_->Discretization()->Comm().MyPID() == 0)
  {
    std::cout << "\n";
    std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    std::cout << "<                                                  >" << std::endl;
    std::cout << "< ScaTra-Coupling with 1D Artery Network activated >" << std::endl;
  }

  // init the mesh tying object, which does all the work
  arttoscatracoupling_ = POROMULTIPHASESCATRA::UTILS::CreateAndInitArteryCouplingStrategy(
      artscatradis_, scatradis_, myscatraparams.sublist("ARTERY COUPLING"), "ArtScatraCouplCon",
      "COUPLEDDOFS_ARTSCATRA", "COUPLEDDOFS_SCATRA");

  return;
}

/*----------------------------------------------------------------------*
 | setup                                               kremheller 04/18 |
 *----------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyArtery::SetupMeshtying()
{
  // Initialize rhs vector
  rhs_ = Teuchos::rcp(new Epetra_Vector(*arttoscatracoupling_->FullMap(), true));

  // Initialize increment vector
  comb_increment_ = Teuchos::rcp(new Epetra_Vector(*arttoscatracoupling_->FullMap(), true));

  // initialize scatra-artery_scatra-systemmatrix_
  comb_systemmatrix_ =
      Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
          *arttoscatracoupling_->GlobalExtractor(), *arttoscatracoupling_->GlobalExtractor(), 81,
          false, true));

  arttoscatracoupling_->Setup();

  return;
}

/*-----------------------------------------------------------------------*
 | return global map of degrees of freedom              kremheller 04/18 |
 *-----------------------------------------------------------------------*/
const Epetra_Map& SCATRA::MeshtyingStrategyArtery::DofRowMap() const
{
  return *arttoscatracoupling_->FullMap();
}

/*-----------------------------------------------------------------------*
 | return global map of degrees of freedom              kremheller 04/18 |
 *-----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> SCATRA::MeshtyingStrategyArtery::ArtScatraDofRowMap() const
{
  return arttoscatracoupling_->ArteryDofRowMap();
}

/*----------------------------------------------------------------------*
 | evaluate                                            kremheller 04/18 |
 *----------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyArtery::EvaluateMeshtying()
{
  // here we just assemble matrix and rhs of artery-scatra problem
  // actual coupling (meshtying) is evaluated in Solve
  // reason for that is that we need the system matrix of the continuous scatra
  // problem with DBCs applied which is performed directly before calling solve

  artscatratimint_->PrepareLinearSolve();

  return;
}


/*----------------------------------------------------------------------------------*
 | initialize system matrix for scatra-artery interface coupling   kremheller 04/18 |
 *------------------------------------------------------------------------------    */
Teuchos::RCP<LINALG::SparseOperator> SCATRA::MeshtyingStrategyArtery::InitSystemMatrix() const
{
  return Teuchos::rcp(
      new LINALG::SparseMatrix(*(scatratimint_->Discretization()->DofRowMap()), 27, false, true));
  ;
}

/*-------------------------------------------------------------------------------*
 | return linear solver for global system of linear equations   kremheller 04/18 |
 *-------------------------------------------------------------------------------*/
const LINALG::Solver& SCATRA::MeshtyingStrategyArtery::Solver() const
{
  if (scatratimint_->Solver() == Teuchos::null) dserror("Invalid linear solver!");

  return *scatratimint_->Solver();
}

/*------------------------------------------------------------------------------*
 | instantiate strategy for Newton-Raphson convergence check   kremheller 04/18 |
 *------------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyArtery::InitConvCheckStrategy()
{
  convcheckstrategy_ = Teuchos::rcp(new SCATRA::ConvCheckStrategyPoroMultiphaseScatraArtMeshTying(
      scatratimint_->ScatraParameterList()->sublist("NONLINEAR")));

  return;
}

/*------------------------------------------------------------------------------------------*
 | solve linear system of equations for scatra-scatra interface coupling   kremheller 04/18 |
 *------------------------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyArtery::Solve(const Teuchos::RCP<LINALG::Solver>& solver,  //!< solver
    const Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,  //!< system matrix
    const Teuchos::RCP<Epetra_Vector>& increment,              //!< increment vector
    const Teuchos::RCP<Epetra_Vector>& residual,               //!< residual vector
    const Teuchos::RCP<Epetra_Vector>& phinp,                  //!< state vector at time n+1
    const int& iteration,  //!< number of current Newton-Raphson iteration
    const Teuchos::RCP<LINALG::KrylovProjector>& projector  //!< Krylov projector
    ) const
{
  // setup the system (evaluate mesh tying)
  // reason for this being done here is that we need the system matrix of the continuous scatra
  // problem with DBCs applied which is performed directly before calling solve

  SetupSystem(systemmatrix, residual);

  comb_systemmatrix_->Complete();

  // solve
  comb_increment_->PutScalar(0.0);
  solver->Solve(
      comb_systemmatrix_->EpetraOperator(), comb_increment_, rhs_, true, iteration == 1, projector);

  // extract increments of scatra and artery-scatra field
  Teuchos::RCP<const Epetra_Vector> artscatrainc;
  Teuchos::RCP<const Epetra_Vector> myinc;
  arttoscatracoupling_->ExtractSingleFieldVectors(comb_increment_, myinc, artscatrainc);

  // update the scatra increment, update iter is performed outside
  increment->Update(1.0, *(myinc), 1.0);
  // update the artery-scatra field
  artscatratimint_->UpdateIter(artscatrainc);

  return;
}

/*------------------------------------------------------------------------------------------*
 | solve linear system of equations for scatra-scatra interface coupling   kremheller 04/18 |
 *------------------------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyArtery::SetupSystem(
    const Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,  //!< system matrix
    const Teuchos::RCP<Epetra_Vector>& residual                //!< residual vector
    ) const
{
  arttoscatracoupling_->SetSolutionVectors(
      scatratimint_->Phinp(), Teuchos::null, artscatratimint_->Phinp());

  arttoscatracoupling_->SetupSystem(comb_systemmatrix_, rhs_,
      Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(systemmatrix),
      artscatratimint_->SystemMatrix(), residual, artscatratimint_->Residual(),
      scatratimint_->DirichMaps(), artscatratimint_->DirichMaps());
}

/*-------------------------------------------------------------------------*
 | set time integrator for scalar transport in arteries   kremheller 04/18 |
 *------------------------------------------------------------------------ */
void SCATRA::MeshtyingStrategyArtery::UpdateArtScatraIter(
    Teuchos::RCP<const Epetra_Vector> combined_inc)
{
  Teuchos::RCP<const Epetra_Vector> artscatrainc;
  Teuchos::RCP<const Epetra_Vector> myinc;
  arttoscatracoupling_->ExtractSingleFieldVectors(combined_inc, myinc, artscatrainc);

  artscatratimint_->UpdateIter(artscatrainc);

  return;
}

/*-------------------------------------------------------------------------*
 | set time integrator for scalar transport in arteries   kremheller 04/18 |
 *------------------------------------------------------------------------ */
void SCATRA::MeshtyingStrategyArtery::SetArteryScatraTimeIntegrator(
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> artscatratimint)
{
  artscatratimint_ = artscatratimint;
  if (artscatratimint_ == Teuchos::null) dserror("could not set artery scatra time integrator");

  return;
}

/*-------------------------------------------------------------------------*
 | set time integrator for artery problems                kremheller 04/18 |
 *------------------------------------------------------------------------ */
void SCATRA::MeshtyingStrategyArtery::SetArteryTimeIntegrator(
    Teuchos::RCP<ADAPTER::ArtNet> arttimint)
{
  arttimint_ = arttimint;
  if (arttimint_ == Teuchos::null) dserror("could not set artery time integrator");

  return;
}

/*-------------------------------------------------------------------------*
 | set element pairs that are close                       kremheller 03/19 |
 *------------------------------------------------------------------------ */
void SCATRA::MeshtyingStrategyArtery::SetNearbyElePairs(
    const std::map<int, std::set<int>>* nearbyelepairs)
{
  arttoscatracoupling_->SetNearbyElePairs(nearbyelepairs);

  return;
}

/*--------------------------------------------------------------------------*
 | setup the coupled matrix                                kremheller 04/18 |
 *--------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyArtery::PrepareTimeStep() const
{
  artscatratimint_->PrepareTimeStep();
  return;
}

/*--------------------------------------------------------------------------*
 | setup the coupled matrix                                kremheller 04/18 |
 *--------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyArtery::SetArteryPressure() const
{
  artscatradis_->SetState(2, "one_d_artery_pressure", arttimint_->Pressurenp());
  return;
}

/*--------------------------------------------------------------------------*
 | apply mesh movement on artery coupling                  kremheller 07/18 |
 *--------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyArtery::ApplyMeshMovement()
{
  arttoscatracoupling_->ApplyMeshMovement();
  return;
}

/*--------------------------------------------------------------------------*
 | check if initial fields match                           kremheller 04/18 |
 *--------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyArtery::CheckInitialFields() const
{
  arttoscatracoupling_->CheckInitialFields(scatratimint_->Phinp(), artscatratimint_->Phinp());

  return;
}
