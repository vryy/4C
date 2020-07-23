/*----------------------------------------------------------------------*/
/*! \file

\brief routines for coupling with artery network

\level 3

*----------------------------------------------------------------------*/

#include "porofluidmultiphase_meshtying_strategy_artery.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_adapter/ad_art_net.H"
#include "../drt_inpar/inpar_bio.H"
#include "../drt_art_net/art_net_utils.H"
#include "../drt_poromultiphase_scatra/poromultiphase_scatra_artery_coupling_base.H"
#include "../drt_poromultiphase_scatra/poromultiphase_scatra_utils.H"
#include "../linalg/linalg_solver.H"
#include "porofluidmultiphase_utils.H"
#include "../linalg/linalg_utils_sparse_algebra_print.H"


#include "../drt_io/io.H"
#include "../drt_io/io_control.H"

#include <Teuchos_TimeMonitor.hpp>


/*----------------------------------------------------------------------*
 | constructor                                (public) kremheller 04/18 |
 *----------------------------------------------------------------------*/
POROFLUIDMULTIPHASE::MeshtyingStrategyArtery::MeshtyingStrategyArtery(
    POROFLUIDMULTIPHASE::TimIntImpl* porofluidmultitimint, const Teuchos::ParameterList& probparams,
    const Teuchos::ParameterList& poroparams)
    : MeshtyingStrategyBase(porofluidmultitimint, probparams, poroparams)
{
  const Teuchos::ParameterList& artdyn = DRT::Problem::Instance()->ArterialDynamicParams();

  arterydis_ = DRT::Problem::Instance()->GetDis("artery");

  if (!arterydis_->Filled()) arterydis_->FillComplete();

  INPAR::ARTDYN::TimeIntegrationScheme timintscheme =
      DRT::INPUT::IntegralValue<INPAR::ARTDYN::TimeIntegrationScheme>(artdyn, "DYNAMICTYP");

  Teuchos::RCP<IO::DiscretizationWriter> artery_output = arterydis_->Writer();
  artery_output->WriteMesh(0, 0.0);

  // build art net time integrator
  artnettimint_ =
      ART::UTILS::CreateAlgorithm(timintscheme, arterydis_, artdyn.get<int>("LINEAR_SOLVER"),
          probparams, artdyn, DRT::Problem::Instance()->ErrorFile()->Handle(), artery_output);

  // set to false
  artnettimint_->SetSolveScatra(false);

  // initialize
  artnettimint_->Init(probparams, artdyn, "artery_scatra");

  // print user info
  if (porofluidmultitimint->Discretization()->Comm().MyPID() == 0)
  {
    std::cout << "\n";
    std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    std::cout << "<                                                  >" << std::endl;
    std::cout << "<    Coupling with 1D Artery Network activated     >" << std::endl;
  }

  // initialize mesh tying object
  arttoporofluidcoupling_ = POROMULTIPHASESCATRA::UTILS::CreateAndInitArteryCouplingStrategy(
      arterydis_, porofluidmultitimint->Discretization(), poroparams.sublist("ARTERY COUPLING"),
      "ArtPorofluidCouplCon", "COUPLEDDOFS_ART", "COUPLEDDOFS_PORO");

  // Initialize rhs vector
  rhs_ = Teuchos::rcp(new Epetra_Vector(*arttoporofluidcoupling_->FullMap(), true));

  // Initialize increment vector
  comb_increment_ = Teuchos::rcp(new Epetra_Vector(*arttoporofluidcoupling_->FullMap(), true));
  // Initialize phinp vector
  comb_phinp_ = Teuchos::rcp(new Epetra_Vector(*arttoporofluidcoupling_->FullMap(), true));

  // initialize Poromultiphase-elasticity-systemmatrix_
  comb_systemmatrix_ =
      Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
          *arttoporofluidcoupling_->GlobalExtractor(), *arttoporofluidcoupling_->GlobalExtractor(),
          81, false, true));

  return;
}


/*----------------------------------------------------------------------*
| Destructor dtor (public)                             kremheller 04/18 |
*-----------------------------------------------------------------------*/
POROFLUIDMULTIPHASE::MeshtyingStrategyArtery::~MeshtyingStrategyArtery() { return; }

/*----------------------------------------------------------------------*
 | prepare time loop                                   kremheller 04/18 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::MeshtyingStrategyArtery::PrepareTimeLoop()
{
  artnettimint_->PrepareTimeLoop();
  return;
}

/*----------------------------------------------------------------------*
 | setup the variables to do a new time step  (public) kremheller 04/18 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::MeshtyingStrategyArtery::PrepareTimeStep()
{
  artnettimint_->PrepareTimeStep();
  return;
}

/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                     kremheller 04/18 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::MeshtyingStrategyArtery::Update()
{
  artnettimint_->TimeUpdate();
  return;
}

/*--------------------------------------------------------------------------*
 | initialize the linear solver                            kremheller 07/20 |
 *--------------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::MeshtyingStrategyArtery::InitializeLinearSolver(
    Teuchos::RCP<LINALG::Solver> solver)
{
  const Teuchos::ParameterList& porofluidparams =
      DRT::Problem::Instance()->PoroFluidMultiPhaseDynamicParams();
  const int linsolvernumber = porofluidparams.get<int>("LINEAR_SOLVER");
  const Teuchos::ParameterList& solverparams =
      DRT::Problem::Instance()->SolverParams(linsolvernumber);
  const int solvertype =
      DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(solverparams, "SOLVER");
  // no need to do the rest for direct solvers
  if (solvertype == INPAR::SOLVER::umfpack or solvertype == INPAR::SOLVER::superlu or
      solvertype == INPAR::SOLVER::amesos_klu_nonsym)
    return;

  if (solvertype != INPAR::SOLVER::aztec_msr && solvertype != INPAR::SOLVER::belos)
    dserror("aztec solver expected");

  const int azprectype =
      DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(solverparams, "AZPREC");

  // plausibility check
  switch (azprectype)
  {
    case INPAR::SOLVER::azprec_AMGnxn:
    {
      // no plausibility checks here
      // if you forget to declare an xml file you will get an error message anyway
    }
    break;
    default:
      dserror("AMGnxn preconditioner expected");
      break;
  }

  // equip smoother for fluid matrix block with empty parameter sublists to trigger null space
  // computation
  Teuchos::ParameterList& blocksmootherparams1 = solver->Params().sublist("Inverse1");
  blocksmootherparams1.sublist("Aztec Parameters");
  blocksmootherparams1.sublist("MueLu Parameters");

  porofluidmultitimint_->Discretization()->ComputeNullSpaceIfNecessary(blocksmootherparams1);

  Teuchos::ParameterList& blocksmootherparams2 = solver->Params().sublist("Inverse2");
  blocksmootherparams2.sublist("Aztec Parameters");
  blocksmootherparams2.sublist("MueLu Parameters");

  arterydis_->ComputeNullSpaceIfNecessary(blocksmootherparams2);
}

/*--------------------------------------------------------------------------*
 | solve linear system of equations                        kremheller 04/18 |
 *--------------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::MeshtyingStrategyArtery::LinearSolve(Teuchos::RCP<LINALG::Solver> solver,
    Teuchos::RCP<LINALG::SparseOperator> sysmat, Teuchos::RCP<Epetra_Vector> increment,
    Teuchos::RCP<Epetra_Vector> residual)
{
  bool matlab = false;
  if (matlab)
  {
    // sparse_matrix
    std::string filename = "../o/mymatrix.dat";
    std::string filename_vc = "../o/myvec.dat";
    LINALG::PrintBlockMatrixInMatlabFormat(filename, *(comb_systemmatrix_));
    LINALG::PrintVectorInMatlabFormat(filename_vc, *rhs_, true);
    dserror("exit");
  }

  comb_systemmatrix_->Complete();

  comb_increment_->PutScalar(0.0);

  // standard solver call
  // system is ready to solve since Dirichlet Boundary conditions have been applied in
  // SetupSystemMatrix or Evaluate
  solver->Solve(comb_systemmatrix_->EpetraOperator(), comb_increment_, rhs_, true, false);

  return;
}

/*----------------------------------------------------------------------*
 | Calculate problem specific norm                     kremheller 03/18 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::MeshtyingStrategyArtery::CalculateNorms(double& preresnorm,
    double& incprenorm, double& prenorm, const Teuchos::RCP<const Epetra_Vector> increment)
{
  preresnorm = UTILS::CalculateVectorNorm(vectornormfres_, rhs_);
  incprenorm = UTILS::CalculateVectorNorm(vectornorminc_, comb_increment_);
  // setup combined solution vector and calculate norm
  arttoporofluidcoupling_->SetupVector(
      comb_phinp_, porofluidmultitimint_->Phinp(), artnettimint_->Pressurenp());
  prenorm = UTILS::CalculateVectorNorm(vectornorminc_, comb_phinp_);

  return;
}

/*----------------------------------------------------------------------*
 | create result test for this field                   kremheller 04/18 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::MeshtyingStrategyArtery::CreateFieldTest()
{
  Teuchos::RCP<DRT::ResultTest> arteryresulttest = artnettimint_->CreateFieldTest();
  DRT::Problem::Instance()->AddFieldTest(arteryresulttest);
  return;
}

/*----------------------------------------------------------------------*
 |  read restart data                                  kremheller 04/18 |
 -----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::MeshtyingStrategyArtery::ReadRestart(const int step)
{
  artnettimint_->ReadRestart(step);

  return;
}

/*----------------------------------------------------------------------*
 | output of solution vector to BINIO                  kremheller 04/18 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::MeshtyingStrategyArtery::Output()
{
  if (porofluidmultitimint_->Step() != 0) artnettimint_->Output(false, Teuchos::null);

  return;
}

/*----------------------------------------------------------------------*
 | evaluate matrix and rhs                             kremheller 04/18 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::MeshtyingStrategyArtery::Evaluate()
{
  // evaluate
  artnettimint_->AssembleMatAndRHS();
  // apply DBC
  artnettimint_->PrepareLinearSolve();

  arttoporofluidcoupling_->SetSolutionVectors(
      porofluidmultitimint_->Phinp(), porofluidmultitimint_->Phin(), artnettimint_->Pressurenp());

  // SetupCoupledArteryPoroFluidSystem();
  arttoporofluidcoupling_->SetupSystem(comb_systemmatrix_, rhs_,
      porofluidmultitimint_->SystemMatrix(), artnettimint_->SystemMatrix(),
      porofluidmultitimint_->RHS(), artnettimint_->RHS(),
      porofluidmultitimint_->GetDBCMapExtractor(), artnettimint_->GetDBCMapExtractor());

  return;
}

/*----------------------------------------------------------------------*
 | extract and update                                  kremheller 04/18 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector>
POROFLUIDMULTIPHASE::MeshtyingStrategyArtery::ExtractAndUpdateIter(
    const Teuchos::RCP<const Epetra_Vector> inc)
{
  Teuchos::RCP<const Epetra_Vector> arterypressinc;
  Teuchos::RCP<const Epetra_Vector> porofluidinc;

  arttoporofluidcoupling_->ExtractSingleFieldVectors(inc, porofluidinc, arterypressinc);

  artnettimint_->UpdateIter(arterypressinc);

  return porofluidinc;
}

/*----------------------------------------------------------------------*
 | artery dof row map                                  kremheller 04/18 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> POROFLUIDMULTIPHASE::MeshtyingStrategyArtery::ArteryDofRowMap() const
{
  return arttoporofluidcoupling_->ArteryDofRowMap();
}

/*-----------------------------------------------------------------------*
 | access to block system matrix of artery poro problem kremheller 04/18 |
 *-----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::BlockSparseMatrixBase>
POROFLUIDMULTIPHASE::MeshtyingStrategyArtery::ArteryPorofluidSysmat() const
{
  return comb_systemmatrix_;
}

/*----------------------------------------------------------------------*
 | return coupled residual                             kremheller 05/18 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> POROFLUIDMULTIPHASE::MeshtyingStrategyArtery::ArteryPorofluidRHS()
    const
{
  return rhs_;
}

/*----------------------------------------------------------------------*
 | extract and update                                  kremheller 04/18 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> POROFLUIDMULTIPHASE::MeshtyingStrategyArtery::CombinedIncrement(
    const Teuchos::RCP<const Epetra_Vector> inc) const
{
  return comb_increment_;
}

/*----------------------------------------------------------------------*
 | check initial fields                                kremheller 06/18 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::MeshtyingStrategyArtery::CheckInitialFields(
    Teuchos::RCP<const Epetra_Vector> vec_cont) const
{
  arttoporofluidcoupling_->CheckInitialFields(vec_cont, artnettimint_->Pressurenp());
  return;
}

/*-------------------------------------------------------------------------*
 | set element pairs that are close                       kremheller 03/19 |
 *------------------------------------------------------------------------ */
void POROFLUIDMULTIPHASE::MeshtyingStrategyArtery::SetNearbyElePairs(
    const std::map<int, std::set<int>>* nearbyelepairs)
{
  arttoporofluidcoupling_->SetNearbyElePairs(nearbyelepairs);
  return;
}

/*-------------------------------------------------------------------------*
 | setup the strategy                                     kremheller 03/19 |
 *------------------------------------------------------------------------ */
void POROFLUIDMULTIPHASE::MeshtyingStrategyArtery::Setup()
{
  arttoporofluidcoupling_->Setup();
  return;
}

/*----------------------------------------------------------------------*
 | apply mesh movement                                 kremheller 06/18 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::MeshtyingStrategyArtery::ApplyMeshMovement() const
{
  arttoporofluidcoupling_->ApplyMeshMovement();
  return;
}

/*----------------------------------------------------------------------*
 | access to blood vessel volume fraction              kremheller 10/19 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector>
POROFLUIDMULTIPHASE::MeshtyingStrategyArtery::BloodVesselVolumeFraction()
{
  return arttoporofluidcoupling_->BloodVesselVolumeFraction();
}
