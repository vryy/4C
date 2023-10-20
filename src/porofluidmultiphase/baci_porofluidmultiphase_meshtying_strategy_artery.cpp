/*----------------------------------------------------------------------*/
/*! \file

\brief routines for coupling with artery network

\level 3

*----------------------------------------------------------------------*/

#include "baci_porofluidmultiphase_meshtying_strategy_artery.H"

#include "baci_adapter_art_net.H"
#include "baci_art_net_utils.H"
#include "baci_inpar_bio.H"
#include "baci_io.H"
#include "baci_io_control.H"
#include "baci_lib_globalproblem.H"
#include "baci_linalg_utils_sparse_algebra_print.H"
#include "baci_linear_solver_method_linalg.H"
#include "baci_porofluidmultiphase_utils.H"
#include "baci_poromultiphase_scatra_artery_coupling_base.H"
#include "baci_poromultiphase_scatra_utils.H"

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

  const bool evaluate_on_lateral_surface = DRT::INPUT::IntegralValue<int>(
      poroparams.sublist("ARTERY COUPLING"), "LATERAL_SURFACE_COUPLING");

  const std::string couplingcondname = std::invoke(
      [&]()
      {
        if (DRT::INPUT::IntegralValue<INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod>(
                DRT::Problem::Instance()->PoroFluidMultiPhaseDynamicParams().sublist(
                    "ARTERY COUPLING"),
                "ARTERY_COUPLING_METHOD") ==
            INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod::ntp)
        {
          return "ArtPorofluidCouplConNodeToPoint";
        }
        else
        {
          return "ArtPorofluidCouplConNodebased";
        }
      });

  // initialize mesh tying object
  arttoporofluidcoupling_ = POROMULTIPHASESCATRA::UTILS::CreateAndInitArteryCouplingStrategy(
      arterydis_, porofluidmultitimint->Discretization(), poroparams.sublist("ARTERY COUPLING"),
      couplingcondname, "COUPLEDDOFS_ART", "COUPLEDDOFS_PORO", evaluate_on_lateral_surface);

  // Initialize rhs vector
  rhs_ = Teuchos::rcp(new Epetra_Vector(*arttoporofluidcoupling_->FullMap(), true));

  // Initialize increment vector
  comb_increment_ = Teuchos::rcp(new Epetra_Vector(*arttoporofluidcoupling_->FullMap(), true));
  // Initialize phinp vector
  comb_phinp_ = Teuchos::rcp(new Epetra_Vector(*arttoporofluidcoupling_->FullMap(), true));

  // initialize Poromultiphase-elasticity-systemmatrix_
  comb_systemmatrix_ =
      Teuchos::rcp(new CORE::LINALG::BlockSparseMatrix<CORE::LINALG::DefaultBlockMatrixStrategy>(
          *arttoporofluidcoupling_->GlobalExtractor(), *arttoporofluidcoupling_->GlobalExtractor(),
          81, false, true));

  return;
}



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
    Teuchos::RCP<CORE::LINALG::Solver> solver)
{
  const Teuchos::ParameterList& porofluidparams =
      DRT::Problem::Instance()->PoroFluidMultiPhaseDynamicParams();
  const int linsolvernumber = porofluidparams.get<int>("LINEAR_SOLVER");
  const Teuchos::ParameterList& solverparams =
      DRT::Problem::Instance()->SolverParams(linsolvernumber);
  const auto solvertype =
      Teuchos::getIntegralValue<INPAR::SOLVER::SolverType>(solverparams, "SOLVER");
  // no need to do the rest for direct solvers
  if (solvertype == INPAR::SOLVER::SolverType::umfpack or
      solvertype == INPAR::SOLVER::SolverType::superlu)
    return;

  if (solvertype != INPAR::SOLVER::SolverType::belos) dserror("Iterative solver expected");

  const auto azprectype =
      Teuchos::getIntegralValue<INPAR::SOLVER::PreconditionerType>(solverparams, "AZPREC");

  // plausibility check
  switch (azprectype)
  {
    case INPAR::SOLVER::PreconditionerType::multigrid_nxn:
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
  blocksmootherparams1.sublist("Belos Parameters");
  blocksmootherparams1.sublist("MueLu Parameters");

  porofluidmultitimint_->Discretization()->ComputeNullSpaceIfNecessary(blocksmootherparams1);

  Teuchos::ParameterList& blocksmootherparams2 = solver->Params().sublist("Inverse2");
  blocksmootherparams2.sublist("Belos Parameters");
  blocksmootherparams2.sublist("MueLu Parameters");

  arterydis_->ComputeNullSpaceIfNecessary(blocksmootherparams2);
}

/*--------------------------------------------------------------------------*
 | solve linear system of equations                        kremheller 04/18 |
 *--------------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::MeshtyingStrategyArtery::LinearSolve(
    Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<CORE::LINALG::SparseOperator> sysmat,
    Teuchos::RCP<Epetra_Vector> increment, Teuchos::RCP<Epetra_Vector> residual)
{
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
void POROFLUIDMULTIPHASE::MeshtyingStrategyArtery::CalculateNorms(std::vector<double>& preresnorm,
    std::vector<double>& incprenorm, std::vector<double>& prenorm,
    const Teuchos::RCP<const Epetra_Vector> increment)
{
  preresnorm.resize(2);
  incprenorm.resize(2);
  prenorm.resize(2);

  prenorm[0] = UTILS::CalculateVectorNorm(vectornorminc_, porofluidmultitimint_->Phinp());
  prenorm[1] = UTILS::CalculateVectorNorm(vectornorminc_, artnettimint_->Pressurenp());

  Teuchos::RCP<const Epetra_Vector> arterypressinc;
  Teuchos::RCP<const Epetra_Vector> porofluidinc;

  arttoporofluidcoupling_->ExtractSingleFieldVectors(comb_increment_, porofluidinc, arterypressinc);

  incprenorm[0] = UTILS::CalculateVectorNorm(vectornorminc_, porofluidinc);
  incprenorm[1] = UTILS::CalculateVectorNorm(vectornorminc_, arterypressinc);

  Teuchos::RCP<const Epetra_Vector> arterypressrhs;
  Teuchos::RCP<const Epetra_Vector> porofluidrhs;

  arttoporofluidcoupling_->ExtractSingleFieldVectors(rhs_, porofluidrhs, arterypressrhs);

  preresnorm[0] = UTILS::CalculateVectorNorm(vectornormfres_, porofluidrhs);
  preresnorm[1] = UTILS::CalculateVectorNorm(vectornormfres_, arterypressrhs);

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
  arttoporofluidcoupling_->SetSolutionVectors(
      porofluidmultitimint_->Phinp(), porofluidmultitimint_->Phin(), artnettimint_->Pressurenp());

  // evaluate the coupling
  arttoporofluidcoupling_->Evaluate(comb_systemmatrix_, rhs_);

  // evaluate artery
  artnettimint_->AssembleMatAndRHS();
  // apply DBC
  artnettimint_->PrepareLinearSolve();

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
Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase>
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
