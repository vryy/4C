/*----------------------------------------------------------------------*/
/*! \file
 \brief routines for coupling between 1D arterial network and 2D/3D
        scatra-algorithm

   \level 3

 *----------------------------------------------------------------------*/

#include "4C_scatra_timint_meshtying_strategy_artery.hpp"

#include "4C_adapter_art_net.hpp"
#include "4C_adapter_scatra_base_algorithm.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_bio.hpp"
#include "4C_linear_solver_method.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_poromultiphase_scatra_artery_coupling_nodebased.hpp"
#include "4C_poromultiphase_scatra_utils.hpp"
#include "4C_scatra_timint_implicit.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | constructor                                         kremheller 04/18 |
 *----------------------------------------------------------------------*/
ScaTra::MeshtyingStrategyArtery::MeshtyingStrategyArtery(
    ScaTra::ScaTraTimIntImpl* scatratimint  //!< scalar transport time integrator
    )
    : MeshtyingStrategyBase(scatratimint)
{
}

/*----------------------------------------------------------------------*
 | init                                                kremheller 04/18 |
 *----------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyArtery::InitMeshtying()
{
  // instantiate strategy for Newton-Raphson convergence check
  init_conv_check_strategy();

  const Teuchos::ParameterList& globaltimeparams =
      Global::Problem::Instance()->poro_multi_phase_scatra_dynamic_params();
  const Teuchos::ParameterList& myscatraparams =
      Global::Problem::Instance()->scalar_transport_dynamic_params();
  if (Core::UTILS::IntegralValue<Inpar::ScaTra::VelocityField>(myscatraparams, "VELOCITYFIELD") !=
      Inpar::ScaTra::velocity_zero)
    FOUR_C_THROW("set your velocity field to zero!");

  // construct artery scatra problem
  Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> art_scatra =
      Teuchos::rcp(new Adapter::ScaTraBaseAlgorithm(globaltimeparams, myscatraparams,
          Global::Problem::Instance()->SolverParams(myscatraparams.get<int>("LINEAR_SOLVER")),
          "artery_scatra", false));

  // initialize the base algo.
  // scatra time integrator is initialized inside.
  art_scatra->Init();

  // only now we must call Setup() on the scatra time integrator.
  // all objects relying on the parallel distribution are
  // created and pointers are set.
  // calls Setup() on the scatra time integrator inside.
  art_scatra->ScaTraField()->Setup();
  Global::Problem::Instance()->AddFieldTest(art_scatra->create_sca_tra_field_test());

  // set the time integrator
  set_artery_scatra_time_integrator(art_scatra->ScaTraField());

  // get the two discretizations
  artscatradis_ = artscatratimint_->discretization();
  scatradis_ = scatratimint_->discretization();

  if (scatratimint_->discretization()->Comm().MyPID() == 0)
  {
    std::cout << "\n";
    std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    std::cout << "<                                                  >" << std::endl;
    std::cout << "< ScaTra-Coupling with 1D Artery Network activated >" << std::endl;
  }

  const bool evaluate_on_lateral_surface = Core::UTILS::IntegralValue<int>(
      Global::Problem::Instance()->poro_fluid_multi_phase_dynamic_params().sublist(
          "ARTERY COUPLING"),
      "LATERAL_SURFACE_COUPLING");

  // set coupling condition name
  const std::string couplingcondname = std::invoke(
      [&]()
      {
        if (Core::UTILS::IntegralValue<
                Inpar::ArteryNetwork::ArteryPoroMultiphaseScatraCouplingMethod>(
                Global::Problem::Instance()->poro_fluid_multi_phase_dynamic_params().sublist(
                    "ARTERY COUPLING"),
                "ARTERY_COUPLING_METHOD") ==
            Inpar::ArteryNetwork::ArteryPoroMultiphaseScatraCouplingMethod::ntp)
        {
          return "ArtScatraCouplConNodeToPoint";
        }
        else
        {
          return "ArtScatraCouplConNodebased";
        }
      });

  // init the mesh tying object, which does all the work
  arttoscatracoupling_ = PoroMultiPhaseScaTra::UTILS::CreateAndInitArteryCouplingStrategy(
      artscatradis_, scatradis_, myscatraparams.sublist("ARTERY COUPLING"), couplingcondname,
      "COUPLEDDOFS_ARTSCATRA", "COUPLEDDOFS_SCATRA", evaluate_on_lateral_surface);

  initialize_linear_solver(myscatraparams);
}

/*----------------------------------------------------------------------*
 | setup                                               kremheller 04/18 |
 *----------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyArtery::setup_meshtying()
{
  // Initialize rhs vector
  rhs_ = Teuchos::rcp(new Epetra_Vector(*arttoscatracoupling_->FullMap(), true));

  // Initialize increment vector
  comb_increment_ = Teuchos::rcp(new Epetra_Vector(*arttoscatracoupling_->FullMap(), true));

  // initialize scatra-artery_scatra-systemmatrix_
  comb_systemmatrix_ =
      Teuchos::rcp(new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
          *arttoscatracoupling_->GlobalExtractor(), *arttoscatracoupling_->GlobalExtractor(), 81,
          false, true));

  arttoscatracoupling_->Setup();

  return;
}

/*----------------------------------------------------------------------*
 | initialize the linear solver                        kremheller 07/20 |
 *----------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyArtery::initialize_linear_solver(
    const Teuchos::ParameterList& scatraparams)
{
  const int linsolvernumber = scatraparams.get<int>("LINEAR_SOLVER");
  const Teuchos::ParameterList& solverparams =
      Global::Problem::Instance()->SolverParams(linsolvernumber);
  const auto solvertype =
      Teuchos::getIntegralValue<Core::LinearSolver::SolverType>(solverparams, "SOLVER");
  // no need to do the rest for direct solvers
  if (solvertype == Core::LinearSolver::SolverType::umfpack or
      solvertype == Core::LinearSolver::SolverType::superlu)
    return;

  if (solvertype != Core::LinearSolver::SolverType::belos)
    FOUR_C_THROW("Iterative solver expected");

  const auto azprectype =
      Teuchos::getIntegralValue<Core::LinearSolver::PreconditionerType>(solverparams, "AZPREC");

  // plausibility check
  switch (azprectype)
  {
    case Core::LinearSolver::PreconditionerType::multigrid_nxn:
    {
      // no plausibility checks here
      // if you forget to declare an xml file you will get an error message anyway
    }
    break;
    default:
      FOUR_C_THROW("AMGnxn preconditioner expected");
      break;
  }

  // equip smoother for fluid matrix block with empty parameter sublists to trigger null space
  // computation
  Teuchos::ParameterList& blocksmootherparams1 = Solver().Params().sublist("Inverse1");
  blocksmootherparams1.sublist("Belos Parameters");
  blocksmootherparams1.sublist("MueLu Parameters");

  scatradis_->compute_null_space_if_necessary(blocksmootherparams1);

  Teuchos::ParameterList& blocksmootherparams2 = Solver().Params().sublist("Inverse2");
  blocksmootherparams2.sublist("Belos Parameters");
  blocksmootherparams2.sublist("MueLu Parameters");

  artscatradis_->compute_null_space_if_necessary(blocksmootherparams2);
}

/*-----------------------------------------------------------------------*
 | return global map of degrees of freedom              kremheller 04/18 |
 *-----------------------------------------------------------------------*/
const Epetra_Map& ScaTra::MeshtyingStrategyArtery::dof_row_map() const
{
  return *arttoscatracoupling_->FullMap();
}

/*-----------------------------------------------------------------------*
 | return global map of degrees of freedom              kremheller 04/18 |
 *-----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ScaTra::MeshtyingStrategyArtery::ArtScatraDofRowMap() const
{
  return arttoscatracoupling_->ArteryDofRowMap();
}

/*-------------------------------------------------------------------------------*
 | return linear solver for global system of linear equations   kremheller 04/18 |
 *-------------------------------------------------------------------------------*/
const Core::LinAlg::Solver& ScaTra::MeshtyingStrategyArtery::Solver() const
{
  if (scatratimint_->Solver() == Teuchos::null) FOUR_C_THROW("Invalid linear solver!");

  return *scatratimint_->Solver();
}

/*------------------------------------------------------------------------------*
 | instantiate strategy for Newton-Raphson convergence check   kremheller 04/18 |
 *------------------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyArtery::init_conv_check_strategy()
{
  convcheckstrategy_ = Teuchos::rcp(new ScaTra::ConvCheckStrategyPoroMultiphaseScatraArtMeshTying(
      scatratimint_->ScatraParameterList()->sublist("NONLINEAR")));

  return;
}

/*------------------------------------------------------------------------------------------*
 | solve linear system of equations for scatra-scatra interface coupling   kremheller 04/18 |
 *------------------------------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyArtery::Solve(
    const Teuchos::RCP<Core::LinAlg::Solver>& solver,                //!< solver
    const Teuchos::RCP<Core::LinAlg::SparseOperator>& systemmatrix,  //!< system matrix
    const Teuchos::RCP<Epetra_Vector>& increment,                    //!< increment vector
    const Teuchos::RCP<Epetra_Vector>& residual,                     //!< residual vector
    const Teuchos::RCP<Epetra_Vector>& phinp,                        //!< state vector at time n+1
    const int iteration,  //!< number of current Newton-Raphson iteration
    Core::LinAlg::SolverParams& solver_params) const
{
  // setup the system (evaluate mesh tying)
  // reason for this being done here is that we need the system matrix of the continuous scatra
  // problem with DBCs applied which is performed directly before calling solve

  SetupSystem(systemmatrix, residual);

  comb_systemmatrix_->Complete();

  // solve
  comb_increment_->PutScalar(0.0);
  solver_params.refactor = true;
  solver_params.reset = iteration == 1;
  solver->Solve(comb_systemmatrix_->EpetraOperator(), comb_increment_, rhs_, solver_params);

  // extract increments of scatra and artery-scatra field
  Teuchos::RCP<const Epetra_Vector> artscatrainc;
  Teuchos::RCP<const Epetra_Vector> myinc;
  extract_single_field_vectors(comb_increment_, myinc, artscatrainc);

  // update the scatra increment, update iter is performed outside
  increment->Update(1.0, *(myinc), 1.0);
  // update the artery-scatra field
  artscatratimint_->UpdateIter(artscatrainc);

  return;
}

/*------------------------------------------------------------------------------------------*
 | solve linear system of equations for scatra-scatra interface coupling   kremheller 04/18 |
 *------------------------------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyArtery::SetupSystem(
    const Teuchos::RCP<Core::LinAlg::SparseOperator>& systemmatrix,  //!< system matrix
    const Teuchos::RCP<Epetra_Vector>& residual                      //!< residual vector
) const
{
  arttoscatracoupling_->SetSolutionVectors(
      scatratimint_->Phinp(), Teuchos::null, artscatratimint_->Phinp());

  // evaluate the 1D-3D coupling
  arttoscatracoupling_->Evaluate(comb_systemmatrix_, rhs_);

  // evaluate 1D sub-problem
  artscatratimint_->PrepareLinearSolve();

  // setup the entire system
  arttoscatracoupling_->SetupSystem(comb_systemmatrix_, rhs_,
      Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(systemmatrix),
      artscatratimint_->SystemMatrix(), residual, artscatratimint_->Residual(),
      scatratimint_->DirichMaps(), artscatratimint_->DirichMaps());
}

/*-------------------------------------------------------------------------*
 | set time integrator for scalar transport in arteries   kremheller 04/18 |
 *------------------------------------------------------------------------ */
void ScaTra::MeshtyingStrategyArtery::UpdateArtScatraIter(
    Teuchos::RCP<const Epetra_Vector> combined_inc)
{
  Teuchos::RCP<const Epetra_Vector> artscatrainc;
  Teuchos::RCP<const Epetra_Vector> myinc;
  extract_single_field_vectors(combined_inc, myinc, artscatrainc);

  artscatratimint_->UpdateIter(artscatrainc);

  return;
}

/*-------------------------------------------------------------------------*
 | extract single field vectors                           kremheller 10/20 |
 *------------------------------------------------------------------------ */
void ScaTra::MeshtyingStrategyArtery::extract_single_field_vectors(
    Teuchos::RCP<const Epetra_Vector> globalvec, Teuchos::RCP<const Epetra_Vector>& vec_cont,
    Teuchos::RCP<const Epetra_Vector>& vec_art) const
{
  arttoscatracoupling_->extract_single_field_vectors(globalvec, vec_cont, vec_art);

  return;
}

/*-------------------------------------------------------------------------*
 | set time integrator for scalar transport in arteries   kremheller 04/18 |
 *------------------------------------------------------------------------ */
void ScaTra::MeshtyingStrategyArtery::set_artery_scatra_time_integrator(
    Teuchos::RCP<ScaTra::ScaTraTimIntImpl> artscatratimint)
{
  artscatratimint_ = artscatratimint;
  if (artscatratimint_ == Teuchos::null)
    FOUR_C_THROW("could not set artery scatra time integrator");

  return;
}

/*-------------------------------------------------------------------------*
 | set time integrator for artery problems                kremheller 04/18 |
 *------------------------------------------------------------------------ */
void ScaTra::MeshtyingStrategyArtery::set_artery_time_integrator(
    Teuchos::RCP<Adapter::ArtNet> arttimint)
{
  arttimint_ = arttimint;
  if (arttimint_ == Teuchos::null) FOUR_C_THROW("could not set artery time integrator");

  return;
}

/*-------------------------------------------------------------------------*
 | set element pairs that are close                       kremheller 03/19 |
 *------------------------------------------------------------------------ */
void ScaTra::MeshtyingStrategyArtery::SetNearbyElePairs(
    const std::map<int, std::set<int>>* nearbyelepairs)
{
  arttoscatracoupling_->SetNearbyElePairs(nearbyelepairs);

  return;
}

/*--------------------------------------------------------------------------*
 | setup the coupled matrix                                kremheller 04/18 |
 *--------------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyArtery::prepare_time_step() const
{
  artscatratimint_->prepare_time_step();
  return;
}

/*--------------------------------------------------------------------------*
 | setup the coupled matrix                                kremheller 04/18 |
 *--------------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyArtery::SetArteryPressure() const
{
  artscatradis_->set_state(2, "one_d_artery_pressure", arttimint_->Pressurenp());
  return;
}

/*--------------------------------------------------------------------------*
 | apply mesh movement on artery coupling                  kremheller 07/18 |
 *--------------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyArtery::ApplyMeshMovement()
{
  arttoscatracoupling_->ApplyMeshMovement();
  return;
}

/*--------------------------------------------------------------------------*
 | check if initial fields match                           kremheller 04/18 |
 *--------------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyArtery::CheckInitialFields() const
{
  arttoscatracoupling_->CheckInitialFields(scatratimint_->Phinp(), artscatratimint_->Phinp());

  return;
}

FOUR_C_NAMESPACE_CLOSE
