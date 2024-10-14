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
#include "4C_linear_solver_method_parameters.hpp"
#include "4C_poromultiphase_scatra_artery_coupling_nodebased.hpp"
#include "4C_poromultiphase_scatra_utils.hpp"
#include "4C_scatra_timint_implicit.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

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
void ScaTra::MeshtyingStrategyArtery::init_meshtying()
{
  // instantiate strategy for Newton-Raphson convergence check
  init_conv_check_strategy();

  const Teuchos::ParameterList& globaltimeparams =
      Global::Problem::instance()->poro_multi_phase_scatra_dynamic_params();
  const Teuchos::ParameterList& myscatraparams =
      Global::Problem::instance()->scalar_transport_dynamic_params();
  if (Teuchos::getIntegralValue<Inpar::ScaTra::VelocityField>(myscatraparams, "VELOCITYFIELD") !=
      Inpar::ScaTra::velocity_zero)
    FOUR_C_THROW("set your velocity field to zero!");

  // construct artery scatra problem
  Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> art_scatra =
      Teuchos::make_rcp<Adapter::ScaTraBaseAlgorithm>(globaltimeparams, myscatraparams,
          Global::Problem::instance()->solver_params(myscatraparams.get<int>("LINEAR_SOLVER")),
          "artery_scatra", false);

  // initialize the base algo.
  // scatra time integrator is initialized inside.
  art_scatra->init();

  // only now we must call setup() on the scatra time integrator.
  // all objects relying on the parallel distribution are
  // created and pointers are set.
  // calls setup() on the scatra time integrator inside.
  art_scatra->scatra_field()->setup();
  Global::Problem::instance()->add_field_test(art_scatra->create_scatra_field_test());

  // set the time integrator
  set_artery_scatra_time_integrator(art_scatra->scatra_field());

  // get the two discretizations
  artscatradis_ = artscatratimint_->discretization();
  scatradis_ = scatratimint_->discretization();

  if (scatratimint_->discretization()->get_comm().MyPID() == 0)
  {
    std::cout << "\n";
    std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    std::cout << "<                                                  >" << std::endl;
    std::cout << "< ScaTra-Coupling with 1D Artery Network activated >" << std::endl;
  }

  const bool evaluate_on_lateral_surface = Global::Problem::instance()
                                               ->poro_fluid_multi_phase_dynamic_params()
                                               .sublist("ARTERY COUPLING")
                                               .get<bool>("LATERAL_SURFACE_COUPLING");

  // set coupling condition name
  const std::string couplingcondname = std::invoke(
      [&]()
      {
        if (Teuchos::getIntegralValue<
                Inpar::ArteryNetwork::ArteryPoroMultiphaseScatraCouplingMethod>(
                Global::Problem::instance()->poro_fluid_multi_phase_dynamic_params().sublist(
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
  arttoscatracoupling_ = PoroMultiPhaseScaTra::Utils::create_and_init_artery_coupling_strategy(
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
  rhs_ = Teuchos::make_rcp<Core::LinAlg::Vector<double>>(*arttoscatracoupling_->full_map(), true);

  // Initialize increment vector
  comb_increment_ =
      Teuchos::make_rcp<Core::LinAlg::Vector<double>>(*arttoscatracoupling_->full_map(), true);

  // initialize scatra-artery_scatra-systemmatrix_
  comb_systemmatrix_ =
      Teuchos::make_rcp<Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>>(
          *arttoscatracoupling_->global_extractor(), *arttoscatracoupling_->global_extractor(), 81,
          false, true);

  arttoscatracoupling_->setup();

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
      Global::Problem::instance()->solver_params(linsolvernumber);
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
    case Core::LinearSolver::PreconditionerType::block_teko:
    {
      // no plausibility checks here
      // if you forget to declare an xml file you will get an error message anyway
    }
    break;
    default:
      FOUR_C_THROW("AMGnxn preconditioner expected");
      break;
  }

  Teuchos::ParameterList& blocksmootherparams1 = solver().params().sublist("Inverse1");
  Core::LinearSolver::Parameters::compute_solver_parameters(*scatradis_, blocksmootherparams1);

  Teuchos::ParameterList& blocksmootherparams2 = solver().params().sublist("Inverse2");
  Core::LinearSolver::Parameters::compute_solver_parameters(*artscatradis_, blocksmootherparams2);
}

/*-----------------------------------------------------------------------*
 | return global map of degrees of freedom              kremheller 04/18 |
 *-----------------------------------------------------------------------*/
const Epetra_Map& ScaTra::MeshtyingStrategyArtery::dof_row_map() const
{
  return *arttoscatracoupling_->full_map();
}

/*-----------------------------------------------------------------------*
 | return global map of degrees of freedom              kremheller 04/18 |
 *-----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ScaTra::MeshtyingStrategyArtery::art_scatra_dof_row_map() const
{
  return arttoscatracoupling_->artery_dof_row_map();
}

/*-------------------------------------------------------------------------------*
 | return linear solver for global system of linear equations   kremheller 04/18 |
 *-------------------------------------------------------------------------------*/
const Core::LinAlg::Solver& ScaTra::MeshtyingStrategyArtery::solver() const
{
  if (scatratimint_->solver() == Teuchos::null) FOUR_C_THROW("Invalid linear solver!");

  return *scatratimint_->solver();
}

/*------------------------------------------------------------------------------*
 | instantiate strategy for Newton-Raphson convergence check   kremheller 04/18 |
 *------------------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyArtery::init_conv_check_strategy()
{
  convcheckstrategy_ = Teuchos::make_rcp<ScaTra::ConvCheckStrategyPoroMultiphaseScatraArtMeshTying>(
      scatratimint_->scatra_parameter_list()->sublist("NONLINEAR"));

  return;
}

/*------------------------------------------------------------------------------------------*
 | solve linear system of equations for scatra-scatra interface coupling   kremheller 04/18 |
 *------------------------------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyArtery::solve(
    const Teuchos::RCP<Core::LinAlg::Solver>& solver,                //!< solver
    const Teuchos::RCP<Core::LinAlg::SparseOperator>& systemmatrix,  //!< system matrix
    const Teuchos::RCP<Core::LinAlg::Vector<double>>& increment,     //!< increment vector
    const Teuchos::RCP<Core::LinAlg::Vector<double>>& residual,      //!< residual vector
    const Teuchos::RCP<Core::LinAlg::Vector<double>>& phinp,         //!< state vector at time n+1
    const int iteration,  //!< number of current Newton-Raphson iteration
    Core::LinAlg::SolverParams& solver_params) const
{
  // setup the system (evaluate mesh tying)
  // reason for this being done here is that we need the system matrix of the continuous scatra
  // problem with DBCs applied which is performed directly before calling solve

  setup_system(systemmatrix, residual);

  comb_systemmatrix_->complete();

  // solve
  comb_increment_->PutScalar(0.0);
  solver_params.refactor = true;
  solver_params.reset = iteration == 1;
  solver->solve(comb_systemmatrix_->epetra_operator(), comb_increment_, rhs_, solver_params);

  // extract increments of scatra and artery-scatra field
  Teuchos::RCP<const Core::LinAlg::Vector<double>> artscatrainc;
  Teuchos::RCP<const Core::LinAlg::Vector<double>> myinc;
  extract_single_field_vectors(comb_increment_, myinc, artscatrainc);

  // update the scatra increment, update iter is performed outside
  increment->Update(1.0, *(myinc), 1.0);
  // update the artery-scatra field
  artscatratimint_->update_iter(*artscatrainc);

  return;
}

/*------------------------------------------------------------------------------------------*
 | solve linear system of equations for scatra-scatra interface coupling   kremheller 04/18 |
 *------------------------------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyArtery::setup_system(
    const Teuchos::RCP<Core::LinAlg::SparseOperator>& systemmatrix,  //!< system matrix
    const Teuchos::RCP<Core::LinAlg::Vector<double>>& residual       //!< residual vector
) const
{
  arttoscatracoupling_->set_solution_vectors(
      scatratimint_->phinp(), Teuchos::null, artscatratimint_->phinp());

  // evaluate the 1D-3D coupling
  arttoscatracoupling_->evaluate(comb_systemmatrix_, rhs_);

  // evaluate 1D sub-problem
  artscatratimint_->prepare_linear_solve();

  // setup the entire system
  arttoscatracoupling_->setup_system(comb_systemmatrix_, rhs_,
      Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(systemmatrix),
      artscatratimint_->system_matrix(), residual, artscatratimint_->residual(),
      scatratimint_->dirich_maps(), artscatratimint_->dirich_maps());
}

/*-------------------------------------------------------------------------*
 | set time integrator for scalar transport in arteries   kremheller 04/18 |
 *------------------------------------------------------------------------ */
void ScaTra::MeshtyingStrategyArtery::update_art_scatra_iter(
    Teuchos::RCP<const Core::LinAlg::Vector<double>> combined_inc)
{
  Teuchos::RCP<const Core::LinAlg::Vector<double>> artscatrainc;
  Teuchos::RCP<const Core::LinAlg::Vector<double>> myinc;
  extract_single_field_vectors(combined_inc, myinc, artscatrainc);

  artscatratimint_->update_iter(*artscatrainc);

  return;
}

/*-------------------------------------------------------------------------*
 | extract single field vectors                           kremheller 10/20 |
 *------------------------------------------------------------------------ */
void ScaTra::MeshtyingStrategyArtery::extract_single_field_vectors(
    Teuchos::RCP<const Core::LinAlg::Vector<double>> globalvec,
    Teuchos::RCP<const Core::LinAlg::Vector<double>>& vec_cont,
    Teuchos::RCP<const Core::LinAlg::Vector<double>>& vec_art) const
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
void ScaTra::MeshtyingStrategyArtery::set_nearby_ele_pairs(
    const std::map<int, std::set<int>>* nearbyelepairs)
{
  arttoscatracoupling_->set_nearby_ele_pairs(nearbyelepairs);

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
void ScaTra::MeshtyingStrategyArtery::set_artery_pressure() const
{
  artscatradis_->set_state(2, "one_d_artery_pressure", arttimint_->pressurenp());
  return;
}

/*--------------------------------------------------------------------------*
 | apply mesh movement on artery coupling                  kremheller 07/18 |
 *--------------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyArtery::apply_mesh_movement()
{
  arttoscatracoupling_->apply_mesh_movement();
  return;
}

/*--------------------------------------------------------------------------*
 | check if initial fields match                           kremheller 04/18 |
 *--------------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyArtery::check_initial_fields() const
{
  arttoscatracoupling_->check_initial_fields(scatratimint_->phinp(), artscatratimint_->phinp());

  return;
}

FOUR_C_NAMESPACE_CLOSE
