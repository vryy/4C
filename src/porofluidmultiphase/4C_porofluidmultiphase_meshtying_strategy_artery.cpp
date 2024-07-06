/*----------------------------------------------------------------------*/
/*! \file

\brief routines for coupling with artery network

\level 3

*----------------------------------------------------------------------*/

#include "4C_porofluidmultiphase_meshtying_strategy_artery.hpp"

#include "4C_adapter_art_net.hpp"
#include "4C_art_net_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_bio.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_utils_sparse_algebra_print.hpp"
#include "4C_linear_solver_method.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_porofluidmultiphase_utils.hpp"
#include "4C_poromultiphase_scatra_artery_coupling_base.hpp"
#include "4C_poromultiphase_scatra_utils.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 | constructor                                (public) kremheller 04/18 |
 *----------------------------------------------------------------------*/
POROFLUIDMULTIPHASE::MeshtyingStrategyArtery::MeshtyingStrategyArtery(
    POROFLUIDMULTIPHASE::TimIntImpl* porofluidmultitimint, const Teuchos::ParameterList& probparams,
    const Teuchos::ParameterList& poroparams)
    : MeshtyingStrategyBase(porofluidmultitimint, probparams, poroparams)
{
  const Teuchos::ParameterList& artdyn = Global::Problem::instance()->arterial_dynamic_params();

  arterydis_ = Global::Problem::instance()->get_dis("artery");

  if (!arterydis_->filled()) arterydis_->fill_complete();

  Inpar::ArtDyn::TimeIntegrationScheme timintscheme =
      Core::UTILS::IntegralValue<Inpar::ArtDyn::TimeIntegrationScheme>(artdyn, "DYNAMICTYP");

  Teuchos::RCP<Core::IO::DiscretizationWriter> artery_output = arterydis_->writer();
  artery_output->write_mesh(0, 0.0);

  // build art net time integrator
  artnettimint_ = Arteries::UTILS::CreateAlgorithm(timintscheme, arterydis_,
      artdyn.get<int>("LINEAR_SOLVER"), probparams, artdyn, artery_output);

  // set to false
  artnettimint_->set_solve_scatra(false);

  // initialize
  artnettimint_->init(probparams, artdyn, "artery_scatra");

  // print user info
  if (porofluidmultitimint->discretization()->get_comm().MyPID() == 0)
  {
    std::cout << "\n";
    std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    std::cout << "<                                                  >" << std::endl;
    std::cout << "<    Coupling with 1D Artery Network activated     >" << std::endl;
  }

  const bool evaluate_on_lateral_surface = Core::UTILS::IntegralValue<int>(
      poroparams.sublist("ARTERY COUPLING"), "LATERAL_SURFACE_COUPLING");

  const std::string couplingcondname = std::invoke(
      [&]()
      {
        if (Core::UTILS::IntegralValue<
                Inpar::ArteryNetwork::ArteryPoroMultiphaseScatraCouplingMethod>(
                Global::Problem::instance()->poro_fluid_multi_phase_dynamic_params().sublist(
                    "ARTERY COUPLING"),
                "ARTERY_COUPLING_METHOD") ==
            Inpar::ArteryNetwork::ArteryPoroMultiphaseScatraCouplingMethod::ntp)
        {
          return "ArtPorofluidCouplConNodeToPoint";
        }
        else
        {
          return "ArtPorofluidCouplConNodebased";
        }
      });

  // initialize mesh tying object
  arttoporofluidcoupling_ = PoroMultiPhaseScaTra::UTILS::CreateAndInitArteryCouplingStrategy(
      arterydis_, porofluidmultitimint->discretization(), poroparams.sublist("ARTERY COUPLING"),
      couplingcondname, "COUPLEDDOFS_ART", "COUPLEDDOFS_PORO", evaluate_on_lateral_surface);

  // Initialize rhs vector
  rhs_ = Teuchos::rcp(new Epetra_Vector(*arttoporofluidcoupling_->full_map(), true));

  // Initialize increment vector
  comb_increment_ = Teuchos::rcp(new Epetra_Vector(*arttoporofluidcoupling_->full_map(), true));
  // Initialize phinp vector
  comb_phinp_ = Teuchos::rcp(new Epetra_Vector(*arttoporofluidcoupling_->full_map(), true));

  // initialize Poromultiphase-elasticity-systemmatrix_
  comb_systemmatrix_ =
      Teuchos::rcp(new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
          *arttoporofluidcoupling_->global_extractor(),
          *arttoporofluidcoupling_->global_extractor(), 81, false, true));

  return;
}



/*----------------------------------------------------------------------*
 | prepare time loop                                   kremheller 04/18 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::MeshtyingStrategyArtery::prepare_time_loop()
{
  artnettimint_->prepare_time_loop();
  return;
}

/*----------------------------------------------------------------------*
 | setup the variables to do a new time step  (public) kremheller 04/18 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::MeshtyingStrategyArtery::prepare_time_step()
{
  artnettimint_->prepare_time_step();
  return;
}

/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                     kremheller 04/18 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::MeshtyingStrategyArtery::update()
{
  artnettimint_->time_update();
  return;
}

/*--------------------------------------------------------------------------*
 | initialize the linear solver                            kremheller 07/20 |
 *--------------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::MeshtyingStrategyArtery::initialize_linear_solver(
    Teuchos::RCP<Core::LinAlg::Solver> solver)
{
  const Teuchos::ParameterList& porofluidparams =
      Global::Problem::instance()->poro_fluid_multi_phase_dynamic_params();
  const int linsolvernumber = porofluidparams.get<int>("LINEAR_SOLVER");
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
  Teuchos::ParameterList& blocksmootherparams1 = solver->params().sublist("Inverse1");
  blocksmootherparams1.sublist("Belos Parameters");
  blocksmootherparams1.sublist("MueLu Parameters");

  porofluidmultitimint_->discretization()->compute_null_space_if_necessary(blocksmootherparams1);

  Teuchos::ParameterList& blocksmootherparams2 = solver->params().sublist("Inverse2");
  blocksmootherparams2.sublist("Belos Parameters");
  blocksmootherparams2.sublist("MueLu Parameters");

  arterydis_->compute_null_space_if_necessary(blocksmootherparams2);
}

/*--------------------------------------------------------------------------*
 | solve linear system of equations                        kremheller 04/18 |
 *--------------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::MeshtyingStrategyArtery::linear_solve(
    Teuchos::RCP<Core::LinAlg::Solver> solver, Teuchos::RCP<Core::LinAlg::SparseOperator> sysmat,
    Teuchos::RCP<Epetra_Vector> increment, Teuchos::RCP<Epetra_Vector> residual,
    Core::LinAlg::SolverParams& solver_params)
{
  comb_systemmatrix_->complete();

  comb_increment_->PutScalar(0.0);

  // standard solver call
  // system is ready to solve since Dirichlet Boundary conditions have been applied in
  // setup_system_matrix or Evaluate
  solver_params.refactor = true;
  solver->solve(comb_systemmatrix_->epetra_operator(), comb_increment_, rhs_, solver_params);

  return;
}

/*----------------------------------------------------------------------*
 | Calculate problem specific norm                     kremheller 03/18 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::MeshtyingStrategyArtery::calculate_norms(std::vector<double>& preresnorm,
    std::vector<double>& incprenorm, std::vector<double>& prenorm,
    const Teuchos::RCP<const Epetra_Vector> increment)
{
  preresnorm.resize(2);
  incprenorm.resize(2);
  prenorm.resize(2);

  prenorm[0] = UTILS::calculate_vector_norm(vectornorminc_, porofluidmultitimint_->phinp());
  prenorm[1] = UTILS::calculate_vector_norm(vectornorminc_, artnettimint_->pressurenp());

  Teuchos::RCP<const Epetra_Vector> arterypressinc;
  Teuchos::RCP<const Epetra_Vector> porofluidinc;

  arttoporofluidcoupling_->extract_single_field_vectors(
      comb_increment_, porofluidinc, arterypressinc);

  incprenorm[0] = UTILS::calculate_vector_norm(vectornorminc_, porofluidinc);
  incprenorm[1] = UTILS::calculate_vector_norm(vectornorminc_, arterypressinc);

  Teuchos::RCP<const Epetra_Vector> arterypressrhs;
  Teuchos::RCP<const Epetra_Vector> porofluidrhs;

  arttoporofluidcoupling_->extract_single_field_vectors(rhs_, porofluidrhs, arterypressrhs);

  preresnorm[0] = UTILS::calculate_vector_norm(vectornormfres_, porofluidrhs);
  preresnorm[1] = UTILS::calculate_vector_norm(vectornormfres_, arterypressrhs);

  return;
}

/*----------------------------------------------------------------------*
 | create result test for this field                   kremheller 04/18 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::MeshtyingStrategyArtery::create_field_test()
{
  Teuchos::RCP<Core::UTILS::ResultTest> arteryresulttest = artnettimint_->create_field_test();
  Global::Problem::instance()->add_field_test(arteryresulttest);
  return;
}

/*----------------------------------------------------------------------*
 |  read restart data                                  kremheller 04/18 |
 -----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::MeshtyingStrategyArtery::read_restart(const int step)
{
  artnettimint_->read_restart(step);

  return;
}

/*----------------------------------------------------------------------*
 | output of solution vector to BINIO                  kremheller 04/18 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::MeshtyingStrategyArtery::output()
{
  if (porofluidmultitimint_->step() != 0) artnettimint_->output(false, Teuchos::null);

  return;
}

/*----------------------------------------------------------------------*
 | evaluate matrix and rhs                             kremheller 04/18 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::MeshtyingStrategyArtery::evaluate()
{
  arttoporofluidcoupling_->set_solution_vectors(
      porofluidmultitimint_->phinp(), porofluidmultitimint_->phin(), artnettimint_->pressurenp());

  // evaluate the coupling
  arttoporofluidcoupling_->evaluate(comb_systemmatrix_, rhs_);

  // evaluate artery
  artnettimint_->assemble_mat_and_rhs();
  // apply DBC
  artnettimint_->prepare_linear_solve();

  // SetupCoupledArteryPoroFluidSystem();
  arttoporofluidcoupling_->setup_system(comb_systemmatrix_, rhs_,
      porofluidmultitimint_->system_matrix(), artnettimint_->system_matrix(),
      porofluidmultitimint_->rhs(), artnettimint_->rhs(),
      porofluidmultitimint_->get_dbc_map_extractor(), artnettimint_->get_dbc_map_extractor());

  return;
}

/*----------------------------------------------------------------------*
 | extract and update                                  kremheller 04/18 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector>
POROFLUIDMULTIPHASE::MeshtyingStrategyArtery::extract_and_update_iter(
    const Teuchos::RCP<const Epetra_Vector> inc)
{
  Teuchos::RCP<const Epetra_Vector> arterypressinc;
  Teuchos::RCP<const Epetra_Vector> porofluidinc;

  arttoporofluidcoupling_->extract_single_field_vectors(inc, porofluidinc, arterypressinc);

  artnettimint_->update_iter(arterypressinc);

  return porofluidinc;
}

/*----------------------------------------------------------------------*
 | artery dof row map                                  kremheller 04/18 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> POROFLUIDMULTIPHASE::MeshtyingStrategyArtery::artery_dof_row_map()
    const
{
  return arttoporofluidcoupling_->artery_dof_row_map();
}

/*-----------------------------------------------------------------------*
 | access to block system matrix of artery poro problem kremheller 04/18 |
 *-----------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase>
POROFLUIDMULTIPHASE::MeshtyingStrategyArtery::artery_porofluid_sysmat() const
{
  return comb_systemmatrix_;
}

/*----------------------------------------------------------------------*
 | return coupled residual                             kremheller 05/18 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector>
POROFLUIDMULTIPHASE::MeshtyingStrategyArtery::artery_porofluid_rhs() const
{
  return rhs_;
}

/*----------------------------------------------------------------------*
 | extract and update                                  kremheller 04/18 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> POROFLUIDMULTIPHASE::MeshtyingStrategyArtery::combined_increment(
    const Teuchos::RCP<const Epetra_Vector> inc) const
{
  return comb_increment_;
}

/*----------------------------------------------------------------------*
 | check initial fields                                kremheller 06/18 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::MeshtyingStrategyArtery::check_initial_fields(
    Teuchos::RCP<const Epetra_Vector> vec_cont) const
{
  arttoporofluidcoupling_->check_initial_fields(vec_cont, artnettimint_->pressurenp());
  return;
}

/*-------------------------------------------------------------------------*
 | set element pairs that are close                       kremheller 03/19 |
 *------------------------------------------------------------------------ */
void POROFLUIDMULTIPHASE::MeshtyingStrategyArtery::set_nearby_ele_pairs(
    const std::map<int, std::set<int>>* nearbyelepairs)
{
  arttoporofluidcoupling_->set_nearby_ele_pairs(nearbyelepairs);
  return;
}

/*-------------------------------------------------------------------------*
 | setup the strategy                                     kremheller 03/19 |
 *------------------------------------------------------------------------ */
void POROFLUIDMULTIPHASE::MeshtyingStrategyArtery::setup()
{
  arttoporofluidcoupling_->setup();
  return;
}

/*----------------------------------------------------------------------*
 | apply mesh movement                                 kremheller 06/18 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::MeshtyingStrategyArtery::apply_mesh_movement() const
{
  arttoporofluidcoupling_->apply_mesh_movement();
  return;
}

/*----------------------------------------------------------------------*
 | access to blood vessel volume fraction              kremheller 10/19 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector>
POROFLUIDMULTIPHASE::MeshtyingStrategyArtery::blood_vessel_volume_fraction()
{
  return arttoporofluidcoupling_->blood_vessel_volume_fraction();
}

FOUR_C_NAMESPACE_CLOSE
