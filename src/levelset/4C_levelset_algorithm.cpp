/*----------------------------------------------------------------------*/
/*! \file

\brief base level-set algorithm

    detailed description in header file levelset_algorithm.H

\level 2


 *------------------------------------------------------------------------------------------------*/


#include "4C_levelset_algorithm.hpp"

#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_levelset_intersection_utils.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_scatra_resulttest.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 | constructor                                          rasthofer 09/13 |
 *----------------------------------------------------------------------*/
ScaTra::LevelSetAlgorithm::LevelSetAlgorithm(Teuchos::RCP<Core::FE::Discretization> dis,
    Teuchos::RCP<Core::LinAlg::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
    Teuchos::RCP<Teuchos::ParameterList> extraparams,
    Teuchos::RCP<Core::IO::DiscretizationWriter> output)
    : ScaTraTimIntImpl(dis, solver, sctratimintparams, extraparams, output),
      levelsetparams_(params),
      reinitaction_(Inpar::ScaTra::reinitaction_none),
      switchreinit_(false),
      pseudostepmax_(0),
      pseudostep_(0),
      dtau_(0.0),
      thetareinit_(0.0),
      initvolminus_(0.0),
      initialphireinit_(Teuchos::null),
      nb_grad_val_(Teuchos::null),
      reinitinterval_(-1),
      reinitband_(false),
      reinitbandwidth_(-1.0),
      reinitcorrector_(true),
      useprojectedreinitvel_(Inpar::ScaTra::vel_reinit_integration_point_based),
      lsdim_(Inpar::ScaTra::ls_3D),
      projection_(true),
      reinit_tol_(-1.0),
      reinitvolcorrection_(false),
      interface_eleq_(Teuchos::null),
      extract_interface_vel_(false),
      convel_layers_(-1),
      cpbc_(false)
{
  // DO NOT DEFINE ANY STATE VECTORS HERE (i.e., vectors based on row or column maps)
  // this is important since we have problems which require an extended ghosting
  // this has to be done before all state vectors are initialized
  return;
}



/*----------------------------------------------------------------------*
 | initialize algorithm                                     rauch 09/16 |
 *----------------------------------------------------------------------*/
void ScaTra::LevelSetAlgorithm::init()
{
  // todo #initsetupissue
  // DO NOT CALL init() IN ScaTraTimIntImpl
  // issue with writeflux and probably scalarhandler_
  // this should not be

  return;
}


/*----------------------------------------------------------------------*
 | setup algorithm                                          rauch 09/16 |
 *----------------------------------------------------------------------*/
void ScaTra::LevelSetAlgorithm::setup()
{
  // todo #initsetupissue
  // DO NOT CALL setup() IN ScaTraTimIntImpl
  // issue with writeflux and probably scalarhandler_
  // this should not be

  // -------------------------------------------------------------------
  //                         setup domains
  // -------------------------------------------------------------------
  get_initial_volume_of_minus_domain(phinp_, discret_, initvolminus_);

  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices: local <-> global dof numbering
  // -------------------------------------------------------------------
  const Epetra_Map* dofrowmap = discret_->dof_row_map();

  // -------------------------------------------------------------------
  //         initialize reinitialization
  // -------------------------------------------------------------------
  // get reinitialization strategy
  reinitaction_ = Core::UTILS::IntegralValue<Inpar::ScaTra::ReInitialAction>(
      levelsetparams_->sublist("REINITIALIZATION"), "REINITIALIZATION");

  if (reinitaction_ != Inpar::ScaTra::reinitaction_none)
  {
    // how often to perform reinitialization
    reinitinterval_ = levelsetparams_->sublist("REINITIALIZATION").get<int>("REINITINTERVAL");

    // define band for reinitialization (may either be used for geometric reinitialization or
    // reinitialization equation
    reinitbandwidth_ = levelsetparams_->sublist("REINITIALIZATION").get<double>("REINITBANDWIDTH");

    // set parameters for geometric reinitialization
    if (reinitaction_ == Inpar::ScaTra::reinitaction_signeddistancefunction)
    {
      // reinitialization within band around interface only
      reinitband_ = Core::UTILS::IntegralValue<int>(
          levelsetparams_->sublist("REINITIALIZATION"), "REINITBAND");
    }

    // set parameters for reinitialization equation
    if (reinitaction_ == Inpar::ScaTra::reinitaction_sussman)
    {
      // vector for initial phi (solution of level-set equation) of reinitialization process
      initialphireinit_ = Core::LinAlg::CreateVector(*dofrowmap, true);

      // get pseudo-time step size
      dtau_ = levelsetparams_->sublist("REINITIALIZATION").get<double>("TIMESTEPREINIT");

      // number of pseudo-time steps to capture gradient around zero level-set
      pseudostepmax_ = levelsetparams_->sublist("REINITIALIZATION").get<int>("NUMSTEPSREINIT");

      // theta of ost time integration
      thetareinit_ = levelsetparams_->sublist("REINITIALIZATION").get<double>("THETAREINIT");

      // tolerance for convergence of reinitialization equation
      reinit_tol_ = levelsetparams_->sublist("REINITIALIZATION").get<double>("CONVTOL_REINIT");

      // flag to activate corrector step
      reinitcorrector_ = Core::UTILS::IntegralValue<int>(
          levelsetparams_->sublist("REINITIALIZATION"), "CORRECTOR_STEP");

      // flag to activate calculation of node-based velocity
      useprojectedreinitvel_ = Core::UTILS::IntegralValue<Inpar::ScaTra::VelReinit>(
          levelsetparams_->sublist("REINITIALIZATION"), "VELREINIT");

      if (useprojectedreinitvel_ == Inpar::ScaTra::vel_reinit_node_based)
      {
        // vector for nodal velocity for reinitialization
        // velocities (always three velocity components per node)
        // (get noderowmap of discretization for creating this multivector)
        const Epetra_Map* noderowmap = discret_->node_row_map();
        nb_grad_val_ = Teuchos::rcp(new Epetra_MultiVector(*noderowmap, 3, true));
      }

      // get dimension
      lsdim_ = Core::UTILS::IntegralValue<Inpar::ScaTra::LSDim>(
          levelsetparams_->sublist("REINITIALIZATION"), "DIMENSION");
    }

    if (reinitaction_ == Inpar::ScaTra::reinitaction_ellipticeq)
    {
      // number of iterations steps to solve nonlinear equation
      pseudostepmax_ = levelsetparams_->sublist("REINITIALIZATION").get<int>("NUMSTEPSREINIT");

      // tolerance for convergence of reinitialization equation
      reinit_tol_ = levelsetparams_->sublist("REINITIALIZATION").get<double>("CONVTOL_REINIT");

      // get dimension
      lsdim_ = Core::UTILS::IntegralValue<Inpar::ScaTra::LSDim>(
          levelsetparams_->sublist("REINITIALIZATION"), "DIMENSION");

      // use L2-projection of grad phi and related quantities
      projection_ = Core::UTILS::IntegralValue<int>(
          levelsetparams_->sublist("REINITIALIZATION"), "PROJECTION");
      if (projection_ == true)
      {
        // vector for nodal level-set gradient for reinitialization
        // gradients (always three gradient components per node)
        // (get noderowmap of discretization for creating this multivector)
        const Epetra_Map* noderowmap = discret_->node_row_map();
        nb_grad_val_ = Teuchos::rcp(new Epetra_MultiVector(*noderowmap, 3, true));
      }
    }

    // flag to correct volume after reinitialization
    reinitvolcorrection_ = Core::UTILS::IntegralValue<int>(
        levelsetparams_->sublist("REINITIALIZATION"), "REINITVOLCORRECTION");

    // initialize level-set to signed distance function if required
    if (Core::UTILS::IntegralValue<int>(
            levelsetparams_->sublist("REINITIALIZATION"), "REINIT_INITIAL"))
    {
      reinitialization();

      // reset phin vector
      // remark: for BDF2, we do not have to set phinm,
      //        since it is not used for the first time step as a backward Euler step is performed
      phin_->Update(1.0, *phinp_, 0.0);

      // reinitialization is done, reset flag
      switchreinit_ = false;
    }
  }

  // -------------------------------------------------------------------
  //       initialize treatment of velocity from Navier-Stokes
  // -------------------------------------------------------------------
  // set potential extraction of interface velocity
  extract_interface_vel_ =
      Core::UTILS::IntegralValue<int>(*levelsetparams_, "EXTRACT_INTERFACE_VEL");
  if (extract_interface_vel_)
  {
    // set number of element layers around interface where velocity field form Navier-Stokes is kept
    convel_layers_ = levelsetparams_->get<int>("NUM_CONVEL_LAYERS");
    if (convel_layers_ < 1)
      FOUR_C_THROW(
          "Set number of element layers around interface where velocity "
          "field form Navier-Stokes should be kept");
  }

  // set flag for modification of convective velocity at contact points
  // check whether there are level-set contact point conditions
  std::vector<Core::Conditions::Condition*> lscontactpoint;
  discret_->get_condition("LsContact", lscontactpoint);

  if (not lscontactpoint.empty()) cpbc_ = true;

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::LevelSetAlgorithm::get_initial_volume_of_minus_domain(
    const Teuchos::RCP<const Epetra_Vector>& phinp,
    const Teuchos::RCP<const Core::FE::Discretization>& scatradis, double& volumedomainminus) const
{
  double volplus = 0.0;
  double surf = 0.0;
  std::map<int, Core::Geo::BoundaryIntCells> interface;
  interface.clear();
  // reconstruct interface and calculate volumes, etc ...
  ScaTra::LevelSet::Intersection intersect;
  intersect.capture_zero_level_set(phinp, scatradis, volumedomainminus, volplus, surf, interface);
}

/*----------------------------------------------------------------------*
 | time loop                                            rasthofer 09/13 |
 *----------------------------------------------------------------------*/
void ScaTra::LevelSetAlgorithm::time_loop()
{
  // safety check
  check_is_init();
  check_is_setup();

  // provide information about initial field (do not do for restarts!)
  if (step() == 0)
  {
    // write out initial state
    check_and_write_output_and_restart();

    // compute error for problems with analytical solution (initial field!)
    evaluate_error_compared_to_analytical_sol();
  }

  while ((step_ < stepmax_) and ((time_ + 1e-12) < maxtime_))
  {
    // -------------------------------------------------------------------
    // prepare time step
    // -------------------------------------------------------------------
    prepare_time_step();

    // -------------------------------------------------------------------
    //                  solve level-set equation
    // -------------------------------------------------------------------
    solve();

    // -----------------------------------------------------------------
    //                     reinitialize level-set
    // -----------------------------------------------------------------
    // will be done only if required
    reinitialization();

    // -------------------------------------------------------------------
    //                         update solution
    //        current solution becomes old solution of next time step
    // -------------------------------------------------------------------
    update_state();

    // -------------------------------------------------------------------
    //       evaluate error for problems with analytical solution
    // -------------------------------------------------------------------
    evaluate_error_compared_to_analytical_sol();

    // -------------------------------------------------------------------
    //                         output of solution
    // -------------------------------------------------------------------
    check_and_write_output_and_restart();

  }  // while

  return;
}


/*----------------------------------------------------------------------*
 | setup the variables to do a new time step            rasthofer 09/13 |
 *----------------------------------------------------------------------*/
void ScaTra::LevelSetAlgorithm::prepare_time_step()
{
  // prepare basic scalar transport solver
  ScaTraTimIntImpl::prepare_time_step();

  return;
}


/*----------------------------------------------------------------------*
 | solve level-set equation and perform correction      rasthofer 09/13 |
 *----------------------------------------------------------------------*/
void ScaTra::LevelSetAlgorithm::solve()
{
  // -----------------------------------------------------------------
  //                    solve level-set equation
  // -----------------------------------------------------------------
  if (solvtype_ == Inpar::ScaTra::solvertype_nonlinear)
    nonlinear_solve();
  else
    linear_solve();

  return;
}


/*----------------------------------------------------------------------*
 | reinitialize level-set                               rasthofer 09/13 |
 *----------------------------------------------------------------------*/
void ScaTra::LevelSetAlgorithm::reinitialization()
{
  // check for reinitialization action first

  if (reinitaction_ != Inpar::ScaTra::reinitaction_none)
  {
    // check if level-set field should be reinitialized in this time step
    if (step_ % reinitinterval_ == 0)
    {
      // -----------------------------------------------------------------
      //             capture interface before reinitialization
      // -----------------------------------------------------------------
      // get volume distribution before reinitialization
      // initalize structure holding the interface
      // potentially required for reinitialization via signed distance to interface
      std::map<int, Core::Geo::BoundaryIntCells> zerolevelset;
      zerolevelset.clear();
      capture_interface(zerolevelset);

      // -----------------------------------------------------------------
      //                    reinitialize level-set
      // -----------------------------------------------------------------
      // select reinitialization method
      switch (reinitaction_)
      {
        case Inpar::ScaTra::reinitaction_signeddistancefunction:
        {
          // reinitialization via signed distance to interface
          reinit_geo(zerolevelset);
          break;
        }
        case Inpar::ScaTra::reinitaction_sussman:
        {
          // reinitialization via solving equation to steady state
          reinit_eq();
          break;
        }
        case Inpar::ScaTra::reinitaction_ellipticeq:
        {
          // reinitialization via elliptic equation
          reinit_elliptic(zerolevelset);
          break;
        }
        default:
        {
          FOUR_C_THROW("Unknown reinitialization method!");
          break;
        }
      }

      // -----------------------------------------------------------------
      //                       correct volume
      // -----------------------------------------------------------------
      // since the reinitialization may shift the interface,
      // this method allows for correcting the volume
      if (reinitvolcorrection_) correct_volume();
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::LevelSetAlgorithm::check_and_write_output_and_restart()
{
  // time measurement: output of solution
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:    + output of solution");

  // -----------------------------------------------------------------
  //      standard paraview and gmsh output for scalar field
  // -----------------------------------------------------------------

  // solution output and potentially restart data and/or flux data
  if (is_result_step())
  {
    // step number and time (only after that data output is possible)
    output_->new_step(step_, time_);

    // write domain decomposition for visualization (only once at step "upres"!)
    if (step_ == upres_) output_->write_element_data(true);

    // write state vectors
    output_state();

    // write output to Gmsh postprocessing files
    if (outputgmsh_) output_to_gmsh(step_, time_);
  }

  // add restart data
  if (is_restart_step()) write_restart();

  // -----------------------------------------------------------------
  //             further level-set specific values
  // -----------------------------------------------------------------
  output_of_level_set_specific_values();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ScaTra::LevelSetAlgorithm::output_of_level_set_specific_values()
{
  // capture interface, evalute mass conservation, write to file
  std::map<int, Core::Geo::BoundaryIntCells> zerolevelset;
  zerolevelset.clear();
  capture_interface(zerolevelset, true);
}

/*----------------------------------------------------------------------*
 | perform result test                                  rasthofer 01/14 |
 *----------------------------------------------------------------------*/
void ScaTra::LevelSetAlgorithm::test_results()
{
  problem_->add_field_test(Teuchos::rcp(new ScaTra::ScaTraResultTest(Teuchos::rcp(this, false))));
  problem_->test_all(discret_->get_comm());

  return;
}


/*----------------------------------------------------------------------*
 | set time and step value                              rasthofer 04/14 |
 *----------------------------------------------------------------------*/
void ScaTra::LevelSetAlgorithm::set_time_step(const double time, const int step)
{
  // call base class function
  ScaTra::ScaTraTimIntImpl::set_time_step(time, step);

  return;
}

FOUR_C_NAMESPACE_CLOSE
