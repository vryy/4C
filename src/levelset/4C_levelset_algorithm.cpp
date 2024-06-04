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
SCATRA::LevelSetAlgorithm::LevelSetAlgorithm(Teuchos::RCP<DRT::Discretization> dis,
    Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
    Teuchos::RCP<Teuchos::ParameterList> extraparams,
    Teuchos::RCP<CORE::IO::DiscretizationWriter> output)
    : ScaTraTimIntImpl(dis, solver, sctratimintparams, extraparams, output),
      levelsetparams_(params),
      reinitaction_(INPAR::SCATRA::reinitaction_none),
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
      useprojectedreinitvel_(INPAR::SCATRA::vel_reinit_integration_point_based),
      lsdim_(INPAR::SCATRA::ls_3D),
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
void SCATRA::LevelSetAlgorithm::Init()
{
  // todo #initsetupissue
  // DO NOT CALL Init() IN ScaTraTimIntImpl
  // issue with writeflux and probably scalarhandler_
  // this should not be

  return;
}


/*----------------------------------------------------------------------*
 | setup algorithm                                          rauch 09/16 |
 *----------------------------------------------------------------------*/
void SCATRA::LevelSetAlgorithm::Setup()
{
  // todo #initsetupissue
  // DO NOT CALL Setup() IN ScaTraTimIntImpl
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
  reinitaction_ = CORE::UTILS::IntegralValue<INPAR::SCATRA::ReInitialAction>(
      levelsetparams_->sublist("REINITIALIZATION"), "REINITIALIZATION");

  if (reinitaction_ != INPAR::SCATRA::reinitaction_none)
  {
    // how often to perform reinitialization
    reinitinterval_ = levelsetparams_->sublist("REINITIALIZATION").get<int>("REINITINTERVAL");

    // define band for reinitialization (may either be used for geometric reinitialization or
    // reinitialization equation
    reinitbandwidth_ = levelsetparams_->sublist("REINITIALIZATION").get<double>("REINITBANDWIDTH");

    // set parameters for geometric reinitialization
    if (reinitaction_ == INPAR::SCATRA::reinitaction_signeddistancefunction)
    {
      // reinitialization within band around interface only
      reinitband_ = CORE::UTILS::IntegralValue<int>(
          levelsetparams_->sublist("REINITIALIZATION"), "REINITBAND");
    }

    // set parameters for reinitialization equation
    if (reinitaction_ == INPAR::SCATRA::reinitaction_sussman)
    {
      // vector for initial phi (solution of level-set equation) of reinitialization process
      initialphireinit_ = CORE::LINALG::CreateVector(*dofrowmap, true);

      // get pseudo-time step size
      dtau_ = levelsetparams_->sublist("REINITIALIZATION").get<double>("TIMESTEPREINIT");

      // number of pseudo-time steps to capture gradient around zero level-set
      pseudostepmax_ = levelsetparams_->sublist("REINITIALIZATION").get<int>("NUMSTEPSREINIT");

      // theta of ost time integration
      thetareinit_ = levelsetparams_->sublist("REINITIALIZATION").get<double>("THETAREINIT");

      // tolerance for convergence of reinitialization equation
      reinit_tol_ = levelsetparams_->sublist("REINITIALIZATION").get<double>("CONVTOL_REINIT");

      // flag to activate corrector step
      reinitcorrector_ = CORE::UTILS::IntegralValue<int>(
          levelsetparams_->sublist("REINITIALIZATION"), "CORRECTOR_STEP");

      // flag to activate calculation of node-based velocity
      useprojectedreinitvel_ = CORE::UTILS::IntegralValue<INPAR::SCATRA::VelReinit>(
          levelsetparams_->sublist("REINITIALIZATION"), "VELREINIT");

      if (useprojectedreinitvel_ == INPAR::SCATRA::vel_reinit_node_based)
      {
        // vector for nodal velocity for reinitialization
        // velocities (always three velocity components per node)
        // (get noderowmap of discretization for creating this multivector)
        const Epetra_Map* noderowmap = discret_->NodeRowMap();
        nb_grad_val_ = Teuchos::rcp(new Epetra_MultiVector(*noderowmap, 3, true));
      }

      // get dimension
      lsdim_ = CORE::UTILS::IntegralValue<INPAR::SCATRA::LSDim>(
          levelsetparams_->sublist("REINITIALIZATION"), "DIMENSION");
    }

    if (reinitaction_ == INPAR::SCATRA::reinitaction_ellipticeq)
    {
      // number of iterations steps to solve nonlinear equation
      pseudostepmax_ = levelsetparams_->sublist("REINITIALIZATION").get<int>("NUMSTEPSREINIT");

      // tolerance for convergence of reinitialization equation
      reinit_tol_ = levelsetparams_->sublist("REINITIALIZATION").get<double>("CONVTOL_REINIT");

      // get dimension
      lsdim_ = CORE::UTILS::IntegralValue<INPAR::SCATRA::LSDim>(
          levelsetparams_->sublist("REINITIALIZATION"), "DIMENSION");

      // use L2-projection of grad phi and related quantities
      projection_ = CORE::UTILS::IntegralValue<int>(
          levelsetparams_->sublist("REINITIALIZATION"), "PROJECTION");
      if (projection_ == true)
      {
        // vector for nodal level-set gradient for reinitialization
        // gradients (always three gradient components per node)
        // (get noderowmap of discretization for creating this multivector)
        const Epetra_Map* noderowmap = discret_->NodeRowMap();
        nb_grad_val_ = Teuchos::rcp(new Epetra_MultiVector(*noderowmap, 3, true));
      }
    }

    // flag to correct volume after reinitialization
    reinitvolcorrection_ = CORE::UTILS::IntegralValue<int>(
        levelsetparams_->sublist("REINITIALIZATION"), "REINITVOLCORRECTION");

    // initialize level-set to signed distance function if required
    if (CORE::UTILS::IntegralValue<int>(
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
      CORE::UTILS::IntegralValue<int>(*levelsetparams_, "EXTRACT_INTERFACE_VEL");
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
  std::vector<CORE::Conditions::Condition*> lscontactpoint;
  discret_->GetCondition("LsContact", lscontactpoint);

  if (not lscontactpoint.empty()) cpbc_ = true;

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SCATRA::LevelSetAlgorithm::get_initial_volume_of_minus_domain(
    const Teuchos::RCP<const Epetra_Vector>& phinp,
    const Teuchos::RCP<const DRT::Discretization>& scatradis, double& volumedomainminus) const
{
  double volplus = 0.0;
  double surf = 0.0;
  std::map<int, CORE::GEO::BoundaryIntCells> interface;
  interface.clear();
  // reconstruct interface and calculate volumes, etc ...
  SCATRA::LEVELSET::Intersection intersect;
  intersect.CaptureZeroLevelSet(phinp, scatradis, volumedomainminus, volplus, surf, interface);
}

/*----------------------------------------------------------------------*
 | time loop                                            rasthofer 09/13 |
 *----------------------------------------------------------------------*/
void SCATRA::LevelSetAlgorithm::TimeLoop()
{
  // safety check
  check_is_init();
  check_is_setup();

  // provide information about initial field (do not do for restarts!)
  if (Step() == 0)
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
    Solve();

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
void SCATRA::LevelSetAlgorithm::prepare_time_step()
{
  // prepare basic scalar transport solver
  ScaTraTimIntImpl::prepare_time_step();

  return;
}


/*----------------------------------------------------------------------*
 | solve level-set equation and perform correction      rasthofer 09/13 |
 *----------------------------------------------------------------------*/
void SCATRA::LevelSetAlgorithm::Solve()
{
  // -----------------------------------------------------------------
  //                    solve level-set equation
  // -----------------------------------------------------------------
  if (solvtype_ == INPAR::SCATRA::solvertype_nonlinear)
    nonlinear_solve();
  else
    linear_solve();

  return;
}


/*----------------------------------------------------------------------*
 | reinitialize level-set                               rasthofer 09/13 |
 *----------------------------------------------------------------------*/
void SCATRA::LevelSetAlgorithm::reinitialization()
{
  // check for reinitialization action first

  if (reinitaction_ != INPAR::SCATRA::reinitaction_none)
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
      std::map<int, CORE::GEO::BoundaryIntCells> zerolevelset;
      zerolevelset.clear();
      capture_interface(zerolevelset);

      // -----------------------------------------------------------------
      //                    reinitialize level-set
      // -----------------------------------------------------------------
      // select reinitialization method
      switch (reinitaction_)
      {
        case INPAR::SCATRA::reinitaction_signeddistancefunction:
        {
          // reinitialization via signed distance to interface
          reinit_geo(zerolevelset);
          break;
        }
        case INPAR::SCATRA::reinitaction_sussman:
        {
          // reinitialization via solving equation to steady state
          reinit_eq();
          break;
        }
        case INPAR::SCATRA::reinitaction_ellipticeq:
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
void SCATRA::LevelSetAlgorithm::check_and_write_output_and_restart()
{
  // time measurement: output of solution
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:    + output of solution");

  // -----------------------------------------------------------------
  //      standard paraview and gmsh output for scalar field
  // -----------------------------------------------------------------

  // solution output and potentially restart data and/or flux data
  if (IsResultStep())
  {
    // step number and time (only after that data output is possible)
    output_->NewStep(step_, time_);

    // write domain decomposition for visualization (only once at step "upres"!)
    if (step_ == upres_) output_->WriteElementData(true);

    // write state vectors
    output_state();

    // write output to Gmsh postprocessing files
    if (outputgmsh_) output_to_gmsh(step_, time_);
  }

  // add restart data
  if (IsRestartStep()) write_restart();

  // -----------------------------------------------------------------
  //             further level-set specific values
  // -----------------------------------------------------------------
  output_of_level_set_specific_values();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void SCATRA::LevelSetAlgorithm::output_of_level_set_specific_values()
{
  // capture interface, evalute mass conservation, write to file
  std::map<int, CORE::GEO::BoundaryIntCells> zerolevelset;
  zerolevelset.clear();
  capture_interface(zerolevelset, true);
}

/*----------------------------------------------------------------------*
 | perform result test                                  rasthofer 01/14 |
 *----------------------------------------------------------------------*/
void SCATRA::LevelSetAlgorithm::TestResults()
{
  problem_->AddFieldTest(Teuchos::rcp(new SCATRA::ScaTraResultTest(Teuchos::rcp(this, false))));
  problem_->TestAll(discret_->Comm());

  return;
}


/*----------------------------------------------------------------------*
 | set time and step value                              rasthofer 04/14 |
 *----------------------------------------------------------------------*/
void SCATRA::LevelSetAlgorithm::SetTimeStep(const double time, const int step)
{
  // call base class function
  SCATRA::ScaTraTimIntImpl::SetTimeStep(time, step);

  return;
}

FOUR_C_NAMESPACE_CLOSE
