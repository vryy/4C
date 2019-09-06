/*----------------------------------------------------------------------*/
/*! \file

\brief base level-set algorithm

    detailed description in header file levelset_algorithm.H

\level 2

\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236

 *------------------------------------------------------------------------------------------------*/


#include "levelset_algorithm.H"

#include "levelset_intersection_utils.H"

#include "../drt_scatra/scatra_resulttest.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_io/io.H"

#include <Teuchos_TimeMonitor.hpp>


/*----------------------------------------------------------------------*
 | constructor                                          rasthofer 09/13 |
 *----------------------------------------------------------------------*/
SCATRA::LevelSetAlgorithm::LevelSetAlgorithm(Teuchos::RCP<DRT::Discretization> dis,
    Teuchos::RCP<LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
    Teuchos::RCP<Teuchos::ParameterList> extraparams, Teuchos::RCP<IO::DiscretizationWriter> output)
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
 | deconstructor                                        rasthofer 09/13 |
 *----------------------------------------------------------------------*/
SCATRA::LevelSetAlgorithm::~LevelSetAlgorithm() {}


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
  GetInitialVolumeOfMinusDomain(phinp_, discret_, initvolminus_);

  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices: local <-> global dof numbering
  // -------------------------------------------------------------------
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // -------------------------------------------------------------------
  //         initialize reinitialization
  // -------------------------------------------------------------------
  // get reinitialization strategy
  reinitaction_ = DRT::INPUT::IntegralValue<INPAR::SCATRA::ReInitialAction>(
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
      reinitband_ = DRT::INPUT::IntegralValue<int>(
          levelsetparams_->sublist("REINITIALIZATION"), "REINITBAND");
    }

    // set parameters for reinitialization equation
    if (reinitaction_ == INPAR::SCATRA::reinitaction_sussman)
    {
      // vector for initial phi (solution of level-set equation) of reinitialization process
      initialphireinit_ = LINALG::CreateVector(*dofrowmap, true);

      // get pseudo-time step size
      dtau_ = levelsetparams_->sublist("REINITIALIZATION").get<double>("TIMESTEPREINIT");

      // number of pseudo-time steps to capture gradient around zero level-set
      pseudostepmax_ = levelsetparams_->sublist("REINITIALIZATION").get<int>("NUMSTEPSREINIT");

      // theta of ost time integration
      thetareinit_ = levelsetparams_->sublist("REINITIALIZATION").get<double>("THETAREINIT");

      // tolerance for convergence of reinitialization equation
      reinit_tol_ = levelsetparams_->sublist("REINITIALIZATION").get<double>("CONVTOL_REINIT");

      // flag to activate corrector step
      reinitcorrector_ = DRT::INPUT::IntegralValue<int>(
          levelsetparams_->sublist("REINITIALIZATION"), "CORRECTOR_STEP");

      // flag to activate calculation of node-based velocity
      useprojectedreinitvel_ = DRT::INPUT::IntegralValue<INPAR::SCATRA::VelReinit>(
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
      lsdim_ = DRT::INPUT::IntegralValue<INPAR::SCATRA::LSDim>(
          levelsetparams_->sublist("REINITIALIZATION"), "DIMENSION");
    }

    if (reinitaction_ == INPAR::SCATRA::reinitaction_ellipticeq)
    {
      // number of iterations steps to solve nonlinear equation
      pseudostepmax_ = levelsetparams_->sublist("REINITIALIZATION").get<int>("NUMSTEPSREINIT");

      // tolerance for convergence of reinitialization equation
      reinit_tol_ = levelsetparams_->sublist("REINITIALIZATION").get<double>("CONVTOL_REINIT");

      // get dimension
      lsdim_ = DRT::INPUT::IntegralValue<INPAR::SCATRA::LSDim>(
          levelsetparams_->sublist("REINITIALIZATION"), "DIMENSION");

      // use L2-projection of grad phi and related quantities
      projection_ = DRT::INPUT::IntegralValue<int>(
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
    reinitvolcorrection_ = DRT::INPUT::IntegralValue<int>(
        levelsetparams_->sublist("REINITIALIZATION"), "REINITVOLCORRECTION");

    // initialize level-set to signed distance function if required
    if (DRT::INPUT::IntegralValue<int>(
            levelsetparams_->sublist("REINITIALIZATION"), "REINIT_INITIAL"))
    {
      Reinitialization();

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
      DRT::INPUT::IntegralValue<int>(*levelsetparams_, "EXTRACT_INTERFACE_VEL");
  if (extract_interface_vel_)
  {
    // set number of element layers around interface where velocity field form Navier-Stokes is kept
    convel_layers_ = levelsetparams_->get<int>("NUM_CONVEL_LAYERS");
    if (convel_layers_ < 1)
      dserror(
          "Set number of element layers around interface where velocity "
          "field form Navier-Stokes should be kept");
  }

  // set flag for modification of convective velocity at contact points
  // check whether there are level-set contact point conditions
  std::vector<DRT::Condition*> lscontactpoint;
  discret_->GetCondition("LsContact", lscontactpoint);

  if (not lscontactpoint.empty()) cpbc_ = true;

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SCATRA::LevelSetAlgorithm::GetInitialVolumeOfMinusDomain(
    const Teuchos::RCP<const Epetra_Vector>& phinp,
    const Teuchos::RCP<const DRT::Discretization>& scatradis, double& volumedomainminus) const
{
  double volplus = 0.0;
  double surf = 0.0;
  std::map<int, GEO::BoundaryIntCells> interface;
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
  CheckIsInit();
  CheckIsSetup();

  // provide information about initial field (do not do for restarts!)
  if (Step() == 0)
  {
    // write out initial state
    Output();

    // compute error for problems with analytical solution (initial field!)
    EvaluateErrorComparedToAnalyticalSol();
  }

  while ((step_ < stepmax_) and ((time_ + EPS12) < maxtime_))
  {
    // -------------------------------------------------------------------
    // prepare time step
    // -------------------------------------------------------------------
    PrepareTimeStep();

    // -------------------------------------------------------------------
    //                  solve level-set equation
    // -------------------------------------------------------------------
    Solve();

    // -----------------------------------------------------------------
    //                     reinitialize level-set
    // -----------------------------------------------------------------
    // will be done only if required
    Reinitialization();

    // -------------------------------------------------------------------
    //                         update solution
    //        current solution becomes old solution of next time step
    // -------------------------------------------------------------------
    UpdateState();

    // -------------------------------------------------------------------
    //       evaluate error for problems with analytical solution
    // -------------------------------------------------------------------
    EvaluateErrorComparedToAnalyticalSol();

    // -------------------------------------------------------------------
    //                         output of solution
    // -------------------------------------------------------------------
    Output();

  }  // while

  return;
}


/*----------------------------------------------------------------------*
 | setup the variables to do a new time step            rasthofer 09/13 |
 *----------------------------------------------------------------------*/
void SCATRA::LevelSetAlgorithm::PrepareTimeStep()
{
  // prepare basic scalar transport solver
  ScaTraTimIntImpl::PrepareTimeStep();

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
    NonlinearSolve();
  else
    LinearSolve();

  return;
}


/*----------------------------------------------------------------------*
 | reinitialize level-set                               rasthofer 09/13 |
 *----------------------------------------------------------------------*/
void SCATRA::LevelSetAlgorithm::Reinitialization()
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
      std::map<int, GEO::BoundaryIntCells> zerolevelset;
      zerolevelset.clear();
      CaptureInterface(zerolevelset);

      // -----------------------------------------------------------------
      //                    reinitialize level-set
      // -----------------------------------------------------------------
      // select reinitialization method
      switch (reinitaction_)
      {
        case INPAR::SCATRA::reinitaction_signeddistancefunction:
        {
          // reinitialization via signed distance to interface
          ReinitGeo(zerolevelset);
          break;
        }
        case INPAR::SCATRA::reinitaction_sussman:
        {
          // reinitialization via solving equation to steady state
          ReinitEq();
          break;
        }
        case INPAR::SCATRA::reinitaction_ellipticeq:
        {
          // reinitialization via elliptic equation
          ReinitElliptic(zerolevelset);
          break;
        }
        default:
        {
          dserror("Unknown reinitialization method!");
          break;
        }
      }

      // -----------------------------------------------------------------
      //                       correct volume
      // -----------------------------------------------------------------
      // since the reinitialization may shift the interface,
      // this method allows for correcting the volume
      if (reinitvolcorrection_) CorrectVolume();
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | output of solution                                   rasthofer 09/13 |
 *----------------------------------------------------------------------*/
void SCATRA::LevelSetAlgorithm::Output(const int num)
{
  // time measurement: output of solution
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:    + output of solution");

  // -----------------------------------------------------------------
  //      standard paraview and gmsh output for scalar field
  // -----------------------------------------------------------------

  // solution output and potentially restart data and/or flux data
  if (DoOutput())
  {
    // step number and time (only after that data output is possible)
    output_->NewStep(step_, time_);

    // write domain decomposition for visualization (only once at step "upres"!)
    if (step_ == upres_) output_->WriteElementData(true);

    // write state vectors
    OutputState();

    // write output to Gmsh postprocessing files
    if (outputgmsh_) OutputToGmsh(step_, time_);

    // add restart data
    if (step_ % uprestart_ == 0) OutputRestart();
  }

  // -----------------------------------------------------------------
  //             further level-set specific values
  // -----------------------------------------------------------------
  OutputOfLevelSetSpecificValues();

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void SCATRA::LevelSetAlgorithm::OutputOfLevelSetSpecificValues()
{
  // capture interface, evalute mass conservation, write to file
  std::map<int, GEO::BoundaryIntCells> zerolevelset;
  zerolevelset.clear();
  CaptureInterface(zerolevelset, true);
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
