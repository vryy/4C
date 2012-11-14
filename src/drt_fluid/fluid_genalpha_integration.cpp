/*!----------------------------------------------------------------------
\file fluid_genalpha_integration.cpp

\brief Time integration according to dis. C. Whiting

Documentation see header.


<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#include "../drt_fluid/fluid_genalpha_integration.H"
#include "../drt_lib/drt_globalproblem.H"
#include "fluid_utils.H"
#include "../drt_io/io_control.H"
#include "../drt_nurbs_discret/drt_nurbs_discret.H"
#include "../drt_nurbs_discret/drt_apply_nurbs_initial_condition.H"
#include "../drt_fluid_ele/fluid_ele_action.H"

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Constructor (public)                                     gammi 06/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
FLD::FluidGenAlphaIntegration::FluidGenAlphaIntegration(
  const Teuchos::RCP<DRT::Discretization>&      actdis,
  const Teuchos::RCP<LINALG::Solver>&           solver,
  const Teuchos::RCP<ParameterList>&            params,
  const Teuchos::RCP<IO::DiscretizationWriter>& output,
  bool                                          alefluid,
  RefCountPtr<map<int,vector<int> > > pbcmapmastertoslave
):TimInt(actdis, solver, params, output),
  dis_   (discret_),
  alefluid_(alefluid),
  density_(1.0),
  writestresses_(params_->get("write stresses", 0)),
  pbcmapmastertoslave_(pbcmapmastertoslave),
  locsysman_(Teuchos::null)
{
  // -------------------------------------------------------------------
  // create timers and time monitor
  // -------------------------------------------------------------------
  timedyntot_           = TimeMonitor::getNewTimer("dynamic routine total"          );
  timedyninit_          = TimeMonitor::getNewTimer(" + initial phase"               );
  timedynloop_          = TimeMonitor::getNewTimer(" + time loop"                   );
  timenlnloop_          = TimeMonitor::getNewTimer("   + nonlinear iteration"       );
  timesparsitypattern_  = TimeMonitor::getNewTimer("      + set up and complete sparsity pattern");
  timeeleloop_          = TimeMonitor::getNewTimer("      + element calls"          );
  timenonlinup_         = TimeMonitor::getNewTimer("      + update and calc. of intermediate sols");
  timeapplydirich_      = TimeMonitor::getNewTimer("      + apply dirich cond."   );
  timeevaldirich_       = TimeMonitor::getNewTimer("      + evaluate dirich cond.");
  timesolver_           = TimeMonitor::getNewTimer("      + solver calls"         );
  timeout_              = TimeMonitor::getNewTimer("      + output and statistics");

  // time measurement --- start TimeMonitor tm0
  tm0_ref_        = Teuchos::rcp(new TimeMonitor(*timedyntot_ ));

  // time measurement --- start TimeMonitor tm7
  tm7_ref_        = Teuchos::rcp(new TimeMonitor(*timedyninit_ ));

  discret_->ComputeNullSpaceIfNecessary(solver_->Params(),true);

  // -------------------------------------------------------------------
  // check whether we have locsys BCs and create LocSysManager if so
  // -------------------------------------------------------------------
  {
    std::vector<DRT::Condition*> locsysconditions(0);
    discret_->GetCondition("Locsys", locsysconditions);
    if (locsysconditions.size())
    {
      locsysman_ = Teuchos::rcp(new DRT::UTILS::LocsysManager(*discret_,true));
    }
  }

  // ensure that degrees of freedom in the discretization have been set
  if (!discret_->Filled() || !actdis->HaveDofs()) discret_->FillComplete();

  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  // -------------------------------------------------------------------
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // -------------------------------------------------------------------
  // init some class variables (parallelism)

  // get the processor ID from the communicator
  myrank_  = discret_->Comm().MyPID();

  //--------------------------------------------------------------------
  // init some class variables (time integration)

  // time step size
  dt_     = params_->get<double>("time step size");
  // maximum number of timesteps
  endstep_= params_->get<int>   ("max number timesteps");
  // maximum simulation time
  endtime_= params_->get<double>("total time");

  // generalized alpha parameters
  // (choice of third parameter necessary but not sufficiant for second
  // order accuracy)
  //           gamma_  = 0.5 + alphaM_ - alphaF_
  alphaM_ = params_->get<double>("alpha_M");
  alphaF_ = params_->get<double>("alpha_F");
  gamma_  = params_->get<double>("gamma");

  // use of predictor
  predictor_ = params_->get<string>("predictor","steady_state_predictor");

  // parameter for linearisation scheme (fixed point like or newton like)
  newton_    = DRT::INPUT::get<INPAR::FLUID::LinearisationAction>(*params_, "Linearisation");

  itenum_    = 0;
  itemax_    = params_->get<int>   ("max nonlin iter steps");

  physicaltype_ = DRT::INPUT::get<INPAR::FLUID::PhysicalType>(*params_, "Physical Type");

  //--------------------------------------------------------------------
  // init some class variables (algorithm)

  numdim_ = params_->get<int>("number of velocity degrees of freedom");

  // -------------------------------------------------------------------
  // create empty system matrix --- stiffness and mass are assembled in
  // one system matrix!
  // -------------------------------------------------------------------

  // This is a first estimate for the number of non zeros in a row of
  // the matrix. Assuming a structured 3d-fluid mesh we have 27 adjacent
  // nodes with 4 dofs each. (27*4=108)
  // We do not need the exact number here, just for performance reasons
  // a 'good' estimate

  // sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,108,false,true));
  sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,500,false,true));

  // sysmat might be singular (if we have a purely Dirichlet constrained
  // problem, the pressure mode is defined only up to a constant)
  // in this case, we need a basis vector for the nullspace/kernel
  vector<DRT::Condition*> KSPcond;
  discret_->GetCondition("KrylovSpaceProjection",KSPcond);
  int numcond = KSPcond.size();
  int numfluid = 0;
  for(int icond = 0; icond < numcond; icond++)
  {
    const std::string* name = KSPcond[icond]->Get<std::string>("discretization");
    if (*name == "fluid") numfluid++;
  }
  if (numfluid == 1)
  {
    project_ = true;
    w_       = LINALG::CreateVector(*dofrowmap,true);
    c_       = LINALG::CreateVector(*dofrowmap,true);
  }
  else if (numfluid == 0)
  {
    project_ = false;
    w_       = Teuchos::null;
    c_       = Teuchos::null;
  }
  else
    dserror("Received more than one KrylovSpaceCondition for fluid field");

  // -------------------------------------------------------------------
  // create empty vectors
  // -------------------------------------------------------------------

  // Vectors passed to the element
  // -----------------------------

  // accelerations at time n+1 and n and n+alpha_M
  accnp_        = LINALG::CreateVector(*dofrowmap,true);
  accn_         = LINALG::CreateVector(*dofrowmap,true);
  accam_        = LINALG::CreateVector(*dofrowmap,true);

  // velocities and pressures at time n+1, n and n+alpha_F
  velnp_        = LINALG::CreateVector(*dofrowmap,true);
  veln_         = LINALG::CreateVector(*dofrowmap,true);
  velaf_        = LINALG::CreateVector(*dofrowmap,true);

  velnm_available_=false;

  if(predictor_=="constant_increment_predictor")
  {
    // velocities and pressure at time n-1
    velnm_        = LINALG::CreateVector(*dofrowmap,true);
  }

  // scalar at time n+1 (required only for Neumann boundary conditions)
  scaaf_       = LINALG::CreateVector(*dofrowmap,true);

  // grid displacements and velocities for the ale case
  if (alefluid_)
  {
    dispnp_     = LINALG::CreateVector(*dofrowmap,true);
    dispn_      = LINALG::CreateVector(*dofrowmap,true);
    dispnm_     = LINALG::CreateVector(*dofrowmap,true);
    gridveln_   = LINALG::CreateVector(*dofrowmap,true);
    gridvelaf_  = LINALG::CreateVector(*dofrowmap,true);
  }

  // Vectors associated to boundary conditions
  // -----------------------------------------

  // a vector of zeros to be used to enforce zero dirichlet boundary conditions
  zeros_        = LINALG::CreateVector(*dofrowmap,true);

  // object holds maps/subsets for DOFs subjected to Dirichlet BCs and otherwise
  dbcmaps_ = Teuchos::rcp(new LINALG::MapExtractor());
  {
    ParameterList eleparams;
    // other parameters needed by the elements
    eleparams.set("total time",time_);
    discret_->EvaluateDirichlet(eleparams, zeros_, Teuchos::null, Teuchos::null,
                                Teuchos::null, dbcmaps_);
    zeros_->PutScalar(0.0); // just in case of change
  }

  // the vector containing body and surface forces
  neumann_loads_= LINALG::CreateVector(*dofrowmap,true);

  // Vectors used for solution process
  // ---------------------------------

  // The residual vector --- more or less the rhs for the incremental
  // formulation!!!
  residual_     = LINALG::CreateVector(*dofrowmap,true);

  // The force vector (a copy of residual_ without Dirichlet
  // conditions applied)
  force_        = LINALG::CreateVector(*dofrowmap,true);

  // Nonlinear iteration increment vector
  increment_    = LINALG::CreateVector(*dofrowmap,true);

  // -------------------------------------------------------------------
  // get a vector layout from the discretization for a vector which only
  // contains the velocity dofs and for one vector which only contains
  // pressure degrees of freedom.
  // -------------------------------------------------------------------

  FLD::UTILS::SetupFluidSplit(*discret_,numdim_,velpressplitter_);

  // -------------------------------------------------------------------
  // initialize turbulence-statistics evaluation
  // -------------------------------------------------------------------
  ParameterList *  modelparams =&(params_->sublist("TURBULENCE MODEL"));

  // flag for special flow: currently channel flow or flow in a lid-driven cavity
  special_flow_ = modelparams->get<string>("CANONICAL_FLOW","no");

  // all averaging is done in this statistics manager
//  statisticsmanager_=Teuchos::rcp(new FLD::TurbulenceStatisticManager(*this));
  statisticsmanager_=Teuchos::null;

  if (special_flow_ != "no")
  {
    // parameters for sampling/dumping period --- used for ad-hoc
    // modification of itemax for turbulent channel flows
    samstart_  = modelparams->get<int>("SAMPLING_START",1);
  }

  // -------------------------------------------------------------------
  // initialize outflow boundary stabilization if required
  // -------------------------------------------------------------------
  ParameterList *  stabparams=&(params_->sublist("STABILIZATION"));

  // flag for potential Neumann-type outflow stabilization
  outflow_stab_ = stabparams->get<string>("OUTFLOW_STAB","yes_outstab");

  // the vector containing potential Neumann-type outflow stabilization
  if(outflow_stab_ == "yes_outstab")
    outflow_stabil_= LINALG::CreateVector(*dofrowmap,true);

  // (fine-scale) subgrid viscosity?
  fssgv_ = params_->sublist("TURBULENCE MODEL").get<string>("FSSUGRVISC","No");

  // -------------------------------------------------------------------
  // necessary only for the AVM3 approach: fine-scale solution vector
  // -------------------------------------------------------------------
  if (fssgv_ != "No") fsvelaf_ = LINALG::CreateVector(*dofrowmap,true);

  // -------------------------------------------------------------------
  // initialise the dynamic Smagorinsky model (the smoothed quantities)
  // -------------------------------------------------------------------
  if (modelparams->get<string>("TURBULENCE_APPROACH","DNS_OR_RESVMM_LES")
      ==
      "CLASSICAL_LES")
  {

    string model = modelparams->get<string>("PHYSICAL_MODEL","no_model");

    if(model=="Dynamic_Smagorinsky")
    {
      // get one instance of the dynamic Smagorinsky class
      DynSmag_=Teuchos::rcp(new FLD::DynSmagFilter(discret_            ,
                                     pbcmapmastertoslave_,
                                     *params_             ));
    }

    if (model=="Dynamic_Smagorinsky"
        ||
        model=="Smagorinsky_with_van_Driest_damping"
        ||
        model=="Smagorinsky")
    {
      if(numdim_!=3)
      {
        dserror("Smagorinsky models only available for 3d!\n");
      }
    }
  }

  // -------------------------------------------------------------------
  // check whether we have a coupling to a turbulent inflow generating
  // computation and initialise the transfer if necessary
  // -------------------------------------------------------------------
  turbulent_inflow_condition_
    = Teuchos::rcp(new TransferTurbulentInflowCondition(discret_,dbcmaps_));

  //--------------------------------------------------------------------
  // do output to screen
  this->GenAlphaEchoToScreen("print start-up info");

  // ---------------------------------------------------------------------
  // set general fluid parameter defined before
  // ---------------------------------------------------------------------
  SetElementGeneralFluidParameter();
  SetElementTurbulenceParameter();

  // end time measurement for timeloop

  tm7_ref_ = null;

  // extra discretisation for mixed/hybrid Dirichlet conditions
  MHD_evaluator_=Teuchos::rcp(new FluidMHDEvaluate(discret_));

  interface_ = Teuchos::rcp(new FLD::UTILS::MapExtractor());
  {
    interface_->Setup(*actdis);

    // build inner velocity map
    // dofs at the interface are excluded
    // we use only velocity dofs and only those without Dirichlet constraint
    const Teuchos::RCP<const LINALG::MapExtractor> dbcmaps = DirichMaps();
    std::vector<Teuchos::RCP<const Epetra_Map> > maps;
    maps.push_back(interface_->OtherMap());
    maps.push_back(dbcmaps->OtherMap());
    innervelmap_ = LINALG::MultiMapExtractor::MergeMaps(maps);

    interfaceforcen_ = Teuchos::rcp(new Epetra_Vector(*(interface_->FSICondMap())));
  }




  return;

} // FluidGenAlphaIntegration::FluidGenAlphaIntegration


/*----------------------------------------------------------------------*
 | Destructor dtor (public)                                  gammi 06/07|
 *----------------------------------------------------------------------*/
FLD::FluidGenAlphaIntegration::~FluidGenAlphaIntegration()
{
  return;
}// FluidGenAlphaIntegration::~FluidGenAlphaIntegration



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Contains the time loop for generalised alpha                         |
 |                                                           gammi 06/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//

void FLD::FluidGenAlphaIntegration::GenAlphaTimeloop()
{

  // start time measurement for timeloop
  tm2_ref_ = Teuchos::rcp(new TimeMonitor(*timedynloop_));

  bool stop_timeloop=false;
  while (stop_timeloop==false)
  {
    // -------------------------------------------------------------------
    //     preparation of time step by performing several procedures
    // -------------------------------------------------------------------
    this->GenAlphaPrepareTimeStep();

    // -------------------------------------------------------------------
    //                     solve nonlinear equation
    // -------------------------------------------------------------------
    this->DoGenAlphaPredictorCorrectorIteration();

    // -------------------------------------------------------------------
    //                         update solution
    // -------------------------------------------------------------------
    this->GenAlphaTimeUpdate();

    // -------------------------------------------------------------------
    //  statistics time sample and output of solution and statistics
    // -------------------------------------------------------------------
    this->GenAlphaStatisticsAndOutput();

    // -------------------------------------------------------------------
    // evaluate error for test flows with analytical solutions
    // -------------------------------------------------------------------
    this->EvaluateErrorComparedToAnalyticalSol();

    // -------------------------------------------------------------------
    //                    stop criterium for timeloop
    // -------------------------------------------------------------------

    if(step_>=endstep_||time_>=endtime_)
    {
	stop_timeloop=true;
    }
  }

  // end time measurement for timeloop
  tm2_ref_ = null;

  tm0_ref_ = null; // end total time measurement
  if(discret_->Comm().MyPID()==0)
  {
    cout<<"\n"<<"\n";
  }
  TimeMonitor::summarize();


  return;
}// FluidGenAlphaIntegration::GenAlphaIntegrateFromTo


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | preparation of time step by performing several procedures   vg 11/08 |
 -----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidGenAlphaIntegration::GenAlphaPrepareTimeStep()
{
  // -------------------------------------------------------------------
  //              set time dependent parameters
  // -------------------------------------------------------------------
  this->GenAlphaIncreaseTimeAndStep();

  // -------------------------------------------------------------------
  //                         out to screen
  // -------------------------------------------------------------------
  this->GenAlphaEchoToScreen("print time algorithm info");

  // -------------------------------------------------------------------
  //     predict new values for velocity and pressure
  // -------------------------------------------------------------------
  this->GenAlphaPredictNewSolutionValues();

  // -------------------------------------------------------------------
  //         evaluate dirichlet and neumann boundary conditions
  // -------------------------------------------------------------------
  // start time measurement for application of dirichlet conditions
  tm1_ref_ = Teuchos::rcp(new TimeMonitor(*timeevaldirich_));

  this->GenAlphaApplyDirichletAndNeumann();

  // end time measurement for application of dirichlet conditions
  tm1_ref_=null;

  // -------------------------------------------------------------------
  //           preparation of AVM3-based scale separation
  // -------------------------------------------------------------------
  if (step_==1 and fssgv_ != "No") this->AVM3Preparation();

  // -------------------------------------------------------------------
  //      calculate initial acceleration according to predicted
  //                  velocities and boundary values
  // -------------------------------------------------------------------
  this->GenAlphaCalcInitialAccelerations();

  return;
} // FluidGenAlphaIntegration::GenAlphaPrepareTimeStep


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Iterative procedure to solve the nonlinear problem resulting from    |
 | the time discrete version                                            |
 |                                                           gammi 06/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//

void FLD::FluidGenAlphaIntegration::DoGenAlphaPredictorCorrectorIteration(
  )
{
  double            tcpu     ;

  // reset iteration counter, initialize max iteration counter
  itenum_  = 0;
  // currently default for turbulent channel flow: only one iteration before sampling
  if (special_flow_ == "channel_flow_of_height_2" && step_<samstart_ )
    itemax_  = 1;
  else itemax_  = params_->get<int>   ("max nonlin iter steps");

  // stop nonlinear iteration when both increment-norms are below this
  // bound
  ittol_     =params_->get<double>("tolerance for nonlin iter");

  // set
  if (params_->get<string>("CONVCHECK","L_2_norm")
      ==
      "L_2_norm_without_residual_at_itemax")
  {
    skiplastelecall_=true;
  }
  else
  {
    skiplastelecall_=false;
  }

  // start time measurement for nonlinear iteration
  tm6_ref_ = Teuchos::rcp(new TimeMonitor(*timenlnloop_));

  // -------------------------------------------------------------------
  //  Evaluate acceleration and velocity at the intermediate time level
  //                     n+alpha_M and n+alpha_F
  //
  //                             -> (0)
  // -------------------------------------------------------------------
  // start time measurement for nonlinear update
  tm9_ref_ = Teuchos::rcp(new TimeMonitor(*timenonlinup_));

  this->GenAlphaComputeIntermediateSol();

  // time measurement --- stop TimeMonitor tm9
  tm9_ref_        = null;

  //--------------------------------------------------------------------
  // do output to screen
  this->GenAlphaEchoToScreen("print first nonlinear iteration info");

  // -------------------------------------------------------------------
  // call elements to calculate residual and matrix for first iteration
  // -------------------------------------------------------------------

  this->GenAlphaAssembleResidualAndMatrix();

  {

    // in the case of local systems we have to rotate back for the
    // convergence check ...
    if (locsysman_ != Teuchos::null)
    {
      locsysman_->RotateLocalToGlobal(residual_);
    }

    Teuchos::RCP<Epetra_Vector> onlyvel = velpressplitter_.ExtractOtherVector(residual_);
    Teuchos::RCP<Epetra_Vector> onlypre = velpressplitter_.ExtractCondVector (residual_);

    // extract velocity and pressure residuals from rhs vector
    LINALG::Export(*residual_,*onlyvel);
    LINALG::Export(*residual_,*onlypre);

    onlypre->Norm2(&L2preresnorm_);
    onlyvel->Norm2(&L2velresnorm_);

    // in the case of local systems we have to forward again ...
    if (locsysman_ != Teuchos::null)
    {
      locsysman_->RotateGlobalToLocal(residual_);
    }
  }

  //--------------------------------------------------------------------
  // do output to screen
  this->GenAlphaEchoToScreen("print residual norms before first iteration");
  this->GenAlphaEchoToScreen("print timing of element calls");

  // -------------------------------------------------------------------
  // now to the iteration
  // -------------------------------------------------------------------
  bool stopnonliniter = false;
  double badestnlnnorm = 1000.0;
  while (stopnonliniter==false)
  {
    itenum_++;

    // -------------------------------------------------------------------
    // solve for increments
    // -------------------------------------------------------------------
    // start time measurement for solver call
    tm5_ref_ = Teuchos::rcp(new TimeMonitor(*timesolver_));

    // get cpu time
    tcpu=Teuchos::Time::wallTime();

    this->GenAlphaCalcIncrement(badestnlnnorm);

    // end time measurement for application of dirichlet conditions
    tm5_ref_=null;
    dtsolve_=Teuchos::Time::wallTime()-tcpu;

    // start time measurement for nonlinear update
    tm9_ref_ = Teuchos::rcp(new TimeMonitor(*timenonlinup_));

    // -------------------------------------------------------------------
    // update estimates by incremental solution
    // -------------------------------------------------------------------
    this->GenAlphaNonlinearUpdate();


    // -------------------------------------------------------------------
    //  Evaluate acceleration and velocity at the intermediate time level
    //                     n+alpha_M and n+alpha_F
    //
    //                          (i)->(i+1)
    // -------------------------------------------------------------------
    this->GenAlphaComputeIntermediateSol();

    // time measurement --- stop TimeMonitor tm9
    tm9_ref_        = null;

    // -------------------------------------------------------------------
    // call elements to calculate residual for convergence check and
    // matrix for the next step
    // skip if the residual check in the last iteration is suppressed
    // -------------------------------------------------------------------
    if(!(itenum_ == itemax_) || !(skiplastelecall_))
    {
      this->GenAlphaAssembleResidualAndMatrix();
    }

    // -------------------------------------------------------------------
    // do convergence check
    // -------------------------------------------------------------------
    stopnonliniter=this->GenAlphaNonlinearConvergenceCheck(badestnlnnorm);
  }

  // end time measurement for nonlinear iteration
  tm6_ref_ = null;

  return;
}// FluidGenAlphaIntegration::DoGenAlphaPredictorCorrectorIteration


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Predict velocity and pressure of the new timestep. Up to now, we    |
 |  use a constant predictor for the velocity and the pressure or one   |
 |  of the following: zero-acceleration, constant-acceleration          |
 |  predictor, constant-increment predictor                             |
 |                                                                      |
 |  Remark: For Dirichlet nodes, no matter what was set here, velnp     |
 |          will be overwritten by the prescribed value. The            |
 |          accelerations are calculated after these Dirichlet values   |
 |          have been set.                                              |
 |                                                                      |
 |                                                           gammi 06/07|
 -----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidGenAlphaIntegration::GenAlphaPredictNewSolutionValues()
{

  // first guess is the same for all predictors
  velnp_->Update(1.0,*veln_ ,0.0);

  if(predictor_=="steady_state_predictor")
  {
    // steady state predictor
    //
    //       n+1    n
    //      u    = u
    //       (0)
    //
    //  and
    //
    //       n+1    n
    //      p    = p
    //       (0)

    // there is nothing more to do
  }
  else if(predictor_=="zero_acceleration_predictor")
  {
    // zero acceleration predictor
    //
    //       n+1    n                   n
    //      u    = u  + (1-gamma)*dt*acc
    //       (0)
    //
    //  and
    //
    //       n+1    n
    //      p    = p
    //       (0)
    //

    // split between acceleration and pressure
    Teuchos::RCP<Epetra_Vector> inc = velpressplitter_.ExtractOtherVector(accn_);
    inc->Scale((1.0-gamma_)*dt_);

    velpressplitter_.AddOtherVector(inc,velnp_);
  }
  else if(predictor_=="constant_acceleration_predictor")
  {
    // constant acceleration predictor
    //
    //       n+1    n         n
    //      u    = u  + dt*acc
    //       (0)
    //
    //  and
    //
    //       n+1    n
    //      p    = p
    //       (0)
    //

    Teuchos::RCP<Epetra_Vector> inc = velpressplitter_.ExtractOtherVector(accn_);
    inc->Scale(dt_);

    velpressplitter_.AddOtherVector(inc,velnp_);
  }
  else if(predictor_=="constant_increment_predictor")
  {
    // constant increment predictor
    //
    //       n+1      n    n-1
    //      u    = 2*u  - u
    //       (0)
    //
    //  and
    //
    //       n+1    n
    //      p    = p
    //       (0)
    //

    if(velnm_available_)
    {
      Teuchos::RCP<Epetra_Vector> onlyveln  = velpressplitter_.ExtractOtherVector(veln_ );
      Teuchos::RCP<Epetra_Vector> onlyvelnm = velpressplitter_.ExtractOtherVector(velnm_);
      onlyvelnm->Scale(-1.0);

      velpressplitter_.AddOtherVector(onlyveln ,velnp_);
      velpressplitter_.AddOtherVector(onlyvelnm,velnp_);
    }
  }

  return;
} // FluidGenAlphaIntegration::GenAlphaPredictNewSolutionValues


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Apply Dirichlet boundary conditions to velocity vector and          |
 |  calculate accelerations according to prescribed Dirichlet values.   |
 |  Apply surface Neumann conditions.                                   |
 |                                                           gammi 06/07|
 -----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidGenAlphaIntegration::GenAlphaApplyDirichletAndNeumann()
{
  // --------------------------------------------------
  // apply Dirichlet conditions to velnp

  // in the case of local systems we have to rotate forward ...
  if (locsysman_ != Teuchos::null)
  {
    locsysman_->RotateGlobalToLocal(velnp_);
  }

  ParameterList eleparams;

  // total time required for Dirichlet conditions
  eleparams.set("total time",time_);
  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("velnp",velnp_);

  // predicted dirichlet values
  // velnp then also holds prescribed new dirichlet values
  discret_->EvaluateDirichlet(eleparams,velnp_,null,null,null,dbcmaps_);
  discret_->ClearState();

  // Transfer of boundary data if necessary
  turbulent_inflow_condition_->Transfer(veln_,velnp_,time_);

  // in the case of local systems we have to rotate back into global
  // Cartesian frame
  if (locsysman_ != Teuchos::null)
  {
    locsysman_->RotateLocalToGlobal(velnp_);
  }

  // -------------------------------------------------------------------
  //  Set time parameter for element call
  // -------------------------------------------------------------------
  SetElementTimeParameter();

  // --------------------------------------------------
  // evaluate Neumann conditions
  // ---------------------------------------------------

  neumann_loads_->PutScalar(0.0);
  discret_->SetState("scaaf",scaaf_);
  discret_->EvaluateNeumann(eleparams,*neumann_loads_);
  discret_->ClearState();

  return;
} // FluidGenAlphaIntegration::GenAlphaApplyDirichletAndNeumann


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  calculate accelerations according to prescribed Dirichlet values    |
 |  and predicted solution values.                                      |
 |                                                           gammi 06/07|
 -----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidGenAlphaIntegration::GenAlphaCalcInitialAccelerations()
{
  // --------------------------------------------------
  // adjust accnp according to Dirichlet values of velnp
  //
  //                                  n+1     n
  //                               vel   - vel
  //       n+1      n  gamma-1.0      (0)
  //    acc    = acc * --------- + ------------
  //       (0)           gamma      gamma * dt
  //

  accnp_->Update(1.0,*velnp_,-1.0,*veln_,0.0);
  accnp_->Update((gamma_-1.0)/gamma_,*accn_,1.0/(gamma_*dt_));

  return;
} // FluidGenAlphaIntegration::GenAlphaCalcInitialAccelerations

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Evaluate acceleration and velocity at the intermediate time level    |
 | n+alpha_M and n+alpha_F                                              |
 |                                                           gammi 06/07|
 -----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidGenAlphaIntegration::GenAlphaComputeIntermediateSol()
{
  //       n+alphaM                n+1                      n
  //    acc         = alpha_M * acc     + (1-alpha_M) *  acc
  //       (i)                     (i)

  accam_->Update((alphaM_),*accnp_,(1.0-alphaM_),*accn_,0.0);

  //       n+alphaF              n+1                   n
  //      u         = alpha_F * u     + (1-alpha_F) * u
  //       (i)                   (i)

  velaf_->Update((alphaF_),*velnp_,(1.0-alphaF_),*veln_,0.0);

  return;
} // FluidGenAlphaIntegration::GenAlphaComputeIntermediateSol


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Update the solution after convergence of the nonlinear iteration.   |
 |  Current solution becomes old solution of next timestep.             |
 |                                                           gammi 06/07|
 -----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidGenAlphaIntegration::GenAlphaTimeUpdate()
{
  //--------------------------------------------------
  // solution of this step becomes most recent solution of the last step

  // store old velos for constant increment predictor
  if(predictor_=="constant_increment_predictor")
  {
    velnm_available_=true;
    // velocities and pressure at time n-1
    velnm_->Update(1.0,*veln_ ,0.0);
  }
  // for velocities and pressure
  veln_->Update(1.0,*velnp_ ,0.0);
  // for the accelerations
  accn_->Update(1.0,*accnp_ ,0.0);

  if(params_->sublist("STABILIZATION").get<string>("TDS")=="time_dependent")
  {
    // create the parameters for the discretization
    ParameterList eleparams;
    // action for elements
    eleparams.set<int>("action",FLD::calc_fluid_genalpha_update_for_subscales);

    // update time paramters
    eleparams.set("gamma"  ,gamma_ );
    eleparams.set("dt"     ,dt_    );

    // call loop over elements
    discret_->Evaluate(eleparams,null,null,null,null,null);
  }

  if (alefluid_)
  {

    /*
       2nd order finite difference approximation of the new grid velocity

                    (equivalent to what is done for BDF2)

                           /                      \
        n+1        1      |     n+1       n    n-1 |
     u_G     = -------- * | 3* d   - 4 * d  + d    |
                2 * dt    |                        |
                           \                      /
    */
    gridveln_->Update(1.5/dt_, *dispnp_, -2.0/dt_, *dispn_, 0.0);
    gridveln_->Update(0.5/dt_, *dispnm_, 1.0);

    //    n+1         n
    //   d      ---> d
    //
    dispnm_  ->Update(1.0,*dispn_,0.0);
    dispn_   ->Update(1.0,*dispnp_,0.0);
  }


  return;
} // FluidGenAlphaIntegration::GenAlphaTimeUpdate


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | statistics time sample and output of solution and statistics         |
 |                                                             vg 11/08 |
 -----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidGenAlphaIntegration::GenAlphaStatisticsAndOutput()
{
  // time measurement --- start TimeMonitor tm8
  tm8_ref_ = Teuchos::rcp(new TimeMonitor(*timeout_ ));

  // -------------------------------------------------------------------
  //   add calculated velocity to mean value calculation (statistics)
  // -------------------------------------------------------------------
//  statisticsmanager_->DoTimeSample(step_,time_,0.0);

  // -------------------------------------------------------------------
  //   calculate lift and drag values
  // -------------------------------------------------------------------
  this->LiftDrag();

  // -------------------------------------------------------------------
  //                         output of solution
  // -------------------------------------------------------------------
  this->GenAlphaOutput();

  // time measurement --- stop TimeMonitor tm8
  tm8_ref_        = null;

  return;
} // FluidGenAlphaIntegration::StatisticsAndOutput


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Write solution to file                                   gammi 06/07|
 -----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidGenAlphaIntegration::GenAlphaOutput()
{
#if 0
  {
    const Epetra_Map* dofrowmap = discret_->DofRowMap();

    vector<double> my_y;
    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = discret_->lRowNode(lnodeid);

      if (abs(lnode->X()[0])<1e-6 && abs(lnode->X()[2])<1e-6)
      {
        my_y.push_back(lnode->X()[1]);
      }

      // the set of degrees of freedom associated with the node
      vector<int> nodedofset = discret_->Dof(lnode);
    }

    vector<int> num_y(discret_->Comm().NumProc(),0);
    num_y[myrank_]=my_y.size();
    vector<int> num_y_all(discret_->Comm().NumProc(),0);
    discret_->Comm().SumAll(&num_y[0],&num_y_all[0],(discret_->Comm().NumProc()));

    int n            =0;
    int n_lower_procs=0;

    for(int rr=0;rr<(discret_->Comm().NumProc());++rr)
    {
      n+=num_y_all[rr];
    }

    for(int rr=0;rr<myrank_;++rr)
    {
      n_lower_procs+=num_y_all[rr];
    }

    vector<int>    rank     (n,0);
    vector<int>    rank_loc (n,0);

    vector<int>    id       (n,0);
    vector<int>    id_loc   (n,0);

    vector<double> u    (n,0.0);
    vector<double> u_loc(n,0.0);

    vector<double> v    (n,0.0);
    vector<double> v_loc(n,0.0);

    vector<double> w    (n,0.0);
    vector<double> w_loc(n,0.0);

    vector<double> p    (n,0.0);
    vector<double> p_loc(n,0.0);

    vector<double> y    (n,0.0);
    vector<double> y_loc(n,0.0);


    int count=n_lower_procs;

    int lid=-1;

    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = discret_->lRowNode(lnodeid);

      if (abs(lnode->X()[0])<1e-6 && abs(lnode->X()[2])<1e-6)
      {
        id_loc  [count]=lnode->Id();
        rank_loc[count]=myrank_;

        y_loc[count]=lnode->X()[1];

        // the set of degrees of freedom associated with the node
        vector<int> nodedofset = discret_->Dof(lnode);

        lid = dofrowmap->LID(nodedofset[0]);
        u_loc[count]=(*velnp_)[lid];

        lid = dofrowmap->LID(nodedofset[1]);
        v_loc[count]=(*velnp_)[lid];

        lid = dofrowmap->LID(nodedofset[2]);
        w_loc[count]=(*velnp_)[lid];

        lid = dofrowmap->LID(nodedofset[3]);
        p_loc[count]=(*velnp_)[lid];

        ++count;
      }
    }


    discret_->Comm().SumAll(&rank_loc[0],&rank[0],n);
    discret_->Comm().SumAll(&id_loc  [0],&id[0]  ,n);

    discret_->Comm().SumAll(&y_loc[0],&y[0],n);
    discret_->Comm().SumAll(&u_loc[0],&u[0],n);
    discret_->Comm().SumAll(&v_loc[0],&v[0],n);
    discret_->Comm().SumAll(&w_loc[0],&w[0],n);
    discret_->Comm().SumAll(&p_loc[0],&p[0],n);

    if(myrank_==0)
    {
      for(int rr=0;rr<n;++rr)
      {
        printf("NR  %3d %5d %18.10e  %13.5e  %13.5e  %13.5e  %13.5e  %13.5e\n",rank[rr],id[rr],time_,y[rr],u[rr],v[rr],w[rr],p[rr]);
      }
    }
  }
#endif

  //-------------------------------------------- output of solution
  if (step_%upres_ == 0)  //write solution
  {
    output_->NewStep    (step_,time_);

    output_->WriteVector("velnp"   , velnp_);

    // output real pressure
    Teuchos::RCP<Epetra_Vector> pressure = velpressplitter_.ExtractCondVector(velnp_);
    pressure->Scale(density_);
    output_->WriteVector("pressure", pressure);

    if (alefluid_)
    {
      output_->WriteVector("dispnp", dispnp_);
    }


    if(params_->sublist("TURBULENCE MODEL").get<string>("PHYSICAL_MODEL","no_model")
       ==
       "Dynamic_Smagorinsky")
    {
      const Epetra_Map* dofrowmap = discret_->DofRowMap();

      RefCountPtr<Epetra_Vector> filteredvel = LINALG::CreateVector(*dofrowmap,true);

      DynSmag_->OutputofAveragedVel(filteredvel);
      output_->WriteVector("filteredvel",filteredvel);
    }


    //only perform stress calculation when output is needed
    if (writestresses_)
    {
     RefCountPtr<Epetra_Vector> traction = CalcStresses();
     output_->WriteVector("traction",traction);
    }

    // write domain decomposition for visualization (only once!)
    if (step_==upres_)
     output_->WriteElementData();

    // dumping of turbulence statistics if required
//    statisticsmanager_->DoOutput(*output_,step_);

    // do restart if we have to
    if (step_%uprestart_ == 0)
    {
      output_->WriteVector("veln",  veln_ );
      output_->WriteVector("accn",  accn_ );

      if (alefluid_)
      {
        output_->WriteVector("dispn"   ,dispn_   );
        output_->WriteVector("dispnm"  ,dispnm_  );
        output_->WriteVector("gridveln",gridveln_);
      }

      // write mesh in each restart step --- the elements are required since
      // they contain history variables (the time dependent subscales)
      // But never do this for step 0 (visualization of initial field) since
      // it would lead to writing the mesh twice for step 0
      // (FluidBaseAlgorithm already wrote the mesh) -> HDF5 writer will claim!
      if (step_!=0)
        output_->WriteMesh(step_,time_);
    }

  }
  // write restart also when uprestart_ is not a integer multiple of upres_
  else if (step_%uprestart_ == 0)
  {
    output_->NewStep    (step_,time_);

    output_->WriteVector("velnp", velnp_);
    output_->WriteVector("veln" , veln_ );
    output_->WriteVector("accn" , accn_ );

    if (alefluid_)
    {
      output_->WriteVector("dispnp" , dispnp_  );
      output_->WriteVector("dispn"   ,dispn_   );
      output_->WriteVector("dispnm"  ,dispnm_  );
      output_->WriteVector("gridveln",gridveln_);
    }

    //only perform stress calculation when output is needed
    if (writestresses_)
    {
     RefCountPtr<Epetra_Vector> traction = CalcStresses();
     output_->WriteVector("traction",traction);
    }

    // dumping of turbulence statistics if required
//    statisticsmanager_->DoOutput(*output_,step_);

    // write mesh in each restart step --- the elements are required since
    // they contain history variables (the time dependent subscales)
    // But never do this for step 0 (visualization of initial field) since
    // it would lead to writing the mesh twice for step 0
    // (FluidBaseAlgorithm already wrote the mesh) -> HDF5 writer will claim!
    if (step_!=0)
      output_->WriteMesh(step_,time_);
  }

  return;
} // FluidGenAlphaIntegration::GenAlphaOutput


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Assemble residual and system matrix. Dirichlet conditions applied in |
 | here, the true residual is stored in force_.                         |
 |                                                           gammi 07/07|
 -----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidGenAlphaIntegration::GenAlphaAssembleResidualAndMatrix()
{

  // -------------------------------------------------------------------
  // Filter velocity for dynamic Smagorinsky model --- this provides
  // the necessary dynamic constant
  // -------------------------------------------------------------------
  if (params_->sublist("TURBULENCE MODEL").get<string>("TURBULENCE_APPROACH","DNS_OR_RESVMM_LES")
      ==
      "CLASSICAL_LES")
  {
    if(params_->sublist("TURBULENCE MODEL").get<string>("PHYSICAL_MODEL","no_model")
       ==
       "Dynamic_Smagorinsky"
      )
    {
      this->ApplyFilterForDynamicComputationOfCs();
    }
  }

  // get cpu time
  double tcpu=Teuchos::Time::wallTime();

  // -------------------------------------------------------------------
  // call elements to calculate residual and matrix
  // -------------------------------------------------------------------
  // zero out the stiffness matrix
  // we keep the sparsity pattern throughout the calculation for
  // performance reasons
  // start time measurement for generation of sparsity pattern
  {
    RefCountPtr<TimeMonitor> timesparsitypattern_ref_
      =
      Teuchos::rcp(new TimeMonitor(*timesparsitypattern_));

    sysmat_->Zero();

    if (shapederivatives_ != Teuchos::null)
    {
      shapederivatives_->Zero();
    }

    timesparsitypattern_ref_=null;
  }

  // Neumann loads to residual
  residual_->Update(1.0,*neumann_loads_,0.0);

  // start time measurement for element call
  tm3_ref_ = Teuchos::rcp(new TimeMonitor(*timeeleloop_));

  // add stabilization term at Neumann outflow boundary if required
  if(outflow_stab_ == "yes_outstab")
  {
    // create the parameters for the discretization
    ParameterList condparams;

    discret_->ClearState();
    discret_->SetState("velaf",velaf_);
    discret_->SetState("scaaf",scaaf_);
    condparams.set("thsl",alphaF_*gamma_*dt_);
    condparams.set("rhs time factor",1.0);
    condparams.set<int>("action",FLD::calc_Neumann_inflow);
    condparams.set<int>("Physical Type",physicaltype_);
    condparams.set("using generalized-alpha time integration",true);

    std::string condstring("FluidNeumannInflow");
    discret_->EvaluateCondition(condparams,
                                sysmat_,
                                Teuchos::null,
                                residual_,
                                Teuchos::null,
                                Teuchos::null,
                                condstring);

    discret_->ClearState();
  }

  // create the parameters for the discretization
  ParameterList eleparams;

  // action for elements
  eleparams.set<int>("action",FLD::calc_fluid_genalpha_sysmat_and_residual);
  // choose what to assemble
  if (itenum_ == itemax_)
  {
    eleparams.set("assemble matrix 1",false);
  }
  else
  {
    eleparams.set("assemble matrix 1",true);
  }
  eleparams.set("assemble matrix 2",false);
  eleparams.set("assemble vector 1",true);
  eleparams.set("assemble vector 2",false);
  eleparams.set("assemble vector 3",false);

  // other parameters that might be needed by the elements
  {
    ParameterList& timelist = eleparams.sublist("time integration parameters");

    timelist.set("alpha_M",alphaM_);
    timelist.set("alpha_F",alphaF_);
    timelist.set("gamma"  ,gamma_ );
    timelist.set("time"   ,time_  );
    timelist.set("dt"     ,dt_    );
  }

  // do not compute the element matrix if itmax is reached
  // in this case, only the residual is required for the convergence check
  if (itenum_<itemax_)
  {
    eleparams.set("compute element matrix",true);
  }
  else
  {
    eleparams.set("compute element matrix",false);
  }

  // parameters for nonlinear treatment (linearisation) and low-Mach-number solver
  eleparams.set<int>("Linearisation",newton_);

  // parameters for usage of conservative/convective form
  eleparams.set("CONVFORM",params_->get<string>("form of convective term"));

  // parameters for stabilisation
  {
    eleparams.sublist("STABILIZATION")    = params_->sublist("STABILIZATION");
  }

  // parameters for a turbulence model
  {
    eleparams.sublist("TURBULENCE MODEL") = params_->sublist("TURBULENCE MODEL");
  }

  // (fine-scale) subgrid viscosity flag
  eleparams.set("fs subgrid viscosity",fssgv_);

  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("u and p (n+1      ,trial)",velnp_ );
  discret_->SetState("u and p (n+alpha_F,trial)",velaf_ );
  discret_->SetState("acc     (n+alpha_M,trial)",accam_ );

  if (alefluid_)
  {
    discret_->SetState("dispnp"    , dispnp_   );
    discret_->SetState("gridvelaf" , gridvelaf_);
  }

  //----------------------------------------------------------------------
  // decide whether VM3-based solution approach or standard approach
  //----------------------------------------------------------------------
  if (fssgv_ != "No") AVM3Separation();

  //----------------------------------------------------------------------
  // call loop over elements
  //----------------------------------------------------------------------
  // do not assemble the elemetn matrix if itmax is reached
  // in this case, only the residual is required for the convergence check
  if (itenum_<itemax_)
  {
    discret_->Evaluate(eleparams,sysmat_,shapederivatives_,residual_,Teuchos::null,Teuchos::null);
  }
  else
  {
    discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,residual_,Teuchos::null,Teuchos::null);
  }

  discret_->ClearState();

  // get density
  density_ = eleparams.get<double>("density");

  //----------------------------------------------------------------------
  // extended statistics (plane average of Cs) for dynamic
  // Smagorinsky model --- communication part, store values
  //----------------------------------------------------------------------
//  statisticsmanager_->StoreElementValues(step_);

  //----------------------------------------------------------------------
  // remember force vector for stress computation
  //----------------------------------------------------------------------
  if(abs(density_)<1e-16)
  {
    dserror("zero density not allowed in dynamic computations\n");
  }
  force_->Update(density_,*residual_,0.0);

  // rotate forces from global to local co-ordinate system
  if (locsysman_ != Teuchos::null)
  {
    locsysman_->RotateGlobalToLocal(force_);
  }


  //----------------------------------------------------------------------
  // apply weak Dirichlet boundary conditions to sysmat_ and residual_
  //----------------------------------------------------------------------
  {
    // vector containing weak dirichlet loads

    RefCountPtr<Epetra_Vector> wdbcloads = LINALG::CreateVector(*(discret_->DofRowMap()),true);

    ParameterList weakdbcparams;

    // set action for elements
    weakdbcparams.set<int>("action"    ,FLD::enforce_weak_dbc);
    weakdbcparams.set("gdt"       ,gamma_*dt_        );
    weakdbcparams.set("afgdt"     ,alphaF_*gamma_*dt_);
    weakdbcparams.set("total time",time_             );
    weakdbcparams.set("using p^{n+1} generalized-alpha time integration",true);

    // set the only required state vectors
    discret_->SetState("velaf",velaf_);
    discret_->SetState("velnp",velnp_);
    if (alefluid_)
    {
      discret_->SetState("dispnp"    , dispnp_   );
      discret_->SetState("gridvelaf" , gridvelaf_);
    }

    // evaluate all line weak Dirichlet boundary conditions
    discret_->EvaluateConditionUsingParentData
      (weakdbcparams      ,
       sysmat_            ,
       Teuchos::null      ,
       wdbcloads          ,
       Teuchos::null      ,
       Teuchos::null      ,
       "LineWeakDirichlet");

    // evaluate all surface weak Dirichlet boundary conditions
    discret_->EvaluateConditionUsingParentData
      (weakdbcparams      ,
       sysmat_            ,
       Teuchos::null      ,
       wdbcloads          ,
       Teuchos::null      ,
       Teuchos::null      ,
       "SurfaceWeakDirichlet");

    // clear state
    discret_->ClearState();

    // update the residual
    residual_->Update(1.0,*wdbcloads,1.0);
  }


  //----------------------------------------------------------------------
  // apply mixed/hybrid Dirichlet boundary conditions
  //----------------------------------------------------------------------
  vector<DRT::Condition*> MHDcnd;
  discret_->GetCondition("SurfaceMixHybDirichlet",MHDcnd);

  if(MHDcnd.size()!=0)
  {

    ParameterList mhdbcparams;

    // set action for elements
    mhdbcparams.set<int>("action"    ,FLD::mixed_hybrid_dbc);
    mhdbcparams.set("alpha_F",alphaF_);
    mhdbcparams.set("gamma"  ,gamma_ );
    mhdbcparams.set("dt"     ,dt_    );
    mhdbcparams.set("total time",time_                 );
    mhdbcparams.set("using p^{n+1} generalized-alpha time integration",true);

    // set the required state vectors
    discret_->SetState("velaf"    ,velaf_);
    discret_->SetState("velnp",velnp_);

    discret_->EvaluateConditionUsingParentData
      (mhdbcparams          ,
       sysmat_              ,
       Teuchos::null        ,
       residual_            ,
       Teuchos::null        ,
       Teuchos::null        ,
       "LineMixHybDirichlet");

    // clear state
    discret_->ClearState();

    bool doold=false;

    if(doold)
    {
      discret_->EvaluateConditionUsingParentData
        (mhdbcparams          ,
         sysmat_              ,
         Teuchos::null        ,
         residual_            ,
         Teuchos::null        ,
         Teuchos::null        ,
         "SurfaceMixHybDirichlet");
    }
    else
    {
      MHD_evaluator_->BoundaryElementLoop(
        mhdbcparams   ,
        velaf_        ,
        velnp_        ,
        residual_     ,
        SystemMatrix());
    }
  }

  //----------------------------------------------------------------------
  // apply consistent outflow boundary condition for conservative element
  //     formulation (expressions arising from partial integration)
  //----------------------------------------------------------------------
  if(params_->get<string>("form of convective term")=="conservative")
  {
    ParameterList apply_cons_params;
    // set action for elements
    apply_cons_params.set<int>("action",FLD::conservative_outflow_bc);
    apply_cons_params.set("timefac_mat" ,alphaF_*gamma_*dt_);
    apply_cons_params.set("timefac_rhs" ,               1.0);

    // set required state vectors
    discret_->ClearState();

    discret_->SetState("u and p (trial)",velaf_);
    if (alefluid_)
    {
      discret_->SetState("dispnp", dispnp_);
    }

    // call loop over elements
    discret_->EvaluateCondition(apply_cons_params                      ,
                                sysmat_                                ,
                                Teuchos::null                          ,
                                residual_                              ,
                                Teuchos::null                          ,
                                Teuchos::null                          ,
                                "SurfaceConservativeOutflowConsistency");

    discret_->ClearState();
  }

  // end time measurement for element call
  tm3_ref_=null;

  // start time measurement for generation of sparsity pattern
  {
    RefCountPtr<TimeMonitor> timesparsitypattern_ref_ = Teuchos::rcp(new TimeMonitor(*timesparsitypattern_));
    // finalize the system matrix
    sysmat_->Complete();

    if (shapederivatives_ != Teuchos::null)
    {
      shapederivatives_->Complete();
      // apply Dirichlet conditions to a non-diagonal matrix
      // (The Dirichlet rows will become all zero, no diagonal one.)
      shapederivatives_->ApplyDirichlet(*(dbcmaps_->CondMap()),false);
    }

    timesparsitypattern_ref_ = null;
  }


  // -------------------------------------------------------------------
  // Apply strong Dirichlet boundary conditions to system of equations
  // residuals are supposed to be zero at boundary conditions
  // -------------------------------------------------------------------
  // start time measurement for application of dirichlet conditions
  tm4_ref_ = Teuchos::rcp(new TimeMonitor(*timeapplydirich_));

  {
    // cast EpetraOperator sysmat_ to a LINALG::SparseMatrix in order to
    // use the matrix-multiply capabilities in the locsysmanager
    LINALG::SparseMatrix* A = dynamic_cast<LINALG::SparseMatrix*>(sysmat_.get());
    if (A!=NULL)
    {
      // transform to local co-ordinate systems
      if (locsysman_!=Teuchos::null)
      {
        locsysman_->RotateGlobalToLocal(Teuchos::rcp(A,false),residual_);
      }

      // apply the dirichlet conditions to the (rotated) system
      zeros_->PutScalar(0.0);
      {
        LINALG::ApplyDirichlettoSystem(Teuchos::rcp(A,false)           ,
                                       increment_             ,
                                       residual_              ,
                                       GetLocSysTrafo()       ,
                                       zeros_                 ,
                                       *(dbcmaps_->CondMap()));
      }
    }
    else
    {
      // transform to local co-ordinate systems
      if (locsysman_!=Teuchos::null)
      {
        dserror("expecting LINALG::SparseMatrix as a SparseOperator. Cannot use locsys.\n");
      }

      // apply the dirichlet conditions to the (rotated) system
      zeros_->PutScalar(0.0);
      {
        LINALG::ApplyDirichlettoSystem(sysmat_                ,
                                       increment_             ,
                                       residual_              ,
                                       zeros_                 ,
                                       *(dbcmaps_->CondMap()));
      }
    }
  }

  // -------------------------------------------------------------------
  // For purely Dirichlet bounded problems, compute mean basis function
  // weight vector and define nullspace for matrix
  // -------------------------------------------------------------------
  if (project_)
  {
    DRT::Condition* KSPcond=discret_->GetCondition("KrylovSpaceProjection");

    // in this case, we want to project out some zero pressure modes
    const string* definition = KSPcond->Get<string>("weight vector definition");

    if(*definition == "pointvalues")
    {
      // zero w and c
      w_->PutScalar(0.0);
      c_->PutScalar(0.0);

      // get pressure
      const vector<double>* mode = KSPcond->Get<vector<double> >("mode");

      for(int rr=0;rr<numdim_;++rr)
      {
        if(abs((*mode)[rr]>1e-14))
        {
          dserror("expecting only an undetermined pressure");
        }
      }

      int predof = numdim_;

      Teuchos::RCP<Epetra_Vector> presmode = velpressplitter_.ExtractCondVector(*w_);

      presmode->PutScalar((*mode)[predof]);

      /* export to vector to normalize against
      //
      // Note that in the case of definition pointvalue based,
      // the average pressure will vanish in a pointwise sense
      //
      //    +---+
      //     \
      //      +   p_i  = 0
      //     /
      //    +---+
      */
#if 0
      if(myrank_==0)
      {
        int lid =0;
        double one=1.0;

        presmode->PutScalar(0.0);
        int err = presmode->ReplaceMyValues(1,&one,&lid);
      }
#endif

      LINALG::Export(*presmode,*w_);

      // export to vector of ones
      presmode->PutScalar(1.0);

      LINALG::Export(*presmode,*c_);
    }
    else if(*definition == "integration")
    {
      // zero w and c
      w_->PutScalar(0.0);
      c_->PutScalar(0.0);

      ParameterList mode_params;

      // set action for elements
      mode_params.set<int>("action",FLD::integrate_shape);

      if (alefluid_)
      {
        discret_->SetState("dispnp",dispnp_);
      }

      /* evaluate KrylovSpaceProjection condition in order to get
      // integrated nodal basis functions w_
      // Note that in the case of definition integration based,
      // the average pressure will vanish in an integral sense
      //
      //                    /              /                      /
      //   /    \          |              |  /          \        |  /    \
      //  | w_*p | = p_i * | N_i(x) dx =  | | N_i(x)*p_i | dx =  | | p(x) | dx = 0
      //   \    /          |              |  \          /        |  \    /
      //  	          /  	         /  		        /
      */

      discret_->EvaluateCondition
      (mode_params        ,
       Teuchos::null      ,
       Teuchos::null      ,
       w_                 ,
       Teuchos::null      ,
       Teuchos::null      ,
       "KrylovSpaceProjection");

      // get pressure
      const vector<double>* mode = KSPcond->Get<vector<double> >("mode");

      for(int rr=0;rr<numdim_;++rr)
      {
        if(abs((*mode)[rr]>1e-14))
        {
          dserror("expecting only an undetermined pressure");
        }
      }

      Teuchos::RCP<Epetra_Vector> presmode = velpressplitter_.ExtractCondVector(*w_);

      // export to vector of ones
      presmode->PutScalar(1.0);
      LINALG::Export(*presmode,*c_);
    }
    else
    {
      dserror("unknown definition of weight vector w for restriction of Krylov space");
    }
  }

  // end time measurement for application of dirichlet conditions
  tm4_ref_=null;

  // end measurement element call
  dtele_=Teuchos::Time::wallTime()-tcpu;

  return;
} // FluidGenAlphaIntegration::GenAlphaAssembleResidualAndMatrix

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Solve linear problem                                      gammi 06/07|
 -----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidGenAlphaIntegration::GenAlphaCalcIncrement(const double nlnres)
{
  bool   isadapttol    = params_->get<bool>  ("ADAPTCONV"                ,true  );
  double adaptolbetter = params_->get<double>("ADAPTCONV_BETTER"         ,0.01  );
  double nlntolsol     = params_->get<double>("tolerance for nonlin iter",1.e-10);

  //--------------------------- adapt tolerance  in the convergence limit
  if (isadapttol && itenum_>1)
  {
    solver_->AdaptTolerance(nlntolsol,nlnres,adaptolbetter);
  }

  //-------solve for residual displacements to correct incremental displacements
  increment_->PutScalar(0.0);

  // always refactor the matrix for a new solver call --- we assume that
  // it has changed since the last call
  bool refactor=true;
  // never reset solver from time integration level unless it is the first step
  // the preconditioner does the job on its own according to the AZreuse
  // parameter
  bool reset=false;
  if(step_==1 && itenum_ == 1)
  {
    reset=true;
  }

  solver_->Solve(sysmat_->EpetraOperator(),
                increment_               ,
                residual_                ,
                refactor                 ,
                reset                    ,
                w_                       ,
                c_                       ,
                project_                 );

  solver_->ResetTolerance();

  return;
} // FluidGenAlphaIntegration::GenAlphaCalcIncrement

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | update the current acceleration, velocity and pressure    gammi 06/07|
 -----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidGenAlphaIntegration::GenAlphaNonlinearUpdate()
{

  // -------------------------------------------------------------------
  // split between accelerations and pressure increments
  Teuchos::RCP<Epetra_Vector> accinc = velpressplitter_.ExtractOtherVector(increment_);
  Teuchos::RCP<Epetra_Vector> preinc = velpressplitter_.ExtractCondVector (increment_);

  // ------------------------------------------------------
  // update acceleration
  //
  //        n+1         n+1
  //     acc      =  acc    + dacc
  //        (i+1)       (i)
  //
  velpressplitter_.AddOtherVector(accinc,accnp_);

  // ------------------------------------------------------
  // use updated acceleration to update velocity. Since
  //
  //    n+1         n            n                 +-   n+1       n -+
  // vel      =  vel   + dt * acc   + gamma * dt * | acc     - acc   | =
  //    (i+1)                                      +-   (i+1)       -+
  //
  //                n            n                 +-   n+1       n -+
  //          =  vel   + dt * acc   + gamma * dt * | acc     - acc   | +
  //                                               +-   (i)         -+
  //
  //                                      n+1
  //             + gamma * dt * dacc = vel     +  gamma * dt * dacc =
  //                                      (i)
  //               n+1
  //          = vel     +   dvel
  //               (i)           ,
  //
  // the velocity could be updated incremently.
  //
  // Incremental updates seem to be dangerous in this place. Round-off
  // errors could evolve on velocity and pressure seperately through the
  // nonlinear process. That means that accelerations and velocities at
  // the current iteration level might get decoupled due to numerical
  // errors.
  //
  // For this reason we use a different update formula here, making sure
  // that accelerations and velocities at the current iteration level are
  // alway syncronised:
  //
  //    n+1         n    +-         -+           n                   n+1
  // vel      =  vel   + | 1 - gamma | * dt * acc  + gamma * dt * acc
  //    (i+1)            +-         -+                               (i+1)
  //
  Teuchos::RCP<Epetra_Vector> vel    = velpressplitter_.ExtractOtherVector(veln_ );
  Teuchos::RCP<Epetra_Vector> accold = velpressplitter_.ExtractOtherVector(accn_ );
  Teuchos::RCP<Epetra_Vector> accnew = velpressplitter_.ExtractOtherVector(accnp_);

  vel->Update((1.0-gamma_)*dt_,*accold,gamma_*dt_,*accnew,1.0);
  velpressplitter_.InsertOtherVector(vel,velnp_);

  // ------------------------------------------------------
  // update pressure
  //
  //         n+1          n+1
  //     pres      =  pres    + dpres
  //         (i+1)        (i)
  //

  // rescaled pressure to preserve symmetry of pressure
  // and continuity part in matrix
  preinc->Scale(gamma_*dt_);
  velpressplitter_.AddCondVector(preinc,velnp_);
  preinc->Scale(1.0/gamma_*dt_);

  return;
} // FluidGenAlphaIntegration::GenAlphaNonlinearUpdate



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | check for convergence of nonlinear iteration              gammi 06/07|
 -----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
bool FLD::FluidGenAlphaIntegration::GenAlphaNonlinearConvergenceCheck(double& badestnlnnorm)
{
  bool stopnonliniter = false;

  // extract velocity and pressure increments from increment vector
  Teuchos::RCP<Epetra_Vector> onlyvel = velpressplitter_.ExtractOtherVector(increment_);
  Teuchos::RCP<Epetra_Vector> onlypre = velpressplitter_.ExtractCondVector (increment_);

  // calculate L2_Norm of increments
  double L2incaccnorm;
  onlyvel->Norm2(&L2incaccnorm);
  onlypre->Norm2(&L2incprenorm_);

  L2incvelnorm_= L2incaccnorm*gamma_*dt_;

  // rescaling of pressure to keep symmetry of matrix
  L2incprenorm_*=gamma_*dt_;

  // extract velocity and pressure solutions from solution vector
  onlyvel = velpressplitter_.ExtractOtherVector(velnp_);
  onlypre = velpressplitter_.ExtractCondVector (velnp_);

  // calculate L2_Norm of solution
  double L2velnorm;
  double L2prenorm;
  onlyvel->Norm2(&L2velnorm);
  onlypre->Norm2(&L2prenorm);


  if (L2velnorm<EPS5)
  {
    L2velnorm = 1.0;
  }
  if (L2prenorm<EPS5)
  {
    L2prenorm = 1.0;
  }

  // compute relative increment
  L2incvelnorm_ /= L2velnorm;
  L2incprenorm_ /= L2prenorm;

  // extract velocity and pressure residuals from rhs vector
  onlyvel = velpressplitter_.ExtractOtherVector(residual_);
  onlypre = velpressplitter_.ExtractCondVector (residual_);

  onlypre->Norm2(&L2preresnorm_);
  onlyvel->Norm2(&L2velresnorm_);

  badestnlnnorm = max(L2preresnorm_,L2velresnorm_);
  badestnlnnorm = max(badestnlnnorm,L2incvelnorm_);
  badestnlnnorm = max(badestnlnnorm,L2incprenorm_);

  // out to screen
  this->GenAlphaEchoToScreen("print residual and increment values");
  this->GenAlphaEchoToScreen("print timing of solver and element calls");

  // this is the convergence check
  if((L2incvelnorm_ <= ittol_ && L2incprenorm_ <= ittol_)
     &&
     (L2velresnorm_ <= ittol_ && L2preresnorm_ <= ittol_))
  {
    stopnonliniter=true;
    this->GenAlphaEchoToScreen("print nonlin iter converged");
  }

  // currently default: only one iteration before sampling
  if (special_flow_ == "channel_flow_of_height_2" && step_<samstart_ )
  {
    stopnonliniter=true;
    this->GenAlphaEchoToScreen("print warning, only one iteration before sampling");
  }

  // warn if itemax is reached without convergence, but proceed to
  // next timestep...
  if (itenum_ == itemax_ && stopnonliniter!=true)
  {
    stopnonliniter=true;
    this->GenAlphaEchoToScreen("print warning, nonlin iter not converged");
  }

  return stopnonliniter;
} // FluidGenAlphaIntegration::GenAlphaNonlinearConvergenceCheck


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Read restart information                                  gammi 06/07|
 -----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidGenAlphaIntegration::ReadRestart(int step)
{

  IO::DiscretizationReader reader(discret_,step);
  time_ = reader.ReadDouble ("time");
  step_ = reader.ReadInt    ("step");

  if(myrank_==0)
  {
    printf("Reading restart info to restart at step %d, time %5f\n\n",step_,time_);
  }

  reader.ReadVector(velnp_,"velnp");
  reader.ReadVector(veln_ ,"veln" );
  reader.ReadVector(accn_ ,"accn" );

  if (alefluid_)
  {
    reader.ReadVector(dispnp_  ,"dispnp"  );
    reader.ReadVector(dispn_   ,"dispn"   );
    reader.ReadVector(dispnm_  ,"dispnm"  );
    reader.ReadVector(gridveln_,"gridveln");
  }

  // read the previously written elements including the history data

  reader.ReadMesh(step_);

  // read previous averages
//  statisticsmanager_->Restart(reader,step_);

  // since ReadMesh can change the overall dof numbering due
  // to a FillComplete() call performed internally, the underlying maps
  // of state vectors can become illegal - especially in case of
  // multiphysics problems & periodic boundary conditions!
  // It is better to check consistency here:
    if (not (discret_->DofRowMap())->SameAs(velnp_->Map()))
      dserror("Global dof numbering in maps does not match");
    if (not (discret_->DofRowMap())->SameAs(veln_->Map()))
      dserror("Global dof numbering in maps does not match");
    if (not (discret_->DofRowMap())->SameAs(accn_->Map()))
      dserror("Global dof numbering in maps does not match");

  return;
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | prepare AVM3-based scale separation                         vg 10/08 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidGenAlphaIntegration::AVM3Preparation()
{
  // zero matrix
  sysmat_->Zero();

  // add Neumann loads
  residual_->Update(1.0,*neumann_loads_,0.0);

  // create the parameters for the discretization
  ParameterList eleparams;

  // action for elements
  eleparams.set<int>("action",FLD::calc_fluid_genalpha_sysmat_and_residual);
  eleparams.set("assemble matrix 1",true);
  eleparams.set("assemble matrix 2",false);
  eleparams.set("assemble vector 1",true);
  eleparams.set("assemble vector 2",false);
  eleparams.set("assemble vector 3",false);
  eleparams.set("compute element matrix",true);

  // other parameters that might be needed by the elements
  {
    ParameterList& timelist = eleparams.sublist("time integration parameters");

    timelist.set("alpha_M",alphaM_);
    timelist.set("alpha_F",alphaF_);
    timelist.set("gamma"  ,gamma_ );
    timelist.set("time"   ,time_  );
    timelist.set("dt"     ,dt_    );
  }

  // parameters for nonlinear treatment (linearisation) and low-Mach-number solver
  eleparams.set<int>("Linearisation",newton_);

  // parameters for stabilisation
  {
    eleparams.sublist("STABILIZATION")    = params_->sublist("STABILIZATION");
  }

  // parameters for a turbulence model
  {
    eleparams.sublist("TURBULENCE MODEL") = params_->sublist("TURBULENCE MODEL");
  }

  // (fine-scale) subgrid viscosity flag
  eleparams.set("fs subgrid viscosity","No");

  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("u and p (n+1      ,trial)",velnp_ );
  discret_->SetState("u and p (n+alpha_F,trial)",velaf_ );
  discret_->SetState("acc     (n+alpha_M,trial)",accam_ );

  // element evaluation for getting system matrix
  // -> we merely need matrix "structure" below, not the actual contents
  discret_->Evaluate(eleparams,sysmat_,residual_);
  discret_->ClearState();

  // complete system matrix
  sysmat_->Complete();

  // apply DBC to system matrix
  LINALG::ApplyDirichlettoSystem(sysmat_,increment_,residual_,zeros_,*(dbcmaps_->CondMap()));

  // get scale-separation matrix
  {
    // this is important to have!!!
    MLAPI::Init();

    // extract the ML parameters
    ParameterList&  mlparams = solver_->Params().sublist("ML Parameters");;

    // get toggle vector for Dirchlet boundary conditions
    const Epetra_Vector& dbct = *Dirichlet();

    // get nullspace parameters
    double* nullspace = mlparams.get("null space: vectors",(double*)NULL);
    if (!nullspace) dserror("No nullspace supplied in parameter list");
    int nsdim = mlparams.get("null space: dimension",1);

    // modify nullspace to ensure that DBC are fully taken into account
    if (nullspace)
    {
      const int length = SystemMatrix()->OperatorRangeMap().NumMyElements();
      for (int i=0; i<nsdim; ++i)
        for (int j=0; j<length; ++j)
          if (dbct[j]!=0.0) nullspace[i*length+j] = 0.0;
    }

    // get plain aggregation Ptent
    RCP<Epetra_CrsMatrix> crsPtent;
    GetPtent(*SystemMatrix()->EpetraMatrix(),mlparams,nullspace,crsPtent);
    LINALG::SparseMatrix Ptent(crsPtent);

    // compute scale-separation matrix: S = I - Ptent*Ptent^T
    Sep_ = LINALG::Multiply(Ptent,false,Ptent,true);
    Sep_->Scale(-1.0);
    RCP<Epetra_Vector> tmp = LINALG::CreateVector(Sep_->RowMap(),false);
    tmp->PutScalar(1.0);
    RCP<Epetra_Vector> diag = LINALG::CreateVector(Sep_->RowMap(),false);
    Sep_->ExtractDiagonalCopy(*diag);
    diag->Update(1.0,*tmp,1.0);
    Sep_->ReplaceDiagonalValues(*diag);

    //complete scale-separation matrix and check maps
    Sep_->Complete(Sep_->DomainMap(),Sep_->RangeMap());
    if (!Sep_->RowMap().SameAs(SystemMatrix()->RowMap())) dserror("rowmap not equal");
    if (!Sep_->RangeMap().SameAs(SystemMatrix()->RangeMap())) dserror("rangemap not equal");
    if (!Sep_->DomainMap().SameAs(SystemMatrix()->DomainMap())) dserror("domainmap not equal");
  }

  return;
}// FluidGenAlphaIntegration::AVM3Preparation


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | AVM3-based scale separation                                 vg 10/08 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidGenAlphaIntegration::AVM3Separation()
{
  // get fine-scale part of velocity
  Sep_->Multiply(false,*velaf_,*fsvelaf_);

  // set fine-scale vector
  discret_->SetState("fsvelaf",fsvelaf_);

  return;
}// FluidGenAlphaIntegration::AVM3Separation


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  set initial flow field for test cases                    gammi 06/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidGenAlphaIntegration::SetInitialFlowField(
  const INPAR::FLUID::InitialField initfield,
  const int startfuncno
 )
{

  //----------------------------------------------------------------------
  // Initialfield from function
  //----------------------------------------------------------------------
  if(initfield == INPAR::FLUID::initfield_field_by_function or
     initfield == INPAR::FLUID::initfield_disturbed_field_from_function)
  {
    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = discret_->lRowNode(lnodeid);
      // the set of degrees of freedom associated with the node
      vector<int> nodedofset = discret_->Dof(lnode);

      for(int index=0;index<numdim_+1;++index)
      {
        int gid = nodedofset[index];

        double initialval=DRT::Problem::Instance()->Funct(startfuncno-1).Evaluate(index,lnode->X(),time_,NULL);

        velnp_->ReplaceGlobalValues(1,&initialval,&gid);
        veln_ ->ReplaceGlobalValues(1,&initialval,&gid);
      }
    }

    // for NURBS discretizations we have to solve a least squares problem,
    // with high accuracy! (do nothing for Lagrangian polynomials)
    DRT::NURBS::apply_nurbs_initial_condition(
        *discret_  ,
        DRT::Problem::Instance()->ErrorFile()->Handle(),
        DRT::Problem::Instance()->UMFPACKSolverParams(),
        startfuncno,
        velnp_     );

    // initialize veln_ as well. That's what we actually want to do here!
      veln_->Update(1.0,*velnp_ ,0.0);

    //----------------------------------------------------------------------
    // random perturbations for field
    //----------------------------------------------------------------------
    if (initfield == INPAR::FLUID::initfield_disturbed_field_from_function)
    {
      //
      {
        const Epetra_Map* dofrowmap = discret_->DofRowMap();

        int err =0;

        // random noise is perc percent of the initial profile

        double perc = params_->sublist("TURBULENCE MODEL").get<double>("CHAN_AMPL_INIT_DIST",0.1);

        // out to screen
        if (myrank_==0)
        {
          cout << "Disturbed initial profile:   max. " << perc*100 << "% random perturbation\n";
          cout << "\n\n";
        }

        // loop all nodes on the processor
        for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();++lnodeid)
        {
          // get the processor local node
          DRT::Node*  lnode      = discret_->lRowNode(lnodeid);
          // the set of degrees of freedom associated with the node
          vector<int> nodedofset = discret_->Dof(lnode);

          // check whether we have a pbc condition on this node
          vector<DRT::Condition*> mypbc;

          lnode->GetCondition("SurfacePeriodic",mypbc);

          bool is_slave=false;

          // check whether a periodic boundary condition is active on this node
          if (mypbc.size()>0)
          {
            // yes, we have one
            for (unsigned numcond=0;numcond<mypbc.size();++numcond)
            {

              const string* mymasterslavetoggle
                = mypbc[numcond]->Get<string>("Is slave periodic boundary condition");

              if(!(*mymasterslavetoggle=="Master"))
              {
                // the node is a slave --- so don't do anything
                is_slave=true;
              }
            }
          }

          if(is_slave)
          {
            continue;
          }

          // the noise is proportional to the maximum component of the
          // undisturbed initial field in this point
          double initialval=0;

          // direction with max. profile
          int flowdirection = 0;

          for(int index=0;index<numdim_;++index)
          {
            int gid = nodedofset[index];
            int lid = dofrowmap->LID(gid);

            double thisval=(*velnp_)[lid];
            if (initialval*initialval < thisval*thisval)
            {
              initialval=thisval;

              // remember the direction of maximum flow
              flowdirection = index;
            }
          }

          // add random noise on initial function field
          for(int index=0;index<numdim_;++index)
          {
            int gid = nodedofset[index];

            double randomnumber = DRT::Problem::Instance()->Random()->Uni();

            double noise = initialval * randomnumber * perc;

            // full noise only in main flow direction
              // one third noise orthogonal to flow direction
            if (index != flowdirection)
            {
              noise *= 1./3.;
            }

            err += velnp_->SumIntoGlobalValues(1,&noise,&gid);
            err += veln_ ->SumIntoGlobalValues(1,&noise,&gid);

          }

          if(err!=0)
          {
            dserror("dof not on proc");
          }
        }
      }
    }
  }
  //----------------------------------------------------------------------
  // Initialfield for Beltrami flow
  //----------------------------------------------------------------------
  else if(initfield == INPAR::FLUID::initfield_beltrami_flow)
  {
    int gid;
    int lid;

    const Epetra_Map* dofrowmap = discret_->DofRowMap();


    int err =0;

    int numdim  = params_->get<int>("number of velocity degrees of freedom");
    int npredof = numdim;

    double         p;
    vector<double> u  (numdim);
    vector<double> xyz(numdim);


    if(numdim!=3)
    {
      dserror("Beltrami flow is three dimensional flow!");
    }

    // set constants for analytical solution
    double a      = M_PI/4.0;
    double d      = M_PI/2.0;

    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = discret_->lRowNode(lnodeid);

      // the set of degrees of freedom associated with the node
      vector<int> nodedofset = discret_->Dof(lnode);

      // set node coordinates
      for(int dim=0;dim<numdim;dim++)
      {
        xyz[dim]=lnode->X()[dim];
      }

      // compute initial pressure
      p = -a*a/2.0 *
        ( exp(2.0*a*xyz[0])
          + exp(2.0*a*xyz[1])
          + exp(2.0*a*xyz[2])
          + 2.0 * sin(a*xyz[0] + d*xyz[1]) * cos(a*xyz[2] + d*xyz[0]) * exp(a*(xyz[1]+xyz[2]))
          + 2.0 * sin(a*xyz[1] + d*xyz[2]) * cos(a*xyz[0] + d*xyz[1]) * exp(a*(xyz[2]+xyz[0]))
          + 2.0 * sin(a*xyz[2] + d*xyz[0]) * cos(a*xyz[1] + d*xyz[2]) * exp(a*(xyz[0]+xyz[1]))
          );

      // compute initial velocities
      u[0] = -a * ( exp(a*xyz[0]) * sin(a*xyz[1] + d*xyz[2]) +
                    exp(a*xyz[2]) * cos(a*xyz[0] + d*xyz[1]) );
      u[1] = -a * ( exp(a*xyz[1]) * sin(a*xyz[2] + d*xyz[0]) +
                    exp(a*xyz[0]) * cos(a*xyz[1] + d*xyz[2]) );
      u[2] = -a * ( exp(a*xyz[2]) * sin(a*xyz[0] + d*xyz[1]) +
                    exp(a*xyz[1]) * cos(a*xyz[2] + d*xyz[0]) );
      // initial velocities
      for(int nveldof=0;nveldof<numdim;nveldof++)
      {
        gid = nodedofset[nveldof];
        lid = dofrowmap->LID(gid);
        err += velnp_->ReplaceMyValues(1,&(u[nveldof]),&lid);
        err += veln_ ->ReplaceMyValues(1,&(u[nveldof]),&lid);
     }

      // initial pressure
      gid = nodedofset[npredof];
      lid = dofrowmap->LID(gid);
      err += velnp_->ReplaceMyValues(1,&p,&lid);
      err += veln_ ->ReplaceMyValues(1,&p,&lid);

    } // end loop nodes lnodeid
    if(err!=0)
    {
      dserror("dof not on proc");
    }

  }
  else
  {
    dserror("no other initial fields than zero, function and beltrami are available up to now");
  }

  // for safety purposes --- make the initial flow field compatible
  // to the initial Dirichlet boundary condition
  {
    ParameterList eleparams;

    // total time required for Dirichlet conditions
    eleparams.set("total time",time_);
    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("velnp",velnp_);
    // predicted dirichlet values
    // velnp then also holds prescribed new dirichlet values
    discret_->EvaluateDirichlet(eleparams,velnp_,null,null,null,dbcmaps_);
    discret_->ClearState();

    // set vector values needed by elements
    discret_->SetState("velnp",veln_);
    // predicted dirichlet values
    // velnp then also holds prescribed new dirichlet values
    discret_->EvaluateDirichlet(eleparams,veln_,null,null,null,dbcmaps_);
    discret_->ClearState();
  }

  return;
} // end FluidGenAlphaIntegration::SetInitialFlowField


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | evaluate error for test cases with analytical solutions   gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidGenAlphaIntegration::EvaluateErrorComparedToAnalyticalSol()
{
  INPAR::FLUID::CalcError calcerr = DRT::INPUT::get<INPAR::FLUID::CalcError>(*params_,"calculate error");

  //------------------------------------------------------- beltrami flow
  switch (calcerr)
  {
  case INPAR::FLUID::initfield_zero_field:
  case INPAR::FLUID::initfield_field_by_function:
  case INPAR::FLUID::initfield_disturbed_field_from_function:
  case INPAR::FLUID::initfield_flame_vortex_interaction:
    // do nothing --- no analytical solution available
    break;
  case INPAR::FLUID::initfield_beltrami_flow:
  {
    // create the parameters for the discretization
    ParameterList eleparams;

    eleparams.set<double>("L2 integrated velocity error",0.0);
    eleparams.set<double>("L2 integrated pressure error",0.0);

    // action for elements
    eleparams.set("action","calc_fluid_beltrami_error");
    // actual time for elements
    eleparams.set("total time",time_);

    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("u and p at time n+1 (converged)",velnp_);

    // call loop over elements
    discret_->Evaluate(eleparams,null,null,null,null,null);
    discret_->ClearState();


    double locvelerr = eleparams.get<double>("L2 integrated velocity error");
    double locpreerr = eleparams.get<double>("L2 integrated pressure error");

    double velerr = 0;
    double preerr = 0;

    discret_->Comm().SumAll(&locvelerr,&velerr,1);
    discret_->Comm().SumAll(&locpreerr,&preerr,1);

    // for the L2 norm, we need the square root
    velerr = sqrt(velerr);
    preerr = sqrt(preerr);


    if (myrank_ == 0)
    {
      printf("\n  L2_err for beltrami flow:  velocity %15.8e  pressure %15.8e\n\n",
             velerr,preerr);
    }
  }
  break;
  default:
    dserror("Cannot calculate error. Unknown type of analytical test problem");
  break;
  }

  return;
} // end EvaluateErrorComparedToAnalyticalSol


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | filter quantities for dynamic Smagorinsky model. Compute averaged    |
 | values for LijMij and MijMij.                             gammi 12/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidGenAlphaIntegration::ApplyFilterForDynamicComputationOfCs()
{
  TEUCHOS_FUNC_TIME_MONITOR("FLD::FluidGenAlphaIntegration::ApplyFilterForDynamicComputationOfCs");

  // time measurement
  double tcpu=Teuchos::Time::wallTime();

  // perform filtering and computation of Cs
  {
    const Teuchos::RCP<const Epetra_Vector> dirichtoggle = Dirichlet();
    DynSmag_->ApplyFilterForDynamicComputationOfCs(velaf_,scaaf_,0.0,dirichtoggle);
  }

  dtfilter_=Teuchos::Time::wallTime()-tcpu;

  return;
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Print info to shell                                      gammi 01/08|
 -----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidGenAlphaIntegration::GenAlphaEchoToScreen(
  const string& what_to_print
  )
{
  if (myrank_==0)
  {
    if(what_to_print == "print start-up info")
    {
      //--------------------------------------------------------------------
      /* output of initial time integration related data */
      {
        // alpha_F
        cout << "Generalized Alpha parameter: alpha_F = ";
        cout << params_->get<double>("alpha_F");
        cout << "\n";

        // alpha_M
        cout << "                             alpha_M = ";
        cout << params_->get<double>("alpha_M");
        cout << "\n";

        // gamma
        cout << "                             gamma   = ";
        cout << gamma_;
        cout << "\n";

        if(abs(gamma_  - (0.5 + alphaM_ - alphaF_)>1e-6))
        {
          cout << "Definition of gamma NOT second order accurate, should be " << (0.5 + alphaM_ - alphaF_);
        }
        cout << "\n";

        // linearisation
	cout << "Linearisation (1=fixed_point; 2=Newton; 3=minimal): ";
	cout << DRT::INPUT::get<INPAR::FLUID::LinearisationAction>(*params_, "Linearisation");
        cout << endl;

        // predictor
	cout << "                             ";
        cout << "using " << predictor_ << " as an initial guess" << endl << endl;
      }

      //--------------------------------------------------------------------
      /* output of stabilisation details */
      {
        ParameterList *  stabparams=&(params_->sublist("STABILIZATION"));

        // general

        string stabtype =stabparams->get<string>("STABTYPE");

        cout << "Stabilization type         : ";
        cout << stabtype;
        cout << endl;

        if (stabtype == "residual_based"||stabtype == "inconsistent")
        {
          cout << "                             ";
          cout << "tau according to ";
          cout << stabparams->get<string>("DEFINITION_TAU");
          cout << endl;
        }

        // time-dependent subgrid scales
        cout << "                             ";
        cout << stabparams->get<string>("TDS")<< endl;
        cout << endl;

        // transient term with consistency check for quasistatic subgrid scales
        if(stabparams->get<string>("TDS") == "quasistatic")
        {
          if(stabparams->get<string>("TRANSIENT")=="yes_transient")
          {
            dserror("The quasistatic version of the residual-based stabilization currently does not support the incorporation of the transient term.");
          }
        }
        cout << "                             ";
        cout << "TRANSIENT       = ";
        cout << stabparams->get<string>("TRANSIENT");
        cout << endl;

        // supg stabilization?
        cout <<  "                             ";
        cout << "SUPG            = ";
        cout << stabparams->get<string>("SUPG")           ;
        cout << endl;

        // pspg stabilization?
        cout <<  "                             ";
        cout << "PSPG            = ";
        cout << stabparams->get<string>("PSPG")           ;
        cout << endl;

        // viscous stabilization?
        cout <<  "                             ";
        cout << "VSTAB           = ";
        cout << stabparams->get<string>("VSTAB")          ;
        cout << endl;

        // least-squares continuity stabilization?
        cout <<  "                             ";
        cout << "CSTAB           = ";
        cout << stabparams->get<string>("CSTAB")          ;
        cout << endl;

        // resvmm turbulence modeling, cross-stress part?
        cout <<  "                             ";
        cout << "CROSS-STRESS    = ";
        cout << stabparams->get<string>("CROSS-STRESS")   ;
        cout << endl;

        // resvmm turbulence modeling, Reynolds-stress part?
        cout <<  "                             ";
        cout << "REYNOLDS-STRESS = ";
        cout << stabparams->get<string>("REYNOLDS-STRESS");
        cout << endl;
        cout << endl;

        // do we use a conservative approach?
        cout << "                             ";
        cout << "Choosing cross-stress stabilisation based on ";
        cout << params_->get<string>("form of convective term");
        cout << " equation";
        cout << endl;
        cout << endl;

      }

      //--------------------------------------------------------------------
      /* output of turbulence model if any */
      {
        ParameterList *  modelparams =&(params_->sublist("TURBULENCE MODEL"));

        if (modelparams->get<string>("TURBULENCE_APPROACH", "none")
            !=
            "none")
        {
          // a canonical flow with homogeneous directions would allow a
          // spatial averaging of data
          string special_flow_
            =
            modelparams->get<string>("CANONICAL_FLOW","no");

          // get homogeneous directions
          string homdir = modelparams->get<string>("HOMDIR","not_specified");

          // Underresolved DNS, traditional LES (Smagorinsky type), RANS?
          // including consistecy checks
          cout << "Turbulence approach        : ";
          cout << modelparams->get<string>("TURBULENCE_APPROACH");
          cout << endl << endl;

          if(modelparams->get<string>("TURBULENCE_APPROACH")
             ==
             "RANS")
          {
            dserror("RANS approaches not implemented yet\n");
          }
          else if(modelparams->get<string>("TURBULENCE_APPROACH")
                  ==
                  "CLASSICAL_LES")
          {
            string physmodel = modelparams->get<string>("PHYSICAL_MODEL");

            cout << "                             ";
            cout << physmodel;
            cout << endl;

            if (physmodel == "Smagorinsky")
            {
              cout << "                             " ;
              cout << "with Smagorinsky constant Cs= ";
              cout << params_->sublist("SUBGRID VISCOSITY").get<double>("C_SMAGORINSKY") ;
            }
            else if(physmodel == "Smagorinsky_with_van_Driest_damping")
            {
              if (special_flow_ != "channel_flow_of_height_2"
                  ||
                  homdir != "xz")
              {
                dserror("The van Driest damping is only implemented for a channel flow with wall \nnormal direction y");
              }

              cout << "                             "          ;
              cout << "- Smagorinsky constant:   Cs   = "      ;
              cout << params_->sublist("SUBGRID VISCOSITY").get<double>("C_SMAGORINSKY");
              cout << endl;

              cout << "                             "          ;
              cout << "- viscous length      :   l_tau= "      ;
              cout << params_->sublist("SUBGRID VISCOSITY").get<double>("CHANNEL_L_TAU");
              cout << endl;
            }
            else if(physmodel == "Dynamic_Smagorinsky")
            {
              if (special_flow_ != "channel_flow_of_height_2"
                  ||
                  homdir != "xz")
              {
                cout << "                             ";
                cout << "+ clipping of negative values\n";
              }
              else
              {
                cout << "                             ";
                cout << "+ clipping of negative values\n";
                cout << "                             ";
                cout << "+ in plane averaging of Cs_delta_sq\n";
              }
            }
            cout << endl;
          }
        }

        //--------------------------------------------------------------------
        /* output of fine-scale subgrid-vicosity approach if any */
        if (fssgv_ != "No")
        {
          cout << "Fine-scale subgrid-viscosity approach based on AVM3: ";
          cout << endl << endl;
          cout << params_->get<string>("fs subgrid viscosity");
          cout << " with Smagorinsky constant Cs= ";
          cout << params_->sublist("SUBGRID VISCOSITY").get<double>("C_SMAGORINSKY") ;
          cout << endl << endl;
        }

      }

    }
    else if (what_to_print == "print time algorithm info")
    {
      //--------------------------------------------------------------------
      /* output of current time integration data */

      printf("TIME: %11.4E/%11.4E  DT = %11.4E     GenAlpha     STEP = %4d/%4d \n",
             time_,endtime_,dt_,step_,endstep_);
    }
    else if (what_to_print == "print first nonlinear iteration info")
    {
      printf("+------------+");
      printf("-------------------+");
      printf("--------------+");
      printf("--------------+");
      printf("--------------+");
      printf("--------------+\n");
      printf("|- step/max -|");
      printf("- tol      [norm] -|");
      printf("- vel-error --|");
      printf("- pre-error --|");
      printf("- vres-norm --|");
      printf("- pres-norm --|\n");

    }
    else if (what_to_print == "print residual norms before first iteration")
    {
      //--------------------------------------------------------------------
      /* output of residuals before first solution */

      printf(      "|     ---    |");
      printf("           [L_2 ]  |");
      printf(     "      ---     |");
      printf(     "      ---     |");
      printf(         " %10.3E   |",L2velresnorm_);
      printf(         " %10.3E   |",L2preresnorm_);

      fflush(stdout);
    }
    else if (what_to_print == "print timing of element calls")
    {
      //--------------------------------------------------------------------
      /* output of timing for first element call (without solver) */

      printf("                 (te=%10.3E)",dtele_);
      // additional output for dynamic Smagorinsky model --- the time spent on
      // filtering
      if (params_->sublist("TURBULENCE MODEL").get<string>("TURBULENCE_APPROACH",
                                                          "DNS_OR_RESVMM_LES")
          ==
          "CLASSICAL_LES")
      {
        if(params_->sublist("TURBULENCE MODEL").get<string>("PHYSICAL_MODEL",
                                                           "no_model")
           ==
           "Dynamic_Smagorinsky"
          )
        {
          printf("(tf=%10.3E)",dtfilter_);
        }
      }
      cout << endl;

      // reset time increment after output
      dtele_    = 0;
      dtfilter_ = 0;
    }
    else if (what_to_print == "print timing of solver and element calls")
    {
      //--------------------------------------------------------------------
      /* output of timing for further element calls (including solver) */

      printf("  (ts=%10.3E)(te=%10.3E)",dtsolve_,dtele_);
      // additional output for dynamic Smagorinsky model --- the time spent on
      // filtering
      if (params_->sublist("TURBULENCE MODEL").get<string>("TURBULENCE_APPROACH",
                                                          "DNS_OR_RESVMM_LES")
          ==
          "CLASSICAL_LES")
      {
        if(params_->sublist("TURBULENCE MODEL").get<string>("PHYSICAL_MODEL",
                                                           "no_model")
           ==
           "Dynamic_Smagorinsky"
          )
        {
          printf("(tf=%10.3E)",dtfilter_);
        }
      }
      cout << endl;

      // reset time increments after output
      dtsolve_  = 0;
      dtele_    = 0;
      dtfilter_ = 0;
    }
    else if(what_to_print == "print residual and increment values")
    {
      if ((itenum_ == itemax_) && skiplastelecall_ )
      {
        printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   |",
               itenum_,itemax_,ittol_,L2incvelnorm_,L2incprenorm_);
        printf("     ---      |      ---     |");
      }
      else
      {
        printf("|  %3d/%3d   |" ,itenum_,itemax_);
        printf(" %10.3E[L_2 ]  |",ittol_         );
        printf(" %10.3E   |"     ,L2incvelnorm_  );
        printf(" %10.3E   |"     ,L2incprenorm_  );
        printf(" %10.3E   |"     ,L2velresnorm_  );
        printf(" %10.3E   |"     ,L2preresnorm_  );
      }
    }
    else if (what_to_print == "print warning, nonlin iter not converged")
    {
      //--------------------------------------------------------------------
      /* output of warning if no convergence was achieved */

      printf("+--------------------------------------------------------------------------------------------+\n");
      printf("| >>>>>> not converged in itemax steps! matrix of last step not recomputed (invalid)         |\n");
      printf("+--------------------------------------------------------------------------------------------+\n");
    }
    else if (what_to_print == "print warning, only one iteration before sampling")
    {
      //--------------------------------------------------------------------
      /* output of warning if no convergence was achieved */

      printf("+--------------------------------------------------------------------------------------------+\n");
      printf("| >>>>>> only one iteration before sampling! matrix of last step not recomputed (invalid)    |\n");
      printf("+--------------------------------------------------------------------------------------------+\n");
    }
    else if(what_to_print == "print nonlin iter converged")
    {
      //--------------------------------------------------------------------
      /* close box in case of convergence */

      printf("+------------+-------------------+--------------+--------------+--------------+--------------+ \n");
    }
    else
    {
      dserror("Don't know what to print\n");
    }

    fflush(stdout);

  } // if (myrank_==0)


  return;
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |                                                           gammi 02/08|
 -----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidGenAlphaIntegration::UpdateGridv()
{
  /*
       2nd order finite difference approximation of the new grid velocity

                    (equivalent to what is done for BDF2)

                           /       \
        n+1        3      |  n+1  n |    1       n
     u_G     = -------- * | d   -d  | - --- * u_G
                2 * dt    |         |    2
                           \       /
  */
  gridvelaf_->Update(-0.5,*gridveln_,0.0);
  gridvelaf_->Update( 3.0/(2.0*dt_),*dispnp_,
                     -3.0/(2.0*dt_),*dispn_,1.0);

  /*

    now do a linear interpolation to the intermediate time level

            n+af      n+1               n    /          \
         u_G     = u_G    * alphaF + u_G  * | 1 - alphaF |
                                             \          /
  */
  gridvelaf_->Update((1-alphaF_),*gridveln_,alphaF_);

#if 0
  /*
                             /       \     /            \
        n+af     alphaM     |  n+1  n |   |       alphaM |   .n   . n+alphaM
     u_G     = ---------- * | d   -d  | + | 1.0 - ------ | * d  = d
               gamma * dt   |         |   |        gamma |
                             \       /     \            /

  */

  gridvelaf_->Update(1.0-alphaM_/gamma_,*gridveln_,0.0);

  gridvelaf_->Update( alphaM_/(gamma_*dt_),*dispnp_,
                     -alphaM_/(gamma_*dt_),*dispn_,1.0);



  /*
                shouldn't it be something like

                             /       \     /            \
        n+af     alphaF     |  n+1  n |   |       alphaF |   .n
     u_G     = ---------- * | d   -d  | + | 1.0 - ------ | * d
               gamma * dt   |         |   |        gamma |
                             \       /     \            /


      ??????????????????????????????????????????????????????

  gridvelaf_->Update(1.0-alphaF_/gamma_,*gridveln_,0.0);

  gridvelaf_->Update( alphaF_/(gamma_*dt_),*dispnp_,
                     -alphaF_/(gamma_*dt_),*dispn_,1.0);
  */
#endif
  return;
}



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | LiftDrag                                                  chfoe 11/07|
 |                                    copied and modified by gammi 02/08|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*!
\brief calculate lift&drag forces and angular momenta

Lift and drag forces are based upon the right hand side (force) entities
of the corresponding nodes. The contribution of the end node of a line
is entirely added to a present L&D force.

Notice: Angular moments obtained from lift&drag forces currently refere
        to the initial configuration, i.e. are built with the coordinates
        X of a particular node irrespective of its current position.
*/
void FLD::FluidGenAlphaIntegration::LiftDrag() const
{
  // in this map, the results of the lift drag calculation are stored
  RCP<map<int,vector<double> > > liftdragvals;

  FLD::UTILS::LiftDrag(*discret_,*force_,*params_,liftdragvals);


}//FluidGenAlphaIntegration::LiftDrag


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  calculate traction vector at (Dirichlet) boundary (public) gjb 05/08|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
Teuchos::RCP<Epetra_Vector> FLD::FluidGenAlphaIntegration::CalcStresses()
{
  string condstring("FluidStressCalc");
  Teuchos::RCP<Epetra_Vector> integratedshapefunc = IntegrateInterfaceShape(condstring);

  // compute traction values at specified nodes; otherwise do not touch the zero values
  for (int i=0;i<integratedshapefunc->MyLength();i++)
  {
    if ((*integratedshapefunc)[i] != 0.0)
    {
      // overwrite integratedshapefunc values with the calculated traction coefficients,
      // which are reconstructed out of the nodal forces (trueresidual_) using the
      // same shape functions on the boundary as for velocity and pressure.
      (*integratedshapefunc)[i] = (*force_)[i]/(*integratedshapefunc)[i];
    }
  }

  return integratedshapefunc;

} // FluidGenAlphaIntegration::CalcStresses()


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FLD::FluidGenAlphaIntegration::IntegrateInterfaceShape(
  std::string condname)
{
  ParameterList eleparams;
  // set action for elements
  eleparams.set<int>("action",FLD::integrate_Shapefunction);

  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // create vector (+ initialization with zeros)
  Teuchos::RCP<Epetra_Vector> integratedshapefunc
    =
    LINALG::CreateVector(*dofrowmap,true);

  // call loop over elements
  discret_->ClearState();
  if (alefluid_)
  {
    discret_->SetState("dispnp", dispnp_);
  }
  discret_->EvaluateCondition(eleparams,integratedshapefunc,condname);
  discret_->ClearState();

  return integratedshapefunc;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::FluidGenAlphaIntegration::UseBlockMatrix(Teuchos::RCP<std::set<int> >     condelements,
                                                   const LINALG::MultiMapExtractor& domainmaps  ,
                                                   const LINALG::MultiMapExtractor& rangemaps   ,
                                                   bool                             splitmatrix )
{
  Teuchos::RCP<LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy> > mat;

  if (splitmatrix)
  {
    mat = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
                         domainmaps,
                         rangemaps,
                         500,
                         false,
                         true));
    //    mat->SetCondElements(condelements);
    sysmat_ = mat;
  }

  // if we never build the matrix nothing will be done
  if (params_->get<bool>("shape derivatives"))
  {
    // allocate special mesh moving matrix
    mat = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
                         domainmaps,
                         rangemaps,
                         500,
                         false,
                         true));
    //    mat->SetCondElements(condelements);
    shapederivatives_ = mat;
  }

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::FluidGenAlphaIntegration::LinearRelaxationSolve(
  Teuchos::RCP<Epetra_Vector> relax
  )
{
  dserror("No steepest descent for genalpha fluid\n");

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::FluidGenAlphaIntegration::AddDirichCond(const Teuchos::RCP<const Epetra_Map> maptoadd)
{
  std::vector<Teuchos::RCP<const Epetra_Map> > condmaps;
  condmaps.push_back(maptoadd);
  condmaps.push_back(dbcmaps_->CondMap());
  Teuchos::RCP<Epetra_Map> condmerged = LINALG::MultiMapExtractor::MergeMaps(condmaps);
  *dbcmaps_ = LINALG::MapExtractor(*(discret_->DofRowMap()), condmerged);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::FluidGenAlphaIntegration::RemoveDirichCond(const Teuchos::RCP<const Epetra_Map> maptoremove)
{
  std::vector<Teuchos::RCP<const Epetra_Map> > othermaps;
  othermaps.push_back(maptoremove);
  othermaps.push_back(dbcmaps_->OtherMap());
  Teuchos::RCP<Epetra_Map> othermerged = LINALG::MultiMapExtractor::MergeMaps(othermaps);
  *dbcmaps_ = LINALG::MapExtractor(*(discret_->DofRowMap()), othermerged, false);
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Teuchos::RCP<const Epetra_Vector> FLD::FluidGenAlphaIntegration::Dirichlet()
{
  if (dbcmaps_ == Teuchos::null)
    dserror("Dirichlet map has not been allocated");
  Teuchos::RCP<Epetra_Vector> dirichones = LINALG::CreateVector(*(dbcmaps_->CondMap()),false);
  dirichones->PutScalar(1.0);
  Teuchos::RCP<Epetra_Vector> dirichtoggle = LINALG::CreateVector(*(discret_->DofRowMap()),true);
  dbcmaps_->InsertCondVector(dirichones, dirichtoggle);
  return dirichtoggle;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Teuchos::RCP<const Epetra_Vector> FLD::FluidGenAlphaIntegration::InvDirichlet()
{
  if (dbcmaps_ == Teuchos::null)
    dserror("Dirichlet map has not been allocated");
  Teuchos::RCP<Epetra_Vector> dirichzeros = LINALG::CreateVector(*(dbcmaps_->CondMap()),true);
  Teuchos::RCP<Epetra_Vector> invtoggle = LINALG::CreateVector(*(discret_->DofRowMap()),false);
  invtoggle->PutScalar(1.0);
  dbcmaps_->InsertCondVector(dirichzeros, invtoggle);
  return invtoggle;
}

/*----------------------------------------------------------------------*
 | set scalar field at end of time step                       gjb 01/11 |
 *----------------------------------------------------------------------*/
void FLD::FluidGenAlphaIntegration::SetScalarField(
   RCP<const Epetra_Vector> scalarnp,
   Teuchos::RCP<DRT::Discretization> scatradis,
   const int whichscalar)
{
  // initializations
  int err(0);
  double value(0.0);
  vector<int> nodedofs;

  //--------------------------------------------------------------------------
  // Filling the scaaf-vector with scalar at time n+1 at pressure dofs
  //--------------------------------------------------------------------------
  // loop all nodes on the processor
  for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
  {
    // get the processor's local scatra node
    DRT::Node* lscatranode = scatradis->lRowNode(lnodeid);

    // find out the global dof id of the last(!) dof at the scatra node
    const int numscatradof = scatradis->NumDof(lscatranode);
    int globalscatradofid(-1);
    if (whichscalar == (-1))
    {
      // default: always take the last scatra dof for each node
      globalscatradofid = scatradis->Dof(lscatranode,numscatradof-1);
    }
    else
    {
      // respect the wish of the user
      globalscatradofid = scatradis->Dof(lscatranode,whichscalar);
    }
    const int localscatradofid = scalarnp->Map().LID(globalscatradofid);
    if (localscatradofid < 0)
      dserror("localdofid not found in map for given globaldofid");

    // get the processor's local fluid node
    DRT::Node* lnode = discret_->lRowNode(lnodeid);
    // get the global ids of degrees of freedom associated with this node
    nodedofs = discret_->Dof(lnode);
    // get global and processor's local pressure dof id (using the map!)
    const int globaldofid = nodedofs[numdim_];
    const int localdofid = scaaf_->Map().LID(globaldofid);
    if (localdofid < 0)
      dserror("localdofid not found in map for given globaldofid");

    value = (*scalarnp)[localscatradofid];
    err = scaaf_->ReplaceMyValue(localdofid,0,value);
    if (err != 0) dserror("error while inserting value into scaaf_");
  }

  return;

} // FluidGenAlphaIntegration::SetScalarField

// -------------------------------------------------------------------
// set general fluid parameter (AE 01/2011)
// -------------------------------------------------------------------

void FLD::FluidGenAlphaIntegration::SetElementGeneralFluidParameter()
{
  ParameterList eleparams;

  eleparams.set<int>("action",FLD::set_general_fluid_parameter);

  // set general element parameters
  eleparams.set("form of convective term","convective");
  eleparams.set<int>("Linearisation",newton_);
  eleparams.set<int>("Physical Type", physicaltype_);

  // parameter for stabilization
  eleparams.sublist("STABILIZATION") = params_->sublist("STABILIZATION");

  // timealgorithm is gen_alpha
  eleparams.set<int>("TimeIntegrationScheme", INPAR::FLUID::timeint_gen_alpha);

  // call standard loop over elements
  discret_->Evaluate(eleparams,null,null,null,null,null);
  return;
}

// -------------------------------------------------------------------
// set turbulence parameters                         rasthofer 11/2011
// -------------------------------------------------------------------
void FLD::FluidGenAlphaIntegration::SetElementTurbulenceParameter()
{
  ParameterList eleparams;

  eleparams.set<int>("action",FLD::set_turbulence_parameter);

  // set general parameters for turbulent flow
  eleparams.sublist("TURBULENCE MODEL") = params_->sublist("TURBULENCE MODEL");

  // set model-dependent parameters
  eleparams.sublist("SUBGRID VISCOSITY") = params_->sublist("SUBGRID VISCOSITY");
  eleparams.sublist("MULTIFRACTAL SUBGRID SCALES") = params_->sublist("MULTIFRACTAL SUBGRID SCALES");

  // call standard loop over elements
  discret_->Evaluate(eleparams,null,null,null,null,null);
  return;
}


// -------------------------------------------------------------------
// set general time parameter (AE 01/2011)
// -------------------------------------------------------------------

void FLD::FluidGenAlphaIntegration::SetElementTimeParameter()
{
  ParameterList eleparams;

  eleparams.set<int>("action",FLD::set_time_parameter);

  // set general element parameters
  eleparams.set("dt",dt_);
  // There is only one time integration scheme in FluidGenAlphaIntegration
  eleparams.set("total time",time_-(1-alphaF_)*dt_);
  eleparams.set("theta",alphaF_*gamma_/alphaM_);
  eleparams.set("omtheta",1-(alphaF_*gamma_/alphaM_));
  eleparams.set("alphaF",alphaF_);
  eleparams.set("alphaM",alphaM_);
  eleparams.set("gamma",gamma_);

  // call standard loop over elements
  discret_->Evaluate(eleparams,null,null,null,null,null);
  return;
}

// -------------------------------------------------------------------
// provide access to turbulence statistics manager (gjb 06/2011)
// -------------------------------------------------------------------
Teuchos::RCP<FLD::TurbulenceStatisticManager> FLD::FluidGenAlphaIntegration::TurbulenceStatisticManager()
  {return statisticsmanager_;}

// -------------------------------------------------------------------
// provide access to box filter for dynamic Smagorinsky model rasthofer
// -------------------------------------------------------------------
Teuchos::RCP<FLD::DynSmagFilter> FLD::FluidGenAlphaIntegration::DynSmagFilter() {return DynSmag_; }

// -------------------------------------------------------------------
// extrapolate from time mid-point to end-point         (mayr 12/2011)
// -------------------------------------------------------------------
Teuchos::RCP<Epetra_Vector> FLD::FluidGenAlphaIntegration::ExtrapolateEndPoint
(
  Teuchos::RCP<Epetra_Vector> vecn,
  Teuchos::RCP<Epetra_Vector> vecm
)
{
  Teuchos::RCP<Epetra_Vector> vecnp = Teuchos::rcp(new Epetra_Vector(*vecm));

  // For gen-alpha extrapolate mid-point quantities to end-point.
  vecnp->Update((alphaF_-1.0)/alphaF_,*vecn,1.0/alphaF_);

  return vecnp;
}

// -------------------------------------------------------------------
// return time integration factor                    (mayr.mt 03/2012)
// -------------------------------------------------------------------
double FLD::FluidGenAlphaIntegration::TimIntParam() const
{
  double retval = 0.0;
  switch (TimIntScheme())
  {
  case INPAR::FLUID::timeint_afgenalpha:
  case INPAR::FLUID::timeint_gen_alpha:
  case INPAR::FLUID::timeint_npgenalpha:
    // this is the interpolation weight for quantities from last time step
    retval = 1.0 - alphaF_;
  break;
  case INPAR::FLUID::timeint_one_step_theta:
  case INPAR::FLUID::timeint_bdf2:
  case INPAR::FLUID::timeint_stationary:
    dserror("OST, BDF2 and stationary time integration parameters are not defined in gen-alpha.");
  break;
  default:
    dserror("Unknown time integration scheme");
  break;
  }
  return retval;
}

// -------------------------------------------------------------------
// -------------------------------------------------------------------
void FLD::FluidGenAlphaIntegration::DisplacementToVelocity(
    Teuchos::RCP<Epetra_Vector> fcx,
    Teuchos::RCP<Epetra_Vector> ddgpre,
    Teuchos::RCP<Epetra_Vector> dugpre
)
{
#ifdef DEBUG
  // check, whether maps are the same
  if (! fcx->Map().SameAs(ddgpre->Map())) { dserror("Maps do not match, but they have to."); }
  if (! fcx->Map().SameAs(dugpre->Map())) { dserror("Maps do not match, but they have to."); }
#endif

  // get interface velocity at t(n)
  const Teuchos::RCP<Epetra_Vector> veln = Interface()->ExtractFSICondVector(Veln());

  /*
   * Delta u(n+1,i+1) = fac * [ Delta d(n+1,i+1) - dt * u_fluid(n) + Delta d_(predicted) ]
   *
   *                  - Delta u_fluid(predicted)
   *
   *             / = 2 / dt   if interface time integration is second order
   * with fac = |
   *             \ = 1 / dt   if interface time integration is first order
   */
  const double ts = TimeScaling();
  const double dt = Dt();
  fcx->Update(-dt*ts,*veln,ts,*ddgpre,ts);
  fcx->Update(-1.0,*dugpre,1.0);
}

// -------------------------------------------------------------------
// -------------------------------------------------------------------
void FLD::FluidGenAlphaIntegration::VelocityToDisplacement(
    Teuchos::RCP<Epetra_Vector> fcx,
    Teuchos::RCP<Epetra_Vector> ddgpre,
    Teuchos::RCP<Epetra_Vector> dugpre
)
{
#ifdef DEBUG
  // check, whether maps are the same
  if (! fcx->Map().SameAs(ddgpre->Map())) { dserror("Maps do not match, but they have to."); }
  if (! fcx->Map().SameAs(dugpre->Map())) { dserror("Maps do not match, but they have to."); }
#endif

  // get interface velocity at t(n)
  const Teuchos::RCP<Epetra_Vector> veln = Interface()->ExtractFSICondVector(Veln());

  /*
   * Delta d(n+1,i+1) = fac * [ Delta u(n+1,i+1) + Delta u(predicted)]
   *
   *                  + dt * u(n) - Delta d_structure(predicted)
   *
   *             / = dt / 2   if interface time integration is second order
   * with fac = |
   *             \ = dt       if interface time integration is first order
   */
  const double ts = 1.0/TimeScaling();
  fcx->Update(Dt(), *veln, ts, *dugpre, ts);
  fcx->Update(-1.0, *ddgpre, 1.0);
}
