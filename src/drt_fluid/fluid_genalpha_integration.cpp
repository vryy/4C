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
#ifdef CCADISCRET


#include "../drt_fluid/fluid_genalpha_integration.H"
#include "../drt_lib/drt_globalproblem.H"
#include "vm3_solver.H"
#include "fluid_utils.H"

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
  RefCountPtr<DRT::Discretization> actdis,
  LINALG::Solver&                  solver,
  ParameterList&                   params,
  IO::DiscretizationWriter&        output,
  bool                             alefluid
  )
  :
  // call constructor for "nontrivial" objects
  discret_(actdis),
  solver_ (solver),
  params_ (params),
  output_ (output),
  alefluid_(alefluid),
  density_(1.0),
  time_(0.0),
  step_(0),
  uprestart_(params.get("write restart every", -1)),
  upres_(params.get("write solution every", -1)),
  writestresses_(params.get("write stresses", 0))
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
  tm0_ref_        = rcp(new TimeMonitor(*timedyntot_ ));

  // time measurement --- start TimeMonitor tm7
  tm7_ref_        = rcp(new TimeMonitor(*timedyninit_ ));

  numdim_ = params_.get<int>("number of velocity degrees of freedom");

  // type of solver: low-Mach-number or incompressible solver
  loma_ = params_.get<string>("low-Mach-number solver","No");

  // -------------------------------------------------------------------
  // connect degrees of freedom for periodic boundary conditions
  // -------------------------------------------------------------------
  PeriodicBoundaryConditions::PeriodicBoundaryConditions pbc(discret_);
  pbc.UpdateDofsForPeriodicBoundaryConditions();

  pbcmapmastertoslave_ = pbc.ReturnAllCoupledNodesOnThisProc();

  // ensure that degrees of freedom in the discretization have been set
  if (!discret_->Filled()) discret_->FillComplete();

  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  // -------------------------------------------------------------------
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // -------------------------------------------------------------------
  // get a vector layout from the discretization for a vector which only
  // contains the velocity dofs and for one vector which only contains
  // pressure degrees of freedom.
  // -------------------------------------------------------------------

  FLD::UTILS::SetupFluidSplit(*discret_,numdim_,velpressplitter_);

  // -------------------------------------------------------------------
  // get the processor ID from the communicator
  // -------------------------------------------------------------------
  myrank_  = discret_->Comm().MyPID();

  // -------------------------------------------------------------------
  // create empty system matrix --- stiffness and mass are assembled in
  // one system matrix!
  // -------------------------------------------------------------------

  // This is a first estimate for the number of non zeros in a row of
  // the matrix. Assuming a structured 3d-fluid mesh we have 27 adjacent
  // nodes with 4 dofs each. (27*4=108)
  // We do not need the exact number here, just for performance reasons
  // a 'good' estimate

  sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,108,false,true));

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

  // velocity/density at time n+1
  vedeaf_       = LINALG::CreateVector(*dofrowmap,true);

  if (loma_ != "No")
  {
    // velocity/density at time n and n-1
    vedenp_     = LINALG::CreateVector(*dofrowmap,true);
    veden_      = LINALG::CreateVector(*dofrowmap,true);
  }

  // grid displacements and velocities for the ale case
  if (alefluid_)
  {
    dispnp_     = LINALG::CreateVector(*dofrowmap,true);
    dispn_      = LINALG::CreateVector(*dofrowmap,true);
    gridveln_   = LINALG::CreateVector(*dofrowmap,true);
    gridvelaf_  = LINALG::CreateVector(*dofrowmap,true);
  }


  // Vectors associated to boundary conditions
  // -----------------------------------------

  // toggle vector indicating which dofs have Dirichlet BCs
  dirichtoggle_ = LINALG::CreateVector(*dofrowmap,true);
  invdirtoggle_ = LINALG::CreateVector(*dofrowmap,true);

  // a vector of zeros to be used to enforce zero dirichlet boundary conditions
  zeros_        = LINALG::CreateVector(*dofrowmap,true);

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


  //--------------------------------------------------------------------
  // init some class variables (time integration)

  // time step size
  dt_     = params_.get<double>("time step size");
  // maximum number of timesteps
  endstep_= params_.get<int>   ("max number timesteps");
  // maximum simulation time
  endtime_= params_.get<double>("total time");

  // generalized alpha parameters
  // (choice of third parameter necessary but not sufficiant for second
  // order accuracy)
  //           gamma_  = 0.5 + alphaM_ - alphaF_
  alphaM_ = params_.get<double>("alpha_M");
  alphaF_ = params_.get<double>("alpha_F");
  gamma_  = params_.get<double>("gamma");

  // parameter for linearisation scheme (fixed point like or newton like)
  newton_   = params_.get<string>("Linearisation");

  itenum_   = 0;
  itemax_   = params_.get<int>   ("max nonlin iter steps");

  // -------------------------------------------------------------------
  // initialize turbulence-statistics evaluation
  // -------------------------------------------------------------------
  ParameterList *  modelparams =&(params_.sublist("TURBULENCE MODEL"));

  // flag for special flow: currently channel flow or flow in a lid-driven cavity
  special_flow_ = modelparams->get<string>("CANONICAL_FLOW","no");

  // all averaging is done in this statistics manager
  statisticsmanager_=rcp(new TurbulenceStatisticManager(*this));

  if (special_flow_ != "no")
  {
    // parameters for sampling/dumping period --- used for ad-hoc
    // modification of itemax for turbulent channel flows
    samstart_  = modelparams->get<int>("SAMPLING_START",1);
  }

  // -------------------------------------------------------------------
  // initialize outflow boundary stabilization if required
  // -------------------------------------------------------------------
  ParameterList *  stabparams=&(params_.sublist("STABILIZATION"));

  // flag for potential Neumann-type outflow stabilization
  outflow_stab_ = stabparams->get<string>("OUTFLOW_STAB","no_outstab");

  // the vector containing potential Neumann-type outflow stabilization
  if(outflow_stab_ == "yes_outstab")
    outflow_stabil_= LINALG::CreateVector(*dofrowmap,true);

  // (fine-scale) subgrid viscosity?
  fssgv_ = params_.get<string>("fs subgrid viscosity","No");

  // -------------------------------------------------------------------
  // necessary only for the VM3 approach:
  // coarse- and fine-scale solution vectors + respective ouptput
  // -------------------------------------------------------------------
  if (fssgv_ != "No")
  {
    csvelaf_  = LINALG::CreateVector(*dofrowmap,true);
    fsvelaf_  = LINALG::CreateVector(*dofrowmap,true);
    convaf_   = LINALG::CreateVector(*dofrowmap,true);
    csconvaf_ = LINALG::CreateVector(*dofrowmap,true);
    fsconvaf_ = LINALG::CreateVector(*dofrowmap,true);
  }
  // -------------------------------------------------------------------
  // initialise vectors for dynamic Smagorinsky model
  // (the smoothed quantities)
  // -------------------------------------------------------------------
  if (modelparams->get<string>("TURBULENCE_APPROACH","DNS_OR_RESVMM_LES")
      ==
      "CLASSICAL_LES")
  {
    if(modelparams->get<string>("PHYSICAL_MODEL","no_model")
       ==
       "Dynamic_Smagorinsky"
      )
    {
      
      // ---------------------------------------------------------------
      // get a vector layout from the discretization to construct
      // ---------------------------------------------------------------
      const Epetra_Map* noderowmap = discret_->NodeRowMap();

      filtered_vel_                   = rcp(new Epetra_MultiVector(*noderowmap,3,true));
      filtered_reynoldsstress_        = rcp(new Epetra_MultiVector(*noderowmap,9,true));
      filtered_modeled_subgrid_stress_= rcp(new Epetra_MultiVector(*noderowmap,9,true));

      averaged_LijMij_                = rcp(new vector<double>);
      averaged_MijMij_                = rcp(new vector<double>);
    }
  }

  //--------------------------------------------------------------------
  // do output to screen
  this->GenAlphaEchoToScreen("print start-up info");

  // end time measurement for timeloop

  // set vedenp-vector values to 1.0 for incompressible flow, for the time being
  if (loma_ == "No") vedeaf_->PutScalar(1.0);

  tm7_ref_ = null;

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
  tm2_ref_ = rcp(new TimeMonitor(*timedynloop_));

  bool stop_timeloop=false;
  while (stop_timeloop==false)
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
    tm1_ref_ = rcp(new TimeMonitor(*timeevaldirich_));

    this->GenAlphaApplyDirichletAndNeumann();

    // end time measurement for application of dirichlet conditions
    tm1_ref_=null;

    // -------------------------------------------------------------------
    //      calculate initial acceleration according to predicted
    //                  velocities and boundary values
    // -------------------------------------------------------------------
    this->GenAlphaCalcInitialAccelerations();

    // -------------------------------------------------------------------
    //                     solve nonlinear equation
    // -------------------------------------------------------------------
    this->DoGenAlphaPredictorCorrectorIteration();

    // -------------------------------------------------------------------
    //                         update solution
    // -------------------------------------------------------------------
    this->GenAlphaTimeUpdate();

    // -------------------------------------------------------------------
    // evaluate error for test flows with analytical solutions
    // -------------------------------------------------------------------
    this->EvaluateErrorComparedToAnalyticalSol();


    // time measurement --- start TimeMonitor tm8
    tm8_ref_        = rcp(new TimeMonitor(*timeout_ ));

    // -------------------------------------------------------------------
    // add calculated velocity to mean value calculation (statistics)
    // -------------------------------------------------------------------
    statisticsmanager_->DoTimeSample(step_,time_);

    // -------------------------------------------------------------------
    //                         output of solution
    // -------------------------------------------------------------------
    this->GenAlphaOutput();

    // time measurement --- stop TimeMonitor tm8
    tm8_ref_        = null;

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
  else itemax_  = params_.get<int>   ("max nonlin iter steps");

  // stop nonlinear iteration when both increment-norms are below this
  // bound
  ittol_     =params_.get<double>("tolerance for nonlin iter");

  // set
  if (params_.get<string>("CONVCHECK","L_2_norm")
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
  tm6_ref_ = rcp(new TimeMonitor(*timenlnloop_));

  // -------------------------------------------------------------------
  //  Evaluate acceleration and velocity at the intermediate time level
  //                     n+alpha_M and n+alpha_F
  //
  //                             -> (0)
  // -------------------------------------------------------------------
  // start time measurement for nonlinear update
  tm9_ref_ = rcp(new TimeMonitor(*timenonlinup_));

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
    Teuchos::RCP<Epetra_Vector> onlyvel = velpressplitter_.ExtractOtherVector(residual_);
    Teuchos::RCP<Epetra_Vector> onlypre = velpressplitter_.ExtractCondVector (residual_);

    // extract velocity and pressure residuals from rhs vector
    LINALG::Export(*residual_,*onlyvel);
    LINALG::Export(*residual_,*onlypre);

    onlypre->Norm2(&L2preresnorm_);
    onlyvel->Norm2(&L2velresnorm_);
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
    tm5_ref_ = rcp(new TimeMonitor(*timesolver_));

    // get cpu time
    tcpu=ds_cputime();

    this->GenAlphaCalcIncrement(badestnlnnorm);

    // end time measurement for application of dirichlet conditions
    tm5_ref_=null;
    dtsolve_=ds_cputime()-tcpu;

    // start time measurement for nonlinear update
    tm9_ref_ = rcp(new TimeMonitor(*timenonlinup_));

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
 |  use a constant predictor for the velocity and the pressure.         |
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

  //       n+1    n
  //      u    = u
  //       (0)
  //
  //  and
  //
  //       n+1    n
  //      p    = p
  //       (0)

  velnp_->Update(1.0,*veln_ ,0.0);


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

  ParameterList eleparams;

  // other parameters needed by the elements
  eleparams.set("total time",time_);
  eleparams.set("delta time",dt_  );
  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("velnp",velnp_);
  // predicted dirichlet values
  // velnp then also holds prescribed new dirichlet values
  // dirichtoggle is 1 for dirichlet dofs, 0 elsewhere
  discret_->EvaluateDirichlet(eleparams,velnp_,null,null,dirichtoggle_);
  discret_->ClearState();

  // --------------------------------------------------
  // evaluate Neumann conditions
  eleparams.set("total time",time_-(1-alphaF_)*dt_);
  eleparams.set("thsl",1.);

  discret_->SetState("vedenp",vedeaf_);
  neumann_loads_->PutScalar(0.0);
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

  // for velocities and pressure
  veln_->Update(1.0,*velnp_ ,0.0);
  // for the accelerations
  accn_->Update(1.0,*accnp_ ,0.0);


  if(params_.sublist("STABILIZATION").get<string>("TDS")=="time_dependent")
  {
    // create the parameters for the discretization
    ParameterList eleparams;
    // action for elements
    eleparams.set("action","time update for subscales");

    // update time paramters
    eleparams.set("gamma"  ,gamma_ );
    eleparams.set("dt"     ,dt_    );

    // call loop over elements
    discret_->Evaluate(eleparams,null,null,null,null,null);
  }

  if (alefluid_)
  {
    //            n+1   n
    //   .n+1    d   - d     gamma - 1   .n        .n
    //   d    = ---------- + --------- * d    ---> d
    //          gamma * dt     gamma
    //

    double gdtinv = 1.0/(gamma_*dt_);
    gridveln_->Update(gdtinv,*dispnp_,-gdtinv,*dispn_,(gamma_-1.0)/gamma_);

    //    n+1         n
    //   d      ---> d
    //

    dispn_   ->Update(1.0,*dispnp_,0.0);

  }


  return;
} // FluidGenAlphaIntegration::GenAlphaTimeUpdate


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
  //-------------------------------------------- output of solution
  if (step_%upres_ == 0)  //write solution
  {
    output_.NewStep    (step_,time_);

    output_.WriteVector("velnp"   , velnp_);

    // output real pressure
    Teuchos::RCP<Epetra_Vector> pressure = velpressplitter_.ExtractCondVector(velnp_);
    pressure->Scale(density_);
    output_.WriteVector("pressure", pressure);

    if (alefluid_)
    {
      output_.WriteVector("dispnp", dispnp_);
    }

    //only perform stress calculation when output is needed
    if (writestresses_)
    {
     RefCountPtr<Epetra_Vector> traction = CalcStresses();
     output_.WriteVector("traction",traction);
    }

    // do restart if we have to
    if (step_%uprestart_ == 0)
    {
      output_.WriteVector("veln",  veln_ );
      output_.WriteVector("accnp", accnp_);
      output_.WriteVector("accn",  accn_ );

      if (alefluid_)
      {
        output_.WriteVector("dispn"   ,dispn_   );
        output_.WriteVector("gridveln",gridveln_);
      }

      // write mesh in each restart step --- the elements are required since
      // they contain history variables (the time dependent subscales)
      output_.WriteMesh(step_,time_);
    }
  }
  // write restart also when uprestart_ is not a integer multiple of upres_
  else if (step_%uprestart_ == 0)
  {
    output_.NewStep    (step_,time_);

    output_.WriteVector("velnp", velnp_);
    output_.WriteVector("veln" , veln_ );
    output_.WriteVector("accnp", accnp_);
    output_.WriteVector("accn" , accn_ );

    if (alefluid_)
    {
      output_.WriteVector("dispnp", dispnp_);
      output_.WriteVector("dispn"   ,dispn_   );
      output_.WriteVector("gridveln",gridveln_);
    }

    //only perform stress calculation when output is needed
    if (writestresses_)
    {
     RefCountPtr<Epetra_Vector> traction = CalcStresses();
     output_.WriteVector("traction",traction);
    }

    // write mesh in each restart step --- the elements are required since
    // they contain history variables (the time dependent subscales)
    output_.WriteMesh(step_,time_);
  }

  // dumping of turbulence statistics if required
  statisticsmanager_->DoOutput(step_);

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
  if (params_.sublist("TURBULENCE MODEL").get<string>("TURBULENCE_APPROACH","DNS_OR_RESVMM_LES")
      ==
      "CLASSICAL_LES")
  {
    if(params_.sublist("TURBULENCE MODEL").get<string>("PHYSICAL_MODEL","no_model")
       ==
       "Dynamic_Smagorinsky"
      )
    {
      this->ApplyFilterForDynamicComputationOfCs();
    }
  }

  // get cpu time
  double tcpu=ds_cputime();

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
      rcp(new TimeMonitor(*timesparsitypattern_));

    sysmat_->Zero();

    timesparsitypattern_ref_=null;
  }

  // Neumann loads to residual
  residual_->Update(1.0,*neumann_loads_,0.0);

  // start time measurement for element call
  tm3_ref_ = rcp(new TimeMonitor(*timeeleloop_));

  // create the parameters for the discretization
  ParameterList eleparams;

  // add stabilization term at Neumann outflow boundary if required
  if(outflow_stab_ == "yes_outstab")
  {
    discret_->ClearState();
    discret_->SetState("velnp",velnp_);
    eleparams.set("thsl",1.);
    eleparams.set("outflow stabilization",outflow_stab_);

    discret_->SetState("vedeaf",vedeaf_);
    outflow_stabil_->PutScalar(0.0);
    discret_->EvaluateNeumann(eleparams,*outflow_stabil_);
    discret_->ClearState();

    // add Neumann-type stabilization term to residual vector
    residual_->Update(1.0,*outflow_stabil_,1.0);
  }

  // compute convective stresses when scale-similarity models are used
  if (fssgv_ == "scale_similarity" ||
      fssgv_ == "mixed_Smagorinsky_all" || fssgv_ == "mixed_Smagorinsky_small")
  {
    // action for elements
    eleparams.set("action","calc_convective_stresses");

    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("velaf",velaf_);
    discret_->SetState("vedeaf",vedeaf_);

    // element evaluation for getting convective stresses at nodes
    discret_->Evaluate(eleparams,null,convaf_);
    discret_->ClearState();
  }

  // action for elements
  eleparams.set("action","calc_fluid_genalpha_sysmat_and_residual");
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

  // parameters for nonlinear treatment (linearisation)
  eleparams.set("Linearisation",newton_);

  // parameters for stabilisation
  {
    eleparams.sublist("STABILIZATION")    = params_.sublist("STABILIZATION");
  }

  // parameters for a turbulence model
  {
    eleparams.sublist("TURBULENCE MODEL") = params_.sublist("TURBULENCE MODEL");
  }

  // (fine-scale) subgrid viscosity flag
  eleparams.set("fs subgrid viscosity",fssgv_);

  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("u and p (n+1      ,trial)",velnp_ );
  discret_->SetState("u and p (n+alpha_F,trial)",velaf_ );
  discret_->SetState("acc     (n+alpha_M,trial)",accam_ );
  discret_->SetState("vedeaf"                   ,vedeaf_);

  if (alefluid_)
  {
    discret_->SetState("dispnp"    , dispnp_   );
    discret_->SetState("gridvelaf" , gridvelaf_);
  }

  //----------------------------------------------------------------------
  // decide whether VM3-based solution approach or standard approach
  //----------------------------------------------------------------------
  if (fssgv_ != "No")
  {
    // extract the ML parameters
    ParameterList&  mllist = solver_.Params().sublist("ML Parameters");

    if (step_ == 1)
    {
      // zero fine-scale vector
      fsvelaf_->PutScalar(0.0);

      // set coarse- and fine-scale vectors
      discret_->SetState("fsvelaf",fsvelaf_);
      discret_->SetState("csvelaf",csvelaf_);
      discret_->SetState("csconvaf",csconvaf_);

      // element evaluation for getting system matrix
      discret_->Evaluate(eleparams,sysmat_,residual_);

      // complete system matrix
      sysmat_->Complete();

      // apply DBC to system matrix
      LINALG::ApplyDirichlettoSystem(sysmat_,increment_,residual_,zeros_,dirichtoggle_);

      // call VM3 constructor with system matrix
      vm3_solver_ = rcp(new VM3_Solver(SystemMatrix(),dirichtoggle_,mllist,true,true) );

      // zero system matrix again
      sysmat_->Zero();

      // add Neumann loads and potential Neumann-type outflow stabilization again
      residual_->Update(1.0,*neumann_loads_,0.0);
      if(outflow_stab_ == "yes_outstab")
        residual_->Update(1.0,*outflow_stabil_,1.0);
    }

    // check whether VM3 solver exists
    if (vm3_solver_ == null) dserror("vm3_solver not allocated");

    // call VM3 scale separation to get coarse- and fine-scale part of solution
    vm3_solver_->Separate(csvelaf_,fsvelaf_,velaf_);

    // call VM3 scale separation to get coarse-scale part of convective stresses
    if (fssgv_ == "scale_similarity" ||
        fssgv_ == "mixed_Smagorinsky_all" || fssgv_ == "mixed_Smagorinsky_small")
    {
      vm3_solver_->Separate(csconvaf_,fsconvaf_,convaf_);
      discret_->SetState("csconvaf",csconvaf_);
      discret_->SetState("csvelaf",csvelaf_);
    }

    // set fine-scale vector
    discret_->SetState("fsvelaf",fsvelaf_);

  }

  //----------------------------------------------------------------------
  // call loop over elements
  //----------------------------------------------------------------------
  // do not assemble the elemetn matrix if itmax is reached
  // in this case, only the residual is required for the convergence check
  if (itenum_<itemax_)
  {
    discret_->Evaluate(eleparams,sysmat_,residual_);
  }
  else
  {
    discret_->Evaluate(eleparams,Teuchos::null,residual_);
  }

  discret_->ClearState();

  // get density
  density_ = eleparams.get<double>("density");

  //----------------------------------------------------------------------
  // extended statistics (plane average of Cs) for dynamic 
  // Smagorinsky model --- communication part, store values
  //----------------------------------------------------------------------
  statisticsmanager_->StoreElementValues(step_);

  //----------------------------------------------------------------------
  // apply weak Dirichlet boundary conditions to sysmat_ and residual_
  //----------------------------------------------------------------------
  {
    ParameterList weakdbcparams;

    // set action for elements
    weakdbcparams.set("action","enforce_weak_dbc");
    weakdbcparams.set("afgdt",alphaF_*gamma_*dt_);
    weakdbcparams.set("total time",time_);

    // set the only required state vector
    discret_->SetState("u and p (n+alpha_F,trial)",velaf_);

    // evaluate
    discret_->EvaluateConditionUsingParentData
      (weakdbcparams      ,
       sysmat_            ,
       Teuchos::null      ,
       residual_          ,
       Teuchos::null      ,
       Teuchos::null      ,
       "LineWeakDirichlet");

    // clear state
    discret_->ClearState();
  }

  // end time measurement for element call
  tm3_ref_=null;
  
  // remember force vector for stress computation
  *force_=Epetra_Vector(*residual_);
  force_->Scale(density_);

  // start time measurement for generation of sparsity pattern
  {
    RefCountPtr<TimeMonitor> timesparsitypattern_ref_ = rcp(new TimeMonitor(*timesparsitypattern_));
    // finalize the system matrix
    sysmat_->Complete();

    timesparsitypattern_ref_ = null;
  }

  // -------------------------------------------------------------------
  // Apply dirichlet boundary conditions to system of equations residual
  // discplacements are supposed to be zero at boundary conditions
  // -------------------------------------------------------------------
  // start time measurement for application of dirichlet conditions
  tm4_ref_ = rcp(new TimeMonitor(*timeapplydirich_));

  zeros_->PutScalar(0.0);
  {
    LINALG::ApplyDirichlettoSystem(sysmat_,increment_,residual_,
                                   zeros_,dirichtoggle_);
  }

  // end time measurement for application of dirichlet conditions
  tm4_ref_=null;

  // end measurement element call
  dtele_=ds_cputime()-tcpu;


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
  bool   isadapttol    = params_.get<bool>  ("ADAPTCONV"                ,true  );
  double adaptolbetter = params_.get<double>("ADAPTCONV_BETTER"         ,0.01  );
  double nlntolsol     = params_.get<double>("tolerance for nonlin iter",1.e-10);

  //--------------------------- adapt tolerance  in the convergence limit
  if (isadapttol && itenum_>1)
  {
    solver_.AdaptTolerance(nlntolsol,nlnres,adaptolbetter);
  }

  //-------solve for residual displacements to correct incremental displacements
  increment_->PutScalar(0.0);

  solver_.Solve(sysmat_->EpetraOperator(),
                increment_               ,
                residual_                ,
                true                     ,
                itenum_==-1              );

  solver_.ResetTolerance();

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
  //               (i)
  accinc->Scale(gamma_*dt_);
  velpressplitter_.AddOtherVector(accinc,velnp_);

  // ------------------------------------------------------
  // update pressure
  //
  //         n+1          n+1
  //     pres      =  pres    + dpres
  //         (i+1)        (i)
  //
  if(numdim_==3)
  {
    // rescaled pressure to preserve symmetry of pressure
    // and continuity part in matrix
    preinc->Scale(gamma_*dt_);
    velpressplitter_.AddCondVector(preinc,velnp_);
    preinc->Scale(1.0/gamma_*dt_);
  }
  else
  {
    velpressplitter_.AddCondVector(preinc,velnp_);
  }
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
  if(numdim_==3)
  {
    L2incprenorm_*=gamma_*dt_;
  }

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

  reader.ReadVector(velnp_,"velnp");
  reader.ReadVector(veln_ ,"veln" );
  reader.ReadVector(accnp_,"accnp");
  reader.ReadVector(accn_ ,"accn" );

  if (alefluid_)
  {
       reader.ReadVector(dispnp_  ,"dispnp"  );
       reader.ReadVector(dispn_   ,"dispn"   );
       reader.ReadVector(gridveln_,"gridveln");
  }

  // read the previously written elements including the history data
  reader.ReadMesh(step_);

  return;
}


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
  int whichinitialfield,
  int startfuncno
 )
{

  //----------------------------------------------------------------------
  // Initialfield from function
  //----------------------------------------------------------------------
  if(whichinitialfield==2 ||whichinitialfield==3)
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

        double initialval=DRT::UTILS::FunctionManager::Instance().Funct(startfuncno-1).Evaluate(index,lnode->X());

        velnp_->ReplaceGlobalValues(1,&initialval,&gid);
        veln_ ->ReplaceGlobalValues(1,&initialval,&gid);
      }
    }
    //----------------------------------------------------------------------
    // random perturbations for field
    //----------------------------------------------------------------------
    if (whichinitialfield==3)
    {
      //
      {
        const Epetra_Map* dofrowmap = discret_->DofRowMap();

        int err =0;

        // random noise is perc percent of the initial profile

        double perc = params_.sublist("TURBULENCE MODEL").get<double>("CHAN_AMPL_INIT_DIST",0.1);

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

          // check whether a periodic boundary condition is active on this node
          if (mypbc.size()>0)
          {
            // yes, we have one

            // get the list of all his slavenodes
            map<int, vector<int> >::iterator master = pbcmapmastertoslave_.find(lnode->Id());

            // slavenodes are ignored
            if(master == pbcmapmastertoslave_.end())
            {
              // the node is a slave --- so don't do anything
              continue;
            }
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

            double randomnumber = 2*((double)rand()-((double) RAND_MAX)/2.)/((double) RAND_MAX);

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
  else if(whichinitialfield==8)
  {
    int gid;
    int lid;

    const Epetra_Map* dofrowmap = discret_->DofRowMap();


    int err =0;

    int numdim  = params_.get<int>("number of velocity degrees of freedom");
    int npredof = numdim;

    double         p;
    vector<double> u  (numdim);
    vector<double> xyz(numdim);


    if(numdim!=3)
    {
      dserror("Beltrami flow is three dimensional flow!");
    }

    // set constants for analytical solution
    double a      = PI/4.0;
    double d      = PI/2.0;

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

  int calcerr = params_.get<int>("eval err for analyt sol");

  //------------------------------------------------------- beltrami flow
  switch (calcerr)
  {
  case 0:
    // do nothing --- no analytical solution available
    break;
  case 2:
    // do nothing --- no analytical solution available
    break;
  case 3:
    // do nothing --- no analytical solution available
    break;
  case 8:
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
  double tcpu=ds_cputime();

{
    TEUCHOS_FUNC_TIME_MONITOR("FLD::FluidGenAlphaIntegration::BoxFilter");

    // generate a parameterlist for communication and control
    ParameterList filterparams;
    // action for elements
    filterparams.set("action","calc_fluid_box_filter");
    
    // set state vector to pass distributed vector to the element
    discret_->ClearState();
    discret_->SetState("u and p (trial)",velaf_);
      
    // define element matrices and vectors --- they are used to
    // transfer information into the element routine and back
    Epetra_SerialDenseMatrix ep_reystress_hat(3,3);
    Epetra_SerialDenseMatrix ep_modeled_stress_grid_scale_hat(3,3);
    Epetra_SerialDenseVector ep_velaf_hat (3);
    Epetra_SerialDenseVector dummy1;
    Epetra_SerialDenseVector dummy2;
    
    // ---------------------------------------------------------------
    // get a vector layout from the discretization to construct
    const Epetra_Map* noderowmap = discret_->NodeRowMap();
    
    // alloc an additional vector to store/add up the patch volume
    RCP<Epetra_Vector> patchvol     = rcp(new Epetra_Vector(*noderowmap,true));

    // zero out vecs
    filtered_vel_                   = rcp(new Epetra_MultiVector(*noderowmap,3,true));
    filtered_reynoldsstress_        = rcp(new Epetra_MultiVector(*noderowmap,9,true));
    filtered_modeled_subgrid_stress_= rcp(new Epetra_MultiVector(*noderowmap,9,true));

    // ---------------------------------------------------------------
    // loop all elements on this proc (including ghosted ones)
    for (int nele=0;nele<discret_->NumMyColElements();++nele)
    {
      // get the element
      DRT::Element* ele = discret_->lColElement(nele);
  
      // reset element matrices and vectors --- they are used to
      // transfer information into the element routine and back
      memset(ep_reystress_hat.A()                ,0,9*sizeof(double));
      memset(ep_modeled_stress_grid_scale_hat.A(),0,9*sizeof(double));
      memset(ep_velaf_hat.Values()               ,0,3*sizeof(double));
      
      // get element location vector, dirichlet flags and ownerships
      vector<int> lm;
      vector<int> lmowner;
      ele->LocationVector(*discret_,lm,lmowner);

      // call the element evaluate method to integrate functions
      // against heaviside function element
      int err = ele->Evaluate(filterparams,
                              *discret_,
                              lm,
                              ep_reystress_hat,
                              ep_modeled_stress_grid_scale_hat,
                              ep_velaf_hat,dummy1,dummy2);
      if (err) dserror("Proc %d: Element %d returned err=%d",
                       discret_->Comm().MyPID(),ele->Id(),err);

      // get contribution to patch volume of this element. Add it up.
      double volume_contribution =filterparams.get<double>("volume_contribution");

      // loop all nodes of this element, add values to the global vectors
      DRT::Node** elenodes=ele->Nodes();
      for(int nn=0;nn<ele->NumNode();++nn)
      {
        DRT::Node* node = (elenodes[nn]);
        
        // we are interested only in  row nodes
        if(node->Owner() == myrank_)
        {

          // now assemble the computed values into the global vector
          int    id = (node->Id());

          double val = volume_contribution;
          patchvol->SumIntoGlobalValues(1,&val,&id);

          for (int idim =0;idim<3;++idim)
          {
            val = ep_velaf_hat(idim);
            ((*filtered_vel_)(idim))->SumIntoGlobalValues(1,&val,&id);
            
            for (int jdim =0;jdim<3;++jdim)
            {
              val = ep_reystress_hat (idim,jdim);
              ((*filtered_reynoldsstress_ )       (3*idim+jdim))->SumIntoGlobalValues(1,&val,&id);
              
              val = ep_modeled_stress_grid_scale_hat(idim,jdim);
              ((*filtered_modeled_subgrid_stress_)(3*idim+jdim))->SumIntoGlobalValues(1,&val,&id);
            }
          }
        }
      }
    } // end elementloop

    // ---------------------------------------------------------------
    // send add values from masters and slaves
    {
      map<int, vector<int> >::iterator masternode;

      double val;
      double velaf_val[3];
      double reystress_val[3][3];
      double modeled_subgrid_stress_val[3][3];


      // loop all masternodes on this proc
      for(masternode =pbcmapmastertoslave_.begin();
          masternode!=pbcmapmastertoslave_.end();
          ++masternode)
      {
        // add all slave values to mastervalue
        vector<int>::iterator slavenode;
        
        int lid = noderowmap->LID(masternode->first);
        
        val = (*patchvol)[lid];
        
        for (int idim =0;idim<3;++idim)
        {
          velaf_val[idim]=((*((*filtered_vel_)(idim)))[lid]);
          
          for (int jdim =0;jdim<3;++jdim)
          {
            reystress_val             [idim][jdim]=(*((*filtered_reynoldsstress_         ) (3*idim+jdim)))[lid];
            modeled_subgrid_stress_val[idim][jdim]=(*((*filtered_modeled_subgrid_stress_ ) (3*idim+jdim)))[lid];
          }
        }
        
        // loop all this masters slaves
        for(slavenode=(masternode->second).begin();slavenode!=(masternode->second).end();++slavenode)
        {
          lid = noderowmap->LID(*slavenode);
          val += (*patchvol)[lid];
          for (int idim =0;idim<3;++idim)
          {
            velaf_val[idim]+=((*((*filtered_vel_)(idim)))[lid]);
            
            for (int jdim =0;jdim<3;++jdim)
            {
              reystress_val             [idim][jdim]+=(*((*filtered_reynoldsstress_         ) (3*idim+jdim)))[lid];
              modeled_subgrid_stress_val[idim][jdim]+=(*((*filtered_modeled_subgrid_stress_ ) (3*idim+jdim)))[lid];
            } // end loop jdim
          } // end loop idim
        }  // end loop slaves
        
        // replace value by sum
        lid = noderowmap->LID(masternode->first);
        patchvol->ReplaceMyValues(1,&val,&lid);
        
        for (int idim =0;idim<3;++idim)
        {
          ((*filtered_vel_)(idim))->ReplaceMyValues(1,&(velaf_val[idim]),&lid);
          
          for (int jdim =0;jdim<3;++jdim)
          {
            ((*filtered_reynoldsstress_)(3*idim+jdim))->ReplaceMyValues(1,&(reystress_val[idim][jdim]),&lid);
            ((*filtered_modeled_subgrid_stress_)(3*idim+jdim))->ReplaceMyValues(1,&(modeled_subgrid_stress_val[idim][jdim]),&lid);
          } // end loop jdim
        } // end loop idim

        // loop all this masters slaves
        for(slavenode=(masternode->second).begin();slavenode!=(masternode->second).end();++slavenode)
        {
          lid = noderowmap->LID(*slavenode);
          patchvol->ReplaceMyValues(1,&val,&lid);
          
          for (int idim =0;idim<3;++idim)
          {
            ((*filtered_vel_)(idim))->ReplaceMyValues(1,&(velaf_val[idim]),&lid);
            
            for (int jdim =0;jdim<3;++jdim)
            {
              ((*filtered_reynoldsstress_)(3*idim+jdim))->ReplaceMyValues(1,&(reystress_val[idim][jdim]),&lid);
              ((*filtered_modeled_subgrid_stress_)(3*idim+jdim))->ReplaceMyValues(1,&(modeled_subgrid_stress_val[idim][jdim]),&lid);
            } // end loop jdim
          } // end loop idim
        } // end loop slaves
      } // end loop masters
    }

    // ---------------------------------------------------------------
    // zero out dirichlet nodes
    {
      // get a rowmap for the dofs
      const Epetra_Map* dofrowmap = discret_->DofRowMap();
      
      // loop all nodes on the processor
      for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();++lnodeid)
      {
        // get the processor local node
        DRT::Node*  lnode       = discret_->lRowNode(lnodeid);
        
        // the set of degrees of freedom associated with the node
        vector<int> nodedofset = discret_->Dof(lnode);
        
        // check whether the node is on a wall, i.e. all velocity dofs
        // are Dirichlet constrained
        int is_no_slip_node =0;
        for(int index=0;index<numdim_;++index)
        {
          int gid = nodedofset[index];
          int lid = dofrowmap->LID(gid);
          
          if ((*dirichtoggle_)[lid]==1)
          {
            is_no_slip_node++;
          }
        }
        
        // this node is on a wall
        if (is_no_slip_node == numdim_)
        {
          for (int idim =0;idim<3;++idim)
          {
            double val = 0.0;
            ((*filtered_vel_)(idim))->ReplaceMyValues(1,&val,&lnodeid);
            
            for (int jdim =0;jdim<3;++jdim)
            {
              val = 0.0;
              ((*filtered_reynoldsstress_ ) (3*idim+jdim))->ReplaceMyValues(1,&val,&lnodeid);
              
              val = 0.0;
              ((*filtered_modeled_subgrid_stress_ ) (3*idim+jdim))->ReplaceMyValues(1,&val,&lnodeid);
            } // end loop jdim
          } // end loop idim
        }
      } // end loop all nodes
    }

    // ---------------------------------------------------------------
    // scale vectors by element patch sizes

    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();++lnodeid)
    {
      double thisvol = (*patchvol)[lnodeid];

      for (int idim =0;idim<3;++idim)
      {
        double val = ((*((*filtered_vel_)(idim)))[lnodeid])/thisvol;
        ((*filtered_vel_)(idim))->ReplaceMyValues(1,&val,&lnodeid);
        
        for (int jdim =0;jdim<3;++jdim)
        {
          val = ((*((*filtered_reynoldsstress_ ) (3*idim+jdim)))[lnodeid])/thisvol;
          ((*filtered_reynoldsstress_ ) (3*idim+jdim))->ReplaceMyValues(1,&val,&lnodeid);

          val = ((*((*filtered_modeled_subgrid_stress_ ) (3*idim+jdim)))[lnodeid])/thisvol;
          ((*filtered_modeled_subgrid_stress_ ) (3*idim+jdim))->ReplaceMyValues(1,&val,&lnodeid);
        } // end loop jdim
      } // end loop idim
    } // end loop nodes

    // clean up
    discret_->ClearState();
  }


  {
    TEUCHOS_FUNC_TIME_MONITOR("FLD::FluidGenAlphaIntegration::ComputeCs");

    // get the column map in order to communicate the result to all ghosted nodes
    const Epetra_Map* nodecolmap = discret_->NodeColMap();

    // allocate distributed vectors in col map format to have the filtered
    // quantities available on ghosted nodes
    col_filtered_vel_                    = rcp(new Epetra_MultiVector(*nodecolmap,3,true));
    col_filtered_reynoldsstress_         = rcp(new Epetra_MultiVector(*nodecolmap,9,true));
    col_filtered_modeled_subgrid_stress_ = rcp(new Epetra_MultiVector(*nodecolmap,9,true));

    // export filtered vectors in rowmap to columnmap format
    LINALG::Export(*filtered_vel_                   ,*col_filtered_vel_);
    LINALG::Export(*filtered_reynoldsstress_        ,*col_filtered_reynoldsstress_);
    LINALG::Export(*filtered_modeled_subgrid_stress_,*col_filtered_modeled_subgrid_stress_);
    
    // get ordered layers of elements in which LijMij and MijMij are averaged
    // the layers are obtained from the statistics evaluation via the
    // parameterlist (should be available after construction of the 
    // statistics manager)
    if (planecoords_ == null)
    {
      ParameterList *  modelparams =&(params_.sublist("TURBULENCE MODEL"));
      
      planecoords_=modelparams->get<RefCountPtr<vector<double> > >("planecoords_");
    }
    
    averaged_LijMij_->resize((*planecoords_).size()-1);
    averaged_MijMij_->resize((*planecoords_).size()-1);
    
    vector<int>  count_for_average      ((*planecoords_).size()-1);
    vector<int>  local_count_for_average((*planecoords_).size()-1);
    
    vector <double> local_ele_sum_LijMij((*planecoords_).size()-1);
    vector <double> local_ele_sum_MijMij((*planecoords_).size()-1);

    
    for (unsigned rr=0;rr<(*planecoords_).size()-1;++rr)
    {
      (*averaged_LijMij_)    [rr]=0;
      (*averaged_MijMij_)    [rr]=0;
      local_ele_sum_LijMij   [rr]=0;
      local_ele_sum_MijMij   [rr]=0;
      local_count_for_average[rr]=0;
    }

    // generate a parameterlist for communication and control
    ParameterList calc_smag_const_params;
    // action for elements
    calc_smag_const_params.set("action","calc_smagorinsky_const");
    
    // hand filtered global vectors down to the element
    calc_smag_const_params.set("col_filtered_vel"                   ,col_filtered_vel_);
    calc_smag_const_params.set("col_filtered_reynoldsstress"        ,col_filtered_reynoldsstress_);
    calc_smag_const_params.set("col_filtered_modeled_subgrid_stress",col_filtered_modeled_subgrid_stress_);
    
    // dummy matrices and vectors for element call
    Epetra_SerialDenseMatrix dummym1;
    Epetra_SerialDenseMatrix dummym2;
    Epetra_SerialDenseVector dummyv1;
    Epetra_SerialDenseVector dummyv2;
    Epetra_SerialDenseVector dummyv3;
    
    // loop all elements on this proc (excluding ghosted ones)
    for (int nele=0;nele<discret_->NumMyRowElements();++nele)
    {
      // get the element
      DRT::Element* ele = discret_->lRowElement(nele);
      
      // get element location vector, dirichlet flags and ownerships
      vector<int> lm;
      vector<int> lmowner;
      ele->LocationVector(*discret_,lm,lmowner);
      
      // call the element evaluate method to integrate functions
      // against heaviside function element
      int err = ele->Evaluate(calc_smag_const_params,
                              *discret_,
                              lm,
                              dummym1,dummym2,
                              dummyv1,dummyv2,dummyv3);
      if (err) dserror("Proc %d: Element %d returned err=%d",
                       discret_->Comm().MyPID(),ele->Id(),err);
      
    
      // get the result from the element call
      double LijMij = calc_smag_const_params.get<double>("LijMij");
      double MijMij = calc_smag_const_params.get<double>("MijMij");
      double center = calc_smag_const_params.get<double>("center");

      // add result into result vetor

      // for this purpose, determine the layer (the plane for average)
      int  nlayer;
      bool found = false;
      for (nlayer=0;nlayer<(int)(*planecoords_).size()-1;)
      {
        if(center<(*planecoords_)[nlayer+1])
        {
          found = true;
          break;
        }
        nlayer++;
      }
      if (found ==false)
      {
        dserror("could not determine element layer");
      }
      
      // add it up
      local_ele_sum_LijMij[nlayer] += LijMij;
      local_ele_sum_MijMij[nlayer] += MijMij;
      
      local_count_for_average[nlayer]++;
    } // end loop over elements

    // now add all the stuff from the different processors

    for (unsigned rr=0;rr<(*planecoords_).size()-1;++rr)
    {
      discret_->Comm().SumAll(&(local_count_for_average[rr]),&(count_for_average[rr])  ,1);
      discret_->Comm().SumAll(&(local_ele_sum_LijMij[rr])   ,&((*averaged_LijMij_)[rr]),1);
      discret_->Comm().SumAll(&(local_ele_sum_MijMij[rr])   ,&((*averaged_MijMij_)[rr]),1);
    }

    // do averaging
    for (unsigned rr=0;rr<(*planecoords_).size()-1;++rr)
    {
      (*averaged_LijMij_)[rr]/=count_for_average[rr];
      (*averaged_MijMij_)[rr]/=count_for_average[rr];
    }
    // provide necessary information for the elements
    {
      ParameterList *  modelparams =&(params_.sublist("TURBULENCE MODEL"));

      modelparams->set<RefCountPtr<vector<double> > >("averaged_LijMij_",averaged_LijMij_);
      modelparams->set<RefCountPtr<vector<double> > >("averaged_MijMij_",averaged_MijMij_);
      modelparams->set<RefCountPtr<vector<double> > >("planecoords_"    ,planecoords_    );
    }
  }

  dtfilter_=ds_cputime()-tcpu;

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
        cout << params_.get<double>("alpha_F");
        cout << "\n";

        // alpha_M
        cout << "                             alpha_M = ";
        cout << params_.get<double>("alpha_M");
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
	cout << "Linearisation              : ";
	cout << params_.get<string>("Linearisation");
        cout << endl;
      }

      //--------------------------------------------------------------------
      /* output of stabilisation details */
      {
        ParameterList *  stabparams=&(params_.sublist("STABILIZATION"));

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
      }

      //--------------------------------------------------------------------
      /* output of turbulence model if any */
      {
        ParameterList *  modelparams =&(params_.sublist("TURBULENCE MODEL"));

        if (modelparams->get<string>("TURBULENCE_APPROACH", "none")
            !=
            "none")
        {
          // a canonical flow with homogeneous directions would allow a
          // spatial averaging of data
          string special_flow_
            =
            modelparams->get<string>("CANONICAL_FLOW","no");

          string hom_plane
            =
            modelparams->get<string>("CHANNEL_HOMPLANE","not specified");

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
              cout << modelparams->get<double>("C_SMAGORINSKY") ;
            }
            else if(physmodel == "Smagorinsky_with_van_Driest_damping")
            {
              if (special_flow_ != "channel_flow_of_height_2"
                  ||
                  hom_plane != "xz")
              {
                dserror("The van Driest damping is only implemented for a channel flow with wall \nnormal direction y");
              }

              cout << "                             "          ;
              cout << "- Smagorinsky constant:   Cs   = "      ;
              cout << modelparams->get<double>("C_SMAGORINSKY");
              cout << endl;

              cout << "                             "          ;
              cout << "- viscous length      :   l_tau= "      ;
              cout << modelparams->get<double>("CHANNEL_L_TAU");
              cout << endl;
            }
            else if(physmodel == "Dynamic_Smagorinsky")
            {
              if (special_flow_ != "channel_flow_of_height_2"
                  ||
                  hom_plane != "xz")
              {
                dserror("The dynamic Smagorinsky model is only implemented for a channel flow with \nwall normal direction y");
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
          cout << params_.get<string>("fs subgrid viscosity");

          if (fssgv_ == "Smagorinsky_all" || fssgv_ == "Smagorinsky_small" ||
              fssgv_ == "mixed_Smagorinsky_all" || fssgv_ == "mixed_Smagorinsky_small")
          {
            cout << " with Smagorinsky constant Cs= ";
            cout << modelparams->get<double>("C_SMAGORINSKY") ;
          }
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
      if (params_.sublist("TURBULENCE MODEL").get<string>("TURBULENCE_APPROACH",
                                                          "DNS_OR_RESVMM_LES")
          ==
          "CLASSICAL_LES")
      {
        if(params_.sublist("TURBULENCE MODEL").get<string>("PHYSICAL_MODEL",
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
      if (params_.sublist("TURBULENCE MODEL").get<string>("TURBULENCE_APPROACH",
                                                          "DNS_OR_RESVMM_LES")
          ==
          "CLASSICAL_LES")
      {
        if(params_.sublist("TURBULENCE MODEL").get<string>("PHYSICAL_MODEL",
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
                             /       \     /            \
        n+af     alphaM     |  n+1  n |   |       alphaM |   .n
     u_G     = ---------- * | d   -d  | + | 1.0 - ------ | * d
               gamma * dt   |         |   |        gamma |
                             \       /     \            /

  */

  gridvelaf_->Update(1.0-alphaM_/gamma_,*gridveln_,0.0);

  gridvelaf_->Update( alphaM_/(gamma_*dt_),*dispnp_,
                     -alphaM_/(gamma_*dt_),*dispn_,1.0);

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

  FLD::UTILS::LiftDrag(*discret_,*force_,params_,liftdragvals);


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
  eleparams.set("action","integrate_Shapefunction");

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
void FLD::FluidGenAlphaIntegration::UseBlockMatrix(Teuchos::RCP<std::set<int> > condelements,
                                                   const LINALG::MultiMapExtractor& domainmaps,
                                                   const LINALG::MultiMapExtractor& rangemaps,
                                                   bool splitmatrix)
{
  if (splitmatrix)
  {
    Teuchos::RCP<LINALG::BlockSparseMatrix<FLD::UTILS::InterfaceSplitStrategy> > mat =
      Teuchos::rcp(new LINALG::BlockSparseMatrix<FLD::UTILS::InterfaceSplitStrategy>(domainmaps,rangemaps,108,false,true));
    mat->SetCondElements(condelements);
    sysmat_ = mat;
  }
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


#endif
