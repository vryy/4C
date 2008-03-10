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

#include "fluid_genalpha_integration.H"
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
FluidGenAlphaIntegration::FluidGenAlphaIntegration(
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
  restartstep_(0),
  upres_(params.get("write solution every", -1)),
  writestep_(0)
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
#ifdef PERF
  timeeleassemble_      = TimeMonitor::getNewTimer("       + time for assembly"     );
  timeelegetdoflocation_= TimeMonitor::getNewTimer("       + time to get lm"        );

  timeelegetvelnp_      = TimeMonitor::getNewTimer("        + get global vectors an set node values");
  timeeleinitsmag_      = TimeMonitor::getNewTimer("        + initialize Smagorinsky model");
  timeeleinitstab_      = TimeMonitor::getNewTimer("        + initialize stabilization flags");
  timeelesysmat_        = TimeMonitor::getNewTimer("        + complete call to sysmat");

  timeelederxy2_        = TimeMonitor::getNewTimer("         + second derivatives");
  timeelederxy_         = TimeMonitor::getNewTimer("         + first derivatives");
  timeeletau_           = TimeMonitor::getNewTimer("         + computation of tau");
  timeelegalerkin_      = TimeMonitor::getNewTimer("         + Galerkin part (matrix and rhs)");
  timeelepspg_          = TimeMonitor::getNewTimer("         + pspg part (matrix and rhs)");
  timeelesupg_          = TimeMonitor::getNewTimer("         + supg part (matrix and rhs)");
  timeelecstab_         = TimeMonitor::getNewTimer("         + cstab part (matrix and rhs)");
  timeelevstab_         = TimeMonitor::getNewTimer("         + vstab part (matrix and rhs)");
  timeelecrossrey_      = TimeMonitor::getNewTimer("         + cross and reynolds stresses");
  timeeleintertogp_     = TimeMonitor::getNewTimer("         + interpolation to gausspoint");
  timeeleseteledata_    = TimeMonitor::getNewTimer("         + blitz init, get bodyforce and xyze");
  timeeletdextras_      = TimeMonitor::getNewTimer("         + time dependent fine scale specials");
#endif

  // time measurement --- start TimeMonitor tm0
  tm0_ref_        = rcp(new TimeMonitor(*timedyntot_ ));

  // time measurement --- start TimeMonitor tm7
  tm7_ref_        = rcp(new TimeMonitor(*timedyninit_ ));

  numdim_ = params_.get<int>("number of velocity degrees of freedom");

  // -------------------------------------------------------------------
  // connect degrees of freedom for periodic boundary conditions
  // -------------------------------------------------------------------
  PeriodicBoundaryConditions::PeriodicBoundaryConditions pbc(discret_);
  pbc.UpdateDofsForPeriodicBoundaryConditions();

  mapmastertoslave_ = pbc.ReturnAllCoupledNodesOnThisProc();

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

  FLUID_UTILS::SetupFluidSplit(*discret_,numdim_,velpressplitter_);

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

  // grid displacements and velocities for the ale case
  if (alefluid_)
  {
    dispnp_     = LINALG::CreateVector(*dofrowmap,true);;
    dispn_      = LINALG::CreateVector(*dofrowmap,true);;
    gridveln_   = LINALG::CreateVector(*dofrowmap,true);;
    gridvelaf_  = LINALG::CreateVector(*dofrowmap,true);;
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

  // -------------------------------------------------------------------
  // initialize turbulence-statistics evaluation
  // -------------------------------------------------------------------
  ParameterList *  modelparams =&(params_.sublist("TURBULENCE MODEL"));

  // flag for special flow: currently channel flow or flow in a lid-driven cavity
  special_flow_ = modelparams->get<string>("CANONICAL_FLOW","no");

  if (special_flow_ != "no")
  {
    // parameters for sampling/dumping period
    samstart_  = modelparams->get<int>("SAMPLING_START",1);
    samstop_   = modelparams->get<int>("SAMPLING_STOP",1);
    dumperiod_ = modelparams->get<int>("DUMPING_PERIOD",1);

    if (special_flow_ == "lid_driven_cavity")
      turbulencestatistics_ldc_=rcp(new TurbulenceStatisticsLdc(discret_,params_));
    else if (special_flow_ == "channel_flow_of_height_2")
      turbulencestatistics_=rcp(new TurbulenceStatistics(discret_,params_));
  }

  // (fine-scale) subgrid viscosity?
  fssgv_ = params_.get<string>("fs subgrid viscosity","No");

  // -------------------------------------------------------------------
  // necessary only for the VM3 approach:
  // coarse- and fine-scale solution vectors + respective ouptput
  // -------------------------------------------------------------------
  if (fssgv_ != "No")
  {
    csvelaf_ = LINALG::CreateVector(*dofrowmap,true);
    fsvelaf_ = LINALG::CreateVector(*dofrowmap,true);
  }

  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct
  // -------------------------------------------------------------------
  const Epetra_Map* noderowmap = discret_->NodeRowMap();

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
      filtered_vel_                   = rcp(new Epetra_MultiVector(*noderowmap,3,true));
      filtered_reynoldsstress_        = rcp(new Epetra_MultiVector(*noderowmap,9,true));
      filtered_modeled_subgrid_stress_= rcp(new Epetra_MultiVector(*noderowmap,9,true));

      averaged_LijMij_                = rcp(new vector<double>);
      averaged_MijMij_                = rcp(new vector<double>);
    }
  }

  this->GenAlphaEchoToScreen("print start-up info");

  // end time measurement for timeloop


  //--------------------------------------------------------------------
  // init some class variables

  dt_     = params_.get<double>("time step size");

  alphaM_ = params_.get<double>("alpha_M");
  alphaF_ = params_.get<double>("alpha_F");

  // choice of third parameter necessary but not sufficiant for second
  // order accuracy
  gamma_  = 0.5 + alphaM_ - alphaF_;

  // parameter for linearisation scheme (fixed point like or newton like)
  newton_ = params_.get<bool>("Use reaction terms for linearisation",false);

  // maximum number of timesteps
  endstep_  = params_.get<int>   ("max number timesteps");
  // maximum simulation time
  endtime_  = params_.get<double>("total time");


  tm7_ref_ = null;

  return;

} // FluidGenAlphaIntegration::FluidGenAlphaIntegration


/*----------------------------------------------------------------------*
 | Destructor dtor (public)                                  gammi 06/07|
 *----------------------------------------------------------------------*/
FluidGenAlphaIntegration::~FluidGenAlphaIntegration()
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

void FluidGenAlphaIntegration::GenAlphaTimeloop()
{
  //--------------------------------------------------------------------
  // do output to screen
  this->GenAlphaEchoToScreen("print start-up info");


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
    if(special_flow_ != "no" && step_>=samstart_ && step_<=samstop_)
    {
      if(special_flow_ == "lid_driven_cavity")
        turbulencestatistics_ldc_->DoTimeSample(velnp_);
      else if(special_flow_ == "channel_flow_of_height_2")
        this->GenAlphaTakeSample();
    }

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

void FluidGenAlphaIntegration::DoGenAlphaPredictorCorrectorIteration(
  )
{
  double            tcpu     ;

  // reset iteration counter, initialize max iteration counter
  itenum_  = 0;
  itemax_  = params_.get<int>   ("max nonlin iter steps");

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
void FluidGenAlphaIntegration::GenAlphaPredictNewSolutionValues()
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
void FluidGenAlphaIntegration::GenAlphaApplyDirichletAndNeumann()
{
  // --------------------------------------------------
  // apply Dirichlet conditions to velnp

  ParameterList eleparams;

  // choose what to assemble
  eleparams.set("assemble matrix 1",false);
  eleparams.set("assemble matrix 2",false);
  eleparams.set("assemble vector 1",true);
  eleparams.set("assemble vector 2",false);
  eleparams.set("assemble vector 3",false);
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

  // not implemented yet

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
void FluidGenAlphaIntegration::GenAlphaCalcInitialAccelerations()
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
void FluidGenAlphaIntegration::GenAlphaComputeIntermediateSol()
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
void FluidGenAlphaIntegration::GenAlphaTimeUpdate()
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

    // choose what to assemble --- nothing
    eleparams.set("assemble matrix 1",false);
    eleparams.set("assemble matrix 2",false);
    eleparams.set("assemble vector 1",false);
    eleparams.set("assemble vector 2",false);
    eleparams.set("assemble vector 3",false);

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
void FluidGenAlphaIntegration::GenAlphaOutput()
{
  //-------------------------------------------- output of solution
  restartstep_ += 1;
  writestep_ += 1;

  if (writestep_ == upres_)  //write solution
  {
    writestep_= 0;
    output_.NewStep    (step_,time_);

    output_.WriteVector("velnp"   , velnp_);

    if (alefluid_)
    {
      output_.WriteVector("dispnp", dispnp_);
    }

    // do restart if we have to
    if (restartstep_ == uprestart_)
    {
      restartstep_ = 0;

      output_.WriteVector("veln ", veln_ );
      output_.WriteVector("accnp", accnp_);
      output_.WriteVector("accn ", accn_ );

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
  if ((restartstep_ == uprestart_) && (writestep_ > 0))
  {
    restartstep_ = 0;

    output_.NewStep    (step_,time_);

    output_.WriteVector("velnp", velnp_);
    output_.WriteVector("veln ", veln_ );
    output_.WriteVector("accnp", accnp_);
    output_.WriteVector("accn ", accn_ );

    if (alefluid_)
    {
      output_.WriteVector("dispnp", dispnp_);
      output_.WriteVector("dispn"   ,dispn_   );
      output_.WriteVector("gridveln",gridveln_);
    }


    // write mesh in each restart step --- the elements are required since
    // they contain history variables (the time dependent subscales)
    output_.WriteMesh(step_,time_);

    if(special_flow_ == "channel_flow_of_height_2"  && dumperiod_ == 0)
    {
      turbulencestatistics_->TimeAverageMeansAndOutputOfStatistics(step_);
      turbulencestatistics_->ClearStatistics();
    }
  }

  // dumping of turbulence statistics if required
  if (special_flow_ != "no"  && dumperiod_ != 0)
  {
    int samstep = step_-samstart_+1;
    double dsamstep=samstep;
    double ddumperiod=dumperiod_;

    if (fmod(dsamstep,ddumperiod)==0)
    {
      if (special_flow_ == "lid_driven_cavity")
        turbulencestatistics_ldc_->DumpStatistics(step_);
      else if (special_flow_ == "channel_flow_of_height_2")
        turbulencestatistics_->DumpStatistics(step_);
    }
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
void FluidGenAlphaIntegration::GenAlphaAssembleResidualAndMatrix()
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

  // zero out residual
  residual_->PutScalar(0.0);

  // add Neumann loads to residual
  residual_->Update(1.0,*neumann_loads_,0.0);


  // start time measurement for element call
  tm3_ref_ = rcp(new TimeMonitor(*timeeleloop_));

  // create the parameters for the discretization
  ParameterList eleparams;

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

  // do not compute the elemetn matrix if itmax is reached
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
  eleparams.set("include reactive terms for linearisation"    ,newton_);

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

#ifdef PERF
  // pass time monitor to element for detailed time measurement

  eleparams.set<RefCountPtr<Teuchos::Time> >("time for assembly",timeeleassemble_      );
  eleparams.set<RefCountPtr<Teuchos::Time> >("time to get lm"   ,timeelegetdoflocation_);

  eleparams.set<RefCountPtr<Teuchos::Time> >("get global vectors an set node values"    ,timeelegetvelnp_  );
  eleparams.set<RefCountPtr<Teuchos::Time> >("initialise Smagorinsky model"             ,timeeleinitsmag_  );
  eleparams.set<RefCountPtr<Teuchos::Time> >("initialise stabilization flags"           ,timeeleinitstab_  );
  eleparams.set<RefCountPtr<Teuchos::Time> >("complete call to sysmat"                  ,timeelesysmat_    );

  eleparams.set<RefCountPtr<Teuchos::Time> >("time used for second derivatives"         ,timeelederxy2_    );
  eleparams.set<RefCountPtr<Teuchos::Time> >("time used for secondfirst"                ,timeelederxy_     );
  eleparams.set<RefCountPtr<Teuchos::Time> >("time used for computation of tau"         ,timeeletau_       );
  eleparams.set<RefCountPtr<Teuchos::Time> >("time used for galerkin loops"             ,timeelegalerkin_  );
  eleparams.set<RefCountPtr<Teuchos::Time> >("time used for pspg loop"                  ,timeelepspg_      );
  eleparams.set<RefCountPtr<Teuchos::Time> >("time used for supg loop"                  ,timeelesupg_      );
  eleparams.set<RefCountPtr<Teuchos::Time> >("time used for cstab loop"                 ,timeelecstab_     );
  eleparams.set<RefCountPtr<Teuchos::Time> >("time used for vstab loop"                 ,timeelevstab_     );
  eleparams.set<RefCountPtr<Teuchos::Time> >("time used for cross and reynolds loop"    ,timeelecrossrey_  );
  eleparams.set<RefCountPtr<Teuchos::Time> >("time used to interpolate to gauss points" ,timeeleintertogp_ );
  eleparams.set<RefCountPtr<Teuchos::Time> >("set basic element data"                   ,timeeleseteledata_);
  eleparams.set<RefCountPtr<Teuchos::Time> >("time dependent fine scale specials"       ,timeeletdextras_  );
#endif


  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("u and p (n+1      ,trial)",velnp_);
  discret_->SetState("u and p (n+alpha_F,trial)",velaf_);
  discret_->SetState("acc     (n+alpha_M,trial)",accam_);

  if (alefluid_)
  {
    discret_->SetState("dispnp"    , dispnp_   );
    discret_->SetState("gridvelaf" , gridvelaf_);
  }


  // extended statistics (plane average of Cs, (Cs_delta)^2, visceff)
  // for dynamic Smagorinsky model

  RefCountPtr<vector<double> > global_incr_Cs_sum;
  RefCountPtr<vector<double> > local_Cs_sum;
  global_incr_Cs_sum =  rcp(new vector<double> );
  local_Cs_sum       =  rcp(new vector<double> );


  RefCountPtr<vector<double> > global_incr_Cs_delta_sq_sum;
  RefCountPtr<vector<double> > local_Cs_delta_sq_sum;
  global_incr_Cs_delta_sq_sum =  rcp(new vector<double> );
  local_Cs_delta_sq_sum       =  rcp(new vector<double> );


  RefCountPtr<vector<double> > global_incr_visceff_sum;
  RefCountPtr<vector<double> > local_visceff_sum;
  global_incr_visceff_sum =  rcp(new vector<double> );
  local_visceff_sum       =  rcp(new vector<double> );


  if (params_.sublist("TURBULENCE MODEL").get<string>("TURBULENCE_APPROACH","DNS_OR_RESVMM_LES")
      ==
      "CLASSICAL_LES")
  {
    if(params_.sublist("TURBULENCE MODEL").get<string>("PHYSICAL_MODEL","no_model")
       ==
       "Dynamic_Smagorinsky"
       ||
       params_.sublist("TURBULENCE MODEL").get<string>("PHYSICAL_MODEL","no_model")
       ==
       "Smagorinsky_with_van_Driest_damping"
      )
    {
      // get ordered layers of elements in which LijMij and MijMij are averaged
      if (planecoords_ == null)
      {
        planecoords_ = rcp( new vector<double>((turbulencestatistics_->ReturnNodePlaneCoords()).size()));
        (*planecoords_) = turbulencestatistics_->ReturnNodePlaneCoords();
      }

      local_Cs_sum->resize      (planecoords_->size()-1,0.0);
      global_incr_Cs_sum->resize(planecoords_->size()-1,0.0);

      local_Cs_delta_sq_sum->resize      (planecoords_->size()-1,0.0);
      global_incr_Cs_delta_sq_sum->resize(planecoords_->size()-1,0.0);

      local_visceff_sum->resize      (planecoords_->size()-1,0.0);
      global_incr_visceff_sum->resize(planecoords_->size()-1,0.0);

      eleparams.sublist("TURBULENCE MODEL").set<RefCountPtr<vector<double> > >("planecoords_",planecoords_);
      eleparams.sublist("TURBULENCE MODEL").set<RefCountPtr<vector<double> > >("local_Cs_sum",local_Cs_sum);
      eleparams.sublist("TURBULENCE MODEL").set<RefCountPtr<vector<double> > >("local_Cs_delta_sq_sum",local_Cs_delta_sq_sum);
      eleparams.sublist("TURBULENCE MODEL").set<RefCountPtr<vector<double> > >("local_visceff_sum",local_visceff_sum);
    }
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

      // set fine-scale vector
      discret_->SetState("fsvelaf",fsvelaf_);

      // element evaluation for getting system matrix
      discret_->Evaluate(eleparams,sysmat_,residual_);

      // complete system matrix
      sysmat_->Complete();

      // apply DBC to system matrix
      LINALG::ApplyDirichlettoSystem(sysmat_,increment_,residual_,zeros_,dirichtoggle_);

      // call VM3 constructor with system matrix
      vm3_solver_ = rcp(new VM3_Solver(sysmat_,dirichtoggle_,mllist,true,true) );

      // zero system matrix again
      sysmat_->Zero();

      // add Neumann loads again
      residual_->Update(1.0,*neumann_loads_,0.0);
    }

    // check whether VM3 solver exists
    if (vm3_solver_ == null) dserror("vm3_solver not allocated");

    // call the VM3 scale separation to get fine-scale part of solution
    vm3_solver_->Separate(fsvelaf_,velaf_);

    // set fine-scale vector
    discret_->SetState("fsvelaf",fsvelaf_);

  }

  //----------------------------------------------------------------------
  // call loop over elements
  //----------------------------------------------------------------------
  discret_->Evaluate(eleparams,sysmat_,residual_);
  discret_->ClearState();

  // get density
  density_ = eleparams.get<double>("density");

  // end time measurement for element call
  tm3_ref_=null;

  // extended statistics (plane average of Cs) for dynamic Smagorinsky model --- communication part
  if (params_.sublist("TURBULENCE MODEL").get<string>("TURBULENCE_APPROACH","DNS_OR_RESVMM_LES")
      ==
      "CLASSICAL_LES")
  {
    if(params_.sublist("TURBULENCE MODEL").get<string>("PHYSICAL_MODEL","no_model")
       ==
       "Dynamic_Smagorinsky"
       ||
       params_.sublist("TURBULENCE MODEL").get<string>("PHYSICAL_MODEL","no_model")
       ==
       "Smagorinsky_with_van_Driest_damping"
      )
    {
      // now add all the stuff from the different processors
      discret_->Comm().SumAll(&((*local_Cs_sum               )[0]),
                              &((*global_incr_Cs_sum         )[0]),
                              local_Cs_sum->size());
      discret_->Comm().SumAll(&((*local_Cs_delta_sq_sum      )[0]),
                              &((*global_incr_Cs_delta_sq_sum)[0]),
                              local_Cs_delta_sq_sum->size());
      discret_->Comm().SumAll(&((*local_visceff_sum          )[0]),
                              &((*global_incr_visceff_sum    )[0]),
                              local_visceff_sum->size());
    }
    turbulencestatistics_->ReplaceCsIncrement(global_incr_Cs_sum,
                                              global_incr_Cs_delta_sq_sum,
                                              global_incr_visceff_sum);
  }

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
void FluidGenAlphaIntegration::GenAlphaCalcIncrement(const double nlnres)
{
  bool   isadapttol    = params_.get<bool>("ADAPTCONV",true);
  double adaptolbetter = params_.get<double>("ADAPTCONV_BETTER",0.01);
  double nlntolsol     = params_.get<double>("tolerance for nonlin iter",1.e-10);

  //--------------------------- adapt tolerance  in the convergence limit
  if (isadapttol && itenum_>1) solver_.AdaptTolerance(nlntolsol,nlnres,adaptolbetter);

  //-------solve for residual displacements to correct incremental displacements
  increment_->PutScalar(0.0);
  solver_.Solve(sysmat_->EpetraMatrix(),increment_,residual_,true,itenum_==-1);
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
void FluidGenAlphaIntegration::GenAlphaNonlinearUpdate()
{
  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  // -------------------------------------------------------------------
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  int numlocdofs = dofrowmap->NumMyElements();

  int*   dof = dofrowmap->MyGlobalElements ();


  int predof = numdim_+1;

  // loop all dofs on this proc
  for (int lid=0; lid<numlocdofs; ++lid)
  {
    int gid = dof[lid];

    // if the dof is belonging to an acceleration/velocity
    if ((gid+1)%predof != 0)
    {
      // ------------------------------------------------------
      // update acceleration
      //
      //        n+1         n+1
      //     acc      =  acc    + dacc
      //        (i+1)       (i)
      //
      double dacc = (*increment_)[lid];

      accnp_->SumIntoGlobalValues(1,&dacc,&gid);

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

      //
      double dvel = gamma_*dt_*dacc;

      velnp_->SumIntoGlobalValues(1,&dvel,&gid);
    }
    else
    {
      // ------------------------------------------------------
      // update pressure
      //
      //         n+1          n+1
      //     pres      =  pres    + dpres
      //         (i+1)        (i)
      //

      double dpres = (*increment_)[lid];
      velnp_->SumIntoGlobalValues(1,&dpres,&gid);

    }
  }

#if 0  // DEBUG IO  --- rhs of linear system
  {
    const Epetra_Map* dofrowmap       = discret_->DofRowMap();

    int rr;


    double* data = residual_->Values();
    for(rr=0;rr<residual_->MyLength();rr++)
    {
      int  gid = dofrowmap->GID(rr);

      if(gid%4==0)
      printf("%4d vel %22.15e\n",gid,data[rr]);
    }
  }

#endif


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
bool FluidGenAlphaIntegration::GenAlphaNonlinearConvergenceCheck(double& badestnlnnorm)
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
void FluidGenAlphaIntegration::ReadRestart(int step)
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
void FluidGenAlphaIntegration::SetInitialFlowField(
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

        double perc = 0.1;

        if (params_.sublist("TURBULENCE MODEL").get<string>("CANONICAL_FLOW","no")
            ==
            "channel_flow_of_height_2")
        {
          perc = params_.sublist("TURBULENCE MODEL").get<double>("CHAN_AMPL_INIT_DIST");
        }

        // out to screen
        if (myrank_==0)
        {
          cout << "Disturbed initial profile:   max. " << perc << "\% random perturbation\n";
          cout << "\n\n";
        }

        // loop all nodes on the processor
        for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();++lnodeid)
        {
          // get the processor local node
          DRT::Node*  lnode      = discret_->lRowNode(lnodeid);
          // the set of degrees of freedom associated with the node
          vector<int> nodedofset = discret_->Dof(lnode);

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
void FluidGenAlphaIntegration::EvaluateErrorComparedToAnalyticalSol()
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
    // choose what to assemble --- nothing
    eleparams.set("assemble matrix 1",false);
    eleparams.set("assemble matrix 2",false);
    eleparams.set("assemble vector 1",false);
    eleparams.set("assemble vector 2",false);
    eleparams.set("assemble vector 3",false);
    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("u and p at time n+1 (converged)",velnp_);

    // call loop over elements
    discret_->Evaluate(eleparams,sysmat_,null,residual_,null,null);
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
 | evaluate error for test cases with analytical solutions   gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FluidGenAlphaIntegration::GenAlphaTakeSample()
{
  // --------------------------------------------------------------------
  // add up X, X^2 for velocities and pressure
  turbulencestatistics_->DoTimeSample(velnp_,*force_);

  // create the parameters for the discretization
  ParameterList eleparams;

  // --------------------------------------------------------------------
  // do time averaging for subscales and residual --- this is of interest
  // only for time dependent subscales

  // action for elements
  eleparams.set("action","time average for subscales and residual");

      eleparams.set("assemble matrix 1",false);
      eleparams.set("assemble matrix 2",false);
      eleparams.set("assemble vector 1",false);
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

      // parameters for stabilisation
      {
        eleparams.sublist("STABILIZATION")    = params_.sublist("STABILIZATION");
      }

      // set vector values needed by elements
      discret_->ClearState();
      discret_->SetState("u and p (n+1      ,trial)",velnp_);
      discret_->SetState("u and p (n+alpha_F,trial)",velaf_);
      discret_->SetState("acc     (n+alpha_M,trial)",accam_);

      // get ordered layers of elements in which LijMij and MijMij are averaged
      if (planecoords_ == null)
      {
        planecoords_ = rcp( new vector<double>((turbulencestatistics_->ReturnNodePlaneCoords()).size()));
      }

      (*planecoords_) = turbulencestatistics_->ReturnNodePlaneCoords();

      RefCountPtr<vector<double> > local_incrres;
      local_incrres=  rcp(new vector<double> );
      local_incrres->resize(3*(planecoords_->size()-1),0.0);

      RefCountPtr<vector<double> > global_incrres;
      global_incrres=  rcp(new vector<double> );
      global_incrres->resize(3*(planecoords_->size()-1),0.0);

      RefCountPtr<vector<double> > local_incrres_sq;
      local_incrres_sq=  rcp(new vector<double> );
      local_incrres_sq->resize(3*(planecoords_->size()-1),0.0);

      RefCountPtr<vector<double> > global_incrres_sq;
      global_incrres_sq=  rcp(new vector<double> );
      global_incrres_sq->resize(3*(planecoords_->size()-1),0.0);

      RefCountPtr<vector<double> > local_incrsacc;
      local_incrsacc=  rcp(new vector<double> );
      local_incrsacc->resize(3*(planecoords_->size()-1),0.0);

      RefCountPtr<vector<double> > global_incrsacc;
      global_incrsacc=  rcp(new vector<double> );
      global_incrsacc->resize(3*(planecoords_->size()-1),0.0);

      RefCountPtr<vector<double> > local_incrsacc_sq;
      local_incrsacc_sq=  rcp(new vector<double> );
      local_incrsacc_sq->resize(3*(planecoords_->size()-1),0.0);

      RefCountPtr<vector<double> > global_incrsacc_sq;
      global_incrsacc_sq=  rcp(new vector<double> );
      global_incrsacc_sq->resize(3*(planecoords_->size()-1),0.0);

      eleparams.set<RefCountPtr<vector<double> > >("planecoords_",planecoords_);
      eleparams.set<RefCountPtr<vector<double> > >("incrres"    ,local_incrres);
      eleparams.set<RefCountPtr<vector<double> > >("incrres_sq" ,local_incrres_sq);
      eleparams.set<RefCountPtr<vector<double> > >("incrsacc"   ,local_incrsacc);
      eleparams.set<RefCountPtr<vector<double> > >("incrsacc_sq",local_incrsacc_sq);

      // call loop over elements
      {
        discret_->Evaluate(eleparams,null,null,null,null,null);
        discret_->ClearState();
      }

      discret_->Comm().SumAll(&((*local_incrres)[0])    ,&((*global_incrres)[0])    ,3*(planecoords_->size()-1));
      discret_->Comm().SumAll(&((*local_incrres_sq)[0]) ,&((*global_incrres_sq)[0]) ,3*(planecoords_->size()-1));
      discret_->Comm().SumAll(&((*local_incrsacc)[0])   ,&((*global_incrsacc)[0])   ,3*(planecoords_->size()-1));
      discret_->Comm().SumAll(&((*local_incrsacc_sq)[0]),&((*global_incrsacc_sq)[0]),3*(planecoords_->size()-1));

      turbulencestatistics_->AddToResAverage(global_incrres,
                                             global_incrres_sq,
                                             global_incrsacc,
                                             global_incrsacc_sq);
      return;
}



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
void FluidGenAlphaIntegration::ApplyFilterForDynamicComputationOfCs()
{

  // time measurement
  double tcpu=ds_cputime();

  // get the dofrowmap for access of velocity dofs
  const Epetra_Map* dofrowmap       = discret_->DofRowMap();

  // generate a parameterlist for communication and control
  ParameterList filterparams;
  // action for elements
  filterparams.set("action","calc_fluid_box_filter");

  // set state vector to pass distributed vector to the element
  discret_->ClearState();
  discret_->SetState("u and p (trial)",velaf_);

  // loop all nodes on the processor
  for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
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

    // skip this node if it is on a wall
    if (is_no_slip_node == numdim_)
    {
      continue;
    }

    // now check whether we have a pbc condition on this node
    vector<DRT::Condition*> mypbc;

    lnode->GetCondition("SurfacePeriodic",mypbc);

    // check whether a periodic boundary condition is active on this node
    bool ispbcmaster = false;

    if (mypbc.size()>0)
    {
      // get the list of all his slavenodes
      map<int, vector<int> >::iterator master = mapmastertoslave_.find(lnode->Id());

      // slavenodes are ignored
      if(master == mapmastertoslave_.end())
      {
        // the node is a slave --- so don't do anything
        continue;
      }

      // we have a master. Remember this cause we have to extend the patch
      // to the other side...
      ispbcmaster = true;
    }

    // generate a vector of all adjacent elements
    vector <DRT::Element*> patcheles;
    for(int rr=0;rr<lnode->NumElement();++rr)
    {
      patcheles.push_back(lnode->Elements()[rr]);
    }

    // add the elements connected to the slavenodes --- master and
    // slavenodes are treated like they were identical!
    if(ispbcmaster == true)
    {
      for (unsigned slavecount = 0;slavecount<mapmastertoslave_[lnode->Id()].size();++slavecount)
      {
        // get the corresponding slavenodes
        DRT::Node*  slavenode = discret_->gNode(mapmastertoslave_[lnode->Id()][slavecount]);

        // add the elements
        for(int rr=0;rr<slavenode->NumElement();++rr)
        {
          patcheles.push_back(slavenode->Elements()[rr]);
        }
      }
    }

    // define element matrices and vectors --- they are used to
    // transfer information into the element routine and back
    Epetra_SerialDenseMatrix ep_reystress_hat(3,3);
    Epetra_SerialDenseMatrix ep_modeled_stress_grid_scale_hat(3,3);
    Epetra_SerialDenseVector ep_velaf_hat (3);
    Epetra_SerialDenseVector dummy1;
    Epetra_SerialDenseVector dummy2;

    // the patch volume has to be initialised to zero
    double patchvolume = 0;

    // loop all adjacent elements to this node
    for (unsigned nele=0;nele<patcheles.size();++nele)
    {
      // get the neighbouring element
      DRT::Element* nbele = (patcheles[nele]);

      // get element location vector, dirichlet flags and ownerships
      vector<int> lm;
      vector<int> lmowner;
      nbele->LocationVector(*discret_,lm,lmowner);

      // call the element evaluate method to integrate functions
      // against heaviside function element
      int err = nbele->Evaluate(filterparams,
                                *discret_,
                                lm,
                                ep_reystress_hat,
                                ep_modeled_stress_grid_scale_hat,
                                ep_velaf_hat,dummy1,dummy2);
      if (err) dserror("Proc %d: Element %d returned err=%d",
                       discret_->Comm().MyPID(),nbele->Id(),err);

      // get contribution to patch volume of this element. Add it up.
      double volume_contribution =filterparams.get<double>("volume_contribution");

      patchvolume+=volume_contribution;
    }

    // wrap Epetra Object in Blitz array
    blitz::Array<double, 1> velaf_hat(ep_velaf_hat.Values(),
                                      blitz::shape(ep_velaf_hat.Length()),
                                      blitz::neverDeleteData);

    blitz::Array<double, 2> reystress_hat(ep_reystress_hat.A(),
                                          blitz::shape(ep_reystress_hat.M(),ep_reystress_hat.N()),
                                          blitz::neverDeleteData,
                                          blitz::ColumnMajorArray<2>());

    blitz::Array<double, 2> modeled_stress_grid_scale_hat(ep_modeled_stress_grid_scale_hat.A(),
                                                          blitz::shape(ep_modeled_stress_grid_scale_hat.M(),ep_modeled_stress_grid_scale_hat.N()),
                                                          blitz::neverDeleteData,
                                                          blitz::ColumnMajorArray<2>());

    // normalize the computed convolution products by the complete patchvolume
    reystress_hat                /=patchvolume;
    modeled_stress_grid_scale_hat/=patchvolume;
    velaf_hat                    /=patchvolume;

    // now assemble the computed values into the global vector
    double val = 0;
    int    id  = (lnode->Id());

    for (int idim =0;idim<3;++idim)
    {
      val = velaf_hat(idim);
      ((*filtered_vel_)(idim))->ReplaceGlobalValues(1,&val,&id);

      for (int jdim =0;jdim<3;++jdim)
      {
        val = reystress_hat (idim,jdim);
        ((*filtered_reynoldsstress_ )       (3*idim+jdim))->ReplaceGlobalValues(1,&val,&id);

        val = modeled_stress_grid_scale_hat(idim,jdim);
        ((*filtered_modeled_subgrid_stress_)(3*idim+jdim))->ReplaceGlobalValues(1,&val,&id);
      }
    }

    // for masternodes, all slavenodes get the same values
    if (ispbcmaster == true)
    {
      for (unsigned slavecount = 0;slavecount<mapmastertoslave_[lnode->Id()].size();++slavecount)
      {
        // get the corresponding slavenodes
        DRT::Node*  slavenode = discret_->gNode(mapmastertoslave_[lnode->Id()][slavecount]);

        int    slaveid  = (slavenode->Id());

        for (int idim =0;idim<3;++idim)
        {
          val = (velaf_hat(idim));
          ((*filtered_vel_)(idim))->ReplaceGlobalValues(1,&val,&slaveid);

          for (int jdim =0;jdim<3;++jdim)
          {
            val = reystress_hat (idim,jdim);
            ((*filtered_reynoldsstress_ )       (3*idim+jdim))->ReplaceGlobalValues(1,&val,&slaveid);

            val = modeled_stress_grid_scale_hat(idim,jdim);
            ((*filtered_modeled_subgrid_stress_)(3*idim+jdim))->ReplaceGlobalValues(1,&val,&slaveid);
          }
        }
      }
    }
  }

  // clean up
  discret_->ClearState();

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
  if (planecoords_ == null)
  {
    planecoords_ = rcp( new vector<double>((turbulencestatistics_->ReturnNodePlaneCoords()).size()));
  }

  (*planecoords_) = turbulencestatistics_->ReturnNodePlaneCoords();

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
  }

  // now add all the stuff from the different processors

  for (unsigned rr=0;rr<(*planecoords_).size()-1;++rr)
  {
    discret_->Comm().SumAll(&(local_count_for_average[rr]),&(count_for_average[rr]),1);
    discret_->Comm().SumAll(&(local_ele_sum_LijMij[rr]),&((*averaged_LijMij_)[rr]),1);
    discret_->Comm().SumAll(&(local_ele_sum_MijMij[rr]),&((*averaged_MijMij_)[rr]),1);
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
    modelparams->set<RefCountPtr<vector<double> > >("planecoords_",planecoords_);
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
void FluidGenAlphaIntegration::GenAlphaEchoToScreen(
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
        cout << &endl;

        // alpha_M
        cout << "                             alpha_M = ";
        cout << params_.get<double>("alpha_M");
        cout << &endl;

        // gamma
        cout << "                             gamma   = ";
        cout << 0.5
                +
                params_.get<double>("alpha_M")
                -
                params_.get<double>("alpha_F");
        cout << &endl <<&endl;

        // linearisation
        if(params_.get<bool>("Use reaction terms for linearisation",false)
           ==
           true)
        {
          cout << "Linearisation              : ";
          cout << "Including reactive terms (Newton-like)" <<&endl;
        }
        else
        {
          cout << "Linearisation              : ";
          cout << "Without reactive terms (fixed-point-like)" <<&endl;
        }
        cout << &endl;
      }

      //--------------------------------------------------------------------
      /* output of stabilisation details */
      {
        ParameterList *  stabparams=&(params_.sublist("STABILIZATION"));

        // general
        cout << "Stabilization type         : ";
        cout << stabparams->get<string>("STABTYPE");
        cout << &endl;

        // time-dependent subgrid scales
        cout << "                             ";
        cout << stabparams->get<string>("TDS")<< &endl;
        cout << &endl;

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
        cout << &endl;

        // supg stabilization?
        cout <<  "                             ";
        cout << "SUPG            = ";
        cout << stabparams->get<string>("SUPG")           ;
        cout << &endl;

        // pspg stabilization?
        cout <<  "                             ";
        cout << "PSPG            = ";
        cout << stabparams->get<string>("PSPG")           ;
        cout << &endl;

        // viscous stabilization?
        cout <<  "                             ";
        cout << "VSTAB           = ";
        cout << stabparams->get<string>("VSTAB")          ;
        cout << &endl;

        // least-squares continuity stabilization?
        cout <<  "                             ";
        cout << "CSTAB           = ";
        cout << stabparams->get<string>("CSTAB")          ;
        cout << &endl;

        // resvmm turbulence modeling, cross-stress part?
        cout <<  "                             ";
        cout << "CROSS-STRESS    = ";
        cout << stabparams->get<string>("CROSS-STRESS")   ;
        cout << &endl;

        // resvmm turbulence modeling, Reynolds-stress part?
        cout <<  "                             ";
        cout << "REYNOLDS-STRESS = ";
        cout << stabparams->get<string>("REYNOLDS-STRESS");
        cout << &endl;
        cout << &endl;
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
          cout << &endl << &endl;

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
            cout << &endl;

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
              cout << &endl;

              cout << "                             "          ;
              cout << "- viscous length      :   l_tau= "      ;
              cout << modelparams->get<double>("CHANNEL_L_TAU");
              cout << &endl;
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
            cout << &endl;
          }

          if (special_flow_ == "channel_flow_of_height_2")
          {
            cout << "                             " ;
            cout << "Turbulence statistics are evaluated ";
            cout << "for a turbulent channel flow.\n";
            cout << "                             " ;
            cout << "The solution is averaged over the homogeneous ";
            cout << hom_plane;
            cout << " plane and over time.\n";
          }
          cout << &endl;
          cout << &endl;
        }

        //--------------------------------------------------------------------
        /* output of fine-scale subgrid-vicosity approach if any */
        if (fssgv_ != "No")
        {
          cout << "Fine-scale subgrid-viscosity approach based on AVM3: ";
          cout << &endl << &endl;
          cout << params_.get<string>("fs subgrid viscosity");

          if (fssgv_ == "Smagorinsky_all" || fssgv_ == "Smagorinsky_small" ||
              fssgv_ == "mixed_Smagorinsky_all" || fssgv_ == "mixed_Smagorinsky_small")
          {
            cout << " with Smagorinsky constant Cs= ";
            cout << modelparams->get<double>("C_SMAGORINSKY") ;
          }
          cout << &endl << &endl;
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
      cout << &endl;

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
      cout << &endl;

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
    else if("print nonlin iter converged")
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
void FluidGenAlphaIntegration::UpdateGridv()
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

Idea of this routine:

create

map< label, set<DRT::Node*> >

which is a set of nodes to each L&D Id
nodal forces of all the nodes within one set are added to one L&D force

Notice: Angular moments obtained from lift&drag forces currently refere
        to the initial configuration, i.e. are built with the coordinates
        X of a particular node irrespective of its current position.
*/
void FluidGenAlphaIntegration::LiftDrag() const
{
  std::map< const int, std::set<DRT::Node* > > ldnodemap;
  std::map< const int, const std::vector<double>* > ldcoordmap;

  // allocate and initialise LiftDrag conditions
  std::vector<DRT::Condition*> ldconds;
  discret_->GetCondition("LIFTDRAG",ldconds);

  // space dimension of the problem
  const int ndim = params_.get<int>("number of velocity degrees of freedom");

  // there is an L&D condition if it has a size
  if( ldconds.size() )
  {

    // prepare output
    if (myrank_==0)
    {
      cout << "Lift and drag calculation:" << "\n";
      if (ndim == 2)
      {
        cout << "lift'n'drag Id  ";
        cout << "    F_x         ";
        cout << "    F_y         ";
        cout << "    M_z :\n"     ;
      }
      if (ndim == 3)
      {
        cout << "lift'n'drag Id  ";
        cout << "    F_x         ";
        cout << "    F_y         ";
        cout << "    F_z         ";
        cout << "    M_x         ";
        cout << "    M_y         ";
        cout << "    M_z : \n"    ;
      }
    }

    // ---------------------------------------------------------------
    // sort data

    // loop L&D conditions (i.e. lines in .dat file)
    for( unsigned i=0; i<ldconds.size(); ++i)
    {
      /* get label of present LiftDrag condition  */
      const unsigned int label = ldconds[i]->Getint("label");
      /* get new nodeset for new label OR:
         return pointer to nodeset for known label ... */
      std::set<DRT::Node*>& nodes = ldnodemap[label];

      // centre coordinates to present label
      ldcoordmap[label] = ldconds[i]->Get<vector<double> >("centerCoord");

      /* get pointer to its nodal Ids*/
      const vector<int>* ids = ldconds[i]->Get<vector<int> >("Node Ids");

      /* put all nodes belonging to the L&D line or surface into
         'nodes' which are associated with the present label */
      for (unsigned j=0; j<ids->size(); ++j)
      {
        // give me present node Id
        const int node_id = (*ids)[j];
        // put it into nodeset of actual label if node is new and mine
        if( discret_->HaveGlobalNode(node_id) && discret_->gNode(node_id)->Owner()==myrank_ )
	  nodes.insert(discret_->gNode(node_id));
      }
    } // end loop over conditions


    // now step the label map
    for( std::map< const int, std::set<DRT::Node*> >::const_iterator labelit = ldnodemap.begin();
         labelit != ldnodemap.end(); ++labelit )
    {
      // pointer to nodeset of present label
      const std::set<DRT::Node*>& nodes = labelit->second;
      // the present label
      const int label = labelit->first;
      // vector with lift&drag forces
      std::vector<double> values(6,0.0);
      // vector with lift&drag forces after communication
      std::vector<double> resultvec(6,0.0);

      // get also pointer to centre coordinates
      const std::vector<double>* centerCoord = ldcoordmap[label];

      // loop all nodes within my set
      for( std::set<DRT::Node*>::const_iterator actnode = nodes.begin();
           actnode != nodes.end();
           ++actnode)
      {
        // pointer to nodal coordinates
        const double* x = (*actnode)->X();
        const Epetra_BlockMap& rowdofmap = force_->Map();
        const std::vector<int> dof = discret_->Dof(*actnode);

        std::vector<double> distances (3);
        for (unsigned j=0; j<3; ++j)
        {
          distances[j]= x[j]-(*centerCoord)[j];
        }
        // get nodal forces
        const double fx = (*force_)[rowdofmap.LID(dof[0])];
        const double fy = (*force_)[rowdofmap.LID(dof[1])];
        const double fz = (*force_)[rowdofmap.LID(dof[2])];
        values[0] += fx;
        values[1] += fy;
        values[2] += fz;

        // calculate nodal angular momenta
        values[3] += distances[1]*fz-distances[2]*fy;
        values[4] += distances[2]*fx-distances[0]*fz;
        values[5] += distances[0]*fy-distances[1]*fx;
      } // end: loop over nodes

      // care for the fact that we are (most likely) parallel
      force_->Comm().SumAll (&(values[0]), &(resultvec[0]), 6);

      // do the output
      if (myrank_==0)
      {
        if (ndim == 2)
	{
	  cout << "     " << label << "         ";
      cout << std::scientific << resultvec[0] << "    ";
	  cout << std::scientific << resultvec[1] << "    ";
	  cout << std::scientific << resultvec[5];
	  cout << "\n";
        }
        if (ndim == 3)
	{
	  cout << "     " << label << "         ";
      cout << std::scientific << resultvec[0] << "    ";
	  cout << std::scientific << resultvec[1] << "    ";
	  cout << std::scientific << resultvec[2] << "    ";
	  cout << std::scientific << resultvec[3] << "    ";
	  cout << std::scientific << resultvec[4] << "    ";
	  cout << std::scientific << resultvec[5];
	  cout << "\n";
	}
      }
    } // end: loop over L&D labels
    if (myrank_== 0)
    {
      cout << "\n";
    }
  }
}//FluidGenAlphaIntegration::LiftDrag



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FluidGenAlphaIntegration::IntegrateInterfaceShape(
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
  discret_->SetState("dispnp", dispnp_);
  discret_->EvaluateCondition(eleparams,integratedshapefunc,condname);
  discret_->ClearState();

  return integratedshapefunc;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FluidGenAlphaIntegration::LinearRelaxationSolve(
  Teuchos::RCP<Epetra_Vector> relax
  )
{
  dserror("No steepest descent for genalpha fluid\n");

  return;
}


#endif
