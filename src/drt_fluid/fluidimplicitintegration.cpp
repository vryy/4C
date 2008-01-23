/*!----------------------------------------------------------------------
\file fluidimplicitintegration.cpp
\brief Control routine for fluid time integration. Includes

     o Single step one-step-theta time integration

     o Two step BDF2 Gear's methode (with optional one-step-theta
                                     start step)



<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "fluidimplicitintegration.H"
#include "../drt_lib/drt_nodematchingoctree.H"
#include "../drt_lib/drt_periodicbc.H"
#include "../drt_lib/drt_function.H"


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Constructor (public)                                     gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
FluidImplicitTimeInt::FluidImplicitTimeInt(RefCountPtr<DRT::Discretization> actdis,
                                           LINALG::Solver&       solver,
                                           ParameterList&        params,
                                           IO::DiscretizationWriter& output,
                                           bool alefluid) :
  // call constructor for "nontrivial" objects
  discret_(actdis),
  solver_ (solver),
  params_ (params),
  output_ (output),
  alefluid_(alefluid),
  time_(0.0),
  step_(0),
  restartstep_(0),
  uprestart_(params.get("write restart every", -1)),
  writestep_(0),
  upres_(params.get("write solution every", -1)),
  writestresses_(params.get<int>("write stresses", 0))
{

  // -------------------------------------------------------------------
  // get the basic parameters first
  // -------------------------------------------------------------------

  // bound for the number of timesteps
  stepmax_                   =params_.get<int>   ("max number timesteps");
  // max. sim. time
  maxtime_                   =params_.get<double>("total time");

  // parameter for time-integration
  theta_                     =params_.get<double>("theta");
  // which kind of time-integration
  timealgo_                  =params_.get<FLUID_TIMEINTTYPE>("time int algo");

  // parameter for linearisation scheme (fixed point like or newton like)
  newton_ = params_.get<bool>("Use reaction terms for linearisation",false);

  dtp_ = dta_ = params_.get<double>("time step size");

  // (fine-scale) subgrid viscosity?
  fssgv_ = params_.get<int>("fs subgrid viscosity",0);

  // Smagorinsky model parameter from turbulence model sublist
  ParameterList *  turbmodelparams =&(params_.sublist("TURBULENCE MODEL"));
  Cs_fs_ = turbmodelparams->get<double>("C_SMAGORINSKY",0.0);

  // -------------------------------------------------------------------
  // connect degrees of freedom for periodic boundary conditions
  // -------------------------------------------------------------------
  PeriodicBoundaryConditions::PeriodicBoundaryConditions pbc(discret_);
  pbc.UpdateDofsForPeriodicBoundaryConditions();

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
  //
  // The maps are designed assuming that every node has pressure and
  // velocity degrees of freedom --- this won't work for inf-sup stable
  // elements at the moment!
  // -------------------------------------------------------------------
  {
    const int numdim = params_.get<int>("number of velocity degrees of freedom");

    // Allocate integer vectors which will hold the dof number of the
    // velocity or pressure dofs
    vector<int> velmapdata;
    vector<int> premapdata;

    velmapdata.reserve(discret_->NumMyRowNodes()*numdim);
    premapdata.reserve(discret_->NumMyRowNodes());

    int countveldofs = 0;
    int countpredofs = 0;

    for (int i=0; i<discret_->NumMyRowNodes(); ++i)
    {
      DRT::Node* node = discret_->lRowNode(i);
      vector<int> dof = discret_->Dof(node);

      int numdofs = dof.size();
      if (numdofs==numdim+1)
      {
        for (int j=0; j<numdofs-1; ++j)
        {
	  // add this velocity dof to the velmapdata vector
	  velmapdata.push_back(dof[j]);
          countveldofs+=1;
        }

        // add this pressure dof to the premapdata vector
        premapdata.push_back(dof[numdofs-1]);
        countpredofs += 1;
      }
      else
      {
        dserror("up to now fluid expects numdim vel + one pre dofs");
      }
    }

    // the rowmaps are generated according to the pattern provided by
    // the data vectors
    velrowmap_ = rcp(new Epetra_Map(-1,
                                    velmapdata.size(),&velmapdata[0],0,
                                    discret_->Comm()));
    prerowmap_ = rcp(new Epetra_Map(-1,
                                    premapdata.size(),&premapdata[0],0,
                                    discret_->Comm()));
  }
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
  maxentriesperrow_ = 108;

  // initialize standard (stabilized) system matrix
  sysmat_ = null;

  // -------------------------------------------------------------------
  // create empty vectors
  // -------------------------------------------------------------------

  // Vectors passed to the element
  // -----------------------------

  // accelerations at time n and n-1
  accn_         = LINALG::CreateVector(*dofrowmap,true);
  accnm_        = LINALG::CreateVector(*dofrowmap,true);

  // velocities and pressures at time n+1, n and n-1
  velnp_        = LINALG::CreateVector(*dofrowmap,true);
  veln_         = LINALG::CreateVector(*dofrowmap,true);
  velnm_        = LINALG::CreateVector(*dofrowmap,true);

  if (alefluid_)
  {
    dispnp_       = LINALG::CreateVector(*dofrowmap,true);
    dispn_        = LINALG::CreateVector(*dofrowmap,true);
    gridv_        = LINALG::CreateVector(*dofrowmap,true);
  }

  // histvector --- a linear combination of velnm, veln (BDF)
  //                or veln, accn (One-Step-Theta)
  hist_         = LINALG::CreateVector(*dofrowmap,true);

  // Vectors associated to boundary conditions
  // -----------------------------------------

  // toggle vector indicating which dofs have Dirichlet BCs
  dirichtoggle_ = LINALG::CreateVector(*dofrowmap,true);
  // opposite of dirichtoggle vector, ie for each component
  invtoggle_    = LINALG::CreateVector(*dofrowmap,false);

  // a vector of zeros to be used to enforce zero dirichlet boundary conditions
  zeros_        = LINALG::CreateVector(*dofrowmap,true);

  // the vector containing body and surface forces
  neumann_loads_= LINALG::CreateVector(*dofrowmap,true);

  // Vectors used for solution process
  // ---------------------------------

  // rhs: standard (stabilized) residual vector (rhs for the incremental form)
  residual_     = LINALG::CreateVector(*dofrowmap,true);
  trueresidual_ = LINALG::CreateVector(*dofrowmap,true);

  // Nonlinear iteration increment vector
  incvel_       = LINALG::CreateVector(*dofrowmap,true);

  // necessary only for the VM3 approach
  if (fssgv_ > 0)
  {
    // initialize (fine-scale) subgrid-viscosity system matrix
    sysmat_sv_ = null;

    // residual vector containing (fine-scale) subgrid-viscosity residual
    residual_sv_  = LINALG::CreateVector(*dofrowmap,true);
  }

  // -------------------------------------------------------------------
  // create timers and time monitor
  // -------------------------------------------------------------------
  timedyntot_     = TimeMonitor::getNewTimer("dynamic routine total"     );
  timedyninit_    = TimeMonitor::getNewTimer(" + initial phase"          );
  timedynloop_    = TimeMonitor::getNewTimer(" + time loop"              );
  timenlnloop_    = TimeMonitor::getNewTimer("   + nonlinear iteration"  );
  timeeleloop_    = TimeMonitor::getNewTimer("      + element calls"     );
  timeapplydirich_= TimeMonitor::getNewTimer("      + apply dirich cond.");
  timesolver_     = TimeMonitor::getNewTimer("      + solver calls"      );

  return;

} // FluidImplicitTimeInt::FluidImplicitTimeInt


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Start the time integration. Allows                                   |
 |                                                                      |
 |  o starting steps with different algorithms                          |
 |  o the "standard" time integration                                   |
 |                                                           gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FluidImplicitTimeInt::Integrate()
{
  // time measurement --- start TimeMonitor tm0 and tm1
  tm0_ref_        = rcp(new TimeMonitor(*timedyntot_ ));

  // bound for the number of startsteps
  int    numstasteps         =params_.get<int>   ("number of start steps");

  if (timealgo_==timeint_stationary)
    // stationary case
    this->SolveStationaryProblem();

  else  // instationary case
  {
    // start procedure
    if (step_<numstasteps)
    {
      if (numstasteps>stepmax_)
      {
        dserror("more startsteps than steps");
      }

      dserror("no starting steps supported");
    }

    // continue with the final time integration
    this->TimeLoop();
  }

  // print the results of time measurements

  tm0_ref_ = null; // end total time measurement
  //cout<<endl<<endl;
  TimeMonitor::summarize();

  return;
} // FluidImplicitTimeInt::Integrate



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | contains the time loop                                    gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FluidImplicitTimeInt::TimeLoop()
{
  // start time measurement for timeloop
  tm2_ref_ = rcp(new TimeMonitor(*timedynloop_));

  while (step_<stepmax_ and time_<maxtime_)
  {
    PrepareTimeStep();

    // -------------------------------------------------------------------
    //                     solve nonlinear equation
    // -------------------------------------------------------------------
    this->NonlinearSolve();

    // -------------------------------------------------------------------
    //                         update solution
    //        current solution becomes old solution of next timestep
    //
    // One-step-Theta: (step>1)
    //
    //  accn_  = (velnp_-veln_) / (Theta * dt) - (1/Theta -1) * accn_
    //  "(n+1)"
    //
    //  velnm_ =veln_
    //  veln_  =velnp_
    //
    // BDF2:           (step>1)
    //
    //               2*dt(n)+dt(n-1)		  dt(n)+dt(n-1)
    //  accn_   = --------------------- velnp_ - --------------- veln_
    //             dt(n)*[dt(n)+dt(n-1)]	  dt(n)*dt(n-1)
    //
    //                     dt(n)
    //           + ----------------------- velnm_
    //             dt(n-1)*[dt(n)+dt(n-1)]
    //
    //
    //  velnm_ =veln_
    //  veln_  =velnp_
    //
    // BDF2 and  One-step-Theta: (step==1)
    //
    // The given formulas are only valid from the second timestep. In the
    // first step, the acceleration is calculated simply by
    //
    //  accn_  = (velnp_-veln_) / (dt)
    //
    // -------------------------------------------------------------------
    this->TimeUpdate();

    // -------------------------------------------------------------------
    // evaluate error for test flows with analytical solutions
    // -------------------------------------------------------------------
    this->EvaluateErrorComparedToAnalyticalSol();

    // -------------------------------------------------------------------
    //                         output of solution
    // -------------------------------------------------------------------
    this->Output();

    // -------------------------------------------------------------------
    //                    calculate lift'n'drag forces
    // -------------------------------------------------------------------
    this->LiftDrag();

    // -------------------------------------------------------------------
    //                       update time step sizes
    // -------------------------------------------------------------------
    dtp_ = dta_;

    // -------------------------------------------------------------------
    //                    stop criterium for timeloop
    // -------------------------------------------------------------------
  }

  // end time measurement for timeloop
  tm2_ref_ = null;

  return;
} // FluidImplicitTimeInt::TimeIntegrateFromTo


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | set part of the residual vector belonging to the old timestep        |
 |                                                           gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FluidImplicitTimeInt::SetOldPartOfRighthandside()
{
  /*

  One-step-Theta:

                 hist_ = veln_ + dt*(1-Theta)*accn_


  BDF2: for constant time step:

                 hist_ = 4/3 veln_ - 1/3 velnm_

  */
  switch (timealgo_)
  {
  case timeint_one_step_theta: /* One step Theta time integration */
    hist_->Update(1.0, *veln_, dta_*(1.0-theta_), *accn_, 0.0);
    break;

  case timeint_bdf2:	/* 2nd order backward differencing BDF2	*/
    hist_->Update(4./3., *veln_, -1./3., *velnm_, 0.0);
    break;

  default:
    dserror("Time integration scheme unknown!");
  }
  return;
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  do explicit predictor step to start nonlinear iteration from a      |
 |  better value                                             gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FluidImplicitTimeInt::ExplicitPredictor()
{
  const double fact1 = dta_*(1.0+dta_/dtp_);
  const double fact2 = DSQR(dta_/dtp_);

  velnp_->Update( fact1,*accn_ ,1.0);
  velnp_->Update(-fact2,*veln_ ,1.0);
  velnp_->Update( fact2,*velnm_,1.0);

  return;
} // FluidImplicitTimeInt::ExplicitPredictor


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | setup the variables to do a new time step                 u.kue 06/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FluidImplicitTimeInt::PrepareTimeStep()
{
  // -------------------------------------------------------------------
  //              set time dependent parameters
  // -------------------------------------------------------------------
  step_ += 1;

  time_ += dta_;

  // for bdf2 theta is set  by the timestepsizes, 2/3 for const. dt
  if (timealgo_==timeint_bdf2)
  {
    theta_ = (dta_+dtp_)/(2.0*dta_ + dtp_);
  }

  // -------------------------------------------------------------------
  //                         out to screen
  // -------------------------------------------------------------------
  if (myrank_==0)
  {
    switch (timealgo_)
    {
    case timeint_one_step_theta:
      printf("TIME: %11.4E/%11.4E  DT = %11.4E  One-Step-Theta  STEP = %4d/%4d \n",
             time_,maxtime_,dta_,step_,stepmax_);
      break;
    case timeint_bdf2:
      printf("TIME: %11.4E/%11.4E  DT = %11.4E     BDF2         STEP = %4d/%4d \n",
             time_,maxtime_,dta_,step_,stepmax_);
      break;
    default:
      dserror("parameter out of range: IOP\n");
    } /* end of switch(timealgo) */
  }

  // -------------------------------------------------------------------
  // set part of the rhs vector belonging to the old timestep
  //
  //
  //         One-step-Theta:
  //
  //                 hist_ = veln_ + dta*(1-Theta)*accn_
  //
  //
  //         BDF2: for constant time step:
  //
  //                   hist_ = 4/3 veln_ - 1/3 velnm_
  //
  // -------------------------------------------------------------------
  this->SetOldPartOfRighthandside();

  // -------------------------------------------------------------------
  //                     do explicit predictor step
  //
  //                     +-                                      -+
  //                     | /     dta \          dta  veln_-velnm_ |
  // velnp_ =veln_ + dta | | 1 + --- | accn_ - ----- ------------ |
  //                     | \     dtp /          dtp     dtp       |
  //                     +-                                      -+
  //
  // -------------------------------------------------------------------
  if (step_>1)
  {
    this->ExplicitPredictor();
  }

  // -------------------------------------------------------------------
  //         evaluate dirichlet and neumann boundary conditions
  // -------------------------------------------------------------------
  {
    ParameterList eleparams;
    // action for elements
    eleparams.set("action","calc_fluid_eleload");
    // choose what to assemble
    eleparams.set("assemble matrix 1",false);
    eleparams.set("assemble matrix 2",false);
    eleparams.set("assemble vector 1",true);
    eleparams.set("assemble vector 2",false);
    eleparams.set("assemble vector 3",false);
    // other parameters needed by the elements
    eleparams.set("total time",time_);
    eleparams.set("delta time",dta_);
    eleparams.set("thsl",theta_*dta_);

    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("velnp",velnp_);
    // predicted dirichlet values
    // velnp then also holds prescribed new dirichlet values
    // dirichtoggle is 1 for dirichlet dofs, 0 elsewhere
    discret_->EvaluateDirichlet(eleparams,*velnp_,*dirichtoggle_);
    discret_->ClearState();

    // evaluate Neumann conditions
    eleparams.set("total time",time_);
    eleparams.set("thsl",theta_*dta_);

    neumann_loads_->PutScalar(0.0);
    discret_->EvaluateNeumann(eleparams,*neumann_loads_);
    discret_->ClearState();
  }

  //----------------------- compute an inverse of the dirichtoggle vector
  invtoggle_->PutScalar(1.0);
  invtoggle_->Update(-1.0,*dirichtoggle_,1.0);
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | contains the nonlinear iteration loop                     gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FluidImplicitTimeInt::NonlinearSolve()
{
  // start time measurement for nonlinear iteration
  tm6_ref_ = rcp(new TimeMonitor(*timenlnloop_));

  const Epetra_Map* dofrowmap       = discret_->DofRowMap();

  // ---------------------------------------------- nonlinear iteration
  // maximum number of nonlinear iteration steps
  const int     itemax    =params_.get<int>   ("max nonlin iter steps");

  // ------------------------------- stop nonlinear iteration when both
  //                                 increment-norms are below this bound
  const double  ittol     =params_.get<double>("tolerance for nonlin iter");

  int               itnum         = 0;
  bool              stopnonliniter = false;

  double            dtsolve = 0;
  double            dtele   = 0;
  double            tcpu   ;

  if (myrank_ == 0)
  {
    printf("+------------+-------------------+--------------+--------------+--------------+--------------+\n");
    printf("|- step/max -|- tol      [norm] -|-- vel-res ---|-- pre-res ---|-- vel-inc ---|-- pre-inc ---|\n");
  }

  while (stopnonliniter==false)
  {
    itnum++;

    double density = 1;

    // -------------------------------------------------------------------
    // call elements to calculate system matrix
    // -------------------------------------------------------------------
    {
      // start time measurement for element call
      tm3_ref_ = rcp(new TimeMonitor(*timeeleloop_));
      // get cpu time
      tcpu=ds_cputime();

      sysmat_ = LINALG::CreateMatrix(*dofrowmap,maxentriesperrow_);

      // add Neumann loads
      residual_->Update(1.0,*neumann_loads_,0.0);

      // create the parameters for the discretization
      ParameterList eleparams;

      // action for elements
      if (timealgo_==timeint_stationary)
          eleparams.set("action","calc_fluid_stationary_systemmat_and_residual");
      else
          eleparams.set("action","calc_fluid_systemmat_and_residual");

      // other parameters that might be needed by the elements
      eleparams.set("total time",time_);
      eleparams.set("thsl",theta_*dta_);
      eleparams.set("fs subgrid viscosity",fssgv_);
      eleparams.set("fs Smagorinsky parameter",Cs_fs_);
      eleparams.set("include reactive terms for linearisation",newton_);

      // set vector values needed by elements
      discret_->ClearState();
      discret_->SetState("velnp",velnp_);
      discret_->SetState("hist"  ,hist_ );
      if (alefluid_)
      {
        discret_->SetState("dispnp", dispnp_);
        discret_->SetState("gridv", gridv_);
      }

      // decide whether VM3-based solution approach or standard approach
      if (fssgv_ > 0)
      {
        // extract the ML parameters
        ParameterList&  mllist = solver_.Params().sublist("ML Parameters");

        // subgrid-viscosity-scaling vector
        sugrvisc_     = LINALG::CreateVector(*dofrowmap,true);

        // define flag for computation of matrices (true only in first time step)
        bool compute;

        if (step_ == 1) 
        {
          compute=true;

          // create scale-separation matrix
          scalesep_ = LINALG::CreateMatrix(*dofrowmap,maxentriesperrow_);

          // create (fine-scale) subgrid-viscosity matrix
          sysmat_sv_ = LINALG::CreateMatrix(*dofrowmap,maxentriesperrow_); 

          // call loop over elements (two matrices + subgr.-visc.-scal. vector)
          discret_->Evaluate(eleparams,sysmat_,sysmat_sv_,residual_,sugrvisc_);
          discret_->ClearState();

          // finalize the (fine-scale) subgrid-viscosity matrix
          LINALG::Complete(*sysmat_sv_);

          // apply DBC to (fine-scale) subgrid-viscosity matrix
          LINALG::ApplyDirichlettoSystem(sysmat_sv_,incvel_,residual_sv_,zeros_,dirichtoggle_);

          // extract the ML parameters
          ParameterList&  mllist = solver_.Params().sublist("ML Parameters");

          // call the VM3 constructor
          RCP<VM3_Solver> vm3_solver = rcp(new VM3_Solver::VM3_Solver(scalesep_,sysmat_sv_,sysmat_,sugrvisc_,residual_sv_,residual_,velnp_,dirichtoggle_,mllist,compute) );

          // call the VM3 scale separator for incremental formulation:
          // precomputation of unscaled S^T*M*S
          vm3_solver-> VM3_Solver::Separate(scalesep_,sysmat_sv_);
        }
        else
        {
          compute=false;

          // call loop over elements (one matrix + subgr.-visc.-scal. vector)
          discret_->Evaluate(eleparams,sysmat_,residual_,sugrvisc_);
          discret_->ClearState();
        }
        // call the VM3 constructor
        RCP<VM3_Solver> vm3_solver = rcp(new VM3_Solver::VM3_Solver(scalesep_,sysmat_sv_,sysmat_,sugrvisc_,residual_sv_,residual_,velnp_,dirichtoggle_,mllist,compute) );

        residual_sv_->PutScalar(0.0);
        // call the VM3 scaling:
        // scale precomputed matrix product by subgrid-viscosity-scaling vector
        vm3_solver->VM3_Solver::Scale(sysmat_sv_,sysmat_,sugrvisc_,residual_sv_,residual_,velnp_,true);
      }
      else
      {
        // call standard loop over elements
        discret_->Evaluate(eleparams,sysmat_,residual_);
        discret_->ClearState();
      }

      density = eleparams.get("density", 0.0);

      // finalize the complete matrix
      LINALG::Complete(*sysmat_);
      maxentriesperrow_ = sysmat_->MaxNumEntries();

      // end time measurement for element call
      tm3_ref_=null;
      dtele=ds_cputime()-tcpu;
    }

    // How to extract the density from the fluid material?
    trueresidual_->Update(density/dta_/theta_,*residual_,0.0);

    // blank residual DOFs which are on Dirichlet BC
    // We can do this because the values at the dirichlet positions
    // are not used anyway.
    // We could avoid this though, if velrowmap_ and prerowmap_ would
    // not include the dirichlet values as well. But it is expensive
    // to avoid that.
    {
      Epetra_Vector residual(*residual_);
      residual_->Multiply(1.0,*invtoggle_,residual,0.0);
    }

    double incvelnorm_L2;
    double incprenorm_L2;

    double velnorm_L2;
    double prenorm_L2;

    double vresnorm;
    double presnorm;

    {
      Epetra_Vector onlyvel(*velrowmap_);
      Epetra_Import importer(*velrowmap_,residual_->Map());

      int err = onlyvel.Import(*residual_,importer,Insert);
      if (err) dserror("Import using importer returned err=%d",err);
      onlyvel.Norm2(&vresnorm);

      err = onlyvel.Import(*incvel_,importer,Insert);
      if (err) dserror("Import using importer returned err=%d",err);
      onlyvel.Norm2(&incvelnorm_L2);

      err = onlyvel.Import(*velnp_,importer,Insert);
      if (err) dserror("Import using importer returned err=%d",err);
      onlyvel.Norm2(&velnorm_L2);
    }
    {
      Epetra_Vector onlypre(*prerowmap_);
      Epetra_Import importer(*prerowmap_,residual_->Map());

      int err = onlypre.Import(*residual_,importer,Insert);
      if (err) dserror("Import using importer returned err=%d",err);
      onlypre.Norm2(&presnorm);

      err = onlypre.Import(*incvel_,importer,Insert);
      if (err) dserror("Import using importer returned err=%d",err);
      onlypre.Norm2(&incprenorm_L2);

      err = onlypre.Import(*velnp_,importer,Insert);
      if (err) dserror("Import using importer returned err=%d",err);
      onlypre.Norm2(&prenorm_L2);
    }

    // care for the case that nothing really happens in the velocity
    // or pressure field
    if (velnorm_L2 < 1e-5)
    {
      velnorm_L2 = 1.0;
    }
    if (prenorm_L2 < 1e-5)
    {
      prenorm_L2 = 1.0;
    }

    //-------------------------------------------------- output to screen
    /* special case of very first iteration step:
        - solution increment is not yet available
        - convergence check is not required (we solve at least once!)    */
    if (itnum == 1)
    {
      if (myrank_ == 0)
      {
        printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   |      --      |      --      |",
               itnum,itemax,ittol,vresnorm,presnorm);
        printf(" (      --     ,te=%10.3E)\n",dtele);
      }
    }
    /* ordinary case later iteration steps:
        - solution increment can be printed
        - convergence check should be done*/
    else
    {
    // this is the convergence check
    // We always require at least one solve. Otherwise the
    // perturbation at the FSI interface might get by unnoticed.
      if (vresnorm <= ittol and presnorm <= ittol and
          incvelnorm_L2/velnorm_L2 <= ittol and incprenorm_L2/prenorm_L2 <= ittol)
      {
        stopnonliniter=true;
        if (myrank_ == 0)
        {
          printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   | %10.3E   | %10.3E   |",
                 itnum,itemax,ittol,vresnorm,presnorm,
                 incvelnorm_L2/velnorm_L2,incprenorm_L2/prenorm_L2);
          printf(" (ts=%10.3E,te=%10.3E)\n",dtsolve,dtele);
          printf("+------------+-------------------+--------------+--------------+--------------+--------------+\n");
        }
        break;
      }
      else // if not yet converged
        if (myrank_ == 0)
        {
          printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   | %10.3E   | %10.3E   |",
                 itnum,itemax,ittol,vresnorm,presnorm,
                 incvelnorm_L2/velnorm_L2,incprenorm_L2/prenorm_L2);
          printf(" (ts=%10.3E,te=%10.3E)\n",dtsolve,dtele);
        }
    }

    // warn if itemax is reached without convergence, but proceed to
    // next timestep...
    if ((itnum == itemax) and (vresnorm > ittol or presnorm > ittol or
                             incvelnorm_L2/velnorm_L2 > ittol or
                             incprenorm_L2/prenorm_L2 > ittol))
    {
      stopnonliniter=true;
      if (myrank_ == 0)
      {
        printf("+---------------------------------------------------------------+\n");
        printf("|            >>>>>> not converged in itemax steps!              |\n");
        printf("+---------------------------------------------------------------+\n");
      }
      break;
    }

    //--------- Apply dirichlet boundary conditions to system of equations
    //          residual discplacements are supposed to be zero at
    //          boundary conditions
    incvel_->PutScalar(0.0);
    {
      // start time measurement for application of dirichlet conditions
      tm4_ref_ = rcp(new TimeMonitor(*timeapplydirich_));

      LINALG::ApplyDirichlettoSystem(sysmat_,incvel_,residual_,zeros_,dirichtoggle_);

      // end time measurement for application of dirichlet conditions
      tm4_ref_=null;
    }

    //-------solve for residual displacements to correct incremental displacements
    {
      // start time measurement for solver call
      tm5_ref_ = rcp(new TimeMonitor(*timesolver_));
      // get cpu time
      tcpu=ds_cputime();

      solver_.Solve(sysmat_,incvel_,residual_,true,itnum==1);

      // end time measurement for application of dirichlet conditions
      tm5_ref_=null;
      dtsolve=ds_cputime()-tcpu;
    }

    //------------------------------------------------ update (u,p) trial
    velnp_->Update(1.0,*incvel_,1.0);

    //------------------------------------------------- check convergence
    //this->NonlinearConvCheck(stopnonliniter,itnum,dtele,dtsolve);
  }

  // end time measurement for nonlinear iteration
  tm6_ref_ = null;

  return;
} // FluidImplicitTimeInt::NonlinearSolve


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | build linear system matrix and rhs                        u.kue 11/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FluidImplicitTimeInt::Evaluate(Teuchos::RCP<const Epetra_Vector> vel)
{
  const Epetra_Map* dofrowmap = discret_->DofRowMap();
  sysmat_ = LINALG::CreateMatrix(*dofrowmap,maxentriesperrow_);

  // set the new solution we just got
  if (vel!=Teuchos::null)
  {
    incvel_->Update(1.0, *vel, 0.);

    //------------------------------------------------ update (u,p) trial
    velnp_->Update(1.0,*incvel_,1.0);
  }

  // add Neumann loads
  residual_->Update(1.0,*neumann_loads_,0.0);

  // create the parameters for the discretization
  ParameterList eleparams;

  // action for elements
  if (timealgo_==timeint_stationary)
    eleparams.set("action","calc_fluid_stationary_systemmat_and_residual");
  else
    eleparams.set("action","calc_fluid_systemmat_and_residual");

  // other parameters that might be needed by the elements
  eleparams.set("total time",time_);
  eleparams.set("thsl",theta_*dta_);
  eleparams.set("include reactive terms for linearisation",newton_);

  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("velnp",velnp_);
  discret_->SetState("hist"  ,hist_ );
  if (alefluid_)
  {
    discret_->SetState("dispnp", dispnp_);
    discret_->SetState("gridv", gridv_);
  }

  // call loop over elements
  discret_->Evaluate(eleparams,sysmat_,residual_);
  discret_->ClearState();

  density_ = eleparams.get("density", 0.0);

  // finalize the system matrix
  LINALG::Complete(*sysmat_);
  maxentriesperrow_ = sysmat_->MaxNumEntries();

  trueresidual_->Update(density_/dta_/theta_,*residual_,0.0);

  // Apply dirichlet boundary conditions to system of equations
  // residual displacements are supposed to be zero at boundary
  // conditions
  incvel_->PutScalar(0.0);
  LINALG::ApplyDirichlettoSystem(sysmat_,incvel_,residual_,zeros_,dirichtoggle_);
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | convergence check for nonlinear iteration                 gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FluidImplicitTimeInt::NonlinearConvCheck(
  bool&  stopnonliniter,
  int    itnum         ,
  double dtele         ,
  double dtsolve
  )
{

  RefCountPtr<Epetra_Vector> onlyvel_ = LINALG::CreateVector(*velrowmap_,true);
  RefCountPtr<Epetra_Vector> onlypre_ = LINALG::CreateVector(*prerowmap_,true);

  // ---------------------------------------------- nonlinear iteration
  // maximum number of nonlinear iteration steps
  const int     itemax    =params_.get<int>   ("max nonlin iter steps");

  // ------------------------------- stop nonlinear iteration when both
  //                                 increment-norms are below this bound
  const double  ittol     =params_.get<double>("tolerance for nonlin iter");


  // extract velocity and pressure increments from increment vector
  LINALG::Export(*incvel_,*onlyvel_);
  LINALG::Export(*incvel_,*onlypre_);
  // calculate L2_Norm of increments
  double incvelnorm_L2;
  double incprenorm_L2;
  onlyvel_->Norm2(&incvelnorm_L2);
  onlypre_->Norm2(&incprenorm_L2);

  // extract velocity and pressure solutions from solution vector
  LINALG::Export(*velnp_,*onlyvel_);
  LINALG::Export(*velnp_,*onlypre_);
  // calculate L2_Norm of solution
  double velnorm_L2;
  double prenorm_L2;
  onlyvel_->Norm2(&velnorm_L2);
  onlypre_->Norm2(&prenorm_L2);

  if (velnorm_L2<EPS5)
  {
    velnorm_L2 = 1.0;
  }
  if (prenorm_L2<EPS5)
  {
    prenorm_L2 = 1.0;
  }

  // out to screen
  if(myrank_ == 0)
  {
    printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   |",
           itnum,itemax,ittol,incvelnorm_L2/velnorm_L2,incprenorm_L2/prenorm_L2);
    printf(" (te=%10.3E,ts=%10.3E)\n",dtele,dtsolve);
  }

  // this is the convergence check
  if(incvelnorm_L2/velnorm_L2 <= ittol && incprenorm_L2/prenorm_L2 <= ittol)
  {
    stopnonliniter=true;
    if(myrank_ == 0)
    {
      printf("+------------+-------------------+--------------+--------------+ \n");
    }
  }

  // warn if itemax is reached without convergence, but proceed to
  // next timestep...
  if (itnum == itemax
      &&
      (incvelnorm_L2/velnorm_L2 > ittol
       ||
       incprenorm_L2/prenorm_L2 > ittol))
  {
    stopnonliniter=true;
    if(myrank_ == 0)
    {
      printf("+---------------------------------------------------------------+\n");
      printf("|            >>>>>> not converged in itemax steps!              |\n");
      printf("+---------------------------------------------------------------+\n");
    }
  }

}// FluidImplicitTimeInt::NonlinearConvCheck

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                           gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FluidImplicitTimeInt::TimeUpdate()
{

  // update acceleration
  if (step_ == 1)
  {
    accnm_->PutScalar(0.0);

    // do just a linear interpolation within the first timestep
    accn_->Update( 1.0/dta_,*velnp_,1.0);

    accn_->Update(-1.0/dta_,*veln_ ,1.0);

    // ???
    accnm_->Update(1.0,*accn_,0.0);

  }
  else
  {
    // prev. acceleration becomes (n-1)-accel. of next time step
    accnm_->Update(1.0,*accn_,0.0);

    /*

    One-step-Theta:

    acc(n+1) = (vel(n+1)-vel(n)) / (Theta * dt(n)) - (1/Theta -1) * acc(n)


    BDF2:

                   2*dt(n)+dt(n-1)		    dt(n)+dt(n-1)
      acc(n+1) = --------------------- vel(n+1) - --------------- vel(n)
                 dt(n)*[dt(n)+dt(n-1)]	            dt(n)*dt(n-1)

                         dt(n)
               + ----------------------- vel(n-1)
                 dt(n-1)*[dt(n)+dt(n-1)]

      */

      switch (timealgo_)
      {
          case timeint_one_step_theta: /* One step Theta time integration */
          {
            const double fact1 = 1.0/(theta_*dta_);
            const double fact2 =-1.0/theta_ +1.0;	/* = -1/Theta + 1 */

            accn_->Update( fact1,*velnp_,0.0);
            accn_->Update(-fact1,*veln_ ,1.0);
            accn_->Update( fact2,*accnm_ ,1.0);

            break;
          }
          case timeint_bdf2:	/* 2nd order backward differencing BDF2	*/
          {
            if (dta_*dtp_ < EPS15)
              dserror("Zero time step size!!!!!");
            const double sum = dta_ + dtp_;

            accn_->Update((2.0*dta_+dtp_)/(dta_*sum),*velnp_,
                          - sum /(dta_*dtp_),*veln_ ,0.0);
            accn_->Update(dta_/(dtp_*sum),*velnm_,1.0);
          }
          break;
          default:
            dserror("Time integration scheme unknown for mass rhs!");
      }
    }

  // solution of this step becomes most recent solution of the last step
  velnm_->Update(1.0,*veln_ ,0.0);
  veln_ ->Update(1.0,*velnp_,0.0);

  if (alefluid_)
    dispn_->Update(1.0,*dispnp_,0.0);

  return;
}// FluidImplicitTimeInt::TimeUpdate

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | output of solution vector to binio                        gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FluidImplicitTimeInt::Output()
{

  //-------------------------------------------- output of solution

  //increase counters
    restartstep_ += 1;
    writestep_ += 1;

    #if 0
    // write domain decomposition for visualization
    const Epetra_Map* elerowmap = discret_->ElementRowMap(); 
    RCP<Epetra_Vector> domain_decomp = LINALG::CreateVector(*elerowmap,true);
    for (int lid=0;lid<elerowmap->NumMyElements();++lid)
    {
    	double eleowner=(double) (discret_->lRowElement(lid)->Owner());
    	domain_decomp->ReplaceMyValues(1, &eleowner, &lid);
    	// global element ids
    	//double gid=(double) (discret_->lRowElement(i)->Id());
    	//myelevalues->ReplaceMyValues(1, &gid, &lid);
    }
    #endif

  if (writestep_ == upres_)  //write solution
    {
      writestep_= 0;

      output_.NewStep    (step_,time_);
      output_.WriteVector("velnp", velnp_);
      //output_.WriteVector("residual", trueresidual_);
      //output_.WriteVector("domain_decomp",domain_decomp);
      if (alefluid_)
        output_.WriteVector("dispnp", dispnp_);

      //only perform stress calculation when output is needed
      if (writestresses_)
      {
       RefCountPtr<Epetra_Vector> traction = CalcStresses();
       output_.WriteVector("traction",traction);
      }

      if (restartstep_ == uprestart_) //add restart data
      {
        restartstep_ = 0;

        output_.WriteVector("accn", accn_);
    	output_.WriteVector("veln", veln_);
    	output_.WriteVector("velnm", velnm_);
      }
    }

  // write restart also when uprestart_ is not a integer multiple of upres_
  if ((restartstep_ == uprestart_) && (writestep_ > 0))
  {
    restartstep_ = 0;

    output_.NewStep    (step_,time_);
    output_.WriteVector("velnp", velnp_);
    //output_.WriteVector("residual", trueresidual_);
    //output_.WriteVector("domain_decomp",domain_decomp);
    if (alefluid_)
      output_.WriteVector("dispnp", dispnp_);

    //only perform stress calculation when output is needed
    if (writestresses_)
    {
     RefCountPtr<Epetra_Vector> traction = CalcStresses();
     output_.WriteVector("traction",traction);
    }

    output_.WriteVector("accn", accn_);
    output_.WriteVector("veln", veln_);
    output_.WriteVector("velnm", velnm_);
  }


#if 0  // DEBUG IO --- the whole systemmatrix
      {
	  int rr;
	  int mm;
	  for(rr=0;rr<residual_->MyLength();rr++)
	  {
	      int NumEntries;

	      vector<double> Values(maxentriesperrow);
	      vector<int> Indices(maxentriesperrow);

	      sysmat_->ExtractGlobalRowCopy(rr,
					    maxentriesperrow,
					    NumEntries,
					    &Values[0],&Indices[0]);
	      printf("Row %4d\n",rr);

	      for(mm=0;mm<NumEntries;mm++)
	      {
		  printf("mat [%4d] [%4d] %26.19e\n",rr,Indices[mm],Values[mm]);
	      }
	  }
      }
#endif

#if 0  // DEBUG IO  --- rhs of linear system
      {
	  int rr;
	  double* data = residual_->Values();
	  for(rr=0;rr<residual_->MyLength();rr++)
	  {
	      printf("global %22.15e\n",data[rr]);
	  }
      }

#endif

#if 0  // DEBUG IO --- neumann_loads
      {
	  int rr;
	  double* data = neumann_loads_->Values();
	  for(rr=0;rr<incvel_->MyLength();rr++)
	  {
	      printf("neum[%4d] %26.19e\n",rr,data[rr]);
	  }
      }

#endif


#if 0  // DEBUG IO --- incremental solution
      if (myrank_==0)
      {
        int rr;
        double* data = incvel_->Values();
        for(rr=0;rr<incvel_->MyLength();rr++)
        {
          printf("sol[%4d] %26.19e\n",rr,data[rr]);
        }
      }

#endif


#if 0  // DEBUG IO --- the solution vector after convergence
      if (myrank_==0)
      {
        int rr;
        double* data = velnp_->Values();
        for(rr=0;rr<velnp_->MyLength();rr++)
        {
          printf("velnp %22.15e\n",data[rr]);
        }
      }
#endif

     #if 0  // DEBUG IO  --- traction
     {
	  int rr;
	  double* data = traction->Values();
	  for(rr=0;rr<traction->MyLength();rr++)
	  {
	      printf("traction %22.15e\n",data[rr]);
	  }
     }
     #endif

  return;
} // FluidImplicitTimeInt::Output


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |                                                             kue 04/07|
 -----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FluidImplicitTimeInt::ReadRestart(int step)
{
  IO::DiscretizationReader reader(discret_,step);
  time_ = reader.ReadDouble("time");
  step_ = reader.ReadInt("step");

  reader.ReadVector(velnp_,"velnp");
  reader.ReadVector(veln_, "veln");
  reader.ReadVector(velnm_,"velnm");
  reader.ReadVector(accn_ ,"accn");
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  set initial flow field for test cases                    gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FluidImplicitTimeInt::SetInitialFlowField(
  int whichinitialfield,
  int startfuncno
  )
{
  //------------------------------------------------------- beltrami flow
  if(whichinitialfield == 8)
  {
    const Epetra_Map* dofrowmap = discret_->DofRowMap();


    int err =0;

    const int numdim  = params_.get<int>("number of velocity degrees of freedom");
    const int npredof = numdim;

    double         p;
    vector<double> u  (numdim);
    vector<double> xyz(numdim);


    if(numdim!=3)
    {
      dserror("Beltrami flow is three dimensional flow!");
    }

    // set constants for analytical solution
    const double a      = PI/4.0;
    const double d      = PI/2.0;

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
        const int gid = nodedofset[nveldof];
        int lid = dofrowmap->LID(gid);
        err += velnp_->ReplaceMyValues(1,&(u[nveldof]),&lid);
        err += veln_ ->ReplaceMyValues(1,&(u[nveldof]),&lid);
        err += velnm_->ReplaceMyValues(1,&(u[nveldof]),&lid);
     }

      // initial pressure
      const int gid = nodedofset[npredof];
      int lid = dofrowmap->LID(gid);
      err += velnp_->ReplaceMyValues(1,&p,&lid);
      err += veln_ ->ReplaceMyValues(1,&p,&lid);
      err += velnm_->ReplaceMyValues(1,&p,&lid);

    } // end loop nodes lnodeid
    if(err!=0)
    {
      dserror("dof not on proc");
    }
  }
  else if(whichinitialfield==2 ||whichinitialfield==3)
  {
    const int numdim = params_.get<int>("number of velocity degrees of freedom");


    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = discret_->lRowNode(lnodeid);
      // the set of degrees of freedom associated with the node
      const vector<int> nodedofset = discret_->Dof(lnode);

      for(int index=0;index<numdim+1;++index)
      {
        int gid = nodedofset[index];

        double initialval=DRT::UTILS::FunctionManager::Instance().Funct(startfuncno-1).Evaluate(index,lnode->X());

        velnp_->ReplaceGlobalValues(1,&initialval,&gid);
        veln_ ->ReplaceGlobalValues(1,&initialval,&gid);
      }
    }
  }
  else
  {
    dserror("no other initial fields than zero, function and beltrami are available up to now");
  }

  return;
} // end SetInitialFlowField


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | evaluate error for test cases with analytical solutions   gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FluidImplicitTimeInt::EvaluateErrorComparedToAnalyticalSol()
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
 |  set default parameter list (static/public)               gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FluidImplicitTimeInt::SetDefaults(ParameterList& params)
{
  // number of degrees of freedom
  // if this one is not set afterwards, the algo hopefully crashes
  params.set<int> ("number of velocity degrees of freedom" ,-1);

  // evaluate error compared to analytical solution
  params.set<int>("eval err for analyt sol",0);

  // -------------------------------------------------- time integration
  // timestepsize
  params.set<double>           ("time step size"           ,0.01);
  // max. sim. time
  params.set<double>           ("total time"               ,0.0);
  // parameter for time-integration
  params.set<double>           ("theta"                    ,0.6667);
  // which kind of time-integration
  params.set<FLUID_TIMEINTTYPE>("time int algo"            ,timeint_one_step_theta);
  // bound for the number of timesteps
  params.set<int>              ("max number timesteps"     ,0);
  // number of steps with start algorithm
  params.set<int>              ("number of start steps"    ,0);
  // parameter for start algo
  params.set<double>           ("start theta"              ,0.6667);

  // ---------------------------------------------- nonlinear iteration
  // maximum number of nonlinear iteration steps
  params.set<int>             ("max nonlin iter steps"     ,5);
  // stop nonlinear iteration when both incr-norms are below this bound
  params.set<double>          ("tolerance for nonlin iter" ,1E-6);

  return;
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | solve stationary fluid problem   			               gjb 10/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FluidImplicitTimeInt::SolveStationaryProblem()
{
  // start time measurement for timeloop
  tm2_ref_ = rcp(new TimeMonitor(*timedynloop_));

  // set theta to one in order to avoid misuse
  theta_ = 1.0;


  // -------------------------------------------------------------------
  // pseudo time loop (continuation loop)
  // -------------------------------------------------------------------
  // slightly increasing b.c. values by given (pseudo-)timecurves to reach
  // convergence also for higher Reynolds number flows
  // as a side effect, you can do parameter studies for different Reynolds
  // numbers within only ONE simulation when you apply a proper
  // (pseudo-)timecurve

  while (step_< stepmax_)
  {
    // -------------------------------------------------------------------
    //              set (pseudo-)time dependent parameters
    // -------------------------------------------------------------------
       step_ += 1;
       time_ += dta_;

	// -------------------------------------------------------------------
	//                         out to screen
	// -------------------------------------------------------------------
	   if (myrank_==0)
	   {
	       printf("Stationary Fluid Solver - STEP = %4d/%4d \n",step_,stepmax_);
	   }

	// -------------------------------------------------------------------
	//         evaluate dirichlet and neumann boundary conditions
	// -------------------------------------------------------------------
	   {
	     ParameterList eleparams;
	     // action for elements
	     eleparams.set("action","calc_fluid_eleload");
	     // choose what to assemble
	     eleparams.set("assemble matrix 1",false);
	     eleparams.set("assemble matrix 2",false);
	     eleparams.set("assemble vector 1",true);
	     eleparams.set("assemble vector 2",false);
	     eleparams.set("assemble vector 3",false);
	     // other parameters needed by the elements
	     eleparams.set("total time",time_);
	     eleparams.set("delta time",dta_);
	     eleparams.set("thsl",1.0); // no timefac in stationary case
         eleparams.set("fs subgrid viscosity",fssgv_);

	     // set vector values needed by elements
	     discret_->ClearState();
	     discret_->SetState("velnp",velnp_);
	     // predicted dirichlet values
	     // velnp then also holds prescribed new dirichlet values
	     // dirichtoggle is 1 for dirichlet dofs, 0 elsewhere
	     discret_->EvaluateDirichlet(eleparams,*velnp_,*dirichtoggle_);
	     discret_->ClearState();

	     // evaluate Neumann b.c.
	     neumann_loads_->PutScalar(0.0);
	     discret_->EvaluateNeumann(eleparams,*neumann_loads_);
	     discret_->ClearState();
	   }

	   // compute an inverse of the dirichtoggle vector
	   invtoggle_->PutScalar(1.0);
	   invtoggle_->Update(-1.0,*dirichtoggle_,1.0);


    // -------------------------------------------------------------------
    //                     solve nonlinear equation system
    // -------------------------------------------------------------------
    this->NonlinearSolve();


    // -------------------------------------------------------------------
    //                         output of solution
    // -------------------------------------------------------------------
    this->Output();

  } // end of time loop

  // end time measurement for timeloop
  tm2_ref_ = null;

  return;
} // FluidImplicitTimeInt::SolveStationaryProblem


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  calculate (wall) stresses (public)                         gjb 07/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
RefCountPtr<Epetra_Vector> FluidImplicitTimeInt::CalcStresses()
{
     ParameterList eleparams;
     // set action for elements
     eleparams.set("action","integrate_Shapefunction");

     // get a vector layout from the discretization to construct matching
     // vectors and matrices
     //                 local <-> global dof numbering
     const Epetra_Map* dofrowmap = discret_->DofRowMap();

     // create vector (+ initialization with zeros)
     RefCountPtr<Epetra_Vector> integratedshapefunc = LINALG::CreateVector(*dofrowmap,true);

     // call loop over elements to evaluate the condition
     discret_->ClearState();
     const string condstring("FluidStressCalc");
     discret_->EvaluateCondition(eleparams,integratedshapefunc,condstring);
     discret_->ClearState();

     // compute traction values at specified nodes; otherwise do not touch the zero values
     for (int i=0;i<integratedshapefunc->MyLength();i++)
	  {
	      if ((*integratedshapefunc)[i] != 0.0)
	      {
	         // overwerite integratedshapefunc values with the calculated traction coefficients
	        (*integratedshapefunc)[i] = (*trueresidual_)[i]/(*integratedshapefunc)[i];
	      }
	  }

     // compute traction values at nodes   ---not used any more-----
     // component-wise calculation:  traction := 0.0*traction + 1.0*(trueresidual_ / integratedshapefunc_)
     // int err = traction->ReciprocalMultiply(1.0,*integratedshapefunc,*trueresidual_,0.0);
     // if (err > 0) dserror("error in traction->ReciprocalMultiply");


     return integratedshapefunc;

} // FluidImplicitTimeInt::CalcStresses()


/*----------------------------------------------------------------------*
 | Destructor dtor (public)                                  gammi 04/07|
 *----------------------------------------------------------------------*/
FluidImplicitTimeInt::~FluidImplicitTimeInt()
{
  return;
}



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | LiftDrag                                                  chfoe 11/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*!
\brief calculate lift&drag forces and angular momenta

Lift and drag forces are based upon the right hand side true-residual entities
of the corresponding nodes. The contribution of the end node of a line is entirely
added to a present L&D force.

Idea of this routine:

create

map< label, set<DRT::Node*> >

which is a set of nodes to each L&D Id
nodal forces of all the nodes within one set are added to one L&D force

Notice: Angular moments obtained from lift&drag forces currently refere to the
        initial configuration, i.e. are built with the coordinates X of a particular
        node irrespective of its current position.
*/
void FluidImplicitTimeInt::LiftDrag() const
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
        cout << "lift'n'drag Id      F_x             F_y             M_z :" << "\n";
      }
      if (ndim == 3)
      {
        cout << "lift'n'drag Id      F_x             F_y             F_z           ";
        cout << "M_x             M_y             M_z :" << "\n";
      }
    }

    // sort data
    for( unsigned i=0; i<ldconds.size(); ++i) // loop L&D conditions (i.e. lines in .dat file)
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

      /* put all nodes belonging to the L&D line or surface into 'nodes' which are
         associated with the present label */
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
      const std::set<DRT::Node*>& nodes = labelit->second; // pointer to nodeset of present label
      const int label = labelit->first;                    // the present label
      std::vector<double> values(6,0.0);             // vector with lift&drag forces
      std::vector<double> resultvec(6,0.0);          // vector with lift&drag forces after communication

      // get also pointer to centre coordinates
      const std::vector<double>* centerCoord = ldcoordmap[label];

      // loop all nodes within my set
      for( std::set<DRT::Node*>::const_iterator actnode = nodes.begin(); actnode != nodes.end(); ++actnode)
      {
        const double* x = (*actnode)->X(); // pointer to nodal coordinates
        const Epetra_BlockMap& rowdofmap = trueresidual_->Map();
        const std::vector<int> dof = discret_->Dof(*actnode);

        std::vector<double> distances (3);
        for (unsigned j=0; j<3; ++j)
        {
          distances[j]= x[j]-(*centerCoord)[j];
        }
        // get nodal forces
        const double fx = (*trueresidual_)[rowdofmap.LID(dof[0])];
        const double fy = (*trueresidual_)[rowdofmap.LID(dof[1])];
        const double fz = (*trueresidual_)[rowdofmap.LID(dof[2])];
        values[0] += fx;
        values[1] += fy;
        values[2] += fz;

        // calculate nodal angular momenta
        values[3] += distances[1]*fz-distances[2]*fy;
        values[4] += distances[2]*fx-distances[0]*fz;
        values[5] += distances[0]*fy-distances[1]*fx;
      } // end: loop over nodes

      // care for the fact that we are (most likely) parallel
      trueresidual_->Comm().SumAll (&(values[0]), &(resultvec[0]), 6);

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
}//FluidImplicitTimeInt::LiftDrag



#endif /* CCADISCRET       */
