/*!----------------------------------------------------------------------
\file
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
#ifdef TRILINOS_PACKAGE
#ifdef D_FLUID

#include "fluidimplicitintegration.H"

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
                                           DiscretizationWriter& output) :
  // call constructor for "nontrivial" objects
  discret_(actdis),
  solver_ (solver),
  params_ (params),
  output_ (output),
  time_(0.0),
  step_(0),
  restartstep_(0),
  uprestart_(params.get("write restart every", -1))
{

  int numdim = params_.get<int>("number of velocity degrees of freedom");

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
    // Allocate integer vectors which will hold the dof number of the
    // velocity or pressure dofs
    vector<int> velmapdata;
    vector<int> premapdata;

    velmapdata.reserve(discret_->NumMyRowNodes()*numdim);
    premapdata.reserve(discret_->NumMyRowNodes());

    for (int i=0; i<discret_->NumMyRowNodes(); ++i)
    {
      DRT::Node* node = discret_->lRowNode(i);
      for (int j=0; j<numdim; ++j)
      {
	  // add this velocity dof to the velmapdata vector
	  velmapdata.push_back(node->Dof()[j]);
      }
      // add this pressure dof to the premapdata vector
      premapdata.push_back(node->Dof()[numdim]);
    }

    // the rowmaps are generated according to the pattern provided by
    // the data vectors
    velrowmap_ = rcp(new Epetra_Map(discret_->NumGlobalNodes()*numdim,
                                    velmapdata.size(),&velmapdata[0],0,
                                    discret_->Comm()));
    prerowmap_ = rcp(new Epetra_Map(discret_->NumGlobalNodes(),
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

  // histvector --- a linear combination of velnm, veln (BDF)
  //                or veln, accn (One-Step-Theta)
  hist_         = LINALG::CreateVector(*dofrowmap,true);

  // Vectors associated to boundary conditions
  // -----------------------------------------

  // toggle vector indicating which dofs have Dirichlet BCs
  dirichtoggle_ = LINALG::CreateVector(*dofrowmap,true);

  // a vector of zeros to be used to enforce zero dirichlet boundary conditions
  zeros_        = LINALG::CreateVector(*dofrowmap,true);

  // the vector containing body and surface forces
  neumann_loads_= LINALG::CreateVector(*dofrowmap,true);

  // Vectors used for solution process 
  // ---------------------------------

  // The residual vector --- more or less the rhs for the incremental
  // formulation!!!
  residual_     = LINALG::CreateVector(*dofrowmap,true);

  // Nonlinear iteration increment vector
  incvel_       = LINALG::CreateVector(*dofrowmap,true);


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
  tm1_ref_        = rcp(new TimeMonitor(*timedyninit_));

  // initialise some variables
  double dta  = 0.0;
  double dtp  = 0.0;


  // ----------------------------------------------------- stop criteria
  // bound for the number of startsteps
  int    numstasteps         =params_.get<int>   ("number of start steps");
  // bound for the number of timesteps
  int    stepmax             =params_.get<int>   ("max number timesteps");
  // max. sim. time
  double maxtime             =params_.get<double>("total time");

  // parameter for time-integration
  double            theta    =params_.get<double>("theta");
  // which kind of time-integration
  FLUID_TIMEINTTYPE timealgo =params_.get<FLUID_TIMEINTTYPE>("time int algo");

  // parameter for start algorithm
  double starttheta          =params_.get<double>("start theta");
  FLUID_TIMEINTTYPE startalgo=timeint_one_step_theta;

  dtp = dta  =params_.get<double>("time step size");

  // time measurement --- this causes the TimeMonitor tm1 to stop here
  //                                                (call of destructor)
  tm1_ref_ = null;

  // start procedure
  if (step_<numstasteps)
  {
    if(numstasteps>stepmax)
    {
      dserror("more startsteps than steps");
    }
    
    this->TimeIntegrateFromTo(
      step_,
      time_,
      dta,
      dtp,
      numstasteps,
      maxtime,
      startalgo,
      starttheta
      );
  }

  // continue with the final time integration
  this->TimeIntegrateFromTo(
    step_,
    time_,
    dta,
    dtp,
    stepmax,
    maxtime,
    timealgo,
    theta
    );

  // print the results of time measurements

  tm0_ref_ = null; // end total time measurement
  cout<<endl<<endl;
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
void FluidImplicitTimeInt::TimeIntegrateFromTo(
  int&               step,
  double&            time,
  double&            dta,
  double&            dtp,
  int                endstep,
  double             endtime,
  FLUID_TIMEINTTYPE  timealgo,
  double             theta
  )
{
  // start time measurement for timeloop
  tm2_ref_ = rcp(new TimeMonitor(*timedynloop_));

  bool stop_timeloop=false;

  while (stop_timeloop==false)
  {
    // -------------------------------------------------------------------
    //              set time dependent parameters
    // -------------------------------------------------------------------
    step++;

    time+=dta;

    // for bdf2 theta is set  by the timestepsizes, 2/3 for const. dt
    if(timealgo==timeint_bdf2)
    {
      theta= (dta+dtp)/(2.0*dta + dtp);
    }
    
    // -------------------------------------------------------------------
    //                         out to screen
    // -------------------------------------------------------------------
    if (myrank_==0)
    {
      switch(timealgo)
      {
          case timeint_one_step_theta:
            printf("TIME: %11.4E/%11.4E  DT = %11.4E  One-Step-Theta  STEP = %4d/%4d \n",
                   time,endtime,dta,step,endstep);
            break;
          case timeint_bdf2:
            printf("TIME: %11.4E/%11.4E  DT = %11.4E     BDF2         STEP = %4d/%4d \n",
                   time,endtime,dta,step,endstep);
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
    this->SetOldPartOfRighthandside(timealgo,dta,theta);

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
    if (step>1)
    {
      this->ExplicitPredictor(dta,dtp);
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
     eleparams.set("total time",time);
     eleparams.set("delta time",dta);
     eleparams.set("time constant for integration",theta*dta);

     // set vector values needed by elements
     discret_->ClearState();
     discret_->SetState("u and p at time n+1 (trial)",velnp_);
     // predicted dirichlet values
     // velnp then also holds prescribed new dirichlet values
     // dirichtoggle is 1 for dirichlet dofs, 0 elsewhere
     discret_->EvaluateDirichlet(eleparams,*velnp_,*dirichtoggle_);
     discret_->ClearState();

     // evaluate Neumann conditions
     eleparams.set("total time",time);
     eleparams.set("time constant for integration",theta*dta);

     neumann_loads_->PutScalar(0.0);
     discret_->EvaluateNeumann(eleparams,*neumann_loads_);
     discret_->ClearState();
    }

    // -------------------------------------------------------------------
    //                     solve nonlinear equation
    // -------------------------------------------------------------------
    this->NonlinearSolve(dta,theta);


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
    this->TimeUpdate(timealgo,step,dta,dtp,theta);

    
    // -------------------------------------------------------------------
    // evaluate error for test flows with analytical solutions
    // -------------------------------------------------------------------
    this->EvaluateErrorComparedToAnalyticalSol(time);

    
    // -------------------------------------------------------------------
    //                         output of solution
    // -------------------------------------------------------------------
    this->Output(step,time);


    // -------------------------------------------------------------------
    //                       update time step sizes
    // -------------------------------------------------------------------
    dtp = dta;


    // -------------------------------------------------------------------
    //                    stop criterium for timeloop
    // -------------------------------------------------------------------

    if(step==endstep||time>=endtime)
    {
	stop_timeloop=true;
    }

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
void FluidImplicitTimeInt::SetOldPartOfRighthandside(
  FLUID_TIMEINTTYPE time_algo,
  double            dta,
  double            theta
  )
{
  /*

  One-step-Theta:

                 hist_ = veln_ + dt*(1-Theta)*accn_


  BDF2: for constant time step:

                 hist_ = 4/3 veln_ - 1/3 velnm_

  */
  switch (time_algo)
  {
      case timeint_one_step_theta: /* One step Theta time integration */

        hist_->PutScalar(0.0);
        hist_->Update(               1.0,*veln_,1.0);
        hist_->Update(dta * (1.0 -theta),*accn_,1.0);
        break;

      case timeint_bdf2:	/* 2nd order backward differencing BDF2	*/
        hist_->Update(4./3.,*veln_,-1./3.,*velnm_,0.0);
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
void FluidImplicitTimeInt::ExplicitPredictor(
  double dta,
  double dtp
  )
{
  double fact1 = dta*(1.0+dta/dtp);
  double fact2 = DSQR(dta/dtp);

  velnp_->Update( fact1,*accn_ ,1.0);
  velnp_->Update(-fact2,*veln_ ,1.0);
  velnp_->Update( fact2,*velnm_,1.0);

  return;
} // FluidImplicitTimeInt::ExplicitPredictor


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | contains the nonlinear iteration loop                     gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FluidImplicitTimeInt::NonlinearSolve(
  double dta,
  double theta
  )
{
  // start time measurement for nonlinear iteration
  tm6_ref_ = rcp(new TimeMonitor(*timenlnloop_));

  const Epetra_Map* dofrowmap       = discret_->DofRowMap();

  int               itnum         = 0;
  bool              stopnonliniter = false;

  double            dtsolve;
  double            dtele  ;
  double            tcpu   ;

  if(myrank_ == 0)
  {
    printf("+------------+-------------------+--------------+--------------+ \n");
    printf("|- step/max -|- tol      [norm] -|- vel-error --|- pre-error --| \n");
  }

  while (stopnonliniter==false)
  {
    itnum++;

    // -------------------------------------------------------------------
    // call elements to calculate system matrix
    // -------------------------------------------------------------------
    {
      // start time measurement for element call
      tm3_ref_ = rcp(new TimeMonitor(*timeeleloop_));
      // get cpu time
      tcpu=ds_cputime();

      // zero out the stiffness matrix
      sysmat_ = LINALG::CreateMatrix(*dofrowmap,maxentriesperrow_);
      // zero out residual
      residual_->PutScalar(0.0);
      residual_->Update(1.0,*neumann_loads_,0.0);

      // create the parameters for the discretization
      ParameterList eleparams;

      // action for elements
      eleparams.set("action","calc_fluid_systemmat_and_residual");
      // choose what to assemble
      eleparams.set("assemble matrix 1",true);
      eleparams.set("assemble matrix 2",false);
      eleparams.set("assemble vector 1",true);
      eleparams.set("assemble vector 2",false);
      eleparams.set("assemble vector 3",false);
      // other parameters that might be needed by the elements
      eleparams.set("time constant for integration",theta*dta);
      // set vector values needed by elements
      discret_->ClearState();
      discret_->SetState("u and p at time n+1 (trial)",velnp_);
      discret_->SetState("old solution data for rhs"  ,hist_ );

      // call loop over elements
      {
        discret_->Evaluate(eleparams,sysmat_,null,residual_,null,null);
        discret_->ClearState();
      }

      // finalize the system matrix
      LINALG::Complete(*sysmat_);
      maxentriesperrow_ = sysmat_->MaxNumEntries();

      // end time measurement for element call
      tm3_ref_=null;
      dtele=ds_cputime()-tcpu;
    }

    //--------- Apply dirichlet boundary conditions to system of equations
    //          residual discplacements are supposed to be zero at
    //          boundary conditions
    incvel_->PutScalar(0.0);
    zeros_->PutScalar(0.0);
    {
      // start time measurement for application of dirichlet conditions
      tm4_ref_ = rcp(new TimeMonitor(*timeapplydirich_));

      LINALG::ApplyDirichlettoSystem(sysmat_,incvel_,residual_,
				     zeros_,dirichtoggle_);

      // end time measurement for application of dirichlet conditions
      tm4_ref_=null;
    }

    //-------solve for residual displacements to correct incremental displacements
    {
      // start time measurement for element call
      tm5_ref_ = rcp(new TimeMonitor(*timesolver_));
      // get cpu time
      tcpu=ds_cputime();

      bool initsolver = false;
      if (itnum==1) // init solver in first iteration only
      {
        initsolver = true;
      }

      solver_.Solve(sysmat_,incvel_,residual_,true,initsolver);

      // end time measurement for application of dirichlet conditions
      tm5_ref_=null;
      dtsolve=ds_cputime()-tcpu;
    }
    //------------------------------------------------ update (u,p) trial
    velnp_->Update(1.0,*incvel_,1.0);

    //------------------------------------------------- check convergence
    this->NonlinearConvCheck(stopnonliniter,itnum,dtele,dtsolve);
  }

  // end time measurement for nonlinear iteration
  tm6_ref_ = null;

  return;
} // FluidImplicitTimeInt::NonlinearSolve

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
  int     itemax    =params_.get<int>   ("max nonlin iter steps");

  // ------------------------------- stop nonlinear iteration when both
  //                                 increment-norms are below this bound
  double  ittol     =params_.get<double>("tolerance for nonlin iter");


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
void FluidImplicitTimeInt::TimeUpdate(
  FLUID_TIMEINTTYPE time_algo,
  int               step,
  double            dta ,
  double            dtp ,
  double            theta
  )
{

  // update acceleration
  if (step == 1)
  {
    accnm_->PutScalar(0.0);

    // do just a linear interpolation within the first timestep
    accn_->Update( 1.0/dta,*velnp_,1.0);

    accn_->Update(-1.0/dta,*veln_ ,1.0);

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

      switch (time_algo)
      {
          case timeint_one_step_theta: /* One step Theta time integration */
          {
            double fact1 = 1.0/(theta*dta);
            double fact2 =-1.0/theta +1.0;	/* = -1/Theta + 1		*/

            accn_->Update( fact1,*velnp_,0.0);
            accn_->Update(-fact1,*veln_ ,1.0);
            accn_->Update( fact2,*accnm_ ,1.0);

            break;
          }
          case timeint_bdf2:	/* 2nd order backward differencing BDF2	*/
          {
            if (dta*dtp < EPS15)
              dserror("Zero time step size!!!!!");
            double sum = dta + dtp;

            accn_->Update((2.0*dta+dtp)/(dta*sum),*velnp_,
                                 - sum /(dta*dtp),*veln_ ,0.0);
            accn_->Update(dta/(dtp*sum),*velnm_,1.0);
          }
          break;
          default:
            dserror("Time integration scheme unknown for mass rhs!");
      }
    }

  // solution of this step becomes most recent solution of the last step
  velnm_->Update(1.0,*veln_ ,0.0);
  veln_ ->Update(1.0,*velnp_,0.0);

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
void FluidImplicitTimeInt::Output(
  int    step,
  double time
  )
{

  //-------------------------------------------- output of solution
  output_.NewStep    (step,time);
  output_.WriteVector("velnp", velnp_);

  // do restart if we have to
  restartstep_ += 1;
  if (restartstep_ == uprestart_)
  {
    restartstep_ = 0;

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
  DiscretizationReader reader(discret_,step);
  time_ = reader.ReadDouble("time");
  step_ = reader.ReadInt("step");

  reader.ReadVector(velnp_,"velnp");
  reader.ReadVector(veln_, "veln");
  reader.ReadVector(velnm_,"velnm");
  reader.ReadVector(accn_, "accn");
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
  int whichinitialfield
  )
{
  //------------------------------------------------------- beltrami flow
  if(whichinitialfield == 8)
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
      DRT::DofSet nodedofset = lnode->Dof();

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
        err += velnm_->ReplaceMyValues(1,&(u[nveldof]),&lid);
     }
      
      // initial pressure
      gid = nodedofset[npredof];
      lid = dofrowmap->LID(gid);
      err += velnp_->ReplaceMyValues(1,&p,&lid);
      err += veln_ ->ReplaceMyValues(1,&p,&lid);
      err += velnm_->ReplaceMyValues(1,&p,&lid);

    } // end loop nodes lnodeid
    if(err!=0)
    {
      dserror("dof not on proc");
    }
  }
  else
  {
    dserror("no other initial fields than zero and beltrami are available up to now");
  }
  
  return;
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | evaluate error for test cases with analytical solutions   gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FluidImplicitTimeInt::EvaluateErrorComparedToAnalyticalSol(
  double time
  )
{

  int calcerr = params_.get<int>("eval err for analyt sol");
  
  //------------------------------------------------------- beltrami flow
  switch (calcerr)
  {
      case 0:
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
        eleparams.set("total time",time);
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
}


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


/*----------------------------------------------------------------------*
 | Destructor dtor (public)                                  gammi 04/07|
 *----------------------------------------------------------------------*/
FluidImplicitTimeInt::~FluidImplicitTimeInt()
{
  return;
}// FluidImplicitTimeInt::~FluidImplicitTimeInt::


#endif /* D_FLUID          */
#endif /* TRILINOS_PACKAGE */
#endif /* CCADISCRET       */
