/*!----------------------------------------------------------------------
\file condifimplicitintegration.cpp
\brief Control routine for convection-diffusion (in)stationary solvers,

     including instationary solvers based on

     o one-step-theta time-integration scheme

     o two-step BDF2 time-integration scheme 
       (with potential one-step-theta start algorithm)

     and stationary solver.

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "condifimplicitintegration.H"
#include "../drt_lib/drt_periodicbc.H"


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Constructor (public)                                        vg 05/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
CondifImplicitTimeInt::CondifImplicitTimeInt(RefCountPtr<DRT::Discretization> actdis,
                                             LINALG::Solver&                  solver,
                                             ParameterList&                   params,
                                             IO::DiscretizationWriter&        output) :
  // call constructor for "nontrivial" objects
  discret_(actdis),
  solver_ (solver),
  params_ (params),
  output_ (output),
  time_(0.0),
  step_(0),
  restartstep_(0),
  uprestart_(params.get("write restart every", -1)),
  writestep_(0),
  upres_(params.get("write solution every", -1))
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

  dtp_ = dta_  =params_.get<double>("time step size");

  // (fine-scale) subgrid diffusivity?
  fssgd_  =params_.get<int>("fs subgrid viscosity",0);

  // type of convective velocity field
  cdvel_  =params_.get<int>("condif velocity field");

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
  // get the processor ID from the communicator
  // -------------------------------------------------------------------
  myrank_  = discret_->Comm().MyPID();

  // -------------------------------------------------------------------
  // create empty system matrix --- stiffness and mass are assembled in
  // one system matrix!
  // -------------------------------------------------------------------

  // This is a first estimate for the number of non zeros in a row of
  // the matrix. Assuming a structured 3d-condif mesh we have 27 adjacent
  // nodes with 1 dof each.
  // We do not need the exact number here, just for performance reasons
  // a 'good' estimate
  maxentriesperrow_ = 27;

  // initialize standard (stabilized) system matrix
  sysmat_ = null;

  // -------------------------------------------------------------------
  // create empty vectors
  // -------------------------------------------------------------------

  // Vectors passed to the element
  // -----------------------------

  // temporal solution derivatives at time n and n-1
  phidtn_       = LINALG::CreateVector(*dofrowmap,true);
  phidtnm_      = LINALG::CreateVector(*dofrowmap,true);

  // solutions at time n+1, n and n-1
  phinp_        = LINALG::CreateVector(*dofrowmap,true);
  phin_         = LINALG::CreateVector(*dofrowmap,true);
  phinm_        = LINALG::CreateVector(*dofrowmap,true);

  // histvector --- a linear combination of phinm, phin (BDF)
  //                or phin, phidtn (One-Step-Theta)
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

  // The residual vector --- more or less the rhs
  residual_     = LINALG::CreateVector(*dofrowmap,true);

  // necessary only for the VM3 approach
  if (fssgd_ > 0)
  {
    // initialize (fine-scale) diffusivity matrix
    sysmat_sd_ = null;
  }

  // -------------------------------------------------------------------
  // create timers and time monitor
  // -------------------------------------------------------------------
  timedyntot_     = TimeMonitor::getNewTimer("dynamic routine total"     );
  timedyninit_    = TimeMonitor::getNewTimer(" + initial phase"          );
  timedynloop_    = TimeMonitor::getNewTimer(" + time loop"              );
  timeeleloop_    = TimeMonitor::getNewTimer("      + element calls"     );
  timeapplydirich_= TimeMonitor::getNewTimer("      + apply dirich cond.");
  timesolver_     = TimeMonitor::getNewTimer("      + solver calls"      );

  return;

} // CondifImplicitTimeInt::CondifImplicitTimeInt


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Start the time integration. Allows                                   |
 |                                                                      |
 |  o starting steps with different algorithms                          |
 |  o the "standard" time integration                                   |
 |                                                              vg 05/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void CondifImplicitTimeInt::Integrate()
{
  // time measurement --- start TimeMonitor tm0 and tm1
  tm0_ref_        = rcp(new TimeMonitor(*timedyntot_ ));

  // bound for the number of startsteps
  int    numstasteps         =params_.get<int>   ("number of start steps");

  
  if (timealgo_==timeint_stationary) // stationary case
    SolveStationaryProblem();
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
  cout<<endl<<endl;
  TimeMonitor::summarize();

  return;
} // CondifImplicitTimeInt::Integrate



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | contains the time loop                                       vg 05/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void CondifImplicitTimeInt::TimeLoop()
{
  // start time measurement for timeloop
  tm2_ref_ = rcp(new TimeMonitor(*timedynloop_));

  while (step_<stepmax_ and time_<maxtime_)
  {
    PrepareTimeStep();

    // -------------------------------------------------------------------
    //                     solve nonlinear equation
    // -------------------------------------------------------------------
    Solve();


    // -------------------------------------------------------------------
    //                         update solution
    //        current solution becomes old solution of next timestep
    //
    // One-step-Theta: (step>1)
    //
    //  phidtn_  = (phinp_-phin_) / (Theta * dt) - (1/Theta -1) * phin_
    //  "(n+1)"
    //
    //  phinm_ =phin_
    //  phin_  =phinp_
    //
    // BDF2:           (step>1)
    //
    //               2*dt(n)+dt(n-1)		  dt(n)+dt(n-1)
    //  phidtn_   = --------------------- phinp_ - --------------- phin_
    //             dt(n)*[dt(n)+dt(n-1)]	  dt(n)*dt(n-1)
    //
    //                     dt(n)
    //           + ----------------------- phinm_
    //             dt(n-1)*[dt(n)+dt(n-1)]
    //
    //
    //  phinm_ =phin_
    //  phin_  =phinp_
    //
    // BDF2 and  One-step-Theta: (step==1)
    //
    // The given formulas are only valid from the second timestep. In the
    // first step, the time derivative of phi is simply calculated as
    //
    //  phidtn_  = (phinp_-phin_) / (dt)
    //
    // -------------------------------------------------------------------
    TimeUpdate();


    // -------------------------------------------------------------------
    //                         output of solution
    // -------------------------------------------------------------------
    Output();


    // -------------------------------------------------------------------
    //                       update time step sizes
    // -------------------------------------------------------------------
    dtp_ = dta_;

  }

  // end time measurement for timeloop
  tm2_ref_ = null;

  return;
} // CondifImplicitTimeInt::TimeLoop


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | set part of the residual vector belonging to the old timestep        |
 |                                                              vg 05/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void CondifImplicitTimeInt::SetOldPartOfRighthandside()
{
  /*

  One-step-Theta:

                 hist_ = phin_ + dt*(1-Theta)*phidtn_


  BDF2: for constant time step:

                 hist_ = 4/3 phin_ - 1/3 phinm_

  */
  switch (timealgo_)
  {
  case timeint_one_step_theta: /* One step Theta time integration */
    hist_->Update(1.0, *phin_, dta_*(1.0-theta_), *phidtn_, 0.0);
    break;

  case timeint_bdf2:	/* 2nd order backward differencing BDF2	*/
    hist_->Update(4./3., *phin_, -1./3., *phinm_, 0.0);
    break;

  default:
    dserror("Time integration scheme unknown!");
  }
  return;
} // CondifImplicitTimeInt::SetOldPartOfRighthandside


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | setup the variables to do a new time step                    vg 08/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void CondifImplicitTimeInt::PrepareTimeStep()
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
    } /* end of switch(timealgo_) */
  }

  // -------------------------------------------------------------------
  // set part of the rhs vector belonging to the old timestep
  //
  //
  //         One-step-Theta:
  //
  //                 hist_ = phin_ + dta*(1-Theta)*phidtn_
  //
  //
  //         BDF2: for constant time step:
  //
  //                   hist_ = 4/3 phin_ - 1/3 phinm_
  //
  // -------------------------------------------------------------------
  SetOldPartOfRighthandside();

  // -------------------------------------------------------------------
  //         evaluate dirichlet and neumann boundary conditions
  // -------------------------------------------------------------------
  {
    ParameterList eleparams;
    // action for elements
    eleparams.set("action","calc_condif_eleload");
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
    eleparams.set("condif velocity field",cdvel_);

    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("phinp",phinp_);
    // predicted dirichlet values
     // phinp then also holds prescribed new dirichlet values
     // dirichtoggle is 1 for dirichlet dofs, 0 elsewhere
    discret_->EvaluateDirichlet(eleparams,*phinp_,*dirichtoggle_);
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
} // CondifImplicitTimeInt::PrepareTimeStep


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | contains the solver                                          vg 05/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void CondifImplicitTimeInt::Solve(
  bool is_stat //if true, stationary formulations are used in the element
  )
{
  const Epetra_Map* dofrowmap       = discret_->DofRowMap();

  double            dtsolve = 0;
  double            dtele   = 0;
  double            tcpu   ;

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
    eleparams.set("action","calc_condif_systemmat_and_residual");

    // other parameters that might be needed by the elements
    eleparams.set("total time",time_);
    eleparams.set("thsl",theta_*dta_);
    eleparams.set("condif velocity field",cdvel_);
    eleparams.set("fs subgrid viscosity",fssgd_);
    eleparams.set("using stationary formulation",is_stat);

    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("phinp",phinp_);
    discret_->SetState("hist"  ,hist_ );

    // decide whether VM3-based solution approach or standard approach
    if (fssgd_ > 0)
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

        // create (fine-scale) subgrid-diffusivity matrix
        sysmat_sd_ = LINALG::CreateMatrix(*dofrowmap,maxentriesperrow_); 

        // call loop over elements (two matrices + subgr.-visc.-scal. vector)
        discret_->Evaluate(eleparams,sysmat_,sysmat_sd_,residual_,sugrvisc_);
        discret_->ClearState();

        // finalize the (fine-scale) subgrid-diffusivity matrix
        LINALG::Complete(*sysmat_sd_);

        // apply DBC to (fine-scale) subgrid-diffusivity matrix
        LINALG::ApplyDirichlettoSystem(sysmat_sd_,phinp_,residual_,phinp_,dirichtoggle_);

        // call the VM3 constructor
        RCP<VM3_Solver> vm3_solver = rcp(new VM3_Solver(scalesep_,sysmat_sd_,sysmat_,sugrvisc_,zeros_,zeros_,zeros_,dirichtoggle_,mllist,compute) );

        // call the VM3 scale separator: precomputation of unscaled S^T*M*S
        vm3_solver->Separate(scalesep_,sysmat_sd_);
      }
      else
      {
        compute=false;

        // call loop over elements (one matrix + subgr.-visc.-scal. vector)
        discret_->Evaluate(eleparams,sysmat_,residual_,sugrvisc_);
        discret_->ClearState();
      }
      // call the VM3 constructor
      RCP<VM3_Solver> vm3_solver = rcp(new VM3_Solver(scalesep_,sysmat_sd_,sysmat_,sugrvisc_,zeros_,zeros_,zeros_,dirichtoggle_,mllist,compute) );

      // call the VM3 scaling:
      // scale precomputed matrix product by subgrid-viscosity-scaling vector
      vm3_solver->Scale(sysmat_sd_,sysmat_,sugrvisc_,zeros_,zeros_,zeros_,false );
    }
    else
    {
      // call standard loop over elements
      discret_->Evaluate(eleparams,sysmat_,residual_);
      discret_->ClearState();
    }

    // finalize the complete matrix
    LINALG::Complete(*sysmat_);

    // end time measurement for element call
    tm3_ref_=null;
    dtele=ds_cputime()-tcpu;
  }

  // Apply dirichlet boundary conditions to standard (stabilized) and complete 
  // system matrix
  {
    // start time measurement for application of dirichlet conditions
    tm4_ref_ = rcp(new TimeMonitor(*timeapplydirich_));

    LINALG::ApplyDirichlettoSystem(sysmat_,phinp_,residual_,phinp_,dirichtoggle_);

    // end time measurement for application of dirichlet conditions
    tm4_ref_=null;
  }

  //-------solve
  {
    // start time measurement for solver call
    tm5_ref_ = rcp(new TimeMonitor(*timesolver_));
    // get cpu time
    tcpu=ds_cputime();

    solver_.Solve(sysmat_,phinp_,residual_,true,true); 

    // end time measurement for solver call
    tm5_ref_=null;
    dtsolve=ds_cputime()-tcpu;
  }

  return;
} // CondifImplicitTimeInt::Solve

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                              vg 05/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void CondifImplicitTimeInt::TimeUpdate()
{

  // update time derivative of phi
  if (step_ == 1)
  {
    phidtnm_->PutScalar(0.0);

    // do just a linear interpolation within the first timestep
    phidtn_->Update( 1.0/dta_,*phinp_,1.0);

    phidtn_->Update(-1.0/dta_,*phin_ ,1.0);

    // ???
    phidtnm_->Update(1.0,*phidtn_,0.0);

  }
  else
  {
    // prev. time derivative of phi becomes (n-1)-time derivative of phi 
    // of next time step
    phidtnm_->Update(1.0,*phidtn_,0.0);

    /*

    One-step-Theta:

    phi(n+1) = (phi(n+1)-phi(n)) / (Theta * dt(n)) - (1/Theta -1) * phidt(n)


    BDF2:

                   2*dt(n)+dt(n-1)		    dt(n)+dt(n-1)
      phidt(n+1) = --------------------- phi(n+1) - --------------- phi(n)
                 dt(n)*[dt(n)+dt(n-1)]	            dt(n)*dt(n-1)

                         dt(n)
               + ----------------------- phi(n-1)
                 dt(n-1)*[dt(n)+dt(n-1)]

      */

      switch (timealgo_)
      {
          case timeint_one_step_theta: /* One step Theta time integration */
          {
            double fact1 = 1.0/(theta_*dta_);
            double fact2 =-1.0/theta_ +1.0;	/* = -1/Theta + 1		*/

            phidtn_->Update( fact1,*phinp_,0.0);
            phidtn_->Update(-fact1,*phin_ ,1.0);
            phidtn_->Update( fact2,*phidtnm_ ,1.0);

            break;
          }
          case timeint_bdf2:	/* 2nd order backward differencing BDF2	*/
          {
            if (dta_*dtp_ < EPS15)
              dserror("Zero time step size!!!!!");
            double sum = dta_ + dtp_;

            phidtn_->Update((2.0*dta_+dtp_)/(dta_*sum),*phinp_,
                                 - sum /(dta_*dtp_),*phin_ ,0.0);
            phidtn_->Update(dta_/(dtp_*sum),*phinm_,1.0);
          }
          break;
          default:
            dserror("Time integration scheme unknown for mass rhs!");
      }
    }

  // solution of this step becomes most recent solution of the last step
  phinm_->Update(1.0,*phin_ ,0.0);
  phin_ ->Update(1.0,*phinp_,0.0);

  return;
}// CondifImplicitTimeInt::TimeUpdate

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | output of solution vector to binio                           vg 05/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void CondifImplicitTimeInt::Output()
{

  //-------------------------------------------- output of solution
  //increase counters
  restartstep_ += 1;
  writestep_ += 1;

  if (writestep_ == upres_)  //write solution
    {
      writestep_= 0;

      output_.NewStep    (step_,time_);
      output_.WriteVector("phinp", phinp_);
      //output_.WriteVector("residual", residual_);

      if (restartstep_ == uprestart_) //add restart data
      {
        restartstep_ = 0;

        output_.WriteVector("phidtn", phidtn_);
        output_.WriteVector("phin", phin_);
        output_.WriteVector("phinm", phinm_);
      }
    }

  // write restart also when uprestart_ is not a integer multiple of upres_
  if ((restartstep_ == uprestart_) && (writestep_ > 0))
  {
    restartstep_ = 0;

    output_.NewStep    (step_,time_);
    output_.WriteVector("phinp", phinp_);
    //output_.WriteVector("residual", residual_);

    output_.WriteVector("phidtn", phidtn_);
    output_.WriteVector("phin", phin_);
    output_.WriteVector("phinm", phinm_);
  }


#if 0  // DEBUG IO --- the whole systemmatrix
      {
	  int rr;
	  int mm;
	  for(rr=0;rr<residual_->MyLength();rr++)
	  {
	      int NumEntries;

	      vector<double> Values(27);
	      vector<int> Indices(27);

	      sysmat_->ExtractGlobalRowCopy(rr,
					    27,
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
	      printf("rhs [%4d] %22.15e\n",rr,data[rr]);
	  }
      }
#endif

#if 0  // DEBUG IO --- neumann_loads
      {
	  int rr;
	  double* data = neumann_loads_->Values();
	  for(rr=0;rr<phinp_->MyLength();rr++)
	  {
	      printf("neum[%4d] %26.19e\n",rr,data[rr]);
	  }
      }

#endif


#if 1  //DEBUG IO --- incremental solution
      for (int proc=0; proc<phinp_->Comm().NumProc(); ++proc)
      {
        if (proc==myrank_)
        {
          printf("Proc %d\n",myrank_); fflush(stdout);
          for(int rr=0;rr<phinp_->MyLength();++rr)
          {
            printf("sol[%4d] %26.19e\n",phinp_->Map().GID(rr),(*phinp_)[rr]);
          }
        }
        fflush(stdout);
        phinp_->Comm().Barrier();
      }
#endif
//cout << *phinp_;

  return;
} // CondifImplicitTimeInt::Output


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |                                                             kue 04/07|
 -----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void CondifImplicitTimeInt::ReadRestart(int step)
{
  IO::DiscretizationReader reader(discret_,step);
  time_ = reader.ReadDouble("time");
  step_ = reader.ReadInt("step");

  reader.ReadVector(phinp_,"phinp");
  reader.ReadVector(phin_, "phin");
  reader.ReadVector(phinm_,"phinm");
  reader.ReadVector(phidtn_, "phidtn");
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | solve stationary convection-diffusion problem                vg 05/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void CondifImplicitTimeInt::SolveStationaryProblem()
{
  // start time measurement for stationary solver
  tm2_ref_ = rcp(new TimeMonitor(*timedynloop_));

  // -------------------------------------------------------------------
  //                         out to screen
  // -------------------------------------------------------------------
  if (myrank_==0) printf("Stationary Solver\n");
  step_=1;

  // -------------------------------------------------------------------
  //         evaluate dirichlet and neumann boundary conditions
  // -------------------------------------------------------------------
  {
     ParameterList eleparams;
     // action for elements
     eleparams.set("action","calc_condif_eleload");
     // choose what to assemble
     eleparams.set("assemble matrix 1",false);
     eleparams.set("assemble matrix 2",false);
     eleparams.set("assemble vector 1",true);
     eleparams.set("assemble vector 2",false);
     eleparams.set("assemble vector 3",false);
     // other parameter needed by the elements
     eleparams.set("condif velocity field",cdvel_);
     eleparams.set("fs subgrid viscosity",fssgd_);
     eleparams.set("using stationary formulation",true);

     // set vector values needed by elements
     discret_->ClearState();
     discret_->SetState("phinp",phinp_);
     // predicted dirichlet values
     // phinp then also holds prescribed new dirichlet values
     // dirichtoggle is 1 for dirichlet dofs, 0 elsewhere
     discret_->EvaluateDirichlet(eleparams,*phinp_,*dirichtoggle_);
     discret_->ClearState();

     // evaluate Neumann conditions
     neumann_loads_->PutScalar(0.0);
     discret_->EvaluateNeumann(eleparams,*neumann_loads_);
     discret_->ClearState();
  }

  //----------------------- compute an inverse of the dirichtoggle vector
  invtoggle_->PutScalar(1.0);
  invtoggle_->Update(-1.0,*dirichtoggle_,1.0);

  // -------------------------------------------------------------------
  //                     solve nonlinear equation
  // -------------------------------------------------------------------
  Solve(true);

  // -------------------------------------------------------------------
  //                         output of solution
  // -------------------------------------------------------------------
  Output();


  // end time measurement for timeloop
  tm2_ref_ = null;

  return;
} // CondifImplicitTimeInt::SolveStationaryProblem


/*----------------------------------------------------------------------*
 | Destructor dtor (public)                                     vg 05/07|
 *----------------------------------------------------------------------*/
CondifImplicitTimeInt::~CondifImplicitTimeInt()
{
  return;
}// CondifImplicitTimeInt::~CondifImplicitTimeInt::


#endif /* CCADISCRET       */
