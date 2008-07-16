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

#include "../drt_lib/drt_globalproblem.H"
#include "scatra_timint_implicit.H"
#include "../drt_fluid/drt_periodicbc.H"
#include "../drt_fluid/vm3_solver.H"
#include "../drt_lib/drt_function.H"

// for the condition writer output
/*
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/adapter_utils.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io.H"
*/

/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01    |
  | general problem data                                                 |
  | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Constructor (public)                                        vg 05/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
ScaTraImplicitTimeInt::ScaTraImplicitTimeInt(RCP<DRT::Discretization>      actdis,
                                             RCP<LINALG::Solver>           solver,
                                             RCP<ParameterList>            params,
                                             RCP<IO::DiscretizationWriter> output) :
  // call constructor for "nontrivial" objects
  discret_(actdis),
  solver_ (solver),
  params_ (params),
  output_ (output),
  time_(0.0),
  step_(0),
  restartstep_(0),
  uprestart_(params->get("write restart every", -1)),
  writestep_(0),
  upres_(params->get("write solution every", -1))
{

  // -------------------------------------------------------------------
  // create timers and time monitor
  // -------------------------------------------------------------------
  timetotal_    = TimeMonitor::getNewTimer("dynamic routine total");
  timeinit_     = TimeMonitor::getNewTimer(" + initialization"    );
  timetimeloop_ = TimeMonitor::getNewTimer(" + time loop"         );
  timeelement_  = TimeMonitor::getNewTimer("      + element calls");
  timeavm3_     = TimeMonitor::getNewTimer("           + avm3"    );
  timeapplydbc_ = TimeMonitor::getNewTimer("      + apply DBC"    );
  timesolver_   = TimeMonitor::getNewTimer("      + solver calls" );

  // time measurement: total --- start TimeMonitor tm0
  tm0_ref_        = rcp(new TimeMonitor(*timetotal_ ));

  // time measurement: initialization --- start TimeMonitor tm7
  tm1_ref_        = rcp(new TimeMonitor(*timeinit_ ));

  // -------------------------------------------------------------------
  // get the basic parameters first
  // -------------------------------------------------------------------

  // bound for the number of timesteps
  stepmax_                   =params_->get<int>   ("max number timesteps");
  // max. sim. time
  maxtime_                   =params_->get<double>("total time");

  // parameter for time-integration
  theta_                     =params_->get<double>("theta");
  // which kind of time-integration
  timealgo_                  =params_->get<INPUTPARAMS::ScaTraTimeIntegrationScheme>("time int algo");

  dtp_ = dta_  =params_->get<double>("time step size");

  // (fine-scale) subgrid diffusivity?
  fssgd_  =params_->get<string>("fs subgrid diffusivity","No");

  // type of convective velocity field
  cdvel_  =params_->get<int>("velocity field");

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
  // the matrix. Assuming a structured 3d mesh we have 27 adjacent
  // nodes with 1 dof each.
  // We do not need the exact number here, just for performance reasons
  // a 'good' estimate

  // initialize standard (stabilized) system matrix
  sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,27));

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

  // get noderowmap of discretization
  const Epetra_Map* noderowmap = discret_->NodeRowMap();
  /// convective velocity (always three velocity components per node)
  convel_         = rcp(new Epetra_MultiVector(*noderowmap,3,true));

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

  // -------------------------------------------------------------------
  // necessary only for the VM3 approach:
  // initialize subgrid-diffusivity matrix + respective ouptput
  // -------------------------------------------------------------------
  if (fssgd_ != "No")
  {
    sysmat_sd_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,27));

    if (myrank_ == 0)
    {
      // Output
      cout << "Fine-scale subgrid-diffusivity approach based on AVM3: ";
      cout << &endl << &endl;
      cout << params_->get<string>("fs subgrid diffusivity");
      cout << &endl << &endl;
    }
  }

  // set initial field
  SetInitialField(params_->get<int>("scalar initial field"), params_->get<int>("scalar initial field func number"));

  // end time measurement for initialization
  tm1_ref_ = null;

  return;

} // ScaTraImplicitTimeInt::ScaTraImplicitTimeInt


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
void ScaTraImplicitTimeInt::Integrate()
{
  // bound for the number of startsteps
  int    numstasteps         =params_->get<int>   ("number of start steps");

  if (timealgo_==INPUTPARAMS::timeint_stationary) // stationary case
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
    TimeLoop();
  }

  // end total time measurement
  tm0_ref_ = null;

  // print the results of time measurements
  //cout<<endl<<endl;
  TimeMonitor::summarize();

  return;
} // ScaTraImplicitTimeInt::Integrate



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | contains the time loop                                       vg 05/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void ScaTraImplicitTimeInt::TimeLoop()
{
  // time measurement: time loop --- start TimeMonitor tm2
  tm2_ref_ = rcp(new TimeMonitor(*timetimeloop_));

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
    Update();


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
} // ScaTraImplicitTimeInt::TimeLoop


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
void ScaTraImplicitTimeInt::SetOldPartOfRighthandside()
{
  /*

  One-step-Theta:

                 hist_ = phin_ + dt*(1-Theta)*phidtn_


  BDF2: for constant time step:

                 hist_ = 4/3 phin_ - 1/3 phinm_

  */
  switch (timealgo_)
  {
  case INPUTPARAMS::timeint_one_step_theta: /* One step Theta time integration */
    hist_->Update(1.0, *phin_, dta_*(1.0-theta_), *phidtn_, 0.0);
    break;

  case INPUTPARAMS::timeint_bdf2:	/* 2nd order backward differencing BDF2	*/
    hist_->Update(4./3., *phin_, -1./3., *phinm_, 0.0);
    break;

  default:
    dserror("Time integration scheme unknown!");
  }
  return;
} // ScaTraImplicitTimeInt::SetOldPartOfRighthandside


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | setup the variables to do a new time step                    vg 08/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void ScaTraImplicitTimeInt::PrepareTimeStep()
{
  // -------------------------------------------------------------------
  //              set time dependent parameters
  // -------------------------------------------------------------------
  step_ += 1;

  time_ += dta_;

  // for bdf2 theta is set  by the timestepsizes, 2/3 for const. dt
  if (timealgo_==INPUTPARAMS::timeint_bdf2)
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
    case INPUTPARAMS::timeint_one_step_theta:
      printf("TIME: %11.4E/%11.4E  DT = %11.4E  One-Step-Theta  STEP = %4d/%4d \n",
             time_,maxtime_,dta_,step_,stepmax_);
      break;
    case INPUTPARAMS::timeint_bdf2:
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
    discret_->EvaluateDirichlet(eleparams,phinp_,null,null,dirichtoggle_);
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
} // ScaTraImplicitTimeInt::PrepareTimeStep


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | contains the solver                                          vg 05/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void ScaTraImplicitTimeInt::Solve(
  bool is_stat //if true, stationary formulations are used in the element
  )
{
  double tcpu;
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // -------------------------------------------------------------------
  // call elements to calculate system matrix
  // -------------------------------------------------------------------
  {
    // time measurement: element --- start TimeMonitor tm3
    tm3_ref_ = rcp(new TimeMonitor(*timeelement_));
    // get cpu time
    tcpu=ds_cputime();

    sysmat_->Zero();

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
    eleparams.set("fs subgrid diffusivity",fssgd_);
    eleparams.set("using stationary formulation",is_stat);

    //provide velocity field (export to column map necessary for parallel evaluation)
    const Epetra_Map* nodecolmap = discret_->NodeColMap();
    RefCountPtr<Epetra_MultiVector> tmp = rcp(new Epetra_MultiVector(*nodecolmap,3));
    LINALG::Export(*convel_,*tmp);
    eleparams.set("velocity field",tmp);

    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("phinp",phinp_);
    discret_->SetState("hist" ,hist_ );

#if 0
    // AVMS solver: VM3 solution approach with extended matrix system
    // (may replace direct algebraic VM3 solution approach)
    // begin first encapsulation of AVMS solution approach
    if (fssgd_ != "No")
    {
      // create all-scale subgrid-diffusivity matrix
      sysmat_sd_->Zero();

      // call loop over elements (two matrices)
      discret_->Evaluate(eleparams,sysmat_,sysmat_sd_,residual_);
      discret_->ClearState();
    }
    // end first encapsulation of AVMS solution approach
#endif
    // direct algebraic VM3 solution approach
    // begin encapsulation of direct algebraic VM3 solution approach
    if (fssgd_ != "No")
    {
      // time measurement: avm3 --- start TimeMonitor tm4
      tm4_ref_ = rcp(new TimeMonitor(*timeavm3_));

      // extract the ML parameters
      ParameterList&  mllist = solver_->Params().sublist("ML Parameters");

      // subgrid-viscosity-scaling vector
      sugrvisc_ = LINALG::CreateVector(*dofrowmap,true);

      if (step_ == 1)
      {
        // create normalized all-scale subgrid-diffusivity matrix
        sysmat_sd_->Zero();

        // end time measurement for avm3 evaluation
        tm4_ref_=null;

        // call loop over elements (two matrices + subgr.-visc.-scal. vector)
        discret_->Evaluate(eleparams,sysmat_,sysmat_sd_,residual_,sugrvisc_);
        discret_->ClearState();

        // time measurement: avm3 --- start TimeMonitor tm4
        tm4_ref_ = rcp(new TimeMonitor(*timeavm3_));

        // finalize the normalized all-scale subgrid-diffusivity matrix
        sysmat_sd_->Complete();

        // apply DBC to normalized all-scale subgrid-diffusivity matrix
        LINALG::ApplyDirichlettoSystem(sysmat_sd_,phinp_,residual_,phinp_,dirichtoggle_);

        // call the VM3 constructor
        vm3_solver_ = rcp(new VM3_Solver(sysmat_sd_,dirichtoggle_,mllist,true,false) );
      }
      else
      {
        // end time measurement for avm3 evaluation
        tm4_ref_=null;

        // call loop over elements (one matrix + subgr.-visc.-scal. vector)
        discret_->Evaluate(eleparams,sysmat_,null,residual_,sugrvisc_);
        discret_->ClearState();

        // time measurement: avm3 --- start TimeMonitor tm4
        tm4_ref_ = rcp(new TimeMonitor(*timeavm3_));
      }
      // check whether VM3 solver exists
      if (vm3_solver_ == null) dserror("vm3_solver not allocated");

      // call the VM3 scaling:
      // scale precomputed matrix product by subgrid-viscosity-scaling vector
      vm3_solver_->Scale(sysmat_sd_,sysmat_,zeros_,zeros_,sugrvisc_,zeros_,false );

      // end time measurement for avm3 evaluation
      tm4_ref_=null;
    }
    // end encapsulation of direct algebraic VM3 solution approach
    else
    {
      // call standard loop over elements
      discret_->Evaluate(eleparams,sysmat_,residual_);
      discret_->ClearState();
    }

    // finalize the complete matrix
    sysmat_->Complete();

    // end time measurement for element
    tm3_ref_=null;
    dtele_=ds_cputime()-tcpu;
  }

  // Apply dirichlet boundary conditions to system matrix
  {
    // time measurement: application of dbc --- start TimeMonitor tm5
    tm5_ref_ = rcp(new TimeMonitor(*timeapplydbc_));

    LINALG::ApplyDirichlettoSystem(sysmat_,phinp_,residual_,phinp_,dirichtoggle_);

    // end time measurement for application of dbc
    tm5_ref_=null;
  }

  //-------solve
  {
    // time measurement: solver --- start TimeMonitor tm6
    tm6_ref_ = rcp(new TimeMonitor(*timesolver_));
    // get cpu time
    tcpu=ds_cputime();

#if 0
    // AVMS solver: VM3 solution approach with extended matrix system
    // (may replace direct algebraic VM3 solution approach)
    // begin second encapsulation of AVMS solution approach
    if (fssgd_ != "No")
    {
      // add standard matrix to subgrid-diffusivity matrix: fine-scale matrix
      sysmat_sd_->Add(*sysmat_,false,1.0,1.0);

      // finalize fine-scale matrix
      sysmat_sd_->Complete();

      // apply DBC to fine-scale matrix
      LINALG::ApplyDirichlettoSystem(sysmat_sd_,phinp_,residual_,phinp_,dirichtoggle_);

      // extract the ML parameters
      ParameterList&  mllist = solver_.Params().sublist("ML Parameters");

      // call the AVMS constructor
      avms_solver_ = rcp(new AVMS_Solver(sysmat_sd_,sysmat_,dirichtoggle_,mllist) );

      // Apply the AVMS solver
      avms_solver_->Solve(*residual_,*phinp_,solver_.Params());
    }
    else
    // end second encapsulation of AVMS solution approach
#endif
      solver_->Solve(sysmat_->EpetraMatrix(),phinp_,residual_,true,true);

    // end time measurement for solver
    tm6_ref_=null;
    dtsolve_=ds_cputime()-tcpu;
  }

  return;
} // ScaTraImplicitTimeInt::Solve

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
void ScaTraImplicitTimeInt::Update()
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
          case INPUTPARAMS::timeint_one_step_theta: /* One step Theta time integration */
          {
            double fact1 = 1.0/(theta_*dta_);
            double fact2 =-1.0/theta_ +1.0;	/* = -1/Theta + 1		*/

            phidtn_->Update( fact1,*phinp_,0.0);
            phidtn_->Update(-fact1,*phin_ ,1.0);
            phidtn_->Update( fact2,*phidtnm_ ,1.0);

            break;
          }
          case INPUTPARAMS::timeint_bdf2:	/* 2nd order backward differencing BDF2	*/
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
}// ScaTraImplicitTimeInt::TimeUpdate

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | output of solution vector to binio                           vg 05/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void ScaTraImplicitTimeInt::Output()
{
  //-------------------------------------------- output of solution
  //increase counters
  restartstep_ += 1;
  writestep_ += 1;

  if (writestep_ == upres_)  //write solution
    {
      writestep_= 0;

      output_->NewStep    (step_,time_);
      output_->WriteVector("phinp", phinp_);
      output_->WriteVector("convec_velocity", convel_,IO::DiscretizationWriter::nodevector);
      //output_->WriteVector("residual", residual_);

      RCP<Epetra_MultiVector> flux = CalcFlux();
      output_->WriteVector("flux", flux, IO::DiscretizationWriter::nodevector);

      // write domain decomposition for visualization (only once!)
      if (step_==upres_)
       output_->WriteElementData();

      if (restartstep_ == uprestart_) //add restart data
      {
        restartstep_ = 0;

        output_->WriteVector("phidtn", phidtn_);
        output_->WriteVector("phin", phin_);
        output_->WriteVector("phinm", phinm_);
      }
    }

  // write restart also when uprestart_ is not a integer multiple of upres_
  if ((restartstep_ == uprestart_) && (writestep_ > 0))
  {
    restartstep_ = 0;

    output_->NewStep    (step_,time_);
    output_->WriteVector("phinp", phinp_);
    output_->WriteVector("velocity", convel_,IO::DiscretizationWriter::nodevector);
    //output_->WriteVector("residual", residual_);

    output_->WriteVector("phidtn", phidtn_);
    output_->WriteVector("phin", phin_);
    output_->WriteVector("phinm", phinm_);
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


#if 0  //DEBUG IO --- incremental solution
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

#if 0  //DEBUG IO --- convective velocity
convel_->Print(cout);
#endif

  return;
} // ScaTraImplicitTimeInt::Output


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |                                                             kue 04/07|
 -----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void ScaTraImplicitTimeInt::ReadRestart(int step)
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
void ScaTraImplicitTimeInt::SolveStationaryProblem()
{
  // time measurement: time loop (stationary) --- start TimeMonitor tm2
  tm2_ref_ = rcp(new TimeMonitor(*timetimeloop_));

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

     // other parameter needed by the elements
     eleparams.set("condif velocity field",cdvel_);
     eleparams.set("fs subgrid diffusivity",fssgd_);
     eleparams.set("using stationary formulation",true);

     // set vector values needed by elements
     discret_->ClearState();
     discret_->SetState("phinp",phinp_);
     // predicted dirichlet values
     // phinp then also holds prescribed new dirichlet values
     // dirichtoggle is 1 for dirichlet dofs, 0 elsewhere
     discret_->EvaluateDirichlet(eleparams,phinp_,null,null,dirichtoggle_);
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


  // end time measurement for time loop (stationary)
  tm2_ref_ = null;

  return;
} // ScaTraImplicitTimeInt::SolveStationaryProblem



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | update the velocity field                                   gjb 04/08|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void ScaTraImplicitTimeInt::SetVelocityField(int veltype, int velfuncno)
{
    if (veltype != cdvel_)
        dserror("velocity field type does not match: got %d, but expected %d!",veltype,cdvel_);

    if (veltype == 0) // zero
        convel_->PutScalar(0); // just to be sure!
    else if (veltype == 1)  // function
    {
    int numdim =3;
    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = discret_->lRowNode(lnodeid);
      for(int index=0;index<numdim;++index)
      {
        double value=DRT::UTILS::FunctionManager::Instance().Funct(velfuncno-1).Evaluate(index,lnode->X());
        convel_->ReplaceMyValue (lnodeid, index, value);
      }
    }
    }
    else
        dserror("error in setVelocityField");

    return;

} // ScaTraImplicitTimeInt::SetVelocityField


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | update the velocity field                                   gjb 04/08|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void ScaTraImplicitTimeInt::SetVelocityField(int veltype, Teuchos::RCP<const Epetra_Vector> extvel)
{
  if (veltype != cdvel_)
    dserror("velocity field type does not match: got %d, but expected %d!",veltype,cdvel_);

  // check vector compatibility and determine space dimension
  int numdim =-1;
  if (extvel->MyLength()== (2* convel_->MyLength()))
    numdim = 2;
  else if (extvel->MyLength()== (3* convel_->MyLength()))
    numdim = 3;
  else
    dserror("velocity vectors do not match in size");

  if ((numdim == 3) or (numdim == 2))
  {
    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      for(int index=0;index<numdim;++index)
      {
        double value = (*extvel)[lnodeid*numdim + index];
        //printf("myvelocityvalue[%d][%d] = %3.16lf\n",lnodeid,index,value);
        convel_->ReplaceMyValue(lnodeid, index, value);
      }
    }
  }

  return;

} // ScaTraImplicitTimeInt::SetVelocityField


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  set initial field for phi                                gjb   04/08|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void ScaTraImplicitTimeInt::SetInitialField(int init, int startfuncno)
{
  if (init == 0) // zero_field
  { // just to be sure!
    phinp_->PutScalar(0);
    phinm_->PutScalar(0); 
    phin_->PutScalar(0);
  }
  else if (init == 1)  // field_by_function
  {
    const Epetra_Map* dofrowmap = discret_->DofRowMap();

    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = discret_->lRowNode(lnodeid);
      // the set of degrees of freedom associated with the node
      vector<int> nodedofset = discret_->Dof(lnode);

      const int gid = nodedofset[0];
      int lid = dofrowmap->LID(gid);
      double phi0=DRT::UTILS::FunctionManager::Instance().Funct(startfuncno-1).Evaluate(0,lnode->X());

      phinp_->ReplaceMyValues(1,&phi0,&lid);
      phin_->ReplaceMyValues(1,&phi0,&lid);
      phinm_->ReplaceMyValues(1,&phi0,&lid);
    }
  }
  else if (init==2) // field_by_condition
  {
    dserror("Initialfield by condition not finished yet;");
    // access the initial field condition
    vector<DRT::Condition*> cond;
    discret_->GetCondition("InitialField", cond);

    const Epetra_Map* dofrowmap = discret_->DofRowMap();

    for (unsigned i=0; i<cond.size(); ++i)
    {
      cout<<"Applied InitialField Condition "<<i<<endl; 

      // loop all nodes on the processor
      for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
      {
        // get the processor local node
        DRT::Node*  lnode      = discret_->lRowNode(lnodeid);

        vector<DRT::Condition*> mycond;
        lnode->GetCondition("InitialField",mycond);

        if (mycond.size()>0)
        {
          // the set of degrees of freedom associated with the node
          vector<int> nodedofset = discret_->Dof(lnode);

          // get initial value from condition
          double phi0 = 2.0;

          // set initial value
          const int gid = nodedofset[0];
          int lid = dofrowmap->LID(gid);
          phinp_->ReplaceMyValues(1,&phi0,&lid);
          phin_->ReplaceMyValues(1,&phi0,&lid);
          phinm_->ReplaceMyValues(1,&phi0,&lid);
        }
      }
    }
  }
  else
    dserror("unknown option for condif initial field: %d", init);

  return;
} // ScaTraImplicitTimeInt::SetInitialField


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  calculate mass / heat flux vector                        gjb   04/08|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
Teuchos::RCP<Epetra_MultiVector> ScaTraImplicitTimeInt::CalcFlux()
{
  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  const Epetra_Map* dofrowmap = discret_->DofRowMap();
  // get the noderowmap
  const Epetra_Map* noderowmap = discret_->NodeRowMap();

  // empty vector for (normal) mass or heat fluxes
  Teuchos::RCP<Epetra_MultiVector> flux = rcp(new Epetra_MultiVector(*noderowmap,3,true));

  // we have only 1 dof per node, so we have to treat each spatial direction separately
  Teuchos::RCP<Epetra_Vector> fluxx = LINALG::CreateVector(*dofrowmap,true);
  Teuchos::RCP<Epetra_Vector> fluxy = LINALG::CreateVector(*dofrowmap,true);
  Teuchos::RCP<Epetra_Vector> fluxz = LINALG::CreateVector(*dofrowmap,true);

  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("phinp",phinp_);
  // set action for elements
  ParameterList eleparams;
  eleparams.set("action","calc_condif_flux");

  //provide velocity field (export to column map necessary for parallel evaluation)
  const Epetra_Map* nodecolmap = discret_->NodeColMap();
  RefCountPtr<Epetra_MultiVector> vel = rcp(new Epetra_MultiVector(*nodecolmap,3));
  LINALG::Export(*convel_,*vel);
  eleparams.set("velocity field",vel);

  // control parameters (not yet connected to input file)
  string fluxcomputation("domain"); // domain/condition
  string fluxtype("totalflux"); // calculate total/diffusive/convective flux vectors

  // visualization of total flux vector is default at the moment
  eleparams.set("fluxtype",fluxtype); // noflux / totalflux / diffusiveflux

  if (fluxcomputation=="domain")
  {
    // evaluate fluxes in the whole computational domain (e.g., for visualization of particle path-lines)
    discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,fluxx,fluxy,fluxz);
  }
  else if (fluxcomputation=="condition")
  {
    // evaluate fluxes on surface condition only
    // if restriction to normal(!) fluxes is needed put it here
    string condstring("FluxCalculation");
    discret_->EvaluateCondition(eleparams,Teuchos::null,Teuchos::null,fluxx,fluxy,fluxz,condstring);
  }
  else
    dserror("Unknown parameter for flux calculation.");

  // clean up
  discret_->ClearState();
  
  // insert values into final flux vector for visualization
  for (int i = 0;i<flux->MyLength();++i)
  {
    flux->ReplaceMyValue(i,0,(*fluxx)[i]);
    flux->ReplaceMyValue(i,1,(*fluxy)[i]);
    flux->ReplaceMyValue(i,2,(*fluxz)[i]);
  }

  return flux;
} // ScaTraImplicitTimeInt::CalcNormalFlux

/*----------------------------------------------------------------------*
 | Destructor dtor (public)                                     vg 05/07|
 *----------------------------------------------------------------------*/
ScaTraImplicitTimeInt::~ScaTraImplicitTimeInt()
{
  return;
}// ScaTraImplicitTimeInt::~ScaTraImplicitTimeInt::


#endif /* CCADISCRET       */
