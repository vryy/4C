/*!----------------------------------------------------------------------
\file
\brief Control routine for convection-diffusion time integration. Includes

     o Single step one-step-theta time integration

     o Two step BDF2 Gear's methode (with optional one-step-theta
                                     start step)



<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "condifimplicitintegration.H"
#include "../drt_lib/drt_nodematchingoctree.H"
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
                                             LINALG::Solver&       solver,
                                             ParameterList&        params,
                                             IO::DiscretizationWriter& output) :
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
  // get a vector layout from the discretization
  // -------------------------------------------------------------------
  {
    // Allocate integer vectors which will hold the number of the dofs
    vector<int> phimapdata;

    phimapdata.reserve(discret_->NumMyRowNodes());

    int countdofs = 0;
    for (int i=0; i<discret_->NumMyRowNodes(); ++i)
    {
      DRT::Node* node = discret_->lRowNode(i);
      vector<int> dof = discret_->Dof(node);

      int numdofs=discret_->Dof(node).size();
      if(numdofs==1)
      {
         // add this dof to the premapdata vector
        phimapdata.push_back(dof[numdofs-1]);
        countdofs += 1;
      }
      else
      {
        dserror("condif expects one dof");
      }
    }

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
  // the matrix. Assuming a structured 3d-condif mesh we have 27 adjacent
  // nodes with 1 dof each.
  // We do not need the exact number here, just for performance reasons
  // a 'good' estimate
  maxentriesperrow_ = 27;

  sysmat_ = null;

  // -------------------------------------------------------------------
  // create empty vectors
  // -------------------------------------------------------------------

  // Vectors passed to the element
  // -----------------------------
  // solutions at time n+1, n and n-1
  phinp_        = LINALG::CreateVector(*dofrowmap,true);
  phin_         = LINALG::CreateVector(*dofrowmap,true);
  phinm_        = LINALG::CreateVector(*dofrowmap,true);

  // Vectors associated to boundary conditions
  // -----------------------------------------

  // toggle vector indicating which dofs have Dirichlet BCs
  dirichtoggle_ = LINALG::CreateVector(*dofrowmap,true);

  // the vector containing body and surface forces
  neumann_loads_= LINALG::CreateVector(*dofrowmap,true);

  // The residual vector --- more or less the rhs
  residual_     = LINALG::CreateVector(*dofrowmap,true);

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

  // velocity field
  int cdvel  =params_.get<int>("condif velocity field");

  // time measurement --- this causes the TimeMonitor tm1 to stop here
  //                                                (call of destructor)
  tm1_ref_ = null;

  if (timealgo==timeint_stationary)
    // stationary case
    this->SolveStationaryProblem(cdvel);
  else  // instationary case
  {
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
      cdvel,
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
      cdvel,
      stepmax,
      maxtime,
      timealgo,
      theta
      );
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
void CondifImplicitTimeInt::TimeIntegrateFromTo(
  int&               step,
  double&            time,
  double&            dta,
  double&            dtp,
  int                cdvel,
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
     eleparams.set("total time",time);
     eleparams.set("theta",theta);
     eleparams.set("time step size",dta);
     eleparams.set("condif velocity field",cdvel);

     // set vector values needed by elements
     discret_->ClearState();
     discret_->SetState("phi at time n+1 (trial)",phinp_);
     discret_->SetState("phi at time n (trial)",phin_);
     // predicted dirichlet values
     // phinp then also holds prescribed new dirichlet values
     // dirichtoggle is 1 for dirichlet dofs, 0 elsewhere
     discret_->EvaluateDirichlet(eleparams,*phinp_,*dirichtoggle_);
     discret_->ClearState();

     // evaluate Neumann conditions
     eleparams.set("total time",time);
     eleparams.set("time constant for integration",theta);

     neumann_loads_->PutScalar(0.0);
     discret_->EvaluateNeumann(eleparams,*neumann_loads_);
     discret_->ClearState();
    }

    // -------------------------------------------------------------------
    //                     solve equation
    // -------------------------------------------------------------------
    this->Solve(dta,theta,cdvel,false);


    // -------------------------------------------------------------------
    //                         update solution
    //        current solution becomes old solution of next timestep
    //
    //  phinm_ =phin_
    //  phin_  =phinp_
    //
    // -------------------------------------------------------------------
    phinm_->Update(1.0,*phin_ ,0.0);
    phin_ ->Update(1.0,*phinp_,0.0);


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
} // CondifImplicitTimeInt::TimeIntegrateFromTo


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
  double dta,
  double theta,
  int cdvel,
  bool is_stat //if true, stationary formulations are used in the element
  )
{
  const Epetra_Map* dofrowmap       = discret_->DofRowMap();

  double            dtsolve;
  double            dtele  ;
  double            tcpu   ;

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
      eleparams.set("action","calc_condif_systemmat_and_residual");
      // choose what to assemble
      eleparams.set("assemble matrix 1",true);
      eleparams.set("assemble matrix 2",false);
      eleparams.set("assemble vector 1",true);
      eleparams.set("assemble vector 2",false);
      eleparams.set("assemble vector 3",false);
      // other parameters that might be needed by the elements
      eleparams.set("theta",theta);
      eleparams.set("time step size",dta);
      eleparams.set("condif velocity field",cdvel);
      eleparams.set("using stationary formulation",is_stat);
      // set vector values needed by elements
      discret_->ClearState();
      discret_->SetState("phi at time n+1 (trial)",phinp_);
      discret_->SetState("phi at time n (trial)",phin_);

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
    {
      // start time measurement for application of dirichlet conditions
      tm4_ref_ = rcp(new TimeMonitor(*timeapplydirich_));

      LINALG::ApplyDirichlettoSystem(sysmat_,phinp_,residual_,
				     phinp_,dirichtoggle_);

      // end time measurement for application of dirichlet conditions
      tm4_ref_=null;
    }

    //-------solve for residual displacements to correct incremental displacements
    {
      // start time measurement for element call
      tm5_ref_ = rcp(new TimeMonitor(*timesolver_));
      // get cpu time
      tcpu=ds_cputime();

      bool initsolver = true;
      solver_.Solve(sysmat_,phinp_,residual_,true,initsolver);

      // end time measurement for application of dirichlet conditions
      tm5_ref_=null;
      dtsolve=ds_cputime()-tcpu;
    }

  return;
} // CondifImplicitTimeInt::Solve

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | output of solution vector to binio                           vg 05/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void CondifImplicitTimeInt::Output(
  int    step,
  double time
  )
{

  //-------------------------------------------- output of solution
  output_.NewStep    (step,time);
  output_.WriteVector("phinp", phinp_);
  output_.WriteVector("residual", residual_);

  // do restart if we have to
  restartstep_ += 1;
  if (restartstep_ == uprestart_)
  {
    restartstep_ = 0;

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

	      vector<double> Values(10);
	      vector<int> Indices(10);

	      sysmat_->ExtractGlobalRowCopy(rr,
					    10,
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


#if 0  //if 0  //DEBUG IO --- incremental solution
      if (myrank_==0)
      {
        int rr;
        double* data = phinp_->Values();
        for(rr=0;rr<phinp_->MyLength();rr++)
        {
          printf("sol[%4d] %26.19e\n",rr,data[rr]);
        }
      }

#endif  //endif


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
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  set default parameter list (static/public)                  vg 05/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void CondifImplicitTimeInt::SetDefaults(ParameterList& params)
{
  // number of degrees of freedom
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

  return;
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
void CondifImplicitTimeInt::SolveStationaryProblem(int cdvel)
{
  // start time measurement for stationary solver
  tm2_ref_ = rcp(new TimeMonitor(*timedynloop_));

  // -------------------------------------------------------------------
  //                         out to screen
  // -------------------------------------------------------------------
  if (myrank_==0) printf("Stationary Solver");

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
     eleparams.set("condif velocity field",cdvel);
     eleparams.set("using stationary formulation",true);

     // set vector values needed by elements
     discret_->ClearState();
     discret_->SetState("phi at time n+1 (trial)",phinp_);
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

  // -------------------------------------------------------------------
  //                     solve nonlinear equation
  // -------------------------------------------------------------------
  this->Solve(1.0,1.0,cdvel,true);

  // -------------------------------------------------------------------
  //                         output of solution
  // -------------------------------------------------------------------
  this->Output(1,0.0);


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


#endif /* TRILINOS_PACKAGE */
#endif /* CCADISCRET       */
