/*!----------------------------------------------------------------------
\file fluid_genalpha_integration.cpp

\brief Time integration according to dis. C. Whiting




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

#include "fluid_genalpha_integration.H"

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
  IO::DiscretizationWriter&        output)
  :
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

  numdim_ = params_.get<int>("number of velocity degrees of freedom");

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
    // Allocate integer vectors which will hold the dof number of the
    // velocity or pressure dofs
    vector<int> velmapdata;
    vector<int> premapdata;

    velmapdata.reserve(discret_->NumMyRowNodes()*numdim_);
    premapdata.reserve(discret_->NumMyRowNodes());

    int countveldofs = 0;
    int countpredofs = 0;
    
    for (int i=0; i<discret_->NumMyRowNodes(); ++i)
    {
      DRT::Node* node = discret_->lRowNode(i);
      vector<int> dof = discret_->Dof(node);

      int numdofs=discret_->Dof(node).size();
      if(numdofs==numdim_+1)
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

  sysmat_ = null;

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
  increment_    = LINALG::CreateVector(*dofrowmap,true);


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

void FluidGenAlphaIntegration::GenAlphaIntegrateTo(
  int                endstep,
  double             endtime
  )
{

  bool stop_timeloop=false;

  dt_     = params_.get<double>("time step size");

  alphaM_ = params_.get<double>("alpha_M");
  alphaF_ = params_.get<double>("alpha_F");

  // choice of third parameter necessary but not sufficiant for second
  // order accuracy
  gamma_  = 0.5 + alphaM_ - alphaF_;

  
  cout << "Generalized Alpha parameter: alpha_F = " << alphaF_ << &endl;
  cout << "                             alpha_M = " << alphaM_ << &endl;
  cout << "                             gamma   = " << gamma_  << &endl <<&endl;

  
  while (stop_timeloop==false)
  {
    // -------------------------------------------------------------------
    //              set time dependent parameters
    // -------------------------------------------------------------------
    step_++;
    time_+=dt_;

    // -------------------------------------------------------------------
    //                         out to screen
    // -------------------------------------------------------------------
    if (myrank_==0)
    {
      printf("TIME: %11.4E/%11.4E  DT = %11.4E     GenAlpha     STEP = %4d/%4d \n",
             time_,endtime,dt_,step_,endstep);
 
    }

    // -------------------------------------------------------------------
    //     predict new values for velocity and pressure
    // -------------------------------------------------------------------
    this->GenAlphaPredictNewSolutionValues();

    // -------------------------------------------------------------------
    //         evaluate dirichlet and neumann boundary conditions
    // -------------------------------------------------------------------
    this->GenAlphaApplyDirichletAndNeumann();

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
    
    // -------------------------------------------------------------------
    //                         output of solution
    // -------------------------------------------------------------------
    this->GenAlphaOutput();

    // -------------------------------------------------------------------
    //                    stop criterium for timeloop
    // -------------------------------------------------------------------
    
    if(step_==endstep||time_>=endtime)
    {
	stop_timeloop=true;
    }

  }
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

  int               itnum         = 0;

  
  if(myrank_ == 0)
  {
    printf("+------------+-------------------+--------------+--------------+--------------+--------------+ \n");
    printf("|- step/max -|- tol      [norm] -|- vel-error --|- pre-error --|- vres-norm --|- pres-norm --| \n");
  }


  // -------------------------------------------------------------------
  //  Evaluate acceleration and velocity at the intermediate time level
  //                     n+alpha_M and n+alpha_F
  //
  //                             -> (0)
  // -------------------------------------------------------------------
  this->GenAlphaComputeIntermediateSol();
  
    
  // -------------------------------------------------------------------
  // call elements to calculate residual for first iteration
  // -------------------------------------------------------------------
  this->GenAlphaAssembleResidual();

  
  bool              stopnonliniter = false;
  while (stopnonliniter==false)
  {
    itnum++;

    // -------------------------------------------------------------------
    // call elements to calculate systemmatrix. Dirichlet conditions are
    // applied here to the system. 
    // -------------------------------------------------------------------
    this->GenAlphaAssembleMatrix();
    
  
    // -------------------------------------------------------------------
    // solve for increments
    // -------------------------------------------------------------------
    this->GenAlphaCalcIncrement(itnum);


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

    
    // -------------------------------------------------------------------
    // call elements to calculate residual for convergence check and next
    // step
    // -------------------------------------------------------------------
    this->GenAlphaAssembleResidual();

    
    // -------------------------------------------------------------------
    // do convergence check
    // -------------------------------------------------------------------
    stopnonliniter=this->GenAlphaNonlinearConvergenceCheck(itnum);
  }
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
  eleparams.set("delta time",dt_  );
  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("u and p at time n+1 (trial)",velnp_);
  // predicted dirichlet values
  // velnp then also holds prescribed new dirichlet values
  // dirichtoggle is 1 for dirichlet dofs, 0 elsewhere
  discret_->EvaluateDirichlet(eleparams,*velnp_,*dirichtoggle_);
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
  output_.NewStep    (step_,time_);
  output_.WriteVector("velnp"   , velnp_);
  output_.WriteVector("residual", residual_);
  
  // do restart if we have to
  restartstep_ += 1;
  if (restartstep_ == uprestart_)
  {
    restartstep_ = 0;

    output_.WriteVector("veln ", veln_ );
    output_.WriteVector("accnp", accnp_);
    output_.WriteVector("accn ", accn_ );
  }

  return;
} // FluidGenAlphaIntegration::GenAlphaOutput



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Assemble residual. Dirichlet conditions applied in here.  gammi 06/07|
 -----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FluidGenAlphaIntegration::GenAlphaAssembleResidual()
{
  // -------------------------------------------------------------------
  // call elements to calculate residual
  // -------------------------------------------------------------------
  // zero out residual
  residual_->PutScalar(0.0);

  // add Neumann loads to residual
  residual_->Update(1.0,*neumann_loads_,0.0);

  // create the parameters for the discretization
  ParameterList eleparams;
  
  // action for elements
  eleparams.set("action","calc_fluid_genalpha_residual");
  // choose what to assemble
  eleparams.set("assemble matrix 1",false);
  eleparams.set("assemble matrix 2",false);
  eleparams.set("assemble vector 1",true);
  eleparams.set("assemble vector 2",false);
  eleparams.set("assemble vector 3",false);
  // other parameters that might be needed by the elements
  eleparams.set("alpha_M",alphaM_);
  eleparams.set("alpha_F",alphaF_);
  eleparams.set("gamma"  ,gamma_ );
  eleparams.set("time"   ,time_  );
  eleparams.set("dt"     ,dt_    );
  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("u and p (n+1      ,trial)",velnp_);
  discret_->SetState("u and p (n+alpha_F,trial)",velaf_);
  discret_->SetState("acc     (n+alpha_M,trial)",accam_);

  // call loop over elements
  {
    discret_->Evaluate(eleparams,null,null,residual_,null,null);
    discret_->ClearState();
  }

  // -------------------------------------------------------------------
  // apply the Dirichlet boundary conditions to the residual
  // -------------------------------------------------------------------  
  LINALG::ApplyDirichlettoSystem(increment_,residual_,
                                 zeros_,dirichtoggle_);
  
  return;
} // FluidGenAlphaIntegration::GenAlphaAssembleResidual


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Assemble systemmatrix                                     gammi 06/07|
 -----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void  FluidGenAlphaIntegration::GenAlphaAssembleMatrix()
{
  const Epetra_Map* dofrowmap       = discret_->DofRowMap();

  
  // -------------------------------------------------------------------
  // call elements to calculate system matrix
  // -------------------------------------------------------------------
  // zero out the stiffness matrix
  sysmat_ = LINALG::CreateMatrix(*dofrowmap,maxentriesperrow_);

  // create the parameters for the discretization
  ParameterList eleparams;

  // action for elements
  eleparams.set("action","calc_fluid_genalpha_sysmat");
  // choose what to assemble
  eleparams.set("assemble matrix 1",true);
  eleparams.set("assemble matrix 2",false);
  eleparams.set("assemble vector 1",false);
  eleparams.set("assemble vector 2",false);
  eleparams.set("assemble vector 3",false);
  // other parameters that might be needed by the elements
  eleparams.set("alpha_M",alphaM_);
  eleparams.set("alpha_F",alphaF_);
  eleparams.set("gamma"  ,gamma_ );
  eleparams.set("time"   ,time_  );
  eleparams.set("dt"     ,dt_    );
  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("u and p (n+1      ,trial)",velnp_);
  discret_->SetState("u and p (n+alpha_F,trial)",velaf_);
  discret_->SetState("acc     (n+alpha_M,trial)",accam_);

  // call loop over elements
  {
    discret_->Evaluate(eleparams,sysmat_,null,null,null,null);
    discret_->ClearState();
  }

  // finalize the system matrix
  LINALG::Complete(*sysmat_);
  maxentriesperrow_ = sysmat_->MaxNumEntries();

  
  // -------------------------------------------------------------------
  // Apply dirichlet boundary conditions to system of equations residual
  // discplacements are supposed to be zero at boundary conditions
  // -------------------------------------------------------------------
  increment_->PutScalar(0.0);
  zeros_->PutScalar(0.0);
  {
    LINALG::ApplyDirichlettoSystem(sysmat_,increment_,residual_,
                                   zeros_,dirichtoggle_);
  }
return;
}//FluidGenAlphaIntegration::GenAlphaAssembleMatrix




//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Solve linear problem                                      gammi 06/07|
 -----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FluidGenAlphaIntegration::GenAlphaCalcIncrement(
  int itnum
  )
{
 
  //-------solve for residual displacements to correct incremental displacements
  bool initsolver = false;

  if (itnum==1) // init solver in first iteration only
  {
    initsolver = true;
  }
  solver_.Solve(sysmat_,increment_,residual_,true,initsolver);
  
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
      double dacc = (*increment_)[gid];
      
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
      
      double dpres = (*increment_)[gid];
      velnp_->SumIntoGlobalValues(1,&dpres,&gid);

    }
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
bool FluidGenAlphaIntegration::GenAlphaNonlinearConvergenceCheck(
  int          itnum  
  )
{
  bool stopnonliniter = false;

  RefCountPtr<Epetra_Vector> onlyvel_ = LINALG::CreateVector(*velrowmap_,true);
  RefCountPtr<Epetra_Vector> onlypre_ = LINALG::CreateVector(*prerowmap_,true);

  // ---------------------------------------------- nonlinear iteration
  // maximum number of nonlinear iteration steps
  int     itemax    =params_.get<int>   ("max nonlin iter steps");

  // ------------------------------- stop nonlinear iteration when both
  //                                 increment-norms are below this bound
  double  ittol     =params_.get<double>("tolerance for nonlin iter");


  // extract velocity and pressure increments from increment vector
  LINALG::Export(*increment_,*onlyvel_);
  LINALG::Export(*increment_,*onlypre_);
  // calculate L2_Norm of increments
  double incaccnorm_L2;
  double incprenorm_L2;
  onlyvel_->Norm2(&incaccnorm_L2);
  onlypre_->Norm2(&incprenorm_L2);

  double incvelnorm_L2=incaccnorm_L2*gamma_*dt_;

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

  // extract velocity and pressure residuals from rhs vector
  LINALG::Export(*residual_,*onlyvel_);
  LINALG::Export(*residual_,*onlypre_);
  
  double preresnorm_L2;
  onlypre_->Norm2(&preresnorm_L2);

  double velresnorm_L2;
  onlyvel_->Norm2(&velresnorm_L2);
  
  // out to screen
  if(myrank_ == 0)
  {
    printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   | %10.3E   | %10.3E   |\n",
           itnum,itemax,ittol,incvelnorm_L2/velnorm_L2,
           incprenorm_L2/prenorm_L2, velresnorm_L2,preresnorm_L2);
  }

  // this is the convergence check
  if((incvelnorm_L2/velnorm_L2 <= ittol && incprenorm_L2/prenorm_L2 <= ittol)
     ||
     (velresnorm_L2 <= ittol && preresnorm_L2 <= ittol))
  {
    stopnonliniter=true;
    if(myrank_ == 0)
    {
      printf("+------------+-------------------+--------------+--------------+--------------+--------------+ \n");
    }
  }

  // warn if itemax is reached without convergence, but proceed to
  // next timestep...
  if (itnum == itemax
      &&
      ((incvelnorm_L2/velnorm_L2 > ittol
        ||
        incprenorm_L2/prenorm_L2 > ittol
       )
       &&
       (
         velresnorm_L2 > ittol
         ||
         preresnorm_L2 > ittol
       )
      )
    )
  {
    stopnonliniter=true;
    if(myrank_ == 0)
    {
      printf("+---------------------------------------------------------------------------------------------+\n");
      printf("|            >>>>>> not converged in itemax steps!                                            |\n");
      printf("+---------------------------------------------------------------------------------------------+\n");
    }
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
        
        double initialval=DRT::Utils::FunctionManager::Instance().Funct(startfuncno-1).Evaluate(index,lnode->X());

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
        int err =0;
        double perc = 0.1;
        
        // loop all nodes on the processor
        for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
        {
          // get the processor local node
          DRT::Node*  lnode      = discret_->lRowNode(lnodeid);
          // the set of degrees of freedom associated with the node
          vector<int> nodedofset = discret_->Dof(lnode);

          double initialval=0;
          for(int index=0;index<numdim_;++index)
          {
            int gid = nodedofset[index];
        
            initialval=(*velnp_)[gid];
            if (initialval > EPS9)
            {
              break;
            }
          }
          
          for(int index=0;index<numdim_;++index)
          {
            int gid = nodedofset[index];
        
            double randomnumber = ((double)rand())/((double) RAND_MAX);

            double noise = initialval * randomnumber * perc;
            
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



#endif /* D_FLUID          */
#endif /* TRILINOS_PACKAGE */
#endif /* CCADISCRET       */
