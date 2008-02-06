/*!----------------------------------------------------------------------
\file condif_genalpha_integration.cpp
\brief Control routine for convection-diffusion instationary solver
       based on generalized-alpha time-integration scheme

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "condif_genalpha_integration.H"

extern Teuchos::RefCountPtr<Teuchos::ParameterList> globalparameterlist;

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Constructor (public)                                        vg 11/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
CondifGenAlphaIntegration::CondifGenAlphaIntegration(
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
  uprestart_(params.get("write restart every", -1)),
  writestep_(0),
  upres_(params.get("write solution every", -1))
{

  // -------------------------------------------------------------------
  // create timers and time monitor
  // -------------------------------------------------------------------
  timedyntot_     = TimeMonitor::getNewTimer("dynamic routine total"        );
  timedyninit_    = TimeMonitor::getNewTimer(" + initial phase"             );
  timedynloop_    = TimeMonitor::getNewTimer(" + time loop"                 );
  timepcloop_     = TimeMonitor::getNewTimer("   + predictor-corrector"     );
  timeeleloop_    = TimeMonitor::getNewTimer("      + element calls"        );
  timesolupdate_  = TimeMonitor::getNewTimer("      + update and calc. of intermediate sols");
  timeapplydirich_= TimeMonitor::getNewTimer("      + apply dirich cond."   );
  timeevaldirich_ = TimeMonitor::getNewTimer("      + evaluate dirich cond.");
  timesolver_     = TimeMonitor::getNewTimer("      + solver calls"         );
  timeout_        = TimeMonitor::getNewTimer("      + output and statistics");

  // time measurement --- start TimeMonitor tm0
  tm0_ref_        = rcp(new TimeMonitor(*timedyntot_ ));

  // time measurement --- start TimeMonitor tm7
  tm7_ref_        = rcp(new TimeMonitor(*timedyninit_ ));

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
  // initialize standard (stabilized) + discontinuity-capturing system matrix
  sysmat_dc_ = null;

  // -------------------------------------------------------------------
  // create empty vectors
  // -------------------------------------------------------------------

  // Vectors passed to the element
  // -----------------------------

  // temporal solution derivatives at time n and n-1
  phidtnp_      = LINALG::CreateVector(*dofrowmap,true);
  phidtn_       = LINALG::CreateVector(*dofrowmap,true);
  phidtam_      = LINALG::CreateVector(*dofrowmap,true);

  // solutions at time n+1, n and n-1
  phinp_        = LINALG::CreateVector(*dofrowmap,true);
  phin_         = LINALG::CreateVector(*dofrowmap,true);
  phiaf_        = LINALG::CreateVector(*dofrowmap,true);

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

  // The force vector (a copy of residual_ without Dirichlet
  // conditions applied)
  force_        = LINALG::CreateVector(*dofrowmap,true);

  // Increment vector
  increment_    = LINALG::CreateVector(*dofrowmap,true);


  // end time measurement for timeloop
  tm7_ref_ = null;

  return;

} // CondifGenAlphaIntegration::CondifGenAlphaIntegration


/*----------------------------------------------------------------------*
 | Destructor dtor (public)                                     vg 11/07|
 *----------------------------------------------------------------------*/
CondifGenAlphaIntegration::~CondifGenAlphaIntegration()
{
  return;
}// CondifGenAlphaIntegration::~CondifGenAlphaIntegration



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Contains the time loop for the generalized-alpha scheme              |
 |                                                              vg 11/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//

void CondifGenAlphaIntegration::GenAlphaIntegrateTo(
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

  if (myrank_==0)
  {
    cout << "Generalized Alpha parameter: alpha_F = " << alphaF_ << &endl;
    cout << "                             alpha_M = " << alphaM_ << &endl;
    cout << "                             gamma   = " << gamma_  << &endl <<&endl;

    cout <<  "                             " << "INERTIA         = " << ((*globalparameterlist).sublist("FluidStabilisation")).get<string>("INERTIA")        <<&endl;
    cout <<  "                             " << "SUPG            = " << ((*globalparameterlist).sublist("FluidStabilisation")).get<string>("SUPG")           <<&endl;
    cout <<  "                             " << "VSTAB           = " << ((*globalparameterlist).sublist("FluidStabilisation")).get<string>("VSTAB")          <<&endl;
    cout << &endl;

  }

  // start time measurement for timeloop
  tm2_ref_ = rcp(new TimeMonitor(*timedynloop_));

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
    //  predict new values for phi, using a constant predictor so far
    //  (Remark: For Dirichlet nodes, no matter what was set here,
    //   it will be overwritten by the prescribed value. The dphi/dt
    //   will calculated after the Dirichlet values will have been set.)
    // -------------------------------------------------------------------
    //       n+1      n
    //    phi    = phi
    //
    phinp_->Update(1.0,*phin_ ,0.0);

    // -------------------------------------------------------------------
    //         evaluate Dirichlet and Neumann boundary conditions
    // -------------------------------------------------------------------
    // start time measurement for application of Dirichlet conditions
    tm1_ref_ = rcp(new TimeMonitor(*timeevaldirich_));

    this->GenAlphaApplyDirichletAndNeumann();

    // end time measurement for application of Dirichlet conditions
    tm1_ref_=null;

    // -------------------------------------------------------------------
    //      calculate initial time-derivatives according to predicted
    //                  values of phi
    // -------------------------------------------------------------------
    this->GenAlphaCalcInitialTimeDeriv();

    // -------------------------------------------------------------------
    //     solve equation (linear problem: predictor + 1 corrector)
    // -------------------------------------------------------------------
    this->DoGenAlphaPredictorCorrector();

    // -------------------------------------------------------------------
    //                         update solution
    // -------------------------------------------------------------------
    this->GenAlphaTimeUpdate();


    // time measurement --- start TimeMonitor tm8
    tm8_ref_        = rcp(new TimeMonitor(*timeout_ ));

    // -------------------------------------------------------------------
    //                         output of solution
    // -------------------------------------------------------------------
    this->GenAlphaOutput();

    // time measurement --- stop TimeMonitor tm8
    tm8_ref_        = null;

    // -------------------------------------------------------------------
    //                    stop criterium for timeloop
    // -------------------------------------------------------------------

    if(step_>=endstep||time_>=endtime)
    {
	stop_timeloop=true;
    }

  }

  // end time measurement for timeloop
  tm2_ref_ = null;

  tm0_ref_ = null; // end total time measurement
  if(discret_->Comm().MyPID()==0)
  {
    cout<<endl<<endl;
  }
  TimeMonitor::summarize();


  return;
}// CondifGenAlphaIntegration::GenAlphaIntegrateFromTo


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Apply Dirichlet boundary conditions to solution vector              |
 |  Apply surface Neumann conditions (not yet implemented)              |
 |                                                              vg 11/07|
 -----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void CondifGenAlphaIntegration::GenAlphaApplyDirichletAndNeumann()
{
  // --------------------------------------------------
  // apply Dirichlet conditions to phinp

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
  eleparams.set("delta time",dt_  );
  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("phinp",phinp_);
  // predicted dirichlet values
  // phinp then also holds prescribed new dirichlet values
  // dirichtoggle is 1 for dirichlet dofs, 0 elsewhere
  discret_->EvaluateDirichlet(eleparams,*phinp_,*dirichtoggle_);
  discret_->ClearState();



  // --------------------------------------------------
  // evaluate Neumann conditions

  // not yet implemented

  return;
} // CondifGenAlphaIntegration::GenAlphaApplyDirichletAndNeumann


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  calculate dphi/dt according to prescribed Dirichlet values and      |
 |  predicted solution values.                                          |
 |                                                              vg 11/07|
 -----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void CondifGenAlphaIntegration::GenAlphaCalcInitialTimeDeriv()
{
  // ---------------------------------------------------------
  // adjust dphi/dt_np according to Dirichlet values of phi_np
  //
  //                                       n+1     n
  //        n+1          n  gamma-1.0   phi   - phi
  // dphi/dt    = dphi/dt * --------- + ------------
  //                          gamma      gamma * dt
  //

  phidtnp_->Update(1.0,*phinp_,-1.0,*phin_,0.0);
  phidtnp_->Update((gamma_-1.0)/gamma_,*phidtn_,1.0/(gamma_*dt_));

  return;
} // CondifGenAlphaIntegration::GenAlphaCalcInitialTimeDeriv


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Solve the linear problem resulting from the time-discrete version    |
 |                                                              vg 11/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//

void CondifGenAlphaIntegration::DoGenAlphaPredictorCorrector(
  )
{

  double            tcpu   ;

  dtsolve_ = 0;
  dtele_   = 0;


  // start time measurement for predictor-corrector
  tm6_ref_ = rcp(new TimeMonitor(*timepcloop_));

  // -------------------------------------------------------------------
  //  Evaluate dphi/dt and phi at the intermediate time level
  //                     n+alpha_M and n+alpha_F
  //
  //                             -> (0)
  // -------------------------------------------------------------------
  // start time measurement for nonlinear update
  tm9_ref_ = rcp(new TimeMonitor(*timesolupdate_));

  this->GenAlphaComputeIntermediateSol();

  // time measurement --- stop TimeMonitor tm9
  tm9_ref_        = null;

  // -------------------------------------------------------------------
  // call elements to calculate residual and matrix for first iteration
  // -------------------------------------------------------------------
  // start time measurement for element call
  tm3_ref_ = rcp(new TimeMonitor(*timeeleloop_));

  // get cpu time
  tcpu=ds_cputime();

  this->GenAlphaAssembleResidualAndMatrix();

  // end time measurement for element call
  tm3_ref_=null;
  dtele_=ds_cputime()-tcpu;

  // -------------------------------------------------------------------
  // solve for increments
  // -------------------------------------------------------------------
  // start time measurement for solver call
  tm5_ref_ = rcp(new TimeMonitor(*timesolver_));

  // get cpu time
  tcpu=ds_cputime();

  increment_->PutScalar(0.0);
  solver_.Solve(sysmat_,increment_,residual_,true,true);

  // end time measurement for solver call
  tm5_ref_=null;
  dtsolve_=ds_cputime()-tcpu;

  // start time measurement for nonlinear update
  tm9_ref_ = rcp(new TimeMonitor(*timesolupdate_));

  // -------------------------------------------------------------------
  // update estimates by incremental solution
  // -------------------------------------------------------------------
  this->GenAlphaSolUpdate();


  // end time measurement for nonlinear iteration
  tm6_ref_ = null;

  return;
}// CondifGenAlphaIntegration::DoGenAlphaPredictorCorrector


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Evaluate dphi/dt and phi at the intermediate time levels n+alpha_M    |
 | and n+alpha_F                                                        |
 |                                                              vg 11/07|
 -----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void CondifGenAlphaIntegration::GenAlphaComputeIntermediateSol()
{
  //         n+alphaM                    n+1                          n
  //  dphi/dt         = alpha_M * dphi/dt     + (1-alpha_M) *  dphi/dt

  phidtam_->Update((alphaM_),*phidtnp_,(1.0-alphaM_),*phidtn_,0.0);

  //       n+alphaF                n+1                     n
  //    phi         = alpha_F * phi     + (1-alpha_F) * phi

  phiaf_->Update((alphaF_),*phinp_,(1.0-alphaF_),*phin_,0.0);


  return;
} // CondifGenAlphaIntegration::GenAlphaComputeIntermediateSol


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Assemble residual and system matrix. Dirichlet conditions applied in |
 | here, the true residual is stored in force_.                         |
 |                                                              vg 11/07|
 -----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void CondifGenAlphaIntegration::GenAlphaAssembleResidualAndMatrix()
{
  const Epetra_Map* dofrowmap       = discret_->DofRowMap();

  // -------------------------------------------------------------------
  // call elements to calculate residual and matrix
  // -------------------------------------------------------------------
  // zero out the stiffness matrix
  // we keep the sparsity pattern throughout the calculation for
  // performance reasons
  //if (sysmat_==null)
    sysmat_ = LINALG::CreateMatrix(*dofrowmap,maxentriesperrow_);
  //else
  //sysmat_->PutScalar(0.0);

  // zero out residual
  residual_->PutScalar(0.0);

  // add Neumann loads to residual
  residual_->Update(1.0,*neumann_loads_,0.0);

  // create the parameters for the discretization
  ParameterList eleparams;

  // action for elements
  eleparams.set("action","calc_condif_genalpha_sysmat_and_residual");
  eleparams.set("assemble matrix 1",true);
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

  // parameters for stabilisation
  {
    eleparams.sublist("stabilisation") = (*globalparameterlist).sublist("FluidStabilisation");
  }

  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("phi (n+1      ,trial)", phinp_);
  discret_->SetState("phi (n+alpha_F,trial)", phiaf_);
  discret_->SetState("phidt(n+alpha_M,trial)",phidtam_);

  // call loop over elements
  {
    discret_->Evaluate(eleparams,sysmat_,null,residual_,null,null);
    discret_->ClearState();
  }

#if 0  // DEBUG IO --- the whole systemmatrix
      {
	  int rr;
	  int mm;
	  for(rr=0;rr<residual_->MyLength();rr++)
	  {
	      int NumEntries;

	      vector<double> Values(maxentriesperrow_);
	      vector<int> Indices(maxentriesperrow_);

	      sysmat_->ExtractGlobalRowCopy(rr,
					    maxentriesperrow_,
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
    const Epetra_Map* dofrowmap       = discret_->DofRowMap();

    int rr;


    double* data = residual_->Values();
    for(rr=0;rr<residual_->MyLength();rr++)
    {
      int  gid = dofrowmap->GID(rr);

      if(gid%4==0)
      printf("%4d rhs %22.15e\n",gid,data[rr]);
    }
  }

#endif

  // remember force vector for stress computation
  *force_=Epetra_Vector(*residual_);

  // finalize the system matrix
  LINALG::Complete(*sysmat_);
  maxentriesperrow_ = sysmat_->MaxNumEntries();

  // -------------------------------------------------------------------
  // Apply dirichlet boundary conditions to system of equations residual
  // phi values are supposed to be zero at boundary conditions
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

  return;
} // CondifGenAlphaIntegration::GenAlphaAssembleResidualAndMatrix


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | update the current acceleration, velocity and pressure    gammi 06/07|
 -----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void CondifGenAlphaIntegration::GenAlphaSolUpdate()
{
  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  // -------------------------------------------------------------------
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  int numlocdofs = dofrowmap->NumMyElements();

  int*   dof = dofrowmap->MyGlobalElements ();

  // loop all dofs on this proc
  for (int lid=0; lid<numlocdofs; ++lid)
  {
    int gid = dof[lid];

    // ------------------------------------------------------
    // update dphi/dt
    //
    //        n+1             n+1
    // dphi/dt      =  dphi/dt    + ddphi/dt
    //
    double dphidt = (*increment_)[lid];

    phidtnp_->SumIntoGlobalValues(1,&dphidt,&gid);

    // ------------------------------------------------------
    // use ddphi/dt to update phi:
    //
    //    n+1         n
    // phi      =  phi   +  gamma * dt * ddphi/dt
    //
      //
      double dphi = gamma_*dt_*dphidt;

      phinp_->SumIntoGlobalValues(1,&dphi,&gid);

  }

#if 0  // DEBUG IO  --- solution
  {
    const Epetra_Map* dofrowmap       = discret_->DofRowMap();

    int rr;


    double* data = phinp_->Values();
    for(rr=0;rr<phinp_->MyLength();rr++)
    {
      int  gid = dofrowmap->GID(rr);

      if(gid%4==0)
      printf("%4d phi %22.15e\n",gid,data[rr]);
    }
  }
#endif


  return;
} // CondifGenAlphaIntegration::GenAlphaSolUpdate


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Update: current solution becomes old solution of next timestep.     |
 |                                                              vg 11/07|
 -----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void CondifGenAlphaIntegration::GenAlphaTimeUpdate()
{
  //--------------------------------------------------
  // solution of this step becomes most recent solution of the last step

  // for phi
  phin_->Update(1.0,*phinp_ ,0.0);
  // for dphi/dt
  phidtn_->Update(1.0,*phidtnp_ ,0.0);

  return;
} // CondifGenAlphaIntegration::GenAlphaTimeUpdate


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Write solution to file                                   vg 06/07|
 -----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void CondifGenAlphaIntegration::GenAlphaOutput()
{
  //-------------------------------------------- output of solution
  restartstep_ += 1;
  writestep_ += 1;

  if (writestep_ == upres_)  //write solution
  {
    writestep_= 0;
    output_.NewStep    (step_,time_);

    output_.WriteVector("phinp"   , phinp_);

    // do restart if we have to
    if (restartstep_ == uprestart_)
    {
      restartstep_ = 0;

      output_.WriteVector("phin ",   phin_ );
      output_.WriteVector("phidtnp", phidtnp_);
      output_.WriteVector("phidtn ", phidtn_ );
    }
  }

  // write restart also when uprestart_ is not a integer multiple of upres_
  if ((restartstep_ == uprestart_) && (writestep_ > 0))
  {
    restartstep_ = 0;

    output_.NewStep    (step_,time_);

    output_.WriteVector("phinp",   phinp_);
    output_.WriteVector("phin ",   phin_ );
    output_.WriteVector("phidtnp", phidtnp_);
    output_.WriteVector("phidtn ", phidtn_ );
  }

  return;
} // CondifGenAlphaIntegration::GenAlphaOutput


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Read restart information                                     vg 11/07|
 -----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void CondifGenAlphaIntegration::ReadRestart(int step)
{
  IO::DiscretizationReader reader(discret_,step);
  time_ = reader.ReadDouble("time");
  step_ = reader.ReadInt("step");

  reader.ReadVector(phinp_,  "phinp");
  reader.ReadVector(phin_,   "phin");
  reader.ReadVector(phidtnp_,"phidtnp");
  reader.ReadVector(phidtn_, "phidtn");
}


#endif /* CCADISCRET       */
