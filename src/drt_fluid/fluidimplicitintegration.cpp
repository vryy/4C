/*!----------------------------------------------------------------------
\file fluidimplicitintegration.cpp
\brief Control routine for fluid (in)stationary solvers,

     including instationary solvers based on

     o one-step-theta time-integration scheme

     o two-step BDF2 time-integration scheme
       (with potential one-step-theta start algorithm)

     and stationary solver.

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <stdio.h>

#include "fluidimplicitintegration.H"
#include "time_integration_scheme.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/linalg_ana.H"
#include "../drt_lib/drt_nodematchingoctree.H"
#include "drt_periodicbc.H"
#include "../drt_lib/drt_function.H"
#include "fluid_utils.H"
#include "vm3_solver.H"
#include "fluidimpedancecondition.H"


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Constructor (public)                                     gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
FLD::FluidImplicitTimeInt::FluidImplicitTimeInt(RefCountPtr<DRT::Discretization> actdis,
                                                LINALG::Solver&       solver,
                                                ParameterList&        params,
                                                IO::DiscretizationWriter& output,
                                                bool alefluid)
  :
  // call constructor for "nontrivial" objects
  discret_(actdis),
  solver_ (solver),
  params_ (params),
  output_ (output),
  alefluid_(alefluid),
  time_(0.0),
  step_(0),
  extrapolationpredictor_(params.get("do explicit predictor",true)),
  uprestart_(params.get("write restart every", -1)),
  upres_(params.get("write solution every", -1)),
  writestresses_(params.get<int>("write stresses", 0)),
  freesurface_(NULL),
  fsisurface_(NULL),
  inrelaxation_(false)
{

  // -------------------------------------------------------------------
  // get the processor ID from the communicator
  // -------------------------------------------------------------------
  myrank_  = discret_->Comm().MyPID();

  // time measurement: initialization
  TEUCHOS_FUNC_TIME_MONITOR(" + initialization");

  // -------------------------------------------------------------------
  // get the basic parameters first
  // -------------------------------------------------------------------
  // type of solver: low-Mach-number or incompressible solver
  loma_ = params_.get<string>("low-Mach-number solver","No");
  // type of time-integration
  timealgo_ = params_.get<FLUID_TIMEINTTYPE>("time int algo");
  // time-step size
  dtp_ = dta_ = params_.get<double>("time step size");
  // maximum number of timesteps
  stepmax_  = params_.get<int>   ("max number timesteps");
  // maximum simulation time
  maxtime_  = params_.get<double>("total time");
  // parameter theta for time-integration schemes
  theta_    = params_.get<double>("theta");

  // parameter for linearization scheme (fixed-point-like or Newton)
  newton_ = params_.get<string>("Linearisation");

  // (fine-scale) subgrid viscosity?
  fssgv_ = params_.get<string>("fs subgrid viscosity","No");

  // -------------------------------------------------------------------
  // connect degrees of freedom for periodic boundary conditions
  // -------------------------------------------------------------------
  pbc_ = rcp(new PeriodicBoundaryConditions (discret_));
  pbc_->UpdateDofsForPeriodicBoundaryConditions();

  pbcmapmastertoslave_ = pbc_->ReturnAllCoupledNodesOnThisProc();

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

  const int numdim = params_.get<int>("number of velocity degrees of freedom");

  FLD::UTILS::SetupFluidSplit(*discret_,numdim,velpressplitter_);

  // -------------------------------------------------------------------
  // create empty system matrix --- stiffness and mass are assembled in
  // one system matrix!
  // -------------------------------------------------------------------

  // This is a first estimate for the number of non zeros in a row of
  // the matrix. Assuming a structured 3d-fluid mesh we have 27 adjacent
  // nodes with 4 dofs each. (27*4=108)
  // We do not need the exact number here, just for performance reasons
  // a 'good' estimate

  if (not params_.get<int>("Simple Preconditioner",0))
  {
    // initialize standard (stabilized) system matrix
    sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,108,false,true));
  }
  else
  {
    Teuchos::RCP<LINALG::BlockSparseMatrix<FLD::UTILS::VelPressSplitStrategy> > blocksysmat =
      Teuchos::rcp(new LINALG::BlockSparseMatrix<FLD::UTILS::VelPressSplitStrategy>(velpressplitter_,velpressplitter_,108,false,true));
    blocksysmat->SetNumdim(numdim);
    sysmat_ = blocksysmat;
  }

  // -------------------------------------------------------------------
  // create empty vectors
  // -------------------------------------------------------------------

  // additional rhs vector for robin-BC and vector for copying the residual
  robinrhs_      = LINALG::CreateVector(*dofrowmap,true);

  // Vectors passed to the element
  // -----------------------------
  // velocity/pressure at time n+1
  velnp_        = LINALG::CreateVector(*dofrowmap,true);

  // velocity/pressure at time n and n-1
  veln_         = LINALG::CreateVector(*dofrowmap,true);
  velnm_        = LINALG::CreateVector(*dofrowmap,true);

  // acceleration at time n and n-1
  accn_         = LINALG::CreateVector(*dofrowmap,true);
  accnm_        = LINALG::CreateVector(*dofrowmap,true);

  // velocity/density at time n+1
  vedenp_       = LINALG::CreateVector(*dofrowmap,true);

  if (loma_ != "No")
  {
    // velocity/density at time n and n-1
    veden_        = LINALG::CreateVector(*dofrowmap,true);
    vedenm_       = LINALG::CreateVector(*dofrowmap,true);
  }

  // history vector
  hist_           = LINALG::CreateVector(*dofrowmap,true);

  if (alefluid_)
  {
    dispnp_       = LINALG::CreateVector(*dofrowmap,true);
    dispn_        = LINALG::CreateVector(*dofrowmap,true);
    dispnm_       = LINALG::CreateVector(*dofrowmap,true);
    gridv_        = LINALG::CreateVector(*dofrowmap,true);
  }

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

  // right hand side vector for linearised solution;
  rhs_ = LINALG::CreateVector(*dofrowmap,true);

  // Nonlinear iteration increment vector
  incvel_       = LINALG::CreateVector(*dofrowmap,true);

  // initialise vectors and flags for (dynamic) Smagorinsky model
  // ------------------------------------------------------------
  //
  // (the smoothed quantities)
  //
  dynamic_smagorinsky_ = false;

  ParameterList *  modelparams =&(params_.sublist("TURBULENCE MODEL"));

  // flag for special flow: currently channel flow or flow in a lid-driven cavity
  special_flow_ = modelparams->get<string>("CANONICAL_FLOW","no");

  if (modelparams->get<string>("TURBULENCE_APPROACH","DNS_OR_RESVMM_LES")
      ==
      "CLASSICAL_LES")
  {
    // a canonical flow with homogeneous directions would allow a
    // spatial averaging of data
    string hom_plane = modelparams->get<string>("CHANNEL_HOMPLANE","not specified");

    if (myrank_ == 0)
    {

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
            cout << "      no homogeneous directions specified --- so we just use pointwise clipping for Cs\n";
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

    if(modelparams->get<string>("PHYSICAL_MODEL","no_model")
       ==
       "Dynamic_Smagorinsky"
      )
    {
      dynamic_smagorinsky_ = true;

      // get a vector layout from the discretization to construct

      const Epetra_Map* noderowmap    = discret_->NodeRowMap();

      filtered_vel_                   = rcp(new Epetra_MultiVector(*noderowmap,3,true));
      filtered_reynoldsstress_        = rcp(new Epetra_MultiVector(*noderowmap,9,true));
      filtered_modeled_subgrid_stress_= rcp(new Epetra_MultiVector(*noderowmap,9,true));

      averaged_LijMij_                = rcp(new vector<double>);
      averaged_MijMij_                = rcp(new vector<double>);
    }
  }

  // -------------------------------------------------------------------
  // initialize turbulence-statistics evaluation
  // -------------------------------------------------------------------
  //if (special_flow_ == "lid_driven_cavity")
  if (special_flow_ != "no")
  {
    // parameters for sampling/dumping period
    samstart_  = modelparams->get<int>("SAMPLING_START",1);
    samstop_   = modelparams->get<int>("SAMPLING_STOP",stepmax_);
    dumperiod_ = modelparams->get<int>("DUMPING_PERIOD",1);

    if (special_flow_ == "lid_driven_cavity")
      turbulencestatistics_ldc_=rcp(new TurbulenceStatisticsLdc(discret_,params_));
    else if (special_flow_ == "channel_flow_of_height_2")
      turbulencestatistics_=rcp(new TurbulenceStatistics(discret_,params_));
    else if (special_flow_ == "square_cylinder")
      turbulencestatistics_sqc_=rcp(new TurbulenceStatisticsSqc(discret_,params_));
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

  // -------------------------------------------------------------------
  // necessary only for the VM3 approach:
  // coarse- and fine-scale solution vectors + respective ouptput
  // -------------------------------------------------------------------
  if (fssgv_ != "No")
  {
    csvelnp_  = LINALG::CreateVector(*dofrowmap,true);
    fsvelnp_  = LINALG::CreateVector(*dofrowmap,true);
    convnp_   = LINALG::CreateVector(*dofrowmap,true);
    csconvnp_ = LINALG::CreateVector(*dofrowmap,true);
    fsconvnp_ = LINALG::CreateVector(*dofrowmap,true);

    if (myrank_ == 0)
    {
      // Output
      cout << "Fine-scale subgrid-viscosity approach based on AVM3: ";
      cout << &endl << &endl;
      cout << params_.get<string>("fs subgrid viscosity");

      if (fssgv_ == "Smagorinsky_all" || fssgv_ == "Smagorinsky_small" ||
          fssgv_ == "mixed_Smagorinsky_all" || fssgv_ == "mixed_Smagorinsky_small")
      {
        cout << " with Smagorinsky constant Cs= ";
        cout << modelparams->get<double>("C_SMAGORINSKY") ;
      }
      cout << &endl << &endl << &endl;
    }
  }

  // construct impedance bc wrapper
  impedancebc_ = rcp(new UTILS::FluidImpedanceWrapper(discret_, output_, dta_) );

  // get constant density variable for incompressible flow
  // set vedenp-vector values to 1.0 for incompressible flow
  // set density variable to 1.0 for low-Mach-number flow
  if (loma_ == "No")
  {
    ParameterList eleparams;
    eleparams.set("action","get_density");
    discret_->Evaluate(eleparams,null,null,null,null,null);
    density_ = eleparams.get("density", 1.0);
    if (density_ <= 0.0) dserror("received illegal density value");
    vedenp_->PutScalar(1.0);
  }
  else density_ = 1.0;

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
void FLD::FluidImplicitTimeInt::Integrate()
{
  // bound for the number of startsteps
  const int    numstasteps         =params_.get<int>   ("number of start steps");

  // output of stabilization details
  if (myrank_==0)
  {
    ParameterList *  stabparams=&(params_.sublist("STABILIZATION"));

    cout << "Stabilization type         : " << stabparams->get<string>("STABTYPE") << "\n";
    cout << "                             " << stabparams->get<string>("TDS")<< "\n";
    cout << "\n";

    if(stabparams->get<string>("TDS") == "quasistatic")
    {
      if(stabparams->get<string>("TRANSIENT")=="yes_transient")
      {
        dserror("The quasistatic version of the residual-based stabilization currently does not support the incorporation of the transient term.");
      }
    }
    cout <<  "                             " << "TRANSIENT       = " << stabparams->get<string>("TRANSIENT")      <<"\n";
    cout <<  "                             " << "SUPG            = " << stabparams->get<string>("SUPG")           <<"\n";
    cout <<  "                             " << "PSPG            = " << stabparams->get<string>("PSPG")           <<"\n";
    cout <<  "                             " << "VSTAB           = " << stabparams->get<string>("VSTAB")          <<"\n";
    cout <<  "                             " << "CSTAB           = " << stabparams->get<string>("CSTAB")          <<"\n";
    cout <<  "                             " << "CROSS-STRESS    = " << stabparams->get<string>("CROSS-STRESS")   <<"\n";
    cout <<  "                             " << "REYNOLDS-STRESS = " << stabparams->get<string>("REYNOLDS-STRESS")<<"\n";
    cout << "\n";
  }

  if (timealgo_==timeint_stationary)
    // stationary case
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

  // print the results of time measurements
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
void FLD::FluidImplicitTimeInt::TimeLoop()
{
  // time measurement: time loop
  TEUCHOS_FUNC_TIME_MONITOR(" + time loop");

  // how do we want to solve or fluid equations?
  const int dyntype    =params_.get<int>   ("type of nonlinear solve");

  if (dyntype==1)
  {
    if (alefluid_)
      dserror("no ALE possible with linearised fluid");
    if (fssgv_ != "No")
      dserror("no fine-scale solution implemented with linearised fluid");
    /* additionally it remains to mention that for the linearised
       fluid the stbilisation is hard coded to be SUPG/PSPG */
  }

  while (step_<stepmax_ and time_<maxtime_)
  {
    PrepareTimeStep();
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

    switch (dyntype)
    {
    case 0:
      // -----------------------------------------------------------------
      //                     solve nonlinear equation
      // -----------------------------------------------------------------
      NonlinearSolve();
      break;
    case 1:
      // -----------------------------------------------------------------
      //                     solve linearised equation
      // -----------------------------------------------------------------
      LinearSolve();
      break;
    default:
      dserror("Type of dynamics unknown!!");
    }

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
    // For low-Mach-number flow, the same is done for density values,
    // which are located at the "pressure dofs" of "vede"-vectors.
    // -------------------------------------------------------------------

    TimeUpdate();

    // time measurement: output and statistics
    TEUCHOS_FUNC_TIME_MONITOR("      + output and statistics");

    // -------------------------------------------------------------------
    //                    calculate lift'n'drag forces
    // -------------------------------------------------------------------
    const int liftdrag = params_.get<int>("liftdrag");

    if (liftdrag == 0); // do nothing, we don't want lift & drag
    if (liftdrag == 1)
      dserror("how did you manage to get here???");
    if (liftdrag == 2)
      LiftDrag();

    // -------------------------------------------------------------------
    // add calculated velocity to mean value calculation (statistics)
    // -------------------------------------------------------------------
    if(special_flow_ != "no" && step_>=samstart_ && step_<=samstop_)
    //if(special_flow_ == "lid_driven_cavity" && step_>=samstart_ && step_<=samstop_)
    {
      if(special_flow_ == "lid_driven_cavity")
        turbulencestatistics_ldc_->DoTimeSample(velnp_);
      else if(special_flow_ == "channel_flow_of_height_2")
        turbulencestatistics_->DoTimeSample(velnp_,*trueresidual_);
      else if(special_flow_ == "square_cylinder")
        turbulencestatistics_sqc_->DoTimeSample(velnp_);
    }

    // -------------------------------------------------------------------
    // evaluate error for test flows with analytical solutions
    // -------------------------------------------------------------------
    EvaluateErrorComparedToAnalyticalSol();

    // -------------------------------------------------------------------
    //                         output of solution
    // -------------------------------------------------------------------
    Output();

    // -------------------------------------------------------------------
    //                       update time step sizes
    // -------------------------------------------------------------------
    dtp_ = dta_;

    // -------------------------------------------------------------------
    //                    stop criterium for timeloop
    // -------------------------------------------------------------------
  }
} // FluidImplicitTimeInt::TimeLoop


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | setup the variables to do a new time step                 u.kue 06/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidImplicitTimeInt::PrepareTimeStep()
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
  // set part(s) of the rhs vector(s) belonging to the old timestep
  // for low-Mach-number flow: distinguish momentum and continuity part
  // (continuity part only meaningful for low-Mach-number flow)
  //
  // Stationary:
  //
  //               mom: hist_ = 0.0
  //              (con: hist_ = 0.0)
  //
  // One-step-Theta:
  //
  //               mom: hist_ = veln_  + dt*(1-Theta)*accn_
  //              (con: hist_ = densn_ + dt*(1-Theta)*densdtn_)
  //
  //
  // BDF2: for constant time step:
  //
  //               mom: hist_ = 4/3 veln_  - 1/3 velnm_
  //              (con: hist_ = 4/3 densn_ - 1/3 densnm_)
  //
  // -------------------------------------------------------------------
  if (loma_ != "No")
  {
    TIMEINT_THETA_BDF2::SetOldPartOfRighthandside(
      veden_, vedenm_, accn_,
      timealgo_, dta_, theta_,
      hist_);
  }
  else
  {
    TIMEINT_THETA_BDF2::SetOldPartOfRighthandside(
      veln_, velnm_, accn_,
      timealgo_, dta_, theta_,
      hist_);
  }

  // -------------------------------------------------------------------
  //                     do explicit predictor step
  //
  //                      +-                                      -+
  //                      | /     dta \          dta  veln_-velnm_ |
  // velnp_ = veln_ + dta | | 1 + --- | accn_ - ----- ------------ |
  //                      | \     dtp /          dtp     dtp       |
  //                      +-                                      -+
  //
  // -------------------------------------------------------------------
  //
  // We cannot have a predictor in case of monolithic FSI here. There needs to
  // be a way to turn this off.
  if (extrapolationpredictor_)
  {
    if (step_>1)
    {
      TIMEINT_THETA_BDF2::ExplicitPredictor(
          veln_, velnm_, accn_,
          dta_, dtp_,
          velnp_);
    }
  }

  // -------------------------------------------------------------------
  //         evaluate dirichlet and neumann boundary conditions
  // -------------------------------------------------------------------
  {
    ParameterList eleparams;

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
    discret_->EvaluateDirichlet(eleparams,velnp_,null,null,dirichtoggle_);
    discret_->ClearState();

    // evaluate Neumann conditions
    eleparams.set("total time",time_);
    eleparams.set("thsl",theta_*dta_);
    //eleparams.set("inc_density",density_);

    discret_->SetState("vedenp",vedenp_);
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
void FLD::FluidImplicitTimeInt::NonlinearSolve()
{
  inrelaxation_ = false;
  dirichletlines_ = Teuchos::null;
  // Do not remove meshmatrix_ here as we want to reuse its graph.
  // (We pay for the memory anyway if we use it, we might as well keep it.)
  //meshmatrix_ = Teuchos::null;

  // time measurement: nonlinear iteration
  TEUCHOS_FUNC_TIME_MONITOR("   + nonlin. iteration/lin. solve");

  // ---------------------------------------------- nonlinear iteration
  // ------------------------------- stop nonlinear iteration when both
  //                                 increment-norms are below this bound
  const double  ittol     =params_.get<double>("tolerance for nonlin iter");

  //------------------------------ turn adaptive solver tolerance on/off
  const bool   isadapttol    = params_.get<bool>("ADAPTCONV",true);
  const double adaptolbetter = params_.get<double>("ADAPTCONV_BETTER",0.01);

  const bool fluidrobin = params_.get<bool>("fluidrobin", false);

  int               itnum = 0;
  int               itemax = 0;
  bool              stopnonliniter = false;

  // currently default for turbulent channel flow: only one iteration before sampling
  if (special_flow_ == "channel_flow_of_height_2" && step_<samstart_ )
    itemax  = 2;
  else itemax  = params_.get<int>   ("max nonlin iter steps");

  double dtsolve = 0.0;
  double dtele   = 0.0;
  double dtfilter = 0.0;

  if (myrank_ == 0)
  {
    printf("+------------+-------------------+--------------+--------------+--------------+--------------+\n");
    printf("|- step/max -|- tol      [norm] -|-- vel-res ---|-- pre-res ---|-- vel-inc ---|-- pre-inc ---|\n");
  }

  while (stopnonliniter==false)
  {
    itnum++;

    // -------------------------------------------------------------------
    // call elements to calculate system matrix
    // -------------------------------------------------------------------
    {
      // time measurement: element
      TEUCHOS_FUNC_TIME_MONITOR("      + element calls");

      // get cpu time
      const double tcpu=ds_cputime();

      sysmat_->Zero();

      // add Neumann loads
      residual_->Update(1.0,*neumann_loads_,0.0);

      // update impedance boundary condition
      impedancebc_->UpdateResidual(residual_);

      // create the parameters for the discretization
      ParameterList eleparams;

      // add stabilization term at Neumann outflow boundary if required
      if(outflow_stab_ == "yes_outstab")
      {
        discret_->ClearState();
        discret_->SetState("velnp",velnp_);
        eleparams.set("thsl",theta_*dta_);
        //eleparams.set("inc_density",density_);
        eleparams.set("outflow stabilization",outflow_stab_);

        discret_->SetState("vedenp",vedenp_);
        outflow_stabil_->PutScalar(0.0);
        discret_->EvaluateNeumann(eleparams,*outflow_stabil_);
        discret_->ClearState();

        // add Neumann-type stabilization term to residual vector
        residual_->Update(1.0,*outflow_stabil_,1.0);
      }

      // Filter velocity for dynamic Smagorinsky model --- this provides
      // the necessary dynamic constant
      // //
      // convergence check at itemax is skipped for speedup if
      // CONVCHECK is set to L_2_norm_without_residual_at_itemax
      if ((itnum != itemax)
          ||
          (params_.get<string>("CONVCHECK","L_2_norm")
           !=
           "L_2_norm_without_residual_at_itemax"))
      {
        if (dynamic_smagorinsky_)
        {
          // time measurement
          const double tcpufilter=ds_cputime();
          this->ApplyFilterForDynamicComputationOfCs();
          dtfilter=ds_cputime()-tcpufilter;
        }
      }

      // compute convective stresses when scale-similarity models are used
      if (fssgv_ == "scale_similarity" ||
          fssgv_ == "mixed_Smagorinsky_all" ||
          fssgv_ == "mixed_Smagorinsky_small")
      {
        // action for elements
        eleparams.set("action","calc_convective_stresses");

        // set vector values needed by elements
        discret_->ClearState();
        discret_->SetState("velnp",velnp_);
        discret_->SetState("vedenp",vedenp_);

        // element evaluation for getting convective stresses at nodes
        discret_->Evaluate(eleparams,null,convnp_);
        discret_->ClearState();
      }

      // action for elements
      if (timealgo_==timeint_stationary)
        eleparams.set("action","calc_fluid_stationary_systemmat_and_residual");
      else
        eleparams.set("action","calc_fluid_systemmat_and_residual");

      // other parameters that might be needed by the elements
      eleparams.set("total time",time_);
      eleparams.set("thsl",theta_*dta_);
      eleparams.set("dt",dta_);
      eleparams.set("fs subgrid viscosity",fssgv_);
      eleparams.set("Linearisation",newton_);
      eleparams.set("low-Mach-number solver",loma_);

      // parameters for stabilization
      eleparams.sublist("STABILIZATION") = params_.sublist("STABILIZATION");

      // parameters for stabilization
      eleparams.sublist("TURBULENCE MODEL") = params_.sublist("TURBULENCE MODEL");

      // set vector values needed by elements
      discret_->ClearState();
      discret_->SetState("velnp",velnp_);
      discret_->SetState("vedenp",vedenp_);

      discret_->SetState("hist"  ,hist_ );
      if (alefluid_)
      {
        discret_->SetState("dispnp", dispnp_);
        discret_->SetState("gridv", gridv_);
      }

      //----------------------------------------------------------------------
      // decide whether VM3-based solution approach or standard approach
      //----------------------------------------------------------------------
      if (fssgv_ != "No")
      {
        // time measurement: avm3
        TEUCHOS_FUNC_TIME_MONITOR("           + avm3");

        // call the VM3 constructor (only in the first time step)
        if (step_ == 1)
        {
          // zero fine-scale vector
          fsvelnp_->PutScalar(0.0);

          // set coarse- and fine-scale vectors
          discret_->SetState("fsvelnp",fsvelnp_);
          discret_->SetState("csvelnp",csvelnp_);
          discret_->SetState("csconvnp",csconvnp_);

          // element evaluation for getting system matrix
          discret_->Evaluate(eleparams,sysmat_,residual_);


	  // complete system matrix
          sysmat_->Complete();

          // apply DBC to system matrix
          LINALG::ApplyDirichlettoSystem(sysmat_,incvel_,residual_,zeros_,dirichtoggle_);

          // call VM3 constructor with system matrix
          // extract the ML parameters
          ParameterList&  mllist = solver_.Params().sublist("ML Parameters");
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
        vm3_solver_->Separate(csvelnp_,fsvelnp_,velnp_);

        // call VM3 scale separation to get coarse-scale part of convective stresses
        if (fssgv_ == "scale_similarity" ||
            fssgv_ == "mixed_Smagorinsky_all" || fssgv_ == "mixed_Smagorinsky_small")
        {
          vm3_solver_->Separate(csconvnp_,fsconvnp_,convnp_);
          discret_->SetState("csvelnp",csvelnp_);
          discret_->SetState("csconvnp",csconvnp_);
        }

        // set coarse- and fine-scale vectors
        discret_->SetState("fsvelnp",fsvelnp_);
      }

      // convergence check at itemax is skipped for speedup if
      // CONVCHECK is set to L_2_norm_without_residual_at_itemax
      if ((itnum != itemax)
          ||
          (params_.get<string>("CONVCHECK","L_2_norm")
           !=
           "L_2_norm_without_residual_at_itemax"))
      {
        // call standard loop over elements
        discret_->Evaluate(eleparams,sysmat_,residual_);

        discret_->ClearState();

#if 0

        //---------------------------surface tension update
        if (alefluid_ and freesurface_->Relevant())
        {
          // select free surface elements
          std::string condname = "FREESURFCoupling";

          ParameterList eleparams;
          // set action for elements, calc_surface_tension uses node normals
          eleparams.set("action","calc_node_normal");

          // get a vector layout from the discretization to construct matching
          // vectors and matrices
          //                 local <-> global dof numbering
          const Epetra_Map* dofrowmap = discret_->DofRowMap();

          //vector ndnorm0 with pressure-entries is needed for EvaluateCondition
          Teuchos::RCP<Epetra_Vector> ndnorm0 = LINALG::CreateVector(*dofrowmap,true);

          //call loop over elements, note: normal vectors do not yet have length = 1.0
          discret_->ClearState();
          discret_->SetState("dispnp", dispnp_);
          discret_->EvaluateCondition(eleparams,ndnorm0,condname);

          // set action for elements
          eleparams.set("action","calc_surface_tension");
          eleparams.set("thsl",theta_*dta_);
          eleparams.set("dta",dta_);
          discret_->SetState("normals", ndnorm0);
          discret_->EvaluateCondition(eleparams,sysmat_,Teuchos::null,residual_,Teuchos::null,Teuchos::null,condname);
          discret_->ClearState();
        }
        //---------------------------end of surface tension update

#endif

        trueresidual_->Update(ResidualScaling(),*residual_,0.0);

        // finalize the complete matrix
        sysmat_->Complete();

        // If we have a robin condition we need to modify both the rhs and the
        // matrix diagonal corresponding to the dofs at the robin interface.
        if (fluidrobin)
        {
          // Add structral part of Robin force
          // (combination of structral force and velocity)
          residual_->Update(theta_*dta_,*robinrhs_,1.0);

          double alphaf = params_.get<double>("alpharobinf",-1.0);
          double scale = alphaf*theta_*dta_;

          // Add fluid part of Robin force
          // (scaled fluid velocity)
          fsisurface_->AddCondVector(-1.*scale,
                                     fsisurface_->ExtractCondVector(velnp_),
                                     residual_);

          // Note: It is the right thing to test the robin enhanced residual_
          // for convergence, since the velocity terms are to vanish and the
          // structural forces are to cancel with the internal forces.
          //
          // Note: We do not add any external (robin) loads to
          // trueresidual_. This way we get the unbalanced forces at the
          // interface, which can be applied to the structure later on.

          const Epetra_Map& robinmap = *fsisurface_->CondMap();
          int numrdofs = robinmap.NumMyElements();
          int* rdofs = robinmap.MyGlobalElements();
          for (int lid=0; lid<numrdofs; ++lid)
          {
            int gid = rdofs[lid];
            // We assemble with a global id into a filled matrix here. This is
            // fine as we know we do not add new entries but just add to the
            // diagonal.
            //
            // Note: The matrix lives in the full fluid map whereas
            // we loop the robin interface map here, so our local ids are very
            // different from the matrix local ids.
            //
            // Note: This assemble might fail if we have a block matrix here.
            // (No, it won't since the matrix is already filled. :] )
            sysmat_->Assemble(scale,gid,gid);
          }
        }
      }

      // end time measurement for element
      dtele=ds_cputime()-tcpu;
    }

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

    Teuchos::RCP<Epetra_Vector> onlyvel = velpressplitter_.ExtractOtherVector(residual_);
    onlyvel->Norm2(&vresnorm);

    velpressplitter_.ExtractOtherVector(incvel_,onlyvel);
    onlyvel->Norm2(&incvelnorm_L2);

    velpressplitter_.ExtractOtherVector(velnp_,onlyvel);
    onlyvel->Norm2(&velnorm_L2);

    Teuchos::RCP<Epetra_Vector> onlypre = velpressplitter_.ExtractCondVector(residual_);
    onlypre->Norm2(&presnorm);

    velpressplitter_.ExtractCondVector(incvel_,onlypre);
    onlypre->Norm2(&incprenorm_L2);

    velpressplitter_.ExtractCondVector(velnp_,onlypre);
    onlypre->Norm2(&prenorm_L2);

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
        printf(" (      --     ,te=%10.3E",dtele);
        if (dynamic_smagorinsky_)
        {
          printf(",tf=%10.3E",dtfilter);
        }
        printf(")\n");
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
          printf(" (ts=%10.3E,te=%10.3E",dtsolve,dtele);
          if (dynamic_smagorinsky_)
          {
            printf(",tf=%10.3E",dtfilter);
          }
          printf(")\n");
          printf("+------------+-------------------+--------------+--------------+--------------+--------------+\n");

          FILE* errfile = params_.get<FILE*>("err file",NULL);
          if (errfile!=NULL)
          {
            fprintf(errfile,"fluid solve:   %3d/%3d  tol=%10.3E[L_2 ]  vres=%10.3E  pres=%10.3E  vinc=%10.3E  pinc=%10.3E\n",
                    itnum,itemax,ittol,vresnorm,presnorm,
                    incvelnorm_L2/velnorm_L2,incprenorm_L2/prenorm_L2);
          }
        }
        break;
      }
      else // if not yet converged
        if (myrank_ == 0)
        {
          printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   | %10.3E   | %10.3E   |",
                 itnum,itemax,ittol,vresnorm,presnorm,
                 incvelnorm_L2/velnorm_L2,incprenorm_L2/prenorm_L2);
          printf(" (ts=%10.3E,te=%10.3E",dtsolve,dtele);
          if (dynamic_smagorinsky_)
          {
            printf(",tf=%10.3E",dtfilter);
          }
          printf(")\n");
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

        FILE* errfile = params_.get<FILE*>("err file",NULL);
        if (errfile!=NULL)
        {
          fprintf(errfile,"fluid unconverged solve:   %3d/%3d  tol=%10.3E[L_2 ]  vres=%10.3E  pres=%10.3E  vinc=%10.3E  pinc=%10.3E\n",
                  itnum,itemax,ittol,vresnorm,presnorm,
                  incvelnorm_L2/velnorm_L2,incprenorm_L2/prenorm_L2);
        }
      }
      break;
    }

    //--------- Apply dirichlet boundary conditions to system of equations
    //          residual discplacements are supposed to be zero at
    //          boundary conditions
    incvel_->PutScalar(0.0);
    {
      // time measurement: application of dbc
      TEUCHOS_FUNC_TIME_MONITOR("      + apply DBC");
      LINALG::ApplyDirichlettoSystem(sysmat_,incvel_,residual_,zeros_,dirichtoggle_);
    }

    //-------solve for residual displacements to correct incremental displacements
    {
      // time measurement: solver
      TEUCHOS_FUNC_TIME_MONITOR("      + solver calls");

      // get cpu time
      const double tcpusolve=ds_cputime();

      // do adaptive linear solver tolerance (not in first solve)
      if (isadapttol && itnum>1)
      {
        double currresidual = max(vresnorm,presnorm);
        currresidual = max(currresidual,incvelnorm_L2/velnorm_L2);
        currresidual = max(currresidual,incprenorm_L2/prenorm_L2);
        solver_.AdaptTolerance(ittol,currresidual,adaptolbetter);
      }
      solver_.Solve(sysmat_->EpetraOperator(),incvel_,residual_,true,itnum==1);
      solver_.ResetTolerance();

      // end time measurement for solver
      dtsolve = ds_cputime()-tcpusolve;
    }

    //------------------------------------------------ update (u,p) trial
    velnp_->Update(1.0,*incvel_,1.0);

    // free surface update
    if (alefluid_ and freesurface_->Relevant())
    {
      Teuchos::RefCountPtr<Epetra_Vector> fsvelnp = freesurface_->ExtractCondVector(velnp_);
      Teuchos::RefCountPtr<Epetra_Vector> fsdisp = freesurface_->ExtractCondVector(dispn_);
      Teuchos::RefCountPtr<Epetra_Vector> fsdispnp = Teuchos::rcp(new Epetra_Vector(*freesurface_->CondMap()));

      // select free surface elements
      std::string condname = "FREESURFCoupling";

      std::vector<DRT::Condition*> conds;
      discret_->GetCondition(condname, conds);

      // select only heightfunction conditions here
      std::vector<DRT::Condition*> hfconds;
      for (unsigned i=0; i<conds.size(); ++i)
      {
        if (*conds[i]->Get<std::string>("coupling")=="heightfunction")
          hfconds.push_back(conds[i]);
      }

      conds.clear();

      if (hfconds.size()>0)
      {
        ParameterList eleparams;
        // set action for elements
        eleparams.set("action","calc_node_normal");

        // get a vector layout from the discretization to construct matching
        // vectors and matrices
        //                 local <-> global dof numbering
        const Epetra_Map* dofrowmap = discret_->DofRowMap();

        //vector ndnorm0 with pressure-entries is needed for EvaluateCondition
        Teuchos::RCP<Epetra_Vector> ndnorm0 = LINALG::CreateVector(*dofrowmap,true);

        //call loop over elements, note: normal vectors do not yet have length = 1.0
        discret_->ClearState();
        discret_->SetState("dispnp", dispnp_);
        discret_->EvaluateCondition(eleparams,ndnorm0,condname);
        discret_->ClearState();

        //ndnorm contains fsnodes' normal vectors (with arbitrary length). no pressure-entries
        Teuchos::RefCountPtr<Epetra_Vector> ndnorm = freesurface_->ExtractCondVector(ndnorm0);

        std::vector< int > GIDfsnodes;  //GIDs of free surface nodes
        std::vector< int > GIDdof;      //GIDs of current fsnode's dofs
        std::vector< int > rfs;         //local indices for ndnorm and fsvelnp for current node

        //get GIDs of free surface nodes for this processor
        DRT::UTILS::FindConditionedNodes(*discret_,hfconds,GIDfsnodes);

        for (unsigned int node=0; node<(GIDfsnodes.size()); node++)
        {
          //get LID for this node
          int ndLID = (discret_->NodeRowMap())->LID(GIDfsnodes[node]);
          if (ndLID == -1) dserror("No LID for free surface node");

          //get vector of this node's dof GIDs
          GIDdof.clear();
          discret_->Dof(discret_->lRowNode(ndLID), GIDdof);
          GIDdof.pop_back();  //free surface nodes: pop pressure dof

          //numdof = dim, no pressure
          int numdof = GIDdof.size();
          rfs.clear();
          rfs.resize(numdof);

          //get local indices for dofs in ndnorm and fsvelnp
          for (int i=0; i<numdof; i++)
          {
            int rgid = GIDdof[i];
            if (!ndnorm->Map().MyGID(rgid) or !fsvelnp->Map().MyGID(rgid))
              dserror("Sparse vector does not have global row %d",rgid);
            rfs[i] = ndnorm->Map().LID(rgid);
          }

          double length = 0.0;
          for (int i=0; i<numdof; i++)
            length += (*ndnorm)[rfs[i]] * (*ndnorm)[rfs[i]];
          length = sqrt(length);

          double pointproduct = 0.0;
          for (int i=0; i<numdof; i++)
          {
            (*ndnorm)[rfs[i]] = (1.0/length) * (*ndnorm)[rfs[i]];
            //height function approach: left hand side of eq. 15
            pointproduct += (*ndnorm)[rfs[i]] * (*fsvelnp)[rfs[i]];
          }

          for (int i=0; i<numdof; i++)
          {
            //height function approach: last entry of u_G is delta_phi/delta_t,
            //the other entries are zero
            if (i == numdof-1)
              (*fsvelnp)[rfs[i]]  = pointproduct / (*ndnorm)[rfs[i]];
            else
              (*fsvelnp)[rfs[i]] = 0.0;
          }
        }
      }

      fsdispnp->Update(1.0,*fsdisp,dta_,*fsvelnp,0.0);

      freesurface_->InsertCondVector(fsdispnp,dispnp_);
      freesurface_->InsertCondVector(fsvelnp,gridv_);
    }
  }
} // FluidImplicitTimeInt::NonlinearSolve


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | the time step of a linearised fluid                      chfoe 02/08 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*
This fluid implementation is designed to be quick(er) but has a couple of
drawbacks:
- currently it is incapable of ALE fluid solutions
- the order of accuracy in time is fixed to 1, i.e. some more steps may be required
- some effort has to be made if correct nodal forces are required as this
  implementation does a total solve rather than an incremental one.
*/
void FLD::FluidImplicitTimeInt::LinearSolve()
{
  // time measurement: linearised fluid
  TEUCHOS_FUNC_TIME_MONITOR("   + nonlin. iteration/lin. solve");

  if (myrank_ == 0)
    cout << "solution of linearised fluid   ";

  // -------------------------------------------------------------------
  // call elements to calculate system matrix
  // -------------------------------------------------------------------

  // get cpu time
  const double tcpuele = ds_cputime();
  {
    // time measurement: element
    TEUCHOS_FUNC_TIME_MONITOR("      + element calls");

    sysmat_->Zero();

    // add Neumann loads
    rhs_->Update(1.0,*neumann_loads_,0.0);

    // create the parameters for the discretization
    ParameterList eleparams;

    // action for elements
    if (timealgo_==timeint_stationary)
      dserror("no stationary solution with linearised fluid!!!");
    else
      eleparams.set("action","calc_linear_fluid");

    // other parameters that might be needed by the elements
    eleparams.set("total time",time_);
    eleparams.set("thsl",theta_*dta_);
    eleparams.set("low-Mach-number solver",loma_);

    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("velnp",velnp_);
    discret_->SetState("vedenp",vedenp_);
    discret_->SetState("hist"  ,hist_ );

    // call standard loop over linear elements
    discret_->Evaluate(eleparams,sysmat_,rhs_);
    discret_->ClearState();

    // finalize the complete matrix
    sysmat_->Complete();
  }
  // end time measurement for element
  const double dtele = ds_cputime() - tcpuele;

  //--------- Apply dirichlet boundary conditions to system of equations
  //          residual velocities (and pressures) are supposed to be zero at
  //          boundary conditions
  //velnp_->PutScalar(0.0);

  {
    // time measurement: application of dbc
    TEUCHOS_FUNC_TIME_MONITOR("      + apply DBC");

    LINALG::ApplyDirichlettoSystem(sysmat_,velnp_,rhs_,velnp_,dirichtoggle_);
  }

  //-------solve for total new velocities and pressures
  // get cpu time
  const double tcpusolve = ds_cputime();
  {
    // time measurement: solver
    TEUCHOS_FUNC_TIME_MONITOR("      + solver calls");

    /* possibly we could accelerate it if the reset variable
       is true only every fifth step, i.e. set the last argument to false
       for 4 of 5 timesteps or so. */
    solver_.Solve(sysmat_->EpetraOperator(),velnp_,rhs_,true,true);
  }
  // end time measurement for solver
  const double dtsolve = ds_cputime() - tcpusolve;

  if (myrank_ == 0)
    cout << "te=" << dtele << ", ts=" << dtsolve << "\n\n" ;

} // FluidImplicitTimeInt::LinearSolve


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | build linear system matrix and rhs                        u.kue 11/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidImplicitTimeInt::Evaluate(Teuchos::RCP<const Epetra_Vector> vel)
{
  sysmat_->Zero();
  if (meshmovematrix_ != Teuchos::null)
    meshmovematrix_->Zero();

  // set the new solution we just got
  if (vel!=Teuchos::null)
  {
    const int len = vel->MyLength();

    // Take Dirichlet values from velnp and add vel to veln for non-Dirichlet
    // values.
    //
    // There is no epetra operation for this! Maybe we could have such a beast
    // in ANA?

    double* veln  = &(*veln_)[0];
    double* velnp = &(*velnp_)[0];
    double* dt    = &(*dirichtoggle_)[0];
    double* idv   = &(*invtoggle_)[0];
    const double* incvel = &(*vel)[0];

    //------------------------------------------------ update (u,p) trial
    for (int i=0; i<len; ++i)
    {
      velnp[i] = velnp[i]*dt[i] + (veln[i] + incvel[i])*idv[i];
    }
  }

  // add Neumann loads
  residual_->Update(1.0,*neumann_loads_,0.0);

  // update impedance boundary condition
  impedancebc_->UpdateResidual(residual_);

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
  eleparams.set("dt",dta_);
  eleparams.set("fs subgrid viscosity",fssgv_);
  eleparams.set("Linearisation",newton_);
  eleparams.set("low-Mach-number solver",loma_);

  // parameters for stabilization
  eleparams.sublist("STABILIZATION") = params_.sublist("STABILIZATION");

  // parameters for stabilization
  eleparams.sublist("TURBULENCE MODEL") = params_.sublist("TURBULENCE MODEL");

  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("velnp",velnp_);
  discret_->SetState("vedenp",vedenp_);

  discret_->SetState("hist",hist_ );
  if (alefluid_)
  {
    discret_->SetState("dispnp", dispnp_);
    discret_->SetState("gridv", gridv_);
  }

  // call loop over elements
  discret_->Evaluate(eleparams,sysmat_,meshmovematrix_,residual_,Teuchos::null,Teuchos::null);
  discret_->ClearState();

  // finalize the system matrix
  sysmat_->Complete();

  if (meshmovematrix_ != Teuchos::null)
  {
    meshmovematrix_->Complete();
    // apply Dirichlet conditions to a non-diagonal matrix
    // (The Dirichlet rows will become all zero, no diagonal one.)
    meshmovematrix_->ApplyDirichlet(dirichtoggle_,false);
  }

  trueresidual_->Update(ResidualScaling(),*residual_,0.0);

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
 | current solution becomes most recent solution of next timestep       |
 |                                                           gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidImplicitTimeInt::TimeUpdate()
{
  // prev. acceleration becomes (n-1)-accel. of next time step
  accnm_->Update(1.0,*accn_,0.0);

  if (loma_ != "No")
  {
    // compute momentum/density time derivative
    TIMEINT_THETA_BDF2::CalculateAcceleration(
        vedenp_, veden_, vedenm_, accnm_,
        timealgo_, step_, theta_, dta_, dtp_,
        accn_);
  }
  else
  {
    // compute acceleration
    // note a(n+1) is directly stored in a(n),
    // hence we use a(n-1) as a(n) (see line above)
    TIMEINT_THETA_BDF2::CalculateAcceleration(
        velnp_, veln_, velnm_, accnm_,
        timealgo_, step_, theta_, dta_, dtp_,
        accn_);
  }

  // solution of this step becomes most recent solution of the last step
  velnm_->Update(1.0,*veln_ ,0.0);
  veln_ ->Update(1.0,*velnp_,0.0);

  if (alefluid_)
  {
    dispnm_->Update(1.0,*dispn_,0.0);
    dispn_ ->Update(1.0,*dispnp_,0.0);
  }

  // -------------------------------------------------------------------
  // treat impedance BC
  // note: these methods return without action, if the problem does not
  //       have impedance boundary conditions
  // -------------------------------------------------------------------
  discret_->ClearState();
  discret_->SetState("velnp",velnp_);
  discret_->SetState("hist",hist_);

  impedancebc_->FlowRateCalculation(time_,dta_);
  impedancebc_->OutflowBoundary(time_,dta_,theta_);

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
void FLD::FluidImplicitTimeInt::Output()
{

  //-------------------------------------------- output of solution

  if (step_%upres_ == 0)  //write solution
  {
    output_.NewStep    (step_,time_);
    output_.WriteVector("velnp", velnp_);

    // output (hydrodynamic) pressure
    Teuchos::RCP<Epetra_Vector> pressure = velpressplitter_.ExtractCondVector(velnp_);
    pressure->Scale(density_);
    output_.WriteVector("pressure", pressure);

    //output_.WriteVector("residual", trueresidual_);
    if (alefluid_)
      output_.WriteVector("dispnp", dispnp_);

    //only perform stress calculation when output is needed
    if (writestresses_)
    {
      RCP<Epetra_Vector> traction = CalcStresses();
      output_.WriteVector("traction",traction);
    }

    // write domain decomposition for visualization (only once!)
    if (step_==upres_)
     output_.WriteElementData();

    if (uprestart_ != 0 && step_%uprestart_ == 0) //add restart data
    {
      output_.WriteVector("accn", accn_);
      output_.WriteVector("veln", veln_);
      output_.WriteVector("velnm", velnm_);

      if (alefluid_)
      {
        output_.WriteVector("dispn", dispn_);
        output_.WriteVector("dispnm",dispnm_);
      }
      // also write impedance bc information if required
      // Note: this method acts only if there is an impedance BC
      impedancebc_->WriteRestart(output_);
    }
  }

  // write restart also when uprestart_ is not a integer multiple of upres_
  else if (uprestart_ != 0 && step_%uprestart_ == 0)
  {
    output_.NewStep    (step_,time_);
    output_.WriteVector("velnp", velnp_);
    //output_.WriteVector("residual", trueresidual_);
    if (alefluid_)
    {
      output_.WriteVector("dispnp", dispnp_);
      output_.WriteVector("dispn", dispn_);
      output_.WriteVector("dispnm",dispnm_);
    }

    //only perform stress calculation when output is needed
    if (writestresses_)
    {
      RCP<Epetra_Vector> traction = CalcStresses();
      output_.WriteVector("traction",traction);
    }

    output_.WriteVector("accn", accn_);
    output_.WriteVector("veln", veln_);
    output_.WriteVector("velnm", velnm_);

    // also write impedance bc information if required
    // Note: this method acts only if there is an impedance BC
    impedancebc_->WriteRestart(output_);
  }

  // dumping of turbulence statistics if required
  //if (special_flow_ == "lid_driven_cavity" && step_>=samstart_ && step_<=samstop_)
  if (special_flow_ != "no" && step_>=samstart_ && step_<=samstop_)
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
      else if (special_flow_ == "square_cylinder")
        turbulencestatistics_sqc_->DumpStatistics(step_);
    }
  }
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
void FLD::FluidImplicitTimeInt::ReadRestart(int step)
{
  IO::DiscretizationReader reader(discret_,step);
  time_ = reader.ReadDouble("time");
  step_ = reader.ReadInt("step");

  reader.ReadVector(velnp_,"velnp");
  reader.ReadVector(veln_, "veln");
  reader.ReadVector(velnm_,"velnm");
  reader.ReadVector(accn_ ,"accn");

  if (alefluid_)
  {
    reader.ReadVector(dispnp_,"dispnp");
    reader.ReadVector(dispn_ , "dispn");
    reader.ReadVector(dispnm_,"dispnm");
  }
  // also read impedance bc information if required
  // Note: this method acts only if there is an impedance BC
  impedancebc_->ReadRestart(reader);
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |                                                           chfoe 01/08|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidImplicitTimeInt::UpdateGridv()
{
  // get order of accuracy of grid velocity determination
  // from input file data
  const int order  = params_.get<int>("order gridvel");

  switch (order)
  {
    case 1:
      /* get gridvelocity from BE time discretisation of mesh motion:
           -> cheap
           -> easy
           -> limits FSI algorithm to first order accuracy in time

                  x^n+1 - x^n
             uG = -----------
                    Delta t                        */
      gridv_->Update(1/dta_, *dispnp_, -1/dta_, *dispn_, 0.0);
    break;
    case 2:
      /* get gridvelocity from BDF2 time discretisation of mesh motion:
           -> requires one more previous mesh position or displacemnt
           -> somewhat more complicated
           -> allows second order accuracy for the overall flow solution  */
      gridv_->Update(1.5/dta_, *dispnp_, -2.0, *dispn_, 0.0);
      gridv_->Update(0.5, *dispnm_, 1.0);
    break;
  }
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
void FLD::FluidImplicitTimeInt::SetInitialFlowField(
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
  else if(whichinitialfield==2 || whichinitialfield==3)
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

    // add random perturbation
    if(whichinitialfield==3)
    {
      const int numdim = params_.get<int>("number of velocity degrees of freedom");

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

      double bmvel=0;
      double mybmvel=0;
      double thisvel=0;
      // loop all nodes on the processor
      for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();++lnodeid)
      {
        // get the processor local node
        DRT::Node*  lnode      = discret_->lRowNode(lnodeid);
        // the set of degrees of freedom associated with the node
        vector<int> nodedofset = discret_->Dof(lnode);

        for(int index=0;index<numdim;++index)
        {
          int gid = nodedofset[index];
          int lid = dofrowmap->LID(gid);

          thisvel=(*velnp_)[lid];
          if (mybmvel*mybmvel < thisvel*thisvel) mybmvel=thisvel;
        }
      }

      // the noise is proportional to the bulk mean velocity of the
      // undisturbed initial field (=2/3*maximum velocity)
      mybmvel=2*mybmvel/3;
      discret_->Comm().MaxAll(&mybmvel,&bmvel,1);

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
          if(master == pbcmapmastertoslave_.end()) continue;
        }

        // add random noise on initial function field
        for(int index=0;index<numdim;++index)
        {
          int gid = nodedofset[index];

          double randomnumber = 2*((double)rand()-((double) RAND_MAX)/2.)/((double) RAND_MAX);

          double noise = perc * bmvel * randomnumber;

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
  else
  {
    dserror("no other initial fields than zero, function and beltrami are available up to now");
  }

  return;
} // end SetInitialFlowField


/*----------------------------------------------------------------------*
 | set time-step-related fields for low-Mach-number flow       vg 08/08 |
 *----------------------------------------------------------------------*/
void FLD::FluidImplicitTimeInt::SetTimeLomaFields(
   RCP<const Epetra_Vector> densnp,
   RCP<const Epetra_Vector> densn,
   RCP<const Epetra_Vector> densnm)
{
  // check vector compatibility and determine space dimension
  int numdim =-1;
  int numdof =-1;
  if (velnp_->MyLength()== (3* densnp->MyLength()))
  {
    numdim = 2;
    numdof = 3;
  }
  else if (velnp_->MyLength()== (4* densnp->MyLength()))
  {
    numdim = 3;
    numdof = 4;
  }
  else
    dserror("velocity/pressure and density vectors do not match in size");

  // get velocity dofs for vede-vectors as copies from vel-vectors
  vedenp_->Update(1.0,*velnp_,0.0);
  veden_->Update(1.0,*veln_,0.0);
  vedenm_->Update(1.0,*velnm_,0.0);

  vector<int>    Indices(numdof);
  vector<double> Values(numdof);

  // insert density values in vede-vectors
  // loop all nodes on the processor
  for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
  {
    Indices[numdim] = lnodeid*numdof + numdim;

    Values[numdim] = (*densnp)[lnodeid];
    vedenp_->ReplaceMyValues(1,&Values[numdim],&Indices[numdim]);

    Values[numdim] = (*densn)[lnodeid];
    veden_->ReplaceMyValues(1,&Values[numdim],&Indices[numdim]);

    Values[numdim] = (*densnm)[lnodeid];
    vedenm_->ReplaceMyValues(1,&Values[numdim],&Indices[numdim]);
  }

  return;

} // ScaTraTimIntImpl::SetTimeLomaFields


/*----------------------------------------------------------------------*
 | set outer-iteration-related fields for low-Mach-number flow vg 08/08 |
 *----------------------------------------------------------------------*/
void FLD::FluidImplicitTimeInt::SetIterLomaFields(
   RCP<const Epetra_Vector> densnp)
{
  // check vector compatibility and determine space dimension
  int numdim =-1;
  int numdof =-1;
  if (velnp_->MyLength()== (3* densnp->MyLength()))
  {
    numdim = 2;
    numdof = 3;
  }
  else if (velnp_->MyLength()== (4* densnp->MyLength()))
  {
    numdim = 3;
    numdof = 4;
  }
  else
    dserror("velocity/pressure and density vectors do not match in size");

  // get velocity dofs for vedenp-vector as copies from velnp-vector
  vedenp_->Update(1.0,*velnp_,0.0);

  vector<int>    Indices(numdof);
  vector<double> Values(numdof);

  // insert density values in vedenp-vector
  // loop all nodes on the processor
  for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
  {
    Indices[numdim] = lnodeid*numdof + numdim;

    Values[numdim] = (*densnp)[lnodeid];
    vedenp_->ReplaceMyValues(1,&Values[numdim],&Indices[numdim]);
  }

  return;

} // ScaTraTimIntImpl::SetIterLomaFields


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | evaluate error for test cases with analytical solutions   gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidImplicitTimeInt::EvaluateErrorComparedToAnalyticalSol()
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

    // call loop over elements (assemble nothing)
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
 | solve stationary fluid problem                              gjb 10/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidImplicitTimeInt::SolveStationaryProblem()
{
  // time measurement: time loop (stationary) --- start TimeMonitor tm2
  TEUCHOS_FUNC_TIME_MONITOR(" + time loop");

  // override time integration parameters in order to avoid misuse in
  // NonLinearSolve method below
  const double origdta = dta_;
  dta_= 1.0;
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
   time_ += origdta;
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

     // other parameters needed by the elements
     eleparams.set("total time",time_);
     eleparams.set("delta time",origdta);
     eleparams.set("thsl",1.0); // no timefac in stationary case
     eleparams.set("fs subgrid viscosity",fssgv_);

     // set vector values needed by elements
     discret_->ClearState();
     discret_->SetState("velnp",velnp_);
     // predicted dirichlet values
     // velnp then also holds prescribed new dirichlet values
     // dirichtoggle is 1 for dirichlet dofs, 0 elsewhere
     discret_->EvaluateDirichlet(eleparams,velnp_,null,null,dirichtoggle_);
     discret_->ClearState();

     // evaluate Neumann b.c.
     //eleparams.set("inc_density",density_);

     discret_->SetState("vedenp",vedenp_);
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
    NonlinearSolve();

    // -------------------------------------------------------------------
    //                    calculate lift'n'drag forces
    // -------------------------------------------------------------------
    int liftdrag = params_.get<int>("liftdrag");

    if(liftdrag == 0); // do nothing, we don't want lift & drag
    if(liftdrag == 1)
      dserror("how did you manage to get here???");
    if(liftdrag == 2)
      LiftDrag();

    // -------------------------------------------------------------------
    //                         output of solution
    // -------------------------------------------------------------------
    Output();

  } // end of time loop

} // FluidImplicitTimeInt::SolveStationaryProblem


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  calculate traction vector at (Dirichlet) boundary (public) gjb 07/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
Teuchos::RCP<Epetra_Vector> FLD::FluidImplicitTimeInt::CalcStresses()
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
      (*integratedshapefunc)[i] = (*trueresidual_)[i]/(*integratedshapefunc)[i];
    }
  }

  return integratedshapefunc;

} // FluidImplicitTimeInt::CalcStresses()


/*----------------------------------------------------------------------*
 | Destructor dtor (public)                                  gammi 04/07|
 *----------------------------------------------------------------------*/
FLD::FluidImplicitTimeInt::~FluidImplicitTimeInt()
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
void FLD::FluidImplicitTimeInt::LiftDrag() const
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

      // -------------------------------------------------------------------
      // add calculated lift and drag to mean value calculation (statistics)
      // -------------------------------------------------------------------
      if(special_flow_ == "square_cylinder" && step_>=samstart_ && step_<=samstop_)
      {
        double dragforce = resultvec[0];
        double liftforce = resultvec[1];
        turbulencestatistics_sqc_->DoLiftDragTimeSample(dragforce,liftforce);
      }

    } // end: loop over L&D labels
    if (myrank_== 0)
    {
      cout << "\n";
    }
  }
}//FluidImplicitTimeInt::LiftDrag

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | filter quantities for dynamic Smagorinsky model. Compute averaged    |
 | values for LijMij and MijMij.                             gammi 02/08|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidImplicitTimeInt::ApplyFilterForDynamicComputationOfCs()
{
  // check plausiblity of dimension of problem --- only 3d is
  // valid for turbulence calculations
  const int numdim  = params_.get<int>("number of velocity degrees of freedom");
  if(numdim!=3)
  {
    dserror("Only 3d problems are allowed for dynamic Smagorinsky");
  }

  // get the dofrowmap for access of velocity dofs
  const Epetra_Map* dofrowmap       = discret_->DofRowMap();

  // generate a parameterlist for communication and control
  ParameterList filterparams;

  // set filter action for elements
  filterparams.set("action","calc_fluid_box_filter");

  // set state vector to pass distributed vector to the element
  discret_->ClearState();
  discret_->SetState("u and p (trial)",velnp_);

  // ----------------------------------------------------------
  // ----------------------------------------------------------
  // compute smoothed (averaged) velocities and stresses for
  // all nodes
  // ----------------------------------------------------------
  // ----------------------------------------------------------

  // get additional connectivity for periodic boundary conditions
  map<int,vector<int> > mapmastertoslave
    =
    pbc_->ReturnAllCoupledNodesOnThisProc();


  // loop all nodes on the processor
  for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
  {
    // ----------------------------------------
    // get the processor local node
    DRT::Node*  lnode       = discret_->lRowNode(lnodeid);

    // the set of degrees of freedom associated with the node
    vector<int> nodedofset = discret_->Dof(lnode);

    // ----------------------------------------
    // determine a patch of all elements adjacent to this nodex

    // check whether the node is on a wall, i.e. all velocity
    // dofs are Dirichlet constrained
    int is_no_slip_node =0;
    for(int index=0;index<numdim;++index)
    {
      int gid = nodedofset[index];
      int lid = dofrowmap->LID(gid);

      if ((*dirichtoggle_)[lid]==1)
      {
        is_no_slip_node++;
      }
    }

    // skip this node if it is on a wall
    if (is_no_slip_node == numdim)
    {
      continue;
    }

    // now check whether we have a pbc condition on this node
    vector<DRT::Condition*> mypbc;

    lnode->GetCondition("SurfacePeriodic",mypbc);

    // check whether a periodic boundary condition is active
    // on this node
    bool ispbcmaster = false;

    if (mypbc.size()>0)
    {
      // get the list of all his slavenodes
      map<int, vector<int> >::iterator master
        =
        mapmastertoslave.find(lnode->Id());

      // slavenodes are ignored
      if(master == mapmastertoslave.end())
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
      for (unsigned slavecount = 0;
           slavecount<mapmastertoslave[lnode->Id()].size();
           ++slavecount)
      {
        // get the corresponding slavenodes
        DRT::Node*  slavenode = discret_->gNode(mapmastertoslave[lnode->Id()][slavecount]);

        // add the elements
        for(int rr=0;rr<slavenode->NumElement();++rr)
        {
          patcheles.push_back(slavenode->Elements()[rr]);
        }
      }
    }

    // ----------------------------------------
    // the patch is determined right now --- what follows
    // now is the averaging over this patch


    // define element matrices and vectors --- they are used to
    // transfer information into the element routine and back
    Epetra_SerialDenseMatrix ep_reystress_hat(3,3);
    Epetra_SerialDenseMatrix ep_modeled_stress_grid_scale_hat(3,3);
    Epetra_SerialDenseVector ep_velnp_hat (3);
    Epetra_SerialDenseVector dummy1;
    Epetra_SerialDenseVector dummy2;

    // the patch volume has to be initialised to zero
    double patchvolume = 0;

    // loop all adjacent elements to this node
    for (unsigned nele=0;nele<patcheles.size();++nele)
    {
      // get the adjacent element
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
                                ep_velnp_hat,dummy1,dummy2);
      if (err) dserror("Proc %d: Element %d returned err=%d",
                       discret_->Comm().MyPID(),nbele->Id(),err);

      // get contribution to patch volume of this element. Add it up.
      double volume_contribution =filterparams.get<double>("volume_contribution");

      patchvolume+=volume_contribution;
    }

    // wrap Epetra Object in Blitz array
    blitz::Array<double, 1> velnp_hat(
      ep_velnp_hat.Values(),
      blitz::shape(ep_velnp_hat.Length()),
      blitz::neverDeleteData);

    blitz::Array<double, 2> reystress_hat(
      ep_reystress_hat.A(),
      blitz::shape(ep_reystress_hat.M(),
                   ep_reystress_hat.N()),
      blitz::neverDeleteData,
      blitz::ColumnMajorArray<2>());

    blitz::Array<double, 2> modeled_stress_grid_scale_hat(
      ep_modeled_stress_grid_scale_hat.A(),
      blitz::shape(ep_modeled_stress_grid_scale_hat.M(),
                   ep_modeled_stress_grid_scale_hat.N()),
      blitz::neverDeleteData,
      blitz::ColumnMajorArray<2>());


    // normalize the computed convolution products by the complete patchvolume
    reystress_hat                /=patchvolume;
    modeled_stress_grid_scale_hat/=patchvolume;
    velnp_hat                    /=patchvolume;

    // now assemble the computed values into the global vector
    double val = 0;
    int    id  = (lnode->Id());

    for (int idim =0;idim<3;++idim)
    {
      val = velnp_hat(idim);
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
      for (unsigned slavecount = 0;slavecount<mapmastertoslave[lnode->Id()].size();++slavecount)
      {
        // get the corresponding slavenodes
        DRT::Node*  slavenode = discret_->gNode(mapmastertoslave[lnode->Id()][slavecount]);

        int    slaveid  = (slavenode->Id());

        for (int idim =0;idim<3;++idim)
        {
          val = (velnp_hat(idim));
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


  // -----------------------------------------------------------
  // -----------------------------------------------------------
  // computation of LijMij, MijMij and Cs_delta_sq

  // up to now, only averaging in homogeneous planes is
  // implemented. Otherwise, no aberaging will be applied.
  //
  //
  // -----------------------------------------------------------
  // -----------------------------------------------------------

  // -----------------------------------------------------------
  // initialise plane averaging of Cs for turbulent channel flow
  vector<int>  count_for_average      ;
  vector<int>  local_count_for_average;

  vector <double> local_ele_sum_LijMij;
  vector <double> local_ele_sum_MijMij;

  if (params_.sublist("TURBULENCE MODEL").get<string>("CANONICAL_FLOW","no")
      ==
      "channel_flow_of_height_2")
  {
    // get ordered layers of elements in which LijMij and MijMij are averaged
    if (planecoords_ == null)
    {
      planecoords_ = rcp( new vector<double>((turbulencestatistics_->ReturnNodePlaneCoords()).size()));
    }

    (*planecoords_) = turbulencestatistics_->ReturnNodePlaneCoords();

    averaged_LijMij_->resize((*planecoords_).size()-1);
    averaged_MijMij_->resize((*planecoords_).size()-1);

    count_for_average      .resize((*planecoords_).size()-1);
    local_count_for_average.resize((*planecoords_).size()-1);

    local_ele_sum_LijMij.resize((*planecoords_).size()-1);
    local_ele_sum_MijMij.resize((*planecoords_).size()-1);

    for (unsigned rr=0;rr<(*planecoords_).size()-1;++rr)
    {
      (*averaged_LijMij_)    [rr]=0;
      (*averaged_MijMij_)    [rr]=0;
      local_ele_sum_LijMij   [rr]=0;
      local_ele_sum_MijMij   [rr]=0;
      local_count_for_average[rr]=0;
    }
  }

  // -----------------------------------------------------------
  // Compute Cs_delta_sq on elements, return LijMij and MijMij
  // for plane averaging in case of turbulent channel flows

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

  // loop all elements on this proc (including ghosted ones)
  for (int nele=0;nele<discret_->NumMyColElements();++nele)
  {
    // -----------------------------------------------------------
    // Compute LijMij, MijMij and Cs_delta_sq_

    // get the element
    //DRT::Element* ele = discret_->lRowElement(nele);
    DRT::Element* ele = discret_->lColElement(nele);

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



    // -----------------------------------------------------------
    // do averaging of Cs for turbulent channel flow

    // initialise plane averaging of Cs for turbulent channel flow
    if (params_.sublist("TURBULENCE MODEL").get<string>("CANONICAL_FLOW","no")
        ==
        "channel_flow_of_height_2")
    {
      // only row elements are included in the global averaging
      if(ele->Owner() == myrank_)
      {
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
    }
  }

  if (params_.sublist("TURBULENCE MODEL").get<string>("CANONICAL_FLOW","no")
      ==
      "channel_flow_of_height_2")
  {
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
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FLD::FluidImplicitTimeInt::IntegrateInterfaceShape(std::string condname)
{
  ParameterList eleparams;
  // set action for elements
  eleparams.set("action","integrate_Shapefunction");

  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // create vector (+ initialization with zeros)
  Teuchos::RCP<Epetra_Vector> integratedshapefunc = LINALG::CreateVector(*dofrowmap,true);

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
void FLD::FluidImplicitTimeInt::UseBlockMatrix(Teuchos::RCP<std::set<int> > condelements,
                                               const LINALG::MultiMapExtractor& domainmaps,
                                               const LINALG::MultiMapExtractor& rangemaps,
                                               bool splitmatrix)
{
  Teuchos::RCP<LINALG::BlockSparseMatrix<FLD::UTILS::InterfaceSplitStrategy> > mat;

  if (splitmatrix)
  {
    // (re)allocate system matrix
    mat = Teuchos::rcp(new LINALG::BlockSparseMatrix<FLD::UTILS::InterfaceSplitStrategy>(domainmaps,rangemaps,108,false,true));
    mat->SetCondElements(condelements);
    sysmat_ = mat;
  }

  // if we never build the matrix nothing will be done
  if (params_.get<bool>("shape derivatives"))
  {
    // allocate special mesh moving matrix
    mat = Teuchos::rcp(new LINALG::BlockSparseMatrix<FLD::UTILS::InterfaceSplitStrategy>(domainmaps,rangemaps,108,false,true));
    mat->SetCondElements(condelements);
    meshmovematrix_ = mat;
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::FluidImplicitTimeInt::LinearRelaxationSolve(Teuchos::RCP<Epetra_Vector> relax)
{
  TEUCHOS_FUNC_TIME_MONITOR("FluidImplicitTimeInt::LinearRelaxationSolve");

  //
  // Special linear solve used for steepest descent relaxation as well as
  // Jacobian-free Newton-Krylov on the FSI interface equations. The later one
  // presents a special challenge, as we have to solve the same linear system
  // repeatedly for different rhs. That is why we need the inrelaxation_ flag.
  //
  // Additionally we might want to include the mesh derivatives to get optimal
  // convergance in the Newton loop.
  //
  // This adds even more state to the fluid algorithm class, which is a bad
  // thing. And the explicit storage of the Dirichlet lines is
  // required. However, we do not need any special element code to perform the
  // steepest descent calculation. This is quite a benefit as the special code
  // in the old discretization was a real nightmare.
  //

  if (not inrelaxation_)
  {
    // setup relaxation matrices just once
    //
    // We use these matrices for several solves in Jacobian-free Newton-Krylov
    // solves of the FSI interface equations.

    const Epetra_Map* dofrowmap = discret_->DofRowMap();
    Teuchos::RCP<Epetra_Vector> griddisp = LINALG::CreateVector(*dofrowmap,false);

    // set the grid displacement independent of the trial value at the
    // interface
    griddisp->Update(1., *dispnp_, -1., *dispn_, 0.);

    // dirichtoggle_ has already been set up

    // zero out the stiffness matrix
    sysmat_->Zero();

    // zero out residual, no neumann bc
    residual_->PutScalar(0.0);

    // Get matrix for mesh derivatives. This is not meant to be efficient.
    if (params_.get<bool>("shape derivatives"))
    {
      if (meshmatrix_==Teuchos::null)
      {
        meshmatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*SystemMatrix()));
      }
      else
      {
        meshmatrix_->Zero();
      }
    }

    ParameterList eleparams;
    eleparams.set("action","calc_fluid_systemmat_and_residual");
    eleparams.set("total time",time_);
    eleparams.set("thsl",theta_*dta_);
    eleparams.set("dt",dta_);
    eleparams.set("Linearisation",newton_);
    eleparams.set("low-Mach-number solver",loma_);

    // parameters for stabilization
    eleparams.sublist("STABILIZATION") = params_.sublist("STABILIZATION");

    // parameters for stabilization
    eleparams.sublist("TURBULENCE MODEL") = params_.sublist("TURBULENCE MODEL");

    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("velnp",velnp_);
    discret_->SetState("vedenp",vedenp_);
    discret_->SetState("hist",zeros_ );
    discret_->SetState("dispnp", griddisp);
    discret_->SetState("gridv", zeros_);

    // call loop over elements
    discret_->Evaluate(eleparams,sysmat_,meshmovematrix_,residual_);
    discret_->ClearState();

    // finalize the system matrix
    sysmat_->Complete();

    if (meshmatrix_!=Teuchos::null)
    {
      meshmatrix_->Complete();
    }

    //--------- Apply dirichlet boundary conditions to system of equations
    //          residual displacements are supposed to be zero at
    //          boundary conditions
    dirichletlines_ = SystemMatrix()->ExtractDirichletLines(dirichtoggle_);
    sysmat_->ApplyDirichlet(dirichtoggle_);
  }

  // No, we do not want to have any rhs. There cannot be any.
  residual_->PutScalar(0.0);

  if (meshmatrix_!=Teuchos::null)
  {
    // Calculate rhs due to mesh movement induced by the interface
    // displacement.
    meshmatrix_->Apply(*gridv_,*residual_);
    residual_->Scale(-dta_);
  }

  //--------- Apply dirichlet boundary conditions to system of equations
  //          residual displacements are supposed to be zero at
  //          boundary conditions
  incvel_->PutScalar(0.0);
  LINALG::ApplyDirichlettoSystem(incvel_,residual_,relax,dirichtoggle_);

  //-------solve for residual displacements to correct incremental displacements
  solver_.Solve(sysmat_->EpetraOperator(),incvel_,residual_,not inrelaxation_,not inrelaxation_);

  // and now we need the reaction forces

  if (dirichletlines_->Apply(*incvel_, *trueresidual_)!=0)
    dserror("dirichletlines_->Apply() failed");

  if (meshmatrix_!=Teuchos::null)
  {
    // Calculate rhs due to mesh movement induced by the interface
    // displacement.
    meshmatrix_->Apply(*gridv_,*residual_);
    trueresidual_->Update(dta_,*residual_,1.0);
  }

  trueresidual_->Scale(-ResidualScaling());

  if (not inrelaxation_)
    inrelaxation_ = true;
}


#endif /* CCADISCRET       */
