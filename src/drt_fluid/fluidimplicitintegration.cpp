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
  // af-generalized-alpha parameters: gamma_ = 0.5 + alphaM_ - alphaF_
  // (may be reset below when starting algorithm is used)
  alphaM_   = params_.get<double>("alpha_M");
  alphaF_   = params_.get<double>("alpha_F");
  gamma_    = params_.get<double>("gamma");

  // number of steps for starting algorithm
  numstasteps_ = params_.get<int> ("number of start steps");
  // starting algorithm only for af-generalized-alpha so far
  // -> check for time-integration scheme and reasonability of number of steps
  startalgo_ = false;
  if (numstasteps_ > 0)
  {
    if (timealgo_ != timeint_afgenalpha)
      dserror("no starting algorithm supported for schemes other than af-gen-alpha");
    else startalgo_= true;
    if (numstasteps_>stepmax_)
      dserror("more steps for starting algorithm than steps overall");
  }

  // parameter for linearization scheme (fixed-point-like or Newton)
  newton_ = params_.get<string>("Linearisation");

  // use of predictor
  // (might be used for af-generalized-alpha, but not yet activated)
  //predictor_ = params_.get<string>("predictor","steady_state_predictor");

  // form of convective term
  convform_ = params_.get<string>("form of convective term","convective");

  // fine-scale subgrid viscosity?
  fssgv_ = params_.get<string>("fs subgrid viscosity","No");

  // -------------------------------------------------------------------
  // connect degrees of freedom for periodic boundary conditions
  // -------------------------------------------------------------------
  pbc_ = rcp(new PeriodicBoundaryConditions (discret_));
  pbc_->UpdateDofsForPeriodicBoundaryConditions();

  pbcmapmastertoslave_ = pbc_->ReturnAllCoupledNodesOnThisProc();

  discret_->ComputeNullSpaceIfNecessary(solver_.Params(),true);

  // ensure that degrees of freedom in the discretization have been set
  if (!discret_->Filled() || !actdis->HaveDofs()) discret_->FillComplete();

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
  // velocity/pressure at time n+1, n and n-1
  velnp_        = LINALG::CreateVector(*dofrowmap,true);
  veln_         = LINALG::CreateVector(*dofrowmap,true);
  velnm_        = LINALG::CreateVector(*dofrowmap,true);

  // acceleration at time n+1 and n
  accnp_        = LINALG::CreateVector(*dofrowmap,true);
  accn_         = LINALG::CreateVector(*dofrowmap,true);

  // velocity/density at time n+1
  vedenp_       = LINALG::CreateVector(*dofrowmap,true);

  // vectors only required for af-generalized-alpha scheme
  if (timealgo_==timeint_afgenalpha)
  {
    // velocity/pressure at time n+alpha_F
    velaf_        = LINALG::CreateVector(*dofrowmap,true);

    // acceleration/(density time derivative) at time n+alpha_M
    accam_        = LINALG::CreateVector(*dofrowmap,true);

    // velocity/density at time n+alpha_M and n+alpha_F
    vedeam_       = LINALG::CreateVector(*dofrowmap,true);
    vedeaf_       = LINALG::CreateVector(*dofrowmap,true);
  }

  // vectors only required for low-Mach-number flow
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

  // a vector of zeros to be used to enforce zero dirichlet boundary conditions
  zeros_   = LINALG::CreateVector(*dofrowmap,true);

  // object holds maps/subsets for DOFs subjected to Dirichlet BCs and otherwise
  dbcmaps_ = Teuchos::rcp(new LINALG::MapExtractor());
  {
    ParameterList eleparams;
    // other parameters needed by the elements
    eleparams.set("total time",time_);
    discret_->EvaluateDirichlet(eleparams, zeros_, Teuchos::null, Teuchos::null,
                                Teuchos::null, dbcmaps_);
    zeros_->PutScalar(0.0); // just in case of change
  }

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
  incvel_ = LINALG::CreateVector(*dofrowmap,true);

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
    string homdir = modelparams->get<string>("HOMDIR","not_specified");

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
              homdir != "xz")
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
              homdir != "xz")
          {
            cout << "      no homogeneous directions specified --- so we just use pointwise clipping for Cs\n";
          }
        }
        cout << &endl;
      }

      if (special_flow_ == "channel_flow_of_height_2" or
          special_flow_ == "loma_channel_flow_of_height_2")
      {
        cout << "                             " ;
        cout << "Turbulence statistics are evaluated ";
        cout << "for a turbulent channel flow.\n";
        cout << "                             " ;
        cout << "The solution is averaged over the homogeneous ";
        cout << homdir;
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

      // get one instance of the dynamic Smagorinsky class
      DynSmag_=rcp(new DynSmagFilter(discret_            ,
                                     pbcmapmastertoslave_,
                                     params_             ));

    }
  }

  // -------------------------------------------------------------------
  // initialize turbulence-statistics evaluation
  // -------------------------------------------------------------------
  //
  statisticsmanager_=rcp(new TurbulenceStatisticManager(*this));

  if (special_flow_ != "no")
  {
    // parameters for sampling/dumping period
    samstart_  = modelparams->get<int>("SAMPLING_START",1);
  }

  // -------------------------------------------------------------------
  // initialize outflow boundary stabilization if required
  // -------------------------------------------------------------------
  ParameterList *  stabparams=&(params_.sublist("STABILIZATION"));

  // flag for potential Neumann-type outflow stabilization
  outflow_stab_ = stabparams->get<string>("OUTFLOW_STAB","no_outstab");

  // vector containing potential Neumann-type outflow stabilization term
  if (outflow_stab_ == "yes_outstab")
    outflow_= LINALG::CreateVector(*dofrowmap,true);

  // -------------------------------------------------------------------
  // necessary only for the VM3 approach:
  // fine-scale solution vector + respective ouptput
  // -------------------------------------------------------------------
  if (fssgv_ != "No")
  {
    fsvelnp_  = LINALG::CreateVector(*dofrowmap,true);

    if (myrank_ == 0)
    {
      // Output
      cout << "Fine-scale subgrid-viscosity approach based on AVM3: ";
      cout << &endl << &endl;
      cout << params_.get<string>("fs subgrid viscosity");
      cout << " with Smagorinsky constant Cs= ";
      cout << modelparams->get<double>("C_SMAGORINSKY") ;
      cout << &endl << &endl << &endl;
    }
  }

  // construct impedance bc wrapper
  impedancebc_ = rcp(new UTILS::FluidImpedanceWrapper(discret_, output_, dta_) );

  // get constant density variable for incompressible flow
  // set vedenp/af/am-vector values to 1.0 for incompressible flow
  // set density variable to 1.0 for low-Mach-number flow
  if (loma_ == "No")
  {
    ParameterList eleparams;
    eleparams.set("action","get_density");
    discret_->Evaluate(eleparams,null,null,null,null,null);
    density_ = eleparams.get("density", 1.0);
    if (density_ <= 0.0) dserror("received illegal density value");
    vedenp_->PutScalar(1.0);
    if (timealgo_==timeint_afgenalpha)
    {
      vedeam_->PutScalar(1.0);
      vedeaf_->PutScalar(1.0);
    }
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

  // distinguish stationary and instationary case
  if (timealgo_==timeint_stationary) SolveStationaryProblem();
  else TimeLoop();

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
  const int dyntype = params_.get<int>("type of nonlinear solve");

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
    //                       output to screen
    // -------------------------------------------------------------------
    if (myrank_==0)
    {
      switch (timealgo_)
      {
      case timeint_one_step_theta:
        printf("TIME: %11.4E/%11.4E  DT = %11.4E   One-Step-Theta    STEP = %4d/%4d \n",
              time_,maxtime_,dta_,step_,stepmax_);
        break;
      case timeint_afgenalpha:
        printf("TIME: %11.4E/%11.4E  DT = %11.4E  Generalized-Alpha  STEP = %4d/%4d \n",
               time_,maxtime_,dta_,step_,stepmax_);
        break;
      case timeint_bdf2:
        printf("TIME: %11.4E/%11.4E  DT = %11.4E       BDF2          STEP = %4d/%4d \n",
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
    // -------------------------------------------------------------------
    TimeUpdate();

    // -------------------------------------------------------------------
    //  lift'n'drag forces, statistics time sample and output of solution
    //  and statistics
    // -------------------------------------------------------------------
    StatisticsAndOutput();

    // -------------------------------------------------------------------
    // evaluate error for test flows with analytical solutions
    // -------------------------------------------------------------------
    EvaluateErrorComparedToAnalyticalSol();

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

  // for BDF2, theta is set by the time-step sizes, 2/3 for const. dt
  if (timealgo_==timeint_bdf2) theta_ = (dta_+dtp_)/(2.0*dta_ + dtp_);

  // -------------------------------------------------------------------
  // set part(s) of the rhs vector(s) belonging to the old timestep
  // for low-Mach-number flow: distinguish momentum and continuity part
  // (continuity part only meaningful for low-Mach-number flow)
  // For low-Mach-number flow, all velocity values are multiplied by the
  // respective density values.
  //
  // Stationary/af-generalized-alpha:
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
          timealgo_, dta_, dtp_,
          velnp_);
    }
  }

  // -------------------------------------------------------------------
  //         evaluate Dirichlet and Neumann boundary conditions
  // -------------------------------------------------------------------
  {
    ParameterList eleparams;

    // total time required for Dirichlet conditions
    eleparams.set("total time",time_);

    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("velnp",velnp_);

    // predicted dirichlet values
    // velnp then also holds prescribed new dirichlet values
    discret_->EvaluateDirichlet(eleparams,velnp_,null,null,null);
    discret_->ClearState();

    // set all parameters and states required for Neumann conditions
    if (timealgo_==timeint_afgenalpha)
    {
      eleparams.set("total time",time_-(1-alphaF_)*dta_);
      eleparams.set("thsl",1.0);
    }
    else
    {
      eleparams.set("total time",time_);
      eleparams.set("thsl",theta_*dta_);
    }
    // For generalized-alpha, the following is an approximation,
    // since actually vedeaf would be required, which is not yet
    // available, though.
    discret_->SetState("vedenp",vedenp_);
    neumann_loads_->PutScalar(0.0);

    // evaluate Neumann conditions
    discret_->EvaluateNeumann(eleparams,*neumann_loads_);
    discret_->ClearState();
  }

  // -------------------------------------------------------------------
  //  For af-generalized-alpha time-integration scheme:
  //  set "pseudo-theta", calculate initial accelerations according to
  //  prescribed Dirichlet values for generalized-alpha time
  //  integration as well as velocities, pressures, densities,
  //  accelerations and density time derivatives at intermediate time
  //  steps n+alpha_F and n+alpha_M, respectively, for first iteration.
  // -------------------------------------------------------------------
  if (timealgo_==timeint_afgenalpha)
  {
    // starting algorithm
    if (startalgo_)
    {
      // use backward-Euler-type parameter combination
      if (step_<=numstasteps_)
      {
        alphaM_ = 1.0;
        alphaF_ = 1.0;
        gamma_  = 1.0;
      }
      else
      {
        alphaM_ = params_.get<double>("alpha_M");
        alphaF_ = params_.get<double>("alpha_F");
        gamma_  = params_.get<double>("gamma");
        startalgo_ = false;
      }
    }

    // set "pseudo-theta" for af-generalized-alpha scheme
    theta_ = alphaF_*gamma_/alphaM_;

    // --------------------------------------------------
    // adjust accnp according to Dirichlet values of velnp
    //
    //                                  n+1     n
    //                               vel   - vel
    //       n+1      n  gamma-1.0      (0)
    //    acc    = acc * --------- + ------------
    //       (0)           gamma      gamma * dt
    //
    // in case of conservative form: velocity*density
    if (convform_ == "conservative" and loma_ != "No")
    {
      accnp_->Multiply(1.0,*velnp_,*vedenp_,0.0);
      accnp_->Multiply(-1.0,*veln_,*veden_,1.0);
    }
    else accnp_->Update(1.0,*velnp_,-1.0,*veln_,0.0);
    accnp_->Update((gamma_-1.0)/gamma_,*accn_,1.0/(gamma_*dta_));

    //       n+alphaM                n+1                      n
    //    acc         = alpha_M * acc     + (1-alpha_M) *  acc
    //       (i)                     (i)
    accam_->Update((alphaM_),*accnp_,(1.0-alphaM_),*accn_,0.0);

    //       n+alphaF              n+1                   n
    //      u         = alpha_F * u     + (1-alpha_F) * u
    //       (i)                   (i)
    velaf_->Update((alphaF_),*velnp_,(1.0-alphaF_),*veln_,0.0);

    // values only required for low-Mach-number flow
    if (loma_ != "No")
    {
      vedeaf_->Update((alphaF_),*vedenp_,(1.0-alphaF_),*veden_,0.0);
      vedeam_->Update((alphaM_),*vedenp_,(1.0-alphaM_),*veden_,0.0);
    }
  }

  // -------------------------------------------------------------------
  //           preparation of AVM3-based scale separation
  // -------------------------------------------------------------------
  if (step_==1 and fssgv_ != "No") AVM3Preparation();

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

  int  itnum = 0;
  int  itemax = 0;
  bool stopnonliniter = false;

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

      // add potential Neumann-type outflow stabilization term
      if (outflow_stab_ == "yes_outstab")
      {
        discret_->ClearState();
        eleparams.set("outflow stabilization",outflow_stab_);
        if (timealgo_==timeint_afgenalpha)
        {
          eleparams.set("thsl",1.0);
          discret_->SetState("velnp",velaf_);
          discret_->SetState("vedenp",vedeaf_);
        }
        else
        {
          eleparams.set("thsl",theta_*dta_);
          discret_->SetState("velnp",velnp_);
          discret_->SetState("vedenp",vedenp_);
        }
        outflow_->PutScalar(0.0);
        discret_->EvaluateNeumann(eleparams,*outflow_);
        discret_->ClearState();

        // add Neumann-type outflow term to residual vector
        residual_->Update(1.0,*outflow_,1.0);
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

      // set general element parameters
      eleparams.set("thsl",theta_*dta_);
      eleparams.set("dt",dta_);
      eleparams.set("form of convective term",convform_);
      eleparams.set("fs subgrid viscosity",fssgv_);
      eleparams.set("Linearisation",newton_);
      eleparams.set("low-Mach-number solver",loma_);
      eleparams.set("eos factor",eosfac_);

      // parameters for stabilization
      eleparams.sublist("STABILIZATION") = params_.sublist("STABILIZATION");

      // parameters for stabilization
      eleparams.sublist("TURBULENCE MODEL") = params_.sublist("TURBULENCE MODEL");

      // set general vector values needed by elements
      discret_->ClearState();
      discret_->SetState("hist"  ,hist_ );
      if (alefluid_)
      {
        discret_->SetState("dispnp", dispnp_);
        discret_->SetState("gridv", gridv_);
      }

      // set scheme-specific element parameters and vector values
      if (timealgo_==timeint_stationary)
      {
        eleparams.set("action","calc_fluid_stationary_systemmat_and_residual");
        eleparams.set("using generalized-alpha time integration",false);
        eleparams.set("total time",time_);

        discret_->SetState("velnp",velnp_);
        discret_->SetState("vedenp",vedenp_);
      }
      else if (timealgo_==timeint_afgenalpha)
      {
        eleparams.set("action","calc_fluid_afgenalpha_systemmat_and_residual");
        eleparams.set("using generalized-alpha time integration",true);
        eleparams.set("total time",time_-(1-alphaF_)*dta_);
        eleparams.set("timefacrhs",alphaM_/(gamma_*dta_));

        discret_->SetState("velnp", velaf_ );
        discret_->SetState("vedenp",vedeaf_);
        discret_->SetState("accam", accam_ );
        if (convform_ == "conservative") discret_->SetState("vedeam",vedenp_);
        else                             discret_->SetState("vedeam",vedeam_);
      }
      else
      {
        eleparams.set("action","calc_fluid_systemmat_and_residual");
        eleparams.set("using generalized-alpha time integration",false);
        eleparams.set("total time",time_);

        discret_->SetState("velnp",velnp_);
        discret_->SetState("vedenp",vedenp_);
      }

      //----------------------------------------------------------------------
      // decide whether AVM3-based solution approach or standard approach
      //----------------------------------------------------------------------
      if (fssgv_ != "No") AVM3Separation();

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


        //---------------------------surface tension update
        if (alefluid_ and freesurface_->Relevant())
        {

// 2D or TET: To possibilities work. FSTENS1 is a direct implementation of Wall et
// al. eq. (25) with node normals obtained by weighted assembly of element
// normals. Because geometric considerations are used to find the normals of
// our flat (!) surface elements no second derivatives appear. FSTENS2 employs the
// divergence theorem acc. to Saksono eq. (24).

#define FSTENS2
#undef FSTENS1


          // select free surface elements
          std::string condname = "FREESURFCoupling";

          // get a vector layout from the discretization to construct matching
          // vectors and matrices
          //                 local <-> global dof numbering
          const Epetra_Map* dofrowmap = discret_->DofRowMap();

          //vector ndnorm0 with pressure-entries is needed for EvaluateCondition
          Teuchos::RCP<Epetra_Vector> ndnorm0 = LINALG::CreateVector(*dofrowmap,true);

          ParameterList eleparams;

#ifdef FSTENS1
          // set action for elements, calc_surface_tension uses node normals
          eleparams.set("action","calc_node_normal");

          //call loop over elements, note: normal vectors do not yet have length = 1.0
          discret_->ClearState();
          discret_->SetState("dispnp", dispnp_);
          discret_->EvaluateCondition(eleparams,ndnorm0,condname);
#endif

          // set action for elements
          eleparams.set("action","calc_surface_tension");
          discret_->ClearState();
          discret_->SetState("dispnp", dispnp_);
#ifdef FSTENS1
          discret_->SetState("normals", ndnorm0);
#endif
          discret_->EvaluateCondition(eleparams,sysmat_,Teuchos::null,residual_,Teuchos::null,Teuchos::null,condname);
          discret_->ClearState();
        }
        //---------------------------end of surface tension update

        if (timealgo_==timeint_afgenalpha)
        {
          // For af-generalized-alpha scheme, we already have the true residual,...
          trueresidual_->Update(1.0,*residual_,0.0);

          // ...but the residual vector for the solution rhs has to be scaled.
          residual_->Scale(gamma_*dta_/alphaM_);
        }
        else
        {
          // scaling to get true residual vector for all other schemes
          trueresidual_->Update(ResidualScaling(),*residual_,0.0);
        }

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
    dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), residual_);

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
    if (velnorm_L2 < 1e-5) velnorm_L2 = 1.0;
    if (prenorm_L2 < 1e-5) prenorm_L2 = 1.0;

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

    //--------- Apply Dirichlet boundary conditions to system of equations
    //          residual displacements are supposed to be zero at
    //          boundary conditions
    incvel_->PutScalar(0.0);
    {
      // time measurement: application of dbc
      TEUCHOS_FUNC_TIME_MONITOR("      + apply DBC");
      LINALG::ApplyDirichlettoSystem(sysmat_,incvel_,residual_,zeros_,*(dbcmaps_->CondMap()));
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

    // -------------------------------------------------------------------
    // update velocity and pressure values by increments
    // -------------------------------------------------------------------
    velnp_->Update(1.0,*incvel_,1.0);

    // -------------------------------------------------------------------
    // For af-generalized-alpha, also update accelerations, but do not
    // update time derivatives of density for low-Mach-number flow.
    // Furthermore, calculate velocities, pressures, densities,
    // accelerations and density time derivatives at intermediate time
    // steps n+alpha_F and n+alpha_M, respectively, for next iteration.
    // This has to be done at the end of the iteration, since we might
    // need the velocities at n+alpha_F in a potential coupling
    // algorithm, for instance.
    // -------------------------------------------------------------------
    if (timealgo_==timeint_afgenalpha)
    {
      // -------------------------------------------------------------------
      // separate velocity from pressure values
      velpressplitter_.ExtractOtherVector(incvel_,onlyvel);

      // ------------------------------------------------------
      // use updated velocity to update acceleration:
      //
      //    n+1         n+1
      // acc      =  acc    + (1/gamma*dt)*dvel
      //    (i+1)       (i)
      //
      onlyvel->Scale(1.0/(gamma_*dta_));
      // in case of conservative form: velocity*density
      if (convform_ == "conservative")
      {
        Teuchos::RCP<Epetra_Vector> onlydens = velpressplitter_.ExtractOtherVector(vedenp_);
        Teuchos::RCP<Epetra_Vector> denstimesvel = velpressplitter_.ExtractOtherVector(vedenp_);
        denstimesvel->Multiply(1.0, *onlyvel, *onlydens, 0.0);
        velpressplitter_.AddOtherVector(denstimesvel,accnp_);
      }
      else velpressplitter_.AddOtherVector(onlyvel,accnp_);

      //       n+alphaM                n+1                      n
      //    acc         = alpha_M * acc     + (1-alpha_M) *  acc
      //       (i)                     (i)
      accam_->Update((alphaM_),*accnp_,(1.0-alphaM_),*accn_,0.0);

      //       n+alphaF              n+1                   n
      //      u         = alpha_F * u     + (1-alpha_F) * u
      //       (i)                   (i)
      velaf_->Update((alphaF_),*velnp_,(1.0-alphaF_),*veln_,0.0);

      // values only required for low-Mach-number flow
      if (loma_ != "No")
      {
        vedeaf_->Update((alphaF_),*vedenp_,(1.0-alphaF_),*veden_,0.0);
        vedeam_->Update((alphaM_),*vedenp_,(1.0-alphaM_),*veden_,0.0);
      }
    }

    //------------------------------------------------ free surface update
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
            if (!ndnorm->Map().MyGID(rgid) or !fsvelnp->Map().MyGID(rgid)
                or ndnorm->Map().MyGID(rgid) != fsvelnp->Map().MyGID(rgid))
              dserror("Sparse vector does not have global row  %d or vectors don't match",rgid);
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
    eleparams.set("action","calc_linear_fluid");

    // other parameters that might be needed by the elements
    eleparams.set("total time",time_);
    eleparams.set("thsl",theta_*dta_);
    eleparams.set("low-Mach-number solver",loma_);
    eleparams.set("eos factor",eosfac_);

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

    LINALG::ApplyDirichlettoSystem(sysmat_,velnp_,rhs_,velnp_,*(dbcmaps_->CondMap()));
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
    // Take Dirichlet values from velnp and add vel to veln for non-Dirichlet
    // values.
    Teuchos::RCP<Epetra_Vector> aux = LINALG::CreateVector(*(discret_->DofRowMap()),true);
    aux->Update(1.0, *veln_, 1.0, *vel, 0.0);
    dbcmaps_->InsertOtherVector(dbcmaps_->ExtractOtherVector(aux), velnp_);
  }

  // add Neumann loads
  residual_->Update(1.0,*neumann_loads_,0.0);

  // update impedance boundary condition
  impedancebc_->UpdateResidual(residual_);

  // create the parameters for the discretization
  ParameterList eleparams;

  // set general element parameters
  eleparams.set("thsl",theta_*dta_);
  eleparams.set("dt",dta_);
  eleparams.set("form of convective term",convform_);
  eleparams.set("fs subgrid viscosity",fssgv_);
  eleparams.set("Linearisation",newton_);
  eleparams.set("low-Mach-number solver",loma_);
  eleparams.set("eos factor",eosfac_);

  // parameters for stabilization
  eleparams.sublist("STABILIZATION") = params_.sublist("STABILIZATION");

  // parameters for stabilization
  eleparams.sublist("TURBULENCE MODEL") = params_.sublist("TURBULENCE MODEL");

  // set general vector values needed by elements
  discret_->ClearState();
  discret_->SetState("hist",hist_ );
  if (alefluid_)
  {
    discret_->SetState("dispnp", dispnp_);
    discret_->SetState("gridv", gridv_);
  }

  // set scheme-specific element parameters and vector values
  if (timealgo_==timeint_afgenalpha)
  {
    eleparams.set("action","calc_fluid_afgenalpha_systemmat_and_residual");
    eleparams.set("using generalized-alpha time integration",true);
    eleparams.set("total time",time_-(1-alphaF_)*dta_);
    eleparams.set("timefacrhs",alphaM_/(gamma_*dta_));

    discret_->SetState("velnp", velaf_ );
    discret_->SetState("vedenp",vedeaf_);
    discret_->SetState("accam", accam_ );
    if (convform_ == "conservative") discret_->SetState("vedeam",vedenp_);
    else                             discret_->SetState("vedeam",vedeam_);
  }
  else
  {
    eleparams.set("action","calc_fluid_systemmat_and_residual");
    eleparams.set("using generalized-alpha time integration",false);
    eleparams.set("total time",time_);

    discret_->SetState("velnp",velnp_);
    discret_->SetState("vedenp",vedenp_);
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
    meshmovematrix_->ApplyDirichlet(*(dbcmaps_->CondMap()),false);
  }

  trueresidual_->Update(ResidualScaling(),*residual_,0.0);

  // Apply dirichlet boundary conditions to system of equations
  // residual displacements are supposed to be zero at boundary
  // conditions
  incvel_->PutScalar(0.0);
  LINALG::ApplyDirichlettoSystem(sysmat_,incvel_,residual_,zeros_,*(dbcmaps_->CondMap()));
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                                      |
 | One-step-Theta: (step>1)                                             |
 |                                                                      |
 |  accn_  = (velnp_-veln_) / (Theta * dt) - (1/Theta -1) * accn_"(n+1) |
 |                                                                      |
 |  velnm_ =veln_                                                       |
 |  veln_  =velnp_                                                      |
 |                                                                      |
 |  BDF2:           (step>1)                                            |
 |                                                                      |
 |               2*dt(n)+dt(n-1)              dt(n)+dt(n-1)             |
 |  accn_   = --------------------- velnp_ - --------------- veln_      |
 |            dt(n)*[dt(n)+dt(n-1)]           dt(n)*dt(n-1)             |
 |                                                                      |
 |                     dt(n)                                            |
 |           + ----------------------- velnm_                           |
 |             dt(n-1)*[dt(n)+dt(n-1)]                                  |
 |                                                                      |
 |  velnm_ =veln_                                                       |
 |  veln_  =velnp_                                                      |
 |                                                                      |
 |  BDF2 and  One-step-Theta: (step==1)                                 |
 |                                                                      |
 |  The given formulas are only valid from the second timestep. In the  |
 |  first step, the acceleration is calculated simply by                |
 |                                                                      |
 |  accn_  = (velnp_-veln_) / (dt)                                      |
 |                                                                      |
 |  For low-Mach-number flow, the same is done for density values,      |
 |  which are located at the "pressure dofs" of "vede"-vectors and all  |
 |  velocity values are multiplied by the respective density values.    |
 |                                                                      |
 |                                                           gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidImplicitTimeInt::TimeUpdate()
{
  if (loma_ != "No")
  {
    // compute accelerations and density time derivatives
    TIMEINT_THETA_BDF2::CalculateAcceleration(
        vedenp_, veden_, vedenm_, accn_,
        timealgo_, step_, theta_, dta_, dtp_,
        accnp_);
  }
  else
  {
    // compute accelerations
    TIMEINT_THETA_BDF2::CalculateAcceleration(
        velnp_, veln_, velnm_, accn_,
        timealgo_, step_, theta_, dta_, dtp_,
        accnp_);
  }

  // update old acceleration
  accn_->Update(1.0,*accnp_,0.0);

  // velocities/pressures of this step become most recent
  // velocities/pressures of the last step
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
 | lift'n'drag forces, statistics time sample and output of solution    |
 | and statistics                                              vg 11/08 |
 -----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidImplicitTimeInt::StatisticsAndOutput()
{
  // time measurement: output and statistics
  TEUCHOS_FUNC_TIME_MONITOR("      + output and statistics");

  // -------------------------------------------------------------------
  //          calculate lift'n'drag forces from the residual
  // -------------------------------------------------------------------
  LiftDrag();

  // -------------------------------------------------------------------
  //   add calculated velocity to mean value calculation (statistics)
  // -------------------------------------------------------------------
  statisticsmanager_->DoTimeSample(step_,time_,eosfac_);

  // -------------------------------------------------------------------
  //                         output of solution
  // -------------------------------------------------------------------
  Output();

  // -------------------------------------------------------------------
  //          dumping of turbulence statistics if required
  // -------------------------------------------------------------------
  statisticsmanager_->DoOutput(output_,step_);

  return;
} // FluidImplicitTimeInt::StatisticsAndOutput


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
 | prepare AVM3-based scale separation                         vg 10/08 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidImplicitTimeInt::AVM3Preparation()
{
  {// time measurement: avm3
  TEUCHOS_FUNC_TIME_MONITOR("           + avm3");

  // zero matrix
  sysmat_->Zero();

  // add Neumann loads
  residual_->Update(1.0,*neumann_loads_,0.0);

  // create the parameters for the discretization
  ParameterList eleparams;

  // set general element parameters
  eleparams.set("thsl",theta_*dta_);
  eleparams.set("dt",dta_);
  eleparams.set("form of convective term",convform_);
  eleparams.set("fs subgrid viscosity",fssgv_);
  eleparams.set("Linearisation",newton_);
  eleparams.set("low-Mach-number solver",loma_);
  eleparams.set("eos factor",eosfac_);

  // parameters for stabilization
  eleparams.sublist("STABILIZATION") = params_.sublist("STABILIZATION");

  // parameters for stabilization
  eleparams.sublist("TURBULENCE MODEL") = params_.sublist("TURBULENCE MODEL");

  // set general vector values needed by elements
  discret_->ClearState();
  discret_->SetState("hist"  ,hist_ );

  // zero and set fine-scale vector required by element routines
  fsvelnp_->PutScalar(0.0);
  discret_->SetState("fsvelnp",fsvelnp_);

  // set scheme-specific element parameters and vector values
  if (timealgo_==timeint_stationary)
  {
    eleparams.set("action","calc_fluid_stationary_systemmat_and_residual");
    eleparams.set("using generalized-alpha time integration",false);
    eleparams.set("total time",time_);

    discret_->SetState("velnp",velnp_);
    discret_->SetState("vedenp",vedenp_);
  }
  else if (timealgo_==timeint_afgenalpha)
  {
    eleparams.set("action","calc_fluid_afgenalpha_systemmat_and_residual");
    eleparams.set("using generalized-alpha time integration",true);
    eleparams.set("total time",time_-(1-alphaF_)*dta_);
    eleparams.set("timefacrhs",alphaM_/(gamma_*dta_));

    discret_->SetState("velnp", velaf_ );
    discret_->SetState("vedenp",vedeaf_);
    discret_->SetState("accam", accam_ );
    if (convform_ == "conservative") discret_->SetState("vedeam",vedenp_);
    else                             discret_->SetState("vedeam",vedeam_);
  }
  else
  {
    eleparams.set("action","calc_fluid_systemmat_and_residual");
    eleparams.set("using generalized-alpha time integration",false);
    eleparams.set("total time",time_);

    discret_->SetState("velnp",velnp_);
    discret_->SetState("vedenp",vedenp_);
  }

  // element evaluation for getting system matrix
  // -> we merely need matrix "structure" below, not the actual contents
  discret_->Evaluate(eleparams,sysmat_,residual_);
  discret_->ClearState();

  // complete system matrix
  sysmat_->Complete();

  // apply DBC to system matrix
  LINALG::ApplyDirichlettoSystem(sysmat_,incvel_,residual_,zeros_,*(dbcmaps_->CondMap()));

  // extract ML parameters
  ParameterList&  mllist = solver_.Params().sublist("ML Parameters");

  // call VM3 constructor with system matrix for generating scale-separating matrix
  {
    const Teuchos::RCP<const Epetra_Vector> dirichtoggle = Dirichlet();
    vm3_solver_ = rcp(new VM3_Solver(SystemMatrix(),dirichtoggle,mllist,true,true));
  }
  }// time measurement: avm3

  return;
}// FluidImplicitTimeInt::AVM3Preparation


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | AVM3-based scale separation                                 vg 10/08 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidImplicitTimeInt::AVM3Separation()
{
  {// time measurement: avm3
  TEUCHOS_FUNC_TIME_MONITOR("           + avm3");

  // check whether VM3 solver exists
  if (vm3_solver_ == null) dserror("vm3_solver not allocated");

  if (timealgo_==timeint_afgenalpha)
  {
    // call VM3 scale separation to get fine-scale part of velocity
    // at time n+alpha_F
    vm3_solver_->Separate(fsvelnp_,velaf_);

  }
  else
  {
    // call VM3 scale separation to get fine-scale part of velocity
    // at time n+1
    vm3_solver_->Separate(fsvelnp_,velnp_);
  }

  // set fine-scale vector
  discret_->SetState("fsvelnp",fsvelnp_);
  }// time measurement: avm3

  return;
}// FluidImplicitTimeInt::AVM3Separation


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
          map<int, vector<int> >::iterator master = pbcmapmastertoslave_->find(lnode->Id());

          // slavenodes are ignored
          if(master == pbcmapmastertoslave_->end()) continue;
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
   RCP<const Epetra_Vector> densnm,
   const double             eosfac)
{
  // get factor for equation of state (i.e., therm. press. / gas constant)
  eosfac_ = eosfac;

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

  vector<int>    Indices(numdof);
  vector<double> Values(numdof);

  // There are three different ways for filling the vedenp/veden-vectors:
  // 1) generalized-alpha/conservative: vel-dofs: density,  pre-dofs: density
  // 2) generalized-alpha/convective:   vel-dofs: velocity, pre-dofs: density
  // 3) one-step-theta/BDF2:    vel-dofs: velocity*density, pre-dofs: density
  // Furthermore, the accn-vector is filled with time derivatives of density
  // at pre-dof locations for generalized-alpha time integration (contained in
  // densnm-vector here).
  // Moreover, for BDF2, the same mentioned above for the vedenp/veden-vectors
  // is done for the vedenm-vector.
  if (timealgo_==timeint_afgenalpha)
  {
    if (convform_ == "conservative")
    {
      // loop all nodes on the processor
      for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
      {
        for(int index=0;index<numdim;++index)
        {
          Indices[index] = lnodeid*numdof + index;

          Values[index] = (*densnp)[lnodeid];
          vedenp_->ReplaceMyValues(1,&Values[index],&Indices[index]);

          Values[index] = (*densn)[lnodeid];
          veden_->ReplaceMyValues(1,&Values[index],&Indices[index]);
        }

        Indices[numdim] = lnodeid*numdof + numdim;

        Values[numdim] = (*densnp)[lnodeid];
        vedenp_->ReplaceMyValues(1,&Values[numdim],&Indices[numdim]);

        Values[numdim] = (*densn)[lnodeid];
        veden_->ReplaceMyValues(1,&Values[numdim],&Indices[numdim]);

        Values[numdim] = (*densnm)[lnodeid];
        accn_->ReplaceMyValues(1,&Values[numdim],&Indices[numdim]);
      }
    }
    else
    {
      // get velocity dofs for vede-vectors as copies from vel-vectors
      vedenp_->Update(1.0,*velnp_,0.0);
      veden_->Update(1.0,*veln_,0.0);

      // loop all nodes on the processor
      for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
      {
        Indices[numdim] = lnodeid*numdof + numdim;

        Values[numdim] = (*densnp)[lnodeid];
        vedenp_->ReplaceMyValues(1,&Values[numdim],&Indices[numdim]);

        Values[numdim] = (*densn)[lnodeid];
        veden_->ReplaceMyValues(1,&Values[numdim],&Indices[numdim]);

        Values[numdim] = (*densnm)[lnodeid];
        accn_->ReplaceMyValues(1,&Values[numdim],&Indices[numdim]);
      }
    }
  }
  else
  {
    // insert density values in vedenp- and veden-vectors
    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
    {
      for(int index=0;index<numdim;++index)
      {
        Indices[index] = lnodeid*numdof + index;

        Values[index] = (*densnp)[lnodeid]*(*velnp_)[Indices[index]];
        vedenp_->ReplaceMyValues(1,&Values[index],&Indices[index]);

        Values[index] = (*densn)[lnodeid]*(*veln_)[Indices[index]];
        veden_->ReplaceMyValues(1,&Values[index],&Indices[index]);
      }

      Indices[numdim] = lnodeid*numdof + numdim;

      Values[numdim] = (*densnp)[lnodeid];
      vedenp_->ReplaceMyValues(1,&Values[numdim],&Indices[numdim]);

      Values[numdim] = (*densn)[lnodeid];
      veden_->ReplaceMyValues(1,&Values[numdim],&Indices[numdim]);
    }

    if (timealgo_==timeint_bdf2)
    {
      // insert density values in vedenm-vector
      // loop all nodes on the processor
      for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
      {
        for(int index=0;index<numdim;++index)
        {
          Indices[index] = lnodeid*numdof + index;

          Values[index] = (*densnm)[lnodeid]*(*velnm_)[Indices[index]];
          vedenm_->ReplaceMyValues(1,&Values[index],&Indices[index]);
        }

        Indices[numdim] = lnodeid*numdof + numdim;

        Values[numdim] = (*densnm)[lnodeid];
        vedenm_->ReplaceMyValues(1,&Values[numdim],&Indices[numdim]);
      }
    }
  }

  return;

} // ScaTraTimIntImpl::SetTimeLomaFields


/*----------------------------------------------------------------------*
 | set outer-iteration-related fields for low-Mach-number flow vg 08/08 |
 *----------------------------------------------------------------------*/
void FLD::FluidImplicitTimeInt::SetIterLomaFields(
   RCP<const Epetra_Vector> densnp,
   RCP<const Epetra_Vector> densdtnp,
   const double             eosfac)
{
  // get factor for equation of state (i.e., therm. press. / gas constant)
  eosfac_ = eosfac;

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

  vector<int>    Indices(numdof);
  vector<double> Values(numdof);

  // There are three different ways for filling the vedenp-vector:
  // 1) generalized-alpha/conservative: vel-dofs: density,  pre-dofs: density
  // 2) generalized-alpha/convective:   vel-dofs: velocity, pre-dofs: density
  // 3) one-step-theta/BDF2:    vel-dofs: velocity*density, pre-dofs: density
  // Furthermore, the accnp-vector is filled with time derivatives of density
  // at pre-dof locations for generalized-alpha time integration
  if (timealgo_==timeint_afgenalpha)
  {
    if (convform_ == "conservative")
    {
      // insert density values in vedenp-vector
      // loop all nodes on the processor
      for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
      {
        for(int index=0;index<numdim;++index)
        {
          Indices[index] = lnodeid*numdof + index;

          Values[index] = (*densnp)[lnodeid];
          vedenp_->ReplaceMyValues(1,&Values[index],&Indices[index]);
        }

        Indices[numdim] = lnodeid*numdof + numdim;

        Values[numdim] = (*densnp)[lnodeid];
        vedenp_->ReplaceMyValues(1,&Values[numdim],&Indices[numdim]);

        Values[numdim] = (*densdtnp)[lnodeid];
        accnp_->ReplaceMyValues(1,&Values[numdim],&Indices[numdim]);
      }
    }
    else
    {
      // get velocity dofs for vedenp-vector as copies from velnp-vector
      vedenp_->Update(1.0,*velnp_,0.0);

      // insert density time derivative values in accn-vector
      // loop all nodes on the processor
      for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
      {
        Indices[numdim] = lnodeid*numdof + numdim;

        Values[numdim] = (*densnp)[lnodeid];
        vedenp_->ReplaceMyValues(1,&Values[numdim],&Indices[numdim]);

        Values[numdim] = (*densdtnp)[lnodeid];
        accnp_->ReplaceMyValues(1,&Values[numdim],&Indices[numdim]);
      }
    }
  }
  else
  {
    // insert density values in vedenp-vector
    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
    {
      for(int index=0;index<numdim;++index)
      {
        Indices[index] = lnodeid*numdof + index;

        Values[index] = (*densnp)[lnodeid]*(*velnp_)[Indices[index]];
        vedenp_->ReplaceMyValues(1,&Values[index],&Indices[index]);
      }

      Indices[numdim] = lnodeid*numdof + numdim;

      Values[numdim] = (*densnp)[lnodeid];
      vedenp_->ReplaceMyValues(1,&Values[numdim],&Indices[numdim]);
    }
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
      eleparams.set("form of convective term",convform_);
      eleparams.set("fs subgrid viscosity",fssgv_);

      // set vector values needed by elements
      discret_->ClearState();
      discret_->SetState("velnp",velnp_);
      // predicted dirichlet values
      // velnp then also holds prescribed new dirichlet values
      discret_->EvaluateDirichlet(eleparams,velnp_,null,null,null);
      discret_->ClearState();

      // evaluate Neumann b.c.
      //eleparams.set("inc_density",density_);

      discret_->SetState("vedenp",vedenp_);
      neumann_loads_->PutScalar(0.0);
      discret_->EvaluateNeumann(eleparams,*neumann_loads_);
      discret_->ClearState();
    }

    // -------------------------------------------------------------------
    //           preparation of AVM3-based scale separation
    // -------------------------------------------------------------------
    if (step_==1 and fssgv_ != "No") AVM3Preparation();

    // -------------------------------------------------------------------
    //                     solve nonlinear equation system
    // -------------------------------------------------------------------
    NonlinearSolve();

    // -------------------------------------------------------------------
    //         calculate lift'n'drag forces from the residual
    // -------------------------------------------------------------------
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

Notice: Angular moments obtained from lift&drag forces currently refere to the
        initial configuration, i.e. are built with the coordinates X of a particular
        node irrespective of its current position.
*/
void FLD::FluidImplicitTimeInt::LiftDrag() const
{
  // in this map, the results of the lift drag calculation are stored
  RCP<map<int,vector<double> > > liftdragvals;

  FLD::UTILS::LiftDrag(*discret_,*trueresidual_,params_,liftdragvals);

  if (liftdragvals!=Teuchos::null)
    FLD::UTILS::WriteLiftDragToFile(time_, step_, *liftdragvals);

  return;
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
  // perform filtering and computation of Cs
  {
    const Teuchos::RCP<const Epetra_Vector> dirichtoggle = Dirichlet();
    DynSmag_->ApplyFilterForDynamicComputationOfCs(velnp_,dirichtoggle);
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

    // dbcmaps_ has already been set up

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

    // set general element parameters
    eleparams.set("thsl",theta_*dta_);
    eleparams.set("dt",dta_);
    eleparams.set("Linearisation",newton_);
    eleparams.set("form of convective term",convform_);
    eleparams.set("low-Mach-number solver",loma_);
    eleparams.set("eos factor",eosfac_);

    // parameters for stabilization
    eleparams.sublist("STABILIZATION") = params_.sublist("STABILIZATION");

    // parameters for stabilization
    eleparams.sublist("TURBULENCE MODEL") = params_.sublist("TURBULENCE MODEL");

    // set general vector values needed by elements
    discret_->ClearState();
    discret_->SetState("hist",zeros_ );
    discret_->SetState("dispnp", griddisp);
    discret_->SetState("gridv", zeros_);

    // set scheme-specific element parameters and vector values
    if (timealgo_==timeint_afgenalpha)
    {
      eleparams.set("action","calc_fluid_afgenalpha_systemmat_and_residual");
      eleparams.set("using generalized-alpha time integration",true);
      eleparams.set("total time",time_-(1-alphaF_)*dta_);
      eleparams.set("timefacrhs",alphaM_/(gamma_*dta_));

      discret_->SetState("velnp", velaf_ );
      discret_->SetState("vedenp",vedeaf_);
      discret_->SetState("accam", accam_ );
      if (convform_ == "conservative") discret_->SetState("vedeam",vedenp_);
      else                             discret_->SetState("vedeam",vedeam_);
    }
    else
    {
      eleparams.set("action","calc_fluid_systemmat_and_residual");
      eleparams.set("using generalized-alpha time integration",false);
      eleparams.set("total time",time_);

      discret_->SetState("velnp",velnp_);
      discret_->SetState("vedenp",vedenp_);
    }

    // call loop over elements
    discret_->Evaluate(eleparams,sysmat_,meshmatrix_,residual_);
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
    dirichletlines_ = Teuchos::null;
    dirichletlines_ = SystemMatrix()->ExtractDirichletLines(*(dbcmaps_->CondMap()));
    sysmat_->ApplyDirichlet(*(dbcmaps_->CondMap()));
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
  LINALG::ApplyDirichlettoSystem(incvel_,residual_,relax,*(dbcmaps_->CondMap()));

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


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::FluidImplicitTimeInt::AddDirichCond(const Teuchos::RCP<const Epetra_Map> maptoadd)
{
  std::vector<Teuchos::RCP<const Epetra_Map> > condmaps;
  condmaps.push_back(maptoadd);
  condmaps.push_back(dbcmaps_->CondMap());
  Teuchos::RCP<Epetra_Map> condmerged = LINALG::MultiMapExtractor::MergeMaps(condmaps);
  *dbcmaps_ = LINALG::MapExtractor(*(discret_->DofRowMap()), condmerged);
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Teuchos::RCP<const Epetra_Vector> FLD::FluidImplicitTimeInt::Dirichlet()
{
  if (dbcmaps_ == Teuchos::null)
    dserror("Dirichlet map has not been allocated");
  Teuchos::RCP<Epetra_Vector> dirichones = LINALG::CreateVector(*(dbcmaps_->CondMap()),false);
  dirichones->PutScalar(1.0);
  Teuchos::RCP<Epetra_Vector> dirichtoggle = LINALG::CreateVector(*(discret_->DofRowMap()),true);
  dbcmaps_->InsertCondVector(dirichones, dirichtoggle);
  return dirichtoggle;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Teuchos::RCP<const Epetra_Vector> FLD::FluidImplicitTimeInt::InvDirichlet()
{
  if (dbcmaps_ == Teuchos::null)
    dserror("Dirichlet map has not been allocated");
  Teuchos::RCP<Epetra_Vector> dirichzeros = LINALG::CreateVector(*(dbcmaps_->CondMap()),true);
  Teuchos::RCP<Epetra_Vector> invtoggle = LINALG::CreateVector(*(discret_->DofRowMap()),false);
  invtoggle->PutScalar(1.0);
  dbcmaps_->InsertCondVector(dirichzeros, invtoggle);
  return invtoggle;
}


#endif /* CCADISCRET       */
