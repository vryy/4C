/*!------------------------------------------------------------------------------------------------*
\file topopt_fluidAdjointImplTimeIntegration.cpp

\brief Control routine for fluid adjoint (in)stationary solvers,

     including instationary solvers based on

     o one-step-theta time-integration scheme

     o two-step BDF2 time-integration scheme
       (with potential one-step-theta start algorithm)

     and stationary solver.

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#ifdef CCADISCRET

// include Gmsh output
#define GMSHOUTPUT

#include "topopt_fluidAdjointImplTimeIntegration.H"
#include "../drt_fluid/time_integration_scheme.H"
#include "../linalg/linalg_ana.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_mapextractor.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_fluid/fluid_utils.H"
#include "../drt_fluid/fluidimpedancecondition.H"
#include "../drt_fluid/fluid_volumetric_surfaceFlow_condition.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"
#include "../drt_adapter/adapter_coupling_mortar.H"
#include "../drt_adapter/adapter_topopt.H"
#include "../drt_nurbs_discret/drt_apply_nurbs_initial_condition.H"
#include "../drt_opti/topopt_optimizer.H"

// print error file for function EvaluateErrorComparedToAnalyticalSol()
#include "../drt_io/io_control.H"

// for AVM3 solver:
#include <MLAPI_Workspace.h>
#include <MLAPI_Aggregation.h>

#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_gmsh.H"

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Constructor (public)                                     gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
TOPOPT::ADJOINT::ImplicitTimeInt::ImplicitTimeInt(RefCountPtr<DRT::Discretization> actdis,
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
  step_(0),
  extrapolationpredictor_(params.get("do explicit predictor",true)),
  uprestart_(params.get("write restart every", -1)),
  upres_(params.get("write solution every", -1)),
  surfacesplitter_(NULL)
{
  // -------------------------------------------------------------------
  // get the processor ID from the communicator
  // -------------------------------------------------------------------
  myrank_  = discret_->Comm().MyPID();

  // time measurement: initialization
  TEUCHOS_FUNC_TIME_MONITOR(" + adjoint initialization");

  // -------------------------------------------------------------------
  // get the basic parameters first
  // -------------------------------------------------------------------

  // type of time-integration
  timealgo_ = DRT::INPUT::get<INPAR::FLUID::TimeIntegrationScheme>(params_, "time int algo");
  //genalpha integration scheme (afgenalpha or npgenalpha)
  if (timealgo_==INPAR::FLUID::timeint_afgenalpha or timealgo_==INPAR::FLUID::timeint_npgenalpha)
    is_genalpha_= true;
  else
    is_genalpha_= false;
  // time-step size
  dtp_ = dta_ = params_.get<double>("time step size");
  // maximum number of timesteps
  stepmax_  = params_.get<int>   ("max number timesteps");
  // maximum simulation time
  maxtime_  = params_.get<double>("total time");
  // initial simulation time
  time_ = maxtime_; // adjoint equations start at endtime
  // parameter theta for time-integration schemes
  theta_    = params_.get<double>("theta");
  // compute or set 1.0 - theta for time-integration schemes
  if (timealgo_ == INPAR::FLUID::timeint_one_step_theta)  omtheta_ = 1.0 - theta_;
  else                                      omtheta_ = 0.0;
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
    if (timealgo_ != INPAR::FLUID::timeint_afgenalpha)
      dserror("no starting algorithm supported for schemes other than af-gen-alpha");
    else startalgo_= true;
    if (numstasteps_>stepmax_)
      dserror("more steps for starting algorithm than steps overall");
  }

  // parameter for linearization scheme (fixed-point-like or Newton)
  newton_ = DRT::INPUT::get<INPAR::FLUID::LinearisationAction>(params_, "Linearisation");

  // use of specific predictor
  // (might be used for af-generalized-alpha, but not yet activated)

  if(params_.get<string>("predictor","disabled") == "disabled")
  {
    if(myrank_==0)
    {
      printf("disabled extrapolation predictor\n\n");
    }
    extrapolationpredictor_=false;
  }

  // form of convective term
  convform_ = params_.get<string>("form of convective term","convective");

  // -------------------------------------------------------------------
  // account for potential Neuman inflow terms if required
  // -------------------------------------------------------------------
  neumanninflow_ = false;
  if (params_.get<string>("Neumann inflow","no") == "yes") neumanninflow_ = true;

  // -------------------------------------------------------------------
  // care for periodic boundary conditions
  // -------------------------------------------------------------------
  pbcmapmastertoslave_ = params_.get<RCP<map<int,vector<int> > > >("periodic bc");
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
  numdim_ = params_.get<int>("number of velocity degrees of freedom");

  FLD::UTILS::SetupFluidSplit(*discret_,numdim_,velpressplitter_);

  // -------------------------------------------------------------------
  // create empty system matrix --- stiffness and mass are assembled in
  // one system matrix!
  // -------------------------------------------------------------------

  // This is a first estimate for the number of non zeros in a row of
  // the matrix. Assuming a structured 3d-fluid mesh we have 27 adjacent
  // nodes with 4 dofs each. (27*4=108)
  // We do not need the exact number here, just for performance reasons
  // a 'good' estimate

  if (not params_.get<int>("Simple Preconditioner",0) &&
      not params_.get<int>("AMG BS Preconditioner",0))
  {
    // initialize standard (stabilized) system matrix
    sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,108,false,true));
  }
  else
  {
    Teuchos::RCP<LINALG::BlockSparseMatrix<FLD::UTILS::VelPressSplitStrategy> > blocksysmat =
      Teuchos::rcp(new LINALG::BlockSparseMatrix<FLD::UTILS::VelPressSplitStrategy>(velpressplitter_,velpressplitter_,108,false,true));
    blocksysmat->SetNumdim(numdim_);
    sysmat_ = blocksysmat;
  }

  // -------------------------------------------------------------------
  // create empty vectors
  // -------------------------------------------------------------------

  // Vectors passed to the element
  // -----------------------------
  // velocity/pressure at time n+1, n and n-1
  velnp_ = LINALG::CreateVector(*dofrowmap,true);
  veln_  = LINALG::CreateVector(*dofrowmap,true);
  velnm_ = LINALG::CreateVector(*dofrowmap,true);

  // acceleration/(scalar time derivative) at time n+1 and n
  accnp_ = LINALG::CreateVector(*dofrowmap,true);
  accn_  = LINALG::CreateVector(*dofrowmap,true);

  // velocity/pressure at time n+alpha_F
  velaf_ = LINALG::CreateVector(*dofrowmap,true);

  // acceleration/(scalar time derivative) at time n+alpha_M/(n+alpha_M/n)
  accam_ = LINALG::CreateVector(*dofrowmap,true);

  // history vector
  hist_ = LINALG::CreateVector(*dofrowmap,true);

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
  }

  // the vector containing body and surface forces
  neumann_loads_= LINALG::CreateVector(*dofrowmap,true);

  // Vectors used for solution process
  // ---------------------------------
  // rhs: standard (stabilized) residual vector (rhs for the incremental form)
  residual_      = LINALG::CreateVector(*dofrowmap,true);
  trueresidual_  = LINALG::CreateVector(*dofrowmap,true);

  // right hand side vector for linearised solution;
  rhs_ = LINALG::CreateVector(*dofrowmap,true);

  // initialize pseudo-porosity vector for topology optimization as null
  topopt_porosity_ = Teuchos::null;

  // ---------------------------------------------------------------------
  // set general fluid parameter defined before
  // ---------------------------------------------------------------------
  SetElementGeneralAdjointParameter();
  SetElementTimeParameter();
} // ImplicitTimeInt::ImplicitTimeInt


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
void TOPOPT::ADJOINT::ImplicitTimeInt::Integrate()
{
  // output of stabilization details
  if (myrank_==0)
  {
    ParameterList *  stabparams=&(params_.sublist("STABILIZATION"));

    cout << "Stabilization type         : " << stabparams->get<string>("STABTYPE") << "\n";
    cout << "                             " << stabparams->get<string>("TDS")<< "\n";
    cout << "\n";
    cout << "                             " << "Tau Type        = " << stabparams->get<string>("DEFINITION_TAU") <<"\n";
    cout << "                             " << "Evaluation Tau  = " << stabparams->get<string>("EVALUATION_TAU") <<"\n";
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
    cout << endl;
    cout << "                             " << "Evaluation Mat  = " << stabparams->get<string>("EVALUATION_MAT") <<"\n";
    cout << "\n";
  }

  // distinguish stationary and instationary case
  if (timealgo_==INPAR::FLUID::timeint_stationary) SolveStationaryProblem();
  else                                             TimeLoop();

  // print the results of time measurements
  TimeMonitor::summarize();

  return;
} // ImplicitTimeInt::Integrate



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | contains the time loop                                    gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void TOPOPT::ADJOINT::ImplicitTimeInt::TimeLoop()
{
  // time measurement: time loop
  TEUCHOS_FUNC_TIME_MONITOR(" + adjoint time loop");

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
      case INPAR::FLUID::timeint_one_step_theta:
        printf("TIME: %11.4E/%11.4E  DT = %11.4E   One-Step-Theta (%0.2f)   STEP = %4d/%4d \n",
              time_,maxtime_,dta_,theta_,step_,stepmax_);
        break;
      case INPAR::FLUID::timeint_afgenalpha:
        printf("TIME: %11.4E/%11.4E  DT = %11.4E  Af-Generalized-Alpha  STEP = %4d/%4d \n",
               time_,maxtime_,dta_,step_,stepmax_);
        break;
      case INPAR::FLUID::timeint_npgenalpha:
        printf("TIME: %11.4E/%11.4E  DT = %11.4E  Np-Generalized-Alpha  STEP = %4d/%4d \n",
               time_,maxtime_,dta_,step_,stepmax_);
        break;
      case INPAR::FLUID::timeint_bdf2:
        printf("TIME: %11.4E/%11.4E  DT = %11.4E       BDF2          STEP = %4d/%4d \n",
               time_,maxtime_,dta_,step_,stepmax_);
        break;
      default:
        dserror("parameter out of range: IOP\n");
      } /* end of switch(timealgo) */
    }

    // -----------------------------------------------------------------
    //                     solve nonlinear equation
    // -----------------------------------------------------------------
    NonLinearSolve();

    // -------------------------------------------------------------------
    //                         update solution
    //        current solution becomes old solution of next timestep
    // -------------------------------------------------------------------
    TimeUpdate();


    // -------------------------------------------------------------------
    //  output
    // -------------------------------------------------------------------
    Output();

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
} // ImplicitTimeInt::TimeLoop


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | setup the variables to do a new time step                 u.kue 06/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void TOPOPT::ADJOINT::ImplicitTimeInt::PrepareTimeStep()
{

  // -------------------------------------------------------------------
  //              set time dependent parameters
  // -------------------------------------------------------------------
  step_ += 1;
  time_ += dta_;

  // for BDF2, theta is set by the time-step sizes, 2/3 for const. dt
  if (timealgo_==INPAR::FLUID::timeint_bdf2) theta_ = (dta_+dtp_)/(2.0*dta_ + dtp_);

  // -------------------------------------------------------------------
  //  For af-generalized-alpha time-integration scheme:
  //  set "pseudo-theta", calculate initial accelerations according to
  //  prescribed Dirichlet values for generalized-alpha time
  //  integration and values at intermediate time steps
  // -------------------------------------------------------------------
  if (is_genalpha_)
  {
    // starting algorithm
    if (startalgo_)
    {
      // use backward-Euler-type parameter combination
      if (step_<=numstasteps_)
      {
        if (myrank_==0)
        {
          cout<<"Starting algorithm for Af_GenAlpha active."
              <<"Performing step "<<step_ <<" of "<<numstasteps_
              <<" Backward Euler starting steps"<<endl;
        }
        alphaM_ = 1.0;
        alphaF_ = 1.0;
        gamma_  = 1.0;
      }
      else
      {
        // recall original user wish
        alphaM_ = params_.get<double>("alpha_M");
        alphaF_ = params_.get<double>("alpha_F");
        gamma_  = params_.get<double>("gamma");
        // do not enter starting algorithm section in the future
        startalgo_ = false;
      }
    }

    // compute "pseudo-theta" for af-generalized-alpha scheme
    theta_ = alphaF_*gamma_/alphaM_;
  }

  // -------------------------------------------------------------------
  //  Set time parameter for element call
  // -------------------------------------------------------------------

  SetElementTimeParameter();

  // -------------------------------------------------------------------
  //  evaluate Dirichlet and Neumann boundary conditions
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

    // evaluate Neumann conditions
    neumann_loads_->PutScalar(0.0);
    discret_->EvaluateNeumann(eleparams,*neumann_loads_);
    discret_->ClearState();
  }

  if (is_genalpha_)
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
    GenAlphaUpdateAcceleration();

    // ----------------------------------------------------------------
    // compute values at intermediate time steps
    // ----------------------------------------------------------------
    GenAlphaIntermediateValues();
  }

  return;
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
void TOPOPT::ADJOINT::ImplicitTimeInt::NonLinearSolve()
{
  dirichletlines_ = Teuchos::null;
  // Do not remove meshmatrix_ here as we want to reuse its graph.
  // (We pay for the memory anyway if we use it, we might as well keep it.)
  //meshmatrix_ = Teuchos::null;

  // time measurement: nonlinear iteration
  TEUCHOS_FUNC_TIME_MONITOR("   + lin. iteration/lin. adjoint solve");

  dtsolve_  = 0.0;
  dtele_    = 0.0;

  if (myrank_ == 0)
  {
    printf("+------------+-------------------+--------------+--------------+--------------+--------------+\n");
    printf("|- step/max -|- tol      [norm] -|-- vel-res ---|-- pre-res ---|-- vel-inc ---|-- pre-inc ---|\n");
  }

  // -------------------------------------------------------------------
  // call elements to calculate system matrix and RHS
  // -------------------------------------------------------------------
  {
    // time measurement: element
    TEUCHOS_FUNC_TIME_MONITOR("      + adjoint element calls");

    // get cpu time
    const double tcpu=Teuchos::Time::wallTime();

    sysmat_->Zero();

    // create the parameters for the discretization
    ParameterList eleparams;

    eleparams.set("action","calc_fluid_systemmat_and_residual");

    // add Neumann loads
    residual_->Update(1.0,*neumann_loads_,0.0);


    discret_->ClearState();
    discret_->SetState("velnp",velnp_);
    discret_->SetState("hist",hist_);

    //set additional pseudo-porosity field for topology optimization
    eleparams.set("topopt_porosity",topopt_porosity_);

    // set general vector values needed by elements
    discret_->ClearState();
    discret_->SetState("hist" ,hist_ );
    discret_->SetState("accam",accam_);

    // set scheme-specific element parameters and vector values
    if (is_genalpha_)
    {
      discret_->SetState("velaf",velaf_);
      if (timealgo_==INPAR::FLUID::timeint_npgenalpha)
        discret_->SetState("velnp",velnp_);
    }
    else discret_->SetState("velaf",velnp_);

    // call standard loop over elements
    discret_->Evaluate(eleparams,sysmat_,null,residual_,null,null);

    discret_->ClearState();

    //----------------------------------------------------------------------
    // account for potential Neumann inflow terms
    //----------------------------------------------------------------------
    if (neumanninflow_)
    {
      // create parameter list
      ParameterList condparams;

      // action for elements
      condparams.set("action","calc_Neumann_inflow");

      // set vector values needed by elements
      discret_->ClearState();
      // set scheme-specific element parameters and vector values
      if (is_genalpha_)
        discret_->SetState("velaf",velaf_);
      // there is't any contribution of the pressure or the continuity equation
      // in the case of Neumann inflow
      // -> there is no difference between af_genalpha and np_genalpha
      else discret_->SetState("velaf",velnp_);

      std::string condstring("FluidNeumannInflow");
      discret_->EvaluateCondition(condparams,sysmat_,Teuchos::null,residual_,Teuchos::null,Teuchos::null,condstring);
      discret_->ClearState();
    }

    // scaling to get true residual vector
    trueresidual_->Update(ResidualScaling(),*residual_,0.0);

    // finalize the complete matrix
    sysmat_->Complete();

    // end time measurement for element
    dtele_=Teuchos::Time::wallTime()-tcpu;
  }

  // blank residual DOFs which are on Dirichlet BC
  // We can do this because the values at the dirichlet positions
  // are not used anyway.
  // We could avoid this though, if velrowmap_ and prerowmap_ would
  // not include the dirichlet values as well. But it is expensive
  // to avoid that.
  dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), residual_);

  Teuchos::RCP<Epetra_Vector> onlyvel = velpressplitter_.ExtractOtherVector(residual_);
  onlyvel->Norm2(&vresnorm_);

  velpressplitter_.ExtractOtherVector(velnp_,onlyvel);
  onlyvel->Norm2(&velnorm_L2_);

  Teuchos::RCP<Epetra_Vector> onlypre = velpressplitter_.ExtractCondVector(residual_);
  onlypre->Norm2(&presnorm_);

  velpressplitter_.ExtractCondVector(velnp_,onlypre);
  onlypre->Norm2(&prenorm_L2_);

  // care for the case that nothing really happens in the velocity
  // or pressure field
  if (velnorm_L2_ < 1e-5) velnorm_L2_ = 1.0;
  if (prenorm_L2_ < 1e-5) prenorm_L2_ = 1.0;

  //-------------------------------------------------- output to screen
  printf("| %10.3E   | %10.3E   | %10.3E   | %10.3E   |",
      vresnorm_,presnorm_,
      incvelnorm_L2_/velnorm_L2_,incprenorm_L2_/prenorm_L2_);
  printf(" (      --     ,te=%10.3E",dtele_);
  printf("+--------------+--------------+--------------+--------------+\n");

  FILE* errfile = params_.get<FILE*>("err file",NULL);
  if (errfile!=NULL)
  {
    fprintf(errfile,"fluid solve: vres=%10.3E  pres=%10.3E  vinc=%10.3E  pinc=%10.3E\n",
        vresnorm_,presnorm_,
        incvelnorm_L2_/velnorm_L2_,incprenorm_L2_/prenorm_L2_);
  }

  //--------- Apply Dirichlet boundary conditions to system of equations
  //          residual displacements are supposed to be zero at
  //          boundary conditions
  {
    // time measurement: application of dbc
    TEUCHOS_FUNC_TIME_MONITOR("      + apply adjoint DBC");
    LINALG::ApplyDirichlettoSystem(sysmat_,velnp_,residual_,zeros_,*(dbcmaps_->CondMap()));
  }

  //-------solve for residual displacements to correct incremental displacements
  {
    // time measurement: solver
    TEUCHOS_FUNC_TIME_MONITOR("      + adjoint solver calls");

    // get cpu time
    const double tcpusolve=Teuchos::Time::wallTime();

    solver_.Solve(sysmat_->EpetraOperator(),velnp_,residual_,true,1, w_, c_, project_);

    solver_.ResetTolerance();

    // end time measurement for solver
    dtsolve_ = Teuchos::Time::wallTime()-tcpusolve;
  }


  // -------------------------------------------------------------------
  // For af-generalized-alpha: update accelerations
  // Furthermore, calculate velocities, pressures, scalars and
  // accelerations at intermediate time steps n+alpha_F and n+alpha_M,
  // respectively, for next iteration.
  // This has to be done at the end of the iteration, since we might
  // need the velocities at n+alpha_F in a potential coupling
  // algorithm, for instance.
  // -------------------------------------------------------------------
  if (is_genalpha_)
  {
    GenAlphaUpdateAcceleration();

    GenAlphaIntermediateValues();
  }
} // ImplicitTimeInt::LinearSolve



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | predictor                                                   vg 02/09 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void TOPOPT::ADJOINT::ImplicitTimeInt::Predictor()
{
  // -------------------------------------------------------------------
  // time measurement: nonlinear iteration
  // -------------------------------------------------------------------
  TEUCHOS_FUNC_TIME_MONITOR("   + adjoint predictor");

  // -------------------------------------------------------------------
  // call elements to calculate system matrix and rhs and assemble
  // -------------------------------------------------------------------
  AssembleMatAndRHS();

  // -------------------------------------------------------------------
  // calculate and print out residual norms
  // (blank residual DOFs which are on Dirichlet BC
  // We can do this because the values at the dirichlet positions
  // are not used anyway.
  // We could avoid this though, if velrowmap_ and prerowmap_ would
  // not include the dirichlet values as well. But it is expensive
  // to avoid that.)
  // -------------------------------------------------------------------
  dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_),residual_);

  // -------------------------------------------------------------------
  // take surface volumetric flow rate into account
  //    RCP<Epetra_Vector> temp_vec = rcp(new Epetra_Vector(*vol_surf_flow_bcmaps_,true));
  //    vol_surf_flow_bc_->InsertCondVector( *temp_vec , *residual_);
  // -------------------------------------------------------------------
  vol_flow_rates_bc_extractor_->InsertVolumetricSurfaceFlowCondVector(
    vol_flow_rates_bc_extractor_->ExtractVolumetricSurfaceFlowCondVector(zeros_),
    residual_);

  // -------------------------------------------------------------------
  // extract velocity and pressure and compute respective norms
  // -------------------------------------------------------------------
  Teuchos::RCP<Epetra_Vector> onlyvel = velpressplitter_.ExtractOtherVector(residual_);
  Teuchos::RCP<Epetra_Vector> onlypre = velpressplitter_.ExtractCondVector(residual_);

  onlyvel->Norm2(&vresnorm_);
  onlypre->Norm2(&presnorm_);

  // -------------------------------------------------------------------
  // print out norms (only first processor)
  // -------------------------------------------------------------------
  if (myrank_ == 0)
  {
    printf("+------------+-------------------+--------------+--------------+--------------+--------------+\n");
    printf("| predictor  |  vel. | pre. res. | %10.3E   | %10.3E   |      --      |      --      |",vresnorm_,presnorm_);
    printf(" (      --     ,te=%10.3E)\n",dtele_);
    printf("+------------+-------------------+--------------+--------------+--------------+--------------+\n");
  }

} // ImplicitTimeInt::Predictor


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | preparatives for solver                                     vg 09/11 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void TOPOPT::ADJOINT::ImplicitTimeInt::PrepareSolve()
{
  // call elements to calculate system matrix and rhs and assemble
  AssembleMatAndRHS();

  // apply Dirichlet boundary conditions to system of equations
  ApplyDirichletToSystem();

} // ImplicitTimeInt::PrepareSolve


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | call elements to calculate system matrix/rhs and assemble   vg 02/09 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void TOPOPT::ADJOINT::ImplicitTimeInt::AssembleMatAndRHS()
{
  dtele_    = 0.0;

  // time measurement: element
  TEUCHOS_FUNC_TIME_MONITOR("      + adjoint element calls");

  // get cpu time
  const double tcpu=Teuchos::Time::wallTime();

  sysmat_->Zero();

  // create the parameters for the discretization
  ParameterList eleparams;

  // add Neumann loads
  residual_->Update(1.0,*neumann_loads_,0.0);

  discret_->ClearState();
  discret_->SetState("velnp",velnp_);
  discret_->SetState("hist",hist_);

  discret_->ClearState();


  // set action type
  eleparams.set("action","calc_fluid_systemmat_and_residual");

  // set general vector values needed by elements
  discret_->ClearState();
  discret_->SetState("hist" ,hist_ );
  discret_->SetState("accam",accam_);

  // set scheme-specific element parameters and vector values
  if (is_genalpha_)
  {
    discret_->SetState("velaf",velaf_);
    if (timealgo_==INPAR::FLUID::timeint_npgenalpha)
      discret_->SetState("velnp",velnp_);
  }
  else discret_->SetState("velaf",velnp_);

  // call standard loop over elements
  discret_->Evaluate(eleparams,sysmat_,null,residual_,null,null);
  discret_->ClearState();

  // account for potential Neumann inflow terms
  if (neumanninflow_)
  {
    // create parameter list
    ParameterList condparams;

    // action for elements
    condparams.set("action","calc_Neumann_inflow");

    // set vector values needed by elements
    discret_->ClearState();

    // set scheme-specific element parameters and vector values
    if (is_genalpha_)
      discret_->SetState("velaf",velaf_);
      // there is't any contribution of the pressure or the continuity equation
      // in the case of Neumann inflow
      // -> there is no difference between af_genalpha and np_genalpha
    else discret_->SetState("velaf",velnp_);

    std::string condstring("FluidNeumannInflow");
    discret_->EvaluateCondition(condparams,sysmat_,Teuchos::null,residual_,Teuchos::null,Teuchos::null,condstring);
    discret_->ClearState();
  }

  // scaling to get true residual vector
  trueresidual_->Update(ResidualScaling(),*residual_,0.0);

  // finalize the complete matrix
  sysmat_->Complete();

  // end time measurement for element
  dtele_=Teuchos::Time::wallTime()-tcpu;

} // ImplicitTimeInt::AssembleMatAndRHS


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | application of Dirichlet boundary conditions to system      vg 09/11 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void TOPOPT::ADJOINT::ImplicitTimeInt::ApplyDirichletToSystem()
{
  // -------------------------------------------------------------------
  // apply Dirichlet boundary conditions to system of equations:
  // - Residual displacements are supposed to be zero for resp. dofs.
  // - Time for applying Dirichlet boundary conditions is measured.
  // -------------------------------------------------------------------
  velnp_->PutScalar(0.0);
  {
    TEUCHOS_FUNC_TIME_MONITOR("      + apply DBC");
    LINALG::ApplyDirichlettoSystem(sysmat_,velnp_,residual_,zeros_,*(dbcmaps_->CondMap()));
  }
} // ImplicitTimeInt::ApplyDirichletToSystem


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | update within iteration                                     vg 09/11 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void TOPOPT::ADJOINT::ImplicitTimeInt::IterUpdate()
{

  // -------------------------------------------------------------------
  // For af-generalized-alpha: update accelerations
  // Furthermore, calculate velocities, pressures, scalars and
  // accelerations at intermediate time steps n+alpha_F and n+alpha_M,
  // respectively, for next iteration.
  // This has to be done at the end of the iteration, since we might
  // need the velocities at n+alpha_F in a potential coupling
  // algorithm, for instance.
  // -------------------------------------------------------------------
  if (is_genalpha_)
  {
    GenAlphaUpdateAcceleration();

    GenAlphaIntermediateValues();
  }

} // ImplicitTimeInt::IterUpdate


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | update acceleration for generalized-alpha time integration  vg 02/09 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void TOPOPT::ADJOINT::ImplicitTimeInt::GenAlphaUpdateAcceleration()
{

  //                                  n+1     n
  //                               vel   - vel
  //       n+1      n  gamma-1.0      (i)
  //    acc    = acc * --------- + ------------
  //       (i)           gamma      gamma * dt
  //

  // extract the degrees of freedom associated with velocities
  // only these are allowed to be updated, otherwise you will
  // run into trouble in loma, where the 'pressure' component
  // is used to store the acceleration of the temperature
  Teuchos::RCP<Epetra_Vector> onlyaccn  = velpressplitter_.ExtractOtherVector(accn_ );
  Teuchos::RCP<Epetra_Vector> onlyveln  = velpressplitter_.ExtractOtherVector(veln_ );
  Teuchos::RCP<Epetra_Vector> onlyvelnp = velpressplitter_.ExtractOtherVector(velnp_);

  Teuchos::RCP<Epetra_Vector> onlyaccnp = rcp(new Epetra_Vector(onlyaccn->Map()));

  const double fact1 = 1.0/(gamma_*dta_);
  const double fact2 = 1.0 - (1.0/gamma_);
  onlyaccnp->Update(fact2,*onlyaccn,0.0);
  onlyaccnp->Update(fact1,*onlyvelnp,-fact1,*onlyveln,1.0);

  // copy back into global vector
  LINALG::Export(*onlyaccnp,*accnp_);

} // ImplicitTimeInt::GenAlphaUpdateAcceleration


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | compute values at intermediate time steps for gen.-alpha    vg 02/09 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void TOPOPT::ADJOINT::ImplicitTimeInt::GenAlphaIntermediateValues()
{
  //       n+alphaM                n+1                      n
  //    acc         = alpha_M * acc     + (1-alpha_M) *  acc
  //       (i)                     (i)
  {
    // extract the degrees of freedom associated with velocities
    // only these are allowed to be updated, otherwise you will
    // run into trouble in loma, where the 'pressure' component
    // is used to store the acceleration of the temperature
    Teuchos::RCP<Epetra_Vector> onlyaccn  = velpressplitter_.ExtractOtherVector(accn_ );
    Teuchos::RCP<Epetra_Vector> onlyaccnp = velpressplitter_.ExtractOtherVector(accnp_);

    Teuchos::RCP<Epetra_Vector> onlyaccam = rcp(new Epetra_Vector(onlyaccnp->Map()));

    onlyaccam->Update((alphaM_),*onlyaccnp,(1.0-alphaM_),*onlyaccn,0.0);

    // copy back into global vector
    LINALG::Export(*onlyaccam,*accam_);
  }

  // set intermediate values for velocity
  //
  //       n+alphaF              n+1                   n
  //      u         = alpha_F * u     + (1-alpha_F) * u
  //       (i)                   (i)
  //
  // and pressure
  //
  //       n+alphaF              n+1                   n
  //      p         = alpha_F * p     + (1-alpha_F) * p
  //       (i)                   (i)
  //
  // note that its af-genalpha with mid-point treatment of the pressure,
  // not implicit treatment as for the genalpha according to Whiting
  velaf_->Update((alphaF_),*velnp_,(1.0-alphaF_),*veln_,0.0);

  if (timealgo_==INPAR::FLUID::timeint_npgenalpha)
  {
    // set intermediate values for velocity
    //
    //       n+alphaF              n+1                   n
    //      u         = alpha_F * u     + (1-alpha_F) * u
    //       (i)                   (i)
    //
    // and pressure
    //
    //       n+1
    //      p
    //       (i)
    //
    // note that its af-genalpha with mid-point treatment of the pressure,
    // not implicit treatment as for the genalpha according to Whiting
    Teuchos::RCP<Epetra_Vector> onlypre = velpressplitter_.ExtractCondVector(velnp_);

    LINALG::Export(*onlypre, *velaf_);
  }

} // ImplicitTimeInt::GenAlphaIntermediateValues


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | convergence check                                           vg 09/11 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void TOPOPT::ADJOINT::ImplicitTimeInt::ConvergenceCheck()
{
  // -------------------------------------------------------------------
  // calculate and print out norms for convergence check
  // (blank residual DOFs which are on Dirichlet BC
  // We can do this because the values at the dirichlet positions
  // are not used anyway.
  // We could avoid this though, if velrowmap_ and prerowmap_ would
  // not include the dirichlet values as well. But it is expensive
  // to avoid that.)
  // -------------------------------------------------------------------
  dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_),residual_);

  // -------------------------------------------------------------------
  // take surface volumetric flow rate into account
  //    RCP<Epetra_Vector> temp_vec = rcp(new Epetra_Vector(*vol_surf_flow_bcmaps_,true));
  //    vol_surf_flow_bc_->InsertCondVector( *temp_vec , *residual_);
  // -------------------------------------------------------------------
  vol_flow_rates_bc_extractor_->InsertVolumetricSurfaceFlowCondVector(
    vol_flow_rates_bc_extractor_->ExtractVolumetricSurfaceFlowCondVector(zeros_),
      residual_);

  Teuchos::RCP<Epetra_Vector> onlyvel = velpressplitter_.ExtractOtherVector(residual_);
  onlyvel->Norm2(&vresnorm_);

  velpressplitter_.ExtractOtherVector(velnp_,onlyvel);
  onlyvel->Norm2(&velnorm_L2_);

  Teuchos::RCP<Epetra_Vector> onlypre = velpressplitter_.ExtractCondVector(residual_);
  onlypre->Norm2(&presnorm_);

  velpressplitter_.ExtractCondVector(velnp_,onlypre);
  onlypre->Norm2(&prenorm_L2_);

  // check for any INF's and NaN's
  if (std::isnan(vresnorm_) or
      std::isnan(incvelnorm_L2_) or
      std::isnan(velnorm_L2_) or
      std::isnan(presnorm_) or
      std::isnan(incprenorm_L2_) or
      std::isnan(prenorm_L2_))
    dserror("At least one of the calculated vector norms is NaN.");

  if (abs(std::isinf(vresnorm_)) or
      abs(std::isinf(incvelnorm_L2_)) or
      abs(std::isinf(velnorm_L2_)) or
      abs(std::isinf(presnorm_)) or
      abs(std::isinf(incprenorm_L2_)) or
      abs(std::isinf(prenorm_L2_)))
    dserror("At least one of the calculated vector norms is INF.");

  // care for the case that nothing really happens in velocity
  // or pressure field
  if (velnorm_L2_ < 1e-5) velnorm_L2_ = 1.0;
  if (prenorm_L2_ < 1e-5) prenorm_L2_ = 1.0;

  if (myrank_ == 0)
  {
    printf("| %10.3E   | %10.3E   | %10.3E   | %10.3E   |",
      vresnorm_,presnorm_,incvelnorm_L2_/velnorm_L2_,incprenorm_L2_/prenorm_L2_);
    printf(" (ts=%10.3E,te=%10.3E)\n",dtsolve_,dtele_);
    printf("+--------------+--------------+--------------+--------------+\n");
    FILE* errfile = params_.get<FILE*>("err file",NULL);
    if (errfile!=NULL)
    {
      fprintf(errfile,"fluid solve:   vres=%10.3E  pres=%10.3E  vinc=%10.3E  pinc=%10.3E\n",
        vresnorm_,presnorm_,incvelnorm_L2_/velnorm_L2_,incprenorm_L2_/prenorm_L2_);
    }
  }
} // ImplicitTimeInt::ConvergenceCheck


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | build linear system matrix and rhs                        u.kue 11/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void TOPOPT::ADJOINT::ImplicitTimeInt::Evaluate(Teuchos::RCP<const Epetra_Vector> vel)
{
  sysmat_->Zero();

  // set the new solution we just got
  if (vel!=Teuchos::null)
  {
    // Take Dirichlet values from velnp and add vel to veln for non-Dirichlet
    // values.
    Teuchos::RCP<Epetra_Vector> aux = LINALG::CreateVector(*(discret_->DofRowMap()),true);
    aux->Update(1.0, *veln_, 1.0, *vel, 0.0);
    //    dbcmaps_->InsertOtherVector(dbcmaps_->ExtractOtherVector(aux), velnp_);
    dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(velnp_), aux);

    *velnp_ = *aux;

  }

  // add Neumann loads
  residual_->Update(1.0,*neumann_loads_,0.0);

  discret_->ClearState();
  discret_->SetState("velnp",velnp_);
  discret_->SetState("hist",hist_);

  discret_->ClearState();

  // create the parameters for the discretization
  ParameterList eleparams;

  // set general vector values needed by elements
  discret_->ClearState();
  discret_->SetState("hist" ,hist_ );
  discret_->SetState("accam",accam_);

  // set the only required state vectors
  if (is_genalpha_)
  {
    discret_->SetState("velaf",velaf_);
    if (timealgo_==INPAR::FLUID::timeint_npgenalpha)
      discret_->SetState("velnp",velnp_);
  }
  else discret_->SetState("velaf",velnp_);

  // call loop over elements
  discret_->Evaluate(eleparams,sysmat_,residual_);
  discret_->ClearState();

  // finalize the system matrix
  sysmat_->Complete();

  trueresidual_->Update(ResidualScaling(),*residual_,0.0);

  // Apply dirichlet boundary conditions to system of equations
  // residual displacements are supposed to be zero at boundary
  // conditions
  velnp_->PutScalar(0.0);
  LINALG::ApplyDirichlettoSystem(sysmat_,velnp_,residual_,zeros_,*(dbcmaps_->CondMap()));
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
 |                                                           gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void TOPOPT::ADJOINT::ImplicitTimeInt::TimeUpdate()
{

  ParameterList *  stabparams=&(params_.sublist("STABILIZATION"));

  if(stabparams->get<string>("TDS") == "time_dependent")
  {
    const double tcpu=Teuchos::Time::wallTime();

    if(myrank_==0)
    {
      cout << "time update for subscales";
    }

    // call elements to calculate system matrix and rhs and assemble
    // this is required for the time update of the subgrid scales and
    // makes sure that the current subgrid scales correspond to the
    // current residual
    AssembleMatAndRHS();

    // create the parameters for the discretization
    ParameterList eleparams;
    // action for elements
    eleparams.set("action","time update for subscales");

    // update time paramters
    if (is_genalpha_)
    {
      eleparams.set("gamma"  ,gamma_);
    }
    else if (timealgo_==INPAR::FLUID::timeint_one_step_theta)
    {
      eleparams.set("gamma"  ,theta_);
    }
    else if((timealgo_==INPAR::FLUID::timeint_bdf2))
    {
      eleparams.set("gamma"  ,1.0);
    }

    eleparams.set("dt"     ,dta_    );

    // call loop over elements to update subgrid scales
    discret_->Evaluate(eleparams,null,null,null,null,null);

    if(myrank_==0)
    {
      cout << "("<<Teuchos::Time::wallTime()-tcpu<<")\n";
    }
  }

  // compute accelerations
  {
    Teuchos::RCP<Epetra_Vector> onlyaccn  = velpressplitter_.ExtractOtherVector(accn_ );
    Teuchos::RCP<Epetra_Vector> onlyaccnp = velpressplitter_.ExtractOtherVector(accnp_);
    Teuchos::RCP<Epetra_Vector> onlyvelnm = velpressplitter_.ExtractOtherVector(velnm_);
    Teuchos::RCP<Epetra_Vector> onlyveln  = velpressplitter_.ExtractOtherVector(veln_ );
    Teuchos::RCP<Epetra_Vector> onlyvelnp = velpressplitter_.ExtractOtherVector(velnp_);

    FLD::TIMEINT_THETA_BDF2::CalculateAcceleration(onlyvelnp,
                                              onlyveln ,
                                              onlyvelnm,
                                              onlyaccn ,
                                              timealgo_,
                                              step_    ,
                                              theta_   ,
                                              dta_     ,
                                              dtp_     ,
                                              onlyaccnp);

    // copy back into global vector
    LINALG::Export(*onlyaccnp,*accnp_);
  }

  // update old acceleration
  accn_->Update(1.0,*accnp_,0.0);

  // velocities/pressures of this step become most recent
  // velocities/pressures of the last step
  velnm_->Update(1.0,*veln_ ,0.0);
  veln_ ->Update(1.0,*velnp_,0.0);

  // -------------------------------------------------------------------
  // treat impedance BC
  // note: these methods return without action, if the problem does not
  //       have impedance boundary conditions
  // -------------------------------------------------------------------
  discret_->ClearState();
  discret_->SetState("velnp",velnp_);
  discret_->SetState("hist",hist_);

  discret_->ClearState();

  return;
}// ImplicitTimeInt::TimeUpdate



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | output of solution vector to binio                        gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void TOPOPT::ADJOINT::ImplicitTimeInt::Output()
{

  //  ART_exp_timeInt_->Output();
  // output of solution
  if (step_%upres_ == 0)
  {
    // step number and time
    output_.NewStep(step_,time_);

    // velocity/pressure vector
    output_.WriteVector("velnp",velnp_);
    output_.WriteVector("tract_resid",trac_residual_);
    output_.WriteVector("neumann_loads",neumann_loads_);
    // (hydrodynamic) pressure
    Teuchos::RCP<Epetra_Vector> pressure = velpressplitter_.ExtractCondVector(velnp_);
    output_.WriteVector("pressure", pressure);

#ifdef GMSHOUTPUT
    OutputToGmsh(step_, time_,false);
#endif

    // write domain decomposition for visualization (only once!)
    if (step_==upres_) output_.WriteElementData();

    if (uprestart_ != 0 && step_%uprestart_ == 0) //add restart data
    {
      // acceleration vector at time n+1 and n, velocity/pressure vector at time n and n-1
      output_.WriteVector("accnp",accnp_);
      output_.WriteVector("accn", accn_);
      output_.WriteVector("veln", veln_);
      output_.WriteVector("velnm",velnm_);
    }
  }
  // write restart also when uprestart_ is not a integer multiple of upres_
  else if (uprestart_ > 0 && step_%uprestart_ == 0)
  {
    // step number and time
    output_.NewStep(step_,time_);

    // velocity/pressure vector
    output_.WriteVector("velnp",velnp_);

    // acceleration vector at time n+1 and n, velocity/pressure vector at time n and n-1
    output_.WriteVector("accnp",accnp_);
    output_.WriteVector("accn", accn_);
    output_.WriteVector("veln", veln_);
    output_.WriteVector("velnm",velnm_);
    output_.WriteVector("neumann_loads",neumann_loads_);
  }

  if (topopt_porosity_!=Teuchos::null)
    optimizer_->ImportAdjointFluidData(velnp_,step_);

  return;
} // ImplicitTimeInt::Output


void TOPOPT::ADJOINT::ImplicitTimeInt::OutputToGmsh(
    const int step,
    const double time,
    const bool inflow
    ) const
{
  // turn on/off screen output for writing process of Gmsh postprocessing file
  const bool screen_out = true;

  // create Gmsh postprocessing file
  // 20 steps are kept
  std::string filename = "dummy";
  if (inflow)
  {
    filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("solution_velpres_inflow", step, 20, screen_out, discret_->Comm().MyPID());
    //std::ofstream gmshfilecontent(filename.c_str());
  }
  else
  {
    filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("solution_velpres", step, 20, screen_out, discret_->Comm().MyPID());
    //std::ofstream gmshfilecontent(filename.c_str());
  }
  std::ofstream gmshfilecontent(filename.c_str());

  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" " << "velocity solution \" {" << endl;
    IO::GMSH::VelocityPressureFieldDofBasedToGmsh(discret_, velnp_ , "velocity", gmshfilecontent);
    gmshfilecontent << "};" << endl;
  }
  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" " << "pressure solution\" {" << endl;
    IO::GMSH::VelocityPressureFieldDofBasedToGmsh(discret_, velnp_, "pressure",gmshfilecontent);
    gmshfilecontent << "};" << endl;
  }

  gmshfilecontent.close();
  if (screen_out) std::cout << " done" << endl;

 return;
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |                                                             kue 04/07|
 -----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void TOPOPT::ADJOINT::ImplicitTimeInt::ReadRestart(int step)
{
  //  ART_exp_timeInt_->ReadRestart(step);
  IO::DiscretizationReader reader(discret_,step);
  time_ = reader.ReadDouble("time");
  step_ = reader.ReadInt("step");

  reader.ReadVector(velnp_,"velnp");
  reader.ReadVector(veln_, "veln");
  reader.ReadVector(velnm_,"velnm");
  reader.ReadVector(accnp_,"accnp");
  reader.ReadVector(accn_ ,"accn");

  // set element time parameter after restart:
  // Here it is already needed by AVM3 and impedance boundary condition!!
  SetElementTimeParameter();
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |set restart values (turbulent inflow only)             rasthofer 06/11|
 -----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void TOPOPT::ADJOINT::ImplicitTimeInt::SetRestart(
  const int step,
  const double time,
  Teuchos::RCP<const Epetra_Vector> readvelnp,
  Teuchos::RCP<const Epetra_Vector> readveln,
  Teuchos::RCP<const Epetra_Vector> readvelnm,
  Teuchos::RCP<const Epetra_Vector> readaccnp,
  Teuchos::RCP<const Epetra_Vector> readaccn)
{
  time_ = time;
  step_ = step;

  velnp_->Update(1.0,*readvelnp,0.0);
  veln_->Update(1.0,*readveln,0.0);
  velnm_->Update(1.0,*readvelnm,0.0);
  accnp_->Update(1.0,*readaccnp,0.0);
  accn_->Update(1.0,*readaccn,0.0);

  SetElementTimeParameter();
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
void TOPOPT::ADJOINT::ImplicitTimeInt::SetInitialFlowField(
  const INPAR::FLUID::InitialField initfield,
  const int startfuncno
  )
{
  // initial field by (undisturbed) function (init==2)
  // or disturbed function (init==3)
  if (initfield == INPAR::FLUID::initfield_field_by_function or
      initfield == INPAR::FLUID::initfield_disturbed_field_from_function)
  {
    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = discret_->lRowNode(lnodeid);
      // the set of degrees of freedom associated with the node
      const vector<int> nodedofset = discret_->Dof(lnode);

      for(int index=0;index<numdim_+1;++index)
      {
        int gid = nodedofset[index];

        double initialval=DRT::Problem::Instance()->Funct(startfuncno-1).Evaluate(index,lnode->X(),0.0,NULL);

        velnp_->ReplaceGlobalValues(1,&initialval,&gid);
      }
    }

    // for NURBS discretizations we have to solve a least squares problem,
    // with high accuracy! (do nothing for Lagrangian polynomials)
      DRT::NURBS::apply_nurbs_initial_condition(
        *discret_  ,
        DRT::Problem::Instance()->ErrorFile()->Handle(),
        DRT::Problem::Instance()->FluidSolverParams(),
        startfuncno,
        velnp_     );

    // initialize veln_ as well. That's what we actually want to do here!
    veln_->Update(1.0,*velnp_ ,0.0);

    // add random perturbation of certain percentage to function
    if (initfield == INPAR::FLUID::initfield_disturbed_field_from_function)
    {
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

        for(int index=0;index<numdim_;++index)
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
        for(int index=0;index<numdim_;++index)
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
  // special initial function: two counter-rotating vortices (2-D) and flame front
  // for flame-vortex interaction problem
  else if (initfield == INPAR::FLUID::initfield_flame_vortex_interaction)
  {
    const Epetra_Map* dofrowmap = discret_->DofRowMap();

    int err = 0;

    // define vectors for velocity field, node coordinates and coordinates
    // of left and right vortex
    vector<double> u(numdim_);
    vector<double> xy(numdim_);
    vector<double> xy0_left(numdim_);
    vector<double> xy0_right(numdim_);

    // check whether present flow is indeed two-dimensional
    if (numdim_!=2) dserror("Counter-rotating vortices are a two-dimensional flow!");

    // set laminar burning velocity, vortex strength C (scaled by laminar
    // burning velocity and (squared) vortex radius R
    const double sl = 1.0;
    const double C = 70.0*sl;
    const double R_squared = 16.0;

    // set density in unburnt and burnt phase and initialize actual density
    const double densu = 1.161;
    // -> for "pure fluid" computation: rhob = rhou = 1.161
    //const double densb = 1.161;
    const double densb = 0.157;
    double dens = 1.161;

    // initialize progress variable
    double pv = 0.0;

    // variables for evaluation of progress-variable profile
    // locations separating region 1 from region 2 and region 2 from region 3
    const double loc12 = 98.5;
    const double loc23 = 103.0;

    // define parameters for region 1 (exponential function for curve fitting)
    const double beta1  = 1.65;
    const double delta1 = 1.0;
    const double trans1 = 100.0;

    // define parameters for region 2 (linear function for curve fitting)
    const double abs2 = 0.0879;
    const double fac2 = 0.139309333;
    const double trans2 = 98.5;

    // define parameters for region 3 (exponential function for curve fitting)
    const double beta3  = 3.506209;
    const double delta3 = 4.28875;
    const double trans3 = 103.0;

    // set (scaled) vortex strength C, (squared) vortex radius R and define variables
    double r_squared_left;
    double r_squared_right;

    // set initial locations of vortices
    xy0_left[0] = 37.5;
    xy0_left[1] = 75.0;
    xy0_right[0] = 62.5;
    xy0_right[1] = 75.0;

    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = discret_->lRowNode(lnodeid);

      // the set of degrees of freedom associated with the node
      vector<int> nodedofset = discret_->Dof(lnode);

      // set node coordinates
      for(int dim=0;dim<numdim_;dim++)
      {
        xy[dim]=lnode->X()[dim];
      }

      // compute preliminary values for both vortices
      r_squared_left  = ((xy[0]-xy0_left[0])*(xy[0]-xy0_left[0])
                        +(xy[1]-xy0_left[1])*(xy[1]-xy0_left[1]))/R_squared;
      r_squared_right = ((xy[0]-xy0_right[0])*(xy[0]-xy0_right[0])
                        +(xy[1]-xy0_right[1])*(xy[1]-xy0_right[1]))/R_squared;

      // compute value of progress variable
      if (xy[1] < loc12-EPS10)
        pv = (1.0-(1.0/beta1))*exp((xy[1]-trans1)/delta1);
      else if (xy[1] > loc23+EPS10)
        pv = 1.0-(exp((1.0-beta3)*(xy[1]-trans3)/delta3)/beta3);
      else
        pv = fac2*(xy[1]-trans2) + abs2;

      // compute current density
      dens = densu+(densb-densu)*pv;

      // compute initial velocity components
      // including initial velocity distribution velocity in x2-direction
      u[0] = (C/R_squared)*(-(xy[1]-xy0_left[1])*exp(-r_squared_left/2.0)
                            +(xy[1]-xy0_right[1])*exp(-r_squared_right/2.0));
      u[1] = (C/R_squared)*( (xy[0]-xy0_left[0])*exp(-r_squared_left/2.0)
                            -(xy[0]-xy0_right[0])*exp(-r_squared_right/2.0))
                            + sl*densu/dens;

      // velocity profile due to flame without vortices:
      //u[1] = sl*densu/dens;

      // set initial velocity components
      for(int nveldof=0;nveldof<numdim_;nveldof++)
      {
        const int gid = nodedofset[nveldof];
        int lid = dofrowmap->LID(gid);
        err += velnp_->ReplaceMyValues(1,&(u[nveldof]),&lid);
        err += veln_ ->ReplaceMyValues(1,&(u[nveldof]),&lid);
        err += velnm_->ReplaceMyValues(1,&(u[nveldof]),&lid);
      }
    } // end loop nodes lnodeid

    if (err!=0) dserror("dof not on proc");
  }
  // special initial function: Beltrami flow (3-D)
  else if (initfield == INPAR::FLUID::initfield_beltrami_flow)
  {
    const Epetra_Map* dofrowmap = discret_->DofRowMap();

    int err = 0;

    const int npredof = numdim_;

    double         p;
    vector<double> u  (numdim_);
    vector<double> xyz(numdim_);

    // check whether present flow is indeed three-dimensional
    if (numdim_!=3) dserror("Beltrami flow is a three-dimensional flow!");

    // set constants for analytical solution
    const double a = M_PI/4.0;
    const double d = M_PI/2.0;

    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = discret_->lRowNode(lnodeid);

      // the set of degrees of freedom associated with the node
      vector<int> nodedofset = discret_->Dof(lnode);

      // set node coordinates
      for(int dim=0;dim<numdim_;dim++)
      {
        xyz[dim]=lnode->X()[dim];
      }

      // compute initial velocity components
      u[0] = -a * ( exp(a*xyz[0]) * sin(a*xyz[1] + d*xyz[2]) +
                    exp(a*xyz[2]) * cos(a*xyz[0] + d*xyz[1]) );
      u[1] = -a * ( exp(a*xyz[1]) * sin(a*xyz[2] + d*xyz[0]) +
                    exp(a*xyz[0]) * cos(a*xyz[1] + d*xyz[2]) );
      u[2] = -a * ( exp(a*xyz[2]) * sin(a*xyz[0] + d*xyz[1]) +
                    exp(a*xyz[1]) * cos(a*xyz[2] + d*xyz[0]) );

      // compute initial pressure
      p = -a*a/2.0 *
        ( exp(2.0*a*xyz[0])
          + exp(2.0*a*xyz[1])
          + exp(2.0*a*xyz[2])
          + 2.0 * sin(a*xyz[0] + d*xyz[1]) * cos(a*xyz[2] + d*xyz[0]) * exp(a*(xyz[1]+xyz[2]))
          + 2.0 * sin(a*xyz[1] + d*xyz[2]) * cos(a*xyz[0] + d*xyz[1]) * exp(a*(xyz[2]+xyz[0]))
          + 2.0 * sin(a*xyz[2] + d*xyz[0]) * cos(a*xyz[1] + d*xyz[2]) * exp(a*(xyz[0]+xyz[1]))
          );

      // set initial velocity components
      for(int nveldof=0;nveldof<numdim_;nveldof++)
      {
        const int gid = nodedofset[nveldof];
        int lid = dofrowmap->LID(gid);
        err += velnp_->ReplaceMyValues(1,&(u[nveldof]),&lid);
        err += veln_ ->ReplaceMyValues(1,&(u[nveldof]),&lid);
        err += velnm_->ReplaceMyValues(1,&(u[nveldof]),&lid);
      }

      // set initial pressure
      const int gid = nodedofset[npredof];
      int lid = dofrowmap->LID(gid);
      err += velnp_->ReplaceMyValues(1,&p,&lid);
      err += veln_ ->ReplaceMyValues(1,&p,&lid);
      err += velnm_->ReplaceMyValues(1,&p,&lid);
    } // end loop nodes lnodeid

    if (err!=0) dserror("dof not on proc");
  }
  // special initial function: test case due to Bochev et al. (2007) (2-D)
  else if (initfield == INPAR::FLUID::initfield_bochev_test)
  {
    const Epetra_Map* dofrowmap = discret_->DofRowMap();

    int err = 0;

    // check whether present flow is indeed two-dimensional
    if (numdim_!=2) dserror("Bochev test case is a two-dimensional flow!");

    // define vectors for velocity and pressure field as well as node coordinates
    vector<double> up(numdim_+1);
    vector<double> xy(numdim_);

    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = discret_->lRowNode(lnodeid);

      // the set of degrees of freedom associated with the node
      vector<int> nodedofset = discret_->Dof(lnode);

      // set node coordinates
      for(int dim=0;dim<numdim_;dim++)
      {
        xy[dim]=lnode->X()[dim];
      }

      // compute initial velocity and pressure components
      up[0] = sin(M_PI*xy[0]-0.7)*sin(M_PI*xy[1]+0.2);
      up[1] = cos(M_PI*xy[0]-0.7)*cos(M_PI*xy[1]+0.2);
      up[2] = sin(xy[0])*cos(xy[1])+(cos(1.0)-1.0)*sin(1.0);

      // set initial velocity and pressure components
      for(int ndof=0;ndof<numdim_+1;ndof++)
      {
        const int gid = nodedofset[ndof];
        int lid = dofrowmap->LID(gid);
        err += velnp_->ReplaceMyValues(1,&(up[ndof]),&lid);
        err += veln_ ->ReplaceMyValues(1,&(up[ndof]),&lid);
        err += velnm_->ReplaceMyValues(1,&(up[ndof]),&lid);
      }
    } // end loop nodes lnodeid

    if (err!=0) dserror("dof not on proc");
  }
  else
  {
    dserror("Only initial fields auch as a zero field, initial fields by (un-)disturbed functions and three special initial fields (counter-rotating vortices, Beltrami flow and Bochev test) are available up to now!");
  }

  return;
} // end SetInitialFlowField



/*----------------------------------------------------------------------*
 | sent density field for topology optimization         winklmaier 12/11|
 *----------------------------------------------------------------------*/
void TOPOPT::ADJOINT::ImplicitTimeInt::SetTopOptData(
    RCP<Epetra_Vector> porosity,
    RCP<TOPOPT::Optimizer> optimizer
)
{
  topopt_porosity_ = porosity;
  optimizer_=optimizer;
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
void TOPOPT::ADJOINT::ImplicitTimeInt::EvaluateErrorComparedToAnalyticalSol()
{

  INPAR::FLUID::CalcError calcerr = DRT::INPUT::get<INPAR::FLUID::CalcError>(params_,"calculate error");

  switch(calcerr)
  {
  case INPAR::FLUID::no_error_calculation:
    // do nothing --- no analytical solution available
    break;
  case INPAR::FLUID::beltrami_flow:
  case INPAR::FLUID::channel2D:
  case INPAR::FLUID::gravitation:
  case INPAR::FLUID::shear_flow:
  {
    // create the parameters for the discretization
    ParameterList eleparams;

    // action for elements
    eleparams.set("action","calc_fluid_error");
    eleparams.set<int>("calculate error",calcerr);

    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("u and p at time n+1 (converged)",velnp_);

    // get (squared) error values
    // 0: vel_mag
    // 1: p
    // 2: u_mag,analytical
    // 3: p_analytic
    // (4: vel_x)
    // (5: vel_y)
    // (6: vel_z)
    Teuchos::RCP<Epetra_SerialDenseVector> errors
      = Teuchos::rcp(new Epetra_SerialDenseVector(2+2));
    //  = Teuchos::rcp(new Epetra_SerialDenseVector(numdim_+2+2))

    // call loop over elements (assemble nothing)
    discret_->EvaluateScalars(eleparams, errors);
    discret_->ClearState();

    double velerr = 0.0;
    double preerr = 0.0;

    // integrated analytic solution in order to compute relative error
    double velint = 0.0;
    double pint = 0.0;

    // error in the single velocity components
    //double velerrx = 0.0;
    //double velerry = 0.0;
    //double velerrz = 0.0;

    // for the L2 norm, we need the square root
    velerr = sqrt((*errors)[0]);
    preerr = sqrt((*errors)[1]);

    // analytical vel_mag and p_mag
    velint= sqrt((*errors)[2]);
    pint = sqrt((*errors)[3]);

    if (myrank_ == 0)
    {
      {
        cout.precision(8);
        cout << endl << "----relative L_2 error norm for analytical solution Nr. " <<
          DRT::INPUT::get<INPAR::FLUID::CalcError>(params_,"calculate error") <<
          " ----------" << endl;
        cout << "| velocity:  " << velerr/velint << endl;
        cout << "| pressure:  " << preerr/pint << endl;
        cout << "--------------------------------------------------------------------" << endl << endl;
      }

      //velerrx = sqrt((*errors)[4]);
      //velerry = sqrt((*errors)[5]);
      //if (numdim_==3)
      //  velerrz = sqrt((*errors)[6]);

      // append error of the last time step to the error file
      if ((step_==stepmax_) or (time_==maxtime_))// write results to file
      {
        ostringstream temp;
        const std::string simulation = DRT::Problem::Instance()->OutputControlFile()->FileName();
        const std::string fname = simulation+".relerror";

        std::ofstream f;
        f.open(fname.c_str(),std::fstream::ate | std::fstream::app);
        f << "#| " << simulation << "\n";
        f << "#| Step | Time | rel. L2-error velocity mag |  rel. L2-error pressure  |\n";
        f << step_ << " " << time_ << " " << velerr/velint << " " << preerr/pint << " "<<"\n";
        f.flush();
        f.close();
      }

      ostringstream temp;
      const std::string simulation = DRT::Problem::Instance()->OutputControlFile()->FileName();
      const std::string fname = simulation+"_time.relerror";

      if(step_==1)
      {
        std::ofstream f;
        f.open(fname.c_str());
        f << "#| Step | Time | rel. L2-error velocity mag |  rel. L2-error pressure  |\n";
        f << step_ << " " << time_ << " " << velerr/velint << " " << preerr/pint << " "<<"\n";
        f.flush();
        f.close();
      }
      else
      {
        std::ofstream f;
        f.open(fname.c_str(),std::fstream::ate | std::fstream::app);
        f << step_ << " " << time_ << " " << velerr/velint << " " << preerr/pint << " "<<"\n";
        f.flush();
        f.close();
      }
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
void TOPOPT::ADJOINT::ImplicitTimeInt::SolveStationaryProblem()
{
  // time measurement: time loop (stationary) --- start TimeMonitor tm2
  TEUCHOS_FUNC_TIME_MONITOR(" + adjoint time loop");

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
    printf("Stationary Fluid Adjoint Solver - STEP = %4d/%4d \n",step_,stepmax_);
   }

    SetElementTimeParameter();

    // -------------------------------------------------------------------
    //         evaluate Dirichlet and Neumann boundary conditions
    // -------------------------------------------------------------------
    {
      ParameterList eleparams;

      // other parameters needed by the elements
      eleparams.set("total time",time_);

      // set vector values needed by elements
      discret_->ClearState();
      discret_->SetState("velaf",velnp_);
      // predicted dirichlet values
      // velnp then also holds prescribed new dirichlet values
      discret_->EvaluateDirichlet(eleparams,velnp_,null,null,null);

      discret_->ClearState();

      neumann_loads_->PutScalar(0.0);
      discret_->EvaluateNeumann(eleparams,*neumann_loads_);
      discret_->ClearState();
    }

    // -------------------------------------------------------------------
    //                     solve nonlinear equation system
    // -------------------------------------------------------------------
    NonLinearSolve();

    // -------------------------------------------------------------------
    //                         output of solution
    // -------------------------------------------------------------------
    Output();
  } // end of time loop
} // ImplicitTimeInt::SolveStationaryProblem



/*----------------------------------------------------------------------*
 | Destructor dtor (public)                                  gammi 04/07|
 *----------------------------------------------------------------------*/
TOPOPT::ADJOINT::ImplicitTimeInt::~ImplicitTimeInt()
{
  return;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void TOPOPT::ADJOINT::ImplicitTimeInt::AddDirichCond(const Teuchos::RCP<const Epetra_Map> maptoadd)
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
const Teuchos::RCP<const Epetra_Vector> TOPOPT::ADJOINT::ImplicitTimeInt::Dirichlet()
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
const Teuchos::RCP<const Epetra_Vector> TOPOPT::ADJOINT::ImplicitTimeInt::InvDirichlet()
{
  if (dbcmaps_ == Teuchos::null)
    dserror("Dirichlet map has not been allocated");
  Teuchos::RCP<Epetra_Vector> dirichzeros = LINALG::CreateVector(*(dbcmaps_->CondMap()),true);
  Teuchos::RCP<Epetra_Vector> invtoggle = LINALG::CreateVector(*(discret_->DofRowMap()),false);
  invtoggle->PutScalar(1.0);
  dbcmaps_->InsertCondVector(dirichzeros, invtoggle);
  return invtoggle;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> TOPOPT::ADJOINT::ImplicitTimeInt::VelocityRowMap()
{ return velpressplitter_.OtherMap(); }

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> TOPOPT::ADJOINT::ImplicitTimeInt::PressureRowMap()
{ return velpressplitter_.CondMap(); }



// -------------------------------------------------------------------
// set general fluid parameter (AE 01/2011)
// -------------------------------------------------------------------

void TOPOPT::ADJOINT::ImplicitTimeInt::SetElementGeneralAdjointParameter()
{
  ParameterList eleparams;

  eleparams.set("action","set_general_fluid_parameter");

  // set general element parameters
  eleparams.set("form of convective term",convform_);

  // parameter for stabilization
  eleparams.sublist("STABILIZATION") = params_.sublist("STABILIZATION");

  //set time integration scheme
  eleparams.set<int>("TimeIntegrationScheme", timealgo_);

  // call standard loop over elements TODO uncomment
//  discret_->Evaluate(eleparams,null,null,null,null,null);
  return;
}


// -------------------------------------------------------------------
// set general time parameter (AE 01/2011)
// -------------------------------------------------------------------

void TOPOPT::ADJOINT::ImplicitTimeInt::SetElementTimeParameter()
{
  ParameterList eleparams;

  eleparams.set("action","set_time_parameter");

  // set general element parameters
  eleparams.set("dt",dta_);
  eleparams.set("theta",theta_);
  eleparams.set("omtheta",omtheta_);

  // set scheme-specific element parameters and vector values
  if (timealgo_==INPAR::FLUID::timeint_stationary)
  {
    eleparams.set("total time",time_);
  }
  else if (is_genalpha_)
  {
    if (time_ < maxtime_)
    {dserror("put here time+... instead of time-... is this correct?");
      eleparams.set("total time",time_+(1-alphaF_)*dta_);
    }
    else
    {
      eleparams.set("total time",time_);
    }
    eleparams.set("alphaF",alphaF_);
    eleparams.set("alphaM",alphaM_);
    eleparams.set("gamma",gamma_);
  }
  else
  {
    eleparams.set("total time",time_);
  }

  // call standard loop over elements
  discret_->Evaluate(eleparams,null,null,null,null,null);
  return;
}
#endif  // #ifdef CCADISCRET
