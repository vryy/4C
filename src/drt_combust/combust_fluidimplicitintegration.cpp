/*!----------------------------------------------------------------------
\file combust_fluidimplicitintegration.cpp
\brief class holding implicit time integration schemes for combustion problems

This class is a merger of the standard fluid time integration and the XFSI time integration classes.
member functions. Maybe it will not be kept as a stand-alone class until the end of days, but
unified with a generalized XFEM time integration class.

For the time being, the only available time integration scheme for combustion problems is the
One-step-theta scheme.

Since a combustion problem is always a coupled multi-field problem, this class is remote-controlled
by the combustion algorithm. It does not have a TimeLoop() on its own. This class is only in charge
of finding the solution to the fluid field in a nonlinear iterative procedure in the context of a
combustion problem.

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>

*----------------------------------------------------------------------*/

#include "combust_fluidimplicitintegration.H"
#include "combust3_interpolation.H"
#include "combust_defines.H"
#include "combust_flamefront.H"
#include "combust_fluidresulttest.H"
#include "two_phase_defines.H"

#include "../drt_fluid/time_integration_scheme.H"
#include "../drt_fluid/fluid_utils.H"
#include "../drt_fluid/drt_periodicbc.H"
#include "../drt_fluid/drt_transfer_turb_inflow.H"
#include "../drt_fluid/turbulence_statistic_manager.H"
#include "../drt_geometry/integrationcell_coordtrafo.H"
#include "../drt_geometry/position_array.H"
#include "../drt_io/io.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_utils_parallel.H"
#include "../drt_lib/drt_dofset_independent_pbc.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret_xfem.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/matlist.H"
#include "../drt_xfem/dof_distribution_switcher.H"
#include "../drt_xfem/timeInt_std_SemiLagrange.H"
#include "../drt_xfem/timeInt_std_extrapolation.H"
#include "../drt_xfem/timeInt_enr.H"
#include "../linalg/linalg_ana.H"
#include "../linalg/linalg_serialdensevector.H"


/*------------------------------------------------------------------------------------------------*
 | constructor                                                                        henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
FLD::CombustFluidImplicitTimeInt::CombustFluidImplicitTimeInt(
    const Teuchos::RCP<DRT::Discretization>&      actdis,
    const Teuchos::RCP<LINALG::Solver>&           solver,
    const Teuchos::RCP<Teuchos::ParameterList>&   params,
    const Teuchos::RCP<IO::DiscretizationWriter>& output
):TimInt(actdis, solver, params, output),
  // call constructor for "nontrivial" objects
  xparams_(params_->sublist("XFEM")),
  interfacehandle_(Teuchos::null),
  phinp_(Teuchos::null),
  cout0_(discret_->Comm(), std::cout),
  combusttype_(DRT::INPUT::IntegralValue<INPAR::COMBUST::CombustionType>(params_->sublist("COMBUSTION FLUID"),"COMBUSTTYPE")),
  veljumptype_(DRT::INPUT::IntegralValue<INPAR::COMBUST::VelocityJumpType>(params_->sublist("COMBUSTION FLUID"),"VELOCITY_JUMP_TYPE")),
  fluxjumptype_(DRT::INPUT::IntegralValue<INPAR::COMBUST::FluxJumpType>(params_->sublist("COMBUSTION FLUID"),"FLUX_JUMP_TYPE")),
  xfemtimeint_(DRT::INPUT::IntegralValue<INPAR::COMBUST::XFEMTimeIntegration>(params_->sublist("COMBUSTION FLUID"),"XFEMTIMEINT")),
  xfemtimeint_enr_(DRT::INPUT::IntegralValue<INPAR::COMBUST::XFEMTimeIntegrationEnr>(params_->sublist("COMBUSTION FLUID"),"XFEMTIMEINT_ENR")),
  xfemtimeint_enr_comp_(DRT::INPUT::IntegralValue<INPAR::COMBUST::XFEMTimeIntegrationEnrComp>(params_->sublist("COMBUSTION FLUID"),"XFEMTIMEINT_ENR_COMP")),
  flamespeed_(params_->sublist("COMBUSTION FLUID").get<double>("LAMINAR_FLAMESPEED")),
  marksteinlength_(params_->sublist("COMBUSTION FLUID").get<double>("MARKSTEIN_LENGTH")),
  moldiffusivity_(params_->sublist("COMBUSTION FLUID").get<double>("MOL_DIFFUSIVITY")),
  nitschevel_(params_->sublist("COMBUSTION FLUID").get<double>("NITSCHE_VELOCITY")),
  nitschepres_(params_->sublist("COMBUSTION FLUID").get<double>("NITSCHE_PRESSURE")),
  condensation_(xparams_.get<bool>("DLM_condensation")),
  gmshoutput_(params_->get<bool>("GMSH_OUTPUT")),
  surftensapprox_(DRT::INPUT::IntegralValue<INPAR::COMBUST::SurfaceTensionApprox>(params_->sublist("COMBUSTION FLUID"),"SURFTENSAPPROX")),
  connected_interface_(DRT::INPUT::IntegralValue<int>(params_->sublist("COMBUSTION FLUID"),"CONNECTED_INTERFACE")),
  smoothed_boundary_integration_(DRT::INPUT::IntegralValue<int>(params_->sublist("COMBUSTION FLUID"),"SMOOTHED_BOUNDARY_INTEGRATION")),
  smoothgradphi_(DRT::INPUT::IntegralValue<INPAR::COMBUST::SmoothGradPhi>(params_->sublist("COMBUSTION FLUID"),"SMOOTHGRADPHI")),
  dtele_(0.0),
  startsteps_(params_->get<int> ("number of start steps")),
  dtp_     (params_->get<double> ("time step size")),
  theta_   (params_->get<double>("theta")),
  alphaM_(params_->get<double>("alpha_M")),
  alphaF_(params_->get<double>("alpha_F")),
  gamma_(params_->get<double>("gamma")),
  initstatsol_(DRT::INPUT::IntegralValue<int>(params_->sublist("COMBUSTION FLUID"),"INITSTATSOL")),
  itemaxFRS_(params_->sublist("COMBUSTION FLUID").get<int>("ITE_MAX_FRS")),
  totalitnumFRS_(0),
  curritnumFRS_(0),
  extrapolationpredictor_(params_->get("do explicit predictor",false)),
  writestresses_(params_->get<int>("write stresses", 0)),
  project_(false),
  samstart_(-1),
  samstop_(-1)
{
  //------------------------------------------------------------------------------------------------
  // time measurement: initialization
  //------------------------------------------------------------------------------------------------
  TEUCHOS_FUNC_TIME_MONITOR(" + initialization");

  //------------------------------------------------------------------------------------------------
  // set time integration parameters for computation of initial field from stationary problem
  //------------------------------------------------------------------------------------------------
  if(initstatsol_) step_ = -1;

  //------------------------------------------------------------------------------------------------
  // set time integration parameters for stationary simulation
  //------------------------------------------------------------------------------------------------
  if (timealgo_ == INPAR::FLUID::timeint_stationary)
  {
    dta_ = 1.0;
    dtp_ = 1.0;
    theta_ = 1.0;
    cout0_ << "parameters 'theta' and 'time step size' have been set to 1.0 for stationary problem " << endl;
  }

  numdim_ = params_->get<int>("number of velocity degrees of freedom");
  if (numdim_ != 3) dserror("COMBUST only supports 3D problems.");

  //------------------------------------------------------------------------------------------------
  // connect degrees of freedom for periodic boundary conditions
  //------------------------------------------------------------------------------------------------
  {
    // set elements to 'standard mode' so the parallel redistribution (includes call of FillComplete)
    // hidden in the construction of periodic boundary condition runs correctly
    // remark: the 'elementdofmanager' would get lost packing and unpacking the elements
    ParameterList eleparams;
    eleparams.set("action","set_standard_mode");
    discret_->Evaluate(eleparams);

    // we need to keep pbc_ in order to update the PBCDofset without
    // losing the DofGid range that the PBCDofset is assigned here
    pbc_ = Teuchos::rcp(new PeriodicBoundaryConditions(discret_));
    pbc_->UpdateDofsForPeriodicBoundaryConditions();
    pbcmapmastertoslave_ = pbc_->ReturnAllCoupledColNodes();
  }

  //------------------------------------------------------------------------------------------------
  // consistency checks
  //------------------------------------------------------------------------------------------------
  if (combusttype_ == INPAR::COMBUST::combusttype_premixedcombustion and smoothgradphi_ == INPAR::COMBUST::smooth_grad_phi_none)
    dserror("Every COMBUSTTYPE except Two_Phase_Flow_... need smoothgradphi_ set to anything but smooth_grad_phi_none.");
  if (!smoothed_boundary_integration_ and smoothgradphi_ != INPAR::COMBUST::smooth_grad_phi_none)
    dserror("If SMOOTHGRADPHI is not smooth_grad_phi_none, then the SMOOTHED_BOUNDARY_INTEGRATION has to be set to Yes.");
  if (combusttype_ == INPAR::COMBUST::combusttype_twophaseflow_surf or combusttype_ == INPAR::COMBUST::combusttype_twophaseflowjump)
  {
    if ((surftensapprox_ != INPAR::COMBUST::surface_tension_approx_laplacebeltrami
        and surftensapprox_ != INPAR::COMBUST::surface_tension_approx_fixed_curvature)
        and smoothed_boundary_integration_!=true)
      dserror("All surface tension approximations need a smooth gradient field of phi expect laplace-beltrami and fixed-curvature! Read remark!");
    // surface_tension_approx_laplacebeltrami
    // surface_tension_approx_fixed_curvature: we can use a smoothed normal vector based on the smoothed gradient of phi (SMOOTHGRADPHI = Yes)
    //                                         or we simply use the normal vector of the boundary integration cell (SMOOTHGRADPHI = No)
    // surface_tension_approx_divgrad
    // surface_tension_approx_divgrad_normal: we use the smoothed gradient of phi to compute the curvature of the interface, hence we
    //                                        always have to set SMOOTHGRADPHI = Yes, in addition, we need the normal to the interface,
    //                                        version surface_tension_approx_divgrad_normal simply uses the normal vector of the boundary
    //                                        integration cell, while surface_tension_approx_divgrad uses a smoothed one based on phi
    // surface_tension_approx_laplacebeltrami_smoothed: here, we need a smoothed and a non-smoothed normal, hence, SMOOTHGRADPHI = Yes
    if (params_->sublist("COMBUSTION FLUID").get<double>("VARIABLESURFTENS")!=0.0)
    {
      if (myrank_ == 0)
      {
        std::cout << "-> location-dependent surface tension coefficient:" << std::endl;
        std::cout << "only coefficient linear in x available for the time being" << std::endl;
        if (surftensapprox_ == INPAR::COMBUST::surface_tension_approx_laplacebeltrami or
            surftensapprox_ == INPAR::COMBUST::surface_tension_approx_laplacebeltrami_smoothed)
           std::cout << "currently only connected interface possible for laplace-beltrami" << std::endl;
      }
    }
  }
  if (combusttype_ == INPAR::COMBUST::combusttype_twophaseflowjump
      and (veljumptype_ != INPAR::COMBUST::vel_jump_none or fluxjumptype_ != INPAR::COMBUST::flux_jump_surface_tension))
  {
    if (veljumptype_ != INPAR::COMBUST::vel_jump_none)
    {
      veljumptype_ = INPAR::COMBUST::vel_jump_none;
      std::cout << "Velocity jump is set to NONE for two-phase flow with jumps!" << std::endl;
    }
    if (fluxjumptype_ != INPAR::COMBUST::flux_jump_surface_tension)
    {
      fluxjumptype_ = INPAR::COMBUST::flux_jump_surface_tension;
      std::cout << "Flux jump is set to SURFACE TENSION for two-phase flow with jumps!" << std::endl;
    }
  }

  //------------------------------------------------------------------------------------------------
  // Neumann inflow term if required
  //------------------------------------------------------------------------------------------------
  neumanninflow_ = false;
  if (params_->get<string>("Neumann inflow","no") == "yes") neumanninflow_ = true;

  //------------------------------------------------------------------------------------------------
  // prepare XFEM (initial degree of freedom management)
  //------------------------------------------------------------------------------------------------
  // special option for two-phase flows
  // select enrichments for one or both fields or suppress enrichments completely, i.e., do usual FEM with intersected elements
  xparams_.set<int>("SELECTED_ENRICHMENT",DRT::INPUT::IntegralValue<INPAR::COMBUST::SelectedEnrichment>(params_->sublist("COMBUSTION FLUID"),"SELECTED_ENRICHMENT"));
  
  physprob_.xfemfieldset_.clear();
  // declare physical fields of the problem (continuous fields and discontinuous XFEM fields)
  physprob_.xfemfieldset_.insert(XFEM::PHYSICS::Velx);
  physprob_.xfemfieldset_.insert(XFEM::PHYSICS::Vely);
  physprob_.xfemfieldset_.insert(XFEM::PHYSICS::Velz);
#ifdef COMBUST_NORMAL_ENRICHMENT
  physprob_.xfemfieldset_.insert(XFEM::PHYSICS::Veln);
#endif
  physprob_.xfemfieldset_.insert(XFEM::PHYSICS::Pres);
#ifdef COMBUST_STRESS_BASED
#ifdef COMBUST_EPSPRES_BASED
  // define approach for extra stress field (stress-based Lagrange Multiplier approach)
  physprob_.elementAnsatz_ = rcp(new COMBUST::EpsilonPressureAnsatz());
#endif
#ifdef COMBUST_SIGMA_BASED
  // define approach for extra stress field (stress-based Lagrange Multiplier approach)
  physprob_.elementAnsatz_ = rcp(new COMBUST::CauchyStressAnsatz());
#endif
#else
  // for the Nitsche method assign an arbitrary element ansatz to compile
  physprob_.elementAnsatz_ = rcp(new COMBUST::TauPressureAnsatz());
#endif
  // create dummy instance of interfacehandle holding no flamefront and hence no integration cells
  Teuchos::RCP<COMBUST::InterfaceHandleCombust> ihdummy = rcp(new COMBUST::InterfaceHandleCombust(discret_,Teuchos::null));
  // create dummy instance of dof manager assigning standard enrichments to all nodes
  const Teuchos::RCP<XFEM::DofManager> dofmanagerdummy = rcp(new XFEM::DofManager(ihdummy,Teuchos::null,physprob_.xfemfieldset_,*physprob_.elementAnsatz_,xparams_,Teuchos::null));

  // save dofmanager to be able to plot Gmsh stuff in Output()
  dofmanagerForOutput_ = dofmanagerdummy;

  // pass dof information to elements (no enrichments yet, standard FEM!)
  TransferDofInformationToElements(Teuchos::null, dofmanagerdummy);
  // ensure that degrees of freedom in the discretization have been set
  discret_->FillComplete();

  output_->WriteMesh(0,0.0);

  // store a dofset with the complete fluid unknowns
  standarddofset_ = Teuchos::rcp(new DRT::IndependentPBCDofSet(pbcmapmastertoslave_));
  standarddofset_->Reset();
  standarddofset_->AssignDegreesOfFreedom(*discret_,0,0);

  velpressplitterForOutput_ = rcp(new LINALG::MapExtractor());
  // split based on complete fluid field
  FLD::UTILS::SetupFluidSplit(*discret_,*standarddofset_,3,*velpressplitterForOutput_);

  //------------------------------------------------------------------------------------------------
  // get dof layout from the discretization to construct vectors and matrices
  //------------------------------------------------------------------------------------------------

  // parallel dof distribution contained in dofrowmap: local (LID) <-> global (GID) dof numbering
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // get layout of velocity and pressure dofs in a vector
  const int numdim = params_->get<int>("number of velocity degrees of freedom");
  velpressplitter_ = rcp(new LINALG::MapExtractor());
  FLD::UTILS::SetupFluidSplit(*discret_,numdim,*velpressplitter_);

  //------------------------------------------------------------------------------------------------
  // create empty system matrix - stiffness and mass are assembled in one system matrix!
  //------------------------------------------------------------------------------------------------

  // initialize standard (stabilized) system matrix (and save its graph!)
  sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,27,false,true));

  //------------------------------------------------------------------------------------------------
  // create empty vectors - used for different purposes
  //------------------------------------------------------------------------------------------------

  //-------------------------------
  // vectors passed to the element
  //-------------------------------
  // velocity/pressure at time step n+1, n and n-1
  state_.velnp_ = LINALG::CreateVector(*dofrowmap,true);
  state_.veln_  = LINALG::CreateVector(*dofrowmap,true);
  state_.velnm_ = LINALG::CreateVector(*dofrowmap,true);

  // acceleration at time n+1 and n
  state_.accnp_  = LINALG::CreateVector(*dofrowmap,true);
  state_.accn_   = LINALG::CreateVector(*dofrowmap,true);

  if (timealgo_ == INPAR::FLUID::timeint_afgenalpha)
  {
    state_.velaf_ = LINALG::CreateVector(*dofrowmap,true);
    state_.accam_  = LINALG::CreateVector(*dofrowmap,true);
  }
  else
  {
    state_.velaf_ = Teuchos::null;
    state_.accam_  = Teuchos::null;
  }

  state_.nodalDofDistributionMap_.clear();

  dofmanagerForOutput_->fillDofRowDistributionMaps(
      state_.nodalDofDistributionMap_);

  //---------------------------------------------
  // vectors associated with boundary conditions
  //---------------------------------------------

  // a vector of zeros to be used to enforce zero dirichlet boundary conditions
  zeros_ = LINALG::CreateVector(*dofrowmap,true);

  //-----------------------------------
  // vectors used for solution process
  //-----------------------------------

  // rhs: standard (stabilized) residual vector (rhs for the incremental form)
  residual_     = LINALG::CreateVector(*dofrowmap,true);
  trueresidual_ = LINALG::CreateVector(*dofrowmap,true);

  // nonlinear iteration increment vector
  incvel_ = LINALG::CreateVector(*dofrowmap,true);

  //------------------------------------------------------------------------------------------------
  // prepare turbulence specific stuff
  //------------------------------------------------------------------------------------------------

  // ----------------------------------------------------
  // initialize vectors and flags for turbulence approach
  // ----------------------------------------------------

  // get list of turbulence parameters
  ParameterList *  modelparams =&(params_->sublist("TURBULENCE MODEL"));

  // read turbulence model
  string physmodel = modelparams->get<string>("PHYSICAL_MODEL","no_model");
  if (physmodel != "no_model")
    dserror("XFEM does not support any fancy turbulence models");

  // read flag for special flow
  special_flow_ = modelparams->get<string>("CANONICAL_FLOW","no");
  if ( special_flow_ != "no" and
       special_flow_ != "bubbly_channel_flow" and
       special_flow_ != "combust_oracles" and
       special_flow_ != "backward_facing_step_tp")
    dserror("XFEM does not support special flow %s.", special_flow_.c_str());

  // -------------------------------------------------------------------
  // initialize turbulence statistics evaluation
  // -------------------------------------------------------------------
  turbstatisticsmanager_ = rcp(new FLD::TurbulenceStatisticManager(*this));

  // parameter for sampling/dumping period
  if (special_flow_ != "no")
  {
    samstart_ = modelparams->get<int>("SAMPLING_START",1);
    samstop_  = modelparams->get<int>("SAMPLING_STOP",1);
  }
}

/*------------------------------------------------------------------------------------------------*
 | destructor                                                                         henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
FLD::CombustFluidImplicitTimeInt::~CombustFluidImplicitTimeInt()
{
  return;
}

/*------------------------------------------------------------------------------------------------*
 | out of order!                                                                      henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::Integrate()
{
  dserror("Thou shalt not use this function! Switch over transient/stationary scheme is in combust_dyn!");
}

/*------------------------------------------------------------------------------------------------*
 | out of order!                                                                      henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::TimeLoop()
{
  dserror("Thou shalt not use this function! Use COMBUST::Algorithm::TimeLoop() instead");
}

/*------------------------------------------------------------------------------------------------*
 | check whether FRS iteration is finished                                       winklmaier 11/11 |
 *------------------------------------------------------------------------------------------------*/
bool FLD::CombustFluidImplicitTimeInt::FluidRefSolLoopFinished()
{
  bool finished = false;

  if (timealgo_ == INPAR::FLUID::timeint_stationary
      and curritnumFRS_>=1)
    finished = true; // no FRS iterations in stationary case

  if (curritnumFRS_>=itemaxFRS_) finished = true;

  if (curritnumFRS_>=1) // check whether the algorithms require additional iterations
  {
    if (step_<1) // stationary solution just one time
      finished = true;

    switch (xfemtimeint_enr_) // currently no enrichment computation requires >1 iterations
    {
    case INPAR::COMBUST::xfemtimeintenr_donothing:
    case INPAR::COMBUST::xfemtimeintenr_quasistatic:
    case INPAR::COMBUST::xfemtimeintenr_project:
    case INPAR::COMBUST::xfemtimeintenr_project_scalar:
      break;
    default:
      dserror("enrichment recomputation approach in XFEM time integration not implemented");
    }

    switch (xfemtimeint_)
    {
    case INPAR::COMBUST::xfemtimeint_donothing:
    case INPAR::COMBUST::xfemtimeint_extrapolation:
    {
      finished = true; // the above standard value computations don't require -> iterations
      break;
    }
    case INPAR::COMBUST::xfemtimeint_semilagrange:
    case INPAR::COMBUST::xfemtimeint_mixedSLExtrapol:
    case INPAR::COMBUST::xfemtimeint_mixedSLExtrapolNew:
      break;
    default:
      dserror("standard recomputation approach in XFEM time integration not implemented");
    }
  }

  // if not finished, print iteration counter, update false
  if (!finished)
  {
    totalitnumFRS_++;
    curritnumFRS_++;
    if (myrank_==0 and itemaxFRS_!=1)
      printf("FLUID REFERENCE SOLUTION ITERATION %3d/%3d\n",curritnumFRS_,itemaxFRS_);
  }
  else
    curritnumFRS_=0; // for new fgi or time step
  return finished;
}



void FLD::CombustFluidImplicitTimeInt::ClearTimeInt()
{
  timeIntStd_ = Teuchos::null;
  timeIntEnr_ = Teuchos::null;
}



/*------------------------------------------------------------------------------------------------*
 | prepare a fluid time step                                                          henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::PrepareTimeStep()
{
  totalitnumFRS_ = 0;

  // update old acceleration
  if (state_.accn_ != Teuchos::null)
    state_.accn_->Update(1.0,*state_.accnp_,0.0);

  // velocities/pressures of this step become most recent
  // velocities/pressures of the last step
  if (state_.velnm_ != Teuchos::null)
    state_.velnm_->Update(1.0,*state_.veln_ ,0.0);
  if (state_.veln_ != Teuchos::null)
    state_.veln_ ->Update(1.0,*state_.velnp_,0.0);

  //---------------------------- -------------------------------------------------------------------
  // set time dependent parameters
  // -----------------------------------------------------------------------------------------------
  // reset time step size, 'theta' and time integration algorithm
  dta_      = params_->get<double> ("time step size");
  dtp_      = params_->get<double> ("time step size");
  theta_    = params_->get<double>("theta");
  timealgo_ = DRT::INPUT::get<INPAR::FLUID::TimeIntegrationScheme>(*params_, "time int algo");
  itemax_   = params_->get<int>("max nonlin iter steps");

  step_ += 1;
  time_ += dta_;

  switch(timealgo_)
  {
  case INPAR::FLUID::timeint_stationary:
  {
    // for stationary problems PrepareTimeStep() should only be called for output related reasons
    cout0_ << "/!\\ warning: 'time' and 'time step' are set to 1.0 and 1.0 for output control file" << endl;
    step_ = 1;
    time_ = 1.0;
    theta_ = 1.0;
    break;
  }
  case INPAR::FLUID::timeint_one_step_theta:
  {
    // compute initial field from stationary problem
    if ((step_==0) && initstatsol_)
    {
      cout0_ << "/!\\ warning: initial solution is computed by stationary algorithm" << endl;
      cout0_ << "/!\\ warning: 'time' and 'time step' are set to 0.0 and 1.0 for output control file" << endl;
      timealgo_ = INPAR::FLUID::timeint_stationary;
      time_ =  0.0; // only needed for output
      dta_ =   1.0; // for calculation, we reset this value at the end of NonlinearSolve()
      dtp_ =   1.0; // for calculation, we reset this value at the end of NonlinearSolve()
      theta_ = 1.0;
      // set max iterations for initial stationary algorithm
      itemax_ = params_->get<int>("max nonlin iter steps init stat sol");
    }
    // compute first (instationary) time step differently
    // remark: usually backward Euler (theta=1.0) to damp inital pertubations
    else if (step_==1)
    {
      // get starting 'theta' for first time step
      theta_ = params_->get<double>("start theta");
      cout0_ << "/!\\ first time step computed with theta =  " << theta_ << endl;
    }
    // regular time step
    else if (step_ > 1)
    {
      theta_ = params_->get<double>("theta");
    }
    else
      dserror("number of time step is wrong");
    break;
  }
  case INPAR::FLUID::timeint_afgenalpha:
  {
    // use start parameters for generalized alpha scheme
    if (startsteps_>0)
    {
      if (step_<=startsteps_)
      {
        if (startsteps_ > stepmax_)
          dserror("more starting steps than total time steps");

        cout0_ << "/!\\ first " << startsteps_ << " steps are computed with Backward-Euler scheme" << endl;
        // use backward-Euler-type parameter combination
        alphaM_ = 1.0;
        alphaF_ = 1.0;
        gamma_  = 1.0;
      }
      else
      {
        // recall original parameters from input file
        alphaM_ = params_->get<double>("alpha_M");
        alphaF_ = params_->get<double>("alpha_F");
        gamma_  = params_->get<double>("gamma");
      }
    }
    // compute "pseudo-theta" for af-generalized-alpha scheme
    theta_ = alphaF_*gamma_/alphaM_;

    break;
  }
  case INPAR::FLUID::timeint_bdf2:
  {
    // do a backward Euler step for the first time step
    if (step_==1)
    {
      timealgo_ = INPAR::FLUID::timeint_one_step_theta;
      theta_ = params_->get<double>("start theta");
    }
    // regular time step (step_>1)
    else
    {
      // for BDF2, theta is set by the time-step sizes, 2/3 for const dt
      theta_ = (dta_+dtp_)/(2.0*dta_ + dtp_);
    }
    break;
  }
  default:
    dserror("unknown time integration scheme");
  }
}

/*------------------------------------------------------------------------------------------------*
 | prepare a fluid nonlinear iteration                                                henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::PrepareNonlinearSolve()
{

  // -------------------------------------------------------------------
  // set part of the rhs vector belonging to the old timestep
  // -------------------------------------------------------------------
  // TODO Do we need this?
//  TIMEINT_THETA_BDF2::SetOldPartOfRighthandside(
//      state_.veln_, state_.velnm_, state_.accn_,
//          timealgo_, dta_, theta_,
//          hist_);

  // -------------------------------------------------------------------
  //                     do explicit predictor step
  // -------------------------------------------------------------------
  //
  // We cannot have a predictor in case of monolithic FSI here. There needs to
  // be a way to turn this off.
//  if (extrapolationpredictor_)
//  {
//    // TODO Was davon brauchen wir?
//    if (step_>1)
//    {
//      double timealgo_constant=theta_;
//
//      TIMEINT_THETA_BDF2::ExplicitPredictor(
//        "default",
//        state_.veln_,
//        state_.velnm_,
//        state_.accn_,
//        velpressplitter_,
//        timealgo_,
//        timealgo_constant,
//        dta_,
//        dtp_,
//        state_.velnp_,
//        discret_->Comm());
//    }
//  }

  // -------------------------------------------------------------------
  //         evaluate dirichlet and neumann boundary conditions
  // -------------------------------------------------------------------
  {
    ParameterList eleparams;

    // other parameters needed by the elements
    eleparams.set("total time",time_);
    eleparams.set("delta time",dta_);
    eleparams.set("thsl",theta_*dta_);


    if (step_ > 0) // initial stationary solution does not require reference solution update
    {
      TEUCHOS_FUNC_TIME_MONITOR("   + xfem time integration");

      vector<RCP<Epetra_Vector> > newRowVectorsn; // solution vectors of old time step due to new interface
      newRowVectorsn.push_back(state_.veln_);
      newRowVectorsn.push_back(state_.accn_);

      vector<RCP<Epetra_Vector> > newRowVectorsnp; // solution vectors at new time step
      newRowVectorsnp.push_back(state_.velnp_);
      newRowVectorsnp.push_back(state_.accnp_);

      if(xfemtimeint_ != INPAR::COMBUST::xfemtimeint_donothing) // do something for standard values
      {
        cout0_ << "---  XFEM time integration: adopt standard values to new interface... " << std::flush;
        timeIntStd_->type(totalitnumFRS_,itemaxFRS_); // update algorithm handling
        timeIntStd_->compute(newRowVectorsn,newRowVectorsnp); // call computation
        cout0_ << "done" << std::endl;
      }

      if ((xfemtimeint_enr_!=INPAR::COMBUST::xfemtimeintenr_donothing) and
          (xfemtimeint_enr_!=INPAR::COMBUST::xfemtimeintenr_quasistatic)) // do something for enrichment values
      {
        cout0_ << "---  XFEM time integration: adopt enrichment values to new interface... " << std::flush;
        timeIntEnr_->type(totalitnumFRS_,itemaxFRS_); // update algorithm handling
        timeIntEnr_->compute(newRowVectorsn,newRowVectorsnp); // call computation
        cout0_ << "done" << std::endl;
      }
    }

    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("velnp",state_.velnp_);
    // predicted dirichlet values
    // velnp then also holds prescribed new dirichlet values
    RCP<DRT::DiscretizationXFEM> xdiscret = Teuchos::rcp_dynamic_cast<DRT::DiscretizationXFEM>(discret_, true);
    xdiscret->EvaluateDirichletCombust(eleparams,state_.velnp_,null,null,null,dbcmaps_);

    discret_->ClearState();

    // Transfer of boundary data if necessary
    turbulent_inflow_condition_->Transfer(state_.veln_,state_.velnp_,time_);

    // evaluate Neumann conditions
    neumann_loads_->PutScalar(0.0);
    discret_->EvaluateNeumann(eleparams,*neumann_loads_);
    discret_->ClearState();
  }

  if (timealgo_==INPAR::FLUID::timeint_afgenalpha)
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

  // sysmat might be singular (if we have a purely Dirichlet constrained
  // problem, the pressure mode is defined only up to a constant)
  // in this case, we need a basis vector for the nullspace/kernel
  {
    vector<DRT::Condition*> KSPcond;
    discret_->GetCondition("KrylovSpaceProjection",KSPcond);
    int numcond = KSPcond.size();

    int numfluid = 0;

    for(int icond = 0; icond < numcond; icond++)
    {
      const std::string* name = KSPcond[icond]->Get<std::string>("discretization");
      if (*name == "fluid")
        numfluid++;
    }

    if (numfluid == 1)
    {
      project_ = true;
      w_       = LINALG::CreateVector(*discret_->DofRowMap(),true);
      c_       = LINALG::CreateVector(*discret_->DofRowMap(),true);

      kspsplitter_ = rcp(new FLD::UTILS::KSPMapExtractor());
      kspsplitter_->Setup(*discret_);
    }
    else if (numfluid == 0)
    {
      project_ = false;
      w_       = Teuchos::null;
      c_       = Teuchos::null;
    }
    else
      dserror("Received more than one KrylovSpaceCondition for fluid field");
  }

}

/*------------------------------------------------------------------------------------------------*
 | hand over information about (XFEM) degrees of freedom to elements                     ag 04/09 |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::TransferDofInformationToElements(
    const Teuchos::RCP<COMBUST::FlameFront> flamefront,
    const Teuchos::RCP<XFEM::DofManager> dofmanager
    )
{
  ParameterList eleparams;
  eleparams.set("action","store_xfem_info");
  eleparams.set("dofmanager",dofmanager);
  eleparams.set("DLM_condensation",xparams_.get<bool>("DLM_condensation"));
  eleparams.set("flamefront",flamefront);
//  eleparams.set("phinp",phinp);
//  eleparams.set("gradphi",gradphi);
//  eleparams.set("curvature",curvature);

  discret_->Evaluate(eleparams);
}

/*------------------------------------------------------------------------------------------------*
 | import geometrical information about the interface (integration cells) from the combustion     |
 | algorithm and incorporate it into the fluid field                                  henke 03/09 |
 |
 | remark: Within this routine, no parallel re-distribution is allowed to take place. Before and  |
 | after this function, it's ok to do so.                                                    axel |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::IncorporateInterface(const Teuchos::RCP<COMBUST::FlameFront>& flamefront)
{
  // build instance of DofManager with information about the interface from the interfacehandle
  // remark: DofManager is rebuilt in every inter-field iteration step, because number and position
  // of enriched degrees of freedom change in every time step/FG-iteration
  const Teuchos::RCP<XFEM::DofManager> dofmanager = rcp(new XFEM::DofManager(
      interfacehandle_,
      phinp_,
      physprob_.xfemfieldset_,
      *physprob_.elementAnsatz_,
      xparams_,
      pbcmapmastertoslave_)
  );

  // temporarely save old dofmanager
  const RCP<XFEM::DofManager> olddofmanager = dofmanagerForOutput_;

  // save dofmanager to be able to plot Gmsh stuff in Output()
  dofmanagerForOutput_ = dofmanager;

  // print global and element dofmanager to Gmsh
  dofmanager->toGmsh(step_);
  interfacehandle_->toGmsh(step_);

  // get old dofmaps, compute a new one and get the new one, too
  const Epetra_Map olddofrowmap = *discret_->DofRowMap();
  const Epetra_Map olddofcolmap = *discret_->DofColMap();
  map<XFEM::DofKey,XFEM::DofGID> oldNodalDofColDistrib;
  olddofmanager->fillNodalDofColDistributionMap(oldNodalDofColDistrib);

  // -------------------------------------------------------------------
  // connect degrees of freedom for periodic boundary conditions
  // -------------------------------------------------------------------
  // tell elements about the dofs and the integration
  TransferDofInformationToElements(flamefront, dofmanager);
  // assign degrees of freedom
  // remark: - assign degrees of freedom (first slot)
  //         - build geometry for (Neumann) boundary conditions (third slot);
  //           without Neumann boundary conditions Fillcomplete(true,false,false) will also work
  discret_->FillComplete(true,false,true);
  const Epetra_Map& newdofrowmap = *discret_->DofRowMap();

  // remark: 'true' is needed to prevent iterative solver from crashing
  discret_->ComputeNullSpaceIfNecessary(solver_->Params(),true);

  // anonymous namespace for dofswitcher and startvalues
  {
    //const std::map<XFEM::DofKey, XFEM::DofGID> oldNodalDofDistributionMap(state_.nodalDofDistributionMap_);
    std::map<XFEM::DofKey, XFEM::DofGID> oldNodalDofDistributionMap(state_.nodalDofDistributionMap_);
    dofmanager->fillDofRowDistributionMaps(
        state_.nodalDofDistributionMap_);

    // create switcher
    const XFEM::DofDistributionSwitcher dofswitch(
        interfacehandle_, dofmanager,
        olddofrowmap, newdofrowmap,
        oldNodalDofDistributionMap, state_.nodalDofDistributionMap_
    );

    //---------------------------------------------------------------
    // extract old enrichment dofkeys and values before they are lost
    //---------------------------------------------------------------
    vector<RCP<Epetra_Vector> > oldColStateVectors; // same order of vectors as newRowVectorsn combined with newRowVectorsnp
    RCP<Epetra_Vector> veln = rcp(new Epetra_Vector(olddofcolmap,true));
    {
      LINALG::Export(*state_.veln_,*veln);
      oldColStateVectors.push_back(veln);

      RCP<Epetra_Vector> accn = rcp(new Epetra_Vector(olddofcolmap,true));
      LINALG::Export(*state_.accn_,*accn);
      oldColStateVectors.push_back(accn);

      RCP<Epetra_Vector> velnp = rcp(new Epetra_Vector(olddofcolmap,true));
      LINALG::Export(*state_.velnp_,*velnp);
      oldColStateVectors.push_back(velnp);

      RCP<Epetra_Vector> accnp = rcp(new Epetra_Vector(olddofcolmap,true));
      LINALG::Export(*state_.accnp_,*accnp);
      oldColStateVectors.push_back(accnp);
    }

    //---------------------------------------------
    // switch state vectors to new dof distribution
    //---------------------------------------------
    cout0_ << "---  transform state vectors... " << std::flush;
    // quasi-static enrichment strategy for kink enrichments
    // remark: as soon as the XFEM-time-integration works for kinks, this should be removed
    if (xfemtimeint_enr_==INPAR::COMBUST::xfemtimeintenr_quasistatic)
    {
      cout0_ << "quasi-static enrichment for two-phase flow problems... " << std::flush;

      // accelerations at time n+1 and n
      dofswitch.mapVectorToNewDofDistributionCombust(state_.accnp_,true);
      dofswitch.mapVectorToNewDofDistributionCombust(state_.accn_, true);
      if (timealgo_ == INPAR::FLUID::timeint_afgenalpha)
      {
        dofswitch.mapVectorToNewDofDistributionCombust(state_.accam_,true);
        dofswitch.mapVectorToNewDofDistributionCombust(state_.velaf_,true);
      }
      else
        dofswitch.mapVectorToNewDofDistributionCombust(state_.velnm_,true); //TODO: is this needed for genalpha?
      // velocities and pressures at time n+1, n and n-1
      dofswitch.mapVectorToNewDofDistributionCombust(state_.velnp_,true); // use old velocity as start value
      dofswitch.mapVectorToNewDofDistributionCombust(state_.veln_ ,true);
    }
    else
    {
      // accelerations at time n+1 and n
      dofswitch.mapVectorToNewDofDistributionCombust(state_.accnp_,false);
      dofswitch.mapVectorToNewDofDistributionCombust(state_.accn_ ,false);
      if (timealgo_ == INPAR::FLUID::timeint_afgenalpha)
      {
        dofswitch.mapVectorToNewDofDistributionCombust(state_.accam_,false);
        dofswitch.mapVectorToNewDofDistributionCombust(state_.velaf_,false);
      }
      else
        dofswitch.mapVectorToNewDofDistributionCombust(state_.velnm_,false);
      // velocities and pressures at time n+1, n and n-1
      dofswitch.mapVectorToNewDofDistributionCombust(state_.velnp_,false); // use old velocity as start value
      dofswitch.mapVectorToNewDofDistributionCombust(state_.veln_ ,false);
    }
    cout0_ << "done" << endl;

    // all initial values can be set now; including the enrichment values
    if (step_ == 0)
      SetEnrichmentField(dofmanager,newdofrowmap); // set initial (enrichment) field

    if (timealgo_ == INPAR::FLUID::timeint_afgenalpha
        and (xfemtimeint_!=INPAR::COMBUST::xfemtimeint_donothing
        or ((xfemtimeint_enr_ !=INPAR::COMBUST::xfemtimeintenr_donothing) and (xfemtimeint_enr_ !=INPAR::COMBUST::xfemtimeintenr_quasistatic))))
      dserror("Genalpha has not been tested with anything but Timeint DoNothing and TimeintEnr DoNothing/QuasiStatic.");

    if (step_ > 0) // reference solution update not in first update
    {
      if (totalitnumFRS_==0) // construct time int classes once every time step
      {
        if ((xfemtimeint_!=INPAR::COMBUST::xfemtimeint_donothing) or
            ((xfemtimeint_enr_ !=INPAR::COMBUST::xfemtimeintenr_donothing) and (xfemtimeint_enr_ !=INPAR::COMBUST::xfemtimeintenr_quasistatic)))
        {
          // basic time integration data
          RCP<XFEM::TIMEINT> timeIntData = rcp(new XFEM::TIMEINT(
              discret_,
              olddofmanager,
              dofmanager,
              oldColStateVectors,
              flamefront,
              olddofcolmap,
              newdofrowmap,
              oldNodalDofColDistrib,
              state_.nodalDofDistributionMap_,
              pbcmapmastertoslave_));

          switch (xfemtimeint_enr_)
          {
          case INPAR::COMBUST::xfemtimeintenr_donothing:
          case INPAR::COMBUST::xfemtimeintenr_quasistatic:
            break; // nothing to do
          case INPAR::COMBUST::xfemtimeintenr_project:
          case INPAR::COMBUST::xfemtimeintenr_project_scalar:
          {
            INPAR::COMBUST::XFEMTimeIntegrationEnrComp timeIntEnrComp(DRT::INPUT::IntegralValue<INPAR::COMBUST::XFEMTimeIntegrationEnrComp>(params_->sublist("COMBUSTION FLUID"),"XFEMTIMEINT_ENR_COMP"));

            // enrichment time integration data
            timeIntEnr_ = rcp(new XFEM::EnrichmentProjection(
                *timeIntData,
                xfemtimeint_enr_,
                timeIntEnrComp));
            break;
          }
          default:
            dserror("enrichment recomputation approach in XFEM time integration not implemented");
          }

          switch (xfemtimeint_)
          {
          case INPAR::COMBUST::xfemtimeint_donothing:
            break;
          case INPAR::COMBUST::xfemtimeint_semilagrange:
          case INPAR::COMBUST::xfemtimeint_mixedSLExtrapol:
          case INPAR::COMBUST::xfemtimeint_mixedSLExtrapolNew:
          {
            // time integration data for standard dofs, semi-lagrangian approach
            timeIntStd_ = rcp(new XFEM::SemiLagrange(
                *timeIntData,
                xfemtimeint_,
                veln,
                dta_,
                theta_,
                flamefront,
                flamespeed_,
                true));
            break;
          }
          case INPAR::COMBUST::xfemtimeint_extrapolation:
          {
            // time integration data for standard dofs, extrapolation approach
            timeIntStd_ = rcp(new XFEM::Extrapolation(
                *timeIntData,
                xfemtimeint_,
                veln,
                dta_,
                flamefront,
                true));
            break;
          }
          default:
            dserror("standard recomputation approach in XFEM time integration not implemented");
          }
        }
      }

      else if (curritnumFRS_==0) // new FGI requires little new data to check which dofs shall be recomputed
      {
        switch (xfemtimeint_enr_)
        {
        case INPAR::COMBUST::xfemtimeintenr_donothing:
        case INPAR::COMBUST::xfemtimeintenr_quasistatic: break; // nothing to do
        case INPAR::COMBUST::xfemtimeintenr_project:
        case INPAR::COMBUST::xfemtimeintenr_project_scalar:
        {
          timeIntEnr_->importNewFGIData(
              discret_,
              dofmanager,
              flamefront,
              newdofrowmap,
              state_.nodalDofDistributionMap_,
              oldNodalDofColDistrib); // new data due to gamma^n+1,i+1
          break;
        }
        default: dserror("enrichment recomputation approach in XFEM time integration not implemented");
        }

        switch (xfemtimeint_)
        {
        case INPAR::COMBUST::xfemtimeint_donothing: break;
        case INPAR::COMBUST::xfemtimeint_semilagrange:
        case INPAR::COMBUST::xfemtimeint_mixedSLExtrapol:
        case INPAR::COMBUST::xfemtimeint_mixedSLExtrapolNew:
        case INPAR::COMBUST::xfemtimeint_extrapolation:
        {
          timeIntStd_->importNewFGIData(
              discret_,
              dofmanager,
              flamefront,
              newdofrowmap,
              state_.nodalDofDistributionMap_); // new data due to gamma^n+1,i+1
          break;
        }
        default: dserror("standard recomputation approach in XFEM time integration not implemented");
        }
      } // end else if, nothing more to do
    }
  } // anonymous namespace for dofswitcher and startvalues

  // --------------------------------------------------
  // create remaining vectors with new dof distribution
  // --------------------------------------------------

  zeros_        = LINALG::CreateVector(newdofrowmap,true);

  // object holds maps/subsets for DOFs subjected to Dirichlet BCs and otherwise
  dbcmaps_ = Teuchos::rcp(new LINALG::MapExtractor());
  {
    ParameterList eleparams;
    // other parameters needed by the elements
    eleparams.set("total time",time_);

    RCP<DRT::DiscretizationXFEM> xdiscret = Teuchos::rcp_dynamic_cast<DRT::DiscretizationXFEM>(discret_, true);
    xdiscret->EvaluateDirichletCombust(eleparams, zeros_, Teuchos::null, Teuchos::null, Teuchos::null, dbcmaps_);
    zeros_->PutScalar(0.0); // just in case of change
  }

  neumann_loads_= LINALG::CreateVector(newdofrowmap,true);

  // ---------------------------------
  // Vectors used for solution process
  // ---------------------------------
  residual_     = LINALG::CreateVector(newdofrowmap,true);
  trueresidual_ = LINALG::CreateVector(newdofrowmap,true);
  incvel_       = LINALG::CreateVector(newdofrowmap,true);


  // -------------------------------------------------------------------
  // get a vector layout from the discretization for a vector which only
  // contains the velocity dofs and for one vector which only contains
  // pressure degrees of freedom.
  // -------------------------------------------------------------------
  SetupXFluidSplit(*discret_,dofmanager,*velpressplitter_);

  // -------------------------------------------------------------------------------------
  // create empty system matrix --- stiffness and mass are assembled in one system matrix!
  // -------------------------------------------------------------------------------------
  // initialize system matrix
  sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(newdofrowmap,0,false,true));

  // -------------------------------------------------------------------
  // check whether we have a coupling to a turbulent inflow generating
  // computation and initialise the transfer if necessary
  // -------------------------------------------------------------------
  turbulent_inflow_condition_ = Teuchos::rcp(new TransferTurbulentInflowCondition(discret_,dbcmaps_));
}


/*------------------------------------------------------------------------------------------------*
 | import geometrical information about the interface (integration cells) from the combustion     |
 | algorithm and incorporate it into the fluid field                                  henke 03/09 |
 |
 | remark: Within this routine, no parallel re-distribution is allowed to take place. Before and  |
 | after this function, it's ok to do so.                                                    axel |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::StoreFlameFront(
    const Teuchos::RCP<COMBUST::FlameFront> flamefront,
    bool UpdateDofSet)
{
  if (flamefront!=Teuchos::null)
  {
    interfacehandle_ = flamefront->InterfaceHandle();
    phinp_ = flamefront->Phinp();

    if (UpdateDofSet)
      IncorporateInterface(flamefront);
  }
  else
  {
    interfacehandle_ = Teuchos::null;
    phinp_ = Teuchos::null;
  }


  return;
}


/*------------------------------------------------------------------------------------------------*
 | get convection velocity vector for transfer to scalar transport field           wichmann 02/12 |
 *------------------------------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FLD::CombustFluidImplicitTimeInt::StdVeln(void)
{
  // velocity vector has to be transformed from XFEM format to Standard FEM format, because ScaTra
  // cannot handle XFEM dofs at enriched nodes
  // remark: has to include pressure, since ScaTraTimIntImpl::SetVelocityField() expects it!
  std::set<XFEM::PHYSICS::Field> outputfields;
  outputfields.insert(XFEM::PHYSICS::Velx);
  outputfields.insert(XFEM::PHYSICS::Vely);
  outputfields.insert(XFEM::PHYSICS::Velz);
  outputfields.insert(XFEM::PHYSICS::Pres);

  // TODO: check performance time to built convel vector; if this is costly, it could be stored as
  //       a private member variable of the time integration scheme, since it is built in two places
  //       (here in every nonlinear iteration, and in the Output() function after every time step)

  // convection velocity vector
  const Teuchos::RCP<Epetra_Vector> convel = dofmanagerForOutput_->transformXFEMtoStandardVector(*state_.veln_,
                                                               *standarddofset_,
                                                               state_.nodalDofDistributionMap_,
                                                               outputfields);
  return convel;
}


/*------------------------------------------------------------------------------------------------*
 | get convection velocity vector for transfer to scalar transport field           wichmann 02/12 |
 *------------------------------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FLD::CombustFluidImplicitTimeInt::StdVelnp(void)
{
  // velocity vector has to be transformed from XFEM format to Standard FEM format, because ScaTra
  // cannot handle XFEM dofs at enriched nodes
  // remark: has to include pressure, since ScaTraTimIntImpl::SetVelocityField() expects it!
  std::set<XFEM::PHYSICS::Field> outputfields;
  outputfields.insert(XFEM::PHYSICS::Velx);
  outputfields.insert(XFEM::PHYSICS::Vely);
  outputfields.insert(XFEM::PHYSICS::Velz);
  outputfields.insert(XFEM::PHYSICS::Pres);

  // TODO: check performance time to built convel vector; if this is costly, it could be stored as
  //       a private member variable of the time integration scheme, since it is built in two places
  //       (here in every nonlinear iteration, and in the Output() function after every time step)

  // convection velocity vector
  const Teuchos::RCP<Epetra_Vector> convel = dofmanagerForOutput_->transformXFEMtoStandardVector(*state_.velnp_,
                                                               *standarddofset_,
                                                               state_.nodalDofDistributionMap_,
                                                               outputfields);
  return convel;
}

/*------------------------------------------------------------------------------------------------*
 | get convection velocity vector for transfer to scalar transport field           wichmann 02/12 |
 *------------------------------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FLD::CombustFluidImplicitTimeInt::StdVelaf(void)
{
  if (timealgo_ != INPAR::FLUID::timeint_afgenalpha)
    dserror("Velaf is only available for genalpha time int scheme.");

  // velocity vector has to be transformed from XFEM format to Standard FEM format, because ScaTra
  // cannot handle XFEM dofs at enriched nodes
  // remark: has to include pressure, since ScaTraTimIntImpl::SetVelocityField() expects it!
  std::set<XFEM::PHYSICS::Field> outputfields;
  outputfields.insert(XFEM::PHYSICS::Velx);
  outputfields.insert(XFEM::PHYSICS::Vely);
  outputfields.insert(XFEM::PHYSICS::Velz);
  outputfields.insert(XFEM::PHYSICS::Pres);

  // TODO: check performance time to built convel vector; if this is costly, it could be stored as
  //       a private member variable of the time integration scheme, since it is built in two places
  //       (here in every nonlinear iteration, and in the Output() function after every time step)

  // convection velocity vector
  const Teuchos::RCP<Epetra_Vector> convel = dofmanagerForOutput_->transformXFEMtoStandardVector(*state_.velaf_,
                                                               *standarddofset_,
                                                               state_.nodalDofDistributionMap_,
                                                               outputfields);
  return convel;
}

/*------------------------------------------------------------------------------------------------*
 | get history vector for transfer to scalar transport field                      rasthofer 01/10 |
 | needed for subgrid-velocity                                                                    |
 *------------------------------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> FLD::CombustFluidImplicitTimeInt::Hist()
{
  // velocity vector has to be transformed from XFEM format to Standard FEM format, because ScaTra
  // cannot handle XFEM dofs at enriched nodes
  // remark: has to include pressure, since ScaTraTimIntImpl::SetVelocityField() expects it!
  std::set<XFEM::PHYSICS::Field> outputfields;
  outputfields.insert(XFEM::PHYSICS::Velx);
  outputfields.insert(XFEM::PHYSICS::Vely);
  outputfields.insert(XFEM::PHYSICS::Velz);
  outputfields.insert(XFEM::PHYSICS::Pres);

  // TODO: check performance time to built convel vector; if this is costly, it could be stored as
  //       a private member variable of the time integration scheme, since it is built in two places
  //       (here in every nonlinear iteration, and in the Output() function after every time step)

  // convection velocity vector
  Teuchos::RCP<Epetra_Vector> veln = dofmanagerForOutput_->transformXFEMtoStandardVector(
                                         *state_.veln_, *standarddofset_,
                                         state_.nodalDofDistributionMap_, outputfields);
  // acceleration vector
  Teuchos::RCP<Epetra_Vector> accn = dofmanagerForOutput_->transformXFEMtoStandardVector(
                                         *state_.accn_, *standarddofset_,
                                         state_.nodalDofDistributionMap_, outputfields);

  if (veln->MyLength() != accn->MyLength())
    dserror("vectors must have the same length");

  // history vector (OST: linaer combination of veln and accn)
  Teuchos::RCP<Epetra_Vector> hist = LINALG::CreateVector(*standarddofset_->DofRowMap(),true);

  if (hist->MyLength() != accn->MyLength())
    dserror("vectors must have the same length");

  //TODO es ist sowas auch noch in PrepareNonlinearSolve(). Was brauchen wir?
  //stationary case (timealgo_== INPAR::FLUID::timeint_stationary))
  if ( (timealgo_==INPAR::FLUID::timeint_one_step_theta) or
       (timealgo_==INPAR::FLUID::timeint_afgenalpha) )
    FLD::TIMEINT_THETA_BDF2::SetOldPartOfRighthandside(veln,Teuchos::null, accn,timealgo_, dta_, theta_, hist);
  else
    dserror("time integration scheme not supported");

  return hist;
}

/*------------------------------------------------------------------------------------------------*
 | solve the nonlinear fluid problem                                                  henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::NonlinearSolve()
{
  while (!FluidRefSolLoopFinished()) // iterator between NS-solution and recomputation of reference solution
  {
    PrepareNonlinearSolve();

    TEUCHOS_FUNC_TIME_MONITOR("   + nonlin. iteration/lin. solve");

    // ---------------------------------------------- nonlinear iteration
    // ------------------------------- stop nonlinear iteration when both
    //                                 increment-norms are below this bound
    const double  ittol     = params_->get<double>("tolerance for nonlin iter");

    //------------------------------ turn adaptive solver tolerance on/off
    const bool   isadapttol    = params_->get<bool>("ADAPTCONV");
    const double adaptolbetter = params_->get<double>("ADAPTCONV_BETTER");


    // REMARK:
    // commented reduced number of iterations out as it seems that more iterations
    // are necessary before sampling to obtain a converged result
    // for turbulent channel flow: only one iteration before sampling otherwise use the usual itemax_
    //const int itemax = (special_flow_ == "bubbly_channel_flow" and step_ < samstart_) ? 2 : itemax_ ;

    //const bool fluidrobin = params_->get<bool>("fluidrobin", false);

    int               itnum = 0;
    bool              stopnonliniter = false;

    double dtsolve = 0.0;
    dtele_   = 0.0;

    // out to screen
    PrintTimeStepInfo();

    /*
	  {
		std::ofstream f;
		const std::string fname = DRT::Problem::Instance()->OutputControlFile()->FileName()
								+ ".outifacevelnp.txt";
		if (step_ <= 1)
		  f.open(fname.c_str(),std::fstream::trunc);
		else
		  f.open(fname.c_str(),std::fstream::ate | std::fstream::app);

		f << time_ << " " << (*ivelcolnp)[0] << "  " << "\n";

		f.close();
	  }
     */

    if (myrank_ == 0)
    {
      printf("+------------+-------------------+--------------+--------------+--------------+--------------+\n");
      printf("|- step/max -|- tol      [norm] -|-- vel-res ---|-- pre-res ---|-- vel-inc ---|-- pre-inc ---|\n");
    }

    // TODO add comment
    incvel_->PutScalar(0.0);
    residual_->PutScalar(0.0);
    // increment of the old iteration step - used for update of condensed element stresses
    oldinc_ = LINALG::CreateVector(*discret_->DofRowMap(),true);

    while (stopnonliniter==false)
    {
      itnum++;

#ifdef SUGRVEL_OUTPUT
      std::cout << "writing gmsh output" << std::endl;
      const bool screen_out = false;

      const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("SubgridVelocityFluid", step_, 350, screen_out, 0);
      std::ofstream gmshfilecontent(filename.c_str());
      gmshfilecontent << "View \" " << "SubgridVelocity" << " \" {\n";
      gmshfilecontent.close();

      const std::string filename2 = IO::GMSH::GetNewFileNameAndDeleteOldFiles("Residual", step_, 350, screen_out, 0);
      std::ofstream gmshfilecontent2(filename2.c_str());
      gmshfilecontent2 << "View \" " << "Residual" << " \" {\n";
      gmshfilecontent2.close();

      const std::string filename3 = IO::GMSH::GetNewFileNameAndDeleteOldFiles("Tau", step_, 350, screen_out, 0);
      std::ofstream gmshfilecontent3(filename3.c_str());
      gmshfilecontent3 << "View \" " << "Tau" << " \" {\n";
      gmshfilecontent3.close();
#endif

      // -------------------------------------------------------------------
      // call elements to calculate system matrix
      // -------------------------------------------------------------------
      {
        // time measurement: element
        TEUCHOS_FUNC_TIME_MONITOR("      + element calls");

        // get cpu time
        const double tcpu=Teuchos::Time::wallTime();

        sysmat_->Zero();

        // add Neumann loads
        residual_->Update(1.0,*neumann_loads_,0.0);

        // create the parameters for the discretization
        ParameterList eleparams;

        // action for elements
        eleparams.set("action","calc_fluid_systemmat_and_residual");

        // flag for type of combustion problem
        eleparams.set<int>("combusttype",combusttype_);
        eleparams.set<int>("veljumptype",veljumptype_);
        eleparams.set<int>("fluxjumptype",fluxjumptype_);
        eleparams.set("flamespeed",flamespeed_);
        eleparams.set("marksteinlength",marksteinlength_);
        eleparams.set("nitschevel",nitschevel_);
        eleparams.set("nitschepres",nitschepres_);
        eleparams.set("DLM_condensation",condensation_);

        // parameter for suppressing additional enrichment dofs in two-phase flow problems
        eleparams.set<int>("selectedenrichment",DRT::INPUT::IntegralValue<INPAR::COMBUST::SelectedEnrichment>(params_->sublist("COMBUSTION FLUID"),"SELECTED_ENRICHMENT"));

        // parameters for two-phase flow problems with surface tension
        eleparams.set<int>("surftensapprox",surftensapprox_);
        eleparams.set("variablesurftens",params_->sublist("COMBUSTION FLUID").get<double>("VARIABLESURFTENS"));
        eleparams.set("connected_interface",connected_interface_);

        // smoothed normal vectors for boundary integration
        eleparams.set("smoothed_bound_integration",smoothed_boundary_integration_);
        eleparams.set<int>("smoothgradphi",smoothgradphi_);

        // other parameters that might be needed by the elements
        //eleparams.set("total time",time_);
        //eleparams.set("thsl",theta_*dta_);
        eleparams.set<int>("timealgo",timealgo_);
        eleparams.set("time",time_);
        eleparams.set("dt",dta_);
        eleparams.set("theta",theta_);
        eleparams.set("gamma",gamma_);
        eleparams.set("alphaF",alphaF_);
        eleparams.set("alphaM",alphaM_);

#ifdef SUGRVEL_OUTPUT
        //eleparams.set("step",step_);
#endif

        //eleparams.set("include reactive terms for linearisation",params_->get<bool>("Use reaction terms for linearisation"));
        //type of linearisation: include reactive terms for linearisation
        if(DRT::INPUT::get<INPAR::FLUID::LinearisationAction>(*params_, "Linearisation") == INPAR::FLUID::Newton)
          eleparams.set("include reactive terms for linearisation",true);
        else if (DRT::INPUT::get<INPAR::FLUID::LinearisationAction>(*params_, "Linearisation") == INPAR::FLUID::minimal)
          dserror("LinearisationAction minimal is not defined in the combustion formulation");
        else
          eleparams.set("include reactive terms for linearisation",false);

        // parameters for stabilization
        eleparams.sublist("STABILIZATION") = params_->sublist("STABILIZATION");

        // parameters for stabilization
        eleparams.sublist("TURBULENCE MODEL") = params_->sublist("TURBULENCE MODEL");

        // set vector values needed by elements
        discret_->ClearState();

        // set scheme-specific element parameters and vector values
        if (timealgo_==INPAR::FLUID::timeint_afgenalpha)
          discret_->SetState("velaf",state_.velaf_);
        else
          discret_->SetState("velnp",state_.velnp_);

        discret_->SetState("veln" ,state_.veln_);
        discret_->SetState("velnm",state_.velnm_);
        discret_->SetState("accn" ,state_.accn_);

        if (timealgo_==INPAR::FLUID::timeint_afgenalpha)
          discret_->SetState("accam",state_.accam_);

        discret_->SetState("velpres nodal iterinc",oldinc_);

        // convergence check at itemax is skipped for speedup if
        // CONVCHECK is set to L_2_norm_without_residual_at_itemax
        if ((itnum != itemax_) || (params_->get<string>("CONVCHECK") != "L_2_norm_without_residual_at_itemax"))
        {
          // call standard loop over elements
          discret_->Evaluate(eleparams,sysmat_,residual_);

          discret_->ClearState();

                 // account for potential Neumann inflow terms
          if (neumanninflow_)
          {
            // create parameter list
            ParameterList condparams;

            // action for elements
            condparams.set("action","calc_Neumann_inflow");

            // flag for type of combustion problem
            condparams.set<int>("combusttype",combusttype_);
            condparams.set<int>("veljumptype",veljumptype_);
            condparams.set<int>("fluxjumptype",fluxjumptype_);
            condparams.set("nitschevel",nitschevel_);
            condparams.set("nitschepres",nitschepres_);

            // parameter for suppressing additional enrichment dofs in two-phase flow problems
            condparams.set<int>("selectedenrichment",DRT::INPUT::IntegralValue<INPAR::COMBUST::SelectedEnrichment>(params_->sublist("COMBUSTION FLUID"),"SELECTED_ENRICHMENT"));

            // parameters for two-phase flow problems with surface tension
            condparams.set<int>("surftensapprox",surftensapprox_);
            eleparams.set("variablesurftens",params_->sublist("COMBUSTION FLUID").get<double>("VARIABLESURFTENS"));
            condparams.set("connected_interface",connected_interface_);

            // smoothed normal vectors for boundary integration
            condparams.set("smoothed_bound_integration",smoothed_boundary_integration_);
            condparams.set<int>("smoothgradphi",smoothgradphi_);

            // other parameters that might be needed by the elements
            //condparams.set("total time",time_);
            //condparams.set("thsl",theta_*dta_);
            condparams.set<int>("timealgo",timealgo_);
            condparams.set("time",time_);
            condparams.set("dt",dta_);
            condparams.set("theta",theta_);
            condparams.set("gamma",gamma_);
            condparams.set("alphaF",alphaF_);
            condparams.set("alphaM",alphaM_);

            // set vector values needed by elements
            discret_->ClearState();

            // set scheme-specific element parameters and vector values
            if (timealgo_==INPAR::FLUID::timeint_afgenalpha)
              discret_->SetState("velaf",state_.velaf_);
            else discret_->SetState("velnp",state_.velnp_);

            discret_->SetState("veln" ,state_.veln_);
            discret_->SetState("velnm",state_.velnm_);

            std::string condstring("FluidNeumannInflow");
            discret_->EvaluateConditionUsingParentData(condparams,sysmat_,Teuchos::null,residual_,Teuchos::null,Teuchos::null,condstring);
            discret_->ClearState();
          }

          // scaling to get true residual vector for all other schemes
          trueresidual_->Update(ResidualScaling(),*residual_,0.0);

          // finalize the complete matrix
          sysmat_->Complete();
        }

        // end time measurement for element
        dtele_=Teuchos::Time::wallTime()-tcpu;

      } // end of element call

      // blank residual DOFs which are on Dirichlet BC
      // We can do this because the values at the dirichlet positions
      // are not used anyway.
      // We could avoid this though, if velrowmap_ and prerowmap_ would
      // not include the dirichlet values as well. But it is expensive
      // to avoid that.
      dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), residual_);

      // Krylov projection for solver already required in convergence check
      if (project_)
      {
        DRT::Condition* KSPcond=discret_->GetCondition("KrylovSpaceProjection");

        // in this case, we want to project out some zero pressure modes
        const string* definition = KSPcond->Get<string>("weight vector definition");

        if(*definition == "pointvalues")
        {
          // zero w and c
          w_->PutScalar(0.0);
          c_->PutScalar(0.0);

          // get pressure
          const vector<double>* mode = KSPcond->Get<vector<double> >("mode");

          // make sure the dat file specified a pressure mode
          for(int rr=0;rr<numdim_;++rr)
          {
            if(abs((*mode)[rr])>1e-14)
              dserror("expecting only an undetermined pressure");
          }

          int predof = numdim_;

          /* export to vector to normalize against
        //
        // Note that in the case of definition pointvalue based,
        // the average pressure will vanish in a pointwise sense
        //
        //    +---+
        //     \
        //      +   p_i  = 0
        //     /
        //    +---+
           */

          // write the value specified in the .dat-file to every std. pressure dof
          // (regardless of whether this node also uses enriched pressure)
          dofmanagerForOutput_->overwritePhysicalField(w_, state_.nodalDofDistributionMap_, XFEM::PHYSICS::Pres,
              XFEM::Enrichment::typeStandard, (*mode)[predof], false);

          // select only a certain volume to be affected by this projection
          Teuchos::RCP<Epetra_Vector> tmpkspw = kspsplitter_->ExtractKSPCondVector(*w_);
          w_->PutScalar(0.0);
          LINALG::Export(*tmpkspw,*w_);

          // c_ = w_ except we use 1.0 instead of the .dat-file provided value
          dofmanagerForOutput_->overwritePhysicalField(c_, state_.nodalDofDistributionMap_, XFEM::PHYSICS::Pres,
              XFEM::Enrichment::typeStandard, 1.0, false);

          Teuchos::RCP<Epetra_Vector> tmpkspc = kspsplitter_->ExtractKSPCondVector(*c_);
          c_->PutScalar(0.0);
          LINALG::Export(*tmpkspc,*c_);
        }
        else if(*definition == "integration")
        {
          // zero w and c
          w_->PutScalar(0.0);
          c_->PutScalar(0.0);

          //---------------------------------------------------------
          // fill w_
          //---------------------------------------------------------
          ParameterList eleparams;
          // set action for elements
          eleparams.set("action","integrate_Shapefunction");

          /* evaluate KrylovSpaceProjection condition in order to get
        // integrated nodal basis functions w_
        // Note that in the case of definition integration based,
        // the average pressure will vanish in an integral sense
        //
        //                    /              /                      /
        //   /    \          |              |  /          \        |  /    \
        //  | w_*p | = p_i * | N_i(x) dx =  | | N_i(x)*p_i | dx =  | | p(x) | dx = 0
        //   \    /          |              |  \          /        |  \    /
        //                   /              /                      /
           */

          // call loop over elements
          discret_->ClearState();
          discret_->EvaluateCondition(eleparams, w_, "KrylovSpaceProjection");
          discret_->ClearState();
          //----------------------------------------------------------
          // fill c_
          //----------------------------------------------------------
          // same as for pointwise
          dofmanagerForOutput_->overwritePhysicalField(c_, state_.nodalDofDistributionMap_, XFEM::PHYSICS::Pres,
              XFEM::Enrichment::typeStandard, 1.0, false);

          Teuchos::RCP<Epetra_Vector> tmpkspc = kspsplitter_->ExtractKSPCondVector(*c_);
          c_->PutScalar(0.0);
          LINALG::Export(*tmpkspc,*c_);
        }
        else
        {
          dserror("unknown definition of weight vector w for restriction of Krylov space");
        }

        double cTw;
        c_->Dot(*w_,&cTw);

        double cTres;
        c_->Dot(*residual_,&cTres);

        residual_->Update(-cTres/cTw,*w_,1.0);
      }

      double incvelnorm_L2 = 0.0;
      double velnorm_L2 = 0.0;
      double vresnorm = 0.0;

      Teuchos::RCP<Epetra_Vector> onlyvel = velpressplitter_->ExtractOtherVector(residual_);
      onlyvel->Norm2(&vresnorm);

      velpressplitter_->ExtractOtherVector(incvel_,onlyvel);
      onlyvel->Norm2(&incvelnorm_L2);

      velpressplitter_->ExtractOtherVector(state_.velnp_,onlyvel);
      onlyvel->Norm2(&velnorm_L2);

      double incprenorm_L2 = 0.0;
      double prenorm_L2 = 0.0;
      double presnorm = 0.0;

      Teuchos::RCP<Epetra_Vector> onlypre = velpressplitter_->ExtractCondVector(residual_);
      onlypre->Norm2(&presnorm);

      velpressplitter_->ExtractCondVector(incvel_,onlypre);
      onlypre->Norm2(&incprenorm_L2);

      velpressplitter_->ExtractCondVector(state_.velnp_,onlypre);
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
      double min = 1.0e19;
      double max = 0.0;
      discret_->Comm().MinAll(&dtele_,&min,1);
      discret_->Comm().MaxAll(&dtele_,&max,1);

      if (myrank_ == 0)
      {
        printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   |      --      |      --      |",
               itnum,itemax_,ittol,vresnorm,presnorm);
//        printf(" (      --     ,te=%10.3E",dtele);
//        printf(")\n");
        printf(" (      --     ,te_min=%10.3E,te_max=%10.3E)\n", min, max);
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
      if (vresnorm <= ittol and
          presnorm <= ittol and
          incvelnorm_L2/velnorm_L2 <= ittol and
          incprenorm_L2/prenorm_L2 <= ittol)
      {
        stopnonliniter=true;

        double min = 1.0e19;
        double max = 0.0;
        discret_->Comm().MinAll(&dtele_,&min,1);
        discret_->Comm().MaxAll(&dtele_,&max,1);

        if (myrank_ == 0)
        {
          printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   | %10.3E   | %10.3E   |",
                 itnum,itemax_,ittol,vresnorm,presnorm,
                 incvelnorm_L2/velnorm_L2,incprenorm_L2/prenorm_L2);
//          printf(" (ts=%10.3E,te=%10.3E",dtsolve,dtele);
//          printf(")\n");
          printf(" (ts=%10.3E,te_min=%10.3E,te_max=%10.3E)\n", dtsolve, min, max);
          printf("+------------+-------------------+--------------+--------------+--------------+--------------+\n");

          FILE* errfile = params_->get<FILE*>("err file");
          if (errfile!=NULL)
          {
            fprintf(errfile,"fluid solve:   %3d/%3d  tol=%10.3E[L_2 ]  vres=%10.3E  pres=%10.3E  vinc=%10.3E  pinc=%10.3E\n",
                    itnum,itemax_,ittol,vresnorm,presnorm,
                    incvelnorm_L2/velnorm_L2,incprenorm_L2/prenorm_L2);
          }
        }
#ifdef COMBUST_2D
        // for 2-dimensional problems errors occur because of pseudo 3D code!
        // These error shall get smaller with the following modifications
        const int dim = 2; // z-direction is assumed to be the pseudo-dimension
        if (dim == 2)
        {
          const Epetra_Map* dofcolmap = discret_->DofColMap();
          map<XFEM::DofKey,XFEM::DofGID> dofColDistrib;
          dofmanagerForOutput_->fillNodalDofColDistributionMap(dofColDistrib);

          RCP<Epetra_Vector> velnp = rcp(new Epetra_Vector(*dofcolmap,true));
          LINALG::Export(*state_.velnp_,*velnp);


          for (int inode=0;inode<discret_->NumMyColNodes();inode++)
          {
            DRT::Node* frontnode = discret_->lColNode(inode);
            DRT::Node* backnode = NULL;
            if (frontnode->Id()%2==0) // gerade gid -> nachbar hat gid+1
              backnode = discret_->gNode(frontnode->Id()+1);
            else // ungerade gid -> nachbar hat gid-1
              backnode = discret_->gNode(frontnode->Id()-1);

            const std::set<XFEM::FieldEnr>& frontfieldEnrSet(dofmanagerForOutput_->getNodeDofSet(frontnode->Id()));

            { // check if frontnode and backnode is correct
              LINALG::Matrix<3,1> frontcoords(frontnode->X());
              LINALG::Matrix<3,1> backcoords(backnode->X());

              if ((fabs(frontcoords(0) - backcoords(0)) > 1e-12) || (fabs(frontcoords(1) - backcoords(1)) > 1e-12))
              {
                cout << *frontnode << endl;
                cout << *backnode << endl;
                dserror("wrong order of nodes as thought here!");
              }

              // compare both fieldenrsets
              const std::set<XFEM::FieldEnr>& backfieldEnrSet(dofmanagerForOutput_->getNodeDofSet(backnode->Id()));
              int i=0;
              for (set<XFEM::FieldEnr>::const_iterator frontfieldenr = frontfieldEnrSet.begin();
                  frontfieldenr != frontfieldEnrSet.end();frontfieldenr++)
              {
                int j=0;
                for (set<XFEM::FieldEnr>::const_iterator backfieldenr = backfieldEnrSet.begin();
                    backfieldenr != backfieldEnrSet.end();backfieldenr++)
                {
                  if (i==j)
                  {
                    if (*frontfieldenr != *backfieldenr)
                      dserror("fieldenrsets do not fit");
                  }
                  j++;
                }
                i++;
              }
            } // if the compare is successful, all fieldenrichments for this node fit for front and back

            for (set<XFEM::FieldEnr>::const_iterator fieldenr = frontfieldEnrSet.begin();
                fieldenr != frontfieldEnrSet.end();fieldenr++)
            {
              const XFEM::DofKey frontdofkey(frontnode->Id(),*fieldenr);
              const XFEM::DofKey backdofkey(backnode->Id(),*fieldenr);

              const int frontdofpos = dofColDistrib.find(frontdofkey)->second;
              const int backdofpos = dofColDistrib.find(backdofkey)->second;

              if (fieldenr->getField() != XFEM::PHYSICS::Velz)
              {
                double average = 0.5*((*velnp)[(*dofcolmap).LID(frontdofpos)]
                                               +(*velnp)[(*dofcolmap).LID(backdofpos)]);
                (*velnp)[(*dofcolmap).LID(frontdofpos)] = average;
                (*velnp)[(*dofcolmap).LID(backdofpos)] = average;
              }
              else // Velz values shall be set to zero
              {
                (*velnp)[(*dofcolmap).LID(frontdofpos)] = 0;
                (*velnp)[(*dofcolmap).LID(backdofpos)] = 0;
              }
            }
          }
          LINALG::Export(*velnp,*state_.velnp_);
        }
#endif
        break;
      }
      else // if not yet converged
      {
        double min = 1.0e19;
        double max = 0.0;
        discret_->Comm().MinAll(&dtele_,&min,1);
        discret_->Comm().MaxAll(&dtele_,&max,1);

        if (myrank_ == 0)
        {
          printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   | %10.3E   | %10.3E   |",
                 itnum,itemax_,ittol,vresnorm,presnorm,
                 incvelnorm_L2/velnorm_L2,incprenorm_L2/prenorm_L2);
//          printf(" (ts=%10.3E,te=%10.3E",dtsolve,dtele);
//          printf(")\n");
          printf(" (ts=%10.3E,te_min=%10.3E,te_max=%10.3E)\n", dtsolve, min, max);
        }
      }
    }

    // warn if itemax is reached without convergence, but proceed to
    // next timestep...
    if ((itnum == itemax_) and (vresnorm > ittol or presnorm > ittol or
                             incvelnorm_L2/velnorm_L2 > ittol or
                             incprenorm_L2/prenorm_L2 > ittol))
    {
      stopnonliniter=true;
      if (myrank_ == 0)
      {
        printf("+---------------------------------------------------------------+\n");
        printf("|            >>>>>> not converged in itemax steps!              |\n");
        printf("+---------------------------------------------------------------+\n");

        FILE* errfile = params_->get<FILE*>("err file");
        if (errfile!=NULL)
        {
          fprintf(errfile,"fluid unconverged solve:   %3d/%3d  tol=%10.3E[L_2 ]  vres=%10.3E  pres=%10.3E  vinc=%10.3E  pinc=%10.3E\n",
                  itnum,itemax_,ittol,vresnorm,presnorm,
                  incvelnorm_L2/velnorm_L2,incprenorm_L2/prenorm_L2);
        }
      }

#ifdef COMBUST_2D
      // for 2-dimensional problems errors occur because of pseudo 3D code!
      // These error shall be removed with the following modifications
      const int dim = 2; // z-direction is assumed to be the pseudo-dimension
      if (dim == 2)
      {
        const Epetra_Map* dofcolmap = discret_->DofColMap();
        map<XFEM::DofKey,XFEM::DofGID> dofColDistrib;
        dofmanagerForOutput_->fillNodalDofColDistributionMap(dofColDistrib);

        RCP<Epetra_Vector> velnp = rcp(new Epetra_Vector(*dofcolmap,true));
        LINALG::Export(*state_.velnp_,*velnp);


        for (int inode=0;inode<discret_->NumMyColNodes();inode++)
        {
          DRT::Node* frontnode = discret_->lColNode(inode);
          DRT::Node* backnode = NULL;
          if (frontnode->Id()%2==0) // gerade gid -> nachbar hat gid+1
            backnode = discret_->gNode(frontnode->Id()+1);
          else // ungerade gid -> nachbar hat gid-1
            backnode = discret_->gNode(frontnode->Id()-1);

          const std::set<XFEM::FieldEnr>& frontfieldEnrSet(dofmanagerForOutput_->getNodeDofSet(frontnode->Id()));

          { // check if frontnode and backnode is correct
            LINALG::Matrix<3,1> frontcoords(frontnode->X());
            LINALG::Matrix<3,1> backcoords(backnode->X());

            if ((fabs(frontcoords(0) - backcoords(0)) > 1e-12) || (fabs(frontcoords(1) - backcoords(1)) > 1e-12))
            {
              cout << *frontnode << endl;
              cout << *backnode << endl;
              dserror("wrong order of nodes as thought here!");
            }

            // compare both fieldenrsets
            const std::set<XFEM::FieldEnr>& backfieldEnrSet(dofmanagerForOutput_->getNodeDofSet(backnode->Id()));
            int i=0;
            for (set<XFEM::FieldEnr>::const_iterator frontfieldenr = frontfieldEnrSet.begin();
                frontfieldenr != frontfieldEnrSet.end();frontfieldenr++)
            {
              int j=0;
              for (set<XFEM::FieldEnr>::const_iterator backfieldenr = backfieldEnrSet.begin();
                  backfieldenr != backfieldEnrSet.end();backfieldenr++)
              {
                if (i==j)
                {
                  if (*frontfieldenr != *backfieldenr)
                    dserror("fieldenrsets do not fit");
                }
                j++;
              }
              i++;
            }
          } // if the compare is successful, all fieldenrichments for this node fit for front and back

          for (set<XFEM::FieldEnr>::const_iterator fieldenr = frontfieldEnrSet.begin();
              fieldenr != frontfieldEnrSet.end();fieldenr++)
          {
            const XFEM::DofKey frontdofkey(frontnode->Id(),*fieldenr);
            const XFEM::DofKey backdofkey(backnode->Id(),*fieldenr);

            const int frontdofpos = dofColDistrib.find(frontdofkey)->second;
            const int backdofpos = dofColDistrib.find(backdofkey)->second;

            if (fieldenr->getField() != XFEM::PHYSICS::Velz)
            {
              double average = 0.5*((*velnp)[(*dofcolmap).LID(frontdofpos)]
                                             +(*velnp)[(*dofcolmap).LID(backdofpos)]);
              (*velnp)[(*dofcolmap).LID(frontdofpos)] = average;
              (*velnp)[(*dofcolmap).LID(backdofpos)] = average;
            }
            else // Velz values shall be set to zero
            {
              (*velnp)[(*dofcolmap).LID(frontdofpos)] = 0;
              (*velnp)[(*dofcolmap).LID(backdofpos)] = 0;
            }
          }
        }
        LINALG::Export(*velnp,*state_.velnp_);
      }
#endif
      break;
    }

    // stop if NaNs occur
    if (std::isnan(vresnorm) or
        std::isnan(presnorm) or
        std::isnan(incvelnorm_L2/velnorm_L2) or
        std::isnan(incprenorm_L2/prenorm_L2))
    {
      dserror("NaN's detected! Quitting...");
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
      discret_->Comm().Barrier();
      const double tcpusolve=Teuchos::Time::wallTime();

      // do adaptive linear solver tolerance (not in first solve)
      if (isadapttol and itnum>1)
      {
        double currresidual = max(vresnorm,presnorm);
        currresidual = max(currresidual,incvelnorm_L2/velnorm_L2);
        currresidual = max(currresidual,incprenorm_L2/prenorm_L2);
        solver_->AdaptTolerance(ittol,currresidual,adaptolbetter);
      }

//      // print matrix in matlab format
//      // matrix printing options (DEBUGGING!)
//      cout << "print matrix in matlab format to sparsematrix.mtl" << endl;
//      //cast sysmat to SparseMatrix
//      Teuchos::RCP<LINALG::SparseMatrix> A = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(sysmat_);
//      if (A != Teuchos::null)
//      {
//        const std::string fname = "sparsematrix.mtl";
//        cout << "Col:     MinGID " << A->ColMap().MinAllGID() << " MaxGID "   <<   A->ColMap().MaxAllGID()       <<    "   Row:    MinGID    " <<   A->RowMap().MinAllGID() << " MaxGID " << A->RowMap().MaxAllGID()<< endl;
//        // print to  le in matlab format + const std::string fname = "sparsematrix.mtl";
//        LINALG::PrintMatrixInMatlabFormat(fname,*(A->EpetraMatrix()));
//        // print to screen
//        // (A->EpetraMatrix())->Print(cout);
//        // print sparsity pattern to  le
//        //LINALG::PrintSparsityToPostscript( *(A->EpetraMatrix()) );
//      }
//      else
//      {
//        cout << "failure" << endl;
//        //   Teuchos::RCP<LINALG::BlockSparseMatrixBase>            A  =   +//  state_->BlockSystemMatrix();
//        // LINALG::PrintBlockMatrixInMatlabFormat(fname,*(A));
//      }
//      cout << " ...done" << endl;
//      dserror("Fertig");

      solver_->Solve(sysmat_->EpetraOperator(),incvel_,residual_,true,itnum==1,w_,c_,project_); // defaults: weighted_basis_mean_w = null, kernel_c = null, project = false
      solver_->ResetTolerance();

      // end time measurement for solver
      // remark: to get realistic times, this barrier results in the longest time being printed
      discret_->Comm().Barrier();
      dtsolve = Teuchos::Time::wallTime()-tcpusolve;
    }

    //------------------------------------------------ update (u,p) trial
    state_.velnp_->Update(1.0,*incvel_,1.0);

    // -------------------------------------------------------------------
    // For af-generalized-alpha: update accelerations
    // Furthermore, calculate velocities, pressures, scalars and
    // accelerations at intermediate time steps n+alpha_F and n+alpha_M,
    // respectively, for next iteration.
    // This has to be done at the end of the iteration, since we might
    // need the velocities at n+alpha_F in a potential coupling
    // algorithm, for instance.
    // -------------------------------------------------------------------
    if (timealgo_==INPAR::FLUID::timeint_afgenalpha)
    {
      GenAlphaUpdateAcceleration();

      GenAlphaIntermediateValues();
    }

    //------------------- store nodal increment for element stress update
    oldinc_->Update(1.0,*incvel_,0.0);
    }
  }

#ifdef SUGRVEL_OUTPUT
  const bool screen_out = false;

  const std::string filename = IO::GMSH::GetFileName("SubgridVelocityFluid", step_, screen_out, 0);
  std::ofstream gmshfilecontent(filename.c_str(), ios_base::out | ios_base::app);
  gmshfilecontent << "};\n";
  gmshfilecontent.close();

  const std::string filename2 = IO::GMSH::GetFileName("Residual", step_, screen_out, 0);
  std::ofstream gmshfilecontent2(filename2.c_str(), ios_base::out | ios_base::app);
  gmshfilecontent2 << "};\n";
  gmshfilecontent2.close();

  const std::string filename3 = IO::GMSH::GetFileName("Tau", step_, screen_out, 0);
  std::ofstream gmshfilecontent3(filename3.c_str(), ios_base::out | ios_base::app);
  gmshfilecontent3 << "};\n";
  gmshfilecontent3.close();
#endif

  //--------------------
  // compute error norms
  //--------------------
  INPAR::COMBUST::NitscheError errortype = DRT::INPUT::IntegralValue<INPAR::COMBUST::NitscheError>(params_->sublist("COMBUSTION FLUID"),"NITSCHE_ERROR");
  if(errortype != INPAR::COMBUST::nitsche_error_none)
    FLD::CombustFluidImplicitTimeInt::EvaluateErrorComparedToAnalyticalSol_Nitsche(errortype);

  oldinc_ = Teuchos::null;
  incvel_ = Teuchos::null;
  residual_ = Teuchos::null;
  zeros_ = Teuchos::null;
  sysmat_ = Teuchos::null;
}


/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::Evaluate(Teuchos::RCP<const Epetra_Vector> vel)
{
  dserror("no monolithic FSI tested, check first!");
  sysmat_->Zero();
  return;

}

/*------------------------------------------------------------------------------------------------*
 | update a fluid time step                                                           henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::TimeUpdate()
{
  // compute acceleration
  TIMEINT_THETA_BDF2::CalculateAcceleration(
      state_.velnp_, state_.veln_, state_.velnm_, state_.accn_,
          timealgo_, step_, theta_, dta_, dtp_,
          state_.accnp_);

  return;
}// FluidImplicitTimeInt::TimeUpdate


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | update acceleration for generalized-alpha time integration  vg 02/09 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::CombustFluidImplicitTimeInt::GenAlphaUpdateAcceleration()
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
  Teuchos::RCP<Epetra_Vector> onlyaccn  = velpressplitter_->ExtractOtherVector(state_.accn_ );
  Teuchos::RCP<Epetra_Vector> onlyveln  = velpressplitter_->ExtractOtherVector(state_.veln_ );
  Teuchos::RCP<Epetra_Vector> onlyvelnp = velpressplitter_->ExtractOtherVector(state_.velnp_);

  Teuchos::RCP<Epetra_Vector> onlyaccnp = rcp(new Epetra_Vector(onlyaccn->Map()));

  const double fact1 = 1.0/(gamma_*dta_);
  const double fact2 = 1.0 - (1.0/gamma_);
  onlyaccnp->Update(fact2,*onlyaccn,0.0);
  onlyaccnp->Update(fact1,*onlyvelnp,-fact1,*onlyveln,1.0);

  // copy back into global vector
  LINALG::Export(*onlyaccnp,*state_.accnp_);

} // FluidImplicitTimeInt::GenAlphaUpdateAcceleration


/*----------------------------------------------------------------------*
 | compute values at intermediate time steps for gen.-alpha    vg 02/09 |
 *----------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::GenAlphaIntermediateValues()
{
  //       n+alphaM                n+1                      n
  //    acc         = alpha_M * acc     + (1-alpha_M) *  acc
  //       (i)                     (i)
  {
    // extract the degrees of freedom associated with velocities
    // only these are allowed to be updated, otherwise you will
    // run into trouble in loma, where the 'pressure' component
    // is used to store the acceleration of the temperature
    Teuchos::RCP<Epetra_Vector> onlyaccn  = velpressplitter_->ExtractOtherVector(state_.accn_ );
    Teuchos::RCP<Epetra_Vector> onlyaccnp = velpressplitter_->ExtractOtherVector(state_.accnp_);

    Teuchos::RCP<Epetra_Vector> onlyaccam = rcp(new Epetra_Vector(onlyaccnp->Map()));

    onlyaccam->Update((alphaM_),*onlyaccnp,(1.0-alphaM_),*onlyaccn,0.0);

    // copy back into global vector
    LINALG::Export(*onlyaccam,*state_.accam_);
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
  state_.velaf_->Update((alphaF_),*state_.velnp_,(1.0-alphaF_),*state_.veln_,0.0);

} // FluidImplicitTimeInt::GenAlphaIntermediateValues


/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::Output()
{
  const bool write_visualization_data = step_%upres_ == 0;
  const bool write_restart_data = step_!=0 and uprestart_ != 0 and step_%uprestart_ == 0;
  const bool do_time_sample = special_flow_!="no" && step_>=samstart_ && step_<=samstop_;

  // -------------------------------------------------------------------
  //         calculate lift'n'drag forces from the residual
  // -------------------------------------------------------------------
  LiftDrag();

  //-------------------------------------------- output of solution

  if (write_visualization_data or write_restart_data)
  {
    output_->NewStep(step_,time_);
  }

  if (write_visualization_data or do_time_sample)  //write solution for visualization
  {
    std::set<XFEM::PHYSICS::Field> outputfields;
    outputfields.insert(XFEM::PHYSICS::Velx);
    outputfields.insert(XFEM::PHYSICS::Vely);
    outputfields.insert(XFEM::PHYSICS::Velz);
    outputfields.insert(XFEM::PHYSICS::Pres);

    // transform velocity XFEM vector to (standard FEM) output velocity vector
    Teuchos::RCP<Epetra_Vector> velnp_out = dofmanagerForOutput_->transformXFEMtoStandardVector(
        *state_.velnp_, *standarddofset_, state_.nodalDofDistributionMap_, outputfields);

    if (do_time_sample)
    {
      // -------------------------------------------------------------------
      //   add calculated velocity to mean value calculation (statistics)
      // -------------------------------------------------------------------

      // transform residual XFEM vector to (standard FEM) output residual vector
      Teuchos::RCP<Epetra_Vector> trueresidual_out = dofmanagerForOutput_->transformXFEMtoStandardVector(
          *trueresidual_, *standarddofset_, state_.nodalDofDistributionMap_, outputfields);

      turbstatisticsmanager_->DoTimeSample(
          step_,
          velnp_out,
          trueresidual_out,
          phinp_,
          standarddofset_,
          state_.velnp_);
    }

    if (write_visualization_data)
    {

      // write physical fields on full domain including voids etc.
      if (outputfields.find(XFEM::PHYSICS::Velx) != outputfields.end())
     {
        // output velocity field for visualization
        output_->WriteVector("velocity_smoothed", velnp_out);

        // output (hydrodynamic) pressure for visualization
        Teuchos::RCP<Epetra_Vector> pressure = velpressplitterForOutput_->ExtractCondVector(velnp_out);
        output_->WriteVector("pressure_smoothed", pressure);

        //output_->WriteVector("residual", trueresidual_);

        //only perform stress calculation when output is needed
        if (writestresses_)
        {
//          Teuchos::RCP<Epetra_Vector> traction = CalcStresses();
//          output_->WriteVector("traction",traction);
        }
      }
      else if (outputfields.find(XFEM::PHYSICS::Temp) != outputfields.end())
      {
        output_->WriteVector("temperature_smoothed", velnp_out);
      }

      // write domain decomposition for visualization
      output_->WriteElementData();

#if 0
    for (int lnodeid=0; lnodeid < discret_->NumMyRowNodes(); ++lnodeid)
    {
      // get the processor local node
      DRT::Node* lnode = discret_->lRowNode(lnodeid);

      LINALG::Matrix<3,1> nodecoord(true);
      // get physical coordinates of this node
      nodecoord(0) = lnode->X()[0];
      nodecoord(1) = lnode->X()[1];
      nodecoord(2) = lnode->X()[2];

      LINALG::Matrix<3,1> vel(true);
      // get the set of dof IDs for this node (3 x vel + 1 x pressure) from standard FEM dofset
      const std::vector<int> dofids = (*standarddofset_).Dof(lnode);
      std::vector<int> lids(3);
      //cout << "mindist " << std::setw(18) << std::setprecision(12) << std::scientific << mindist << endl;
      if (lnode->Id()==163) //or
          //lnode->Id()==498 or
          //lnode->Id()==726 or
          //lnode->Id()==242)
      {
        cout << "node " << lnode->Id() << " ";
        // extract velocity values (no pressure!) from global velocity vector
        for (int icomp=0; icomp<3; ++icomp)
        {
          lids[icomp] = velnp_out->Map().LID(dofids[icomp]);
          vel(icomp) = (*velnp_out)[lids[icomp]];
          std::cout << std::setw(18)<< std::setprecision(12) <<vel(icomp) << " ";
        }
        cout << "coordinates ";
        for (int icomp=0; icomp<3; ++icomp)
        {
          std::cout << nodecoord(icomp) << " ";
        }
        cout << endl;
      }
    }
        for (int lnodeid=0; lnodeid < discret_->NumMyRowNodes(); ++lnodeid)
    {
      // get the processor local node
      DRT::Node* lnode = discret_->lRowNode(lnodeid);

      LINALG::Matrix<3,1> nodecoord(true);
      // get physical coordinates of this node
      nodecoord(0) = lnode->X()[0];
      nodecoord(1) = lnode->X()[1];
      nodecoord(2) = lnode->X()[2];

      LINALG::Matrix<3,1> vel(true);
      // get the set of dof IDs for this node (3 x vel + 1 x pressure) from standard FEM dofset
      const std::vector<int> dofids = (*standarddofset_).Dof(lnode);
      std::vector<int> lids(3);
      //cout << "mindist " << std::setw(18) << std::setprecision(12) << std::scientific << mindist << endl;
      if (lnode->Id()==149)

      {
        cout << "node " << lnode->Id() << " ";
        // extract velocity values (no pressure!) from global velocity vector
        for (int icomp=0; icomp<3; ++icomp)
        {
          lids[icomp] = velnp_out->Map().LID(dofids[icomp]);
          vel(icomp) = (*velnp_out)[lids[icomp]];
          std::cout << setw(18)<< std::setprecision(12) <<vel(icomp) << " ";
        }
        cout << "coordinates ";
        for (int icomp=0; icomp<3; ++icomp)
        {
          std::cout << nodecoord(icomp) << " ";
        }
        cout << endl;
      }
    }
            for (int lnodeid=0; lnodeid < discret_->NumMyRowNodes(); ++lnodeid)
    {
      // get the processor local node
      DRT::Node* lnode = discret_->lRowNode(lnodeid);

      LINALG::Matrix<3,1> nodecoord(true);
      // get physical coordinates of this node
      nodecoord(0) = lnode->X()[0];
      nodecoord(1) = lnode->X()[1];
      nodecoord(2) = lnode->X()[2];

      LINALG::Matrix<3,1> vel(true);
      // get the set of dof IDs for this node (3 x vel + 1 x pressure) from standard FEM dofset
      const std::vector<int> dofids = (*standarddofset_).Dof(lnode);
      std::vector<int> lids(3);
      //cout << "mindist " << std::setw(18) << std::setprecision(12) << std::scientific << mindist << endl;
      if (lnode->Id()==227)
      {
        cout << "node " << lnode->Id() << " ";
        // extract velocity values (no pressure!) from global velocity vector
        for (int icomp=0; icomp<3; ++icomp)
        {
          lids[icomp] = velnp_out->Map().LID(dofids[icomp]);
          vel(icomp) = (*velnp_out)[lids[icomp]];
          std::cout << setw(18)<< std::setprecision(12) <<vel(icomp) << " ";
        }
        cout << "coordinates ";
        for (int icomp=0; icomp<3; ++icomp)
        {
          std::cout << nodecoord(icomp) << " ";
        }
        cout << endl;
      }
    }
    for (int lnodeid=0; lnodeid < discret_->NumMyRowNodes(); ++lnodeid)
    {
      // get the processor local node
      DRT::Node* lnode = discret_->lRowNode(lnodeid);

      LINALG::Matrix<3,1> nodecoord(true);
      // get physical coordinates of this node
      nodecoord(0) = lnode->X()[0];
      nodecoord(1) = lnode->X()[1];
      nodecoord(2) = lnode->X()[2];

      LINALG::Matrix<3,1> vel(true);
      // get the set of dof IDs for this node (3 x vel + 1 x pressure) from standard FEM dofset
      const std::vector<int> dofids = (*standarddofset_).Dof(lnode);
      std::vector<int> lids(3);
      //cout << "mindist " << std::setw(18) << std::setprecision(12) << std::scientific << mindist << endl;
      if (lnode->Id()==59)
      {
        cout << "node " << lnode->Id() << " ";
        // extract velocity values (no pressure!) from global velocity vector
        for (int icomp=0; icomp<3; ++icomp)
        {
          lids[icomp] = velnp_out->Map().LID(dofids[icomp]);
          vel(icomp) = (*velnp_out)[lids[icomp]];
          std::cout << setw(18)<<  std::setprecision(12) <<vel(icomp) << " ";
        }
        cout << "coordinates ";
        for (int icomp=0; icomp<3; ++icomp)
        {
          std::cout << nodecoord(icomp) << " ";
        }
        cout << endl;
      }
    }

//    for (int lnodeid=0; lnodeid < discret_->NumMyRowNodes(); ++lnodeid)
//    {
//      // get the processor local node
//      DRT::Node* lnode = discret_->lRowNode(lnodeid);
//
//      LINALG::Matrix<3,1> nodecoord(true);
//      // get physical coordinates of this node
//      nodecoord(0) = lnode->X()[0];
//      nodecoord(1) = lnode->X()[1];
//      nodecoord(2) = lnode->X()[2];
//
//      LINALG::Matrix<3,1> vel(true);
//      // get the set of dof IDs for this node (3 x vel + 1 x pressure) from standard FEM dofset
//      const std::vector<int> dofids = (*standarddofset_).Dof(lnode);
//      std::vector<int> lids(3);
//      //cout << "mindist " << std::setw(18) << std::setprecision(12) << std::scientific << mindist << endl;
//      if (lnode->Id()==517 or
//          lnode->Id()==725 or
//          lnode->Id()==451 or
//          lnode->Id()==243)
//      {
//        cout << "node " << lnode->Id() << " ";
//        // extract velocity values (no pressure!) from global velocity vector
//        for (int icomp=0; icomp<3; ++icomp)
//        {
//          lids[icomp] = velnp_out->Map().LID(dofids[icomp]);
//          vel(icomp) = (*velnp_out)[lids[icomp]];
//          std::cout << std::setprecision(12) <<vel(icomp) << " ";
//        }
//        cout << endl;
//      }
//    }
#endif
    }
  }

  // write restart
  if (write_restart_data)
  {
    cout0_ << "---  write restart... " << std::flush;
    //std::cout << state_.velnp_->GlobalLength() << std::endl;
    output_->WriteVector("velnp", state_.velnp_);
    //std::cout << state_.veln_->GlobalLength() << std::endl;
    output_->WriteVector("veln" , state_.veln_);
    //std::cout << state_.velnm_->GlobalLength() << std::endl;
    output_->WriteVector("velnm", state_.velnm_);
    //std::cout << state_.accnp_->GlobalLength() << std::endl;
    output_->WriteVector("accnp", state_.accnp_);
    //std::cout << state_.accn_->GlobalLength() << std::endl;
    output_->WriteVector("accn" , state_.accn_);
    if (timealgo_ == INPAR::FLUID::timeint_afgenalpha)
    {
      //std::cout << state_.velaf_->GlobalLength() << std::endl;
      output_->WriteVector("velaf" , state_.velaf_);
      //std::cout << state_.accam_->GlobalLength() << std::endl;
      output_->WriteVector("accam" , state_.accam_);
    }
    cout0_ << "done" << std::endl;
  }

  //if (step_ % 10 == 0 or step_== 1) //write every 5th time step only
  {
    OutputToGmsh((char*)"solution_field_pressure",(char*)"solution_field_velocity",step_, time_);
#ifdef GMSH_REF_FIELDS
    OutputToGmsh((char*)"mod_start_field_pres",(char*)"mod_start_field_vel",step_, time_);
#endif
  }

  // --------------------------
  // dump turbulence statistics
  // --------------------------
  turbstatisticsmanager_->DoOutput((*output_), step_, 0.0);

//  if (step_%upres_ == 0)  //write solution
//  {
//    output_.NewStep(step_,time_);
//    //output_.WriteVector("velnp", state_.velnp_);
//    std::set<XFEM::PHYSICS::Field> fields_out;
//    fields_out.insert(XFEM::PHYSICS::Velx);
//    fields_out.insert(XFEM::PHYSICS::Vely);
//    fields_out.insert(XFEM::PHYSICS::Velz);
//    fields_out.insert(XFEM::PHYSICS::Pres);
//
//    //----------------------------------------------------------------------------------------------
//    // output velocity vector
//    //----------------------------------------------------------------------------------------------
//    // TODO: check performance time to built convel vector; if this is costly, it could be stored
//    //       as a private member variable of the time integration scheme, since it is built in two
//    //       places (here after every time step and in the ConVelnp() function in every nonlinear
//    //       iteration)
//
//    Teuchos::RCP<Epetra_Vector> velnp_out = Teuchos::null;
//    if (step_ == 0)
//    {
//      std::cout << "output standard velocity vector for time step 0" << std::endl;
//      velnp_out = state_.velnp_;
//    }
//    else
//    {
//      std::cout << "output transformed velocity vector for time step 0" << std::endl;
//      velnp_out = dofmanagerForOutput_->transformXFEMtoStandardVector(
//          *state_.velnp_, *standarddofset_, state_.nodalDofDistributionMap_, fields_out);
//    }
//
//    output_.WriteVector("velnp", velnp_out);
//
//    //----------------------------------------------------------------------------------------------
//    // output (hydrodynamic) pressure vector
//    //----------------------------------------------------------------------------------------------
//    Teuchos::RCP<Epetra_Vector> pressure = velpressplitterForOutput_.ExtractCondVector(velnp_out);
//    output_.WriteVector("pressure", pressure);
//
//    //only perform stress calculation when output is needed
//    if (writestresses_)
//    {
//      dserror("not supported, yet");
//      Teuchos::RCP<Epetra_Vector> traction = CalcStresses();
//      output_.WriteVector("traction",traction);
//    }
//
//    // write domain decomposition for visualization (only once!)
//    if (step_ == upres_)
//     output_.WriteElementData();
//
//    if (uprestart_ != 0 and step_%uprestart_ == 0) //add restart data
//    {
//      output_.WriteVector("accn", state_.accn_);
//      output_.WriteVector("veln", state_.veln_);
//      output_.WriteVector("velnm", state_.velnm_);
//    }
//    else if (physprob_.xfemfieldset_.find(XFEM::PHYSICS::Temp) != physprob_.xfemfieldset_.end())
//    {
//      output_.WriteVector("temperature", velnp_out);
//    }
//  }
//
//  // write restart also when uprestart_ is not a integer multiple of upres_
//  else if (uprestart_ != 0 and step_%uprestart_ == 0)
//  {
//    output_.NewStep    (step_,time_);
//    output_.WriteVector("velnp", state_.velnp_);
//    //output_.WriteVector("residual", trueresidual_);
//
//    //only perform stress calculation when output is needed
//    if (writestresses_)
//    {
//      Teuchos::RCP<Epetra_Vector> traction = CalcStresses();
//      output_.WriteVector("traction",traction);
//    }
//
//    output_.WriteVector("accn", state_.accn_);
//    output_.WriteVector("veln", state_.veln_);
//    output_.WriteVector("velnm", state_.velnm_);
//  }
//
//  if (discret_->Comm().NumProc() == 1)
//  {
//    OutputToGmsh(step_, time_);
//  }
  return;
}

/*------------------------------------------------------------------------------------------------*
 |                                                                                rasthofer 05/10 |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::ReadRestart(int step)
{
  IO::DiscretizationReader reader(discret_,step);
  time_ = reader.ReadDouble("time");
  step_ = reader.ReadInt("step");

  //std::cout << state_.velnp_->GlobalLength() << std::endl;
  //std::cout << state_.veln_->GlobalLength()  << std::endl;
  //std::cout << state_.velnm_->GlobalLength() << std::endl;

  cout0_ << "Read restart" << std::endl;

  reader.ReadVector(state_.velnp_,"velnp");
  //std::cout << state_.velnp_->GlobalLength() << std::endl;
  reader.ReadVector(state_.veln_, "veln");
  //std::cout << state_.veln_->GlobalLength() << std::endl;
  reader.ReadVector(state_.velnm_,"velnm");
  //std::cout << state_.velnm_->GlobalLength() << std::endl;
  reader.ReadVector(state_.accnp_ ,"accnp");
  //std::cout << state_.accnp_->GlobalLength() << std::endl;
  reader.ReadVector(state_.accn_ ,"accn");
  //std::cout << state_.accn_->GlobalLength() << std::endl;
  if (timealgo_ == INPAR::FLUID::timeint_afgenalpha)
  {
    reader.ReadVector(state_.velaf_ ,"velaf");
    //std::cout << state_.velaf_->GlobalLength() << std::endl;
    reader.ReadVector(state_.accam_ ,"accam");
    //std::cout << state_.accam_->GlobalLength() << std::endl;
  }
}

/*------------------------------------------------------------------------------------------------*
 | write output to Gmsh postprocessing files                                          henke 10/09 |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::OutputToGmsh(
    const char* presName,
    const char* velName,
    const int step,
    const double time
    ) const
{
  const bool screen_out = true;

  // get a copy on column parallel distribution
  Teuchos::RCP<const Epetra_Vector> output_col_vel;
  // remark: strcmp returns '0' if 'true'
  if ((strcmp(presName,"reference_field_pressure")==0) ||
      (strcmp(presName,"mod_start_field_pres")==0))// explicit names for old field
  { // 0 means true here!!!
    output_col_vel = DRT::UTILS::GetColVersionOfRowVector(discret_, state_.veln_);
  }
  else // default case
  {
    output_col_vel = DRT::UTILS::GetColVersionOfRowVector(discret_, state_.velnp_);
  }

  //------------------------
  // write pressure solution
  //------------------------
  if (gmshoutput_ and (this->physprob_.xfemfieldset_.find(XFEM::PHYSICS::Pres) != this->physprob_.xfemfieldset_.end()))
  {
    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles(presName, step, 500, screen_out, discret_->Comm().MyPID());
    std::ofstream gmshfilecontent(filename.c_str());
    {
      //---------------------------------
      // write physical pressure solution
      //---------------------------------
      gmshfilecontent << "View \" " << "Pressure Solution (Physical) \" {\n";
      for (int i=0; i<discret_->NumMyRowElements(); ++i)
      {
        const DRT::Element* ele = discret_->lRowElement(i);

        // create local copy of information about dofs
        const XFEM::ElementDofManager eledofman(*ele,physprob_.elementAnsatz_->getElementAnsatz(ele->Shape()),*dofmanagerForOutput_);

        //-------------------------------------------
        // extract pressure values from global vector
        //-------------------------------------------
        vector<int> lm;
        vector<int> lmowner;
        vector<int> lmstride;
        ele->LocationVector(*(discret_), lm, lmowner, lmstride);
        // extract local values from the global vector
        vector<double> myvelnp(lm.size());
        DRT::UTILS::ExtractMyValues(*output_col_vel, myvelnp, lm);

        const int numparam = eledofman.NumDofPerField(XFEM::PHYSICS::Pres);
        const vector<int>& dofpos = eledofman.LocalDofPosPerField(XFEM::PHYSICS::Pres);
        // get pressure values for this element
        LINALG::SerialDenseMatrix presele(1,numparam);
        for (int iparam=0; iparam<numparam; ++iparam)
          presele(0,iparam) = myvelnp[dofpos[iparam]];

        //---------------------------------------------------------------
        // extract local level-set (G-function) values from global vector
        //---------------------------------------------------------------
#ifdef DEBUG
        // get map of this vector
        const Epetra_BlockMap& phimap = phinp_->Map();
        // check, whether this map is still identical with the current node map in the discretization
        if (not phimap.SameAs(*discret_->NodeColMap())) dserror("node column map has changed!");
#endif
        size_t numnode = ele->NumNode();
        vector<double> myphinp(numnode);
        // extract G-function values to element level
        DRT::UTILS::ExtractMyNodeBasedValues(ele, myphinp, *phinp_);

        //----------------------------------------
        // interpolate values from element to cell
        //----------------------------------------
        const GEO::DomainIntCells& domainintcells = interfacehandle_->ElementDomainIntCells(ele->Id());
        for (GEO::DomainIntCells::const_iterator cell = domainintcells.begin(); cell != domainintcells.end(); ++cell)
        {
          // pressure values for this integration cell
          LINALG::SerialDenseMatrix prescell(1,DRT::UTILS::getNumberOfElementNodes(cell->Shape()));

          switch(combusttype_)
          {
          case INPAR::COMBUST::combusttype_premixedcombustion:
          case INPAR::COMBUST::combusttype_twophaseflowjump:
          {
            XFEM::InterpolateCellValuesFromElementValuesLevelSet(*ele, eledofman, *cell, myphinp,
                XFEM::PHYSICS::Pres, presele, prescell);
            break;
          }
          case INPAR::COMBUST::combusttype_twophaseflow:
          {
            XFEM::InterpolateCellValuesFromElementValuesLevelSetKink(*ele, eledofman, *cell, myphinp,
                XFEM::PHYSICS::Pres, presele, prescell);
            break;
          }
          case INPAR::COMBUST::combusttype_twophaseflow_surf:
          {
            // interpolation function for combined enrichments (jumps in pressure and kinks in velocity field)
            XFEM::InterpolateCellValuesFromElementValuesLevelSetKinkJump(*ele, eledofman, *cell, myphinp,
                XFEM::PHYSICS::Pres, presele, prescell);
            break;
          }
          default:
            dserror("unknown type of combustion problem!");
          }

          // copy values from matrix format to vector format
          const size_t numnodes = DRT::UTILS::getNumberOfElementNodes(cell->Shape());
          LINALG::SerialDenseVector prescellvec(numnodes);
          for (size_t inode=0; inode<numnodes; ++inode)
            prescellvec(inode) = prescell(0,inode);

          // write data to Gmsh file
          IO::GMSH::cellWithScalarFieldToStream(cell->Shape(), prescellvec, cell->CellNodalPosXYZ(), gmshfilecontent);
        }
      }
      gmshfilecontent << "};\n";

      //---------------------------------------
      // write pressure jump at Gaussian points
      //---------------------------------------
#ifdef GMSH_WRITE_JUMPS
      if (combusttype_ == INPAR::COMBUST::combusttype_premixedcombustion or
          combusttype_ == INPAR::COMBUST::combusttype_twophaseflowjump)
      {
        gmshfilecontent << "View \" " << "Pressure Jump \" {\n";
#ifdef GMSH_AVERAGE_JUMP
        double presjumpnorm = 0.0;
#endif
        for (int iele=0; iele<discret_->NumMyRowElements(); ++iele)
        {
          const DRT::Element* ele = discret_->lRowElement(iele);

          // output only for bisected elements
          if (interfacehandle_->ElementSplit(ele->Id()))
          {
            //------------------------------------------------------------------------------------------
            // extract local level-set (G-function) values from global vector
            //------------------------------------------------------------------------------------------
#ifdef DEBUG
            // get map of this vector
            const Epetra_BlockMap& phimap = phinp_->Map();
            // check, whether this map is still identical with the current node map in the discretization
            if (not phimap.SameAs(*discret_->NodeColMap())) dserror("node column map has changed!");
#endif
            size_t numnode = ele->NumNode();
            vector<double> myphinp(numnode);
            // extract G-function values to element level
            DRT::UTILS::ExtractMyNodeBasedValues(ele, myphinp, *phinp_);
#ifdef DEBUG
#ifndef COMBUST_HEX20
            if (numnode != 8) dserror("pressure jump output only available for hex8 elements!");
#endif
#endif
#ifndef COMBUST_HEX20
            LINALG::Matrix<8,1> ephi;
#else
            LINALG::Matrix<20,1> ephi;
#endif
            for (size_t iparam=0; iparam<numnode; ++iparam)
              ephi(iparam) = myphinp[iparam];

            // create local copy of information about dofs
            const XFEM::ElementDofManager eledofman(*ele,physprob_.elementAnsatz_->getElementAnsatz(ele->Shape()),*dofmanagerForOutput_);

            //-------------------------------------------
            // extract pressure values from global vector
            //-------------------------------------------
            vector<int> lm;
            vector<int> lmowner;
            vector<int> lmstride;
            ele->LocationVector(*(discret_), lm, lmowner, lmstride);
            // extract local values from the global vector
            vector<double> myvelnp(lm.size());
            DRT::UTILS::ExtractMyValues(*output_col_vel, myvelnp, lm);

            const size_t numparam = eledofman.NumDofPerField(XFEM::PHYSICS::Pres);
            const vector<int>& dofpos = eledofman.LocalDofPosPerField(XFEM::PHYSICS::Pres);
            LINALG::SerialDenseVector epres(numparam);
            for (size_t iparam=0; iparam<numparam; ++iparam)
              epres(iparam) = myvelnp[dofpos[iparam]];

            // get node coordinates for this element
            LINALG::SerialDenseMatrix xyze(3,numnode);
            GEO::fillInitialPositionArray(ele,xyze);

            // evaluate the enrichment function at the interface (boundary integration cells)
            const XFEM::ElementEnrichmentValues enrvals_plus(*ele,eledofman,XFEM::Enrichment::approachFromPlus,ephi);
            const XFEM::ElementEnrichmentValues enrvals_minus(*ele,eledofman,XFEM::Enrichment::approachFromMinus,ephi);

            //--------------------------------------------------------------------
            // interpolate values from element to Gaussian points of boundary cell
            //--------------------------------------------------------------------
            const GEO::BoundaryIntCells& boundaryintcells = interfacehandle_->ElementBoundaryIntCells(ele->Id());
            for (GEO::BoundaryIntCells::const_iterator cell = boundaryintcells.begin(); cell != boundaryintcells.end(); ++cell)
            {
              // choose an (arbitrary) number of Gaussian points for the output
              DRT::UTILS::GaussRule2D intrule2D = DRT::UTILS::intrule2D_undefined;
              if(cell->Shape() == DRT::Element::tri3)
                intrule2D = DRT::UTILS::intrule_tri_3point;
              else if(cell->Shape() == DRT::Element::quad4)
                intrule2D = DRT::UTILS::intrule_quad_4point;
              else
                dserror("Not implemented for this cell type");
              DRT::UTILS::IntegrationPoints2D intpoints(intrule2D);

              // loop over Gaussian points
              for (int iquad=0; iquad<intpoints.nquad; ++iquad)
              {
                // transform coordinates of this Gaussian point
                // coordinates of this integration point in boundary cell coordinates \eta^boundary
                const LINALG::Matrix<2,1> posEtaBoundary(intpoints.qxg[iquad]);
                // coordinates of this integration point in element coordinates \xi^domain
                LINALG::Matrix<3,1> posXiDomain;
                GEO::mapEtaBToXiD(*cell, posEtaBoundary, posXiDomain);

                LINALG::SerialDenseVector funct(numnode,true);
                DRT::UTILS::shape_function_3D(funct,posXiDomain(0),posXiDomain(1),posXiDomain(2),ele->Shape());

                LINALG::SerialDenseVector enrfunct_plus(numparam,true);
                LINALG::SerialDenseVector enrfunct_minus(numparam,true);

                enrvals_plus.ComputeModifiedEnrichedNodalShapefunction(XFEM::PHYSICS::Pres, funct, enrfunct_plus);
                enrvals_minus.ComputeModifiedEnrichedNodalShapefunction(XFEM::PHYSICS::Pres, funct, enrfunct_minus);

                XFEM::ApproxFunc<0,16> shp_jump;
#ifdef DEBUG
#ifndef COMBUST_HEX20
                if (numparam != 16) dserror("pressure jump output only available for fully enriched hex8 elements!");
#endif
#endif
                // fill approximation functions
                for (std::size_t iparam = 0; iparam < 16; ++iparam)
                  shp_jump.d0(iparam) = enrfunct_minus(iparam) - enrfunct_plus(iparam);

                // pressure jump
                double presjump = 0.0;
                for (int iparam = 0; iparam < 16; ++iparam)
                  presjump += epres(iparam)*shp_jump.d0(iparam);

                LINALG::Matrix<3,1> posXYZDomain(true);
                for (size_t inode=0;inode<numnode;++inode)
                {
                  posXYZDomain(0) += funct(inode)*xyze(0,inode);
                  posXYZDomain(1) += funct(inode)*xyze(1,inode);
                  posXYZDomain(2) += funct(inode)*xyze(2,inode);
                }

                // write data to Gmsh file
                IO::GMSH::ScalarToStream(posXYZDomain, presjump, gmshfilecontent);

#ifdef GMSH_AVERAGE_JUMP
                // compute average pressure jump
                {
                  // here, a triangular boundary integration cell is assumed (numvertices = 3)
                  if (cell->Shape() != DRT::Element::tri3) dserror("computation of average pressure jump not implemented for this boundary cell distype");
                  const size_t numvertices = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement;
                  //const size_t numvertices = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement;

                  LINALG::SerialDenseMatrix cellXiDomaintmp = cell->CellNodalPosXiDomain();
                  const LINALG::Matrix<3,numvertices> cellXiDomain(cellXiDomaintmp);

                  const LINALG::Matrix<2,1> gpinEta2D(intpoints.qxg[iquad]);

                  // jacobian for coupled transformation
                  // get derivatives dxi_3D/deta_2D
                  static LINALG::Matrix<2,numvertices> deriv_eta2D;
                  DRT::UTILS::shape_function_2D_deriv1(deriv_eta2D,gpinEta2D(0,0),gpinEta2D(1,0),DRT::Element::tri3);

                  // calculate dxi3Ddeta2D
                  static LINALG::Matrix<3,2> dXi3Ddeta2D;
                  dXi3Ddeta2D.Clear();
                  for (int i = 0; i < 3; i++)   // dimensions
                    for (int j = 0; j < 2; j++) // derivatives
                      for (int k = 0; k < (int)numvertices; k++)
                        dXi3Ddeta2D(i,j) += cellXiDomain(i,k)*deriv_eta2D(j,k);

                  // transform Gauss point to xi3D space (element parameter space)
                  static LINALG::Matrix<3,1> gpinXi3D;
                  gpinXi3D.Clear();
                  // coordinates of this integration point in element coordinates \xi^domain
                  GEO::mapEtaBToXiD(*cell, gpinEta2D, gpinXi3D);

                  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
                  static LINALG::Matrix<3,numnode> deriv_xi3D;
                  DRT::UTILS::shape_function_3D_deriv1(deriv_xi3D,gpinXi3D(0,0), gpinXi3D(1,0), gpinXi3D(2,0), DRT::Element::hex8);


                  // calculate dx3Ddxi3D
                  static LINALG::Matrix<3,3> dX3DdXi3D;
                  dX3DdXi3D.Clear();
                  for (int i = 0; i < 3; i++)   // dimensions
                    for (int j = 0; j < 3; j++) // derivatives
                      for (int k = 0; k < (int)numnode; k++)
                        dX3DdXi3D(i,j) += xyze(i,k)*deriv_xi3D(j,k);

                  // get the coupled Jacobian dx3Ddeta2D
                  static LINALG::Matrix<3,2> dx3Ddeta2D;
                  dx3Ddeta2D.Clear();
                  for (int i = 0; i < 3; i++)   // dimensions
                    for (int j = 0; j < 2; j++) // derivatives
                      for (int k = 0; k < 3; k++)
                        dx3Ddeta2D(i,j) += dX3DdXi3D(i,k) * dXi3Ddeta2D(k,j);

                  // get deformation factor
                  static LINALG::Matrix<2,2> Jac_tmp; // J^T*J
                  Jac_tmp.Clear();
                  Jac_tmp.MultiplyTN(dx3Ddeta2D,dx3Ddeta2D);

                  if(Jac_tmp.Determinant() == 0.0) dserror("deformation factor for boundary integration is zero");

                  const double deform_factor = sqrt(Jac_tmp.Determinant()); // sqrt(det(J^T*J))

                  const double fac = intpoints.qwgt[iquad]*deform_factor;
                  //cout << "fac " << fac << endl;
                  presjumpnorm += fac*presjump*presjump;
                }
#endif // computation average pressure
              }
            }
          }
        }
#ifdef COLLAPSE_FLAME
        cout << endl;
        // get the processor local node
        DRT::Node* lnode = discret_->lRowNode(0);

        LINALG::Matrix<3,1> nodecoord(true);
        // get physical coordinates of this node
        nodecoord(0) = lnode->X()[0];
        nodecoord(1) = lnode->X()[1];
        nodecoord(2) = lnode->X()[2];

        const double zcoord = nodecoord(2);
        const double deltaz = 2.*abs(zcoord);
        cout << "Netz " << 1./deltaz << endl;
        const double pi = atan(1.)*4.;
        const double area = pi*2.*0.25*deltaz;
        const double avpresjump = sqrt(presjumpnorm/area);
        cout << "avpresjump " << avpresjump << endl;
#endif
        gmshfilecontent << "};\n";
      }
#endif
    }
    gmshfilecontent.close();
    if (screen_out) std::cout << " done" << endl;
  }
#if 0
  //---------------------------
  // write temperature solution
  //---------------------------
  if (gmshoutput_ and (this->physprob_.xfemfieldset_.find(XFEM::PHYSICS::Temp) != this->physprob_.xfemfieldset_.end()) )
  {
    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("solution_field_temperature", step, 10, screen_out, discret_->Comm().MyPID());
    std::ofstream gmshfilecontent(filename.c_str());

    {
      const XFEM::PHYSICS::Field field = XFEM::PHYSICS::Temp;
      gmshfilecontent << "View \" " << "Temperature Solution (Physical) \" {\n";
      for (int i=0; i<discret_->NumMyRowElements(); ++i)
      {
        const DRT::Element* ele = discret_->lRowElement(i);

        // create local copy of information about dofs
        const XFEM::ElementDofManager eledofman(*ele,physprob_.elementAnsatz_->getElementAnsatz(ele->Shape()),*dofmanagerForOutput_);

        vector<int> lm;
        vector<int> lmowner;
        vector<int> lmstride;
        ele->LocationVector(*(discret_), lm, lmowner, lmstride);

        // extract local values from the global vector
        vector<double> myvelnp(lm.size());
        DRT::UTILS::ExtractMyValues(*output_col_vel, myvelnp, lm);

        const int numparam = eledofman.NumDofPerField(field);
        const vector<int>& dofpos = eledofman.LocalDofPosPerField(field);

        LINALG::SerialDenseMatrix elementvalues(1,numparam);
        for (int iparam=0; iparam<numparam; ++iparam)
          elementvalues(0,iparam) = myvelnp[dofpos[iparam]];

        //------------------------------------------------------------------------------------------
        // extract local level-set (G-function) values from global vector
        //------------------------------------------------------------------------------------------
#ifdef DEBUG
        // get map of this vector
        const Epetra_BlockMap& phimap = phinp_->Map();
        // check, whether this map is still identical with the current node map in the discretization
        if (not phimap.SameAs(*discret_->NodeColMap())) dserror("node column map has changed!");
#endif

        size_t numnode = ele->NumNode();
        vector<double> myphinp(numnode);
        // extract G-function values to element level
        DRT::UTILS::ExtractMyNodeBasedValues(ele, myphinp, *phinp_);

        const GEO::DomainIntCells& domainintcells =
          dofmanagerForOutput_->getInterfaceHandle()->GetDomainIntCells(ele);
        for (GEO::DomainIntCells::const_iterator cell =
          domainintcells.begin(); cell != domainintcells.end(); ++cell)
        {
          LINALG::SerialDenseMatrix cellvalues(1,DRT::UTILS::getNumberOfElementNodes(cell->Shape()));

          switch(combusttype_)
          {
          case INPAR::COMBUST::combusttype_premixedcombustion:
          case INPAR::COMBUST::combusttype_twophaseflowjump:
          {
            //
            XFEM::InterpolateCellValuesFromElementValuesLevelSet(*ele, eledofman, *cell, myphinp, field,
              elementvalues, cellvalues);
            break;
          }
          case INPAR::COMBUST::combusttype_twophaseflow:
          {
            //
            XFEM::InterpolateCellValuesFromElementValuesLevelSetKink(*ele, eledofman, *cell, myphinp, field,
              elementvalues, cellvalues);
            break;
          }
          case INPAR::COMBUST::combusttype_twophaseflow_surf:
          {
            // schott May 17, 2010
            // there is interpolation function for combined enrichments (jumps in pressure and kinks in velocity field)
            XFEM::InterpolateCellValuesFromElementValuesLevelSetKinkJump(*ele, eledofman, *cell, myphinp, field,
              elementvalues, cellvalues);
            break;
          }
          default:
            dserror("unknown type of combustion problem!");
          }

          const size_t numnodes = DRT::UTILS::getNumberOfElementNodes(cell->Shape());
          LINALG::SerialDenseVector cellvaluesvec(numnodes);
          for (size_t inode=1; inode<numnode; ++inode)
          {
            cellvaluesvec(inode) = cellvalues(0,inode);
          }
          IO::GMSH::cellWithScalarFieldToStream(
              cell->Shape(), cellvaluesvec, cell->CellNodalPosXYZ(), gmshfilecontent);
        }
      }
      gmshfilecontent << "};\n";
    }
    gmshfilecontent.close();
    if (screen_out) std::cout << " done" << endl;
  }
#endif

#if 0
  if (gmshoutput_)
  {
    std::ostringstream filename;
    std::ostringstream filenamedel;
    const std::string filebase = DRT::Problem::Instance()->OutputControlFile()->FileName();
    filename    << filebase << ".solution_field_pressure_disc_" << std::setw(5) << setfill('0') << step   << ".pos";
    filenamedel << filebase << ".solution_field_pressure_disc_" << std::setw(5) << setfill('0') << step-5 << ".pos";
    std::remove(filenamedel.str().c_str());
    if (screen_out) std::cout << "writing " << std::left << std::setw(50) <<filename.str()<<"...";
    std::ofstream gmshfilecontent(filename.str().c_str());

    const XFEM::PHYSICS::Field field = XFEM::PHYSICS::DiscPres;

    {
      gmshfilecontent << "View \" " << "Discontinous Pressure Solution (Physical) \" {\n";
      for (int i=0; i<discret_->NumMyRowElements(); ++i)
      {
        const DRT::Element* ele = discret_->lRowElement(i);

        static LINALG::Matrix<3,27> xyze_xfemElement;
        GEO::fillInitialPositionArray(ele,xyze_xfemElement);

        const map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> element_ansatz(COMBUST::getElementAnsatz(ele->Shape()));

        // create local copy of information about dofs
        const XFEM::ElementDofManager eledofman(*ele,element_ansatz,*dofmanagerForOutput_);

        vector<int> lm;
        vector<int> lmowner;
        ele->LocationVector(*(discret_), lm, lmowner);

        // extract local values from the global vector
        vector<double> myvelnp(lm.size());
        DRT::UTILS::ExtractMyValues(*output_col_vel, myvelnp, lm);

        const int numparam = eledofman.NumDofPerField(field);
        const vector<int>& dofpos = eledofman.LocalDofPosPerField(field);

        LINALG::SerialDenseVector elementvalues(numparam);
        for (int iparam=0; iparam<numparam; ++iparam)
          elementvalues(iparam) = myvelnp[dofpos[iparam]];

        const GEO::DomainIntCells& domainintcells =
          dofmanagerForOutput_->getInterfaceHandle()->GetDomainIntCells(ele);
        for (GEO::DomainIntCells::const_iterator cell =
          domainintcells.begin(); cell != domainintcells.end(); ++cell)
        {
          LINALG::SerialDenseVector cellvalues(DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
          XFEM::computeScalarCellNodeValuesFromElementUnknowns(*ele, dofmanagerForOutput_->getInterfaceHandle(), eledofman,
              *cell, field, elementvalues, cellvalues);
          IO::GMSH::cellWithScalarFieldToStream(
              cell->Shape(), cellvalues, cell->CellNodalPosXYZ(), gmshfilecontent);
        }
      }
      gmshfilecontent << "};\n";
    }
    gmshfilecontent.close();
    if (screen_out) std::cout << " done" << endl;
  }
#endif

#if 0
  if (gmshoutput_)
  {
    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("solution_field_sigma_disc", step, 5, screen_out, discret_->Comm().MyPID());
    cout << endl;
    const std::string filenamexx = IO::GMSH::GetNewFileNameAndDeleteOldFiles("solution_field_sigmaxx_disc", step, 5, screen_out, discret_->Comm().MyPID());
    cout << endl;
    const std::string filenameyy = IO::GMSH::GetNewFileNameAndDeleteOldFiles("solution_field_sigmayy_disc", step, 5, screen_out, discret_->Comm().MyPID());
    cout << endl;
    const std::string filenamezz = IO::GMSH::GetNewFileNameAndDeleteOldFiles("solution_field_sigmazz_disc", step, 5, screen_out, discret_->Comm().MyPID());
    cout << endl;
    const std::string filenamexy = IO::GMSH::GetNewFileNameAndDeleteOldFiles("solution_field_sigmaxy_disc", step, 5, screen_out, discret_->Comm().MyPID());
    cout << endl;
    const std::string filenamexz = IO::GMSH::GetNewFileNameAndDeleteOldFiles("solution_field_sigmaxz_disc", step, 5, screen_out, discret_->Comm().MyPID());
    cout << endl;
    const std::string filenameyz = IO::GMSH::GetNewFileNameAndDeleteOldFiles("solution_field_sigmayz_disc", step, 5, screen_out, discret_->Comm().MyPID());
    std::ofstream gmshfilecontent(  filename.c_str());
    std::ofstream gmshfilecontentxx(filenamexx.c_str());
    std::ofstream gmshfilecontentyy(filenameyy.c_str());
    std::ofstream gmshfilecontentzz(filenamezz.c_str());
    std::ofstream gmshfilecontentxy(filenamexy.c_str());
    std::ofstream gmshfilecontentxz(filenamexz.c_str());
    std::ofstream gmshfilecontentyz(filenameyz.c_str());

    const XFEM::PHYSICS::Field field = XFEM::PHYSICS::Sigmaxx;

    {
      gmshfilecontent   << "View \" " << "Discontinous Stress Solution (Physical) \" {" << endl;
      gmshfilecontentxx << "View \" " << "Discontinous Stress (xx) Solution (Physical) \" {\n";
      gmshfilecontentyy << "View \" " << "Discontinous Stress (yy) Solution (Physical) \" {\n";
      gmshfilecontentzz << "View \" " << "Discontinous Stress (zz) Solution (Physical) \" {\n";
      gmshfilecontentxy << "View \" " << "Discontinous Stress (xy) Solution (Physical) \" {\n";
      gmshfilecontentxz << "View \" " << "Discontinous Stress (xz) Solution (Physical) \" {\n";
      gmshfilecontentyz << "View \" " << "Discontinous Stress (yz) Solution (Physical) \" {\n";
      for (int i=0; i<discret_->NumMyRowElements(); ++i)
      {
        const DRT::Element* ele = discret_->lRowElement(i);

        // create local copy of information about dofs
        const XFEM::ElementDofManager eledofman(*ele,physprob_.elementAnsatz_->getElementAnsatz(ele->Shape()),*dofmanagerForOutput_);

        vector<int> lm;
        vector<int> lmowner;
        ele->LocationVector(*(discret_), lm, lmowner);

        // extract local values from the global vector
        vector<double> myvelnp(lm.size());
        DRT::UTILS::ExtractMyValues(*state_.velnp_, myvelnp, lm);

        const int numparam = eledofman.NumDofPerField(field);
        const vector<int>& dofposxx = eledofman.LocalDofPosPerField(XFEM::PHYSICS::Sigmaxx);
        const vector<int>& dofposyy = eledofman.LocalDofPosPerField(XFEM::PHYSICS::Sigmayy);
        const vector<int>& dofposzz = eledofman.LocalDofPosPerField(XFEM::PHYSICS::Sigmazz);
        const vector<int>& dofposxy = eledofman.LocalDofPosPerField(XFEM::PHYSICS::Sigmaxy);
        const vector<int>& dofposxz = eledofman.LocalDofPosPerField(XFEM::PHYSICS::Sigmaxz);
        const vector<int>& dofposyz = eledofman.LocalDofPosPerField(XFEM::PHYSICS::Sigmayz);

        LINALG::SerialDenseMatrix elementvalues(9,numparam);
        for (int iparam=0; iparam<numparam; ++iparam) elementvalues(0,iparam) = myvelnp[dofposxx[iparam]];
        for (int iparam=0; iparam<numparam; ++iparam) elementvalues(1,iparam) = myvelnp[dofposxy[iparam]];
        for (int iparam=0; iparam<numparam; ++iparam) elementvalues(2,iparam) = myvelnp[dofposxz[iparam]];
        for (int iparam=0; iparam<numparam; ++iparam) elementvalues(3,iparam) = myvelnp[dofposxy[iparam]];
        for (int iparam=0; iparam<numparam; ++iparam) elementvalues(4,iparam) = myvelnp[dofposyy[iparam]];
        for (int iparam=0; iparam<numparam; ++iparam) elementvalues(5,iparam) = myvelnp[dofposyz[iparam]];
        for (int iparam=0; iparam<numparam; ++iparam) elementvalues(6,iparam) = myvelnp[dofposxz[iparam]];
        for (int iparam=0; iparam<numparam; ++iparam) elementvalues(7,iparam) = myvelnp[dofposyz[iparam]];
        for (int iparam=0; iparam<numparam; ++iparam) elementvalues(8,iparam) = myvelnp[dofposzz[iparam]];

        LINALG::SerialDenseVector elementvaluexx(numparam); for (int iparam=0; iparam<numparam; ++iparam) elementvaluexx(iparam) = myvelnp[dofposxx[iparam]];
        LINALG::SerialDenseVector elementvalueyy(numparam); for (int iparam=0; iparam<numparam; ++iparam) elementvalueyy(iparam) = myvelnp[dofposyy[iparam]];
        LINALG::SerialDenseVector elementvaluezz(numparam); for (int iparam=0; iparam<numparam; ++iparam) elementvaluezz(iparam) = myvelnp[dofposzz[iparam]];
        LINALG::SerialDenseVector elementvaluexy(numparam); for (int iparam=0; iparam<numparam; ++iparam) elementvaluexy(iparam) = myvelnp[dofposxy[iparam]];
        LINALG::SerialDenseVector elementvaluexz(numparam); for (int iparam=0; iparam<numparam; ++iparam) elementvaluexz(iparam) = myvelnp[dofposxz[iparam]];
        LINALG::SerialDenseVector elementvalueyz(numparam); for (int iparam=0; iparam<numparam; ++iparam) elementvalueyz(iparam) = myvelnp[dofposyz[iparam]];


        const GEO::DomainIntCells& domainintcells =
          dofmanagerForOutput_->getInterfaceHandle()->GetDomainIntCells(ele);
        for (GEO::DomainIntCells::const_iterator cell =
          domainintcells.begin(); cell != domainintcells.end(); ++cell)
        {
          //LINALG::SerialDenseMatrix xyze_cell(3, cell->NumNode());
          //cell->NodalPosXYZ(*ele, xyze_cell);
          // TODO remove
          const LINALG::SerialDenseMatrix& xyze_cell = cell->CellNodalPosXYZ();

          {
          LINALG::SerialDenseMatrix cellvalues(9,DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
          XFEM::computeTensorCellNodeValuesFromElementUnknowns(*ele, dofmanagerForOutput_->getInterfaceHandle(), eledofman,
              *cell, field, elementvalues, cellvalues);
           IO::GMSH::cellWithTensorFieldToStream(cell->Shape(), cellvalues, xyze_cell, gmshfilecontent);
          }

          {
          LINALG::SerialDenseVector cellvaluexx(DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
          XFEM::computeScalarCellNodeValuesFromElementUnknowns(*ele, dofmanagerForOutput_->getInterfaceHandle(), eledofman, *cell, field, elementvaluexx, cellvaluexx);
          IO::GMSH::cellWithScalarFieldToStream(cell->Shape(), cellvaluexx, xyze_cell, gmshfilecontentxx);
          }
          {
          LINALG::SerialDenseVector cellvalueyy(DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
          XFEM::computeScalarCellNodeValuesFromElementUnknowns(*ele, dofmanagerForOutput_->getInterfaceHandle(), eledofman, *cell, field, elementvalueyy, cellvalueyy);
          IO::GMSH::cellWithScalarFieldToStream(cell->Shape(), cellvalueyy, xyze_cell, gmshfilecontentyy);
          }
          {
          LINALG::SerialDenseVector cellvaluezz(DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
          XFEM::computeScalarCellNodeValuesFromElementUnknowns(*ele, dofmanagerForOutput_->getInterfaceHandle(), eledofman, *cell, field, elementvaluezz, cellvaluezz);
          IO::GMSH::cellWithScalarFieldToStream(cell->Shape(), cellvaluezz, xyze_cell, gmshfilecontentzz);
          }
          {
          LINALG::SerialDenseVector cellvaluexy(DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
          XFEM::computeScalarCellNodeValuesFromElementUnknowns(*ele, dofmanagerForOutput_->getInterfaceHandle(), eledofman, *cell, field, elementvaluexy, cellvaluexy);
          IO::GMSH::cellWithScalarFieldToStream(cell->Shape(), cellvaluexy, xyze_cell, gmshfilecontentxy);
          }
          {
          LINALG::SerialDenseVector cellvaluexz(DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
          XFEM::computeScalarCellNodeValuesFromElementUnknowns(*ele, dofmanagerForOutput_->getInterfaceHandle(), eledofman, *cell, field, elementvaluexz, cellvaluexz);
          IO::GMSH::cellWithScalarFieldToStream(cell->Shape(), cellvaluexz, xyze_cell, gmshfilecontentxz);
          }
          {
          LINALG::SerialDenseVector cellvalueyz(DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
          XFEM::computeScalarCellNodeValuesFromElementUnknowns(*ele, dofmanagerForOutput_->getInterfaceHandle(), eledofman, *cell, field, elementvalueyz, cellvalueyz);
          IO::GMSH::cellWithScalarFieldToStream(cell->Shape(), cellvalueyz, xyze_cell, gmshfilecontentyz);
          }
        }
      }
      gmshfilecontent   << "};\n";
      gmshfilecontentxx << "};\n";
      gmshfilecontentyy << "};\n";
      gmshfilecontentzz << "};\n";
      gmshfilecontentxy << "};\n";
      gmshfilecontentxz << "};\n";
      gmshfilecontentyz << "};\n";
    }
    if (screen_out) std::cout << " done" << endl;
  }
#endif

  //------------------------
  // write velocity solution
  //------------------------
  if (this->physprob_.xfemfieldset_.find(XFEM::PHYSICS::Velx) != this->physprob_.xfemfieldset_.end())
  {
    PlotVectorFieldToGmsh(output_col_vel, velName,"Velocity Solution (Physical) n+1",true, step, time);
    if (timealgo_ != INPAR::FLUID::timeint_stationary)
    {
//      PlotVectorFieldToGmsh(DRT::UTILS::GetColVersionOfRowVector(discret_, state_.velnm_), "solution_field_velocity_nm","Velocity Solution (Physical) n-1",false, step, time);
//      PlotVectorFieldToGmsh(DRT::UTILS::GetColVersionOfRowVector(discret_, state_.veln_), "solution_field_velocity_n","Velocity Solution (Physical) n",false, step, time);
//      PlotVectorFieldToGmsh(DRT::UTILS::GetColVersionOfRowVector(discret_, state_.accn_), "solution_field_acceleration_n","Acceleration Solution (Physical) n",false, step, time);
//      PlotVectorFieldToGmsh(DRT::UTILS::GetColVersionOfRowVector(discret_, state_.accnp_), "solution_field_acceleration_np","Acceleration Solution (Physical) n+1",false, step, time);
    }
  }
}


/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::PlotVectorFieldToGmsh(
    const Teuchos::RCP<const Epetra_Vector>   vectorfield,
    const std::string filestr,
    const std::string name_in_gmsh,
    const bool plot_to_gnuplot,
    const int step,
    const double time
    ) const
{
  const bool screen_out = true;

  if (gmshoutput_)
  {
    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles(filestr, step, 500, screen_out, discret_->Comm().MyPID());
    std::ofstream gmshfilecontent(filename.c_str());

    {
      gmshfilecontent << "View \" " << name_in_gmsh << "\" {\n";
      for (int i=0; i<discret_->NumMyRowElements(); ++i)
      {
        const DRT::Element* ele = discret_->lRowElement(i);

        // create local copy of information about dofs
#ifdef COMBUST_STRESS_BASED
#ifdef COMBUST_EPSPRES_BASED
        const COMBUST::EpsilonPressureAnsatz elementAnsatz;
#endif
#ifdef COMBUST_SIGMA_BASED
        const COMBUST::CauchyStressAnsatz elementAnsatz;
#endif
#else
        // just define a default element ansatz; it is not used anyway
        const COMBUST::TauPressureAnsatz elementAnsatz;
#endif
        const XFEM::ElementDofManager eledofman(*ele,elementAnsatz.getElementAnsatz(ele->Shape()),*dofmanagerForOutput_);

        vector<int> lm;
        vector<int> lmowner;
        vector<int> lmstride;
        ele->LocationVector(*discret_, lm, lmowner, lmstride);

        // extract local values from the global vector
        vector<double> myvelnp(lm.size());
        DRT::UTILS::ExtractMyValues(*vectorfield, myvelnp, lm);

        const vector<int>& dofposvelx = eledofman.LocalDofPosPerField(XFEM::PHYSICS::Velx);
        const vector<int>& dofposvely = eledofman.LocalDofPosPerField(XFEM::PHYSICS::Vely);
        const vector<int>& dofposvelz = eledofman.LocalDofPosPerField(XFEM::PHYSICS::Velz);
#ifdef COMBUST_NORMAL_ENRICHMENT
        const vector<int>& dofposveln = eledofman.LocalDofPosPerField(XFEM::PHYSICS::Veln);
#endif

        const int numparamvelx = eledofman.NumDofPerField(XFEM::PHYSICS::Velx);
        LINALG::SerialDenseMatrix elementvalues(4, numparamvelx,true);
        for (int iparam=0; iparam<numparamvelx; ++iparam)
        {
          elementvalues(0, iparam) = myvelnp[dofposvelx[iparam]];
          elementvalues(1, iparam) = myvelnp[dofposvely[iparam]];
          elementvalues(2, iparam) = myvelnp[dofposvelz[iparam]];
        }
#ifdef COMBUST_NORMAL_ENRICHMENT
        const int numparamveln = eledofman.NumDofPerField(XFEM::PHYSICS::Veln);
        for (int iparam=0; iparam<numparamveln; ++iparam)
        {
          elementvalues(3, iparam) = myvelnp[dofposveln[iparam]];
        }
#endif

        //------------------------------------------------------------------------------------------
        // extract local level-set (G-function) values from global vector
        //------------------------------------------------------------------------------------------
#ifdef DEBUG
        // get map of this vector
        const Epetra_BlockMap& phimap = phinp_->Map();
        // check, whether this map is still identical with the current node map in the discretization
        if (not phimap.SameAs(*discret_->NodeColMap())) dserror("node column map has changed!");
#endif
        size_t numnode = ele->NumNode();
        vector<double> myphinp(numnode);
        // extract G-function values to element level
        DRT::UTILS::ExtractMyNodeBasedValues(ele, myphinp, *phinp_);

#ifdef COMBUST_NORMAL_ENRICHMENT
        // get pointer to vector holding smoothed  G-function gradient values at the fluid nodes
        const Teuchos::RCP<Epetra_MultiVector> gradphinp = flamefront_->GradPhi();
        std::vector<double> mygradphi;
        DRT::UTILS::ExtractMyNodeBasedValues(ele, mygradphi,*gradphinp);

        if (numnode != 8) dserror("only available for hex8 elements!");
        LINALG::Matrix<3,8> egradphi(true);

        unsigned ipos;
        for (size_t inode=0; inode<numnode; ++inode)
        {
          ipos = inode*3;
          egradphi(0, inode) = mygradphi[ipos  ];
          egradphi(1, inode) = mygradphi[ipos+1];
          egradphi(2, inode) = mygradphi[ipos+2];
        }
#endif

        const GEO::DomainIntCells& domainintcells = interfacehandle_->ElementDomainIntCells(ele->Id());
        for (GEO::DomainIntCells::const_iterator cell = domainintcells.begin(); cell != domainintcells.end(); ++cell)
        {
          LINALG::SerialDenseMatrix cellvalues(3, DRT::UTILS::getNumberOfElementNodes(cell->Shape()));

          switch(combusttype_)
          {
          case INPAR::COMBUST::combusttype_premixedcombustion:
          case INPAR::COMBUST::combusttype_twophaseflowjump:
          {
#ifndef COMBUST_NORMAL_ENRICHMENT
            //
            XFEM::InterpolateCellValuesFromElementValuesLevelSet(*ele, eledofman, *cell, myphinp, XFEM::PHYSICS::Velx,
              elementvalues, cellvalues);
#else
            XFEM::InterpolateCellValuesFromElementValuesLevelSetNormal(*ele, eledofman, *cell, myphinp, egradphi, XFEM::PHYSICS::Velx,
              elementvalues, cellvalues);
#endif
            break;
          }
          case INPAR::COMBUST::combusttype_twophaseflow:
          {
            //
            XFEM::InterpolateCellValuesFromElementValuesLevelSetKink(*ele, eledofman, *cell, myphinp, XFEM::PHYSICS::Velx,
              elementvalues, cellvalues);
            break;
          }
          case INPAR::COMBUST::combusttype_twophaseflow_surf:
          {
            // plots for Velx, Vely, Velz ... (Kink enrichment)
            XFEM::InterpolateCellValuesFromElementValuesLevelSetKinkJump(*ele, eledofman, *cell, myphinp, XFEM::PHYSICS::Velx,
              elementvalues, cellvalues);
            break;
          }
          default:
            dserror("unknown type of combustion problem!");
          }

          IO::GMSH::cellWithVectorFieldToStream(cell->Shape(), cellvalues, cell->CellNodalPosXYZ(), gmshfilecontent);
        }
      }
      gmshfilecontent << "};\n";
#ifdef GMSH_WRITE_JUMPS
      //---------------------------------------
      // write velocity jump at Gaussian points
      //---------------------------------------
      if (combusttype_ == INPAR::COMBUST::combusttype_premixedcombustion or
          combusttype_ == INPAR::COMBUST::combusttype_twophaseflowjump)
      {
#ifdef GMSH_AVERAGE_JUMP
        double veljumpnormsquare = 0.0;
#endif
        gmshfilecontent << "View \" " << "Velocity Jump \" {\n";
        for (int iele=0; iele<discret_->NumMyRowElements(); ++iele)
        {
          const DRT::Element* ele = discret_->lRowElement(iele);

          // output only for bisected elements
          if (interfacehandle_->ElementSplit(ele->Id()))
          {
            //------------------------------------------------------------------------------------------
            // extract local level-set (G-function) values from global vector
            //------------------------------------------------------------------------------------------
#ifdef DEBUG
            // get map of this vector
            const Epetra_BlockMap& phimap = phinp_->Map();
            // check, whether this map is still identical with the current node map in the discretization
            if (not phimap.SameAs(*discret_->NodeColMap())) dserror("node column map has changed!");
#endif
            size_t numnode = ele->NumNode();
            vector<double> myphinp(numnode);
            // extract G-function values to element level
            DRT::UTILS::ExtractMyNodeBasedValues(ele, myphinp, *phinp_);

#ifdef DEBUG
#ifndef COMBUST_HEX20
            if (numnode != 8) dserror("velocity jump output only available for hex8 elements!");
#endif
#endif
#ifndef COMBUST_HEX20
            LINALG::Matrix<8,1> ephi;
#else
            LINALG::Matrix<20,1> ephi;
#endif
            for (size_t iparam=0; iparam<numnode; ++iparam)
              ephi(iparam) = myphinp[iparam];

#ifdef COMBUST_NORMAL_ENRICHMENT
            // get pointer to vector holding smoothed  G-function gradient values at the fluid nodes
            const Teuchos::RCP<Epetra_MultiVector> gradphinp = flamefront_->GradPhi();
            std::vector<double> mygradphi;
            DRT::UTILS::ExtractMyNodeBasedValues(ele, mygradphi,*gradphinp);

            if (numnode != 8) dserror("only available for hex8 elements!");
            LINALG::Matrix<3,8> egradphi(true);

            unsigned ipos;
            for (size_t inode=0; inode<numnode; ++inode)
            {
              ipos = inode*3;
              egradphi(0, inode) = mygradphi[ipos  ];
              egradphi(1, inode) = mygradphi[ipos+1];
              egradphi(2, inode) = mygradphi[ipos+2];
            }
#endif
            // create local copy of information about dofs
            const XFEM::ElementDofManager eledofman(*ele,physprob_.elementAnsatz_->getElementAnsatz(ele->Shape()),*dofmanagerForOutput_);

            //-------------------------------------------
            // extract velocity values from global vector
            //-------------------------------------------
            vector<int> lm;
            vector<int> lmowner;
            vector<int> lmstride;
            ele->LocationVector(*(discret_), lm, lmowner, lmstride);
            // extract local values from the global vector
            vector<double> myvelnp(lm.size());
            DRT::UTILS::ExtractMyValues(*vectorfield, myvelnp, lm);

            const vector<int>& dofposvelx = eledofman.LocalDofPosPerField(XFEM::PHYSICS::Velx);
            const vector<int>& dofposvely = eledofman.LocalDofPosPerField(XFEM::PHYSICS::Vely);
            const vector<int>& dofposvelz = eledofman.LocalDofPosPerField(XFEM::PHYSICS::Velz);
#ifdef COMBUST_NORMAL_ENRICHMENT
            const vector<int>& dofposveln = eledofman.LocalDofPosPerField(XFEM::PHYSICS::Veln);
#endif

            const size_t numparamvelx = eledofman.NumDofPerField(XFEM::PHYSICS::Velx);
            LINALG::SerialDenseMatrix evel(4, numparamvelx,true);
            for (size_t iparam=0; iparam<numparamvelx; ++iparam)
            {
              evel(0,iparam) = myvelnp[dofposvelx[iparam]];
              evel(1,iparam) = myvelnp[dofposvely[iparam]];
              evel(2,iparam) = myvelnp[dofposvelz[iparam]];
            }
#ifdef COMBUST_NORMAL_ENRICHMENT
            const size_t numparamveln = eledofman.NumDofPerField(XFEM::PHYSICS::Veln);
            for (size_t iparam=0; iparam<numparamveln; ++iparam)
            {
              evel(3, iparam) = myvelnp[dofposveln[iparam]];
            }
#endif
            // get node coordinates for this element
            LINALG::SerialDenseMatrix xyze(3,numnode);
            GEO::fillInitialPositionArray(ele,xyze);

            // evaluate the enrichment function at the interface (boundary integration cells)
            const XFEM::ElementEnrichmentValues enrvals_plus(*ele,eledofman,XFEM::Enrichment::approachFromPlus,ephi);
            const XFEM::ElementEnrichmentValues enrvals_minus(*ele,eledofman,XFEM::Enrichment::approachFromMinus,ephi);

            //--------------------------------------------------------------------
            // interpolate values from element to Gaussian points of boundary cell
            //--------------------------------------------------------------------
            const GEO::BoundaryIntCells& boundaryintcells = interfacehandle_->ElementBoundaryIntCells(ele->Id());
            for (GEO::BoundaryIntCells::const_iterator cell = boundaryintcells.begin(); cell != boundaryintcells.end(); ++cell)
            {
              // choose an (arbitrary) number of Gaussian points for the output
              DRT::UTILS::GaussRule2D intrule2D = DRT::UTILS::intrule2D_undefined;
              if(cell->Shape() == DRT::Element::tri3)
                intrule2D = DRT::UTILS::intrule_tri_3point;
              else if(cell->Shape() == DRT::Element::quad4)
                intrule2D = DRT::UTILS::intrule_quad_4point;
              else
                dserror("Not implemented for this cell type");
              // choose an (arbitrary) number of Gaussian points for the output
              const DRT::UTILS::IntegrationPoints2D intpoints(intrule2D);

              // loop over Gaussian points
              for (int iquad=0; iquad<intpoints.nquad; ++iquad)
              {
                // transform coordinates of this Gaussian point
                // coordinates of this integration point in boundary cell coordinates \eta^boundary
                const LINALG::Matrix<2,1> posEtaBoundary(intpoints.qxg[iquad]);
                // coordinates of this integration point in element coordinates \xi^domain
                LINALG::Matrix<3,1> posXiDomain;
                GEO::mapEtaBToXiD(*cell, posEtaBoundary, posXiDomain);

                LINALG::SerialDenseVector funct(numnode,true);
                DRT::UTILS::shape_function_3D(funct,posXiDomain(0),posXiDomain(1),posXiDomain(2),ele->Shape());

#ifdef DEBUG
#ifndef COMBUST_NORMAL_ENRICHMENT
                if (numparamvelx != 16) dserror("velocity jump output only available for fully enriched hex8 elements!");
#else
                if (numparamvelx != 8 or numparamveln != 8) dserror("velocity jump output only available for fully enriched hex8 elements!");
#endif
#endif

#ifndef COMBUST_NORMAL_ENRICHMENT
                LINALG::SerialDenseVector enrfunct_plus(numparamvelx,true);
                LINALG::SerialDenseVector enrfunct_minus(numparamvelx,true);

                enrvals_plus.ComputeModifiedEnrichedNodalShapefunction(XFEM::PHYSICS::Velx, funct, enrfunct_plus);
                enrvals_minus.ComputeModifiedEnrichedNodalShapefunction(XFEM::PHYSICS::Velx, funct, enrfunct_minus);

                XFEM::ApproxFunc<0,16> shp_jump;
                // fill approximation functions
                for (std::size_t iparam = 0; iparam < 16; ++iparam)
                  shp_jump.d0(iparam) = enrfunct_minus(iparam) - enrfunct_plus(iparam);

                // velocity jump
                LINALG::Matrix<3,1> veljump(true);
                for (int iparam = 0; iparam < 16; ++iparam)
                {
                  veljump(0) += evel(0,iparam)*shp_jump.d0(iparam);
                  veljump(1) += evel(1,iparam)*shp_jump.d0(iparam);
                  veljump(2) += evel(2,iparam)*shp_jump.d0(iparam);
                }
#else
#ifdef COLLAPSE_FLAME
                LINALG::Matrix<3,1> normal(true);
                for (unsigned i=0;i<numnode;i++)
                {
                  normal(0) += funct(i)*xyze(0,i);
                  normal(1) += funct(i)*xyze(1,i);
                }
                const double norm = normal.Norm2(); // sqrt(normal(0)*normal(0) + normal(1)*normal(1) + normal(2)*normal(2))
                if (norm == 0.0) dserror("norm of normal vector is zero!");
                normal.Scale(1.0/norm);
#endif
                // temporary arrays holding enriched shape functions (N * \Psi) on either side of the interface
                XFEM::ApproxFuncNormalVector<0,8> enrfunct_plus(true);
                XFEM::ApproxFuncNormalVector<0,8> enrfunct_minus(true);

                enrvals_plus.ComputeNormalShapeFunction(funct, egradphi,
#ifdef COLLAPSE_FLAME
                    normal,
#endif
                    enrfunct_plus);
                enrvals_minus.ComputeNormalShapeFunction(funct, egradphi,
#ifdef COLLAPSE_FLAME
                    normal,
#endif
                    enrfunct_minus);

                // remark: initialization is essentiual here, since standard shape functions must be zero!
                XFEM::ApproxFuncNormalVector<0,8> shp_jump(true);     // [[ ]] notation
                for (size_t iparam = 0; iparam < numparamveln; ++iparam)
                {
                  shp_jump.velx.d0.n(iparam) = enrfunct_minus.velx.d0.n(iparam) - enrfunct_plus.velx.d0.n(iparam);
                  shp_jump.vely.d0.n(iparam) = enrfunct_minus.vely.d0.n(iparam) - enrfunct_plus.vely.d0.n(iparam);
                  shp_jump.velz.d0.n(iparam) = enrfunct_minus.velz.d0.n(iparam) - enrfunct_plus.velz.d0.n(iparam);
                }
                // velocity jump
                LINALG::Matrix<3,1> veljump(true);
                veljump = XFEM::interpolateVectorFieldToIntPointNormal(evel, shp_jump, numparamvelx, ele, eledofman);
#endif

                // norm of velocity jump
                const double veljumpnorm = veljump.Norm2();

                LINALG::Matrix<3,1> posXYZDomain(true);
                for (size_t inode=0;inode<numnode;++inode)
                {
                  posXYZDomain(0) += funct(inode)*xyze(0,inode);
                  posXYZDomain(1) += funct(inode)*xyze(1,inode);
                  posXYZDomain(2) += funct(inode)*xyze(2,inode);
                }

                // write data to Gmsh file
                //IO::GMSH::VectorToStream(posXYZDomain, veljump, gmshfilecontent);
                IO::GMSH::ScalarToStream(posXYZDomain, veljumpnorm, gmshfilecontent);

#ifdef GMSH_AVERAGE_JUMP
                // compute average velocity jump
                {
                  // here, a triangular boundary integration cell is assumed (numvertices = 3)
                  if (cell->Shape() != DRT::Element::tri3) dserror("Not implemented for this cell type");
                  const size_t numvertices = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement;
                  //const size_t numvertices = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement;

                  LINALG::SerialDenseMatrix cellXiDomaintmp = cell->CellNodalPosXiDomain();
                  const LINALG::Matrix<3,numvertices> cellXiDomain(cellXiDomaintmp);

                  const LINALG::Matrix<2,1> gpinEta2D(intpoints.qxg[iquad]);

                  // jacobian for coupled transformation
                  // get derivatives dxi_3D/deta_2D
                  static LINALG::Matrix<2,numvertices> deriv_eta2D;
                  DRT::UTILS::shape_function_2D_deriv1(deriv_eta2D,gpinEta2D(0,0),gpinEta2D(1,0),DRT::Element::tri3);

                  // calculate dxi3Ddeta2D
                  static LINALG::Matrix<3,2> dXi3Ddeta2D;
                  dXi3Ddeta2D.Clear();
                  for (int i = 0; i < 3; i++)   // dimensions
                    for (int j = 0; j < 2; j++) // derivatives
                      for (int k = 0; k < (int)numvertices; k++)
                        dXi3Ddeta2D(i,j) += cellXiDomain(i,k)*deriv_eta2D(j,k);

                  // transform Gauss point to xi3D space (element parameter space)
                  static LINALG::Matrix<3,1> gpinXi3D;
                  gpinXi3D.Clear();
                  // coordinates of this integration point in element coordinates \xi^domain
                  GEO::mapEtaBToXiD(*cell, gpinEta2D, gpinXi3D);

                  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
                  static LINALG::Matrix<3,numnode> deriv_xi3D;
                  DRT::UTILS::shape_function_3D_deriv1(deriv_xi3D,gpinXi3D(0,0), gpinXi3D(1,0), gpinXi3D(2,0), DRT::Element::hex8);


                  // calculate dx3Ddxi3D
                  static LINALG::Matrix<3,3> dX3DdXi3D;
                  dX3DdXi3D.Clear();
                  for (int i = 0; i < 3; i++)   // dimensions
                    for (int j = 0; j < 3; j++) // derivatives
                      for (int k = 0; k < (int)numnode; k++)
                        dX3DdXi3D(i,j) += xyze(i,k)*deriv_xi3D(j,k);

                  // get the coupled Jacobian dx3Ddeta2D
                  static LINALG::Matrix<3,2> dx3Ddeta2D;
                  dx3Ddeta2D.Clear();
                  for (int i = 0; i < 3; i++)   // dimensions
                    for (int j = 0; j < 2; j++) // derivatives
                      for (int k = 0; k < 3; k++)
                        dx3Ddeta2D(i,j) += dX3DdXi3D(i,k) * dXi3Ddeta2D(k,j);

                  // get deformation factor
                  static LINALG::Matrix<2,2> Jac_tmp; // J^T*J
                  Jac_tmp.Clear();
                  Jac_tmp.MultiplyTN(dx3Ddeta2D,dx3Ddeta2D);

                  if(Jac_tmp.Determinant() == 0.0) dserror("deformation factor for boundary integration is zero");

                  const double deform_factor = sqrt(Jac_tmp.Determinant()); // sqrt(det(J^T*J))

                  const double fac = intpoints.qwgt[iquad]*deform_factor;

                  veljumpnormsquare += fac*veljumpnorm;
                }
#endif // computation average velocity jump
              }
            }
          }
        }
#ifdef COLLAPSE_FLAME
#ifdef GMSH_AVERAGE_JUMP
        cout << endl;
        // get the processor local node
        DRT::Node* lnode = discret_->lRowNode(0);

        LINALG::Matrix<3,1> nodecoord(true);
        // get physical coordinates of this node
        nodecoord(0) = lnode->X()[0];
        nodecoord(1) = lnode->X()[1];
        nodecoord(2) = lnode->X()[2];

        const double zcoord = nodecoord(2);
        const double deltaz = 2.*abs(zcoord);
        cout << "Netz " << 1./deltaz << endl;
        const double pi = atan(1.)*4.;
        const double area = pi*2.*0.25*deltaz;
        const double avveljump = sqrt(veljumpnormsquare/area);
        cout << "avveljump " << avveljump << endl;
#endif
#endif
        gmshfilecontent << "};\n";
      }
#endif
    }
    gmshfilecontent.close();
    if (screen_out) std::cout << " done" << endl;
  }
}

/*------------------------------------------------------------------------------------------------*
 | various options to initialize the fluid field                                      henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::SetInitialFlowField(
    const INPAR::COMBUST::InitialField initfield,
    const int initfuncno)
{
  //------------------------------------------
  // switch over different initial flow fields
  //------------------------------------------
  switch(initfield)
  {
  case INPAR::FLUID::initfield_zero_field:
  {
    // nothing to do
    break;
  }
  case INPAR::FLUID::initfield_beltrami_flow:
  {
    const Epetra_Map* dofrowmap = discret_->DofRowMap();

    int err =0;

    const int npredof = numdim_;

    double         p;
    vector<double> u  (numdim_);
    vector<double> xyz(numdim_);

    // check whether present flow is indeed three-dimensional
    if (numdim_!=3) dserror("Beltrami flow is a three-dimensional flow!");

    // set constants for analytical solution
    const double a      = M_PI/4.0;
    const double d      = M_PI/2.0;

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
        err += state_.velnp_->ReplaceMyValues(1,&(u[nveldof]),&lid);
        err += state_.veln_ ->ReplaceMyValues(1,&(u[nveldof]),&lid);
        err += state_.velnm_->ReplaceMyValues(1,&(u[nveldof]),&lid);
      }

      // set initial pressure
      const int gid = nodedofset[npredof];
      int lid = dofrowmap->LID(gid);
      err += state_.velnp_->ReplaceMyValues(1,&p,&lid);
      err += state_.veln_ ->ReplaceMyValues(1,&p,&lid);
      err += state_.velnm_->ReplaceMyValues(1,&p,&lid);

    } // end loop nodes lnodeid

    if(err!=0) dserror("dof not on proc");

    break;
  }
  case INPAR::COMBUST::initfield_field_by_function:
  case INPAR::COMBUST::initfield_disturbed_field_by_function:
  {
    const int numdim = params_->get<int>("number of velocity degrees of freedom");

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

        double initialval=DRT::Problem::Instance()->Funct(initfuncno-1).Evaluate(index,lnode->X(),time_,NULL);

        state_.velnp_->ReplaceGlobalValues(1,&initialval,&gid);
        state_.veln_ ->ReplaceGlobalValues(1,&initialval,&gid);
      }
    }

    // add random perturbation
    if(initfield == INPAR::COMBUST::initfield_disturbed_field_by_function)
    {
      const int numdim = params_->get<int>("number of velocity degrees of freedom");

      const Epetra_Map* dofrowmap = discret_->DofRowMap();

      int err =0;

      // random noise is perc percent of the initial profile

      double perc = params_->sublist("TURBULENCE MODEL").get<double>("CHAN_AMPL_INIT_DIST");

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

          thisvel=(*state_.velnp_)[lid];
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
//          map<int, vector<int> >::iterator master = pbcmapmastertoslave_->find(lnode->Id());

          // slavenodes are ignored
//          if(master == pbcmapmastertoslave_->end()) continue;
        }

        // add random noise on initial function field
        for(int index=0;index<numdim;++index)
        {
          int gid = nodedofset[index];

          double randomnumber = 2*((double)rand()-((double) RAND_MAX)/2.)/((double) RAND_MAX);

          double noise = perc * bmvel * randomnumber;

          err += state_.velnp_->SumIntoGlobalValues(1,&noise,&gid);
          err += state_.veln_ ->SumIntoGlobalValues(1,&noise,&gid);
        }

        if(err!=0)
        {
          dserror("dof not on proc");
        }
      }
    }
    break;
  }
  //----------------------------------------------------------------------------------------------
  // flame-vortex interaction problem: two counter-rotating vortices (2-D) moving the  flame front
  //----------------------------------------------------------------------------------------------
  case INPAR::COMBUST::initfield_flame_vortex_interaction:
  {
    // number space dimensions
    const int nsd = 3;
    // error indicator
    int err = 0;

    // define vectors for velocity field, node coordinates and coordinates of left and right vortices
    LINALG::Matrix<nsd,1> vel(true);
    double pres = 0.0;
    LINALG::Matrix<nsd,1> xyz(true);
    LINALG::Matrix<nsd,1> xyz0_left(true);
    LINALG::Matrix<nsd,1> xyz0_right(true);

    // set initial locations of vortices
    xyz0_left(0)  = 37.5;//87.5+0.78125; //37.5; // x-coordinate left vortex
    xyz0_left(1)  = 75.0; // y-coordinate left vortex
    xyz0_left(2)  = 0.0;  // z-coordinate is 0 (2D problem)
    xyz0_right(0) = 62.5;//12.5+0.78125; //62.5; // x-coordinate right vortex
    xyz0_right(1) = 75.0; // y-coordinate right vortex
    xyz0_right(2) = 0.0;  // z-coordinate is 0 (2D problem)

    // get laminar burning velocity (flame speed)
    if (flamespeed_ != 1.0) dserror("flame speed should be 1.0 for the 'flame-vortex-interaction' case");
    // vortex strength C (scaled by laminar burning velocity)
    const double C = 70.0*flamespeed_; // 70.0*flamespeed_;
    // (squared) vortex radius R
    const double R_squared = 16.0;

    //------------------------
    // get material parameters
    //------------------------
    // arbitrarily take first node on this proc
    DRT::Node* lnode = discret_->lRowNode(0);
    // get list of adjacent elements of the first node
    DRT::Element** elelist = lnode->Elements();
    // get material from first (arbitrary!) element adjacent to this node
    const Teuchos::RCP<MAT::Material> material = elelist[0]->Material();
#ifdef DEBUG
    // check if we really got a list of materials
    dsassert(material->MaterialType() == INPAR::MAT::m_matlist, "Material law is not of type m_matlist");
#endif
    // get material list for this element
    const MAT::MatList* matlist = static_cast<const MAT::MatList*>(material.get());

    // get burnt material (first material in material list)
    Teuchos::RCP<const MAT::Material> matptr0 = matlist->MaterialById(matlist->MatID(0));
    // get unburnt material (second material in material list)
    Teuchos::RCP<const MAT::Material> matptr1 = matlist->MaterialById(matlist->MatID(1));
#ifdef DEBUG
    dsassert(matptr0->MaterialType() == INPAR::MAT::m_fluid, "material is not of type m_fluid");
    dsassert(matptr1->MaterialType() == INPAR::MAT::m_fluid, "material is not of type m_fluid");
#endif
    const MAT::NewtonianFluid* mat0 = static_cast<const MAT::NewtonianFluid*>(matptr0.get());
    const MAT::NewtonianFluid* mat1 = static_cast<const MAT::NewtonianFluid*>(matptr1.get());

    // get the densities
    const double dens_b = mat0->Density();
    if (dens_b != 0.157) dserror("burnt density should be 0.157 for the 'flame-vortex-interaction' case");
    const double dens_u = mat1->Density();
    if (dens_u != 1.161) dserror("unburnt density should be 1.161 for the 'flame-vortex-interaction' case");
    double dens = dens_u;
    // for "pure fluid" computation: rhob = rhou = 1.161
    //const double dens_b = dens_u;

    // get map of global velocity vectors (DofRowMap)
    const Epetra_Map* dofrowmap = standarddofset_->DofRowMap();
    //const Epetra_Map* dofrowmap = discret_->DofRowMap();

//cout << *dofrowmap << endl;
//cout << *(standarddofset_->DofRowMap()) << endl;
//cout << (state_.velnp_->Map()) << endl;

    //--------------------------------
    // loop all nodes on the processor
    //--------------------------------
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node* lnode = discret_->lRowNode(lnodeid);

      // get node coordinates
      for(int idim=0;idim<nsd;idim++)
        xyz(idim)=lnode->X()[idim];

      // get phi value for this node
      const int lid = phinp_->Map().LID(lnode->Id());
      const double gfuncval = (*phinp_)[lid];

      // compute preliminary values for both vortices
      double r_squared_left  = ((xyz(0)-xyz0_left(0))*(xyz(0)-xyz0_left(0))
                               +(xyz(1)-xyz0_left(1))*(xyz(1)-xyz0_left(1)))/R_squared;
      double r_squared_right = ((xyz(0)-xyz0_right(0))*(xyz(0)-xyz0_right(0))
                               +(xyz(1)-xyz0_right(1))*(xyz(1)-xyz0_right(1)))/R_squared;

      //----------------------------------------
      // set density with respect to flame front
      //----------------------------------------
      if (gfuncval >= 0.0) // plus/burnt domain -> burnt material
      {
        dens = dens_b;
        pres = 0.0
        -0.5*(C*C/R_squared)*(exp(-r_squared_left) + exp(-r_squared_right));
      }
      else // minus/unburnt domain -> unburnt material
      {
        dens = dens_u;
        pres = flamespeed_*flamespeed_*dens_u*dens_u*(1.0/dens_u - 1.0/dens_b)
        -0.5*(C*C/R_squared)*(exp(-r_squared_left) + exp(-r_squared_right));
      }
      //----------------------------------------------
      // compute components of initial velocity vector
      //----------------------------------------------
      vel(0) = (C/R_squared)*(-(xyz(1)-xyz0_left(1))*exp(-r_squared_left/2.0)
                            +(xyz(1)-xyz0_right(1))*exp(-r_squared_right/2.0));
      vel(1) = (C/R_squared)*( (xyz(0)-xyz0_left(0))*exp(-r_squared_left/2.0)
                            -(xyz(0)-xyz0_right(0))*exp(-r_squared_right/2.0))
                            + flamespeed_*dens_u/dens;
      // 2D problem -> vel_z = 0.0
      vel(2) = 0.0;
      // velocity profile without vortices
      //vel(1) = sl*densu/dens;

      // access standard FEM dofset (3 x vel + 1 x pressure) to get dof IDs for this node
      const vector<int> nodedofs = (*standarddofset_).Dof(lnode);
      //const vector<int> nodedofs = discret_->Dof(lnode);
      //for (int i=0;i<standardnodedofset.size();i++)
      //{
      //  cout << "component " << i << " standarddofset dofid " << stdnodedofset[i] << endl;
      //}

      //-----------------------------------------
      // set components of initial velocity field
      //-----------------------------------------
      for(int idim=0;idim<nsd+1;idim++)
      {
        const int gid = nodedofs[idim];
        //local node id
        int lid = dofrowmap->LID(gid);
        if(idim==3){ // pressure dof
          err += state_.velnp_->ReplaceMyValues(1,&pres,&lid);
          err += state_.veln_ ->ReplaceMyValues(1,&pres,&lid);
          err += state_.velnm_->ReplaceMyValues(1,&pres,&lid);
        }
        else{ // velocity dof
          err += state_.velnp_->ReplaceMyValues(1,&vel(idim),&lid);
          err += state_.veln_ ->ReplaceMyValues(1,&vel(idim),&lid);
          err += state_.velnm_->ReplaceMyValues(1,&vel(idim),&lid);
        }
      }
    } // end loop nodes lnodeid

    if (err!=0) dserror("dof not on proc");

    break;
  }
  default:
    dserror("type of initial field not available");
  }

  return;
}


void FLD::CombustFluidImplicitTimeInt::SetEnrichmentField(
    const Teuchos::RCP<XFEM::DofManager> dofmanager,
    const Epetra_Map dofrowmap)
{
#ifdef FLAME_VORTEX
  cout0_ << "---  set initial enrichment field for flame-vortex interaction example... " << std::flush;

  // initial field modification for flame_vortex_interaction
  for (int nodeid=0;nodeid<discret_->NumMyRowNodes();nodeid++) // loop over element nodes
  {
    DRT::Node* lnode = discret_->lRowNode(nodeid);
    const std::set<XFEM::FieldEnr>& fieldenrset(dofmanager->getNodeDofSet(lnode->Id()));
    for (set<XFEM::FieldEnr>::const_iterator fieldenr = fieldenrset.begin();
        fieldenr != fieldenrset.end();++fieldenr)
    {
      const XFEM::DofKey dofkey(lnode->Id(), *fieldenr);
      const int dofpos = state_.nodalDofDistributionMap_.find(dofkey)->second;
      if (fieldenr->getEnrichment().Type() == XFEM::Enrichment::typeJump)
      {
        if (fieldenr->getField() == XFEM::PHYSICS::Velx);
          // nothing to do in flame_vortex example
        else if (fieldenr->getField() == XFEM::PHYSICS::Vely)
        {
          (*state_.veln_)[dofrowmap.LID(dofpos)] = 3.197452229299363;
          (*state_.velnp_)[dofrowmap.LID(dofpos)] = 3.197452229299363;
        }
        else if (fieldenr->getField() == XFEM::PHYSICS::Velz);
          // nothing to do in flame_vortex example
        else if (fieldenr->getField() == XFEM::PHYSICS::Pres)
        {
          (*state_.veln_)[dofrowmap.LID(dofpos)] = 3.7122420382165605;
          (*state_.velnp_)[dofrowmap.LID(dofpos)] = 3.7122420382165605;
        }
      } // end if jump enrichment
    } // end loop over fieldenr
  } // end loop over element nodes
  cout0_ << "done" << std::endl;
#endif

#ifdef COLLAPSE_FLAME
  cout0_ << "---  set initial enrichment field for collapsing flame example... " << std::flush;

  // initial field modification for collapse_flame
  const int nsd = 3;
  for (int nodeid=0;nodeid<discret_->NumMyRowNodes();nodeid++) // loop over element nodes
  {
    DRT::Node* lnode = discret_->lRowNode(nodeid);

    LINALG::Matrix<nsd,1> coords(lnode->X());
#ifdef COMBUST_2D
    coords(2)=0.0;
#endif
    double coordsnorm = sqrt(coords(0)*coords(0)+coords(1)*coords(1)+coords(2)*coords(2));

    const int lid = phinp_->Map().LID(lnode->Id());
    const double gfuncval = (*phinp_)[lid];

    const double radius = 0.025;
    const double velrad = 1.0;
    const double velgrad = 40.0;
    const double densu = 1.0;

    const std::set<XFEM::FieldEnr>& fieldenrset(dofmanager->getNodeDofSet(lnode->Id()));
    for (set<XFEM::FieldEnr>::const_iterator fieldenr = fieldenrset.begin();
        fieldenr != fieldenrset.end();++fieldenr)
    {
      const XFEM::DofKey dofkey(lnode->Id(), *fieldenr);
      const int dofpos = state_.nodalDofDistributionMap_.find(dofkey)->second;
      if (fieldenr->getEnrichment().Type() == XFEM::Enrichment::typeJump)
      {
#ifndef COMBUST_NORMAL_ENRICHMENT
        if (fieldenr->getField() == XFEM::PHYSICS::Velx)
        {
          (*state_.veln_)[dofrowmap.LID(dofpos)] = (-0.5*velrad + 0.5*gfuncval*velgrad)*coords(0)/coordsnorm;
          (*state_.velnp_)[dofrowmap.LID(dofpos)] = (-0.5*velrad + 0.5*gfuncval*velgrad)*coords(0)/coordsnorm;

        }
        else if (fieldenr->getField() == XFEM::PHYSICS::Vely)
        {
          (*state_.veln_)[dofrowmap.LID(dofpos)] = (-0.5*velrad + 0.5*gfuncval*velgrad)*coords(1)/coordsnorm;
          (*state_.velnp_)[dofrowmap.LID(dofpos)] = (-0.5*velrad + 0.5*gfuncval*velgrad)*coords(1)/coordsnorm;
        }
        else if (fieldenr->getField() == XFEM::PHYSICS::Velz)
        {
          (*state_.veln_)[dofrowmap.LID(dofpos)] = (-0.5*velrad + 0.5*gfuncval*velgrad)*coords(2)/coordsnorm;
          (*state_.velnp_)[dofrowmap.LID(dofpos)] = (-0.5*velrad + 0.5*gfuncval*velgrad)*coords(2)/coordsnorm;
#ifdef COMBUST_2D
          (*state_.veln_)[dofrowmap.LID(dofpos)] = 0.0;
          (*state_.velnp_)[dofrowmap.LID(dofpos)] = 0.0;
#endif
        }
#else
        if (fieldenr->getField() == XFEM::PHYSICS::Veln)
        {
          // -0.5 *jump + 0.5*dist*kink
          (*state_.veln_)[dofrowmap.LID(dofpos)] = -0.5*velrad + 0.5*gfuncval*velgrad;
          (*state_.velnp_)[dofrowmap.LID(dofpos)] = -0.5*velrad + 0.5*gfuncval*velgrad;
        }
#endif
        else if (fieldenr->getField() == XFEM::PHYSICS::Pres)
        {
          // -0.5 *jump + 0.5*dist*kink
          //(*state_.veln_)[dofrowmap.LID(dofpos)] = -0.5*(0.92) - 0.5*gfuncval*40.0; // 0.5*(7.0)
          //(*state_.velnp_)[dofrowmap.LID(dofpos)] = -0.5*(0.92) - 0.5*gfuncval*40.0; // 0.5*(7.0)
        }
      } // end if jump enrichment
      else if (fieldenr->getEnrichment().Type() == XFEM::Enrichment::typeStandard)
      {
        if (XFEM::plusDomain(gfuncval) == true) // Standard dofs innerhalb des Kreises
        {
          if (fieldenr->getField() == XFEM::PHYSICS::Velx or
              fieldenr->getField() == XFEM::PHYSICS::Vely or
              fieldenr->getField() == XFEM::PHYSICS::Velz )
          {
            (*state_.veln_)[dofrowmap.LID(dofpos)] = 0.0;
            (*state_.velnp_)[dofrowmap.LID(dofpos)] = 0.0;
          }
          else if (fieldenr->getField() == XFEM::PHYSICS::Pres)
          {
            //(*state_.veln_)[dofrowmap.LID(dofpos)] = -1.42; //-6.5;
            //(*state_.velnp_)[dofrowmap.LID(dofpos)] = -1.42; //-6.5;
          }
        }
        else // Standard dofs ausshalb des Kreises
        {
          if (fieldenr->getField() == XFEM::PHYSICS::Velx)
          {
            (*state_.veln_)[dofrowmap.LID(dofpos)] = radius*radius*velrad*coords(0)/(coordsnorm*coordsnorm*coordsnorm);
            (*state_.velnp_)[dofrowmap.LID(dofpos)] = radius*radius*velrad*coords(0)/(coordsnorm*coordsnorm*coordsnorm);
#ifdef COMBUST_2D
            (*state_.veln_)[dofrowmap.LID(dofpos)] = radius*velrad*coords(0)/(coordsnorm*coordsnorm);
            (*state_.velnp_)[dofrowmap.LID(dofpos)] = radius*velrad*coords(0)/(coordsnorm*coordsnorm);
#endif
          }
          else if (fieldenr->getField() == XFEM::PHYSICS::Vely)
          {
            (*state_.veln_)[dofrowmap.LID(dofpos)] = radius*radius*velrad*coords(1)/(coordsnorm*coordsnorm*coordsnorm);
            (*state_.velnp_)[dofrowmap.LID(dofpos)] = radius*radius*velrad*coords(1)/(coordsnorm*coordsnorm*coordsnorm);
#ifdef COMBUST_2D
            (*state_.veln_)[dofrowmap.LID(dofpos)] = radius*velrad*coords(1)/(coordsnorm*coordsnorm);
            (*state_.velnp_)[dofrowmap.LID(dofpos)] = radius*velrad*coords(1)/(coordsnorm*coordsnorm);
#endif
          }
          else if (fieldenr->getField() == XFEM::PHYSICS::Velz)
          {
            (*state_.veln_)[dofrowmap.LID(dofpos)] = radius*radius*velrad*coords(2)/(coordsnorm*coordsnorm*coordsnorm);
            (*state_.velnp_)[dofrowmap.LID(dofpos)] = radius*radius*velrad*coords(2)/(coordsnorm*coordsnorm*coordsnorm);
#ifdef COMBUST_2D
            (*state_.veln_)[dofrowmap.LID(dofpos)] = 0.0;
            (*state_.velnp_)[dofrowmap.LID(dofpos)] = 0.0;
#endif
          }
          else if (fieldenr->getField() == XFEM::PHYSICS::Pres)
          {
            //(*state_.veln_)[dofrowmap.LID(dofpos)] = 0.5*densu*radius*velrad*radius*velrad/(coordsnorm*coordsnorm);
            //(*state_.velnp_)[dofrowmap.LID(dofpos)] = 0.5*densu*radius*velrad*radius*velrad/(coordsnorm*coordsnorm);
          }
        }
      }
    } // end loop over fieldenr
  } // end loop over element nodes
  cout0_ << "done" << std::endl;
#endif

#ifdef COMBUST_TWO_FLAME_FRONTS
  cout0_ << "---  set initial enrichment field for two approaching flame fronts example... " << std::flush;

  // initial field modification for 2-flames example
  for (int nodeid=0;nodeid<discret_->NumMyRowNodes();nodeid++) // loop over element nodes
  {
    DRT::Node* lnode = discret_->lRowNode(nodeid);

    LINALG::Matrix<3,1> coords(lnode->X());

    const int lid = phinp_->Map().LID(lnode->Id());
    const double gfuncval = (*phinp_)[lid];

    const std::set<XFEM::FieldEnr>& fieldenrset(dofmanager->getNodeDofSet(lnode->Id()));
    for (set<XFEM::FieldEnr>::const_iterator fieldenr = fieldenrset.begin();
        fieldenr != fieldenrset.end();++fieldenr)
    {
      const XFEM::DofKey dofkey(lnode->Id(), *fieldenr);
      const int dofpos = state_.nodalDofDistributionMap_.find(dofkey)->second;
      if (fieldenr->getEnrichment().Type() == XFEM::Enrichment::typeJump)
      {
        if (fieldenr->getField() == XFEM::PHYSICS::Velx)
        {
          if (coords(0) > 0.0) // right flame
          {
            (*state_.veln_)[dofrowmap.LID(dofpos)] = 0.5;
            (*state_.velnp_)[dofrowmap.LID(dofpos)] = 0.5;
          }
          else // left flame
          {
            (*state_.veln_)[dofrowmap.LID(dofpos)] = -0.5;
            (*state_.velnp_)[dofrowmap.LID(dofpos)] = -0.5;
          }
        }
        else if (fieldenr->getField() == XFEM::PHYSICS::Vely or
                 fieldenr->getField() == XFEM::PHYSICS::Velz)
        {
          (*state_.veln_)[dofrowmap.LID(dofpos)] = 0;
          (*state_.velnp_)[dofrowmap.LID(dofpos)] = 0;
        }
        else if (fieldenr->getField() == XFEM::PHYSICS::Pres)
        {
          (*state_.veln_)[dofrowmap.LID(dofpos)] = 1.0;
          (*state_.velnp_)[dofrowmap.LID(dofpos)] = 1.0;
        }
      } // end if jump enrichment
      else if (fieldenr->getEnrichment().Type() == XFEM::Enrichment::typeStandard)
      {
        if (gfuncval>=0)
        {
          if (fieldenr->getField() == XFEM::PHYSICS::Velx)
          {
            if (coords(0) > 0)
            {
              (*state_.veln_)[dofrowmap.LID(dofpos)] = 1.0;
              (*state_.velnp_)[dofrowmap.LID(dofpos)] = 1.0;
            }
            else
            {
              (*state_.veln_)[dofrowmap.LID(dofpos)] = -1.0;
              (*state_.velnp_)[dofrowmap.LID(dofpos)] = -1.0;
            }
          }
          else if (fieldenr->getField() == XFEM::PHYSICS::Vely or
                   fieldenr->getField() == XFEM::PHYSICS::Velz)
          {
            (*state_.veln_)[dofrowmap.LID(dofpos)] = 0.0;
            (*state_.velnp_)[dofrowmap.LID(dofpos)] = 0.0;
          }
          else if (fieldenr->getField() == XFEM::PHYSICS::Pres)
          {
            (*state_.veln_)[dofrowmap.LID(dofpos)] = 2.0;
            (*state_.velnp_)[dofrowmap.LID(dofpos)] = 2.0;
          }
        }
        else
        {
          (*state_.veln_)[dofrowmap.LID(dofpos)] = 0.0;
          (*state_.velnp_)[dofrowmap.LID(dofpos)] = 0.0;
        }
      }
    } // end loop over fieldenr
  } // end loop over element nodes
  cout0_ << "done" << std::endl;
#endif

//#ifdef GMSH_REF_FIELDS
  OutputToGmsh((char*)"mod_start_field_pres",(char*)"mod_start_field_vel",Step(), Time());
//#endif
}


// -------------------------------------------------------------------
// -------------------------------------------------------------------
void FLD::CombustFluidImplicitTimeInt::SetupXFluidSplit(
        const DRT::Discretization& dis,
        const RCP<XFEM::DofManager> dofman,
        LINALG::MapExtractor& extractor)
{
  // -------------------------------------------------------------------
  // get a vector layout from the discretization for a vector which only
  // contains the velocity dofs and for one vector which only contains
  // pressure degrees of freedom.
  //
  // The maps are designed assuming that every node has pressure and
  // velocity degrees of freedom --- this won't work for inf-sup stable
  // elements at the moment!
  // -------------------------------------------------------------------

  // Allocate integer vectors which will hold the dof number of the
  // velocity or pressure dofs
  vector<int> velmapdata;
  vector<int> premapdata;

  // collect global dofids for velocity and pressure in vectors
  for (int i=0; i<dis.NumMyRowNodes(); ++i) {
    const DRT::Node* node = dis.lRowNode(i);
    const std::set<XFEM::FieldEnr>& enrvarset(dofman->getNodeDofSet(node->Id()));
    const vector<int> dof = dis.Dof(node);
    dsassert(dof.size() == enrvarset.size(), "mismatch in length!");
    std::set<XFEM::FieldEnr>::const_iterator enrvar;
    size_t countdof = 0;
    for (enrvar = enrvarset.begin(); enrvar != enrvarset.end(); ++enrvar)
    {
      switch (enrvar->getField()) {
      case XFEM::PHYSICS::Velx:
      case XFEM::PHYSICS::Vely:
      case XFEM::PHYSICS::Velz:
        velmapdata.push_back(dof[countdof]);
        break;
      case XFEM::PHYSICS::Pres:
        premapdata.push_back(dof[countdof]);
        break;
      default:
        break;
      }
      countdof++;
    }
  }

  // the rowmaps are generated according to the pattern provided by
  // the data vectors
  RCP<Epetra_Map> velrowmap = rcp(new Epetra_Map(-1,
      velmapdata.size(),&velmapdata[0],0,
      dis.Comm()));
  RCP<Epetra_Map> prerowmap = rcp(new Epetra_Map(-1,
      premapdata.size(),&premapdata[0],0,
      dis.Comm()));

  const Epetra_Map* map = dis.DofRowMap();
  extractor.Setup(*map, prerowmap, velrowmap);
}


/*--------------------------------------------------------------------------------------------*
 | Redistribute the fluid discretization and vectors according to nodegraph   rasthofer 07/11 |
 |                                                                            DA wichmann     |
 *--------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::Redistribute(const Teuchos::RCP<Epetra_CrsGraph> nodegraph)
{
  // the fluid dis must be set to standard mode in order to avoid trouble with the dofmanager
  ParameterList eleparams;
  eleparams.set("action","set_standard_mode");
  discret_->Evaluate(eleparams);

  // the rowmap will become the new distribution of nodes
  const Epetra_BlockMap rntmp = nodegraph->RowMap();
  Epetra_Map newnoderowmap(-1,rntmp.NumMyElements(),rntmp.MyGlobalElements(),0,discret_->Comm());

  // the column map will become the new ghosted distribution of nodes
  const Epetra_BlockMap Mcntmp = nodegraph->ColMap();
  Epetra_Map newnodecolmap(-1,Mcntmp.NumMyElements(),Mcntmp.MyGlobalElements(),0,discret_->Comm());

  // do the redistribution
  discret_->Redistribute(newnoderowmap,newnodecolmap, false, false, false);

  // assign the new dofs, make absolutely sure that we always
  // have all slaves to a master
  // the finite edge weights are not a 100% warranty for that...
  // update the PBCs and PBCDofSet
  // includes call to FillComplete()
  pbc_->PutAllSlavesToMastersProc();
  pbcmapmastertoslave_ = pbc_->ReturnAllCoupledColNodes();

  // create dummy instance of interfacehandle holding no flamefront and hence no integration cells
  Teuchos::RCP<COMBUST::InterfaceHandleCombust> ihdummy = rcp(new COMBUST::InterfaceHandleCombust(discret_,Teuchos::null));
  // create dummy instance of dof manager assigning standard enrichments to all nodes
  const Teuchos::RCP<XFEM::DofManager> dofmanagerdummy = rcp(new XFEM::DofManager(ihdummy,Teuchos::null,physprob_.xfemfieldset_,*physprob_.elementAnsatz_,xparams_,Teuchos::null));

  // save dofmanager to be able to plot Gmsh stuff in Output()
  dofmanagerForOutput_ = dofmanagerdummy;

  // pass dof information to elements (no enrichments yet, standard FEM!)
  TransferDofInformationToElements(Teuchos::null, dofmanagerdummy);

  // ensure that degrees of freedom in the discretization have been set
  if ((not discret_->Filled()) or (not discret_->HaveDofs()))
    discret_->FillComplete();

  // update the dofset with the complete fluid unknowns
  standarddofset_->SetCoupledNodes(pbcmapmastertoslave_);
  standarddofset_->Reset();
  standarddofset_->AssignDegreesOfFreedom(*discret_,0,0);

  // split based on complete fluid field
  FLD::UTILS::SetupFluidSplit(*discret_,*standarddofset_,3,*velpressplitterForOutput_);

  return;
}


void FLD::CombustFluidImplicitTimeInt::TransferVectorsToNewDistribution(
    Teuchos::RCP<COMBUST::FlameFront> flamefront)
{
  // build instance of DofManager with information about the interface from the interfacehandle
  // remark: DofManager is rebuilt in every inter-field iteration step, because number and position
  // of enriched degrees of freedom change in every time step/FG-iteration
  const Teuchos::RCP<XFEM::DofManager> dofmanager = Teuchos::rcp(new XFEM::DofManager(
      flamefront->InterfaceHandle(),
      flamefront->Phinp(),
      physprob_.xfemfieldset_,
      *physprob_.elementAnsatz_,
      xparams_,
      pbcmapmastertoslave_)
  );

  // save dofmanager to be able to plot Gmsh stuff in Output()
  dofmanagerForOutput_ = dofmanager;

  // -------------------------------------------------------------------
  // connect degrees of freedom for periodic boundary conditions
  // -------------------------------------------------------------------
  // tell elements about the dofs and the integration
  TransferDofInformationToElements(flamefront, dofmanager);
  // assign degrees of freedom
  // remark: - assign degrees of freedom (first slot)
  //         - build geometry for (Neumann) boundary conditions (third slot);
  //           without Neumann boundary conditions Fillcomplete(true,false,false) will also work
  discret_->FillComplete(true,false,true);

  dofmanager->fillDofRowDistributionMaps(
        state_.nodalDofDistributionMap_);

  //------------------------------------------------------------------------------------------------
  // get dof layout from the discretization to construct vectors and matrices
  //------------------------------------------------------------------------------------------------

  // parallel dof distribution contained in dofrowmap: local (LID) <-> global (GID) dof numbering
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  //------------------------------------------------------------------------------------------------
  // rewrite vectors according to the new distribution
  //------------------------------------------------------------------------------------------------

  Teuchos::RCP<Epetra_Vector> old;

  if (state_.velnp_ != Teuchos::null)
  {
    old = state_.velnp_;
    state_.velnp_ = LINALG::CreateVector(*dofrowmap,true);
    LINALG::Export(*old, *state_.velnp_);
  }

  if (state_.veln_ != Teuchos::null)
  {
    old = state_.veln_;
    state_.veln_ = LINALG::CreateVector(*dofrowmap,true);
    LINALG::Export(*old, *state_.veln_);
  }

  if (state_.velnm_ != Teuchos::null)
  {
    old = state_.velnm_;
    state_.velnm_ = LINALG::CreateVector(*dofrowmap,true);
    LINALG::Export(*old, *state_.velnm_);
  }

  if (state_.velaf_ != Teuchos::null)
  {
    old = state_.velaf_;
    state_.velaf_ = LINALG::CreateVector(*dofrowmap,true);
    LINALG::Export(*old, *state_.velaf_);
  }

  // acceleration at time n+1 and n
  if (state_.accnp_ != Teuchos::null)
  {
    old = state_.accnp_;
    state_.accnp_ = LINALG::CreateVector(*dofrowmap,true);
    LINALG::Export(*old, *state_.accnp_);
  }

  if (state_.accn_ != Teuchos::null)
  {
    old = state_.accn_;
    state_.accn_ = LINALG::CreateVector(*dofrowmap,true);
    LINALG::Export(*old, *state_.accn_);
  }

  if (state_.accam_ != Teuchos::null)
  {
    old = state_.accam_;
    state_.accam_ = LINALG::CreateVector(*dofrowmap,true);
    LINALG::Export(*old, *state_.accam_);
  }

  // -------------------------------------------------------------------
  // initialize turbulence-statistics evaluation
  // -------------------------------------------------------------------
  // currently omitted

  return;
} // FLD::CombustFluidImplicitTimeInt::Redistribute




/*------------------------------------------------------------------------------------------------*
 | create field test
 *------------------------------------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> FLD::CombustFluidImplicitTimeInt::CreateFieldTest()
{
  return Teuchos::rcp(new FLD::CombustFluidResultTest(*this));
}


/*------------------------------------------------------------------------------------------------*
 | compute mesh dependent error norms for Nitsche's method                           schott 05/10 |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::EvaluateErrorComparedToAnalyticalSol_Nitsche(INPAR::COMBUST::NitscheError& NitscheErrorType)
{
  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("velnp",state_.velnp_);
  discret_->SetState("veln" ,state_.veln_);
  discret_->SetState("velnm",state_.velnm_);
  discret_->SetState("accn" ,state_.accn_);
  discret_->SetState("nodal increment",oldinc_);

  // create parameters for discretization
  ParameterList eleparams;

  eleparams.set("action", "calc_nitsche_error");
  eleparams.set<int>("Nitsche_Compare_Analyt", NitscheErrorType);

  // smoothed normal vectors for boundary integration terms
  eleparams.set("smoothed_bound_integration", smoothed_boundary_integration_);
  eleparams.set<int>("smoothgradphi",smoothgradphi_);
  // flag for type of combustion problem
  eleparams.set<int>("combusttype",combusttype_);

  eleparams.set("flamespeed",flamespeed_);
  eleparams.set("time",time_);

  // set parameters for parts of the whole Nitsche-error (mesh-dependent norms), here norms not square rooted
  eleparams.set<double>("L2 integrated velocity domain error", 0.0);
  eleparams.set<double>("L2 integrated grad_velocity domain error", 0.0);
  eleparams.set<double>("H-1/2 integrated viscosity interface error", 0.0);
  eleparams.set<double>("H1/2 integrated velocity jump interface error", 0.0);
  eleparams.set<double>("H-1/2 integrated flux jump interface error",0.0);
  eleparams.set<double>("L2 integrated pressure domain error", 0.0);
  eleparams.set<double>("L2 integrated grad_pressure domain error", 0.0);
  eleparams.set<double>("L2 integrated weighted pressure domain error", 0.0);
  eleparams.set<double>("Nitsche integrated error", 0.0); // the whole Nitsche error
  eleparams.set<double>("L2 Divergence integrated error", 0.0);
  eleparams.set<double>("L2 Divergence error in omega+", 0.0);
  eleparams.set<double>("L2 Divergence error in omega-", 0.0);
  eleparams.set<double>("L2 integrated velocity jump interface error Combustion", 0.0);
  eleparams.set<double>("L2 integrated flux jump interface error Combustion", 0.0);

  // call loop over elements (but do not assemble anything)
  // discret_->Evaluate calls combust3_evaluate for each element, but without assemply
  discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);

  discret_->ClearState();

  double locVelDomErr         = eleparams.get<double>("L2 integrated velocity domain error");
  double locGradVelDomErr     = eleparams.get<double>("L2 integrated grad_velocity domain error");
  double locViscInterfErr     = eleparams.get<double>("H-1/2 integrated viscosity interface error");
  double locVelJumpInterfErr  = eleparams.get<double>("H1/2 integrated velocity jump interface error");
  double locFluxJumpInterfErr = eleparams.get<double>("H-1/2 integrated flux jump interface error");
  double locPresDomErr        = eleparams.get<double>("L2 integrated pressure domain error");
  double locGradPresDomErr    = eleparams.get<double>("L2 integrated grad_pressure domain error");
  double locWeightPresDomErr  = eleparams.get<double>("L2 integrated weighted pressure domain error");
  double locNitscheErr        = eleparams.get<double>("Nitsche integrated error");
  //double locDivErr           = eleparams.get<double>("L2 Divergence integrated error");
  //double locDivErrPlus       = eleparams.get<double>("L2 Divergence error in omega+");
  //double locDivErrMinus      = eleparams.get<double>("L2 Divergence error in omega-");
  //double locVelJumpErr       = eleparams.get<double>("L2 integrated velocity jump interface error Combustion");

  // initialize global errors
  double VelDomErr         = 0.0;
  double GradVelDomErr     = 0.0;
  double ViscInterfErr     = 0.0;
  double VelJumpInterfErr  = 0.0;
  double FluxJumpInterfErr = 0.0;
  double PresDomErr        = 0.0;
  double GradPresDomErr    = 0.0;
  double WeightPresDomErr  = 0.0;
  double NitscheErr        = 0.0;
  //double DivErr           = 0.0;
  //double DivErrPlus       = 0.0;
  //double DivErrMinus      = 0.0;
  //double VelJumpErr       = 0.0;

  // sum over processors, each list (list of a processor) has length 1
  discret_->Comm().SumAll(&locVelDomErr,        &VelDomErr,        1);
  discret_->Comm().SumAll(&locGradVelDomErr,    &GradVelDomErr,    1);
  discret_->Comm().SumAll(&locViscInterfErr,    &ViscInterfErr,    1);
  discret_->Comm().SumAll(&locVelJumpInterfErr, &VelJumpInterfErr, 1);
  discret_->Comm().SumAll(&locFluxJumpInterfErr,&FluxJumpInterfErr,1);
  discret_->Comm().SumAll(&locPresDomErr,       &PresDomErr,       1);
  discret_->Comm().SumAll(&locGradPresDomErr,   &GradPresDomErr,   1);
  discret_->Comm().SumAll(&locWeightPresDomErr, &WeightPresDomErr, 1);
  discret_->Comm().SumAll(&locNitscheErr,       &NitscheErr,       1);
  //discret_->Comm().SumAll(&locDivErr,          &DivErr,           1);
  //discret_->Comm().SumAll(&locDivErrPlus,      &DivErrPlus,       1);
  //discret_->Comm().SumAll(&locDivErrMinus,     &DivErrMinus,      1);
  //discret_->Comm().SumAll(&locVelJumpErr,      &VelJumpErr,       1);

  // for the norms, we need the square roots
  VelDomErr         = sqrt(VelDomErr);
  GradVelDomErr     = sqrt(GradVelDomErr);
  ViscInterfErr     = sqrt(ViscInterfErr);
  VelJumpInterfErr  = sqrt(VelJumpInterfErr);
  FluxJumpInterfErr = sqrt(FluxJumpInterfErr);
  PresDomErr        = sqrt(PresDomErr);
  GradPresDomErr    = sqrt(GradPresDomErr);
  WeightPresDomErr  = sqrt(WeightPresDomErr);
  NitscheErr        = sqrt(NitscheErr);
  //DivErr            = sqrt(DivErr);
  //DivErrPlus        = sqrt(DivErrPlus);
  //DivErrMinus       = sqrt(DivErrMinus);
  //VelJumpErr        = sqrt(VelJumpErr);

  if (myrank_ == 0)
  {
    printf("\n======================================================================="
           "\n======================= absolute Nitsche errors ======================="
           "\n======= compare analytical solution with approximated solution=========");
    printf("\n  || u-u_h ||_L2(Omega)\t\t\t\t\t%15.8e",                  VelDomErr);
    printf("\n  || sqrt(mu)grad(u-u_h) ||_L2(Omega) \t\t\t%15.8e",       GradVelDomErr);
    printf("\n  || p-p_h ||_L2(Omega)\t\t\t\t\t%15.8e",                  PresDomErr);
    printf("\n  || grad(p-p_h) ||_L2(Omega)\t\t\t\t%15.8e",              GradPresDomErr);
    printf("\n  || {2mu* E(u-u_h)*n} ||_H-1/2(Gamma)\t\t\t%15.8e",       ViscInterfErr);
    printf("\n  || [[u]]-[[u_h]] ||_H1/2(Gamma)\t\t\t%15.8e",            VelJumpInterfErr);
    //printf("\n  || 1/sqrt(mu_max) * (p-p_h) ||_L2(Omega)\t\t%15.8e",     WeightPresDomErr);
    //printf("\n  || div(u) ||_L2(Omega)\t\t\t\t%15.8e",                   DivErr);
    //printf("\n  || div(u) ||_L2(Omega+)\t\t\t\t%15.8e",                  DivErrPlus);
    //printf("\n  || div(u) ||_L2(Omega-)\t\t\t\t%15.8e",                  DivErrMinus);
    //printf("\n  || [| u |] - ju*n ||_L2(Gamma)\t\t\t%15.8e",             VelJumpErr);
    printf("\n  || [[sigma*n]]-[[jflux*n]] ||_H-1/2(Gamma)\t\t%15.8e",   FluxJumpInterfErr);
    printf("\n ||| (u-u_h, p-p_h) |||_Nitsche(Omega)\t\t\t%15.8e",       NitscheErr);
    printf("\n======================================================================="
           "\n=======================================================================\n");
  }

  return;
}



/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::EvaluateErrorComparedToAnalyticalSol()
{
  INPAR::FLUID::InitialField calcerr = DRT::INPUT::get<INPAR::FLUID::InitialField>(*params_, "eval err for analyt sol");

  //------------------------------------------------------- beltrami flow
  switch (calcerr)
  {
  case INPAR::FLUID::initfield_zero_field:
  case INPAR::FLUID::initfield_field_by_function:
  case INPAR::FLUID::initfield_disturbed_field_from_function:
  case INPAR::FLUID::initfield_flame_vortex_interaction:
    // do nothing --- no analytical solution available
    break;
  case INPAR::FLUID::initfield_beltrami_flow:
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
    discret_->SetState("u and p at time n+1 (converged)",state_.velnp_);

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
      printf("\n  L2_err for beltrami flow:  velocity %15.8e  pressure %15.8e\n\n", velerr,preerr);
    }
  }
  break;
  default:
    dserror("Cannot calculate error. Unknown type of analytical test problem");
  }
  return;
}

/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::SolveStationaryProblem()
{

/*
 * This function is commented out because it should not be used by any combustion simulation.
 * The algorithm to be called is the one in the class COMBUST::Algorithm. It accesses directly e.g.
 * NonlinearSolve in this class CombustFluidImplicitTimeInt.
 */

  dserror("This is the wrong stationary algorithm! Use COMBUST::Algorithm::SolveStationaryProblem()");

}


Teuchos::RCP<const DRT::DofSet> FLD::CombustFluidImplicitTimeInt::DofSet()
{
  return standarddofset_;
}


Teuchos::RCP<const Epetra_Map> FLD::CombustFluidImplicitTimeInt::VelocityRowMap()
{
  return velpressplitter_->OtherMap();
}


Teuchos::RCP<const Epetra_Map> FLD::CombustFluidImplicitTimeInt::PressureRowMap()
{
  return velpressplitter_->CondMap();
}


/// return time integration factor
double FLD::CombustFluidImplicitTimeInt::TimIntParam() const
{
  double retval = 1.0;
  switch (TimIntScheme())
  {
  case INPAR::FLUID::timeint_afgenalpha:
  case INPAR::FLUID::timeint_gen_alpha:
  case INPAR::FLUID::timeint_npgenalpha:
    retval = alphaF_;
  break;
  case INPAR::FLUID::timeint_one_step_theta:
    // this is the point where OST is evaluated
    retval = 1.0;
  break;
  case INPAR::FLUID::timeint_bdf2:
    // this is the point where bdf2 is evaluated
    retval = 1.0;
  break;
  case INPAR::FLUID::timeint_stationary:
    // this is the point where stat. is evaluated
    retval = 1.0;
  break;
  default:
    dserror("Unknown time integration scheme");
  break;
  }
  return retval;
}


/*------------------------------------------------------------------------------------------------*
 |                                                                                 rasthofer 07/12|
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::LiftDrag() const
{
  // in this map, the results of the lift drag calculation are stored
  RCP<map<int,vector<double> > > liftdragvals;

  FLD::UTILS::LiftDrag(*discret_,*trueresidual_,*params_,liftdragvals);

  if (liftdragvals!=Teuchos::null and discret_->Comm().MyPID() == 0)
    FLD::UTILS::WriteLiftDragToFile(time_, step_, *liftdragvals);

  return;
} // CombustFluidImplicitTimeInt::LiftDrag

///*------------------------------------------------------------------------------------------------*
// | henke 08/08 |
// *------------------------------------------------------------------------------------------------*/
//Teuchos::RCP<DRT::ResultTest> FLD::CombustFluidImplicitTimeInt::CreateFieldTest()
//{
//  return Teuchos::rcp(new FLD::CombustFluidResultTest(fluid_));
//}



/*------------------------------------------------------------------------------------------------*
| returns matching string for each time integration scheme                              gjb 08/08 |
*-------------------------------------------------------------------------------------------------*/
std::string FLD::CombustFluidImplicitTimeInt::MapTimIntEnumToString(const enum INPAR::FLUID::TimeIntegrationScheme term)
{
  // length of return string is 14 due to usage in formated screen output
  switch (term)
  {
  case INPAR::FLUID::timeint_one_step_theta:
    return "One-Step-Theta";
    break;
  case INPAR::FLUID::timeint_bdf2 :
    return "    BDF2      ";
    break;
  case INPAR::FLUID::timeint_stationary :
    return "  Stationary  ";
    break;
  case INPAR::FLUID::timeint_afgenalpha :
    return "  Gen. Alpha  ";
    break;
  default :
    dserror("Cannot cope with name %d", term);
    return "";
    break;
  }
}

