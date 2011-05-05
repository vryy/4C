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

#undef WRITEOUTSTATISTICS

#include "fluidimplicitintegration.H"
#include "time_integration_scheme.H"
#include "../linalg/linalg_ana.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_mapextractor.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_condition_utils.H"
#include "fluid_utils.H"
#include "fluidimpedancecondition.H"
#include "fluid_volumetric_surfaceFlow_condition.H"
#include "dyn_smag.H"
#include "turbulence_statistic_manager.H"
#include "fluid_utils_mapextractor.H"
#include "fluid_windkessel_optimization.H"
#include "fluid_meshtying.H"
#include "../drt_adapter/adapter_coupling_mortar.H"

#ifdef D_ARTNET
#include "../drt_art_net/art_net_dyn_drt.H"
#include "../drt_art_net/artnetexplicitintegration.H"
#include "fluid_coupling_red_models.H"
#endif // D_ARTNET

#ifdef WRITEOUTSTATISTICS
#include "../drt_io/io_control.H"
#endif

#include "../drt_io/io.H"
#include "../drt_io/io_control.H"

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
  write_wall_shear_stresses_(params.get<int>("write wall shear stresses", 0)),
  surfacesplitter_(NULL),
  inrelaxation_(false),
  msht_(INPAR::FLUID::no_meshtying)
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

  // physical type of fluid flow (incompressible, varying density, loma, Boussinesq approximation)
  physicaltype_ = DRT::INPUT::get<INPAR::FLUID::PhysicalType>(params_, "Physical Type");
  // type of time-integration
  timealgo_ = DRT::INPUT::get<INPAR::FLUID::TimeIntegrationScheme>(params_, "time int algo");
  // time-step size
  dtp_ = dta_ = params_.get<double>("time step size");
  // maximum number of timesteps
  stepmax_  = params_.get<int>   ("max number timesteps");
  // maximum simulation time
  maxtime_  = params_.get<double>("total time");
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

  predictor_ = params_.get<string>("predictor","steady_state_predictor");

  // form of convective term
  convform_ = params_.get<string>("form of convective term","convective");

  // conservative formulation currently not supported in low-Mach-number case
  // when using generalized-alpha time-integration scheme
  if (physicaltype_ == INPAR::FLUID::loma and timealgo_==INPAR::FLUID::timeint_afgenalpha and convform_ == "conservative")
     dserror("conservative formulation currently not supported for low-Mach-number flow within generalized-alpha time-integration scheme");

  // fine-scale subgrid viscosity?
  fssgv_ = params_.get<string>("fs subgrid viscosity","No");

  // -------------------------------------------------------------------
  // account for potential Neuman inflow terms if required
  // -------------------------------------------------------------------
  neumanninflow_ = false;
  if (params_.get<string>("Neumann inflow","no") == "yes") neumanninflow_ = true;

  // -------------------------------------------------------------------
  // care for periodic boundary conditions
  // -------------------------------------------------------------------

  // TODO: ??
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

  if (not params_.get<int>("Simple Preconditioner",0) && not params_.get<int>("AMG BS Preconditioner",0)
      && params_.get<int>("MESHTYING")== INPAR::FLUID::no_meshtying)
  {
    // initialize standard (stabilized) system matrix
    sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,108,false,true));
  }
  else if(params_.get<int>("MESHTYING")!= INPAR::FLUID::no_meshtying)
  {
    msht_ = params_.get<int>("MESHTYING");

    // define parameter list for meshtying
    ParameterList mshtparams;
    mshtparams.set("theta",theta_);
    mshtparams.set<int>("mshtoption", msht_);

    meshtying_ = Teuchos::rcp(new Meshtying(discret_, solver_, mshtparams, surfacesplitter_));
    sysmat_ = meshtying_->Setup();
    //meshtying_->OutputSetUp();
  }
  else
  {
    Teuchos::RCP<LINALG::BlockSparseMatrix<FLD::UTILS::VelPressSplitStrategy> > blocksysmat =
      Teuchos::rcp(new LINALG::BlockSparseMatrix<FLD::UTILS::VelPressSplitStrategy>(velpressplitter_,velpressplitter_,108,false,true));
    blocksysmat->SetNumdim(numdim_);
    sysmat_ = blocksysmat;
  }

  // sysmat might be singular (if we have a purely Dirichlet constrained
  // problem, the pressure mode is defined only up to a constant)
  // in this case, we need a basis vector for the nullspace/kernel
  vector<DRT::Condition*> KSPcond;
  discret_->GetCondition("KrylovSpaceProjection",KSPcond);
  int numcond = KSPcond.size();
  int numfluid = 0;
  for(int icond = 0; icond < numcond; icond++)
  {
    const std::string* name = KSPcond[icond]->Get<std::string>("discretization");
    if (*name == "fluid") numfluid++;
  }
  if (numfluid == 1)
  {
    project_ = true;
    w_       = LINALG::CreateVector(*dofrowmap,true);
    c_       = LINALG::CreateVector(*dofrowmap,true);
    kspsplitter_.Setup(*discret_);
  }
  else if (numfluid == 0)
  {
    project_ = false;
    w_       = Teuchos::null;
    c_       = Teuchos::null;
  }
  else
    dserror("Received more than one KrylovSpaceCondition for fluid field");

  // -------------------------------------------------------------------
  // create empty vectors
  // -------------------------------------------------------------------
  // additional rhs vector for robin-BC and vector for copying the residual
  robinrhs_ = LINALG::CreateVector(*dofrowmap,true);

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

  // scalar at time n+alpha_F/n+1 and n+alpha_M/n
  // (only required for low-Mach-number case)
  scaaf_ = LINALG::CreateVector(*dofrowmap,true);
  scaam_ = LINALG::CreateVector(*dofrowmap,true);

  // history vector
  hist_ = LINALG::CreateVector(*dofrowmap,true);

  if (alefluid_)
  {
    dispnp_ = LINALG::CreateVector(*dofrowmap,true);
    dispn_  = LINALG::CreateVector(*dofrowmap,true);
    dispnm_ = LINALG::CreateVector(*dofrowmap,true);
    gridv_  = LINALG::CreateVector(*dofrowmap,true);
  }

  // Vectors associated to boundary conditions
  // -----------------------------------------

  // create the volumetric-surface-flow condition
#ifdef D_ALE_BFLOW
  if (alefluid_)
  {
    discret_->SetState("dispnp", dispn_);
  }
#endif

  vol_surf_flow_bc_     = rcp(new UTILS::FluidVolumetricSurfaceFlowWrapper(discret_, output_, dta_) );

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

    // evaluate the map of te womersley bcs
    vol_surf_flow_bc_ -> EvaluateMapExtractor(vol_flow_rates_bc_extractor_);
    vol_surf_flow_bc_->EvaluateCondMap(vol_surf_flow_bcmaps_);


#ifdef D_ARTNET
    // -----------------------------------------------------------------
    // Initialize the reduced models
    // -----------------------------------------------------------------

    ART_exp_timeInt_ = dyn_art_net_drt(true);
    // Check if one-dimensional artery network problem exist
    if (ART_exp_timeInt_ != Teuchos::null)
    {
      IO::DiscretizationWriter output_redD(ART_exp_timeInt_->Discretization());
      discret_->ClearState();
      discret_->SetState("velnp", zeros_);
      if (alefluid_)
      {
        discret_->SetState("dispnp", dispnp_);
      }
      coupled3D_redDbc_art_=   rcp(new  UTILS::Fluid_couplingWrapper<ART::ArtNetExplicitTimeInt>
                                   ( discret_,
                                     ART_exp_timeInt_->Discretization(),
                                     ART_exp_timeInt_,
                                     output_redD,
                                     dta_,
                                     ART_exp_timeInt_->Dt()));

    }
#endif //D_ARTNET

#ifdef D_RED_AIRWAYS
    airway_imp_timeInt_ = dyn_red_airways_drt(true);
    // Check if one-dimensional artery network problem exist
    if (airway_imp_timeInt_ != Teuchos::null)
    {
      IO::DiscretizationWriter output_redD(airway_imp_timeInt_->Discretization());
      discret_->ClearState();
      discret_->SetState("velnp", zeros_);
      if (alefluid_)
      {
        discret_->SetState("dispnp", dispnp_);
      }
      coupled3D_redDbc_airways_ =   rcp(new  UTILS::Fluid_couplingWrapper<AIRWAY::RedAirwayImplicitTimeInt>
                                   ( discret_,
                                     airway_imp_timeInt_->Discretization(),
                                     airway_imp_timeInt_,
                                     output_redD,
                                     dta_,
                                     airway_imp_timeInt_->Dt()));

    }
#endif // D_RED_AIRWAYS

    // Evaluate the womersley velocities
    vol_surf_flow_bc_->EvaluateVelocities(velnp_,time_);


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
  scale_similarity_ = false;

  ParameterList *  modelparams =&(params_.sublist("TURBULENCE MODEL"));

  string physmodel = modelparams->get<string>("PHYSICAL_MODEL","no_model");

  // flag for special flow: currently channel flow or flow in a lid-driven cavity
  special_flow_ = modelparams->get<string>("CANONICAL_FLOW","no");

  // warning if classical (all-scale) turbulence model and fine-scale
  // subgrid-viscosity approach are intended to be used simultaneously
  if (fssgv_ != "No" and (physmodel == "Smagorinsky" or physmodel == "Dynamic_Smagorinsky" or physmodel == "Mixed_Scale_Similarity_Eddy_Viscosity_Model"))
    dserror("No combination of classical all-scale subgrid-viscosity turbulence model and fine-scale subgrid-viscosity approach currently possible!");

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
        else if(physmodel == "Scale_Similarity")
        {
          cout << "                             "      ;
          cout << "scale similarity model     ";
          cout << "Constant:   Cl   = "      ;
          cout << modelparams->get<double>("C_SCALE_SIMILARITY");
          cout << "\n                             ";
          if(fssgv_ != "No")
          {
             cout << "\n                             "      ;
             cout << "combined with      ";
          }
        }
        else if(physmodel == "Mixed_Scale_Similarity_Eddy_Viscosity_Model")
        {
          cout << "                             "      ;
          cout << "scale similarity model combined with constant Smagorinsy model     ";
          cout << "\n                             "      ;
          cout << "- Constant:   Cl   = "      ;
          cout << modelparams->get<double>("C_SCALE_SIMILARITY");
          cout << "\n                             "      ;
          cout << "- Constant:   Cs   = "      ;
          cout << modelparams->get<double>("C_SMAGORINSKY");
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

    if((modelparams->get<string>("PHYSICAL_MODEL","no_model")
       ==
       "Scale_Similarity") or (modelparams->get<string>("PHYSICAL_MODEL","no_model")
       ==
       "Mixed_Scale_Similarity_Eddy_Viscosity_Model")
      )
    {
      scale_similarity_ = true;

      const Epetra_Map* nodecolmap = discret_->NodeColMap();
      filteredvel_ = rcp(new Epetra_MultiVector(*nodecolmap,3,true));
      filteredreystr_ = rcp(new Epetra_MultiVector(*nodecolmap,9,true));

      // get one instance of the dynamic Smagorinsky class
      DynSmag_=rcp(new DynSmagFilter(discret_            ,
                                     pbcmapmastertoslave_,
                                     params_             ));
    }

  }

  // -------------------------------------------------------------------
  // necessary only for the VM3 approach:
  // fine-scale solution vector + respective ouptput
  // -------------------------------------------------------------------
  if (fssgv_ != "No")
  {
    fsvelaf_  = LINALG::CreateVector(*dofrowmap,true);

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

  // -------------------------------------------------------------------
  // initialize turbulence-statistics evaluation
  // -------------------------------------------------------------------
  //
  statisticsmanager_=rcp(new TurbulenceStatisticManager(*this));
  // parameter for sampling/dumping period
  if (special_flow_ != "no")
    samstart_ = modelparams->get<int>("SAMPLING_START",1);

  // ---------------------------------------------------------------------
  // set density variable to 1.0 and get gas constant for low-Mach-number
  // flow and get constant density variable for incompressible flow
  // ---------------------------------------------------------------------
  if (physicaltype_ == INPAR::FLUID::loma or physicaltype_ == INPAR::FLUID::varying_density)
  {
    // set density variable to 1.0 for low-Mach-number flow
    density_ = 1.0;

    // get gas constant
    ParameterList eleparams;
    eleparams.set("action","get_gas_constant");
    eleparams.set<int>("Physical Type", physicaltype_);
    discret_->Evaluate(eleparams,null,null,null,null,null);
    gasconstant_ = eleparams.get("gas constant", 1.0);
    // potential check here -> currently not executed
    //if (gasconstant_ < EPS15) dserror("received zero or negative gas constant");
  }
  else
  {
    // set gas constant to 1.0 for incompressible flow
    gasconstant_ = 1.0;

    // get constant density variable for incompressible flow
    ParameterList eleparams;
    eleparams.set("action","get_density");
    eleparams.set<int>("Physical Type", physicaltype_);
    discret_->Evaluate(eleparams,null,null,null,null,null);
    density_ = eleparams.get("density", 1.0);
    if (density_ < EPS15) dserror("received zero or negative density value");
  }

  // initialize all thermodynamic pressure values and its time derivative
  // to one or zero, respectively
  // -> they are kept this way for incompressible flow
  thermpressaf_   = 1.0;
  thermpressam_   = 1.0;
  thermpressdtam_ = 0.0;

#ifdef D_ALE_BFLOW
  if (alefluid_)
  {
    discret_->ClearState();
    discret_->SetState("dispnp", dispnp_);
  }
#endif // D_ALE_BFLOW
  // construct impedance bc wrapper
  impedancebc_     = rcp(new UTILS::FluidImpedanceWrapper(discret_, output_, dta_) );

  Wk_optimization_  = rcp(new UTILS::FluidWkOptimizationWrapper(discret_,
                                                               output_,
                                                              impedancebc_,
                                                               dta_) );

  // ---------------------------------------------------------------------
  // set general fluid parameter defined before
  // ---------------------------------------------------------------------
  SetElementGeneralFluidParameter();

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

    if (timealgo_!=INPAR::FLUID::timeint_stationary)
          cout <<  "                             " << "Tau Type        = " << stabparams->get<string>("DEFINITION_TAU") <<"\n";
    else
    {
      if(stabparams->get<string>("DEFINITION_TAU") == "Barrenechea_Franca_Valentin_Wall" or
          stabparams->get<string>("DEFINITION_TAU") == "Barrenechea_Franca_Valentin_Wall_wo_dt")
        cout <<  "                             " << "Tau             = " << "Barrenechea_Franca_Valentin_Wall_wo_dt" << "\n";
      else if (stabparams->get<string>("DEFINITION_TAU") == "Bazilevs_wo_dt" or
          stabparams->get<string>("DEFINITION_TAU") == "Bazilevs")
        cout <<  "                             " << "Tau             = " << "Bazilevs_wo_dt" << "\n";
    }
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
  if (timealgo_==INPAR::FLUID::timeint_stationary) SolveStationaryProblem();
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
        printf("TIME: %11.4E/%11.4E  DT = %11.4E  Generalized-Alpha  STEP = %4d/%4d \n",
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
    NonlinearSolve();

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
  if (timealgo_==INPAR::FLUID::timeint_bdf2) theta_ = (dta_+dtp_)/(2.0*dta_ + dtp_);

  // -------------------------------------------------------------------
  // set part(s) of the rhs vector(s) belonging to the old timestep
  // (only meaningful for momentum part)
  //
  // stationary/af-generalized-alpha: hist_ = 0.0
  //
  // one-step-Theta:                  hist_ = veln_  + dt*(1-Theta)*accn_
  //
  // BDF2: for constant time step:    hist_ = 4/3 veln_  - 1/3 velnm_
  //
  // -------------------------------------------------------------------
  TIMEINT_THETA_BDF2::SetOldPartOfRighthandside(veln_,velnm_, accn_,
                                        timealgo_, dta_, theta_, hist_);

  // -------------------------------------------------------------------
  //                     do explicit predictor step
  //
  // for example
  //
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

  if(extrapolationpredictor_)
  {
    if (step_>1)
    {
      TIMEINT_THETA_BDF2::ExplicitPredictor(
        predictor_,
        veln_,
        velnm_,
        accn_,
        velpressplitter_,
        timealgo_,
        theta_,
        dta_,
        dtp_,
        velnp_,
        discret_->Comm());
    }
  }

  // -------------------------------------------------------------------
  //  For af-generalized-alpha time-integration scheme:
  //  set "pseudo-theta", calculate initial accelerations according to
  //  prescribed Dirichlet values for generalized-alpha time
  //  integration and values at intermediate time steps
  // -------------------------------------------------------------------
  if (timealgo_==INPAR::FLUID::timeint_afgenalpha)
  {
    // starting algorithm
    if (startalgo_)
    {
      // use backward-Euler-type parameter combination
      if (step_<=numstasteps_)
      {
        cout<<"Starting algorithm for Af_GenAlpha active."
            <<"Performing step "<<step_ <<" of "<<numstasteps_<<
            " Backward Euler starting steps"<<endl;
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

#ifdef D_ALE_BFLOW
    if (alefluid_)
    {
      discret_->SetState("dispnp", dispnp_);
    }
#endif // D_ALE_BFLOW

#ifdef D_ARTNET
    // update the 3D-to-reduced_D coupling data
    // Check if one-dimensional artery network problem exist
    if (ART_exp_timeInt_ != Teuchos::null)
    {
      coupled3D_redDbc_art_->EvaluateDirichlet(velnp_, *(dbcmaps_->CondMap()), time_);
    }
#endif //D_ARTNET
#ifdef D_RED_AIRWAYS
    // update the 3D-to-reduced_D coupling data
    // Check if one-dimensional artery network problem exist
    if (airway_imp_timeInt_ != Teuchos::null)
    {
      coupled3D_redDbc_airways_->EvaluateDirichlet(velnp_, *(dbcmaps_->CondMap()), time_);
    }
#endif //D_RED_AIRWAYS

    // Evaluate the womersley velocities
    vol_surf_flow_bc_->EvaluateVelocities(velnp_,time_);

    discret_->ClearState();

    // set thermodynamic pressure
    eleparams.set("thermodynamic pressure",thermpressaf_);

    // evaluate Neumann conditions
    neumann_loads_->PutScalar(0.0);
    discret_->SetState("scaaf",scaaf_);
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

  // -------------------------------------------------------------------
  //           preparation of AVM3-based scale separation
  // -------------------------------------------------------------------
  if (step_==1 and fssgv_ != "No") AVM3Preparation();

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

  dtsolve_  = 0.0;
  dtele_    = 0.0;
  dtfilter_ = 0.0;

  if (myrank_ == 0)
  {
    printf("+------------+-------------------+--------------+--------------+--------------+--------------+\n");
    printf("|- step/max -|- tol      [norm] -|-- vel-res ---|-- pre-res ---|-- vel-inc ---|-- pre-inc ---|\n");
  }

  while (stopnonliniter==false)
  {
    itnum++;

    // -------------------------------------------------------------------
    // call elements to calculate system matrix and RHS
    // -------------------------------------------------------------------
    {
      // time measurement: element
      TEUCHOS_FUNC_TIME_MONITOR("      + element calls");

      // get cpu time
      const double tcpu=Teuchos::Time::wallTime();

      sysmat_->Zero();

      // create the parameters for the discretization
      ParameterList eleparams;

      // add Neumann loads
      residual_->Update(1.0,*neumann_loads_,0.0);

      // update impedance boundary condition
      impedancebc_->UpdateResidual(residual_);

#ifdef D_ARTNET
      // update the 3D-to-reduced_D coupling data
      // Check if one-dimensional artery network problem exist
      if (ART_exp_timeInt_ != Teuchos::null)
      {
        coupled3D_redDbc_art_->UpdateResidual(residual_);
      }
#endif //D_ARTNET
#ifdef D_RED_AIRWAYS
      // update the 3D-to-reduced_D coupling data
      // Check if one-dimensional artery network problem exist
      if (airway_imp_timeInt_ != Teuchos::null)
      {
        coupled3D_redDbc_airways_->UpdateResidual(residual_);
      }
#endif // D_RED_AIRWAYS


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
          const double tcpufilter=Teuchos::Time::wallTime();
          this->ApplyFilterForDynamicComputationOfCs();
          dtfilter_=Teuchos::Time::wallTime()-tcpufilter;
        }

        if (scale_similarity_)
        {
          //compute filtered velocity
          // time measurement
          const double tcpufilter=Teuchos::Time::wallTime();
          this->ApplyFilterForClassicalLES();
          dtfilter_=Teuchos::Time::wallTime()-tcpufilter;
        }
      }

      // Set action type
      eleparams.set("action","calc_fluid_systemmat_and_residual");

      // parameters for turbulent approach
      eleparams.sublist("TURBULENCE MODEL") = params_.sublist("TURBULENCE MODEL");
      if (scale_similarity_)
      {
        eleparams.set("Filtered velocity",filteredvel_);
        eleparams.set("Filtered reynoldsstress",filteredreystr_);
      }

      // set thermodynamic pressures
      eleparams.set("thermpress at n+alpha_F/n+1",thermpressaf_);
      eleparams.set("thermpress at n+alpha_M/n",thermpressam_);
      eleparams.set("thermpressderiv at n+alpha_M/n+1",thermpressdtam_);

      // set general vector values needed by elements
      discret_->ClearState();
      discret_->SetState("hist" ,hist_ );
      discret_->SetState("accam",accam_);
      discret_->SetState("scaaf",scaaf_);
      discret_->SetState("scaam",scaam_);
      if (alefluid_)
      {
        discret_->SetState("dispnp", dispnp_);
        discret_->SetState("gridv", gridv_);
      }

      // set scheme-specific element parameters and vector values
      if (timealgo_==INPAR::FLUID::timeint_afgenalpha)
           discret_->SetState("velaf",velaf_);
      else discret_->SetState("velaf",velnp_);

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
        discret_->Evaluate(eleparams,sysmat_,null,residual_,null,null);

        discret_->ClearState();


        //---------------------------surface tension update
        if (alefluid_ and surfacesplitter_->FSCondRelevant())
        {
          // employs the divergence theorem acc. to Saksono eq. (24) and does
          // not require second derivatives.

          // select free surface elements
          std::string condname = "FREESURFCoupling";

          ParameterList eleparams;

          // set action for elements
          eleparams.set("action","calc_surface_tension");

          discret_->ClearState();
          discret_->SetState("dispnp", dispnp_);
          discret_->EvaluateCondition(eleparams,Teuchos::null,Teuchos::null,residual_,Teuchos::null,Teuchos::null,condname);
          discret_->ClearState();
        }
        //---------------------------end of surface tension update

        //----------------------------------------------------------------------
        // apply mixed/hybrid Dirichlet boundary conditions
        //----------------------------------------------------------------------
        vector<DRT::Condition*> MHDcnd;
        discret_->GetCondition("SurfaceMixHybDirichlet",MHDcnd);

        if(MHDcnd.size()!=0)
        {
          ParameterList mhdbcparams;

          // set action for elements
          mhdbcparams.set("action"    ,"MixedHybridDirichlet");

          // set the only required state vectors
          if (timealgo_==INPAR::FLUID::timeint_afgenalpha)
          {
            discret_->SetState("u and p (trial)",velaf_);
          }
          else discret_->SetState("u and p (trial)",velnp_);

          // evaluate all mixed hybrid Dirichlet boundary conditions
          discret_->EvaluateConditionUsingParentData
            (mhdbcparams          ,
             sysmat_              ,
             Teuchos::null        ,
             residual_            ,
             Teuchos::null        ,
             Teuchos::null        ,
             "LineMixHybDirichlet");

          discret_->EvaluateConditionUsingParentData
            (mhdbcparams          ,
             sysmat_              ,
             Teuchos::null        ,
             residual_            ,
             Teuchos::null        ,
             Teuchos::null        ,
             "SurfaceMixHybDirichlet");

          // clear state
          discret_->ClearState();

        }

        //----------------------------------------------------------------------
        // account for potential Neumann inflow terms
        //----------------------------------------------------------------------
        if (neumanninflow_)
        {
          // create parameter list
          ParameterList condparams;

          // action for elements
          condparams.set("action","calc_Neumann_inflow");

          // set thermodynamic pressure
          condparams.set("thermpress at n+alpha_F/n+1",thermpressaf_);

          // set vector values needed by elements
          discret_->ClearState();
          discret_->SetState("scaaf",scaaf_);
          if (timealgo_==INPAR::FLUID::timeint_afgenalpha)
               discret_->SetState("velaf",velaf_);
          else discret_->SetState("velaf",velnp_);

          std::string condstring("FluidNeumannInflow");
          discret_->EvaluateCondition(condparams,sysmat_,Teuchos::null,residual_,Teuchos::null,Teuchos::null,condstring);
          discret_->ClearState();
        }

        // scaling to get true residual vector
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
          surfacesplitter_->AddFSICondVector(-1.*scale,
                                             surfacesplitter_->ExtractFSICondVector(velnp_),
                                             residual_);

          // Note: It is the right thing to test the robin enhanced residual_
          // for convergence, since the velocity terms are to vanish and the
          // structural forces are to cancel with the internal forces.
          //
          // Note: We do not add any external (robin) loads to
          // trueresidual_. This way we get the unbalanced forces at the
          // interface, which can be applied to the structure later on.

          const Epetra_Map& robinmap = *surfacesplitter_->FSICondMap();
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
      dtele_=Teuchos::Time::wallTime()-tcpu;
    }

    //double dtsysmatsplit_=0.0;
    if(msht_ != INPAR::FLUID::no_meshtying)
      meshtying_->PrepareMeshtyingSystem(sysmat_, residual_);

    // blank residual DOFs which are on Dirichlet BC
    // We can do this because the values at the dirichlet positions
    // are not used anyway.
    // We could avoid this though, if velrowmap_ and prerowmap_ would
    // not include the dirichlet values as well. But it is expensive
    // to avoid that.
    dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), residual_);

    // Treat the surface volumetric flow rate
    //    RCP<Epetra_Vector> temp_vec = rcp(new Epetra_Vector(*vol_surf_flow_bcmaps_,true));
    //    vol_surf_flow_bc_->InsertCondVector( *temp_vec , *residual_);
    vol_flow_rates_bc_extractor_->InsertVolumetricSurfaceFlowCondVector(
      vol_flow_rates_bc_extractor_->ExtractVolumetricSurfaceFlowCondVector(zeros_),
      residual_);


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

        for(int rr=0;rr<numdim_;++rr)
        {
          if(abs((*mode)[rr]>1e-14))
          {
            dserror("expecting only an undetermined pressure");
          }
        }

        int predof = numdim_;

        Teuchos::RCP<Epetra_Vector> presmode = velpressplitter_.ExtractCondVector(*w_);

        presmode->PutScalar((*mode)[predof]);

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
        Teuchos::RCP<Epetra_Vector> tmpw = LINALG::CreateVector(*(discret_->DofRowMap()),true);
        LINALG::Export(*presmode,*tmpw);
        Teuchos::RCP<Epetra_Vector> tmpkspw = kspsplitter_.ExtractKSPCondVector(*tmpw);
        LINALG::Export(*tmpkspw,*w_);

        // export to vector of ones
        presmode->PutScalar(1.0);
        Teuchos::RCP<Epetra_Vector> tmpc = LINALG::CreateVector(*(discret_->DofRowMap()),true);
        LINALG::Export(*presmode,*tmpc);
        Teuchos::RCP<Epetra_Vector> tmpkspc = kspsplitter_.ExtractKSPCondVector(*tmpc);
        LINALG::Export(*tmpkspc,*c_);
      }
      else if(*definition == "integration")
      {
        // zero w and c
        w_->PutScalar(0.0);
        c_->PutScalar(0.0);

        ParameterList mode_params;

        // set action for elements
        mode_params.set("action","integrate_shape");

        if (alefluid_)
        {
          discret_->SetState("dispnp",dispnp_);
        }

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

        discret_->EvaluateCondition
          (mode_params        ,
           Teuchos::null      ,
           Teuchos::null      ,
           w_                 ,
           Teuchos::null      ,
           Teuchos::null      ,
           "KrylovSpaceProjection");

        // get pressure
        const vector<double>* mode = KSPcond->Get<vector<double> >("mode");

        for(int rr=0;rr<numdim_;++rr)
        {
          if(abs((*mode)[rr]>1e-14))
          {
            dserror("expecting only an undetermined pressure");
          }
        }

        Teuchos::RCP<Epetra_Vector> presmode = velpressplitter_.ExtractCondVector(*w_);

        // export to vector of ones
        presmode->PutScalar(1.0);
        Teuchos::RCP<Epetra_Vector> tmpc = LINALG::CreateVector(*(discret_->DofRowMap()),true);
        LINALG::Export(*presmode,*tmpc);
        Teuchos::RCP<Epetra_Vector> tmpkspc = kspsplitter_.ExtractKSPCondVector(*tmpc);
        LINALG::Export(*tmpkspc,*c_);

        if(msht_!= INPAR::FLUID::no_meshtying)
          meshtying_->KrylovProjection(c_);
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
        printf(" (      --     ,te=%10.3E",dtele_);
        if (dynamic_smagorinsky_ or scale_similarity_)
        {
          printf(",tf=%10.3E",dtfilter_);
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
          printf(" (ts=%10.3E,te=%10.3E",dtsolve_,dtele_);
          if (dynamic_smagorinsky_ or scale_similarity_)
          {
            printf(",tf=%10.3E",dtfilter_);
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
          printf(" (ts=%10.3E,te=%10.3E",dtsolve_,dtele_);
          if (dynamic_smagorinsky_ or scale_similarity_)
          {
            printf(",tf=%10.3E",dtfilter_);
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
    {
      // apply the womersley velocity profile as a dirichlet bc
      LINALG::ApplyDirichlettoSystem(sysmat_,incvel_,residual_,zeros_,*(vol_surf_flow_bcmaps_));
    }
    //-------solve for residual displacements to correct incremental displacements
    {
      // time measurement: solver
      TEUCHOS_FUNC_TIME_MONITOR("      + solver calls");

      // get cpu time
      const double tcpusolve=Teuchos::Time::wallTime();

      // do adaptive linear solver tolerance (not in first solve)
      if (isadapttol && itnum>1)
      {
        double currresidual = max(vresnorm,presnorm);
        currresidual = max(currresidual,incvelnorm_L2/velnorm_L2);
        currresidual = max(currresidual,incprenorm_L2/prenorm_L2);
        solver_.AdaptTolerance(ittol,currresidual,adaptolbetter);
      }

#ifdef WRITEOUTSTATISTICS
    FILE* errfile = params_.get<FILE*>("err file",NULL);
    if(errfile!=NULL)
    {
      fprintf(errfile, "TOBI: Proc %i/%i\tTimeStep %i\tNonlinIter %i\t",myrank_,discret_->Comm().NumProc(),step_,itnum);
    }
#endif

      if (msht_!= INPAR::FLUID::no_meshtying)
        meshtying_->SolveMeshtying(solver_, sysmat_, incvel_, residual_, itnum, w_, c_, project_);
      else
        solver_.Solve(sysmat_->EpetraOperator(),incvel_,residual_,true,itnum==1, w_, c_, project_);

      solver_.ResetTolerance();

      // end time measurement for solver
      dtsolve_ = Teuchos::Time::wallTime()-tcpusolve;
    }

    // -------------------------------------------------------------------
    // update velocity and pressure values by increments
    // -------------------------------------------------------------------
    velnp_->Update(1.0,*incvel_,1.0);

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

    //------------------------------------------------ free surface update
    if (alefluid_ and surfacesplitter_->FSCondRelevant())
    {
      Teuchos::RefCountPtr<Epetra_Vector> fsvelnp = surfacesplitter_->ExtractFSCondVector(velnp_);
      Teuchos::RefCountPtr<Epetra_Vector> fsdisp = surfacesplitter_->ExtractFSCondVector(dispn_);
      Teuchos::RefCountPtr<Epetra_Vector> fsdispnp = Teuchos::rcp(new Epetra_Vector(*surfacesplitter_->FSCondMap()));

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

      // ================ HEIGHTFUNCTION =======================================
      // This is the heightfunction implementation for the partitioned algorithm
      // as it is. This is a very basic implementation - and it is not
      // mass-consistent! We need another evaluate call here to find a
      // mass-consistent heightfunction manipulation matrix.
      //
      // local lagrange is mass-consistent of course.

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
        Teuchos::RefCountPtr<Epetra_Vector> ndnorm = surfacesplitter_->ExtractFSCondVector(ndnorm0);

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

      surfacesplitter_->InsertFSCondVector(fsdispnp,dispnp_);
      surfacesplitter_->InsertFSCondVector(fsvelnp,gridv_);
    }
  }
} // FluidImplicitTimeInt::NonlinearSolve



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | predictor                                                   vg 02/09 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidImplicitTimeInt::Predictor()
{
  // -------------------------------------------------------------------
  // time measurement: nonlinear iteration
  // -------------------------------------------------------------------
  TEUCHOS_FUNC_TIME_MONITOR("   + predictor");

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

  // Treat the surface volumetric flow rate
  //    RCP<Epetra_Vector> temp_vec = rcp(new Epetra_Vector(*vol_surf_flow_bcmaps_,true));
  //    vol_surf_flow_bc_->InsertCondVector( *temp_vec , *residual_);
  vol_flow_rates_bc_extractor_->InsertVolumetricSurfaceFlowCondVector(
    vol_flow_rates_bc_extractor_->ExtractVolumetricSurfaceFlowCondVector(zeros_),
    residual_);


  double vresnorm;
  double presnorm;

  Teuchos::RCP<Epetra_Vector> onlyvel = velpressplitter_.ExtractOtherVector(residual_);
  onlyvel->Norm2(&vresnorm);

  Teuchos::RCP<Epetra_Vector> onlypre = velpressplitter_.ExtractCondVector(residual_);
  onlypre->Norm2(&presnorm);

  if (myrank_ == 0)
  {
    printf("+------------+-------------------+--------------+--------------+--------------+--------------+\n");
    printf("| predictor  |  vel. | pre. res. | %10.3E   | %10.3E   |      --      |      --      |",vresnorm,presnorm);
    printf(" (      --     ,te=%10.3E",dtele_);
    if (dynamic_smagorinsky_ or scale_similarity_) printf(",tf=%10.3E",dtfilter_);
    printf(")\n");
    printf("+------------+-------------------+--------------+--------------+--------------+--------------+\n");
  }

} // FluidImplicitTimeInt::Predictor


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | (multiple) corrector                                        vg 02/09 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidImplicitTimeInt::MultiCorrector()
{
  // -------------------------------------------------------------------
  // time measurement: nonlinear iteration
  // -------------------------------------------------------------------
  TEUCHOS_FUNC_TIME_MONITOR("   + corrector");

  dtsolve_  = 0.0;

  // -------------------------------------------------------------------
  // parameters and variables for nonlinear iteration
  // -------------------------------------------------------------------
  const double ittol = params_.get<double>("tolerance for nonlin iter");
  int          itnum = 0;
  int          itemax = 0;
  bool         stopnonliniter = false;
  double       incvelnorm_L2;
  double       incprenorm_L2;
  double       velnorm_L2;
  double       prenorm_L2;
  double       vresnorm;
  double       presnorm;

  // -------------------------------------------------------------------
  // currently default for turbulent channel flow:
  // only one iteration before sampling
  // -------------------------------------------------------------------
  if (special_flow_ == "channel_flow_of_height_2" && step_<samstart_ )
       itemax = 1;
  else itemax = params_.get<int>("max nonlin iter steps");

  // -------------------------------------------------------------------
  // turn adaptive solver tolerance on/off
  // -------------------------------------------------------------------
  const bool   isadapttol    = params_.get<bool>("ADAPTCONV",true);
  const double adaptolbetter = params_.get<double>("ADAPTCONV_BETTER",0.01);

  // -------------------------------------------------------------------
  // prepare print out for (multiple) corrector
  // -------------------------------------------------------------------
  if (myrank_ == 0)
  {
    printf("+------------+-------------------+--------------+--------------+--------------+--------------+\n");
    printf("|- step/max -|- tol      [norm] -|-- vel-res ---|-- pre-res ---|-- vel-inc ---|-- pre-inc ---|\n");
  }

  // -------------------------------------------------------------------
  // nonlinear iteration loop
  // -------------------------------------------------------------------
  while (stopnonliniter==false)
  {
    itnum++;

    // -------------------------------------------------------------------
    // call elements to calculate system matrix and rhs and assemble
    // -------------------------------------------------------------------
    AssembleMatAndRHS();

    // -------------------------------------------------------------------
    // apply Dirichlet boundary conditions to system of equations:
    // - Residual displacements are supposed to be zero for resp. dofs.
    // - Time for applying Dirichlet boundary conditions is measured.
    // -------------------------------------------------------------------
    incvel_->PutScalar(0.0);
    {
      TEUCHOS_FUNC_TIME_MONITOR("      + apply DBC");
      LINALG::ApplyDirichlettoSystem(sysmat_,incvel_,residual_,zeros_,*(dbcmaps_->CondMap()));
    }
    {
      // apply the womersley velocity profile as a dirichlet bc
      LINALG::ApplyDirichlettoSystem(sysmat_,incvel_,residual_,zeros_,*(vol_surf_flow_bcmaps_));
    }

    // -------------------------------------------------------------------
    // solve for velocity and pressure increments
    // - Adaptive linear solver tolerance is used from second
    //   corrector step on.
    // - Time for solver is measured.
    // -------------------------------------------------------------------
    {
      TEUCHOS_FUNC_TIME_MONITOR("      + solver calls");

      const double tcpusolve=Teuchos::Time::wallTime();

      if (isadapttol and itnum>1)
      {
        double currresidual = max(vresnorm,presnorm);
        currresidual = max(currresidual,incvelnorm_L2/velnorm_L2);
        currresidual = max(currresidual,incprenorm_L2/prenorm_L2);
        solver_.AdaptTolerance(ittol,currresidual,adaptolbetter);
      }

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

          for(int rr=0;rr<numdim_;++rr)
          {
            if(abs((*mode)[rr]>1e-14))
            {
              dserror("expecting only an undetermined pressure");
            }
          }

          int predof = numdim_;

          Teuchos::RCP<Epetra_Vector> presmode = velpressplitter_.ExtractCondVector(*w_);

          presmode->PutScalar((*mode)[predof]);

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
          Teuchos::RCP<Epetra_Vector> tmpw = LINALG::CreateVector(*(discret_->DofRowMap()),true);
          LINALG::Export(*presmode,*tmpw);
          Teuchos::RCP<Epetra_Vector> tmpkspw = kspsplitter_.ExtractKSPCondVector(*tmpw);
          LINALG::Export(*tmpkspw,*w_);

          // export to vector of ones
          presmode->PutScalar(1.0);
          Teuchos::RCP<Epetra_Vector> tmpc = LINALG::CreateVector(*(discret_->DofRowMap()),true);
          LINALG::Export(*presmode,*tmpc);
          Teuchos::RCP<Epetra_Vector> tmpkspc = kspsplitter_.ExtractKSPCondVector(*tmpc);
          LINALG::Export(*tmpkspc,*c_);
        }
        else if(*definition == "integration")
        {
          // zero w and c
          w_->PutScalar(0.0);
          c_->PutScalar(0.0);

          ParameterList mode_params;

          // set action for elements
          mode_params.set("action","integrate_shape");

          if (alefluid_)
          {
            discret_->SetState("dispnp",dispnp_);
          }

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


          discret_->EvaluateCondition
          (mode_params        ,
           Teuchos::null      ,
           Teuchos::null      ,
           w_               ,
           Teuchos::null      ,
           Teuchos::null      ,
           "KrylovSpaceProjection");

          // get pressure
          const vector<double>* mode = KSPcond->Get<vector<double> >("mode");

          for(int rr=0;rr<numdim_;++rr)
          {
            if(abs((*mode)[rr]>1e-14))
            {
              dserror("expecting only an undetermined pressure");
            }
          }

          Teuchos::RCP<Epetra_Vector> presmode = velpressplitter_.ExtractCondVector(*w_);

          // export to vector of ones
          presmode->PutScalar(1.0);
          Teuchos::RCP<Epetra_Vector> tmpc = LINALG::CreateVector(*(discret_->DofRowMap()),true);
          LINALG::Export(*presmode,*tmpc);
          Teuchos::RCP<Epetra_Vector> tmpkspc = kspsplitter_.ExtractKSPCondVector(*tmpc);
          LINALG::Export(*tmpkspc,*c_);
        }
        else
        {
          dserror("unknown definition of weight vector w for restriction of Krylov space");
        }
      }

      solver_.Solve(sysmat_->EpetraOperator(),incvel_,residual_,true,itnum==1, w_, c_, project_);

      solver_.ResetTolerance();

      dtsolve_ = Teuchos::Time::wallTime()-tcpusolve;
    }

    // -------------------------------------------------------------------
    // update velocity and pressure values by increments
    // -------------------------------------------------------------------
    velnp_->Update(1.0,*incvel_,1.0);

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

    // Treat the surface volumetric flow rate
    //    RCP<Epetra_Vector> temp_vec = rcp(new Epetra_Vector(*vol_surf_flow_bcmaps_,true));
    //    vol_surf_flow_bc_->InsertCondVector( *temp_vec , *residual_);
    vol_flow_rates_bc_extractor_->InsertVolumetricSurfaceFlowCondVector(
      vol_flow_rates_bc_extractor_->ExtractVolumetricSurfaceFlowCondVector(zeros_),
      residual_);

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

    // care for the case that nothing really happens in velocity
    // or pressure field
    if (velnorm_L2 < 1e-5) velnorm_L2 = 1.0;
    if (prenorm_L2 < 1e-5) prenorm_L2 = 1.0;

    if (myrank_ == 0)
    {
      printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   | %10.3E   | %10.3E   |",itnum,itemax,ittol,vresnorm,presnorm,
                  incvelnorm_L2/velnorm_L2,incprenorm_L2/prenorm_L2);
      printf(" (ts=%10.3E,te=%10.3E",dtsolve_,dtele_);
      if (dynamic_smagorinsky_ or scale_similarity_) printf(",tf=%10.3E",dtfilter_);
      printf(")\n");
    }

    // -------------------------------------------------------------------
    // check convergence and print out respective information:
    // - stop if convergence is achieved
    // - warn if itemax is reached without convergence, but proceed to
    //   next timestep
    // -------------------------------------------------------------------
    if (vresnorm <= ittol and
        presnorm <= ittol and
        incvelnorm_L2/velnorm_L2 <= ittol and
        incprenorm_L2/prenorm_L2 <= ittol)
    {
      stopnonliniter=true;
      if (myrank_ == 0)
      {
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

    if ((itnum == itemax) and (vresnorm > ittol or
                               presnorm > ittol or
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
  }

} // FluidImplicitTimeInt::MultiCorrector


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | compute values at intermediate time steps for gen.-alpha    vg 02/09 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidImplicitTimeInt::GenAlphaIntermediateValues()
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

} // FluidImplicitTimeInt::GenAlphaIntermediateValues


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | call elements to calculate system matrix/rhs and assemble   vg 02/09 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidImplicitTimeInt::AssembleMatAndRHS()
{
  dtele_    = 0.0;
  dtfilter_ = 0.0;

  // time measurement: element
  TEUCHOS_FUNC_TIME_MONITOR("      + element calls");

  // get cpu time
  const double tcpu=Teuchos::Time::wallTime();

  sysmat_->Zero();

  // create the parameters for the discretization
  ParameterList eleparams;

  // add Neumann loads
  residual_->Update(1.0,*neumann_loads_,0.0);

  // update impedance boundary condition
  impedancebc_->UpdateResidual(residual_);

  // update the 3D-to-reduced_D coupling condition
#ifdef D_ARTNET
  // Check if one-dimensional artery network problem exist
  if (ART_exp_timeInt_ != Teuchos::null)
  {
    coupled3D_redDbc_art_->UpdateResidual(residual_);
  }
#endif //D_ARTNET
#ifdef D_RED_AIRWAYS
  // Check if one-dimensional artery network problem exist
  if (airway_imp_timeInt_ != Teuchos::null)
  {
    coupled3D_redDbc_airways_->UpdateResidual(residual_);
  }
#endif // D_RED_AIRWAYS

  if (dynamic_smagorinsky_)
  {
    // time measurement
    const double tcpufilter=Teuchos::Time::wallTime();
    this->ApplyFilterForDynamicComputationOfCs();
    dtfilter_=Teuchos::Time::wallTime()-tcpufilter;
  }

  if (scale_similarity_)
  {
    //compute filtered velocity
    // time measurement
    const double tcpufilter=Teuchos::Time::wallTime();
    this->ApplyFilterForClassicalLES();
    dtfilter_=Teuchos::Time::wallTime()-tcpufilter;
  }

  // set action type
  eleparams.set("action","calc_fluid_systemmat_and_residual");

  // parameters for turbulence model
  eleparams.sublist("TURBULENCE MODEL") = params_.sublist("TURBULENCE MODEL");
  if (scale_similarity_)
  {
    eleparams.set("Filtered velocity",filteredvel_);
    eleparams.set("Filtered reynoldsstress",filteredreystr_);
  }

  // set thermodynamic pressures
  eleparams.set("thermpress at n+alpha_F/n+1",thermpressaf_);
  eleparams.set("thermpress at n+alpha_M/n",thermpressam_);
  eleparams.set("thermpressderiv at n+alpha_M/n+1",thermpressdtam_);

  // set general vector values needed by elements
  discret_->ClearState();
  discret_->SetState("hist" ,hist_ );
  discret_->SetState("accam",accam_);
  discret_->SetState("scaaf",scaaf_);
  discret_->SetState("scaam",scaam_);
  if (alefluid_)
  {
    discret_->SetState("dispnp", dispnp_);
    discret_->SetState("gridv", gridv_);
  }

  // set scheme-specific element parameters and vector values
  if (timealgo_==INPAR::FLUID::timeint_afgenalpha)
       discret_->SetState("velaf",velaf_);
  else discret_->SetState("velaf",velnp_);

  //----------------------------------------------------------------------
  // AVM3-based solution approach if required
  //----------------------------------------------------------------------
  if (fssgv_ != "No") AVM3Separation();

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

    // set thermodynamic pressure
    condparams.set("thermpress at n+alpha_F/n+1",thermpressaf_);

    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("scaaf",scaaf_);
    if (timealgo_==INPAR::FLUID::timeint_afgenalpha)
         discret_->SetState("velaf",velaf_);
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

} // FluidImplicitTimeInt::AssembleMatAndRHS


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | update acceleration for generalized-alpha time integration  vg 02/09 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidImplicitTimeInt::GenAlphaUpdateAcceleration()
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

} // FluidImplicitTimeInt::GenAlphaUpdateAcceleration


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
  if (shapederivatives_ != Teuchos::null)
    shapederivatives_->Zero();

  // set the new solution we just got
  if (vel!=Teuchos::null)
  {
    // Take Dirichlet values from velnp and add vel to veln for non-Dirichlet
    // values.
    Teuchos::RCP<Epetra_Vector> aux = LINALG::CreateVector(*(discret_->DofRowMap()),true);
    aux->Update(1.0, *veln_, 1.0, *vel, 0.0);
    //    dbcmaps_->InsertOtherVector(dbcmaps_->ExtractOtherVector(aux), velnp_);
    dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(velnp_), aux);

    //
    vol_flow_rates_bc_extractor_->InsertVolumetricSurfaceFlowCondVector(
      vol_flow_rates_bc_extractor_->ExtractVolumetricSurfaceFlowCondVector(velnp_),
      aux);

    *velnp_ = *aux;

  }


  // add Neumann loads
  residual_->Update(1.0,*neumann_loads_,0.0);

  // update impedance boundary condition
  impedancebc_->UpdateResidual(residual_);

  // update the 3D-to-reduce_D coupling condition
#ifdef D_ARTNET
  // Check if one-dimensional artery network problem exist
  if (ART_exp_timeInt_ != Teuchos::null)
  {
    coupled3D_redDbc_art_->UpdateResidual(residual_);
  }
#endif // D_ARTNET
#ifdef D_RED_AIRWAYS
  // Check if one-dimensional artery network problem exist
  if (airway_imp_timeInt_ != Teuchos::null)
  {
    coupled3D_redDbc_airways_->UpdateResidual(residual_);
  }
#endif // D_RED_AIRWAYS

  // create the parameters for the discretization
  ParameterList eleparams;

  // set action type
  eleparams.set("action","calc_fluid_systemmat_and_residual");

  // set thermodynamic pressures
  eleparams.set("thermpress at n+alpha_F/n+1",thermpressaf_);
  eleparams.set("thermpress at n+alpha_M/n",thermpressam_);
  eleparams.set("thermpressderiv at n+alpha_M/n+1",thermpressdtam_);

  // set general vector values needed by elements
  discret_->ClearState();
  discret_->SetState("hist" ,hist_ );
  discret_->SetState("accam",accam_);
  discret_->SetState("scaaf",scaaf_);
  discret_->SetState("scaam",scaam_);
  if (alefluid_)
  {
    discret_->SetState("dispnp", dispnp_);
    discret_->SetState("gridv", gridv_);
  }

  // set scheme-specific element parameters and vector values
  if (timealgo_==INPAR::FLUID::timeint_afgenalpha)
       discret_->SetState("velaf",velaf_);
  else discret_->SetState("velaf",velnp_);

  // call loop over elements
  discret_->Evaluate(eleparams,sysmat_,shapederivatives_,residual_,Teuchos::null,Teuchos::null);
  discret_->ClearState();

  //---------------------------surface tension update
  if (alefluid_ and surfacesplitter_->FSCondRelevant())
  {
    // employs the divergence theorem acc. to Saksono eq. (24) and does not
    // require second derivatives.

    // select free surface elements
    std::string condname = "FREESURFCoupling";

    ParameterList eleparams;

    // set action for elements
    eleparams.set("action","calc_surface_tension");

    discret_->ClearState();
    discret_->SetState("dispnp", dispnp_);
    discret_->EvaluateCondition(eleparams,Teuchos::null,Teuchos::null,residual_,Teuchos::null,Teuchos::null,condname);
    discret_->ClearState();
  }
  //---------------------------end of surface tension update

  // finalize the system matrix
  sysmat_->Complete();

  if (shapederivatives_ != Teuchos::null)
  {
    shapederivatives_->Complete();
    // apply Dirichlet conditions to a non-diagonal matrix
    // (The Dirichlet rows will become all zero, no diagonal one.)
    shapederivatives_->ApplyDirichlet(*(dbcmaps_->CondMap()),false);


    // apply the womersley bc as a dirichlet bc
    shapederivatives_->ApplyDirichlet(*(vol_surf_flow_bcmaps_),false);
  }

  trueresidual_->Update(ResidualScaling(),*residual_,0.0);

  // Apply dirichlet boundary conditions to system of equations
  // residual displacements are supposed to be zero at boundary
  // conditions
  incvel_->PutScalar(0.0);

  LINALG::ApplyDirichlettoSystem(sysmat_,incvel_,residual_,zeros_,*(dbcmaps_->CondMap()));

  // apply Womersley as a Dirichlet BC
  LINALG::ApplyDirichlettoSystem(sysmat_,incvel_,residual_,zeros_,*(vol_surf_flow_bcmaps_));

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
void FLD::FluidImplicitTimeInt::TimeUpdate()
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
    if (timealgo_==INPAR::FLUID::timeint_afgenalpha)
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
    else
    {

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

    TIMEINT_THETA_BDF2::CalculateAcceleration(onlyvelnp,
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

  if (msht_== INPAR::FLUID::sps_coupled or msht_== INPAR::FLUID::sps_pc)
    meshtying_->UpdateLag();

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

#ifdef D_ALE_BFLOW
  if (alefluid_)
  {
    discret_->SetState("dispnp", dispn_);
  }
#endif //D_ALE_BFLOW

  impedancebc_->FlowRateCalculation(time_,dta_);
  impedancebc_->OutflowBoundary(time_,dta_,theta_);

  // get the parameters needed to be optimized
  ParameterList WkOpt_params;
  WkOpt_params.set<double> ("total time", time_);
  WkOpt_params.set<double> ("time step size", dta_);
  impedancebc_->getResultsOfAPeriod(WkOpt_params);

  // update wind kessel optimization condition
  Wk_optimization_->Solve(WkOpt_params);

  // -------------------------------------------------------------------
  // treat the 3D-to-reduced_D couplign condition
  // note: these methods return without action, if the problem does not
  //       have any coupling boundary conditions
  // -------------------------------------------------------------------


#ifdef D_ALE_BFLOW
  if (alefluid_)
  {
    discret_->SetState("dispnp", dispnp_);
  }
#endif // D_ALE_BFLOW

#ifdef D_ARTNET
  // Check if one-dimensional artery network problem exist
  if (ART_exp_timeInt_ != Teuchos::null)
  {
    coupled3D_redDbc_art_->FlowRateCalculation(time_,dta_);
    coupled3D_redDbc_art_->ApplyBoundaryConditions(time_, dta_, theta_);
  }
#endif //D_ARTNET

#ifdef D_RED_AIRWAYS

  // Check if one-dimensional artery network problem exist
  if (airway_imp_timeInt_ != Teuchos::null)
  {
    coupled3D_redDbc_airways_->FlowRateCalculation(time_,dta_);
    coupled3D_redDbc_airways_->ApplyBoundaryConditions(time_, dta_, theta_);
  }
#endif // D_RED_AIRWAYS
  discret_->ClearState();

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

  // compute equation-of-state factor
  const double eosfac = thermpressaf_/gasconstant_;
  // store subfilter stresses for additional output
  if (scale_similarity_)
  {
    RCP<Epetra_Vector> stress12 = CalcSFS(1,2);
    statisticsmanager_->StoreNodalValues(step_, stress12);
  }
  // -------------------------------------------------------------------
  //   add calculated velocity to mean value calculation (statistics)
  // -------------------------------------------------------------------
  statisticsmanager_->DoTimeSample(step_,time_,eosfac);

  // -------------------------------------------------------------------
  //                        compute flow rates
  // -------------------------------------------------------------------
  ComputeFlowRates();

  // -------------------------------------------------------------------
  //                         output of solution
  // -------------------------------------------------------------------
  Output();

  // -------------------------------------------------------------------
  //          dumping of turbulence statistics if required
  // -------------------------------------------------------------------
  statisticsmanager_->DoOutput(output_,step_,eosfac);

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

  //  ART_exp_timeInt_->Output();
  // output of solution
  if (step_%upres_ == 0)
  {
    // step number and time
    output_.NewStep(step_,time_);

    // velocity/pressure vector
    output_.WriteVector("velnp",velnp_);

    // (hydrodynamic) pressure
    Teuchos::RCP<Epetra_Vector> pressure = velpressplitter_.ExtractCondVector(velnp_);
    pressure->Scale(density_);
    output_.WriteVector("pressure", pressure);

    //output_.WriteVector("residual", trueresidual_);
    if (alefluid_) output_.WriteVector("dispnp", dispnp_);

    if (physicaltype_ == INPAR::FLUID::varying_density or physicaltype_ == INPAR::FLUID::boussinesq)
    {
      Teuchos::RCP<Epetra_Vector> scalar_field = velpressplitter_.ExtractCondVector(scaaf_);
      output_.WriteVector("scalar_field", scalar_field);
    }

    //only perform stress calculation when output is needed
    if (writestresses_)
    {
      RCP<Epetra_Vector> traction = CalcStresses();
      output_.WriteVector("traction",traction);
      if (myrank_==0)
        cout<<"Writing stresses"<<endl;
      //only perform wall shear stress calculation when output is needed
      if (write_wall_shear_stresses_)
      {
        RCP<Epetra_Vector> wss = CalcWallShearStresses();
        output_.WriteVector("wss",wss);
      }
    }

    if (scale_similarity_ or dynamic_smagorinsky_)
    {
      const Epetra_Map* dofrowmap = discret_->DofRowMap();
      RCP<Epetra_Vector> filteredvel = LINALG::CreateVector(*dofrowmap,true);
      DynSmag_->OutputofAveragedVel(filteredvel);
      output_.WriteVector("filteredvel",filteredvel);
      if (scale_similarity_)
      {
        if (myrank_==0)
           std::cout << "output of subfilter stresses for scale similarity model ..." << std::endl;
        RCP<Epetra_Vector> stress11 = CalcSFS(1,1);
        output_.WriteVector("sfs11",stress11);
        RCP<Epetra_Vector> stress12 = CalcSFS(1,2);
        output_.WriteVector("sfs12",stress12);
        RCP<Epetra_Vector> stress13 = CalcSFS(1,3);
        output_.WriteVector("sfs13",stress13);
        RCP<Epetra_Vector> stress22 = CalcSFS(2,2);
        output_.WriteVector("sfs22",stress22);
        RCP<Epetra_Vector> stress23 = CalcSFS(2,3);
        output_.WriteVector("sfs23",stress23);
        RCP<Epetra_Vector> stress33 = CalcSFS(3,3);
        output_.WriteVector("sfs33",stress33);
      }
    }

    if (fssgv_ != "No")
    {
      output_.WriteVector("fsvelaf",fsvelaf_);
    }

    // write domain decomposition for visualization (only once!)
    if (step_==upres_) output_.WriteElementData();

    if (uprestart_ != 0 && step_%uprestart_ == 0) //add restart data
    {
      // acceleration vector at time n+1 and n, velocity/pressure vector at time n and n-1
      output_.WriteVector("accnp",accnp_);
      output_.WriteVector("accn", accn_);
      output_.WriteVector("veln", veln_);
      output_.WriteVector("velnm",velnm_);

      if (alefluid_)
      {
        output_.WriteVector("dispn", dispn_);
        output_.WriteVector("dispnm",dispnm_);
      }

      // also write impedance bc information if required
      // Note: this method acts only if there is an impedance BC
      impedancebc_->WriteRestart(output_);

      Wk_optimization_->WriteRestart(output_);
    }

    vol_surf_flow_bc_->Output(output_);

  }
  // write restart also when uprestart_ is not a integer multiple of upres_
  else if (uprestart_ != 0 && step_%uprestart_ == 0)
  {
    // step number and time
    output_.NewStep(step_,time_);

    // velocity/pressure vector
    output_.WriteVector("velnp",velnp_);

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
      //only perform wall shear stress calculation when output is needed
      if (write_wall_shear_stresses_)
      {
        RCP<Epetra_Vector> wss = CalcWallShearStresses();
        output_.WriteVector("wss",wss);
      }
    }

    // acceleration vector at time n+1 and n, velocity/pressure vector at time n and n-1
    output_.WriteVector("accnp",accnp_);
    output_.WriteVector("accn", accn_);
    output_.WriteVector("veln", veln_);
    output_.WriteVector("velnm",velnm_);

    // also write impedance bc information if required
    // Note: this method acts only if there is an impedance BC
    impedancebc_->WriteRestart(output_);

    Wk_optimization_->WriteRestart(output_);
    vol_surf_flow_bc_->Output(output_);
  }

  // write reduced model problem
#ifdef D_ARTNET
  // Check if one-dimensional artery network problem exist
  if (ART_exp_timeInt_ != Teuchos::null)
  {
    RCP<ParameterList> redD_export_params;
    redD_export_params = rcp(new ParameterList());

    redD_export_params->set<int>("step",step_);
    redD_export_params->set<int>("upres",upres_);
    redD_export_params->set<int>("uprestart",uprestart_);
    redD_export_params->set<double>("time",time_);

    ART_exp_timeInt_->Output(true, redD_export_params);
  }
#endif // D_ARTNET

#ifdef D_RED_AIRWAYS
  // Check if one-dimensional artery network problem exist
  if (airway_imp_timeInt_ != Teuchos::null)
  {
    RCP<ParameterList> redD_export_params;
    redD_export_params = rcp(new ParameterList());

    redD_export_params->set<int>("step",step_);
    redD_export_params->set<int>("upres",upres_);
    redD_export_params->set<int>("uprestart",uprestart_);
    redD_export_params->set<double>("time",time_);

    airway_imp_timeInt_->Output(true, redD_export_params);
  }
#endif // D_AIRWAYS

//#define PRINTALEDEFORMEDNODECOORDS // flag for printing all ALE nodes and xspatial in current configuration - only works for 1 processor  devaal 02.2011

// output ALE nodes and xspatial in current configuration - devaal 02.2011
#ifdef PRINTALEDEFORMEDNODECOORDS

    if (discret_->Comm().NumProc() != 1)
	dserror("The flag PRINTALEDEFORMEDNODECOORDS has been switched on, and only works for 1 processor");

    cout << "ALE DISCRETIZATION IN THE DEFORMED CONFIGURATIONS" << endl;
    // does discret_ exist here?
    //cout << "discret_->NodeRowMap()" << discret_->NodeRowMap() << endl;

    //RCP<Epetra_Vector> mynoderowmap = rcp(new Epetra_Vector(discret_->NodeRowMap()));
    //RCP<Epetra_Vector> noderowmap_ = rcp(new Epetra_Vector(discret_->NodeRowMap()));
    //dofrowmap_  = rcp(new discret_->DofRowMap());
    const Epetra_Map* noderowmap = discret_->NodeRowMap();
    const Epetra_Map* dofrowmap = discret_->DofRowMap();

    for (int lid=0; lid<noderowmap->NumGlobalPoints(); lid++)
	{
	    int gid;
	    // get global id of a node
	    gid = noderowmap->GID(lid);
	    // get the node
	    DRT::Node * node = discret_->gNode(gid);
	    // get the coordinates of the node
	    const double * X = node->X();
	    // get degrees of freedom of a node
	    vector<int> gdofs = discret_->Dof(node);
	    //cout << "for node:" << *node << endl;
	    //cout << "this is my gdof vector" << gdofs[0] << " " << gdofs[1] << " " << gdofs[2] << endl;

	    // get displacements of a node
	    vector<double> mydisp (3,0.0);
	    for (int ldof = 0; ldof<3; ldof ++)
	    {
	    	int displid = dofrowmap->LID(gdofs[ldof]);
	    	//cout << "displacement local id - in the rowmap" << displid << endl;
	    	mydisp[ldof] = (*dispnp_)[displid];
	    	//make zero if it is too small
	    	if (abs(mydisp[ldof]) < 0.00001)
	    	    {
			mydisp[ldof] = 0.0;
	    	    }
	    }
	    // Export disp, X
	    double newX = mydisp[0]+X[0];
	    double newY = mydisp[1]+X[1];
	    double newZ = mydisp[2]+X[2];
	    //cout << "NODE " << gid << "  COORD  " << newX << " " << newY << " " << newZ << endl;
	    cout << gid << " " << newX << " " << newY << " " << newZ << endl;
	}
#endif //PRINTALEDEFORMEDNODECOORDS

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
  //  ART_exp_timeInt_->ReadRestart(step);
  IO::DiscretizationReader reader(discret_,step);
  time_ = reader.ReadDouble("time");
  step_ = reader.ReadInt("step");

  reader.ReadVector(velnp_,"velnp");
  reader.ReadVector(veln_, "veln");
  reader.ReadVector(velnm_,"velnm");
  reader.ReadVector(accnp_,"accnp");
  reader.ReadVector(accn_ ,"accn");

  statisticsmanager_->Restart(reader,step);

  if (fssgv_ != "No") AVM3Preparation();

  if (alefluid_)
  {
    reader.ReadVector(dispnp_,"dispnp");
    reader.ReadVector(dispn_ , "dispn");
    reader.ReadVector(dispnm_,"dispnm");
  }
  // also read impedance bc information if required
  // Note: this method acts only if there is an impedance BC
  impedancebc_->ReadRestart(reader);

  Wk_optimization_->ReadRestart(reader);

  vol_surf_flow_bc_->ReadRestart(reader);

#ifdef D_ARTNET
  // Check if one-dimensional artery network problem exist
  if (ART_exp_timeInt_ != Teuchos::null)
  {
    ART_exp_timeInt_->ReadRestart(step_);
  }
#endif // D_ARTNET
#ifdef D_RED_AIRWAYS
  // Check if one-dimensional artery network problem exist
  if (airway_imp_timeInt_ != Teuchos::null)
  {
    airway_imp_timeInt_->ReadRestart(step_);
  }
#endif // D_RED_AIRWAYS

  // Read restart of one-dimensional arterial network
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
      gridv_->Update(1.5/dta_, *dispnp_, -2.0/dta_, *dispn_, 0.0);
      gridv_->Update(0.5/dta_, *dispnm_, 1.0);
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
  // time measurement: avm3
  TEUCHOS_FUNC_TIME_MONITOR("           + avm3");

  // zero matrix
  sysmat_->Zero();

  // add Neumann loads
  residual_->Update(1.0,*neumann_loads_,0.0);

  // create the parameters for the discretization
  ParameterList eleparams;

  // set action type
  eleparams.set("action","calc_fluid_systemmat_and_residual");

  // parameters for turbulence approach
  eleparams.sublist("TURBULENCE MODEL") = params_.sublist("TURBULENCE MODEL");
  if (scale_similarity_)
  {
    //compute filtered velocity
    //set filtered velocity
    this->ApplyFilterForClassicalLES();
    eleparams.set("Filtered velocity",filteredvel_);
    eleparams.set("Filtered reynoldsstress",filteredreystr_);
  }

  // set thermodynamic pressures
  eleparams.set("thermpress at n+alpha_F/n+1",thermpressaf_);
  eleparams.set("thermpress at n+alpha_M/n",thermpressam_);
  eleparams.set("thermpressderiv at n+alpha_M/n+1",thermpressdtam_);

  // set general vector values needed by elements
  discret_->ClearState();
  discret_->SetState("hist" ,hist_ );
  discret_->SetState("accam",accam_);
  discret_->SetState("scaaf",scaaf_);
  discret_->SetState("scaam",scaam_);

  // set fine-scale vector: at this
  discret_->SetState("fsvelaf",fsvelaf_);

  if (alefluid_)
  {
    discret_->SetState("dispnp", dispnp_);
    discret_->SetState("gridv", gridv_);
  }

  // set scheme-specific element parameters and vector values
  if (timealgo_==INPAR::FLUID::timeint_afgenalpha)
       discret_->SetState("velaf",velaf_);
  else discret_->SetState("velaf",velnp_);

  // element evaluation for getting system matrix
  // -> we merely need matrix "structure" below, not the actual contents
  discret_->Evaluate(eleparams,sysmat_,null,residual_,null,null);
  discret_->ClearState();

  // complete system matrix
  sysmat_->Complete();

  // apply DBC to system matrix
  LINALG::ApplyDirichlettoSystem(sysmat_,incvel_,residual_,zeros_,*(dbcmaps_->CondMap()));

  // apply Womersley as a Dirichlet BC
  LINALG::ApplyDirichlettoSystem(sysmat_,incvel_,residual_,zeros_,*(vol_surf_flow_bcmaps_));

  // get scale-separation matrix
  {
    // this is important to have!!!
    MLAPI::Init();

    // extract the ML parameters
    ParameterList&  mlparams = solver_.Params().sublist("ML Parameters");;

    // get toggle vector for Dirchlet boundary conditions
    const Epetra_Vector& dbct = *Dirichlet();

    // get nullspace parameters
    double* nullspace = mlparams.get("null space: vectors",(double*)NULL);
    if (!nullspace) dserror("No nullspace supplied in parameter list");
    int nsdim = mlparams.get("null space: dimension",1);

    // modify nullspace to ensure that DBC are fully taken into account
    if (nullspace)
    {
      const int length = SystemMatrix()->OperatorRangeMap().NumMyElements();
      for (int i=0; i<nsdim; ++i)
        for (int j=0; j<length; ++j)
          if (dbct[j]!=0.0) nullspace[i*length+j] = 0.0;
    }

    // get plain aggregation Ptent
    RCP<Epetra_CrsMatrix> crsPtent;
    GetPtent(*SystemMatrix()->EpetraMatrix(),mlparams,nullspace,crsPtent);
    LINALG::SparseMatrix Ptent(crsPtent);

    // compute scale-separation matrix: S = I - Ptent*Ptent^T
    Sep_ = LINALG::Multiply(Ptent,false,Ptent,true);
    Sep_->Scale(-1.0);
    RCP<Epetra_Vector> tmp = LINALG::CreateVector(Sep_->RowMap(),false);
    tmp->PutScalar(1.0);
    RCP<Epetra_Vector> diag = LINALG::CreateVector(Sep_->RowMap(),false);
    Sep_->ExtractDiagonalCopy(*diag);
    diag->Update(1.0,*tmp,1.0);
    Sep_->ReplaceDiagonalValues(*diag);

    //complete scale-separation matrix and check maps
    Sep_->Complete(Sep_->DomainMap(),Sep_->RangeMap());
    if (!Sep_->RowMap().SameAs(SystemMatrix()->RowMap())) dserror("rowmap not equal");
    if (!Sep_->RangeMap().SameAs(SystemMatrix()->RangeMap())) dserror("rangemap not equal");
    if (!Sep_->DomainMap().SameAs(SystemMatrix()->DomainMap())) dserror("domainmap not equal");
  }

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
  // time measurement: avm3
  TEUCHOS_FUNC_TIME_MONITOR("           + avm3");

  // get fine-scale part of velocity at time n+alpha_F or n+1
  if (timealgo_==INPAR::FLUID::timeint_afgenalpha)
    Sep_->Multiply(false,*velaf_,*fsvelaf_);
  else
    Sep_->Multiply(false,*velnp_,*fsvelaf_);

  // set fine-scale vector
  discret_->SetState("fsvelaf",fsvelaf_);

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

    // for isogeometric problems we have more effort...
    DRT::NURBS::apply_nurbs_initial_condition(
      *discret_  ,
      solver_    ,
      startfuncno,
      velnp_     );

    // initialize veln_ as well.
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
 | set fields for low-Mach-number flow within iteration loop   vg 09/09 |
 *----------------------------------------------------------------------*/
void FLD::FluidImplicitTimeInt::SetIterLomaFields(
   RCP<const Epetra_Vector> scalaraf,
   RCP<const Epetra_Vector> scalaram,
   RCP<const Epetra_Vector> scalardtam,
   const double             thermpressaf,
   const double             thermpressam,
   const double             thermpressdtam,
   Teuchos::RCP<DRT::Discretization> scatradis)
{
  // initializations
  int err(0);
  double value(0.0);
  vector<int> nodedofs;

  //--------------------------------------------------------------------------
  // Filling the scaaf-vector and scaam-vector at time n+alpha_F/n+1 and
  // n+alpha_M/n, respectively, with scalar at pressure dofs
  // Additionally, filling the scaam-vector at time n+alpha_M/n with
  // velocity at time n at velocity dofs for OST/BDF2
  // Filling the accam-vector at time n+alpha_M/n+1, respectively, with
  // scalar time derivative values at pressure dofs
  //--------------------------------------------------------------------------
  // get velocity values at time n in scaam-vector as copy from veln-vector
  scaam_->Update(1.0,*veln_,0.0);

  // loop all nodes on the processor
  for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
  {
    // get the processor's local scatra node
    DRT::Node* lscatranode = scatradis->lRowNode(lnodeid);

    // find out the global dof id of the last(!) dof at the scatra node
    const int numscatradof = scatradis->NumDof(lscatranode);
    const int globalscatradofid = scatradis->Dof(lscatranode,numscatradof-1);
    const int localscatradofid = scalaraf->Map().LID(globalscatradofid);
    if (localscatradofid < 0)
      dserror("localdofid not found in map for given globaldofid");

    // get the processor's local fluid node
    DRT::Node* lnode = discret_->lRowNode(lnodeid);
    // get the global ids of degrees of freedom associated with this node
    nodedofs = discret_->Dof(lnode);
    // get global and processor's local pressure dof id (using the map!)
    const int numdof = discret_->NumDof(lnode);
    const int globaldofid = discret_->Dof(lnode,numdof-1);
    const int localdofid = scaam_->Map().LID(globaldofid);
    if (localscatradofid < 0)
      dserror("localdofid not found in map for given globaldofid");

    // now copy the values
    value = (*scalaraf)[localscatradofid];
    err = scaaf_->ReplaceMyValue(localdofid,0,value);
    if (err != 0) dserror("error while inserting value into scaaf_");

    value = (*scalaram)[localscatradofid];
    err = scaam_->ReplaceMyValue(localdofid,0,value);
    if (err != 0) dserror("error while inserting value into scaam_");

    if (scalardtam != Teuchos::null)
    {
      value = (*scalardtam)[localscatradofid];
    }
    else
    {
      value = 0.0; // for safety reasons: set zeros in accam_
    }
    err = accam_->ReplaceMyValue(localdofid,0,value);
    if (err != 0) dserror("error while inserting value into accam_");
  }

  //--------------------------------------------------------------------------
  // get thermodynamic pressure at n+alpha_F/n+1 and n+alpha_M/n
  // and time derivative of thermodynamic pressure at n+alpha_M/n+1
  //--------------------------------------------------------------------------
  thermpressaf_   = thermpressaf;
  thermpressam_   = thermpressam;
  thermpressdtam_ = thermpressdtam;

  return;

} // ScaTraTimIntImpl::SetIterLomaFields


/*----------------------------------------------------------------------*
 | set fields for low-Mach-number flow at end of time step     vg 09/09 |
 *----------------------------------------------------------------------*/
void FLD::FluidImplicitTimeInt::SetTimeLomaFields(
   RCP<const Epetra_Vector> scalarnp,
   const double             thermpressnp,
   RCP<const Epetra_Vector> scatraresidual,
   Teuchos::RCP<DRT::Discretization> scatradis)
{
  // initializations
  int err(0);
  double value(0.0);
  vector<int> nodedofs;

  //--------------------------------------------------------------------------
  // Filling the scaaf-vector with scalar at time n+1 at pressure dofs
  //--------------------------------------------------------------------------
  // loop all nodes on the processor
  for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
  {
    // get the processor's local scatra node
    DRT::Node* lscatranode = scatradis->lRowNode(lnodeid);

    // find out the global dof id of the last(!) dof at the scatra node
    const int numscatradof = scatradis->NumDof(lscatranode);
    const int globalscatradofid = scatradis->Dof(lscatranode,numscatradof-1);
    const int localscatradofid = scalarnp->Map().LID(globalscatradofid);
    if (localscatradofid < 0)
      dserror("localdofid not found in map for given globaldofid");

    // get the processor's local fluid node
    DRT::Node* lnode = discret_->lRowNode(lnodeid);
    // get the global ids of degrees of freedom associated with this node
    nodedofs = discret_->Dof(lnode);
    // get global and processor's local pressure dof id (using the map!)
    const int globaldofid = nodedofs[numdim_];
    const int localdofid = scaam_->Map().LID(globaldofid);
    if (localdofid < 0)
      dserror("localdofid not found in map for given globaldofid");

    value = (*scalarnp)[localscatradofid];
    err = scaaf_->ReplaceMyValue(localdofid,0,value);
    if (err != 0) dserror("error while inserting value into scaaf_");

    //--------------------------------------------------------------------------
    // Filling the trueresidual vector with scatraresidual at pre-dofs
    //--------------------------------------------------------------------------
    if (scatraresidual != Teuchos::null)
    {
      value = (*scatraresidual)[localscatradofid];
      trueresidual_->ReplaceMyValue(localdofid,0,value);
    }

  }

  //--------------------------------------------------------------------------
  // get thermodynamic pressure at n+1
  //--------------------------------------------------------------------------
  thermpressaf_ = thermpressnp;


  return;

} // ScaTraTimIntImpl::SetTimeLomaFields


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
    Teuchos::RCP<Epetra_SerialDenseVector> errors
#if 0
      = Teuchos::rcp(new Epetra_SerialDenseVector(numdim_+2));
#else
      = Teuchos::rcp(new Epetra_SerialDenseVector(2));
#endif

    // call loop over elements (assemble nothing)
    discret_->EvaluateScalars(eleparams, errors);
    discret_->ClearState();

    double velerr = 0.0;
    double preerr = 0.0;

    // for the L2 norm, we need the square root
    velerr = sqrt((*errors)[0]);
    preerr = sqrt((*errors)[1]);

    if (myrank_ == 0)
    {
      printf("\n  L2_err for beltrami flow:  velocity %15.8e  pressure %15.8e\n\n",
             velerr,preerr);

#if 0
      //Write error in a file
      // following headers need to be included
      // #include "../drt_io/io.H"
      // #include "../drt_io/io_control.H"

      double velerrx = 0.0;
      double velerry = 0.0;
      double velerrz = 0.0;

      velerrx = sqrt((*errors)[2]);
      velerry = sqrt((*errors)[3]);
      velerrz = sqrt((*errors)[4]);

      if ((step_==stepmax_) or (time_==maxtime_))// write results to file
      {
        ostringstream temp;
        const std::string simulation = DRT::Problem::Instance()->OutputControlFile()->FileName();
        const std::string fname = "error.txt";

        std::ofstream f;
        f.open(fname.c_str(),std::fstream::ate | std::fstream::app);
        f << "#| " << simulation << "\n";
        f << "#| Step | Time | L2-error velocity x | L2-error velocity y  | L2-error pressure|\n";
        f << step_ << " " << time_ << " " << velerr << " " << velerrx << " " << velerry << " " << velerrz << " " << preerr << " " <<"\n";
        f.flush();
        f.close();
      }

      ostringstream temp;
      const std::string simulation = DRT::Problem::Instance()->OutputControlFile()->FileName();
      const std::string fname = "error_time_"+simulation+".txt";

      if(step_==1)
      {
        std::ofstream f;
        f.open(fname.c_str());
        f << step_ << " " << time_ << " " << velerr << " " << velerrx << " " << velerry << " " << velerrz << " " << preerr << " " <<"\n";
        f.flush();
        f.close();
      }
      else
      {
        std::ofstream f;
        f.open(fname.c_str(),std::fstream::ate | std::fstream::app);
        f << step_ << " " << time_ << " " << velerr << " " << velerrx << " " << velerry << " " << velerrz << " " << preerr << " " <<"\n";
        f.flush();
        f.close();
      }
#endif
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
    printf("Stationary Fluid Solver - STEP = %4d/%4d \n",step_,stepmax_);
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
#ifdef D_ALE_BFLOW
        if (alefluid_)
        {
          discret_->SetState("dispnp", dispnp_);
        }
#endif // D_ALE_BFLOW
#ifdef D_ARTNET
      // update the 3D-to-reduced_D coupling data
      // Check if one-dimensional artery network problem exist
      if (ART_exp_timeInt_ != Teuchos::null)
      {
        coupled3D_redDbc_art_->EvaluateDirichlet(velnp_,*(dbcmaps_->CondMap()), time_);
      }
#endif //D_ARTNET
#ifdef D_RED_AIRWAYS
      // update the 3D-to-reduced_D coupling data
      // Check if one-dimensional artery network problem exist
      if (airway_imp_timeInt_ != Teuchos::null)
      {
        coupled3D_redDbc_airways_->EvaluateDirichlet(velnp_,*(dbcmaps_->CondMap()), time_);
      }
#endif // D_RED_AIRWAYS


      // Evaluate the womersley velocities
      vol_surf_flow_bc_->EvaluateVelocities(velnp_,time_);

      discret_->ClearState();
      // evaluate Neumann b.c.
      //eleparams.set("inc_density",density_);

      // set thermodynamic pressure
      eleparams.set("thermodynamic pressure",thermpressaf_);

      neumann_loads_->PutScalar(0.0);
      discret_->SetState("scaaf",scaaf_);
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
    //                        compute flow rates
    // -------------------------------------------------------------------
    ComputeFlowRates();

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
/*
calculate lift&drag forces and angular moments

Lift and drag forces are based upon the right hand side true-residual entities
of the corresponding nodes. The contribution of the end node of a line is entirely
added to a present L&D force.

Notice: Angular moments obtained from lift&drag forces currently refer to the
        initial configuration, i.e. are built with the coordinates X of a particular
        node irrespective of its current position.
*/
void FLD::FluidImplicitTimeInt::LiftDrag() const
{
  // in this map, the results of the lift drag calculation are stored
  RCP<map<int,vector<double> > > liftdragvals;

  FLD::UTILS::LiftDrag(*discret_,*trueresidual_,params_,liftdragvals);

  if (liftdragvals!=Teuchos::null and discret_->Comm().MyPID() == 0)
    FLD::UTILS::WriteLiftDragToFile(time_, step_, *liftdragvals);

  return;
}


/*----------------------------------------------------------------------*
 | compute flow rates through desired boundary parts        u.may 01/10 |
 *----------------------------------------------------------------------*/
void FLD::FluidImplicitTimeInt::ComputeFlowRates() const
{
  vector<DRT::Condition*> flowratecond;
  string condstring;

  if(numdim_ == 2)
  {
    condstring = "LineFlowRate";
    discret_->GetCondition("LineFlowRate", flowratecond);
    // if no flowrate condition is present we do not compute anything
    if((int) flowratecond.size()== 0)
      return;
  }
  else if (numdim_ == 3)
  {
    condstring = "SurfFlowRate";
    discret_->GetCondition("SurfFlowRate", flowratecond);
    // if no flowrate condition is present we do not compute anything
    if((int) flowratecond.size()== 0)
      return;
  }
  else
    dserror("flow rate computation is not implemented for the 1D case");

  const std::map<int,double> flowrates = FLD::UTILS::ComputeFlowRates(*discret_, velnp_, condstring);

  // write to file
  if(discret_->Comm().MyPID() == 0)
    FLD::UTILS::WriteFlowRatesToFile(time_, step_, flowrates );

  return;
}

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
 | filter quantities for classical LES models                           |
 *----------------------------------------------------------------------*/
void FLD::FluidImplicitTimeInt::ApplyFilterForClassicalLES()
{
  // perform filtering
  {
    const Teuchos::RCP<const Epetra_Vector> dirichtoggle = Dirichlet();
    // call only filtering
    if (timealgo_==INPAR::FLUID::timeint_afgenalpha)
    {
      DynSmag_->ApplyFilter(velaf_,dirichtoggle);
    }
    else
    {
      DynSmag_->ApplyFilter(velnp_,dirichtoggle);
    }
    // get filtered fields
    filteredvel_->PutScalar(0.0);
    filteredreystr_->PutScalar(0.0);
    //std::cout << "get filtered fields" << std::endl;
    DynSmag_->GetFilteredVelocity(filteredvel_);
    //std::cout << "get filtered vel" << std::endl;
    DynSmag_->GetFilteredReynoldsStress(filteredreystr_);
    //std::cout << "get filtered str" << std::endl;
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
void FLD::FluidImplicitTimeInt::UseBlockMatrix(Teuchos::RCP<std::set<int> >     condelements,
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
    shapederivatives_ = mat;
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
  // convergence in the Newton loop.
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

    // general fluid and time parameter are set in PrepareTimeStep()
    ParameterList eleparams;

    // parameters for stabilization
    eleparams.sublist("TURBULENCE MODEL") = params_.sublist("TURBULENCE MODEL");

    // set thermodynamic pressures
    eleparams.set("thermpress at n+alpha_F/n+1",thermpressaf_);
    eleparams.set("thermpress at n+alpha_M/n",thermpressam_);
    eleparams.set("thermpressderiv at n+alpha_M/n+1",thermpressdtam_);

    // set general vector values needed by elements
    discret_->ClearState();
    discret_->SetState("hist" ,hist_ );
    discret_->SetState("accam",accam_);
    discret_->SetState("scaaf",scaaf_);
    discret_->SetState("scaam",scaam_);
    discret_->SetState("dispnp", griddisp);
    discret_->SetState("gridv", zeros_);

    eleparams.set("action","calc_fluid_systemmat_and_residual");
    // set scheme-specific element parameters and vector values
    if (timealgo_==INPAR::FLUID::timeint_afgenalpha)
         discret_->SetState("velaf",velaf_);
    else discret_->SetState("velaf", velnp_);

    // call loop over elements
    discret_->Evaluate(eleparams,sysmat_,meshmatrix_,residual_,null,null);
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
    dirichletlines_ = SystemMatrix()->ExtractDirichletRows(*(dbcmaps_->CondMap()));
    sysmat_->ApplyDirichlet(*(dbcmaps_->CondMap()));

    // apply Womersley as a Dirichlet BC
    sysmat_->ApplyDirichlet(*(vol_surf_flow_bcmaps_));
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

  // apply Womersley as a Dirichlet BC
  LINALG::ApplyDirichlettoSystem(incvel_,residual_,relax,*(vol_surf_flow_bcmaps_));

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

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> FLD::FluidImplicitTimeInt::VelocityRowMap()
{ return velpressplitter_.OtherMap(); }

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> FLD::FluidImplicitTimeInt::PressureRowMap()
{ return velpressplitter_.CondMap(); }


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  calculate wall sheer stress at (Dirichlet) boundary (public)        |
 |                                                          ismail 08/10|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
Teuchos::RCP<Epetra_Vector> FLD::FluidImplicitTimeInt::CalcWallShearStresses()
{
  // -------------------------------------------------------------------
  // first evaluate the normals at the nodes
  // -------------------------------------------------------------------

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
  if (alefluid_)
  {
    discret_->SetState("dispnp", dispnp_);
  }
  // evaluate the normals of the surface
  discret_->EvaluateCondition(eleparams,ndnorm0,"FluidStressCalc");
  discret_->ClearState();

  // -------------------------------------------------------------------
  // normalise the normal vectors
  // -------------------------------------------------------------------
  for (int i = 0; i < ndnorm0->MyLength();i+=numdim_+1)
  {
    // calculate the length of the normal
    double L = 0.0;
    for (int j = 0; j<numdim_; j++)
    {
      L += ((*ndnorm0)[i+j])*((*ndnorm0)[i+j]);
    }
    L = sqrt(L);

    // normalise the normal vector (if present for the current node)
    if (L > EPS15)
    {
      for (int j = 0; j < numdim_; j++)
      {
        (*ndnorm0)[i+j] /=  L;
      }
    }
  }

  // -------------------------------------------------------------------
  // evaluate the wall shear stress from the traction by removing
  // the normal stresses
  // -------------------------------------------------------------------

  // get traction
  RCP<Epetra_Vector> wss = CalcStresses();

  // loop over all entities within the traction vector
  for (int i = 0; i < ndnorm0->MyLength();i+=numdim_+1)
  {
    // evaluate the normal stress = < traction . normal >
    double normal_stress = 0.0;
    for (int j = 0; j<numdim_; j++)
    {
      normal_stress += (*wss)[i+j] * (*ndnorm0)[i+j];
    }

    // subtract the normal stresses from traction
    for (int j = 0; j<numdim_; j++)
    {
      (*wss)[i+j] -= normal_stress * (*ndnorm0)[i+j];
    }
  }


  // -------------------------------------------------------------------
  // return the wall_shear_stress vector
  // -------------------------------------------------------------------
  return wss;
}

Teuchos::RCP<Epetra_Vector> FLD::FluidImplicitTimeInt::CalcSFS(
   const int   i,
   const int   j
)
{
  // compute filtered quantities with the latest velocity field
  this->ApplyFilterForClassicalLES();

  const Epetra_Map* noderowmap = discret_->NodeRowMap();
  // create vector (+ initialization with zeros)
  Teuchos::RCP<Epetra_Vector> tauSFS = LINALG::CreateVector(*noderowmap,true);

  // get filtered velocity
  RCP<Epetra_Vector> filteredvel = LINALG::CreateVector(*noderowmap,true);
  DynSmag_->FilteredVelComp(filteredvel, i, j);

  // get filtered reynoldsstress
  RCP<Epetra_Vector> filteredreystr = LINALG::CreateVector(*noderowmap,true);
  DynSmag_->FilteredReyStrComp(filteredreystr, i , j);

  tauSFS->Update(1.0,*filteredreystr,-1.0, *filteredvel,0.0);

  return tauSFS;
}

// -------------------------------------------------------------------
// set general fluid parameter (AE 01/2011)
// -------------------------------------------------------------------

void FLD::FluidImplicitTimeInt::SetElementGeneralFluidParameter()
{
  ParameterList eleparams;

  eleparams.set("action","set_general_fluid_parameter");

  // set general element parameters
  eleparams.set("form of convective term",convform_);
  eleparams.set("fs subgrid viscosity",fssgv_);
  eleparams.set<int>("Linearisation",newton_);
  eleparams.set<int>("Physical Type", physicaltype_);

  // parameter for stabilization
  eleparams.sublist("STABILIZATION") = params_.sublist("STABILIZATION");

  // parameter for turbulent flow
  eleparams.sublist("TURBULENCE MODEL") = params_.sublist("TURBULENCE MODEL");

  //set time integration scheme
  eleparams.set<int>("TimeIntegrationScheme", timealgo_);

  // call standard loop over elements
  discret_->Evaluate(eleparams,null,null,null,null,null);
  return;
}


// -------------------------------------------------------------------
// set general time parameter (AE 01/2011)
// -------------------------------------------------------------------

void FLD::FluidImplicitTimeInt::SetElementTimeParameter()
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
  else if (timealgo_==INPAR::FLUID::timeint_afgenalpha)
  {
    eleparams.set("total time",time_-(1-alphaF_)*dta_);
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


#endif /* CCADISCRET       */
