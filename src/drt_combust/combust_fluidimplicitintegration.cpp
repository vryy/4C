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

\level 2

<pre>
\maintainer Ursula Rasthofer
            rasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/

#include "combust_fluidimplicitintegration.H"
#include "combust3_interpolation.H"
#include "combust_defines.H"
#include "combust_flamefront.H"
#include "combust_fluidresulttest.H"
#include "two_phase_defines.H"

#include "../drt_combust/combust_utils_time_integration.H"
#include "../drt_fluid/fluid_utils.H"
#include "../drt_lib/drt_periodicbc.H"
#include "../drt_fluid_turbulence/drt_transfer_turb_inflow.H"
#include "../drt_fluid_turbulence/turbulence_statistic_manager.H"
#include "../drt_geometry/integrationcell_coordtrafo.H"
#include "../drt_geometry/position_array.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_utils_parallel.H"
#include "../drt_lib/drt_dofset_independent_pbc.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret_combust.H"
#include "../drt_lib/drt_condition.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/matlist.H"
#include "../drt_xfem/dof_distribution_switcher.H"
#include "../drt_xfem/timeInt_std_SemiLagrange.H"
#include "../drt_xfem/timeInt_std_extrapolation.H"
#include "../drt_xfem/timeInt_enr.H"
#include "../linalg/linalg_ana.H"
#include "../linalg/linalg_krylov_projector.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../drt_comm/comm_utils.H"

// for AVM3-based scale separation
#include <MLAPI_Workspace.h>
#include <MLAPI_Aggregation.h>

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
  combusttype_(DRT::INPUT::IntegralValue<INPAR::COMBUST::CombustionType>(params_->sublist("COMBUSTION FLUID"),"COMBUSTTYPE")),
  veljumptype_(DRT::INPUT::IntegralValue<INPAR::COMBUST::VelocityJumpType>(params_->sublist("COMBUSTION FLUID"),"VELOCITY_JUMP_TYPE")),
  fluxjumptype_(DRT::INPUT::IntegralValue<INPAR::COMBUST::FluxJumpType>(params_->sublist("COMBUSTION FLUID"),"FLUX_JUMP_TYPE")),
  turbmodel_(INPAR::FLUID::no_model),
  xfemtimeint_(DRT::INPUT::IntegralValue<INPAR::COMBUST::XFEMTimeIntegration>(params_->sublist("COMBUSTION FLUID"),"XFEMTIMEINT")),
  xfemtimeint_enr_(DRT::INPUT::IntegralValue<INPAR::COMBUST::XFEMTimeIntegrationEnr>(params_->sublist("COMBUSTION FLUID"),"XFEMTIMEINT_ENR")),
  xfemtimeint_enr_comp_(DRT::INPUT::IntegralValue<INPAR::COMBUST::XFEMTimeIntegrationEnrComp>(params_->sublist("COMBUSTION FLUID"),"XFEMTIMEINT_ENR_COMP")),
  flamespeed_(params_->sublist("COMBUSTION FLUID").get<double>("LAMINAR_FLAMESPEED")),
  marksteinlength_(params_->sublist("COMBUSTION FLUID").get<double>("MARKSTEIN_LENGTH")),
  moldiffusivity_(params_->sublist("COMBUSTION FLUID").get<double>("MOL_DIFFUSIVITY")),
  nitschevel_(params_->sublist("COMBUSTION FLUID").get<double>("NITSCHE_VELOCITY")),
  nitschepres_(params_->sublist("COMBUSTION FLUID").get<double>("NITSCHE_PRESSURE")),
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
  xfemdiscret_(Teuchos::null),
  samstart_(-1),
  samstop_(-1),
  excludeXfem_(false),
  external_loads_(Teuchos::null),
  Sep_(Teuchos::null),
  fsvelafStd_(Teuchos::null),
  fsvelafXfem_(Teuchos::null),
  Cs_(0.0),
  redist_this_step_(true),
  repellant_force_(DRT::INPUT::IntegralValue<int>(params_->sublist("COMBUSTION FLUID"),"REPELLANT_FORCE"))
{
  //------------------------------------------------------------------------------------------------
  // connect degrees of freedom for periodic boundary conditions
  //------------------------------------------------------------------------------------------------
  {
    // set elements to 'standard mode' so the parallel redistribution (includes call of FillComplete)
    // hidden in the construction of periodic boundary condition runs correctly
    // remark: the 'elementdofmanager' would get lost packing and unpacking the elements
    Teuchos::ParameterList eleparams;
    eleparams.set("action","set_standard_mode");
    discret_->Evaluate(eleparams);

    // we need to keep pbc_ in order to update the PBCDofset without
    // losing the DofGid range that the PBCDofset is assigned here
    pbc_ = Teuchos::rcp(new PeriodicBoundaryConditions(discret_));
    pbc_->UpdateDofsForPeriodicBoundaryConditions();
  }

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
  physprob_.xfemfieldset_.insert(XFEM::PHYSICS::Pres);
  // for the Nitsche method assign an arbitrary element ansatz to compile
  physprob_.elementAnsatz_ = Teuchos::rcp(new COMBUST::TauPressureAnsatz());

  // create dummy instance of interfacehandle holding no flamefront and hence no integration cells
  Teuchos::RCP<COMBUST::InterfaceHandleCombust> ihdummy = Teuchos::rcp(new COMBUST::InterfaceHandleCombust(discret_,Teuchos::null));
  // create dummy instance of dof manager assigning standard enrichments to all nodes
  const Teuchos::RCP<XFEM::DofManager> dofmanagerdummy = Teuchos::rcp(new XFEM::DofManager(ihdummy,Teuchos::null,physprob_.xfemfieldset_,xparams_,Teuchos::null));

  // pass dof information to elements (no enrichments yet, standard FEM!)
  TransferDofInformationToElements(Teuchos::null, dofmanagerdummy);
  // ensure that degrees of freedom in the discretization have been set
  discret_->FillComplete();

  output_->WriteMesh(0,0.0);

  return;
}

/*------------------------------------------------------------------------------------------------*
 | Initialization                                                                 rasthofer 04/13 |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::Init()
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
    if (myrank_ == 0)
      IO::cout << "parameters 'theta' and 'time step size' have been set to 1.0 for stationary problem " << IO::endl;
  }

  numdim_ = params_->get<int>("number of velocity degrees of freedom");
  if (numdim_ != 3) dserror("COMBUST only supports 3D problems.");

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
      dserror("All surface tension approximations need a smooth gradient field of phi except laplace-beltrami and fixed-curvature! Read remark!");
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
        IO::cout << "-> location-dependent surface tension coefficient:" << IO::endl;
        IO::cout << "only coefficient linear in x available for the time being" << IO::endl;
        if (surftensapprox_ == INPAR::COMBUST::surface_tension_approx_laplacebeltrami or
            surftensapprox_ == INPAR::COMBUST::surface_tension_approx_laplacebeltrami_smoothed)
           IO::cout << "currently only connected interface possible for laplace-beltrami" << IO::endl;
      }
    }
  }
  if (combusttype_ == INPAR::COMBUST::combusttype_twophaseflowjump
      and (veljumptype_ != INPAR::COMBUST::vel_jump_none or fluxjumptype_ != INPAR::COMBUST::flux_jump_surface_tension))
  {
    if (veljumptype_ != INPAR::COMBUST::vel_jump_none)
    {
      veljumptype_ = INPAR::COMBUST::vel_jump_none;
      IO::cout << "Velocity jump is set to NONE for two-phase flow with jumps!" << IO::endl;
    }
    if (fluxjumptype_ != INPAR::COMBUST::flux_jump_surface_tension)
    {
      fluxjumptype_ = INPAR::COMBUST::flux_jump_surface_tension;
      IO::cout << "Flux jump is set to SURFACE TENSION for two-phase flow with jumps!" << IO::endl;
    }
  }

  // -------------------------------------------------------------------
  // flag for potential nonlinear boundary conditions
  // -------------------------------------------------------------------
  nonlinearbc_ = false;
  if (params_->get<std::string>("Nonlinear boundary conditions","no") == "yes")
    nonlinearbc_ = true;

  // create dummy instance of interfacehandle holding no flamefront and hence no integration cells
  Teuchos::RCP<COMBUST::InterfaceHandleCombust> ihdummy = Teuchos::rcp(new COMBUST::InterfaceHandleCombust(discret_,Teuchos::null));
  // create dummy instance of dof manager assigning standard enrichments to all nodes
  const Teuchos::RCP<XFEM::DofManager> dofmanagerdummy = Teuchos::rcp(new XFEM::DofManager(ihdummy,Teuchos::null,physprob_.xfemfieldset_,xparams_,Teuchos::null));

  // save dofmanager to be able to plot Gmsh stuff in Output()
  dofmanagerForOutput_ = dofmanagerdummy;

  // pass dof information to elements (no enrichments yet, standard FEM!)
  TransferDofInformationToElements(Teuchos::null, dofmanagerdummy);
  // ensure that degrees of freedom in the discretization have been set
  discret_->FillComplete();

  // store a dofset with the complete fluid unknowns
  standarddofset_ = Teuchos::rcp(new DRT::IndependentPBCDofSet(discret_->GetAllPBCCoupledColNodes()));
  standarddofset_->Reset();
  standarddofset_->AssignDegreesOfFreedom(*discret_,0,0);

  velpressplitterForOutput_ = Teuchos::rcp(new LINALG::MapExtractor());
  // split based on complete fluid field
  FLD::UTILS::SetupFluidSplit(*discret_,*standarddofset_,3,*velpressplitterForOutput_);

  //------------------------------------------------------------------------------------------------
  // get dof layout from the discretization to construct vectors and matrices
  //------------------------------------------------------------------------------------------------

  // parallel dof distribution contained in dofrowmap: local (LID) <-> global (GID) dof numbering
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // get layout of velocity and pressure dofs in a vector
  const int numdim = params_->get<int>("number of velocity degrees of freedom");
  velpressplitter_ = Teuchos::rcp(new LINALG::MapExtractor());
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
  Teuchos::ParameterList *  modelparams =&(params_->sublist("TURBULENCE MODEL"));

  // read turbulence model
  std::string physmodel = modelparams->get<std::string>("PHYSICAL_MODEL","no_model");
  if (physmodel != "no_model" and physmodel != "Multifractal_Subgrid_Scales")
    dserror("Turbulence model not yet supported by combustion module!");

  // -------------------------------------------------------------------
  // necessary only for the AVM3 approach:
  // fine-scale solution vector + respective output
  // -------------------------------------------------------------------

  // fine-scale subgrid viscosity for AVM3 approach
  fssgv_ = DRT::INPUT::IntegralValue<INPAR::FLUID::FineSubgridVisc>(params_->sublist("TURBULENCE MODEL"),"FSSUGRVISC");

  if (fssgv_ != INPAR::FLUID::no_fssgv)
  {
    excludeXfem_ = DRT::INPUT::IntegralValue<int>(params_->sublist("TURBULENCE MODEL"),"EXCLUDE_XFEM");
    Cs_=params_->sublist("SUBGRID VISCOSITY").get<double>("C_SMAGORINSKY") ;

    // create dofrowmap of standard dofs of initial field
    fsvelafStd_  = LINALG::CreateVector(*dofrowmap,true);
    fsvelafXfem_  = LINALG::CreateVector(*dofrowmap,true);
    // store dof distribution map of inital unenriched fluid field
    plainnodalDofDistributionMap_.clear();
    for (std::map<XFEM::DofKey, XFEM::DofGID>::const_iterator i=state_.nodalDofDistributionMap_.begin();
         i != state_.nodalDofDistributionMap_.end(); ++i)
      plainnodalDofDistributionMap_[i->first] = i->second;

    if (myrank_ == 0)
    {
      // Output
      std::cout << "Fine-scale subgrid-viscosity approach based on AVM3: ";
      std::cout << std::endl << std::endl;
      std::cout << fssgv_;
      std::cout << " with Smagorinsky constant Cs= ";
      std::cout << Cs_ ;
      std::cout << std::endl << std::endl << std::endl;
    }
  }

  // -------------------------------------------------------------------
  // preparations for multifractal subgrid-scale model
  // -------------------------------------------------------------------
  if (physmodel == "Multifractal_Subgrid_Scales")
  {
    turbmodel_ = INPAR::FLUID::multifractal_subgrid_scales;
    excludeXfem_ = DRT::INPUT::IntegralValue<int>(params_->sublist("TURBULENCE MODEL"),"EXCLUDE_XFEM");

    // create dofrowmap of standard dofs of initial field
    fsvelafStd_  = LINALG::CreateVector(*dofrowmap,true);
    fsvelafXfem_  = LINALG::CreateVector(*dofrowmap,true);
    // store dof distribution map of inital unenriched fluid field
    plainnodalDofDistributionMap_.clear();
    for (std::map<XFEM::DofKey, XFEM::DofGID>::const_iterator i=state_.nodalDofDistributionMap_.begin();
         i != state_.nodalDofDistributionMap_.end(); ++i)
      plainnodalDofDistributionMap_[i->first] = i->second;

    if (myrank_ == 0)
    {
      // Output
      std::cout << "Turbulence model        : ";
      std::cout << physmodel;
      std::cout << &std::endl;

      Teuchos::ParameterList *  mfsmodelparams =&(params_->sublist("MULTIFRACTAL SUBGRID SCALES"));
      std::cout << "                             "      ;
      std::cout << "\n";
      std::cout << "- Csgs:              " << mfsmodelparams->get<double>("CSGS") << "\n";
      std::cout << "- Scale separation:  " << mfsmodelparams->get<std::string>("SCALE_SEPARATION") << "\n";
      if ((DRT::INPUT::IntegralValue<int>(*mfsmodelparams,"CALC_N")))
      {
        std::cout << "- Re_length:         " << mfsmodelparams->get<std::string>("REF_LENGTH") << "\n";
        std::cout << "- Re_vel:            " << mfsmodelparams->get<std::string>("REF_VELOCITY") << "\n";
        std::cout << "- c_nu:              " << mfsmodelparams->get<double>("C_NU") << "\n";
      }
      else
        std::cout << "- N:                 " << mfsmodelparams->get<double>("N") << "\n";
      std::cout << "- near-wall limit:   " << DRT::INPUT::IntegralValue<int>(*mfsmodelparams,"NEAR_WALL_LIMIT") << "\n";
      std::cout << "- beta:              " << mfsmodelparams->get<double>("BETA") << "\n";
      std::cout << "- evaluation B:      " << mfsmodelparams->get<std::string>("EVALUATION_B") << "\n";
      std::cout << "- conservative:      " << mfsmodelparams->get<std::string>("CONVFORM") << "\n";
      std::cout << &std::endl;
    }
  }

  // read flag for special flow
  special_flow_ = modelparams->get<std::string>("CANONICAL_FLOW","no");
  if ( special_flow_ != "no" and
       special_flow_ != "bubbly_channel_flow" and
       special_flow_ != "combust_oracles" and
       special_flow_ != "backward_facing_step_tp")
    dserror("XFEM does not support special flow %s.", special_flow_.c_str());

  // -------------------------------------------------------------------
  // initialize turbulence statistics evaluation
  // -------------------------------------------------------------------
  turbstatisticsmanager_ = Teuchos::rcp(new FLD::TurbulenceStatisticManager(*this));

  // parameter for sampling/dumping period
  if (special_flow_ != "no")
  {
    samstart_ = modelparams->get<int>("SAMPLING_START",1);
    samstop_  = modelparams->get<int>("SAMPLING_STOP",1);
  }

  // build xfem discretization with all its special options
  xfemdiscret_ = Teuchos::rcp_dynamic_cast<DRT::DiscretizationCombust>(discret_, true);
  // build internal faces, if required
  if(params_->sublist("RESIDUAL-BASED STABILIZATION").get<std::string>("STABTYPE")=="edge_based" or
     DRT::INPUT::IntegralValue<bool>(params_->sublist("COMBUSTION FLUID"),"XFEMSTABILIZATION") == true)
  {
    // if the definition of internal faces would be included
    // in the standard discretization, these lines can be removed
    // and CreateInternalFacesExtension() can be called once
    // in the constructor of the fluid time integration
    // since we want to keep the standard discretization as clean as
    // possible, we create internal faces via an enhanced discretization
    // including the faces between elements
    xfemdiscret_->CreateInternalFacesExtension(true);
  }

  if (xfemtimeint_ == INPAR::COMBUST::xfemtimeint_semilagrange)
  {
      // get periodic surface boundary conditions
      std::vector<DRT::Condition*> mysurfpbcs;
      discret_->GetCondition("SurfacePeriodic",mysurfpbcs);
      if(mysurfpbcs.empty())
      {}
      else
      {
        if (myrank_==0)
          std::cout << "WARNING: Semi-Lagrange time integration does not entirely account for PBCS" << std::endl;
      }
  }


#ifdef FLAME_VORTEX
  if (myrank_ == 0)
  {
    //---------------------------------
    // open file for flame surface area
    //---------------------------------
    std::ostringstream tmpfilename;
    const std::string filebase(DRT::Problem::Instance()->OutputControlFile()->FileName());

    tmpfilename << filebase << "." << "flame_area" << ".txt";
    IO::cout << "writing " << left << std::setw(60) <<tmpfilename.str()<<"...";

    const std::string filename = tmpfilename.str();
    std::ofstream gmshfilecontent(filename.c_str());
    gmshfilecontent << "# total flame area for each time step\n";
    gmshfilecontent << "#" << std::setw(6) << "step" << "  " << std::setw(6) << std::setprecision(6) << "time" << "  " << std::setprecision(13) << "flame area" << "\n\n";
  }
#endif

#ifdef DL_INSTAB
  if (myrank_ == 0)
  {
    //---------------------------------
    // open file for flame surface area
    //---------------------------------
    std::ostringstream tmpfilename;
    const std::string filebase(DRT::Problem::Instance()->OutputControlFile()->FileName());

    tmpfilename << filebase << "." << "amplitude" << ".txt";
    IO::cout << "writing " << left << std::setw(60) <<tmpfilename.str()<<"...";

    const std::string filename = tmpfilename.str();
    std::ofstream gmshfilecontent(filename.c_str());
    gmshfilecontent << "# amplitude on middle line for each time step\n";
    gmshfilecontent << "#" << std::setw(6) << "step" << "  " << std::setw(6) << std::setprecision(6) << "time";
    gmshfilecontent << "  " << std::setprecision(13) << "amplitude left";
    gmshfilecontent << "  " << std::setprecision(13) << "amplitude middle";
    gmshfilecontent << "  " << std::setprecision(13) << "amplitude" << "\n\n";
  }
#endif
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
    {
      dserror("enrichment recomputation approach in XFEM time integration not implemented");
      break;
    }
    }

    switch (xfemtimeint_)
    {
    case INPAR::COMBUST::xfemtimeint_donothing:
    case INPAR::COMBUST::xfemtimeint_extrapolationold:
    case INPAR::COMBUST::xfemtimeint_extrapolationnew:
    {
      finished = true; // the above standard value computations don't require -> iterations
      break;
    }
    case INPAR::COMBUST::xfemtimeint_semilagrange:
    case INPAR::COMBUST::xfemtimeint_mixedSLExtrapol:
    case INPAR::COMBUST::xfemtimeint_mixedSLExtrapolNew:
      break;
    default:
    {
      dserror("standard recomputation approach in XFEM time integration not implemented");
      break;
    }
    }
  }

  // if not finished, print iteration counter, update false
  if (!finished)
  {
    totalitnumFRS_++;
    curritnumFRS_++;
    if (myrank_==0 and itemaxFRS_!=1)
      IO::cout << "FLUID REFERENCE SOLUTION ITERATION " << curritnumFRS_ << "/" << itemaxFRS_ << IO::endl;
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
    IO::cout << "/!\\ warning: 'time' and 'time step' are set to 1.0 and 1.0 for output control file" << IO::endl;
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
      IO::cout << "/!\\ warning: initial solution is computed by stationary algorithm" << IO::endl;
      IO::cout << "/!\\ warning: 'time' and 'time step' are set to 0.0 and 1.0 for output control file" << IO::endl;
      timealgo_ = INPAR::FLUID::timeint_stationary;
      time_ =  0.0; // only needed for output
      dta_ =   1.0; // for calculation, we reset this value at the end of Solve()
      dtp_ =   1.0; // for calculation, we reset this value at the end of Solve()
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
      if (myrank_ == 0)
        IO::cout << "/!\\ first time step computed with theta =  " << theta_ << IO::endl;
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

        IO::cout << "/!\\ first " << startsteps_ << " steps are computed with Backward-Euler scheme" << IO::endl;
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
      if (myrank_ == 0)
        IO::cout << "/!\\ first time step computed with theta =  " << theta_ << IO::endl;
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
  {
    dserror("unknown time integration scheme");
    break;
  }
  }
}

/*------------------------------------------------------------------------------------------------*
 | prepare a fluid nonlinear iteration                                                henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::PrepareSolve()
{
  // -------------------------------------------------------------------
  //                     do explicit predictor step
  // -------------------------------------------------------------------
  //
//  {
//    // currently not available
//    if (step_>1)
//    {
//      double timealgo_constant=theta_;
//
//      UTILS::ExplicitPredictor(
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
    Teuchos::ParameterList eleparams;

    // other parameters needed by the elements
    eleparams.set("total time",time_);
    eleparams.set("delta time",dta_);
    eleparams.set("thsl",theta_*dta_);


    if (step_ > 0) // initial stationary solution does not require reference solution update
    {
      TEUCHOS_FUNC_TIME_MONITOR("   + xfem time integration");

      std::vector<Teuchos::RCP<Epetra_Vector> > newRowVectorsn; // solution vectors of old time step due to new interface
      newRowVectorsn.push_back(state_.veln_);
      newRowVectorsn.push_back(state_.accn_);

      std::vector<Teuchos::RCP<Epetra_Vector> > newRowVectorsnp; // solution vectors at new time step
      newRowVectorsnp.push_back(state_.velnp_);
      newRowVectorsnp.push_back(state_.accnp_);

      if(xfemtimeint_ != INPAR::COMBUST::xfemtimeint_donothing) // do something for standard values
      {
        if (myrank_ == 0)
          IO::cout << "---  XFEM time integration: adopt standard values to new interface... ";// << IO::flush;
        timeIntStd_->type(totalitnumFRS_,itemaxFRS_); // update algorithm handling
        timeIntStd_->compute(newRowVectorsn,newRowVectorsnp); // call computation
        if (myrank_ == 0)
          IO::cout << "done" << IO::endl;
      }

      if ((xfemtimeint_enr_!=INPAR::COMBUST::xfemtimeintenr_donothing) and
          (xfemtimeint_enr_!=INPAR::COMBUST::xfemtimeintenr_quasistatic)) // do something for enrichment values
      {
        if (myrank_ == 0)
          IO::cout << "---  XFEM time integration: adopt enrichment values to new interface... ";// << IO::flush;
        timeIntEnr_->type(totalitnumFRS_,itemaxFRS_); // update algorithm handling
        timeIntEnr_->compute(newRowVectorsn,newRowVectorsnp); // call computation
        if (myrank_ == 0)
          IO::cout << "done" << IO::endl;
      }
    }

    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("velnp",state_.velnp_);
    // predicted dirichlet values
    // velnp then also holds prescribed new dirichlet values
    xfemdiscret_->EvaluateDirichletCombust(eleparams,state_.velnp_,Teuchos::null,Teuchos::null,Teuchos::null,dbcmaps_);

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

  // -------------------------------------------------------------------
  // create vectors for Krylov projection if necessary
  // -------------------------------------------------------------------

  // sysmat might be singular (if we have a purely Dirichlet constrained
  // problem, the pressure mode is defined only up to a constant)
  // in this case, we need a basis vector for the nullspace/kernel

  // get condition "KrylovSpaceProjection" from discretization
  std::vector<DRT::Condition*> KSPcond;
  discret_->GetCondition("KrylovSpaceProjection",KSPcond);
  int numcond = KSPcond.size();
  int numfluid = 0;

  DRT::Condition* kspcond = NULL;
  // check if for fluid Krylov projection is required
  for(int icond = 0; icond < numcond; icond++)
  {
    const std::string* name = KSPcond[icond]->Get<std::string>("discretization");
    if (*name == "fluid")
    {
      numfluid++;
      kspcond = KSPcond[icond];
    }
  }

  // initialize variables for Krylov projection if necessary
  if (numfluid == 1)
  {
    SetupKrylovSpaceProjection(kspcond);
    if (myrank_ == 0 and step_ == 1)
    {
      IO::cout << "\nSetup of KrylovSpaceProjection in fluid field\n";
      IO::cout << "#################################################\n";
      IO::cout << "#        WARNING !!!                            #\n";
      IO::cout << "# Krylov-Projection for combustion problems not #\n";
      IO::cout << "# carefully checked!                            #\n";
      IO::cout << "# Read remark!                                  #\n";
      IO::cout << "#################################################" << IO::endl;
      // remark:
      // before using this option, the following points should be carefully checked!
      // 1.) for single-phase flows, combust and fluid should provide the same (!) results, see beltrami test
      //     -> ensured by two nightly tests
      // 2.) for two-phase flows, the krylov projection should also work:
      //     - currently, it does not: this might be traced back to the under-integration of intersected elements,
      //       since this alters the system matrix and hence its kernel
      //     - an appropriate test to check this issue would be the two-phase couette flow with hexahedral integration cells,
      //       since, for this test case, we can ensure an exact integration also in intersected elements
      //     - or channel flow with interface: this is already considered by one of the tests and works fine if the pressure is
      //       not enriched
    }
  }
  else if (numfluid == 0)
  {
    updateprojection_ = false;
    projector_ = Teuchos::null;
  }
  else
    dserror("Received more than one KrylovSpaceCondition for fluid field");

  // -------------------------------------------------------------------
  //                     prepare external load
  // -------------------------------------------------------------------
  // evaluate node-based forces
  // (currently only required for turbulent bubbly channel flow)
  ComputeExternalForces();

  return;
}

/*------------------------------------------------------------------------------------------------*
 | hand over information about (XFEM) degrees of freedom to elements                     ag 04/09 |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::TransferDofInformationToElements(
    const Teuchos::RCP<COMBUST::FlameFront> flamefront,
    const Teuchos::RCP<XFEM::DofManager> dofmanager
    )
{
  Teuchos::ParameterList eleparams;
  eleparams.set("action","store_xfem_info");
  eleparams.set("dofmanager",dofmanager);
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
  const Teuchos::RCP<XFEM::DofManager> dofmanager = Teuchos::rcp(new XFEM::DofManager(
      interfacehandle_,
      phinp_,
      physprob_.xfemfieldset_,
      xparams_,
      discret_->GetAllPBCCoupledColNodes())
  );

  // temporarely save old dofmanager
  const Teuchos::RCP<XFEM::DofManager> olddofmanager = dofmanagerForOutput_;

  // save dofmanager to be able to plot Gmsh stuff in Output()
  dofmanagerForOutput_ = dofmanager;

  // print global and element dofmanager to Gmsh
  dofmanager->toGmsh(step_);
  interfacehandle_->toGmsh(step_);

  // get old dofmaps, compute a new one and get the new one, too
  const Epetra_Map olddofrowmap = *discret_->DofRowMap();
  const Epetra_Map olddofcolmap = *discret_->DofColMap();
  std::map<XFEM::DofKey,XFEM::DofGID> oldNodalDofColDistrib;
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
    std::vector<Teuchos::RCP<Epetra_Vector> > oldColStateVectors; // same order of vectors as newRowVectorsn combined with newRowVectorsnp
    Teuchos::RCP<Epetra_Vector> veln = Teuchos::rcp(new Epetra_Vector(olddofcolmap,true));
    {
      LINALG::Export(*state_.veln_,*veln);
      oldColStateVectors.push_back(veln);

      Teuchos::RCP<Epetra_Vector> accn = Teuchos::rcp(new Epetra_Vector(olddofcolmap,true));
      LINALG::Export(*state_.accn_,*accn);
      oldColStateVectors.push_back(accn);

      Teuchos::RCP<Epetra_Vector> velnp = Teuchos::rcp(new Epetra_Vector(olddofcolmap,true));
      LINALG::Export(*state_.velnp_,*velnp);
      oldColStateVectors.push_back(velnp);

      Teuchos::RCP<Epetra_Vector> accnp = Teuchos::rcp(new Epetra_Vector(olddofcolmap,true));
      LINALG::Export(*state_.accnp_,*accnp);
      oldColStateVectors.push_back(accnp);
    }

    //---------------------------------------------
    // switch state vectors to new dof distribution
    //---------------------------------------------
    if (myrank_ == 0)
      IO::cout << "---  transform state vectors... ";// << IO::flush;
    // quasi-static enrichment strategy for kink enrichments
    // remark: as soon as the XFEM-time-integration works for kinks, this should be removed
    if (xfemtimeint_enr_==INPAR::COMBUST::xfemtimeintenr_quasistatic)
    {
      if (myrank_ == 0)
        IO::cout << "quasi-static enrichment for two-phase flow problems... ";// << IO::flush;

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
    if (myrank_ == 0)
      IO::cout << "done" << IO::endl;

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
#ifdef ORACLES
        const double dens_minus = 1.296;
        const double dens_plus = 0.166;
#endif
#ifdef FLAME_VORTEX
        const double dens_minus = 1.161;
        const double dens_plus = 0.157;
#else
        const double dens_minus = 1.0;
        const double dens_plus = 1.0;
#endif
        const double mflux = flamespeed_*dens_minus;
        const double veljump = -mflux*(1.0/dens_minus - 1.0/dens_plus);
        if ((xfemtimeint_!=INPAR::COMBUST::xfemtimeint_donothing) or
            ((xfemtimeint_enr_ !=INPAR::COMBUST::xfemtimeintenr_donothing) and (xfemtimeint_enr_ !=INPAR::COMBUST::xfemtimeintenr_quasistatic)))
        {
          // basic time integration data
          Teuchos::RCP<XFEM::TIMEINT> timeIntData = Teuchos::rcp(new XFEM::TIMEINT(
              discret_,
              olddofmanager,
              dofmanager,
              oldColStateVectors,
              flamefront,
              olddofcolmap,
              newdofrowmap,
              oldNodalDofColDistrib,
              state_.nodalDofDistributionMap_,
              discret_->GetAllPBCCoupledColNodes()));

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
            timeIntEnr_ = Teuchos::rcp(new XFEM::EnrichmentProjection(
                *timeIntData,
                veljump,
                xfemtimeint_enr_,
                timeIntEnrComp));
            break;
          }
          default:
          {
            dserror("enrichment recomputation approach in XFEM time integration not implemented");
            break;
          }
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
            timeIntStd_ = Teuchos::rcp(new XFEM::SemiLagrange(
                *timeIntData,
                xfemtimeint_,
                veln,
                dta_,
                theta_,
                flamefront,
                veljump,
                true));
            break;
          }
          case INPAR::COMBUST::xfemtimeint_extrapolationold:
          {
            // time integration data for standard dofs, extrapolation approach
            timeIntStd_ = Teuchos::rcp(new XFEM::ExtrapolationOld(
                *timeIntData,
                xfemtimeint_,
                veln,
                dta_,
                flamefront,
                veljump,
                true));
            break;
          }
          case INPAR::COMBUST::xfemtimeint_extrapolationnew:
          {
            // time integration data for standard dofs, extrapolation approach
            timeIntStd_ = Teuchos::rcp(new XFEM::ExtrapolationNew(
                *timeIntData,
                xfemtimeint_,
                veln,
                dta_,
                flamefront,
                true));
            break;
          }
          default:
          {
            dserror("standard recomputation approach in XFEM time integration not implemented");
            break;
          }
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
        default:
        {
          dserror("enrichment recomputation approach in XFEM time integration not implemented");
          break;
        }
        }

        switch (xfemtimeint_)
        {
        case INPAR::COMBUST::xfemtimeint_donothing: break;
        case INPAR::COMBUST::xfemtimeint_semilagrange:
        case INPAR::COMBUST::xfemtimeint_mixedSLExtrapol:
        case INPAR::COMBUST::xfemtimeint_mixedSLExtrapolNew:
        case INPAR::COMBUST::xfemtimeint_extrapolationold:
        case INPAR::COMBUST::xfemtimeint_extrapolationnew:
        {
          timeIntStd_->importNewFGIData(
              discret_,
              dofmanager,
              flamefront,
              newdofrowmap,
              state_.nodalDofDistributionMap_); // new data due to gamma^n+1,i+1
          break;
        }
        default:
        {
          dserror("standard recomputation approach in XFEM time integration not implemented");
          break;
        }
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
    Teuchos::ParameterList eleparams;
    // other parameters needed by the elements
    eleparams.set("total time",time_);

    xfemdiscret_->EvaluateDirichletCombust(eleparams, zeros_, Teuchos::null, Teuchos::null, Teuchos::null, dbcmaps_);
    zeros_->PutScalar(0.0); // just in case of change
  }

  neumann_loads_= LINALG::CreateVector(newdofrowmap,true);

  // ---------------------------------
  // Vectors used for solution process
  // ---------------------------------
  residual_     = LINALG::CreateVector(newdofrowmap,true);
  trueresidual_ = LINALG::CreateVector(newdofrowmap,true);
  incvel_       = LINALG::CreateVector(newdofrowmap,true);

  // -----------------------------------------------
  // fine-scale velocity vector based on new dofset
  // -----------------------------------------------

  if (fssgv_ != INPAR::FLUID::no_fssgv or turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales)
    fsvelafXfem_= LINALG::CreateVector(newdofrowmap,true);

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

  return;
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

    // initial call (UpdateDofSet==false):
    // store flame front without extending the dofset to XFEM-layout
    // remark: this is required to set the initial fluid field relying on an inital interface
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

  //TODO es ist sowas auch noch in PrepareSolve(). Was brauchen wir?
  //stationary case (timealgo_== INPAR::FLUID::timeint_stationary))
  if ( (timealgo_==INPAR::FLUID::timeint_one_step_theta) or
       (timealgo_==INPAR::FLUID::timeint_afgenalpha) )
    COMBUST::UTILS::SetOldPartOfRighthandside(veln,Teuchos::null, accn,timealgo_, dta_, theta_, hist);
  else
    dserror("time integration scheme not supported");

  return hist;
}

/*------------------------------------------------------------------------------------------------*
 | solve the nonlinear fluid problem                                                  henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::Solve()
{
  while (!FluidRefSolLoopFinished()) // iterator between NS-solution and recomputation of reference solution
  {
    PrepareSolve();

    TEUCHOS_FUNC_TIME_MONITOR("   + nonlin. iteration/lin. solve");

    // ---------------------------------------------- nonlinear iteration
    // ------------------------------- stop nonlinear iteration when both
    //                                 increment-norms are below this bound
    const double  ittol     = params_->get<double>("tolerance for nonlin iter");

    //------------------------------ turn adaptive solver tolerance on/off
    const bool   isadapttol    = params_->get<bool>("ADAPTCONV");
    const double adaptolbetter = params_->get<double>("ADAPTCONV_BETTER");

    int               itnum = 0;
    bool              stopnonliniter = false;

    double dtsolve = 0.0;
    dtele_   = 0.0;

    // out to screen
    PrintTimeStepInfo();

    if (myrank_ == 0)
    {
      IO::cout << "+------------+-------------------+--------------+--------------+--------------+--------------+" << IO::endl
               << "|- step/max -|- tol      [norm] -|-- vel-res ---|-- pre-res ---|-- vel-inc ---|-- pre-inc ---|" << IO::endl;
    }

    incvel_->PutScalar(0.0);
    residual_->PutScalar(0.0);

    while (stopnonliniter==false)
    {
      itnum++;

#ifdef SUGRVEL_OUTPUT
      IO::cout << "writing gmsh output" << IO::endl;
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

        // add external loads
        if(external_loads_ != Teuchos::null)
          residual_->Update(1.0/ResidualScaling(),*external_loads_,1.0);

        // create the parameters for the discretization
        Teuchos::ParameterList eleparams;

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

        // parameter for suppressing additional enrichment dofs in two-phase flow problems
        eleparams.set<int>("selectedenrichment",DRT::INPUT::IntegralValue<INPAR::COMBUST::SelectedEnrichment>(params_->sublist("COMBUSTION FLUID"),"SELECTED_ENRICHMENT"));

        // parameters for two-phase flow problems with surface tension
        eleparams.set<int>("surftensapprox",surftensapprox_);
        eleparams.set<double>("variablesurftens",params_->sublist("COMBUSTION FLUID").get<double>("VARIABLESURFTENS"));
        eleparams.set<bool>("connected_interface",connected_interface_);
        eleparams.set<bool>("second_deriv",DRT::INPUT::IntegralValue<int>(params_->sublist("COMBUSTION FLUID"),"L2_PROJECTION_SECOND_DERIVATIVES"));

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


        // additional terms for Nitsche's method (see Diss Florian)
        eleparams.set<bool>("nitsche_convflux",DRT::INPUT::IntegralValue<int>(params_->sublist("COMBUSTION FLUID"),"NITSCHE_CONVFLUX"));
        eleparams.set<bool>("nitsche_convstab",DRT::INPUT::IntegralValue<int>(params_->sublist("COMBUSTION FLUID"),"NITSCHE_CONVSTAB"));
        eleparams.set<bool>("nitsche_convpenalty",DRT::INPUT::IntegralValue<int>(params_->sublist("COMBUSTION FLUID"),"NITSCHE_CONVPENALTY"));
        // further terms and parameters related to Nitsche's method
        eleparams.set<bool>("nitsche_mass",DRT::INPUT::IntegralValue<int>(params_->sublist("COMBUSTION FLUID"),"NITSCHE_MASS"));
        eleparams.set<INPAR::COMBUST::WeightType>("weighttype",DRT::INPUT::IntegralValue<INPAR::COMBUST::WeightType>(params_->sublist("COMBUSTION FLUID"),"NITSCHE_WEIGHT"));


#ifdef SUGRVEL_OUTPUT
        //eleparams.set("step",step_);
#endif

        //eleparams.set("include reactive terms for linearisation",params_->get<bool>("Use reaction terms for linearisation"));
        //type of linearisation: include reactive terms for linearisation
        if(DRT::INPUT::get<INPAR::FLUID::LinearisationAction>(*params_, "Linearisation") == INPAR::FLUID::Newton)
          eleparams.set("include reactive terms for linearisation",true);
        else
          eleparams.set("include reactive terms for linearisation",false);

        // parameters for stabilization
        eleparams.sublist("RESIDUAL-BASED STABILIZATION") = params_->sublist("RESIDUAL-BASED STABILIZATION");

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

        // parameters for turbulence
        eleparams.set<INPAR::FLUID::TurbModelAction>("turbmodel",turbmodel_);
        eleparams.set<INPAR::FLUID::FineSubgridVisc>("fssgv",fssgv_);
        eleparams.set<double>("Cs",Cs_);

        //----------------------------------------------------------------------
        // decide whether AVM3-based solution approach or standard approach
        //----------------------------------------------------------------------
        if (fssgv_ != INPAR::FLUID::no_fssgv or turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales)
        {
          AVM3Separation();

          if (turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales)
            eleparams.sublist("MULTIFRACTAL SUBGRID SCALES") = params_->sublist("MULTIFRACTAL SUBGRID SCALES");
        }

        if ((itnum != itemax_))
        {
          // call standard loop over elements
          discret_->Evaluate(eleparams,sysmat_,residual_);

          discret_->ClearState();

          // add edged-based stabilization, if selected
          if(params_->sublist("RESIDUAL-BASED STABILIZATION").get<std::string>("STABTYPE")=="edge_based" or
             DRT::INPUT::IntegralValue<bool>(params_->sublist("COMBUSTION FLUID"),"XFEMSTABILIZATION") == true)
          {
            // create the parameters for the discretization
            Teuchos::ParameterList faceeleparams;

            // action for elements
            faceeleparams.set("action","calc_edge_based_stab_terms");

            // combustion problem
            faceeleparams.set<int>("combusttype",combusttype_);

            // set xfem stabilization, i.e., ghost penalty
            faceeleparams.set<bool>("xfemstab",DRT::INPUT::IntegralValue<INPAR::COMBUST::SmoothGradPhi>(params_->sublist("COMBUSTION FLUID"),"XFEMSTABILIZATION"));

            // other parameters that might be needed by the elements
            faceeleparams.set<int>("timealgo",timealgo_);
            faceeleparams.set<double>("time",time_);
            faceeleparams.set<double>("dt",dta_);
            faceeleparams.set<double>("theta",theta_);
            faceeleparams.set<double>("gamma",gamma_);
            faceeleparams.set<double>("alphaF",alphaF_);
            faceeleparams.set<double>("alphaM",alphaM_);

            // parameter for suppressing additional enrichment dofs in two-phase flow problems
            faceeleparams.set<int>("selectedenrichment",DRT::INPUT::IntegralValue<INPAR::COMBUST::SelectedEnrichment>(params_->sublist("COMBUSTION FLUID"),"SELECTED_ENRICHMENT"));

            // set phinp
            faceeleparams.set<Teuchos::RCP<Epetra_Vector> >("phinp",phinp_);
            // interface handle
            faceeleparams.set<Teuchos::RCP<COMBUST::InterfaceHandleCombust> >("interface handle",interfacehandle_);

            // parameters for stabilization
            faceeleparams.sublist("RESIDUAL-BASED STABILIZATION") = params_->sublist("RESIDUAL-BASED STABILIZATION");
            faceeleparams.sublist("EDGE-BASED STABILIZATION") = params_->sublist("EDGE-BASED STABILIZATION");

            // set scheme-specific element parameters and vector values
            if (timealgo_==INPAR::FLUID::timeint_afgenalpha)
              discret_->SetState("velaf",state_.velaf_);
            else
            discret_->SetState("velnp",state_.velnp_);

            xfemdiscret_->EvaluateEdgeBasedCombust(faceeleparams,sysmat_,residual_);

            discret_->ClearState();
          }


          // potential Neumann inflow terms
          if (nonlinearbc_)
          {
            // create parameter list
            Teuchos::ParameterList condparams;

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
            discret_->EvaluateCondition(condparams,sysmat_,Teuchos::null,residual_,Teuchos::null,Teuchos::null,condstring);
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
      if (updateprojection_)
      {
        UpdateKrylovSpaceProjection();
      }
      // remove contributions of pressure mode
      // that would not vanish due to the projection
      if (projector_ != Teuchos::null)
        projector_->ApplyPT(*residual_);

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
        IO::cout << "|  " << std::setw(3) << itnum << "/" << std::setw(3) << itemax_ << "   | "
                 << std::setw(10) << std::setprecision(3) << std::scientific << ittol << "[L_2 ]  | "
                 << std::setw(10) << std::setprecision(3) << std::scientific << vresnorm << "   | "
                 << std::setw(10) << std::setprecision(3) << std::scientific << presnorm << "   |      --      |      --      | (      --     ,te_min="
                 << std::setw(10) << std::setprecision(3) << std::scientific << min << ",te_max="
                 << std::setw(10) << std::setprecision(3) << std::scientific << max << ")" << IO::endl;
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
          IO::cout << "|  " << std::setw(3) << itnum << "/" << std::setw(3) << itemax_ << "   | "
                   << std::setw(10) << std::setprecision(3) << std::scientific << ittol << "[L_2 ]  | "
                   << std::setw(10) << std::setprecision(3) << std::scientific << vresnorm << "   | "
                   << std::setw(10) << std::setprecision(3) << std::scientific << presnorm << "   | "
                   << std::setw(10) << std::setprecision(3) << std::scientific << incvelnorm_L2/velnorm_L2 << "   | "
                   << std::setw(10) << std::setprecision(3) << std::scientific << incprenorm_L2/prenorm_L2 << "   | (ts="
                   << std::setw(10) << std::setprecision(3) << std::scientific << dtsolve << ",te_min="
                   << std::setw(10) << std::setprecision(3) << std::scientific << min << ",te_max="
                   << std::setw(10) << std::setprecision(3) << std::scientific << max << ")\n"
                   << "+------------+-------------------+--------------+--------------+--------------+--------------+" << IO::endl;

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
          std::map<XFEM::DofKey,XFEM::DofGID> dofColDistrib;
          dofmanagerForOutput_->fillNodalDofColDistributionMap(dofColDistrib);

          Teuchos::RCP<Epetra_Vector> velnp = Teuchos::rcp(new Epetra_Vector(*dofcolmap,true));
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
                IO::cout << *frontnode << IO::endl;
                IO::cout << *backnode << IO::endl;
                dserror("wrong order of nodes as thought here!");
              }

              // compare both fieldenrsets
              const std::set<XFEM::FieldEnr>& backfieldEnrSet(dofmanagerForOutput_->getNodeDofSet(backnode->Id()));
              int i=0;
              for (std::set<XFEM::FieldEnr>::const_iterator frontfieldenr = frontfieldEnrSet.begin();
                  frontfieldenr != frontfieldEnrSet.end();frontfieldenr++)
              {
                int j=0;
                for (std::set<XFEM::FieldEnr>::const_iterator backfieldenr = backfieldEnrSet.begin();
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

            for (std::set<XFEM::FieldEnr>::const_iterator fieldenr = frontfieldEnrSet.begin();
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
          IO::cout << "|  " << std::setw(3) << itnum << "/" << std::setw(3) << itemax_ << "   | "
                   << std::setw(10) << std::setprecision(3) << std::scientific << ittol << "[L_2 ]  | "
                   << std::setw(10) << std::setprecision(3) << std::scientific << vresnorm << "   | "
                   << std::setw(10) << std::setprecision(3) << std::scientific << presnorm << "   | "
                   << std::setw(10) << std::setprecision(3) << std::scientific << incvelnorm_L2/velnorm_L2 << "   | "
                   << std::setw(10) << std::setprecision(3) << std::scientific << incprenorm_L2/prenorm_L2 << "   | (ts="
                   << std::setw(10) << std::setprecision(3) << std::scientific << dtsolve << ",te_min="
                   << std::setw(10) << std::setprecision(3) << std::scientific << min << ",te_max="
                   << std::setw(10) << std::setprecision(3) << std::scientific << max << ")" << IO::endl;
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
        IO::cout << "+---------------------------------------------------------------+\n"
                 << "|            >>>>>> not converged in itemax steps!              |\n"
                 << "+---------------------------------------------------------------+" << IO::endl;

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
        std::map<XFEM::DofKey,XFEM::DofGID> dofColDistrib;
        dofmanagerForOutput_->fillNodalDofColDistributionMap(dofColDistrib);

        Teuchos::RCP<Epetra_Vector> velnp = Teuchos::rcp(new Epetra_Vector(*dofcolmap,true));
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
              IO::cout << *frontnode << IO::endl;
              IO::cout << *backnode << IO::endl;
              dserror("wrong order of nodes as thought here!");
            }

            // compare both fieldenrsets
            const std::set<XFEM::FieldEnr>& backfieldEnrSet(dofmanagerForOutput_->getNodeDofSet(backnode->Id()));
            int i=0;
            for (std::set<XFEM::FieldEnr>::const_iterator frontfieldenr = frontfieldEnrSet.begin();
                frontfieldenr != frontfieldEnrSet.end();frontfieldenr++)
            {
              int j=0;
              for (std::set<XFEM::FieldEnr>::const_iterator backfieldenr = backfieldEnrSet.begin();
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

          for (std::set<XFEM::FieldEnr>::const_iterator fieldenr = frontfieldEnrSet.begin();
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
        double currresidual = std::max(vresnorm,presnorm);
        currresidual = std::max(currresidual,incvelnorm_L2/velnorm_L2);
        currresidual = std::max(currresidual,incprenorm_L2/prenorm_L2);
        solver_->AdaptTolerance(ittol,currresidual,adaptolbetter);
      }

//      // print matrix in matlab format
//      // matrix printing options (DEBUGGING!)
//      IO::cout << "print matrix in matlab format to sparsematrix.mtl" << IO::endl;
//      //cast sysmat to SparseMatrix
//      Teuchos::RCP<LINALG::SparseMatrix> A = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(sysmat_);
//      if (A != Teuchos::null)
//      {
//        const std::string fname = "sparsematrix.mtl";
//        IO::cout << "Col:     MinGID " << A->ColMap().MinAllGID() << " MaxGID "   <<   A->ColMap().MaxAllGID()       <<    "   Row:    MinGID    " <<   A->RowMap().MinAllGID() << " MaxGID " << A->RowMap().MaxAllGID()<< IO::endl;
//        // print to  le in matlab format + const std::string fname = "sparsematrix.mtl";
//        LINALG::PrintMatrixInMatlabFormat(fname,*(A->EpetraMatrix()));
//        // print to screen
//        // (A->EpetraMatrix())->Print(std::cout);
//        // print sparsity pattern to  le
//        //LINALG::PrintSparsityToPostscript( *(A->EpetraMatrix()) );
//      }
//      else
//      {
//        IO::cout << "failure" << IO::endl;
//        //   Teuchos::RCP<LINALG::BlockSparseMatrixBase>            A  =   +//  state_->BlockSystemMatrix();
//        // LINALG::PrintBlockMatrixInMatlabFormat(fname,*(A));
//      }
//      IO::cout << " ...done" << IO::endl;
//      dserror("Fertig");

      // if Krylov space projection is used, check whether constant pressure
      // is in nullspace of sysmat_
      CheckMatrixNullspace();

      solver_->Solve(sysmat_->EpetraOperator(),incvel_,residual_,true,itnum==1,projector_); // defaults: projector_ = Teuchos::null
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

  incvel_ = Teuchos::null;
  residual_ = Teuchos::null;
  zeros_ = Teuchos::null;
  sysmat_ = Teuchos::null;
}


/*--------------------------------------------------------------------------*
 | setup Krylov projector including first fill                    nis Feb13 |
 *--------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::SetupKrylovSpaceProjection(DRT::Condition* kspcond)
{
  // the Krylov space projection for a combust fluid is questionable for more
  // than one phase, since pressure splitters etc do not seem to account for
  // the different fluids.

  // confirm that mode flags are number of nodal dofs
  const int nummodes = kspcond->GetInt("NUMMODES");
  if (nummodes!=(numdim_+1))
    dserror("Expecting numdim_+1 modes in Krylov projection definition. Check dat-file!");

  // get vector of mode flags as given in dat-file
  const std::vector<int>* modeflags = kspcond->Get<std::vector<int> >("ONOFF");

  // confirm that only the pressure mode is selected for Krylov projection in dat-file
  for(int rr=0;rr<numdim_;++rr)
  {
    if(((*modeflags)[rr])!=0)
    {
      dserror("Expecting only an undetermined pressure. Check dat-file!");
    }
  }
  if(((*modeflags)[numdim_])!=1)
    dserror("Expecting an undetermined pressure. Check dat-file!");
  std::vector<int> activemodeids(1,numdim_);

  // allocate kspsplitter_
  kspsplitter_ = Teuchos::rcp(new FLD::UTILS::KSPMapExtractor());
  // create map of nodes involved in Krylov projection
  kspsplitter_->Setup(*discret_);

  // get from dat-file definition how weights are to be computed
  const std::string* weighttype = kspcond->Get<std::string>("weight vector definition");

  // set flag for projection update true
  // (safe option - normally only necessary for moving interfaces)
  updateprojection_ = true;

  // create the projector
  projector_ = Teuchos::rcp(new LINALG::KrylovProjector(activemodeids,weighttype,discret_->DofRowMap()));

  // update the projector
  UpdateKrylovSpaceProjection();

  return;
} // FLD::CombustFluidImplicitTimeInt::SetupKrylovSpaceProjection


/*--------------------------------------------------------------------------*
 | update projection vectors w_ and c_ for Krylov projection      nis Feb13 |
 *--------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::UpdateKrylovSpaceProjection()
{
  // get Teuchos::RCP to kernel vector of projector
  Teuchos::RCP<Epetra_MultiVector> c = projector_->GetNonConstKernel();
  Teuchos::RCP<Epetra_Vector> c0 = Teuchos::rcp((*c)(0),false);
  c0->PutScalar(0.0);

  // extract vector of pressure-dofs
  Teuchos::RCP<Epetra_Vector> presmode = velpressplitter_->ExtractCondVector(*((*c)(0)));

  const std::string* weighttype = projector_->WeightType();
  // compute w_ as defined in dat-file
  if(*weighttype == "pointvalues")
  {
    /*
    // export to vector to normalize against
    // Note that in the case of definition pointvalue based,
    // the average pressure will vanish in a pointwise sense
    //
    //    +---+
    //     \
    //      +   p_i  = 0
    //     /
    //    +---+
    //
    // (everything is done below)
    */
  }
  else if(*weighttype == "integration")
  {
    // get Teuchos::RCP to weight vector of projector
    Teuchos::RCP<Epetra_MultiVector> w = projector_->GetNonConstWeights();
    Teuchos::RCP<Epetra_Vector> w0 = Teuchos::rcp((*w)(0),false);
    w0->PutScalar(0.0);

    // create parameter list for condition evaluate and ...
    Teuchos::ParameterList mode_params;
    // ... set action for elements to integration of shape functions
    mode_params.set("action","integrate_Shapefunction");

    /*
    // evaluate KrylovSpaceProjection condition in order to get
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

    // compute w_ by evaluating the integrals of all pressure basis functions
    discret_->EvaluateCondition
      (mode_params        ,
       Teuchos::null      ,
       Teuchos::null      ,
       w0                 ,
       Teuchos::null      ,
       Teuchos::null      ,
       "KrylovSpaceProjection");
    //w0->Print(std::cout);

  }
  else
  {
    dserror("unknown definition of weight vector w for restriction of Krylov space");
  }

  // construct c by setting all pressure values to 1.0 and export to c
  presmode->PutScalar(1.0);
  Teuchos::RCP<Epetra_Vector> tmpc = LINALG::CreateVector(*(discret_->DofRowMap()),true);
  LINALG::Export(*presmode,*tmpc);
  Teuchos::RCP<Epetra_Vector> tmpkspc = kspsplitter_->ExtractKSPCondVector(*tmpc);
  LINALG::Export(*tmpkspc,*c0);

  // fillcomplete the projector to compute (w^T c)^(-1)
  projector_->FillComplete();
  //c0->Print(std::cout);

  return;
} // CombustFluidImplicitTimeInt::UpdateKrylovSpaceProjection

/*--------------------------------------------------------------------------*
 | check if constant pressure mode is in kernel of sysmat_     nissen Jan13 |
 *--------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::CheckMatrixNullspace()
{
  //Note: this check is expensive and should only be used in the debug mode
  if (projector_ != Teuchos::null)
  {
    Teuchos::RCP<Epetra_MultiVector> c = projector_->GetNonConstKernel();
    projector_->FillComplete();
    int nsdim = c->NumVectors();
    if (nsdim != 1)
      dserror("Only one mode, namely the constant pressure mode, expected.");

    Epetra_Vector result(c->Map(),false);

    sysmat_->Apply(*c,result);

    double norm=1e9;

    result.Norm2(&norm);

    if(norm>1e-12)
    {
      std::cout << "#####################################################" << std::endl;
      std::cout << "Nullspace check for sysmat_ failed!                  " << std::endl;
      std::cout << "This might be caused by:                             " << std::endl;
      std::cout << " - you don't have pure Dirichlet boundary conditions " << std::endl;
      std::cout << "   or pbcs. pressure level is fixed. -> check datfile" << std::endl;
      std::cout << " - you don't integrate pressure dofs accurately      " << std::endl;
      std::cout << "   enough for sysmat_. constant pressure is not in   " << std::endl;
      std::cout << "   kernel of sysmat_. -> use more gauss points (often" << std::endl;
      std::cout << "   problem with nurbs)                               " << std::endl;
      std::cout << " - unlikely but not impossible: nullspace vector is  " << std::endl;
      std::cout << "   not the constant pressure mode (not totally clear " << std::endl;
      std::cout << "   for xfem, yet). In this case sysmat_ could be     " << std::endl;
      std::cout << "   correct. -> adapt nullspace vector                " << std::endl;
      std::cout << "#####################################################" << std::endl;
      dserror("Nullspace check for sysmat_ failed, Ac returned %12.5e",norm);
    }
  }

  return;
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
  COMBUST::UTILS::CalculateAcceleration(
      state_.velnp_, state_.veln_, state_.velnm_, state_.accn_,
          timealgo_, step_, theta_, dta_, dtp_,
          state_.accnp_);

  return;
}// FluidImplicitTimeInt::TimeUpdate


/*----------------------------------------------------------------------*
 | update acceleration for generalized-alpha time integration  vg 02/09 |
 *----------------------------------------------------------------------*/
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

  Teuchos::RCP<Epetra_Vector> onlyaccnp = Teuchos::rcp(new Epetra_Vector(onlyaccn->Map()));

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

    Teuchos::RCP<Epetra_Vector> onlyaccam = Teuchos::rcp(new Epetra_Vector(onlyaccnp->Map()));

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
  // ##########################################################################################
  // #                        WARNING                                                         #
  // # Be careful when adding new fields to the output section:                               #
  // # All maps are stored in the DiscretizationWriter. If they have been changed by calling  #
  // # FillComplete on the discretization, these new maps are added to the stack while        #
  // # outdated ones are not removed. This has a significant negative impact on the memory    #
  // # requirements, since the memory is not freed and the stored pointers point on nothing,  #
  // # i.e., memory leak. In combustion, FillComplete is called in every time step after the  #
  // # new flame front has been incorporated. However, since velocity and pressure are        #
  // # written on the standard dofset, which does not change, this problem is circumvented.   #
  // # Likewise, it is circumvented for the Owner, which is written only once now. In case of #
  // # redistribution, also the standard dofset has to be rebuild and the new Owner also to   #
  // # be written. Therefore, the map stack of the DiscretizationWriter is erased and all     #
  // # outdated maps are deleted (see below).                                                 #
  // ##########################################################################################
  output_->ClearMapCache();

  const bool write_visualization_data = step_%upres_ == 0;
  const bool write_restart_data = step_!=0 and uprestart_ != 0 and step_%uprestart_ == 0;
  const bool do_time_sample = special_flow_!="no" && step_>=samstart_ && step_<=samstop_;

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
          standarddofset_);
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
      }

      // output finescale velocity field for visualization
      if(fssgv_!= INPAR::FLUID::no_fssgv or turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales)
      {
        // transfer to standard dofset
        Teuchos::RCP<Epetra_Vector> fsvel_out = dofmanagerForOutput_->transformXFEMtoStandardVector(
          *fsvelafXfem_, *standarddofset_, state_.nodalDofDistributionMap_, outputfields);

        output_->WriteVector("fsvelocity", fsvel_out);
      }

      // write domain decomposition for visualization
      output_->WriteElementData(redist_this_step_);
      redist_this_step_ = false;

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
      //std::cout << "mindist " << std::setw(18) << std::setprecision(12) << std::scientific << mindist << IO::endl;
      if (lnode->Id()==163) //or
          //lnode->Id()==498 or
          //lnode->Id()==726 or
          //lnode->Id()==242)
      {
        IO::cout << "node " << lnode->Id() << " ";
        // extract velocity values (no pressure!) from global velocity vector
        for (int icomp=0; icomp<3; ++icomp)
        {
          lids[icomp] = velnp_out->Map().LID(dofids[icomp]);
          vel(icomp) = (*velnp_out)[lids[icomp]];
          IO::cout << std::setw(18)<< std::setprecision(12) <<vel(icomp) << " ";
        }
        IO::cout << "coordinates ";
        for (int icomp=0; icomp<3; ++icomp)
        {
          IO::cout << nodecoord(icomp) << " ";
        }
        IO::cout << IO::endl;
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
      //std::cout << "mindist " << std::setw(18) << std::setprecision(12) << std::scientific << mindist << IO::endl;
      if (lnode->Id()==149)

      {
        IO::cout << "node " << lnode->Id() << " ";
        // extract velocity values (no pressure!) from global velocity vector
        for (int icomp=0; icomp<3; ++icomp)
        {
          lids[icomp] = velnp_out->Map().LID(dofids[icomp]);
          vel(icomp) = (*velnp_out)[lids[icomp]];
          IO::cout << std::setw(18)<< std::setprecision(12) <<vel(icomp) << " ";
        }
        IO::cout << "coordinates ";
        for (int icomp=0; icomp<3; ++icomp)
        {
          IO::cout << nodecoord(icomp) << " ";
        }
        IO::cout << IO::endl;
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
      //std::cout << "mindist " << std::setw(18) << std::setprecision(12) << std::scientific << mindist << IO::endl;
      if (lnode->Id()==227)
      {
        IO::cout << "node " << lnode->Id() << " ";
        // extract velocity values (no pressure!) from global velocity vector
        for (int icomp=0; icomp<3; ++icomp)
        {
          lids[icomp] = velnp_out->Map().LID(dofids[icomp]);
          vel(icomp) = (*velnp_out)[lids[icomp]];
          IO::cout << std::setw(18)<< std::setprecision(12) <<vel(icomp) << " ";
        }
        IO::cout << "coordinates ";
        for (int icomp=0; icomp<3; ++icomp)
        {
          IO::cout << nodecoord(icomp) << " ";
        }
        IO::cout << IO::endl;
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
      //std::cout << "mindist " << std::setw(18) << std::setprecision(12) << std::scientific << mindist << IO::endl;
      if (lnode->Id()==59)
      {
        IO::cout << "node " << lnode->Id() << " ";
        // extract velocity values (no pressure!) from global velocity vector
        for (int icomp=0; icomp<3; ++icomp)
        {
          lids[icomp] = velnp_out->Map().LID(dofids[icomp]);
          vel(icomp) = (*velnp_out)[lids[icomp]];
          IO::cout << std::setw(18)<<  std::setprecision(12) <<vel(icomp) << " ";
        }
        IO::cout << "coordinates ";
        for (int icomp=0; icomp<3; ++icomp)
        {
          IO::cout << nodecoord(icomp) << " ";
        }
        IO::cout << IO::endl;
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
//      //std::cout << "mindist " << std::setw(18) << std::setprecision(12) << std::scientific << mindist << IO::endl;
//      if (lnode->Id()==517 or
//          lnode->Id()==725 or
//          lnode->Id()==451 or
//          lnode->Id()==243)
//      {
//        IO::cout << "node " << lnode->Id() << " ";
//        // extract velocity values (no pressure!) from global velocity vector
//        for (int icomp=0; icomp<3; ++icomp)
//        {
//          lids[icomp] = velnp_out->Map().LID(dofids[icomp]);
//          vel(icomp) = (*velnp_out)[lids[icomp]];
//          IO::cout << std::setprecision(12) <<vel(icomp) << " ";
//        }
//        IO::cout << IO::endl;
//      }
//    }
#endif
    }
  }

  // write restart
  if (write_restart_data)
  {
    if (myrank_ == 0)
      IO::cout << "---  write restart... ";// << IO::flush;
    //IO::cout << state_.velnp_->GlobalLength() << IO::endl;
    output_->WriteVector("velnp", state_.velnp_);
    //IO::cout << state_.veln_->GlobalLength() << IO::endl;
    output_->WriteVector("veln" , state_.veln_);
    //IO::cout << state_.velnm_->GlobalLength() << IO::endl;
    output_->WriteVector("velnm", state_.velnm_);
    //IO::cout << state_.accnp_->GlobalLength() << IO::endl;
    output_->WriteVector("accnp", state_.accnp_);
    //IO::cout << state_.accn_->GlobalLength() << IO::endl;
    output_->WriteVector("accn" , state_.accn_);
    const Teuchos::RCP<Epetra_Vector> phinprow = Teuchos::rcp(new Epetra_Vector(*discret_->NodeRowMap()));
    LINALG::Export(*phinp_,*phinprow);
    output_->WriteVector("phinp" , phinprow);
    if (timealgo_ == INPAR::FLUID::timeint_afgenalpha)
    {
      //IO::cout << state_.velaf_->GlobalLength() << IO::endl;
      output_->WriteVector("velaf" , state_.velaf_);
      //IO::cout << state_.accam_->GlobalLength() << IO::endl;
      output_->WriteVector("accam" , state_.accam_);
    }
    if (myrank_ == 0)
      IO::cout << "done" << IO::endl;
  }

  //if (step_ % 10 == 0 or step_== 1) //write every 5th time step only
  {
    OutputToGmsh((char*)"solution_field_pressure",(char*)"solution_field_velocity",step_, time_);
#ifdef GMSH_REF_FIELDS
    OutputToGmsh((char*)"mod_start_field_pres",(char*)"mod_start_field_vel",step_, time_);
#endif
  }

#ifdef FLAME_VORTEX
  OutputFlameArea(step_,time_);
#endif

#ifdef DL_INSTAB
  OutputInstabAmplitude(step_,time_);
#endif

  // --------------------------
  // dump turbulence statistics
  // --------------------------
  turbstatisticsmanager_->DoOutput((*output_), step_, 0);


  // -------------------------------------------------------------------
  //     calculate and write lift'n'drag forces from the residual
  // -------------------------------------------------------------------
  LiftDrag();

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

  //IO::cout << state_.velnp_->GlobalLength() << IO::endl;
  //IO::cout << state_.veln_->GlobalLength()  << IO::endl;
  //IO::cout << state_.velnm_->GlobalLength() << IO::endl;

  reader.ReadVector(state_.velnp_,"velnp");
  //IO::cout << state_.velnp_->GlobalLength() << IO::endl;
  reader.ReadVector(state_.veln_, "veln");
  //IO::cout << state_.veln_->GlobalLength() << IO::endl;
  reader.ReadVector(state_.velnm_,"velnm");
  //IO::cout << state_.velnm_->GlobalLength() << IO::endl;
  reader.ReadVector(state_.accnp_ ,"accnp");
  //IO::cout << state_.accnp_->GlobalLength() << IO::endl;
  reader.ReadVector(state_.accn_ ,"accn");
  //IO::cout << state_.accn_->GlobalLength() << IO::endl;
  if (timealgo_ == INPAR::FLUID::timeint_afgenalpha)
  {
    reader.ReadVector(state_.velaf_ ,"velaf");
    //IO::cout << state_.velaf_->GlobalLength() << IO::endl;
    reader.ReadVector(state_.accam_ ,"accam");
    //IO::cout << state_.accam_->GlobalLength() << IO::endl;
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

  // check maps before writing fields
#ifdef DEBUG
  // get map of this vector
  const Epetra_BlockMap& phimap = phinp_->Map();
  // check, whether this map is still identical with the current node map in the discretization
  if (not phimap.SameAs(*discret_->NodeColMap())) dserror("node column map has changed!");
#endif

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
        std::vector<int> lm;
        std::vector<int> lmowner;
        std::vector<int> lmstride;
        ele->LocationVector(*(discret_), lm, lmowner, lmstride);
        // extract local values from the global vector
        std::vector<double> myvelnp(lm.size());
        DRT::UTILS::ExtractMyValues(*output_col_vel, myvelnp, lm);

        const int numparam = eledofman.NumDofPerField(XFEM::PHYSICS::Pres);
        const std::vector<int>& dofpos = eledofman.LocalDofPosPerField(XFEM::PHYSICS::Pres);
        // get pressure values for this element
        LINALG::SerialDenseMatrix presele(1,numparam);
        for (int iparam=0; iparam<numparam; ++iparam)
          presele(0,iparam) = myvelnp[dofpos[iparam]];

        //---------------------------------------------------------------
        // extract local level-set (G-function) values from global vector
        //---------------------------------------------------------------
        size_t numnode = ele->NumNode();
        std::vector<double> myphinp(numnode);
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
          {
            dserror("unknown type of combustion problem!");
            break;
          }
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
            size_t numnode = ele->NumNode();
            std::vector<double> myphinp(numnode);
            // extract G-function values to element level
            DRT::UTILS::ExtractMyNodeBasedValues(ele, myphinp, *phinp_);
#ifdef DEBUG
            if (numnode != 8) dserror("pressure jump output only available for hex8 elements!");
#endif
            LINALG::Matrix<8,1> ephi;
            for (size_t iparam=0; iparam<numnode; ++iparam)
              ephi(iparam) = myphinp[iparam];

            // create local copy of information about dofs
            const XFEM::ElementDofManager eledofman(*ele,physprob_.elementAnsatz_->getElementAnsatz(ele->Shape()),*dofmanagerForOutput_);

            //-------------------------------------------
            // extract pressure values from global vector
            //-------------------------------------------
            std::vector<int> lm;
            std::vector<int> lmowner;
            std::vector<int> lmstride;
            ele->LocationVector(*(discret_), lm, lmowner, lmstride);
            // extract local values from the global vector
            std::vector<double> myvelnp(lm.size());
            DRT::UTILS::ExtractMyValues(*output_col_vel, myvelnp, lm);

            const size_t numparam = eledofman.NumDofPerField(XFEM::PHYSICS::Pres);
            const std::vector<int>& dofpos = eledofman.LocalDofPosPerField(XFEM::PHYSICS::Pres);
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
                if (numparam != 16) dserror("pressure jump output only available for fully enriched hex8 elements!");
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
                  //std::cout << "fac " << fac << IO::endl;
                  presjumpnorm += fac*presjump*presjump;
                }
#endif // computation average pressure
              }
            }
          }
        }
#ifdef COLLAPSE_FLAME
        IO::cout << IO::endl;
        // get the processor local node
        DRT::Node* lnode = discret_->lRowNode(0);

        LINALG::Matrix<3,1> nodecoord(true);
        // get physical coordinates of this node
        nodecoord(0) = lnode->X()[0];
        nodecoord(1) = lnode->X()[1];
        nodecoord(2) = lnode->X()[2];

        const double zcoord = nodecoord(2);
        const double deltaz = 2.*abs(zcoord);
        IO::cout << "Netz " << 1./deltaz << IO::endl;
        const double pi = atan(1.)*4.;
        const double area = pi*2.*0.25*deltaz;
        const double avpresjump = sqrt(presjumpnorm/area);
        IO::cout << "avpresjump " << avpresjump << IO::endl;
#endif
        gmshfilecontent << "};\n";
      }

      {
        gmshfilecontent << "View \" " << "Pressure Enrichment \" {\n";

        const Epetra_Map* dofcolmap = discret_->DofColMap();
        std::map<XFEM::DofKey,XFEM::DofGID> dofColDistrib;
        dofmanagerForOutput_->fillNodalDofColDistributionMap(dofColDistrib);

        for (int iele=0; iele<discret_->NumMyRowElements(); ++iele)
        {
          const DRT::Element* ele = discret_->lRowElement(iele);

          // output only for bisected elements (fully enriched elements) and touched elements
          if (interfacehandle_->ElementBisected(ele->Id()) or interfacehandle_->ElementTouched(ele->Id()))
          {
            // get node coordinates for this element
            const size_t numnode = DRT::UTILS::getNumberOfElementNodes(ele->Shape());
            LINALG::SerialDenseMatrix xyze(3,numnode);
            GEO::fillInitialPositionArray(ele,xyze);

            // vector for enrichment values
            LINALG::SerialDenseVector enrichmentval(numnode);

            for (size_t inode = 0; inode<numnode; inode++)
            {
              // node gid
              int nodegid = ele->Nodes()[inode]->Id();
              // setup dof key
              const XFEM::FieldEnr fieldenr(XFEM::PHYSICS::Pres,XFEM::Enrichment(XFEM::Enrichment::typeJump,0));
              const XFEM::DofKey dofkey(nodegid, fieldenr);
              // get gid of correseponding dof of node is enriched
              int dofgid = 0;
              if (dofColDistrib.find(dofkey) != dofColDistrib.end())
              {
                dofgid = dofColDistrib.find(dofkey)->second;
                // store value
                enrichmentval(inode) = (*output_col_vel)[dofcolmap->LID(dofgid)];
              }
              else
                enrichmentval(inode) = 0.0;
            }

            // write to file
            IO::GMSH::cellWithScalarFieldToStream(ele->Shape(), enrichmentval, xyze, gmshfilecontent);
          }// end split
        }// elements
        gmshfilecontent << "};\n";
      }

      }
#endif
    }
    gmshfilecontent.close();
    if (screen_out) IO::cout << " done" << IO::endl;
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

        std::vector<int> lm;
        std::vector<int> lmowner;
        std::vector<int> lmstride;
        ele->LocationVector(*(discret_), lm, lmowner, lmstride);

        // extract local values from the global vector
        std::vector<double> myvelnp(lm.size());
        DRT::UTILS::ExtractMyValues(*output_col_vel, myvelnp, lm);

        const int numparam = eledofman.NumDofPerField(field);
        const std::vector<int>& dofpos = eledofman.LocalDofPosPerField(field);

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
        std::vector<double> myphinp(numnode);
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
    if (screen_out) IO::cout << " done" << IO::endl;
  }
#endif

#if 0
  if (gmshoutput_)
  {
    std::ostringstream filename;
    std::ostringstream filenamedel;
    const std::string filebase = DRT::Problem::Instance()->OutputControlFile()->FileName();
    filename    << filebase << ".solution_field_pressure_disc_" << std::setw(5) << std::setfill('0') << step   << ".pos";
    filenamedel << filebase << ".solution_field_pressure_disc_" << std::setw(5) << std::setfill('0') << step-5 << ".pos";
    std::remove(filenamedel.str().c_str());
    if (screen_out) IO::cout << "writing " << std::left << std::setw(50) <<filename.str()<<"...";
    std::ofstream gmshfilecontent(filename.str().c_str());

    const XFEM::PHYSICS::Field field = XFEM::PHYSICS::DiscPres;

    {
      gmshfilecontent << "View \" " << "Discontinous Pressure Solution (Physical) \" {\n";
      for (int i=0; i<discret_->NumMyRowElements(); ++i)
      {
        const DRT::Element* ele = discret_->lRowElement(i);

        static LINALG::Matrix<3,27> xyze_xfemElement;
        GEO::fillInitialPositionArray(ele,xyze_xfemElement);

        const std::map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> element_ansatz(COMBUST::getElementAnsatz(ele->Shape()));

        // create local copy of information about dofs
        const XFEM::ElementDofManager eledofman(*ele,element_ansatz,*dofmanagerForOutput_);

        std::vector<int> lm;
        std::vector<int> lmowner;
        ele->LocationVector(*(discret_), lm, lmowner);

        // extract local values from the global vector
        std::vector<double> myvelnp(lm.size());
        DRT::UTILS::ExtractMyValues(*output_col_vel, myvelnp, lm);

        const int numparam = eledofman.NumDofPerField(field);
        const std::vector<int>& dofpos = eledofman.LocalDofPosPerField(field);

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
    if (screen_out) IO::cout << " done" << IO::endl;
  }
#endif

#if 0
  if (gmshoutput_)
  {
    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("solution_field_sigma_disc", step, 5, screen_out, discret_->Comm().MyPID());
    IO::cout << IO::endl;
    const std::string filenamexx = IO::GMSH::GetNewFileNameAndDeleteOldFiles("solution_field_sigmaxx_disc", step, 5, screen_out, discret_->Comm().MyPID());
    IO::cout << IO::endl;
    const std::string filenameyy = IO::GMSH::GetNewFileNameAndDeleteOldFiles("solution_field_sigmayy_disc", step, 5, screen_out, discret_->Comm().MyPID());
    IO::cout << IO::endl;
    const std::string filenamezz = IO::GMSH::GetNewFileNameAndDeleteOldFiles("solution_field_sigmazz_disc", step, 5, screen_out, discret_->Comm().MyPID());
    IO::cout << IO::endl;
    const std::string filenamexy = IO::GMSH::GetNewFileNameAndDeleteOldFiles("solution_field_sigmaxy_disc", step, 5, screen_out, discret_->Comm().MyPID());
    IO::cout << IO::endl;
    const std::string filenamexz = IO::GMSH::GetNewFileNameAndDeleteOldFiles("solution_field_sigmaxz_disc", step, 5, screen_out, discret_->Comm().MyPID());
    IO::cout << IO::endl;
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
      gmshfilecontent   << "View \" " << "Discontinous Stress Solution (Physical) \" {" << IO::endl;
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

        std::vector<int> lm;
        std::vector<int> lmowner;
        ele->LocationVector(*(discret_), lm, lmowner);

        // extract local values from the global vector
        std::vector<double> myvelnp(lm.size());
        DRT::UTILS::ExtractMyValues(*state_.velnp_, myvelnp, lm);

        const int numparam = eledofman.NumDofPerField(field);
        const std::vector<int>& dofposxx = eledofman.LocalDofPosPerField(XFEM::PHYSICS::Sigmaxx);
        const std::vector<int>& dofposyy = eledofman.LocalDofPosPerField(XFEM::PHYSICS::Sigmayy);
        const std::vector<int>& dofposzz = eledofman.LocalDofPosPerField(XFEM::PHYSICS::Sigmazz);
        const std::vector<int>& dofposxy = eledofman.LocalDofPosPerField(XFEM::PHYSICS::Sigmaxy);
        const std::vector<int>& dofposxz = eledofman.LocalDofPosPerField(XFEM::PHYSICS::Sigmaxz);
        const std::vector<int>& dofposyz = eledofman.LocalDofPosPerField(XFEM::PHYSICS::Sigmayz);

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
    if (screen_out) IO::cout << " done" << IO::endl;
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
 | write flame area to a file                                                         henke 06/12 |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::OutputFlameArea(
    const int    step,
    const double time
    ) const
{
  //-------------------------
  // write flame surface area
  //-------------------------
  if (combusttype_ == INPAR::COMBUST::combusttype_premixedcombustion)
  {

    double globflamearea = 0.0;
    double locflamearea = 0.0;

    for (int iele=0; iele<discret_->NumMyRowElements(); ++iele)
    {
      const DRT::Element* ele = discret_->lRowElement(iele);

      // output only for bisected elements
      if (interfacehandle_->ElementSplit(ele->Id()))
      {
        size_t numnode = ele->NumNode();
        if (numnode != 8) dserror("flame area output only available for hex8 elements!");

        // get node coordinates for this element
        LINALG::SerialDenseMatrix xyze(3,numnode);
        GEO::fillInitialPositionArray(ele,xyze);

        //--------------------------------------------------------------------
        // interpolate values from element to Gaussian points of boundary cell
        //--------------------------------------------------------------------
        const GEO::BoundaryIntCells& boundaryintcells = interfacehandle_->ElementBoundaryIntCells(ele->Id());
        for (GEO::BoundaryIntCells::const_iterator cell = boundaryintcells.begin(); cell != boundaryintcells.end(); ++cell)
        {
          // here, a triangular boundary integration cell is assumed (numvertices = 3)
          if (cell->Shape() == DRT::Element::tri3)
          {
            const size_t numvertices = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement;
            //else if
            //const size_t numvertices = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement;

            // choose an (arbitrary) number of Gaussian points for the output
            DRT::UTILS::GaussRule2D intrule2D = DRT::UTILS::intrule_tri_3point;
            DRT::UTILS::IntegrationPoints2D intpoints(intrule2D);

            // loop over Gaussian points
            for (int iquad=0; iquad<intpoints.nquad; ++iquad)
            {
              // compute area of boundary integration cell
              {
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
                //std::cout << "fac " << fac << IO::endl;
                locflamearea += fac; // *1.0
              }
            }
          }
          else if (cell->Shape() == DRT::Element::quad4)
          {
            const size_t numvertices = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement;
            //else if
            //const size_t numvertices = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement;

            // choose an (arbitrary) number of Gaussian points for the output
            DRT::UTILS::GaussRule2D intrule2D = DRT::UTILS::intrule_quad_4point;
            DRT::UTILS::IntegrationPoints2D intpoints(intrule2D);

            // loop over Gaussian points
            for (int iquad=0; iquad<intpoints.nquad; ++iquad)
            {
              // compute area of boundary integration cell
              {
                LINALG::SerialDenseMatrix cellXiDomaintmp = cell->CellNodalPosXiDomain();
                const LINALG::Matrix<3,numvertices> cellXiDomain(cellXiDomaintmp);

                const LINALG::Matrix<2,1> gpinEta2D(intpoints.qxg[iquad]);

                // jacobian for coupled transformation
                // get derivatives dxi_3D/deta_2D
                static LINALG::Matrix<2,numvertices> deriv_eta2D;
                DRT::UTILS::shape_function_2D_deriv1(deriv_eta2D,gpinEta2D(0,0),gpinEta2D(1,0),DRT::Element::quad4);

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
                //std::cout << "fac " << fac << IO::endl;
                locflamearea += fac; // *1.0
              }
            }
          }
          else
            dserror("Not implemented for this cell type");
        }
      }
    }

    // sum over all procs
    discret_->Comm().SumAll(&locflamearea,&globflamearea,1);

    if (myrank_ == 0)
    {
      //-------------------------
      // write flame area to file
      //-------------------------
      std::ostringstream tmpfilename;
      const std::string filebase(DRT::Problem::Instance()->OutputControlFile()->FileName());

      tmpfilename << filebase << "." << "flame_area" << ".txt";
      IO::cout << "writing " << std::left << std::setw(60) <<tmpfilename.str()<<"...";

      const std::string filename = tmpfilename.str();

      // output to log-file
      Teuchos::RCP<std::ofstream> log;
      log = Teuchos::rcp(new std::ofstream(filename.c_str(),std::ios::app));
      (*log) <<  " "  << std::setw(6) << step_ << "  " << std::setw(6) << std::setprecision(6) << time_ << "  " << std::setprecision(13) << globflamearea;
      (*log) << &std::endl;
      log->flush();

      IO::cout << " done" << IO::endl;
    }
  }
}



/*------------------------------------------------------------------------------------------------*
 | write amplitude of Darrieus Landau instability to a file                           henke 06/12 |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::OutputInstabAmplitude(
    const int    step,
    const double time
    ) const
{
  //----------------
  // write amplitude
  //----------------
  if (combusttype_ == INPAR::COMBUST::combusttype_premixedcombustion)
  {
    if (myrank_ != 0)
      dserror("output of Darrieus Landau instability only serial!");

    double amplitudemiddle = 0.0;
    double amplitudeleft = 0.0;
    double amplitude;

    const double pi = 3.141592653590;
    LINALG::Matrix<3,1> xyz(true);

    double smallest1Gvalmiddle = 2.0*pi/5.0;
    double smallest2Gvalmiddle  = 2.0*pi/5.0;
    double smallest1Yvalmiddle  = 0.0;
    double smallest2Yvalmiddle  = 0.0;

    double smallest1Gvalleft = 2.0*pi/5.0;
    double smallest2Gvalleft = 2.0*pi/5.0;
    double smallest1Yvalleft = 0.0;
    double smallest2Yvalleft = 0.0;

    double Yvalbottommiddle  = 0.0;
    double Gvalbottommiddle  = 2.0*pi/5.0;
    double Yvaltopmiddle     = 0.0;
    double Gvaltopmiddle     = 2.0*pi/5.0;

    double Yvalbottomleft = 0.0;
    double Gvalbottomleft = 2.0*pi/5.0;
    double Yvaltopleft = 0.0;
    double Gvaltopleft = 2.0*pi/5.0;

    //--------------------------------
    // loop all nodes on the processor
    //--------------------------------
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node* lnode = discret_->lRowNode(lnodeid);

      // get node coordinates
      for(int idim=0;idim<3;idim++)
        xyz(idim)=lnode->X()[idim];

      if ((fabs(xyz(0)-(2.0*pi/10.0)) < 1e-6) and (xyz(2)>0)) // node on middle line
      {
        // get phi value for this node
        const int lid = phinp_->Map().LID(lnode->Id());
        const double gfuncval = (*phinp_)[lid];

        if (fabs(gfuncval) < fabs(smallest1Gvalmiddle))
        {
          smallest2Gvalmiddle  = smallest1Gvalmiddle ;
          smallest2Yvalmiddle  = smallest1Yvalmiddle ;
          smallest1Gvalmiddle  = gfuncval;
          smallest1Yvalmiddle  = xyz(1);
        }
        else if (fabs(gfuncval) < fabs(smallest2Gvalmiddle ))
        {
          smallest2Gvalmiddle  = gfuncval;
          smallest2Yvalmiddle  = xyz(1);
        }
        else
        {
          // do nothing, G-value does not belong to the smallest values
        }
      }
      if ((fabs(xyz(0)-0.0) < 1e-6) and (xyz(2)>0)) // node on left periodic boundary
      {
        // get phi value for this node
        const int lid = phinp_->Map().LID(lnode->Id());
        const double gfuncval = (*phinp_)[lid];

        if (fabs(gfuncval) < fabs(smallest1Gvalleft))
        {
          smallest2Gvalleft = smallest1Gvalleft;
          smallest2Yvalleft = smallest1Yvalleft;
          smallest1Gvalleft = gfuncval;
          smallest1Yvalleft = xyz(1);
        }
        else if (fabs(gfuncval) < fabs(smallest2Gvalleft))
        {
          smallest2Gvalleft = gfuncval;
          smallest2Yvalleft = xyz(1);
        }
        else
        {
          // do nothing, G-value does not belong to the smallest values
        }
      }
    }

    if (smallest2Yvalmiddle  < smallest1Yvalmiddle )
    {
      Yvalbottommiddle = smallest2Yvalmiddle;
      Gvalbottommiddle = smallest2Gvalmiddle;
      Yvaltopmiddle    = smallest1Yvalmiddle;
      Gvaltopmiddle    = smallest1Gvalmiddle;
    }
    else
    {
      Yvalbottommiddle = smallest1Yvalmiddle;
      Gvalbottommiddle = smallest1Gvalmiddle;
      Yvaltopmiddle    = smallest2Yvalmiddle;
      Gvaltopmiddle    = smallest2Gvalmiddle;
    }

    if (smallest2Yvalleft  < smallest1Yvalleft )
    {
      Yvalbottomleft = smallest2Yvalleft;
      Gvalbottomleft = smallest2Gvalleft;
      Yvaltopleft    = smallest1Yvalleft;
      Gvaltopleft    = smallest1Gvalleft;
    }
    else
    {
      Yvalbottomleft = smallest1Yvalleft;
      Gvalbottomleft = smallest1Gvalleft;
      Yvaltopleft    = smallest2Yvalleft;
      Gvaltopleft    = smallest2Gvalleft;
    }

    amplitudemiddle = (Yvaltopmiddle  + Yvalbottommiddle *fabs(Gvaltopmiddle )/fabs(Gvalbottommiddle ))/(1 + fabs(Gvaltopmiddle )/fabs(Gvalbottommiddle ));
    amplitudeleft = (Yvaltopleft + Yvalbottomleft *fabs(Gvaltopleft )/fabs(Gvalbottomleft ))/(1 + fabs(Gvaltopleft )/fabs(Gvalbottomleft ));
    amplitude = amplitudemiddle-amplitudeleft;

    if (myrank_ == 0)
    {
      //-------------------------
      // write flame area to file
      //-------------------------
      std::ostringstream tmpfilename;
      const std::string filebase(DRT::Problem::Instance()->OutputControlFile()->FileName());

      tmpfilename << filebase << "." << "amplitude" << ".txt";
      IO::cout << "writing " << std::left << std::setw(60) <<tmpfilename.str()<<"...";

      const std::string filename = tmpfilename.str();

      // output to log-file
      Teuchos::RCP<std::ofstream> log;
      log = Teuchos::rcp(new std::ofstream(filename.c_str(),std::ios::app));
      (*log) <<  " "  << std::setw(6) << step_ << "  " << std::setw(6) << std::setprecision(6) << time_;
      (*log) << "  " << std::setw(15) << std::setprecision(13) << amplitudeleft;
      (*log) << "  " << std::setw(15) << std::setprecision(13) << amplitudemiddle;
      (*log) << "  " << std::setw(15) << std::setprecision(13) << amplitude;
      (*log) << &std::endl;
      log->flush();

      IO::cout << " done" << IO::endl;
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

#ifdef DEBUG
  // get map of this vector
  const Epetra_BlockMap& phimap = phinp_->Map();
  // check, whether this map is still identical with the current node map in the discretization
  if (not phimap.SameAs(*discret_->NodeColMap())) dserror("node column map has changed!");
#endif

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
        // just define a default element ansatz; it is not used anyway
        const COMBUST::TauPressureAnsatz elementAnsatz;
        const XFEM::ElementDofManager eledofman(*ele,elementAnsatz.getElementAnsatz(ele->Shape()),*dofmanagerForOutput_);

        std::vector<int> lm;
        std::vector<int> lmowner;
        std::vector<int> lmstride;
        ele->LocationVector(*discret_, lm, lmowner, lmstride);

        // extract local values from the global vector
        std::vector<double> myvelnp(lm.size());
        DRT::UTILS::ExtractMyValues(*vectorfield, myvelnp, lm);

        const std::vector<int>& dofposvelx = eledofman.LocalDofPosPerField(XFEM::PHYSICS::Velx);
        const std::vector<int>& dofposvely = eledofman.LocalDofPosPerField(XFEM::PHYSICS::Vely);
        const std::vector<int>& dofposvelz = eledofman.LocalDofPosPerField(XFEM::PHYSICS::Velz);

        const int numparamvelx = eledofman.NumDofPerField(XFEM::PHYSICS::Velx);
        LINALG::SerialDenseMatrix elementvalues(4, numparamvelx,true);
        for (int iparam=0; iparam<numparamvelx; ++iparam)
        {
          elementvalues(0, iparam) = myvelnp[dofposvelx[iparam]];
          elementvalues(1, iparam) = myvelnp[dofposvely[iparam]];
          elementvalues(2, iparam) = myvelnp[dofposvelz[iparam]];
        }

        //------------------------------------------------------------------------------------------
        // extract local level-set (G-function) values from global vector
        //------------------------------------------------------------------------------------------
        size_t numnode = ele->NumNode();
        std::vector<double> myphinp(numnode);
        // extract G-function values to element level
        DRT::UTILS::ExtractMyNodeBasedValues(ele, myphinp, *phinp_);

        const GEO::DomainIntCells& domainintcells = interfacehandle_->ElementDomainIntCells(ele->Id());
        for (GEO::DomainIntCells::const_iterator cell = domainintcells.begin(); cell != domainintcells.end(); ++cell)
        {
          LINALG::SerialDenseMatrix cellvalues(3, DRT::UTILS::getNumberOfElementNodes(cell->Shape()));

          switch(combusttype_)
          {
          case INPAR::COMBUST::combusttype_premixedcombustion:
          case INPAR::COMBUST::combusttype_twophaseflowjump:
          {
            //
            XFEM::InterpolateCellValuesFromElementValuesLevelSet(*ele, eledofman, *cell, myphinp, XFEM::PHYSICS::Velx,
              elementvalues, cellvalues);
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
          {
            dserror("unknown type of combustion problem!");
            break;
          }
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
            size_t numnode = ele->NumNode();
            std::vector<double> myphinp(numnode);
            // extract G-function values to element level
            DRT::UTILS::ExtractMyNodeBasedValues(ele, myphinp, *phinp_);

#ifdef DEBUG
            if (numnode != 8) dserror("velocity jump output only available for hex8 elements!");
#endif
            LINALG::Matrix<8,1> ephi;
            for (size_t iparam=0; iparam<numnode; ++iparam)
              ephi(iparam) = myphinp[iparam];

            // create local copy of information about dofs
            const XFEM::ElementDofManager eledofman(*ele,physprob_.elementAnsatz_->getElementAnsatz(ele->Shape()),*dofmanagerForOutput_);

            //-------------------------------------------
            // extract velocity values from global vector
            //-------------------------------------------
            std::vector<int> lm;
            std::vector<int> lmowner;
            std::vector<int> lmstride;
            ele->LocationVector(*(discret_), lm, lmowner, lmstride);
            // extract local values from the global vector
            std::vector<double> myvelnp(lm.size());
            DRT::UTILS::ExtractMyValues(*vectorfield, myvelnp, lm);

            const std::vector<int>& dofposvelx = eledofman.LocalDofPosPerField(XFEM::PHYSICS::Velx);
            const std::vector<int>& dofposvely = eledofman.LocalDofPosPerField(XFEM::PHYSICS::Vely);
            const std::vector<int>& dofposvelz = eledofman.LocalDofPosPerField(XFEM::PHYSICS::Velz);

            const size_t numparamvelx = eledofman.NumDofPerField(XFEM::PHYSICS::Velx);
            LINALG::SerialDenseMatrix evel(4, numparamvelx,true);
            for (size_t iparam=0; iparam<numparamvelx; ++iparam)
            {
              evel(0,iparam) = myvelnp[dofposvelx[iparam]];
              evel(1,iparam) = myvelnp[dofposvely[iparam]];
              evel(2,iparam) = myvelnp[dofposvelz[iparam]];
            }
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
                if (numparamvelx != 16) dserror("velocity jump output only available for fully enriched hex8 elements!");
#endif

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
        IO::cout << IO::endl;
        // get the processor local node
        DRT::Node* lnode = discret_->lRowNode(0);

        LINALG::Matrix<3,1> nodecoord(true);
        // get physical coordinates of this node
        nodecoord(0) = lnode->X()[0];
        nodecoord(1) = lnode->X()[1];
        nodecoord(2) = lnode->X()[2];

        const double zcoord = nodecoord(2);
        const double deltaz = 2.*abs(zcoord);
        IO::cout << "Netz " << 1./deltaz << IO::endl;
        const double pi = atan(1.)*4.;
        const double area = pi*2.*0.25*deltaz;
        const double avveljump = sqrt(veljumpnormsquare/area);
        IO::cout << "avveljump " << avveljump << IO::endl;
#endif
#endif
        gmshfilecontent << "};\n";
      }

      {
        gmshfilecontent << "View \" " << "Velocity Enrichment \" {\n";

        // get a copy on column parallel distribution
        Teuchos::RCP<const Epetra_Vector> output_col_vel = DRT::UTILS::GetColVersionOfRowVector(discret_, state_.velnp_);

        const Epetra_Map* dofcolmap = discret_->DofColMap();
        std::map<XFEM::DofKey,XFEM::DofGID> dofColDistrib;
        dofmanagerForOutput_->fillNodalDofColDistributionMap(dofColDistrib);

        for (int iele=0; iele<discret_->NumMyRowElements(); ++iele)
        {
          const DRT::Element* ele = discret_->lRowElement(iele);

          // output only for bisected elements(fully enriched elements) and touched elements
          if (interfacehandle_->ElementBisected(ele->Id()) or interfacehandle_->ElementTouched(ele->Id()))
          {
            // get node coordinates for this element
            const size_t numnode = DRT::UTILS::getNumberOfElementNodes(ele->Shape());
            LINALG::SerialDenseMatrix xyze(3,numnode);
            GEO::fillInitialPositionArray(ele,xyze);

            // vector for enrichment values
            LINALG::SerialDenseMatrix enrichmentval(3,numnode);

            for (size_t inode = 0; inode<numnode; inode++)
            {
              // node gid
              int nodegid = ele->Nodes()[inode]->Id();

              // setup dof key
              const XFEM::FieldEnr fieldenrvelx(XFEM::PHYSICS::Velx,XFEM::Enrichment(XFEM::Enrichment::typeJump,0));
              const XFEM::DofKey dofkeyvelx(nodegid, fieldenrvelx);
              // get gid of correseponding dof
              int dofgidx = 0;
              if (dofColDistrib.find(dofkeyvelx) != dofColDistrib.end())
              {
                dofgidx = dofColDistrib.find(dofkeyvelx)->second;
                // store value
                enrichmentval(0,inode) = (*output_col_vel)[dofcolmap->LID(dofgidx)];
              }
              else
                enrichmentval(0,inode) = 0.0;

              // setup dof key
              const XFEM::FieldEnr fieldenrvely(XFEM::PHYSICS::Vely,XFEM::Enrichment(XFEM::Enrichment::typeJump,0));
              const XFEM::DofKey dofkeyvely(nodegid, fieldenrvely);
              // get gid of correseponding dof
              int dofgidy = 0;
              if (dofColDistrib.find(dofkeyvely) != dofColDistrib.end())
              {
                dofgidy = dofColDistrib.find(dofkeyvely)->second;
                // store value
                enrichmentval(1,inode) = (*output_col_vel)[dofcolmap->LID(dofgidy)];
              }
              else
                enrichmentval(1,inode) = 0.0;

              // setup dof key
              const XFEM::FieldEnr fieldenrvelz(XFEM::PHYSICS::Velz,XFEM::Enrichment(XFEM::Enrichment::typeJump,0));
              const XFEM::DofKey dofkeyvelz(nodegid, fieldenrvelz);
              // get gid of correseponding dof
              int dofgidz = 0;
              if (dofColDistrib.find(dofkeyvelz) != dofColDistrib.end())
              {
                dofgidz = dofColDistrib.find(dofkeyvelz)->second;
                // store value
                enrichmentval(2,inode) = (*output_col_vel)[dofcolmap->LID(dofgidz)];
              }
              else
                enrichmentval(2,inode) = 0.0;
            }

              // write to file
              IO::GMSH::cellWithVectorFieldToStream(ele->Shape(), enrichmentval, xyze, gmshfilecontent);
            }// end split
          }// elements
        gmshfilecontent << "};\n";
      }

      }
#endif
    }
    gmshfilecontent.close();
    if (screen_out) IO::cout << " done" << IO::endl;
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
  case INPAR::COMBUST::initfield_zero_field:
  {
    // nothing to do
    break;
  }
  case INPAR::COMBUST::initfield_beltrami_flow:
  {
    const Epetra_Map* dofrowmap = discret_->DofRowMap();

    int err =0;

    const int npredof = numdim_;

    double         p;
    std::vector<double> u  (numdim_);
    std::vector<double> acc(numdim_);
    std::vector<double> xyz(numdim_);

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
      std::vector<int> nodedofset = discret_->Dof(lnode);

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
      int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_fluid);
      if (id==-1) dserror("Newtonian fluid material could not be found");
      const MAT::PAR::Parameter* mat = DRT::Problem::Instance()->Materials()->ParameterById(id);
      const MAT::PAR::NewtonianFluid* actmat = static_cast<const MAT::PAR::NewtonianFluid*>(mat);
      double dens = actmat->density_;
      double visc = actmat->viscosity_;

      p = -a*a/2.0 * dens *
        ( exp(2.0*a*xyz[0])
          + exp(2.0*a*xyz[1])
          + exp(2.0*a*xyz[2])
          + 2.0 * sin(a*xyz[0] + d*xyz[1]) * cos(a*xyz[2] + d*xyz[0]) * exp(a*(xyz[1]+xyz[2]))
          + 2.0 * sin(a*xyz[1] + d*xyz[2]) * cos(a*xyz[0] + d*xyz[1]) * exp(a*(xyz[2]+xyz[0]))
          + 2.0 * sin(a*xyz[2] + d*xyz[0]) * cos(a*xyz[1] + d*xyz[2]) * exp(a*(xyz[0]+xyz[1]))
          );

      // Beltrami is always 3D
      acc[0] = u[0]*(-1.0*d*d*visc/dens);
      acc[1] = u[1]*(-1.0*d*d*visc/dens);
      acc[2] = u[2]*(-1.0*d*d*visc/dens);

      // set initial velocity components
      for(int nveldof=0;nveldof<numdim_;nveldof++)
      {
        const int gid = nodedofset[nveldof];
        int lid = dofrowmap->LID(gid);
        err += state_.velnp_->ReplaceMyValues(1,&(u[nveldof]),&lid);
        err += state_.veln_ ->ReplaceMyValues(1,&(u[nveldof]),&lid);
        err += state_.velnm_->ReplaceMyValues(1,&(u[nveldof]),&lid);

        // set additionally the values for the time derivative to start with an exact acceleration in case of OST (theta!=1.0)
        // set initial acceleration components
        err += state_.accnp_->ReplaceMyValues(1,&(acc[nveldof]),&lid);
        err += state_.accn_ ->ReplaceMyValues(1,&(acc[nveldof]),&lid);
        if (timealgo_ == INPAR::FLUID::timeint_afgenalpha)
         err += state_.accam_->ReplaceMyValues(1,&(acc[nveldof]),&lid);
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
      const std::vector<int> nodedofset = discret_->Dof(lnode);

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
        IO::cout << "Disturbed initial profile:   max. " << perc*100 << "% random perturbation\n";
        IO::cout << "\n\n";
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
        std::vector<int> nodedofset = discret_->Dof(lnode);

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
        std::vector<int> nodedofset = discret_->Dof(lnode);

        // check whether we have a pbc condition on this node
        std::vector<DRT::Condition*> mypbc;

        lnode->GetCondition("SurfacePeriodic",mypbc);

        // check whether a periodic boundary condition is active on this node
        if (mypbc.size()>0)
        {
          // yes, we have one

          // get the list of all his slavenodes
          std::map<int, std::vector<int> >::iterator master = (discret_->GetAllPBCCoupledColNodes())->find(lnode->Id());

          // slavenodes are ignored
          if(master == (discret_->GetAllPBCCoupledColNodes())->end()) continue;
        }

        // add random noise on initial function field
        for(int index=0;index<numdim;++index)
        {
          int gid = nodedofset[index];

          double randomnumber = DRT::Problem::Instance(0)->Random()->Uni();

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

//std::cout << *dofrowmap << IO::endl;
//std::cout << *(standarddofset_->DofRowMap()) << IO::endl;
//std::cout << (state_.velnp_->Map()) << IO::endl;

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
      const std::vector<int> nodedofs = (*standarddofset_).Dof(lnode);
      //const std::vector<int> nodedofs = discret_->Dof(lnode);
      //for (int i=0;i<standardnodedofset.size();i++)
      //{
      //  IO::cout << "component " << i << " standarddofset dofid " << stdnodedofset[i] << IO::endl;
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
  //----------------------------------------------------------------------------------------------
  // flame-vortex interaction problem: two counter-rotating vortices (2-D) moving the  flame front
  //----------------------------------------------------------------------------------------------
  case INPAR::COMBUST::initfield_darrieus_landau_instability:
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

    // get laminar burning velocity (flame speed)
    if (flamespeed_ != 1.0) dserror("flame speed should be 1.0 for the 'Darrieus-Landau-instability' case");

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
    if (dens_b != 0.2) dserror("burnt density should be 0.2 for the 'Darrieus-Landau-instability' case");
    const double dens_u = mat1->Density();
    if (dens_u != 1.0) dserror("unburnt density should be 1.0 for the 'Darrieus-Landau-instability' case");
    double dens = dens_u;
    // for "pure fluid" computation: rhob = rhou = 1.161
    //const double dens_b = dens_u;

    // get map of global velocity vectors (DofRowMap)
    const Epetra_Map* dofrowmap = standarddofset_->DofRowMap();
    //const Epetra_Map* dofrowmap = discret_->DofRowMap();

//std::cout << *dofrowmap << IO::endl;
//std::cout << *(standarddofset_->DofRowMap()) << IO::endl;
//std::cout << (state_.velnp_->Map()) << IO::endl;

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

      //----------------------------------------
      // set density with respect to flame front
      //----------------------------------------
      if (gfuncval >= 0.0) // plus/burnt domain -> burnt material
      {
        dens = dens_b;
        pres = 0.0;
      }
      else // minus/unburnt domain -> unburnt material
      {
        dens = dens_u;
        pres = flamespeed_*flamespeed_*dens_u*dens_u*(1.0/dens_u - 1.0/dens_b);
      }
      //----------------------------------------------
      // compute components of initial velocity vector
      //----------------------------------------------
      vel(0) = 0.0;
      vel(1) = flamespeed_*dens_u/dens;
      // 2D problem -> vel_z = 0.0
      vel(2) = 0.0;
      // velocity profile without vortices
      //vel(1) = sl*densu/dens;

      // access standard FEM dofset (3 x vel + 1 x pressure) to get dof IDs for this node
      const std::vector<int> nodedofs = (*standarddofset_).Dof(lnode);
      //const std::vector<int> nodedofs = discret_->Dof(lnode);
      //for (int i=0;i<standardnodedofset.size();i++)
      //{
      //  IO::cout << "component " << i << " standarddofset dofid " << stdnodedofset[i] << IO::endl;
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
  {
    dserror("type of initial field not available");
    break;
  }
  }

  // -------------------------------------------------------------------
  //           preparation of AVM3-based scale separation
  // -------------------------------------------------------------------

  // remark: we decided to call this function here, since
  //         (i) it has to be called only once
  //         (ii) we do not want to have any XFEM dofs etc already
  // (ii) means that merely the standard field is used for scale-separtion,
  // and the enrichment is thus part for the fine-scales
   if (fssgv_ != INPAR::FLUID::no_fssgv or turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales)
     AVM3Preparation();

  return;
}


void FLD::CombustFluidImplicitTimeInt::SetEnrichmentField(
    const Teuchos::RCP<XFEM::DofManager> dofmanager,
    const Epetra_Map dofrowmap)
{
#ifdef FLAME_VORTEX
  IO::cout << "---  set initial enrichment field for flame-vortex interaction example... " << IO::flush;

  // initial field modification for flame_vortex_interaction
  for (int nodeid=0;nodeid<discret_->NumMyRowNodes();nodeid++) // loop over element nodes
  {
    DRT::Node* lnode = discret_->lRowNode(nodeid);
    const std::set<XFEM::FieldEnr>& fieldenrset(dofmanager->getNodeDofSet(lnode->Id()));
    for (std::set<XFEM::FieldEnr>::const_iterator fieldenr = fieldenrset.begin();
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
  IO::cout << "done" << IO::endl;
#endif

#ifdef DL_INSTAB
  IO::cout << "---  set initial enrichment field for Darrieus-Landau instability example... " << IO::flush;

  // initial field modification for flame_vortex_interaction
  for (int nodeid=0;nodeid<discret_->NumMyRowNodes();nodeid++) // loop over element nodes
  {
    DRT::Node* lnode = discret_->lRowNode(nodeid);
    const std::set<XFEM::FieldEnr>& fieldenrset(dofmanager->getNodeDofSet(lnode->Id()));
    for (std::set<XFEM::FieldEnr>::const_iterator fieldenr = fieldenrset.begin();
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
          (*state_.veln_)[dofrowmap.LID(dofpos)] = 4.5;
          (*state_.velnp_)[dofrowmap.LID(dofpos)] = 4.5;
        }
        else if (fieldenr->getField() == XFEM::PHYSICS::Velz);
          // nothing to do in flame_vortex example
        else if (fieldenr->getField() == XFEM::PHYSICS::Pres)
        {
          (*state_.veln_)[dofrowmap.LID(dofpos)] = 4.5;
          (*state_.velnp_)[dofrowmap.LID(dofpos)] = 4.5;
        }
      } // end if jump enrichment
    } // end loop over fieldenr
  } // end loop over element nodes
  IO::cout << "done" << IO::endl;
#endif

#ifdef COLLAPSE_FLAME
  IO::cout << "---  set initial enrichment field for collapsing flame example... " << IO::flush;

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
    for (std::set<XFEM::FieldEnr>::const_iterator fieldenr = fieldenrset.begin();
        fieldenr != fieldenrset.end();++fieldenr)
    {
      const XFEM::DofKey dofkey(lnode->Id(), *fieldenr);
      const int dofpos = state_.nodalDofDistributionMap_.find(dofkey)->second;
      if (fieldenr->getEnrichment().Type() == XFEM::Enrichment::typeJump)
      {
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
  IO::cout << "done" << IO::endl;
#endif

#ifdef COMBUST_TWO_FLAME_FRONTS
  IO::cout << "---  set initial enrichment field for two approaching flame fronts example... " << IO::flush;

  // initial field modification for 2-flames example
  for (int nodeid=0;nodeid<discret_->NumMyRowNodes();nodeid++) // loop over element nodes
  {
    DRT::Node* lnode = discret_->lRowNode(nodeid);

    LINALG::Matrix<3,1> coords(lnode->X());

    const int lid = phinp_->Map().LID(lnode->Id());
    const double gfuncval = (*phinp_)[lid];

    const std::set<XFEM::FieldEnr>& fieldenrset(dofmanager->getNodeDofSet(lnode->Id()));
    for (std::set<XFEM::FieldEnr>::const_iterator fieldenr = fieldenrset.begin();
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
  IO::cout << "done" << IO::endl;
#endif

//#ifdef GMSH_REF_FIELDS
  OutputToGmsh((char*)"mod_start_field_pres",(char*)"mod_start_field_vel",Step(), Time());
//#endif
}


// -------------------------------------------------------------------
// -------------------------------------------------------------------
void FLD::CombustFluidImplicitTimeInt::SetupXFluidSplit(
        const DRT::Discretization& dis,
        const Teuchos::RCP<XFEM::DofManager> dofman,
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
  std::vector<int> velmapdata;
  std::vector<int> premapdata;

  // collect global dofids for velocity and pressure in vectors
  for (int i=0; i<dis.NumMyRowNodes(); ++i) {
    const DRT::Node* node = dis.lRowNode(i);
    const std::set<XFEM::FieldEnr>& enrvarset(dofman->getNodeDofSet(node->Id()));
    const std::vector<int> dof = dis.Dof(node);
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
  Teuchos::RCP<Epetra_Map> velrowmap = Teuchos::rcp(new Epetra_Map(-1,
      velmapdata.size(),&velmapdata[0],0,
      dis.Comm()));
  Teuchos::RCP<Epetra_Map> prerowmap = Teuchos::rcp(new Epetra_Map(-1,
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
  Teuchos::ParameterList eleparams;
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

  // create dummy instance of interfacehandle holding no flamefront and hence no integration cells
  Teuchos::RCP<COMBUST::InterfaceHandleCombust> ihdummy = Teuchos::rcp(new COMBUST::InterfaceHandleCombust(discret_,Teuchos::null));
  // create dummy instance of dof manager assigning standard enrichments to all nodes
  const Teuchos::RCP<XFEM::DofManager> dofmanagerdummy = Teuchos::rcp(new XFEM::DofManager(ihdummy,Teuchos::null,physprob_.xfemfieldset_,xparams_,Teuchos::null));

  // save dofmanager to be able to plot Gmsh stuff in Output()
  dofmanagerForOutput_ = dofmanagerdummy;

  // pass dof information to elements (no enrichments yet, standard FEM!)
  TransferDofInformationToElements(Teuchos::null, dofmanagerdummy);

  // ensure that degrees of freedom in the discretization have been set
  if ((not discret_->Filled()) or (not discret_->HaveDofs()))
    discret_->FillComplete();

  // update the dofset with the complete fluid unknowns
  standarddofset_->SetCoupledNodes(discret_->GetAllPBCCoupledColNodes());
  standarddofset_->Reset();
  standarddofset_->AssignDegreesOfFreedom(*discret_,0,0);

  // split based on complete fluid field
  FLD::UTILS::SetupFluidSplit(*discret_,*standarddofset_,3,*velpressplitterForOutput_);

  // rebuid internal faces
  if(params_->sublist("RESIDUAL-BASED STABILIZATION").get<std::string>("STABTYPE")=="edge_based" or
     DRT::INPUT::IntegralValue<bool>(params_->sublist("COMBUSTION FLUID"),"XFEMSTABILIZATION") == true)
  {
    // do some checks first
    if(params_->sublist("RESIDUAL-BASED STABILIZATION").get<std::string>("STABTYPE")=="edge_based" and
       DRT::INPUT::IntegralValue<bool>(params_->sublist("COMBUSTION FLUID"),"XFEMSTABILIZATION") == true)
       dserror("Combination of face-based stabilization and ghost-penalty stabilization for XFEM currently not supported!");

    // if the definition of internal faces would be included
    // in the standard discretization, these lines can be removed
    // and CreateInternalFacesExtension() can be called once
    // in the constructor of the fluid time integration
    // since we want to keep the standard discretization as clean as
    // possible, we create interal faces via an enhanced discretization
    // including the faces between elements
    xfemdiscret_->CreateInternalFacesExtension(true);
  }

  // remember that we did a redist
  redist_this_step_ = true;

  return;
}


/*--------------------------------------------------------------------------------------------*
 | Redistribute the fluid discretization and vectors according to nodegraph   rasthofer 07/11 |
 |                                                                            DA wichmann     |
 *--------------------------------------------------------------------------------------------*/
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
      xparams_,
      discret_->GetAllPBCCoupledColNodes())
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

  // create parameters for discretization
  Teuchos::ParameterList eleparams;

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
  INPAR::COMBUST::InitialField calcerr = DRT::INPUT::get<INPAR::COMBUST::InitialField>(*params_, "eval err for analyt sol");

  //------------------------------------------------------- beltrami flow
  switch (calcerr)
  {
  case INPAR::COMBUST::initfield_zero_field:
  case INPAR::COMBUST::initfield_field_by_function:
  case INPAR::COMBUST::initfield_disturbed_field_by_function:
  case INPAR::COMBUST::initfield_flame_vortex_interaction:
  case INPAR::COMBUST::initfield_darrieus_landau_instability:
    // do nothing --- no analytical solution available
    break;
  case INPAR::COMBUST::initfield_beltrami_flow:
  {
    // create the parameters for the discretization
    Teuchos::ParameterList eleparams;

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
    discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
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
  {
    dserror("Cannot calculate error. Unknown type of analytical test problem");
    break;
  }
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
//   dserror("LiftDrag() not yet implemented for combustion problems");
//  // in this map, the results of the lift drag calculation are stored
//  Teuchos::RCP<std::map<int,std::vector<double> > > liftdragvals;
//
//  FLD::UTILS::LiftDrag(*discret_,*trueresidual_,*params_,liftdragvals);
//
//  if (liftdragvals!=Teuchos::null and discret_->Comm().MyPID() == 0)
//    FLD::UTILS::WriteLiftDragToFile(time_, step_, *liftdragvals);

  return;
} // CombustFluidImplicitTimeInt::LiftDrag

/*----------------------------------------------------------------------*
 | prepare AVM3-based scale separation                  rasthofer 03/14 |
 *----------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::AVM3Preparation()
{
  // time measurement: avm3
  TEUCHOS_FUNC_TIME_MONITOR("           + avm3");

  // pass dof information to elements (no enrichments yet, standard FEM!)
  TransferDofInformationToElements(Teuchos::null, dofmanagerForOutput_);

  // zero matrix
  sysmat_->Zero();

  Teuchos::ParameterList eleparams;
  dta_      = params_->get<double> ("time step size");

  if(timealgo_==INPAR::FLUID::timeint_afgenalpha)
  {
    // recall original parameters from input file
    alphaM_ = params_->get<double>("alpha_M");
    alphaF_ = params_->get<double>("alpha_F");
    gamma_  = params_->get<double>("gamma");
    // compute "pseudo-theta" for af-generalized-alpha scheme
    theta_ = alphaF_*gamma_/alphaM_;
  }
  else
  {
    theta_    = params_->get<double>("theta");
  }

  // other parameters needed by the elements
  eleparams.set("time",time_);
  eleparams.set("dt",dta_);
  eleparams.set("thsl",theta_*dta_);
  eleparams.set("theta",theta_);

  // evaluate Neumann conditions
  neumann_loads_= LINALG::CreateVector(*discret_->DofRowMap(),true);
  discret_->EvaluateNeumann(eleparams,*neumann_loads_);
  discret_->ClearState();
  residual_->Update(1.0,*neumann_loads_,0.0);

  // action for elements
  eleparams.set("action","calc_fluid_systemmat_and_residual");

  // flag for type of combustion problem
  eleparams.set<int>("combusttype",combusttype_);
  eleparams.set<int>("veljumptype",veljumptype_);
  eleparams.set<int>("fluxjumptype",fluxjumptype_);
  eleparams.set<double>("flamespeed",flamespeed_);
  eleparams.set<double>("marksteinlength",marksteinlength_);
  eleparams.set<double>("nitschevel",nitschevel_);
  eleparams.set<double>("nitschepres",nitschepres_);

  // parameter for suppressing additional enrichment dofs in two-phase flow problems
  eleparams.set<int>("selectedenrichment",DRT::INPUT::IntegralValue<INPAR::COMBUST::SelectedEnrichment>(params_->sublist("COMBUSTION FLUID"),"SELECTED_ENRICHMENT"));

  // parameters for two-phase flow problems with surface tension
  eleparams.set<int>("surftensapprox",surftensapprox_);
  eleparams.set<double>("variablesurftens",params_->sublist("COMBUSTION FLUID").get<double>("VARIABLESURFTENS"));
  eleparams.set<bool>("connected_interface",connected_interface_);
  eleparams.set<bool>("second_deriv",DRT::INPUT::IntegralValue<int>(params_->sublist("COMBUSTION FLUID"),"L2_PROJECTION_SECOND_DERIVATIVES"));

  // smoothed normal vectors for boundary integration
  eleparams.set("smoothed_bound_integration",smoothed_boundary_integration_);
  eleparams.set<int>("smoothgradphi",smoothgradphi_);

  // other parameters that might be needed by the elements
  eleparams.set<int>("timealgo",timealgo_);
  eleparams.set<double>("gamma",gamma_);
  eleparams.set<double>("alphaF",alphaF_);
  eleparams.set<double>("alphaM",alphaM_);

  // additional terms for Nitsche's method (see Diss Florian)
  eleparams.set<bool>("nitsche_convflux",DRT::INPUT::IntegralValue<int>(params_->sublist("COMBUSTION FLUID"),"NITSCHE_CONVFLUX"));
  eleparams.set<bool>("nitsche_convstab",DRT::INPUT::IntegralValue<int>(params_->sublist("COMBUSTION FLUID"),"NITSCHE_CONVSTAB"));
  eleparams.set<bool>("nitsche_convpenalty",DRT::INPUT::IntegralValue<int>(params_->sublist("COMBUSTION FLUID"),"NITSCHE_CONVPENALTY"));
  // further terms and parameters related to Nitsche's method
  eleparams.set<bool>("nitsche_mass",DRT::INPUT::IntegralValue<int>(params_->sublist("COMBUSTION FLUID"),"NITSCHE_MASS"));
  eleparams.set<INPAR::COMBUST::WeightType>("weighttype",DRT::INPUT::IntegralValue<INPAR::COMBUST::WeightType>(params_->sublist("COMBUSTION FLUID"),"NITSCHE_WEIGHT"));

  //type of linearisation: include reactive terms for linearisation
  if(DRT::INPUT::get<INPAR::FLUID::LinearisationAction>(*params_, "Linearisation") == INPAR::FLUID::Newton)
    eleparams.set<bool>("include reactive terms for linearisation",true);
  else
    eleparams.set<bool>("include reactive terms for linearisation",false);

  // parameters for stabilization
  eleparams.sublist("RESIDUAL-BASED STABILIZATION") = params_->sublist("RESIDUAL-BASED STABILIZATION");

  // parameters for turbulence
  eleparams.set<INPAR::FLUID::TurbModelAction>("turbmodel",turbmodel_);
  eleparams.set<double>("Cs",Cs_);
  eleparams.set<INPAR::FLUID::FineSubgridVisc>("fssgv",fssgv_);

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

  // further level-set values (gradient etc) are not required, since boundary integrals are not evaluated
  // all elements are considered uncut
  eleparams.set<Teuchos::RCP<Epetra_Vector> >("phinpcol",phinp_);
  eleparams.set<Teuchos::RCP<COMBUST::InterfaceHandleCombust> >("interfacehandle_",interfacehandle_);

  // set fine-scale vector
  // dummy vector initialized with zeros
  // Remark:
  // This is necessary because the fssgv_ flag
  // has already been set in SetParameters()
  // Therefore, the function Evaluate() already
  // expects the state vector "fsvelaf"
  if (fssgv_ != INPAR::FLUID::no_fssgv or turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales)
  {
    discret_->SetState("fsvelaf",fsvelafXfem_);

    if (turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales)
      eleparams.sublist("MULTIFRACTAL SUBGRID SCALES") = params_->sublist("MULTIFRACTAL SUBGRID SCALES");
  }

  // element evaluation for getting system matrix
  // -> we merely need matrix "structure" below, not the actual contents
  // call standard loop over elements
  discret_->Evaluate(eleparams,sysmat_,residual_);

  discret_->ClearState();

  // complete system matrix
  sysmat_->Complete();

  // apply DBC to system matrix
  // predicted dirichlet values
  // velnp then also holds prescribed new dirichlet values
  dbcmaps_ = Teuchos::rcp(new LINALG::MapExtractor());
  xfemdiscret_->EvaluateDirichletCombust(eleparams,state_.velnp_,Teuchos::null,Teuchos::null,Teuchos::null,dbcmaps_);
  LINALG::ApplyDirichlettoSystem(sysmat_,incvel_,residual_,zeros_,*(dbcmaps_->CondMap()));

  // get scale-separation matrix
  {
    // this is important to have!!!
    // MLAPI::Init() without arguments uses internally MPI_COMM_WOLRD
    MLAPI::Init();

    // extract the ML parameters:
    Teuchos::ParameterList&  mlparams = solver_->Params().sublist("ML Parameters");
    // remark: we create a new solver with ML preconditioner here, since this allows for also using other solver setups
    // to solve the system of equations
    // get the solver number used form the multifractal subgrid-scale model parameter list
    const int scale_sep_solvernumber = params_->sublist("MULTIFRACTAL SUBGRID SCALES").get<int>("ML_SOLVER");
    if (scale_sep_solvernumber != (-1))    // create a dummy solver
    {
      Teuchos::RCP<LINALG::Solver> solver = Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->SolverParams(scale_sep_solvernumber),
                                            discret_->Comm(),
                                            DRT::Problem::Instance()->ErrorFile()->Handle()));
      // compute the null space,
      discret_->ComputeNullSpaceIfNecessary(solver->Params(),true);
      // and, finally, extract the ML parameters
      mlparams = solver->Params().sublist("ML Parameters");
    }

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
    Teuchos::RCP<Epetra_CrsMatrix> crsPtent;
    MLAPI::GetPtent(*SystemMatrix()->EpetraMatrix(),mlparams,nullspace,crsPtent);
    LINALG::SparseMatrix Ptent(crsPtent,LINALG::View);

    // compute scale-separation matrix: S = I - Ptent*Ptent^T
    Sep_ = LINALG::Multiply(Ptent,false,Ptent,true);
    Sep_->Scale(-1.0);
    Teuchos::RCP<Epetra_Vector> tmp = LINALG::CreateVector(Sep_->RowMap(),false);
    tmp->PutScalar(1.0);
    Teuchos::RCP<Epetra_Vector> diag = LINALG::CreateVector(Sep_->RowMap(),false);
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
}// CombustFluidImplicitTimeInt::AVM3Preparation


/*----------------------------------------------------------------------*
 | AVM3-based scale separation                         rasthofer 03/14  |
 *----------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::AVM3Separation()
{
  // time measurement: avm3
  TEUCHOS_FUNC_TIME_MONITOR("           + avm3");

  // create dofrowmap of standard dofs on standard dofset
  //std::cout<<"AVM3 Separation!"<<std::endl;

  std::set<XFEM::PHYSICS::Field> physicalfields;
  physicalfields.insert(XFEM::PHYSICS::Velx);
  physicalfields.insert(XFEM::PHYSICS::Vely);
  physicalfields.insert(XFEM::PHYSICS::Velz);
  physicalfields.insert(XFEM::PHYSICS::Pres);

  // coarse scale vectors
  // based in standard and XFEM dofset
  Teuchos::RCP<Epetra_Vector> csvelafStd = Teuchos::null;
  Teuchos::RCP<Epetra_Vector> csvelafXFEM = Teuchos::null;

  // -------------------------------------------------------
  // calculate fine scale velocity based on standard dofset
  // -------------------------------------------------------

  if (timealgo_ == INPAR::FLUID::timeint_afgenalpha)
  {
    dserror("Genalpha not supported, yet!");
  }
  else
  {
    // vectors and matrix using initial unenriched dof field

    // from current enriched field to standard field
    Teuchos::RCP<Epetra_Vector> velnpstddof = dofmanagerForOutput_->transformXFEMtoStandardVector(
       *state_.velnp_, *standarddofset_, state_.nodalDofDistributionMap_, physicalfields);

    // from standard field to initial unenriched dof field
    Teuchos::RCP<Epetra_Vector> velnpstd = dofmanagerForOutput_->transformStandardToXFEMVector(
       *velnpstddof,*standarddofset_, plainnodalDofDistributionMap_, physicalfields, &(fsvelafStd_->Map()));

    Sep_->Multiply(false,*velnpstd,*fsvelafStd_);

    // from initial unenriched dof field to standard field
    Teuchos::RCP<Epetra_Vector> fsvelnpstd = dofmanagerForOutput_->transformXFEMtoStandardVector(
       *fsvelafStd_, *standarddofset_, plainnodalDofDistributionMap_, physicalfields);

    // from standard field to current enriched field
    Teuchos::RCP<Epetra_Vector> fsvelnpxfem = dofmanagerForOutput_->transformStandardToXFEMVector(
       *fsvelnpstd,*standarddofset_, state_.nodalDofDistributionMap_, physicalfields,discret_->DofRowMap());

    if (not excludeXfem_)
    {
      Teuchos::RCP<Epetra_Vector> csvelnpstd = Teuchos::rcp(new Epetra_Vector(fsvelafStd_->Map(),true));
      csvelnpstd->Update(1.0,*velnpstd,-1.0,*fsvelafStd_,0.0);

      // from initial unenriched dof field to standard field
      Teuchos::RCP<Epetra_Vector> csvelnpstddof = dofmanagerForOutput_->transformXFEMtoStandardVector(
         *csvelnpstd, *standarddofset_, plainnodalDofDistributionMap_, physicalfields);

      // from standard field to current enriched field
      Teuchos::RCP<Epetra_Vector> csvelnpxfem = dofmanagerForOutput_->transformStandardToXFEMVector(
         *csvelnpstddof,*standarddofset_, state_.nodalDofDistributionMap_, physicalfields,discret_->DofRowMap());

      fsvelafXfem_->Update(1.0,*state_.velnp_,-1.0,*csvelnpxfem,0.0);
    }
    else
      fsvelafXfem_->Update(1.0,*fsvelnpxfem,0.0);
  }

  // set fine-scale vector
  discret_->SetState("fsvelaf",fsvelafXfem_);

  return;
}// CombustFluidImplicitTimeInt::AVM3Separation


/*----------------------------------------------------------------------*
 |                                                     rasthofer 03/14  |
 *----------------------------------------------------------------------*/
const Teuchos::RCP<const Epetra_Vector> FLD::CombustFluidImplicitTimeInt::Dirichlet()
{
  if (dbcmaps_ == Teuchos::null)
    dserror("Dirichlet map has not been allocated");

  Teuchos::RCP<Epetra_Vector> dirichones = LINALG::CreateVector(*(dbcmaps_->CondMap()),false);
  dirichones->PutScalar(1.0);
  Teuchos::RCP<Epetra_Vector> dirichtoggle = LINALG::CreateVector(*(discret_->DofRowMap()),true);
  dbcmaps_->InsertCondVector(dirichones, dirichtoggle);

  return dirichtoggle;
}


/*------------------------------------------------------------------------------------------------*
| returns matching std::string for each time integration scheme                         gjb 08/08 |
*-------------------------------------------------------------------------------------------------*/
std::string FLD::CombustFluidImplicitTimeInt::MapTimIntEnumToString(const enum INPAR::FLUID::TimeIntegrationScheme term)
{
  // length of return std::string is 14 due to usage in formated screen output
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
    break;
  }
  return "";
}


/*------------------------------------------------------------------------------------------------*
| evaluate node-based forces                                                      rasthofer 03/14 |
*-------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::ComputeExternalForces()
{
  // compute repellant bubble-wall interaction force for turbulent bubbly channel flow
  if (special_flow_ == "bubbly_channel_flow" and repellant_force_ == true)
  {
    external_loads_ = LINALG::CreateVector(*discret_->DofRowMap(),true);

    // reference for further details: Bolotnov et al. / International Journal of Multiphase Flow 31 (2011) 647-659
    //
    // F = mu * U * R * ( a_1 / d_w + a_2 / (d_w)^2) * n
    //
    // U: local mean axial fluid velocity -> computed here using a simplified version of inner layer velocity profile

    // band around interface for application of force
    const double epsilon = 0.0175;
    // bubble radius
    const double radius = 0.5;
    // further constants
    const double a1 = 550.0;
    const double a2 = 35.0;
    // friction Reynolds number
    const double Retau = 180.0;
    // friction velocity
    const double utau = 0.065833;
    // parameters for low-law
    const double kappa_inv = 1.0/0.41;
    const double B_log = 5.2;
    // dynamic viscosity of liquid phase
    const double visc = 0.00036574;

    // loop all nodes on this proc
    for (int inode=0; inode<discret_->NumMyRowNodes(); inode++)
    {
      // get node
      DRT::Node* actnode = discret_->lRowNode(inode);
      // get y-coordinate of node
      const double ycoord = actnode->X()[1];

      // check if node is close to the wall, i.e., within band for lubrication force
      const double distance = 1.0-std::abs(ycoord);
      if (distance < (4.0*epsilon))
      {
        // get gid of node
        const int nodegid = actnode->Id();
        // level-set value
        const double phinode = (*phinp_)[(discret_->NodeColMap())->LID(nodegid)];

        // check if node is in band around interface
        if (std::abs(phinode) < epsilon)
        {
          // compute local mean axial fluid velocity
          double U_ax = 0.0;

          // compute approximate axial velocity according to law for viscous sublayer
          const double u_visc = distance * Retau * utau;
          // compute approximate axial velocity according to law for log-layer
          double u_log = 1.0e12;
          if (distance > 1.0e-7)
            u_log = (kappa_inv * log(distance * Retau) + B_log) * utau;

          if (u_visc < u_log) U_ax = u_visc;
          else U_ax = u_log;

          // compute force
          const double force = visc * U_ax * radius * (a1/distance + a2/(distance*distance)) * (-1.0*ycoord/std::abs(ycoord));

          // insert in vector
          const std::set<XFEM::FieldEnr>& fieldenrset(dofmanagerForOutput_->getNodeDofSet(nodegid));
          for (std::set<XFEM::FieldEnr>::const_iterator fieldenr = fieldenrset.begin(); fieldenr != fieldenrset.end();++fieldenr)
          {
            if (fieldenr->getField() == XFEM::PHYSICS::Vely)
            {
              if (fieldenr->getEnrichment().Type() == XFEM::Enrichment::typeStandard)
              {
                const XFEM::DofKey dofkey(nodegid, *fieldenr);
                const int dofpos = state_.nodalDofDistributionMap_.find(dofkey)->second;
                (*external_loads_)[discret_->DofRowMap()->LID(dofpos)] = force;
              }
            }
          }
        } // end if in band
      } // end if distance wall

    } // end loop all nodes
  }
  return;
}

Teuchos::RCP<const Epetra_Vector> FLD::CombustFluidImplicitTimeInt::ReadPhinp(int step)     {
    IO::DiscretizationReader reader(discret_,step);
    const Teuchos::RCP<Epetra_Vector> phinprow = Teuchos::rcp(new Epetra_Vector(*discret_->NodeRowMap()));
    reader.ReadVector(phinprow,"phinp");
    const Teuchos::RCP<Epetra_Vector> phinpcol = Teuchos::rcp(new Epetra_Vector(*discret_->NodeColMap()));
    LINALG::Export(*phinprow,*phinpcol);
    return phinpcol; }
