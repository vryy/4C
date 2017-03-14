/*----------------------------------------------------------------------*/
/*!
\file xfluid.cpp

\brief Control routine for fluid (in)stationary solvers with XFEM,
       including instationary solvers for fluid and fsi problems coupled with an internal embedded interface

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 */
/*----------------------------------------------------------------------*/

#include "xfluid.H"
#include "xfluid_state_creator.H"
#include "xfluid_state.H"
#include "xfluid_outputservice.H"
#include "xfluidresulttest.H"

#include "../drt_lib/drt_discret_xfem.H"
#include "../drt_lib/drt_dofset_transparent_independent.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_assemblestrategy.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_colors.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_dofset_predefineddofnumber.H"

#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_krylov_projector.H"

#include "../drt_cut/cut_elementhandle.H"
#include "../drt_cut/cut_sidehandle.H"
#include "../drt_cut/cut_volumecell.H"
#include "../drt_cut/cut_cutwizard.H"

#include "../drt_io/io.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_io/io_control.H"

#include "../drt_fluid_ele/fluid_ele.H"
#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../drt_fluid_ele/fluid_ele_interface.H"
#include "../drt_fluid_ele/fluid_ele_factory.H"

#include "../drt_fluid/fluid_utils_infnormscaling.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"
#include "../drt_fluid/fluid_utils.H"

#include "../drt_inpar/inpar_parameterlist_utils.H"

#include "../drt_xfem/xfem_condition_manager.H"
#include "../drt_xfem/xfem_dofset.H"
#include "../drt_xfem/xfem_edgestab.H"
#include "../drt_xfem/xfem_neumann.H"

#include "../drt_xfem/xfluid_timeInt_std_SemiLagrange.H"
#include "../drt_xfem/xfluid_timeInt_base.H"
#include "../drt_xfem/xfluid_timeInt.H"
#include "../drt_xfem/xfem_discretization_utils.H"

#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/matlist.H"
#include "../drt_mat/matpar_bundle.H"

#include <Epetra_Time.h>


/*----------------------------------------------------------------------*
 |  Constructor for basic XFluid class                     schott 03/12 |
 *----------------------------------------------------------------------*/
FLD::XFluid::XFluid(
    const Teuchos::RCP<DRT::Discretization>&      actdis,
    const Teuchos::RCP<DRT::Discretization>&      mesh_coupdis,
    const Teuchos::RCP<DRT::Discretization>&      levelset_coupdis,
    const Teuchos::RCP<LINALG::Solver>&           solver,
    const Teuchos::RCP<Teuchos::ParameterList>&   params,
    const Teuchos::RCP<IO::DiscretizationWriter>& output,
    bool                                          alefluid /*= false*/)
  : FluidImplicitTimeInt(actdis, solver, params, output, alefluid),
    xdiscret_(Teuchos::rcp_dynamic_cast<DRT::DiscretizationXFEM>(actdis, true)),
    edgestab_(Teuchos::rcp(new XFEM::XFEM_EdgeStab())),
    turbmodel_(INPAR::FLUID::dynamic_smagorinsky)
{
  //TODO the initialization of coupling objects, dofsets, and so on is not that clear so far, however, strongly
  // depends on the calling algorithms and adapters. Maybe we can improve this at some point.

  // all discretizations which potentially include mesh-based XFEM coupling/boundary conditions
  meshcoupl_dis_.clear();
  levelsetcoupl_dis_.clear();

  if(mesh_coupdis != Teuchos::null)
    meshcoupl_dis_.push_back(mesh_coupdis);

  //TODO: remove this after fixing the SemiLagrangean time integration for multiple mesh coupling objects!
  mc_idx_ = 0; // using this constructor only one mesh coupling discretization is supported so far

  // add the background dis itself for boundary-fitted couplings
  meshcoupl_dis_.push_back(actdis);

  if(levelset_coupdis != Teuchos::null)
    levelsetcoupl_dis_.push_back(levelset_coupdis);


  if(levelsetcoupl_dis_.size()>1)
    dserror("so far the framework is tested just for one level-set coupling object");

  return;
}

void FLD::XFluid::SetDofSetCouplingMap(
    const std::map<std::string, int>& dofset_coupling_map
)
{
  dofset_coupling_map_ = dofset_coupling_map;
}

void FLD::XFluid::AddAdditionalScalarDofsetAndCoupling()
{
  // ensure that dofset with idx=1 in bg_dis carries a dofset with one dof per node to carry the levelset field
  // and to allow to use the bgdis also as a cutterdis (note: cutterdis vectors are based on a dofrowmap and not on a noderowmap...)

  Teuchos::RCP<DRT::DofSetInterface> dofsetaux = Teuchos::rcp(new DRT::DofSetPredefinedDoFNumber(1, 0 ,0, true ));

  // add the dofset to the xfluid dis
  const int dofidx = xdiscret_->AddDofSet(dofsetaux);

  // store the dof index in the dofset_coupling_map_ for right access through the coupling objects
  dofset_coupling_map_.insert(std::pair<std::string, int>("phi_scatra_proxy_in_fluid", dofidx));

  if (dofidx != 1) // the index for the phinp-dofset in the fluid dis we currently expect!!!
    dserror("unexpected dof sets in fluid field - check if the framework works properly also if dofidx != 1?");

  // assign degrees of freedom (as a new dofset has been added!)
  xdiscret_->FillComplete(true, false,false);

  //TODO: check if we can add this dofset and the actdis all the time, even if there is a scatra dis (maybe we would obtain two two-phase conditions?)
  levelsetcoupl_dis_.push_back(xdiscret_);
}

void FLD::XFluid::CheckInitializedDofSetCouplingMap()
{
  if(meshcoupl_dis_.size()>0)
  {
    //TODO: use the dofset_coupling_map_ also for mesh coupling objects!
//    if(dofset_coupling_map_.empty())
//      dserror("you first have to call SetDofSetCouplingMap() if there is a mesh coupling discretization");
  }

  if(levelsetcoupl_dis_.size()>0)
  {
    if(dofset_coupling_map_.empty())
      dserror("you first have to call SetDofSetCouplingMap() if there is a level-set coupling discretization");
    else
    {
      // do not add the additional scalar dofset
    }
  }
  else
  {
    // no scatra discretization is available and therefore also not scatra dofset proxy in the fluid dis
    // this is needed for potential level-set based coupling objects defined on the background discretization

    AddAdditionalScalarDofsetAndCoupling();
  }
}


/*----------------------------------------------------------------------*
 |  initialize algorithm                                   schott 11/14 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::Init(bool createinitialstate)
{

  CheckInitializedDofSetCouplingMap();


  FluidImplicitTimeInt::Init();

  // -------------------------------------------------------------------
  // get input params and print Xfluid specific configurations
  // -------------------------------------------------------------------

  // read xfluid input parameters from list
  SetXFluidParams();

  // check xfluid input parameter combination for consistency & valid choices
  CheckXFluidParams();

  // set element time parameter as ghost penalty solve are called already in the Init for SetInitialFlowField
  SetElementTimeParameter();

  // create internal faces, if not already done in base class init
  if (facediscret_ == Teuchos::null)
  {
    CreateFacesExtension();
  }
  // -------------------------------------------------------------------
  // create a Condition/Coupling Manager
  // -------------------------------------------------------------------

  condition_manager_ = Teuchos::rcp(new XFEM::ConditionManager(
      dofset_coupling_map_,
      discret_,
      meshcoupl_dis_,
      levelsetcoupl_dis_,
      time_,
      step_)
  );

  condition_manager_->Init();

  // build the whole object which then can be used
  condition_manager_->Setup();


  // -------------------------------------------------------------------
  // read restart for all cutter discretizations
  // -------------------------------------------------------------------

  // read the interface displacement and interface velocity for the old timestep which was written in Output
  // we have to do this before ReadRestart() is called to get the right
  // initial CUT corresponding to time t^n at which the last solution was written
  //
  // REMARK: ivelnp_ and idispnp_ will be set again for the new time step in PrepareXFEMSolve()

  const int restart = DRT::Problem::Instance()->Restart();

  if(restart) condition_manager_->ReadRestart(restart);


  //TODO: this has to be removed when different includeinner flags for level-set and mesh cuts can be handled in the cut library
  // -------------------------------------------------------------------
  // set include inner flag
  // -------------------------------------------------------------------

  Teuchos::RCP<XFEM::LevelSetCoupling> two_phase_coupl = condition_manager_->GetLevelSetCoupling("XFEMLevelsetTwophase");
  Teuchos::RCP<XFEM::LevelSetCoupling> combust_coupl   = condition_manager_->GetLevelSetCoupling("XFEMLevelsetCombustion");

  if(two_phase_coupl != Teuchos::null or combust_coupl!= Teuchos::null)
  {
    include_inner_ = true;

    if(condition_manager_->HasMeshCoupling())
    {
      // loop all mesh coupling objects
      for(int mc_idx=0; mc_idx< condition_manager_->NumMeshCoupling(); mc_idx++)
      {
        Teuchos::RCP<XFEM::MeshCoupling> mc_coupl = condition_manager_->GetMeshCoupling(mc_idx);

        if (mc_coupl->CutGeometry()) //Mesh cut and Two-Phase cut not allowed at the same time.
          dserror("two-phase flow coupling and mesh coupling at once is not supported by the cut at the moment, as Node-position and include inner are not handled properly then");
      }
    }
  }
  else
  {
    include_inner_=false;
  }


  // -------------------------------------------------------------------
  // create the state creator
  // -------------------------------------------------------------------
  state_creator_ = Teuchos::rcp(
      new FLD::XFluidStateCreator(
        condition_manager_,
        params_->sublist("XFEM"),
        maxnumdofsets_,
        minnumdofsets_,
        include_inner_));


  // -------------------------------------------------------------------
  // create output dofsets and prepare output for xfluid
  // -------------------------------------------------------------------

  // load GMSH output flags
  if (DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->IOParams(),"OUTPUT_GMSH"))
  {
    output_service_ = Teuchos::rcp(new XFluidOutputServiceGmsh(params_->sublist("XFEM"),xdiscret_,condition_manager_,include_inner_));
  }
  else
  {
    output_service_ = Teuchos::rcp(new XFluidOutputService(xdiscret_,condition_manager_));
  }

  // -------------------------------------------------------------------
  // GMSH discretization output before CUT
  // -------------------------------------------------------------------
  output_service_->GmshOutputDiscretization(eval_eos_, step_);

  // -------------------------------------------------------------------
  // Create velpresssplitter for uncut discretization.
  velpressplitter_std_ = Teuchos::rcp(new LINALG::MapExtractor());
  FLD::UTILS::SetupFluidSplit(*discret_,xdiscret_->InitialDofSet(),numdim_,*velpressplitter_std_);

  // -------------------------------------------------------------------
  // initialize ALE-specific fluid vectors based on the initial dof row map
  // -------------------------------------------------------------------

  if (alefluid_)
  {
    dispnp_  = LINALG::CreateVector(*xdiscret_->InitialDofRowMap(),true);
    dispn_   = LINALG::CreateVector(*xdiscret_->InitialDofRowMap(),true);
    dispnm_  = LINALG::CreateVector(*xdiscret_->InitialDofRowMap(),true);
    gridvnp_ = LINALG::CreateVector(*xdiscret_->InitialDofRowMap(),true);
    gridvn_  = LINALG::CreateVector(*xdiscret_->InitialDofRowMap(),true);
  }


  // -------------------------------------------------------------------
  // create the initial state class
  // -------------------------------------------------------------------
  // note that all vectors w.r.t np have to be set properly

  if (createinitialstate)
    CreateInitialState();

  return;
}// Init()


void FLD::XFluid::SetupFluidDiscretization()
{
  XFEM::UTILS::XFEMDiscretizationBuilder xdisbuilder;

  Teuchos::RCP<DRT::Discretization> xfluiddis;

  //TODO: we should try to resolve this confusing meaning of fluid dis and xfluid dis for xfluid and xfluidfluid!!!

  // XFF-case
  if (DRT::Problem::Instance()->DoesExistDis("xfluid"))
  {
    Teuchos::RCP<DRT::Discretization> fluiddis  = DRT::Problem::Instance()->GetDis("fluid"); // fluid dis is here the embedded mesh (required for XFFSI)
    xfluiddis = DRT::Problem::Instance()->GetDis("xfluid"); // xfluid dis is here the cut mesh
    xdisbuilder.SetupXFEMDiscretization(
      DRT::Problem::Instance()->XFEMGeneralParams(),
      xfluiddis,
      fluiddis,
      "FluidMesh"
    );
  }
  else // standard xfluid case
  {
    xfluiddis = DRT::Problem::Instance()->GetDis("fluid"); // fluid dis is here the cut mesh
    xdisbuilder.SetupXFEMDiscretization(
      DRT::Problem::Instance()->XFEMGeneralParams(),
      xfluiddis
    );
  }
}

/*----------------------------------------------------------------------*
 |  set all xfluid parameters                              schott 02/15 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::SetXFluidParams()
{

  omtheta_      = 1.0 - theta_;

  numdim_       = DRT::Problem::Instance()->NDim();

  Teuchos::ParameterList&   params_xfem    = params_->sublist("XFEM");
  Teuchos::ParameterList&   params_xf_gen  = params_->sublist("XFLUID DYNAMIC/GENERAL");
  Teuchos::ParameterList&   params_xf_stab = params_->sublist("XFLUID DYNAMIC/STABILIZATION");

  // get the maximal number of dofsets that are possible to use
  maxnumdofsets_ = params_->sublist("XFEM").get<int>("MAX_NUM_DOFSETS");

  xfluid_timintapproach_ = DRT::INPUT::IntegralValue<INPAR::XFEM::XFluidTimeIntScheme>(params_xf_gen,"XFLUID_TIMEINT");
  xfluid_timint_check_interfacetips_ = (bool)DRT::INPUT::IntegralValue<int>(params_xf_gen,"XFLUID_TIMEINT_CHECK_INTERFACETIPS");

  // for monolithic problems with xfluid (varying dofrowmaps)
  permutation_map_ = Teuchos::rcp(new std::map<int,int>);
  newton_restart_monolithic_ = false;

  // get interface stabilization specific parameters
  coupling_method_    = DRT::INPUT::IntegralValue<INPAR::XFEM::CouplingMethod>(params_xf_stab,"COUPLING_METHOD");

  // set flag if any edge-based fluid stabilization has to integrated as std or gp stabilization
  {
    bool edge_based = (      params_->sublist("RESIDUAL-BASED STABILIZATION").get<std::string>("STABTYPE")=="edge_based"
                             or params_->sublist("EDGE-BASED STABILIZATION").get<std::string>("EOS_PRES")        != "none"
                             or params_->sublist("EDGE-BASED STABILIZATION").get<std::string>("EOS_CONV_STREAM") != "none"
                             or params_->sublist("EDGE-BASED STABILIZATION").get<std::string>("EOS_CONV_CROSS")  != "none"
                             or params_->sublist("EDGE-BASED STABILIZATION").get<std::string>("EOS_DIV")         != "none");

    // set flag if a viscous or transient (1st or 2nd order) ghost-penalty stabiliation due to Nitsche's method has to be integrated
    bool ghost_penalty= (    (bool)DRT::INPUT::IntegralValue<int>(params_xf_stab,"GHOST_PENALTY_STAB")
                          or (bool)DRT::INPUT::IntegralValue<int>(params_xf_stab,"GHOST_PENALTY_TRANSIENT_STAB")
                          or (bool)DRT::INPUT::IntegralValue<int>(params_xf_stab,"GHOST_PENALTY_2nd_STAB") );

    // determine, whether face-based stabilizing terms are active
    eval_eos_ = edge_based || ghost_penalty;
  }

  if(myrank_ == 0)
  {
    std::cout<<"\nVolume:   Gauss point generating method = "<< params_xfem.get<std::string>("VOLUME_GAUSS_POINTS_BY");
    std::cout<<"\nBoundary: Gauss point generating method = "<< params_xfem.get<std::string>("BOUNDARY_GAUSS_POINTS_BY") << "\n\n";
  }

  // set XFEM-related parameters on element level
  SetElementGeneralFluidXFEMParameter();
  SetFaceGeneralFluidXFEMParameter();
}



// -------------------------------------------------------------------
// set general face fluid parameter (BS 06/2014)
// -------------------------------------------------------------------
void FLD::XFluid::SetElementGeneralFluidXFEMParameter()
{

  Teuchos::ParameterList eleparams;

  eleparams.set<int>("action",FLD::set_general_fluid_xfem_parameter); // do not call another action as then another object of the std-class will be created

  //------------------------------------------------------------------------------------------------------
  // set general element parameters
  eleparams.set("form of convective term",convform_);
  eleparams.set<int>("Linearisation",newton_);
  eleparams.set<int>("Physical Type", physicaltype_);

  // parameter for stabilization
  eleparams.sublist("RESIDUAL-BASED STABILIZATION") = params_->sublist("RESIDUAL-BASED STABILIZATION");

  // get function number of given Oseen advective field if necessary
  if (physicaltype_==INPAR::FLUID::oseen)
    eleparams.set<int>("OSEENFIELDFUNCNO", params_->get<int>("OSEENFIELDFUNCNO"));

  //set time integration scheme
  eleparams.set<int>("TimeIntegrationScheme", timealgo_);

  //------------------------------------------------------------------------------------------------------
  // set general parameters for turbulent flow
  eleparams.sublist("TURBULENCE MODEL") = params_->sublist("TURBULENCE MODEL");

  // set model-dependent parameters
  eleparams.sublist("SUBGRID VISCOSITY") = params_->sublist("SUBGRID VISCOSITY");
  eleparams.sublist("MULTIFRACTAL SUBGRID SCALES") = params_->sublist("MULTIFRACTAL SUBGRID SCALES");


  //------------------------------------------------------------------------------------------------------
  // set general XFEM element parameters

  eleparams.sublist("XFEM")                         = params_->sublist("XFEM");
  eleparams.sublist("XFLUID DYNAMIC/GENERAL")       = params_->sublist("XFLUID DYNAMIC/GENERAL");
  eleparams.sublist("XFLUID DYNAMIC/STABILIZATION") = params_->sublist("XFLUID DYNAMIC/STABILIZATION");


  //------------------------------------------------------------------------------------------------------
  // set the params in the XFEM-parameter-list class
  DRT::ELEMENTS::FluidType::Instance().PreEvaluate(*discret_,eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);

  return;
}

// -------------------------------------------------------------------
// set general face fluid parameter (BS 06/2014)
// -------------------------------------------------------------------
void FLD::XFluid::SetFaceGeneralFluidXFEMParameter()
{

  //------------------------------------------------------------------------------------------------------
  // set general fluid stabilization parameter for faces
  {
    Teuchos::ParameterList faceparams;

    faceparams.set<int>("action",FLD::set_general_face_fluid_parameter);

    faceparams.sublist("EDGE-BASED STABILIZATION")     = params_->sublist("EDGE-BASED STABILIZATION");

    faceparams.set<int>("STABTYPE", DRT::INPUT::IntegralValue<INPAR::FLUID::StabType>( params_->sublist("RESIDUAL-BASED STABILIZATION"), "STABTYPE"));

    faceparams.set<int>("Physical Type", physicaltype_);

    // get function number of given Oseen advective field if necessary
    if (physicaltype_==INPAR::FLUID::oseen) faceparams.set<int>("OSEENFIELDFUNCNO", params_->get<int>("OSEENFIELDFUNCNO"));

    DRT::ELEMENTS::FluidIntFaceType::Instance().PreEvaluate(*discret_,faceparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  }

  //------------------------------------------------------------------------------------------------------
  // set XFEM specific parameter for faces
  {
    Teuchos::ParameterList faceparams;

    faceparams.set<int>("action",FLD::set_general_face_xfem_parameter);

    // set general fluid face parameters are contained in the following two sublists
    faceparams.sublist("XFLUID DYNAMIC/STABILIZATION") = params_->sublist("XFLUID DYNAMIC/STABILIZATION");

    DRT::ELEMENTS::FluidIntFaceType::Instance().PreEvaluate(*discret_,faceparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  }

  return;
}


// -------------------------------------------------------------------
// set general time parameter (AE 01/2011)
// -------------------------------------------------------------------
void FLD::XFluid::SetElementTimeParameter()
{

  Teuchos::ParameterList eleparams;

  // set action
  eleparams.set<int>("action",FLD::set_time_parameter);
  // set time integration scheme
  eleparams.set<int>("TimeIntegrationScheme", timealgo_);
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
    eleparams.set<int>("ost cont and press",params_->get<int>("ost cont and press"));
    eleparams.set<bool>("ost new"          , params_->get<bool>("ost new"));
  }

  // call standard loop over elements
  //discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);

  DRT::ELEMENTS::FluidType::Instance().PreEvaluate(*discret_,eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::XFluid::SetGeneralTurbulenceParameters()
{
  FluidImplicitTimeInt::SetGeneralTurbulenceParameters();
  // mark XFEM fluid in name of statistics outputfile (postfix)
  statistics_outfilename_.append("_xfluid");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::XFluid::CreateInitialState()
{
  // initialize the state class iterator with -1
  // the XFluidState class called from the constructor is then indexed with 0
  // all further first cuts of a new time-step have then index 1 and have to be reset to 0 in PrepareTimeStep()
  state_it_=-1;

  // ---------------------------------------------------------------------
  // create the initial state class
  CreateState();
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::XFluid::CreateState()
{
  discret_->Comm().Barrier();
  TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::CreateState" );

  // ---------------------------------------------------------------------
  // create a new state class
  DestroyState();

  // create new state object
  state_ = GetNewState();

  //--------------------------------------------------------------------------------------
  // initialize the KrylovSpaceProjection
  InitKrylovSpaceProjection();

  //--------------------------------------------------------------------------------------
  if (false/*xfluid_.params_->get<bool>("INFNORMSCALING")*/)
  {
    fluid_infnormscaling_ = Teuchos::rcp(new FLD::UTILS::FluidInfNormScaling(*state_->velpressplitter_));
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::XFluid::DestroyState()
{
  if( state_ != Teuchos::null  and state_.strong_count() > 1)
    dserror("deleting old state class object does not work properly, more than one rcp pointer existent!!!");

  if(state_ != Teuchos::null)
  {
    if(!state_->Destroy())
      dserror("destroying XFluidState object failed");

    // delete the old state object and its content (if no ownership given anymore) not to have two objects in memory at the same time
    state_ = Teuchos::null;
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<FLD::XFluidState> FLD::XFluid::GetNewState()
{
  if(state_!=Teuchos::null)
    dserror("please destroy the old state-class before creating a new one!");

  //-------------------------------------------------------------
  // export background mesh ale displacements
  //-------------------------------------------------------------

  // init col vector holding background ALE displacements for backdis
  Teuchos::RCP<Epetra_Vector> dispnpcol = Teuchos::null;

  if (alefluid_)
  {
    dispnpcol = Teuchos::rcp(new Epetra_Vector(*xdiscret_->InitialDofColMap()));
    LINALG::Export(*dispnp_,*dispnpcol);
  }

  //-------------------------------------------------------------
  // create a temporary state-creator object
  //-------------------------------------------------------------
  // create the new state class (vectors, matrices...)
  state_it_++;

  Teuchos::RCP<FLD::XFluidState> state = state_creator_->Create(
      xdiscret_,
      dispnpcol,     //!< col vector holding background ALE displacements for backdis
      solver_->Params(),
      step_,
      time_
      );

  //--------------------------------------------------------------------------------------
  // initialize ALE state vectors
  if(alefluid_) state->InitALEStateVectors(xdiscret_, dispnp_, gridvnp_);

  return state;
}


/*----------------------------------------------------------------------*
 |  evaluate elements, volumecells and boundary cells      schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::AssembleMatAndRHS( int itnum )
{
  discret_->Comm().Barrier();

  TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Evaluate" );

  //----------------------------------------------------------------------
  // set state vectors for cutter discretization
  condition_manager_->SetState();

  //----------------------------------------------------------------------
  // zero state vectors for interface forces based on cutter discretization
  condition_manager_->ZeroStateVectors_FSI();

  //----------------------------------------------------------------------
  // clear the system matrix and related rhs vectors
  state_->ZeroSystemMatrixAndRhs();

  // clear the coupling matrices and rhs vectors
  state_->ZeroCouplingMatricesAndRhs();

  //----------------------------------------------------------------------
  // set general vector values needed by elements
  discret_->ClearState();

  discret_->SetState("hist" ,state_->hist_ );
  discret_->SetState("veln" ,state_->veln_ );
  discret_->SetState("accam",state_->accam_);
  discret_->SetState("scaaf",state_->scaaf_);
  discret_->SetState("scaam",state_->scaam_);

  if (alefluid_)
  {
    discret_->SetState("dispnp", state_->dispnp_);
    discret_->SetState("gridv", state_->gridvnp_);
  }

  SetStateTimInt();


  //----------------------------------------------------------------------
  if (itnum != itemax_)
  {

    //-------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------
    // Evaluate and Assemble Matrices and rhs vectors
    //-------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------

    //-------------------------------------------------------------------------------
    // evaluate and assemble volume integral based terms
    AssembleMatAndRHS_VolTerms();

    //-------------------------------------------------------------------------------
    // evaluate and assemble face-oriented fluid and ghost penalty stabilizations
    AssembleMatAndRHS_FaceTerms(state_->sysmat_, state_->residual_col_, state_->wizard_);

    //-------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------
    // Finalize Matrices and rhs vectors
    //-------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------

    //-------------------------------------------------------------------------------
    // finalize the complete matrix
    // REMARK: for EpetraFECrs matrices Complete() calls the GlobalAssemble() routine to gather entries from all processors
    // and calls a FillComplete for the first run. For further Newton-steps then the optimized FEAssemble routine is used
    // for speedup.
    state_->sysmat_->Complete();

    //-------------------------------------------------------------------------------
    // finalize the coupling matrices
    state_->CompleteCouplingMatricesAndRhs();

    //-------------------------------------------------------------------------------
    // finalize state vectors based on cutter discretization
    condition_manager_->CompleteStateVectors();

    //-------------------------------------------------------------------------------
    // finalize residual vector
    // need to export residual_col to state_->residual_ (row)
    Epetra_Vector res_tmp(state_->residual_->Map(),true);
    Epetra_Export exporter(state_->residual_col_->Map(),res_tmp.Map());
    int err2 = res_tmp.Export(*state_->residual_col_,exporter,Add);
    if (err2) dserror("Export using exporter returned err=%d",err2);

    // add Neumann loads and contributions from evaluate of volume and face integrals
    state_->residual_->Update(1.0, res_tmp, 1.0, *state_->neumann_loads_, 0.0);

    //-------------------------------------------------------------------------------
    // scaling to get true residual vector
    // negative sign to get forces acting on structural side
    // additional residual-scaling to remove the theta*dt-scaling
    state_->trueresidual_->Update(-1.0*ResidualScaling(),*state_->residual_,0.0);

  }

  //-------------------------------------------------------------------------------
  discret_->ClearState();

  condition_manager_->ClearState();

}

void FLD::XFluid::AssembleMatAndRHS_VolTerms()
{
  //----------------------------------------------------------------------
  //TODO: empty eleparams, could be deleted!
  Teuchos::ParameterList eleparams;

  //------------------------------------------------------------
  DRT::AssembleStrategy strategy(0, 0, state_->sysmat_,Teuchos::null,state_->residual_col_,Teuchos::null,Teuchos::null);

  DRT::Element::LocationArray la( 1 );

  //------------------------------------------------------------
  // call standard loop over elements

  // loop over row elements
  const int numrowele = discret_->NumMyRowElements();

  // REMARK: in this XFEM framework the whole evaluate routine uses only row elements
  // and assembles into EpetraFECrs matrix
  // this is baci-unusual but more efficient in all XFEM applications
  for (int i=0; i<numrowele; ++i)
  {
    DRT::Element* actele = discret_->lRowElement(i);
    //Teuchos::RCP<MAT::Material> mat = actele->Material();

    DRT::ELEMENTS::Fluid * ele = dynamic_cast<DRT::ELEMENTS::Fluid *>( actele );
    if ( ele==NULL )
    {
      dserror( "expect fluid element" );
    }

    DRT::ELEMENTS::FluidEleInterface * impl = DRT::ELEMENTS::FluidFactory::ProvideImplXFEM( actele->Shape(), "xfem");

    GEO::CUT::ElementHandle * e = state_->Wizard()->GetElement( actele );

    if ( e!=NULL)
    {
      std::vector< GEO::CUT::plain_volumecell_set > cell_sets;
      std::vector< std::vector<int> > nds_sets;
      std::vector<std::vector< DRT::UTILS::GaussIntegration > > intpoints_sets;

      bool has_xfem_integration_rule =
          e->GetCellSets_DofSets_GaussPoints( cell_sets, nds_sets, intpoints_sets, include_inner_);

      if(cell_sets.size() != nds_sets.size()) dserror("number of cell_sets and nds_sets not equal!");

      int set_counter = 0;

      for( std::vector< GEO::CUT::plain_volumecell_set>::iterator s=cell_sets.begin();
          s!=cell_sets.end();
          s++)
      {

        GEO::CUT::plain_volumecell_set & cells = *s;
        const std::vector<int> & nds = nds_sets[set_counter];

        //Pointer to material of current volume cell
        //Assumes the plain_volumecell_set are all on the same side of the interface.
        Teuchos::RCP<MAT::Material> mat;
        condition_manager_->GetVolumeCellMaterial(actele,mat,*cells.begin());

        // we have to assemble all volume cells of this set
        // for linear elements, there should be only one volume-cell for each set
        // for quadratic elements, there are some volume-cells with respect to subelements, that have to be assembled at once

        // get element location vector, dirichlet flags and ownerships (discret, nds, la, doDirichlet)
        actele->LocationVector(*discret_,nds,la,false);

        // get dimension of element matrices and vectors
        // Reshape element matrices and vectors and init to zero (rdim, cdim)
        strategy.ClearElementStorage( la[0].Size(), la[0].Size() );


        if(!has_xfem_integration_rule) // use standard integration!!!
        {
          //------------------------------------------------------------
          // Evaluate domain integrals
          TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Evaluate 3) standard domain" );

          // call the element evaluate method
          int err = impl->Evaluate( ele, *discret_, la[0].lm_, eleparams, mat,
              strategy.Elematrix1(),
              strategy.Elematrix2(),
              strategy.Elevector1(),
              strategy.Elevector2(),
              strategy.Elevector3()
          );

          if (err)
            dserror("Proc %d: Element %d returned err=%d",discret_->Comm().MyPID(),actele->Id(),err);
        }
        else
        {
          if(cell_sets.size() != intpoints_sets.size()) dserror("number of cell_sets and intpoints_sets not equal!");

          //------------------------------------------------------------
          // Evaluate domain integrals
          TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Evaluate 1) cut domain" );

          // call the element evaluate method
          int err = impl->EvaluateXFEM( ele, *discret_, la[0].lm_, eleparams, mat,
              strategy.Elematrix1(),
              strategy.Elematrix2(),
              strategy.Elevector1(),
              strategy.Elevector2(),
              strategy.Elevector3(),
              intpoints_sets[set_counter],
              cells);

          if (err)
            dserror("Proc %d: Element %d returned err=%d",discret_->Comm().MyPID(),actele->Id(),err);
        }

        //------------------------------------------------------------
        // Evaluate interface integrals
        // do cut interface condition

        // map of sid and corresponding boundary cells ( for quadratic elements: collected via volumecells of subelements)
        std::map<int, std::vector<GEO::CUT::BoundaryCell*> > element_bcells;

        for ( GEO::CUT::plain_volumecell_set::iterator i=cells.begin(); i!=cells.end(); ++i )
        {
          GEO::CUT::VolumeCell * vc = *i;

            vc->GetBoundaryCellsToBeIntegrated( element_bcells );

        }

        //Set material at interface (Master and Slave side)
        Teuchos::RCP<MAT::Material> matptr_m;
        Teuchos::RCP<MAT::Material> matptr_s; //If not instantiated, it is left as null pointer.

        //Get material pointer for master side (LevelSet: positive side)
        condition_manager_->GetInterfaceMasterMaterial(actele,matptr_m,*cells.begin());

        // split the boundary cells by the different mesh couplings / levelset couplings
        // coupling matrices have to be evaluated for each coupling time separtely and cannot be mixed up
        // e.g. do not mix two-phase flow coupling matrices with XFSI coupling matrices
        std::map<int, std::vector<GEO::CUT::BoundaryCell*> > empty_map;
        empty_map.clear();

        const int num_coupling = condition_manager_->NumCoupling();

        // TODO: use a map instead of a vector, see handling of C_sx... matrices in state-class
        std::vector< std::map<int, std::vector<GEO::CUT::BoundaryCell*> > > coupling_bcells(num_coupling, empty_map);

        for ( std::map<int, std::vector<GEO::CUT::BoundaryCell*> >::const_iterator bc=element_bcells.begin();
            bc!=element_bcells.end(); ++bc )
        {
          int coup_sid = bc->first; // all boundary cells within the current iterator belong to the same side

          const int coup_idx = condition_manager_->GetCouplingIndex(coup_sid, actele->Id());

          std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells = coupling_bcells[coup_idx];

          std::vector<GEO::CUT::BoundaryCell*> & bc_new = bcells[bc->first];
          bc_new.clear();
          std::copy(bc->second.begin(), bc->second.end(), std::inserter(bc_new, bc_new.end()));

          const std::vector<std::pair<int,int> > cloning_information = condition_manager_->GetBCCloneInformation(coup_sid, actele->Id(),coup_idx);
          for (uint clone_id = 0; clone_id < cloning_information.size(); ++clone_id)
          {
//            std::cout << "XFluid - Cloning News: " << coup_idx << " --> " << cloning_information.first << ", " <<
//                coup_sid << " --> " << cloning_information.second << std::endl;
            std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells = coupling_bcells[cloning_information[clone_id].first];

            std::vector<GEO::CUT::BoundaryCell*> & bc_new = bcells[cloning_information[clone_id].second];
            bc_new.clear();
            std::copy(bc->second.begin(), bc->second.end(), std::inserter(bc_new, bc_new.end()));
          }
        }

        // loop all the different couplings
        for(int coupl_idx=0; coupl_idx< num_coupling; coupl_idx++)
        {
          std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells = coupling_bcells[coupl_idx];
          std::map<int, std::vector<DRT::UTILS::GaussIntegration> > bintpoints;

          // for each side that is involved in the cut for this element,
          // the coupling matrices C_fs_, C_sf_ and the rhs_s has to be built
          std::map<int, std::vector<Epetra_SerialDenseMatrix> > side_coupling;

          if ( bcells.size() > 0 )
          {
            TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Evaluate 2) interface" );

            // get boundary cell Gaussian points
            e->BoundaryCellGaussPointsLin( bcells, bintpoints);

            //-----------------------------------------------------------
            // fluid-structure coupling part

            std::map<int, std::vector<int> > patchcouplm; // lm vector for each element/side which couples with the current bg element
            std::vector<int> patchelementslm;             // dofs of all coupling elements which couple with the current bg element

            // initialize the coupling lm vectors for each coupling side
            for ( std::map<int,  std::vector<GEO::CUT::BoundaryCell*> >::const_iterator bc=bcells.begin();
                bc!=bcells.end(); ++bc )
            {
              int coup_sid = bc->first; // all boundary cells within the current iterator belong to the same side

              // Set material for coupling element
              // Get slave material from the condition.
              condition_manager_->GetInterfaceSlaveMaterial(actele,matptr_s,coup_sid);

              // boundary discretization for mesh coupling and background discretization for level-set coupling
              Teuchos::RCP<DRT::Discretization> coupl_dis = condition_manager_->GetCouplingDis( coup_sid );

              std::vector<int> & patchlm = patchcouplm[coup_sid]; // []-operator creates new vector, dofs of current coupling side

              // get dofs for coupling side or coupling element
              if(condition_manager_->IsMeshCoupling(coup_sid))
              {
                // fill patchlm for the element we couple with
                condition_manager_->GetCouplingEleLocationVector(coup_sid,patchlm);
              }
              else if(condition_manager_->IsLevelSetCoupling(coup_sid))
              {
                if(!condition_manager_->IsCoupling( coup_sid, ele->Id() )) continue; // level-set wdbc case

                // get the other nds-set which is connected to the current one via this boundary-cell
                DRT::Element::LocationArray la_other( 1 );

                if(bc->second.empty()) dserror("no boundary cells stored!");

                GEO::CUT::BoundaryCell* boundcell = bc->second[0]; // first boundary-cell
                GEO::CUT::Facet * f = boundcell->GetFacet();

                const GEO::CUT::plain_volumecell_set & vcs = f->Cells();
                if(vcs.size() != 2) dserror("for the given boundary-cells facet, exactly two volume-cells have to be adjacent!");

                std::vector<int> nds_other;

                for(GEO::CUT::plain_volumecell_set::const_iterator it= vcs.begin(); it!=vcs.end(); it++)
                {
                  if((*it)->Position()==GEO::CUT::Point::inside ) // now take the inside volume-cell
                  {
                    nds_other = (*it)->NodalDofSet();
                    break;
                  }
                }

                if((*cells.begin())->Position() == GEO::CUT::Point::inside )
                  dserror("For a two-sided level set coupling, we should not enter here with inside volume-cells!!!");

                // get element location vector, dirichlet flags and ownerships (discret, nds, la, doDirichlet)
                actele->LocationVector(*coupl_dis,nds_other,la_other,false);
                std::copy( la_other[0].lm_.begin(), la_other[0].lm_.end(), std::inserter(patchlm,patchlm.end()));
              }

              // initialize the coupling matrices for each coupling side and the current element
              if(condition_manager_->IsCoupling( coup_sid, ele->Id() ))
              {
                patchelementslm.reserve( patchelementslm.size() + patchlm.size());
                patchelementslm.insert(patchelementslm.end(), patchlm.begin(), patchlm.end());

                const size_t ndof_i = patchlm.size();     // number of dofs of this coupling sides
                const size_t ndof   = la[0].lm_.size();   // number of dofs for background element

                std::vector<Epetra_SerialDenseMatrix> & couplingmatrices = side_coupling[coup_sid]; // the function inserts a new element with that key and returns a reference to its mapped value
                if ( couplingmatrices.size()!=0 )
                  dserror("zero sized vector expected");

                couplingmatrices.resize(3);

                // no coupling for pressure in stress based method, but the coupling matrices include entries for pressure coupling
                couplingmatrices[0].Shape(ndof_i,ndof);  //C_sf = C_uiu
                couplingmatrices[1].Shape(ndof,ndof_i);  //C_fs = C_uui
                couplingmatrices[2].Shape(ndof_i,1);     //rhC_s = rhs_ui
              } // IsCoupling
            } // loop bcs

            const size_t nui = patchelementslm.size(); // sum over number of dofs of all coupling sides
            Epetra_SerialDenseMatrix C_ss(nui,nui);    // coupling matrix for monolithic fluid-structure interaction, struct-struct couplings between different sides

            {
              TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Evaluate 2) interface (only evaluate)" );

              if( CouplingMethod() == INPAR::XFEM::Hybrid_LM_Cauchy_stress or
                  CouplingMethod() == INPAR::XFEM::Hybrid_LM_viscous_stress)
                impl->ElementXfemInterfaceHybridLM(
                    ele,
                    *discret_,
                    la[0].lm_,
                    condition_manager_,
                    intpoints_sets[set_counter],
                    bcells,
                    bintpoints,
                    patchcouplm,
                    side_coupling,
                    eleparams,
                    mat,
                    strategy.Elematrix1(),
                    strategy.Elevector1(),
                    C_ss,
                    cells
                );

              if(CouplingMethod() == INPAR::XFEM::Nitsche)
                impl->ElementXfemInterfaceNIT(
                    ele,
                    *discret_,
                    la[0].lm_,
                    condition_manager_,
                    bcells,
                    bintpoints,
                    patchcouplm,
                    eleparams,
                    matptr_m,
                    matptr_s,
                    strategy.Elematrix1(),
                    strategy.Elevector1(),
                    cells,
                    side_coupling,
                    C_ss
                );

            }

            //------------------------------------------------------------------------------------------
            // Assemble bgele-side coupling matrices for monolithic fluid-structure interaction
            //------------------------------------------------------------------------------------------

            Teuchos::RCP<XFluidState::CouplingState> & coup_state = state_->coup_state_[coupl_idx];

            for ( std::map<int, std::vector<Epetra_SerialDenseMatrix> >::const_iterator sc=side_coupling.begin();
                sc!=side_coupling.end(); ++sc )
            {
              std::vector<Epetra_SerialDenseMatrix>  couplingmatrices = sc->second;
              int coup_sid = sc->first;

              std::vector<int> & patchlm = patchcouplm[coup_sid];

              // assemble C_sf_ = Cuiu
              // create a dummy mypatchlmowner that assembles also non-local rows and communicates the required data
              std::vector<int> mypatchlmowner(patchlm.size(), myrank_);
              {
                TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Evaluate 6) FEAssemble" );
                coup_state->C_sx_->FEAssemble( couplingmatrices[0],patchlm,mypatchlmowner,la[0].lm_);
              }

              // assemble C_fs_ = Cuui
              std::vector<int> mylmowner(la[0].lmowner_.size(), myrank_);
              {
                TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Evaluate 6) FEAssemble" );
                coup_state->C_xs_->FEAssemble( couplingmatrices[1],la[0].lm_,mylmowner, patchlm);
              }

              // assemble rhC_s_col = rhC_ui_col
              Epetra_SerialDenseVector rhC_s_eptvec(::View,couplingmatrices[2].A(),patchlm.size());
              LINALG::Assemble(*(coup_state->rhC_s_col_), rhC_s_eptvec, patchlm, mypatchlmowner);
            }

            if(!side_coupling.empty()) // at least one side contributed to coupling for this element
            {
              // assemble C_ss_ = Cuiui
              std::vector<int> mypatchelementslmowner(patchelementslm.size(), myrank_);
              coup_state->C_ss_->FEAssemble(C_ss, patchelementslm, mypatchelementslmowner, patchelementslm );
            }

          } // bcells.size() > 0
        } // loop coupl index
        //------------------------------------------------------------
        // Assemble matrix and vectors

        // introduce an vector containing the rows for that values have to be communicated
        // REMARK: when assembling row elements also non-row rows have to be communicated
        std::vector<int> myowner(la[0].lmowner_.size(), strategy.Systemvector1()->Comm().MyPID());
        {
          TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Evaluate 6) FEAssemble" );
          // calls the Assemble function for EpetraFECrs matrices including communication of non-row entries
          state_->sysmat_->FEAssemble( strategy.Elematrix1(), la[0].lm_,myowner,la[0].lm_);
        }
        // REMARK:: call Assemble without lmowner
        // to assemble the residual_col vector on only row elements also column nodes have to be assembled
        // do not exclude non-row nodes (modify the real owner to myowner)
        // after assembly the col vector it has to be exported to the row residual_ vector
        // using the 'Add' flag to get the right value for shared nodes
        LINALG::Assemble(*strategy.Systemvector1(),strategy.Elevector1(),la[0].lm_,myowner);

        set_counter += 1;

      } // end of loop over cellsets // end of assembly for each set of cells
    } // end of if(e!=NULL) // assembly for cut elements
    else
    {

      Teuchos::RCP<MAT::Material> mat = actele->Material();

      if(mat->MaterialType()==INPAR::MAT::m_matlist)
        dserror("No matlists allowed here!!");

      // get element location vector, dirichlet flags and ownerships
      actele->LocationVector(*discret_,la,false);

      // get dimension of element matrices and vectors
      // Reshape element matrices and vectors and init to zero
      strategy.ClearElementStorage( la[0].Size(), la[0].Size() );

      {
        TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Evaluate 3) standard domain" );

        // call the element evaluate method
        int err = impl->Evaluate(
            ele, *discret_, la[0].lm_, eleparams, mat,
            strategy.Elematrix1(),
            strategy.Elematrix2(),
            strategy.Elevector1(),
            strategy.Elevector2(),
            strategy.Elevector3() );

        if (err) dserror("Proc %d: Element %d returned err=%d",discret_->Comm().MyPID(),actele->Id(),err);
      }

      // introduce an vector containing the rows for that values have to be communicated
      // REMARK: when assembling row elements also non-row rows have to be communicated
      std::vector<int> myowner(la[0].lmowner_.size(), strategy.Systemvector1()->Comm().MyPID());
      {
        TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Evaluate 6) FEAssemble" );

        // calls the Assemble function for EpetraFECrs matrices including communication of non-row entries
        state_->sysmat_->FEAssemble( strategy.Elematrix1(), la[0].lm_,myowner,la[0].lm_);
      }

      // REMARK:: call Assemble without lmowner
      // to assemble the residual_col vector on only row elements also column nodes have to be assembled
      // do not exclude non-row nodes (modify the real owner to myowner)
      // after assembly the col vector it has to be exported to the row residual_ vector
      // using the 'Add' flag to get the right value for shared nodes
      LINALG::Assemble(*strategy.Systemvector1(),strategy.Elevector1(),la[0].lm_,myowner);

    }
  } // loop row elements

} //AssembleMatAndRHS_VolTerms




void FLD::XFluid::AssembleMatAndRHS_FaceTerms(
    const Teuchos::RCP<LINALG::SparseMatrix> & sysmat,
    const Teuchos::RCP<Epetra_Vector> & residual_col,
    const Teuchos::RCP<GEO::CutWizard> & wizard,
    bool is_ghost_penalty_reconstruct
    )
{
  // call edge stabilization
  if( eval_eos_ || is_ghost_penalty_reconstruct)
  {
    TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Evaluate 4) EOS" );

    Teuchos::ParameterList faceparams;

    // set additional faceparams according to ghost-penalty terms due to Nitsche's method
    faceparams.set("ghost_penalty_reconstruct", is_ghost_penalty_reconstruct); // no XFEM timeintegration reconstruction call

    //------------------------------------------------------------
    // loop over row faces

    Teuchos::RCP<DRT::DiscretizationFaces> xdiscret = Teuchos::rcp_dynamic_cast<DRT::DiscretizationFaces>(discret_, true);

    const int numrowintfaces = xdiscret->NumMyRowFaces();

    // REMARK: in this XFEM framework the whole evaluate routine uses only row internal faces
    // and assembles into EpetraFECrs matrix
    // this is baci-unusual but more efficient in all XFEM applications
    for (int i=0; i<numrowintfaces; ++i)
    {
      DRT::Element* actface = xdiscret->lRowFace(i);

      DRT::ELEMENTS::FluidIntFace * face_ele = dynamic_cast<DRT::ELEMENTS::FluidIntFace *>( actface );
      if ( face_ele==NULL ) dserror( "expect FluidIntFace element" );

      bool gmsh_EOS_out(DRT::INPUT::IntegralValue<int>(params_->sublist("XFEM"),"GMSH_EOS_OUT"));
      edgestab_->EvaluateEdgeStabGhostPenalty(
          faceparams, discret_, face_ele, sysmat, residual_col, wizard, include_inner_, gmsh_EOS_out);
    }
  }
}


/*----------------------------------------------------------------------*
 | integrate shape functions over domain                   schott 12/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::IntegrateShapeFunction(
  Teuchos::ParameterList & eleparams,
  DRT::Discretization & discret,
  Teuchos::RCP<Epetra_Vector> vec
  )
{
  TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::IntegrateShapeFunction" );

  // create an column vector for assembly over row elements that has to be communicated at the end
  Teuchos::RCP<Epetra_Vector> w_col = LINALG::CreateVector(*discret.DofColMap(),true);


  //----------------------------------------------------------------------

  // call standard loop over elements

  DRT::AssembleStrategy strategy(0, 0, Teuchos::null,Teuchos::null,w_col,Teuchos::null,Teuchos::null);

  DRT::Element::LocationArray la( 1 );


  //------------------------------------------------------------
  // loop over row elements
  const int numrowele = discret.NumMyRowElements();

  // REMARK: in this XFEM framework the whole evaluate routine uses only row elements
  // and assembles into EpetraFECrs matrix
  // this is baci-unusual but more efficient in all XFEM applications
  for (int i=0; i<numrowele; ++i)
  {
    DRT::Element* actele = discret.lRowElement(i);
    Teuchos::RCP<MAT::Material> mat = actele->Material();

    DRT::ELEMENTS::Fluid * ele = dynamic_cast<DRT::ELEMENTS::Fluid *>( actele );
    if ( ele==NULL )
    {
      dserror( "expect fluid element" );
    }

    DRT::ELEMENTS::FluidEleInterface * impl = DRT::ELEMENTS::FluidFactory::ProvideImplXFEM( actele->Shape(), "xfem");

    GEO::CUT::ElementHandle * e = state_->Wizard()->GetElement( actele );


    if ( e!=NULL )
    {

      std::vector< GEO::CUT::plain_volumecell_set > cell_sets;
      std::vector< std::vector<int> > nds_sets;
      std::vector<std::vector< DRT::UTILS::GaussIntegration > > intpoints_sets;

      bool has_xfem_integration_rule =
          e->GetCellSets_DofSets_GaussPoints( cell_sets, nds_sets, intpoints_sets, false); //(include_inner=false)

      if(cell_sets.size() != nds_sets.size()) dserror("number of cell_sets and nds_sets not equal!");

      int set_counter = 0;

      for( std::vector< GEO::CUT::plain_volumecell_set>::iterator s=cell_sets.begin();
          s!=cell_sets.end();
          s++)
      {
        GEO::CUT::plain_volumecell_set & cells = *s;
        const std::vector<int> & nds = nds_sets[set_counter];

        // we have to assemble all volume cells of this set
        // for linear elements, there should be only one volumecell for each set
        // for quadratic elements, there are some volumecells with respect to subelements, that have to be assembled at once


        // get element location vector, dirichlet flags and ownerships (discret, nds, la, doDirichlet)
        actele->LocationVector(discret,nds,la,false);

        // get dimension of element matrices and vectors
        // Reshape element matrices and vectors and init to zero (rdim, cdim)
        strategy.ClearElementStorage( la[0].Size(), la[0].Size() );

        if(!has_xfem_integration_rule)
        {
          // call the element evaluate method
          Epetra_SerialDenseMatrix elemat1;
          Epetra_SerialDenseMatrix elemat2;
          Epetra_SerialDenseVector elevec2;
          Epetra_SerialDenseVector elevec3;
          Teuchos::ParameterList params;
          params.set<int>("action",FLD::integrate_shape);
          Teuchos::RCP<MAT::Material> mat = ele->Material();
          int err = impl->EvaluateService( ele, params, mat, discret, la[0].lm_, elemat1, elemat2, strategy.Elevector1(), elevec2, elevec3 );

          if (err)
            dserror("Proc %d: Element %d returned err=%d",discret.Comm().MyPID(),actele->Id(),err);
        }
        else
        {
          if(cell_sets.size() != intpoints_sets.size()) dserror("number of cell_sets and intpoints_sets not equal!");

          //------------------------------------------------------------
          // Evaluate domain integrals
          TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Evaluate 1) cut domain" );

          // call the element evaluate method
          int err = impl->IntegrateShapeFunctionXFEM(
              ele, discret, la[0].lm_, strategy.Elevector1(),
              intpoints_sets[set_counter],
              cells
          );

          if (err)
            dserror("Proc %d: Element %d returned err=%d",discret.Comm().MyPID(),actele->Id(),err);
        }


        //------------------------------------------------------------
        // Assemble vector

        // introduce an vector containing the rows for that values have to be communicated
        // REMARK: when assembling row elements also non-row rows have to be communicated
        std::vector<int> myowner;
        for(size_t index=0; index<la[0].lmowner_.size(); index++)
        {
          myowner.push_back(strategy.Systemvector1()->Comm().MyPID());
        }

        // REMARK:: call Assemble without lmowner
        // to assemble the residual_col vector on only row elements also column nodes have to be assembled
        // do not exclude non-row nodes (modify the real owner to myowner)
        // after assembly the col vector it has to be exported to the row residual_ vector
        // using the 'Add' flag to get the right value for shared nodes
        LINALG::Assemble(*strategy.Systemvector1(),strategy.Elevector1(),la[0].lm_,myowner);

        set_counter += 1;

      } // end of loop over cellsets // end of assembly for each set of cells
    } // end of if(e!=NULL) // assembly for cut elements
    else
    {
      TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Evaluate 3) standard domain" );

      // get element location vector, dirichlet flags and ownerships
      actele->LocationVector(discret,la,false);

      // get dimension of element matrices and vectors
      // Reshape element matrices and vectors and init to zero
      strategy.ClearElementStorage( la[0].Size(), la[0].Size() );

      // call the element evaluate method
      Epetra_SerialDenseMatrix elemat1;
      Epetra_SerialDenseMatrix elemat2;
      Epetra_SerialDenseVector elevec2;
      Epetra_SerialDenseVector elevec3;
      Teuchos::ParameterList params;
      params.set<int>("action",FLD::integrate_shape);
      Teuchos::RCP<MAT::Material> mat = ele->Material();
      int err = impl->EvaluateService( ele, params, mat, discret, la[0].lm_, elemat1, elemat2, strategy.Elevector1(), elevec2, elevec3 );

      if (err) dserror("Proc %d: Element %d returned err=%d",discret.Comm().MyPID(),actele->Id(),err);

      // introduce an vector containing the rows for that values have to be communicated
      // REMARK: when assembling row elements also non-row rows have to be communicated
      std::vector<int> myowner;
      for(size_t index=0; index<la[0].lmowner_.size(); index++)
      {
        myowner.push_back(strategy.Systemvector1()->Comm().MyPID());
      }

      // REMARK:: call Assemble without lmowner
      // to assemble the residual_col vector on only row elements also column nodes have to be assembled
      // do not exclude non-row nodes (modify the real owner to myowner)
      // after assembly the col vector it has to be exported to the row w_ vector
      // using the 'Add' flag to get the right value for shared nodes
      LINALG::Assemble(*strategy.Systemvector1(),strategy.Elevector1(),la[0].lm_,myowner);

    }


  }

  discret.ClearState();


  //-------------------------------------------------------------------------------
  // need to export residual_col to systemvector1 (residual_)
  Epetra_Vector vec_tmp(vec->Map(),false);
  Epetra_Export exporter(strategy.Systemvector1()->Map(),vec_tmp.Map());
  int err2 = vec_tmp.Export(*strategy.Systemvector1(),exporter,Add);
  if (err2) dserror("Export using exporter returned err=%d",err2);
  vec->Scale(1.0,vec_tmp);
}


/*----------------------------------------------------------------------*
 |  evaluate gradient penalty terms to reconstruct ghost values  schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::AssembleMatAndRHS_GradientPenalty(
    Teuchos::RCP<LINALG::MapExtractor>         ghost_penaly_dbcmaps,
    Teuchos::RCP<LINALG::SparseMatrix> sysmat_gp,
    Teuchos::RCP<Epetra_Vector> residual_gp,
    Teuchos::RCP<Epetra_Vector> vec
)
{
  TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::AssembleMatAndRHS_GradientPenalty" );

  // create a new sysmat with reusing the old graph (without the DBC modification) when savegraph-flag is switched on
  // for the first iteration we need to create a new matrix without reusing the graph as the matrix could have been used
  // for another assembly

  //TODO: check if this is necessary or worse!
//  sysmat_gp->Zero()

  residual_gp->PutScalar(0.0);
  Teuchos::RCP<Epetra_Vector> residual_gp_col = LINALG::CreateVector(*state_->xfluiddofcolmap_, true);

  //----------------------------------------------------------------------
  // set general vector values needed by elements
  discret_->ClearState();

  if (alefluid_)
  {
    //dserror("which vectors have to be set for gradient penalty for timeintegration in alefluid?!");
    //In principle we would not need gridv, as tau is anyway set to 1.0 at the end ...
    discret_->SetState("dispnp", state_->dispnp_);
    discret_->SetState("gridv", state_->gridvnp_);
  }

  // set scheme-specific element parameters and vector values
  discret_->SetState("velaf",vec);



  //----------------------------------------------------------------------

  // call loop over face-elements
  AssembleMatAndRHS_FaceTerms(sysmat_gp, residual_gp_col, state_->wizard_, true);

  discret_->ClearState();


  //----------------------------------------------------------------------

  // insert already dummy ones such that Complete does not clear the memory for all rows
  // for which not ghost-penalty term has been assembled
  // for these rows we later have to assemble ones, as we solve for the whole vector

  const Epetra_Map & dbctoggle = *(ghost_penaly_dbcmaps->CondMap());

  bool diagonalblock = true;

  for (int i=0; i<sysmat_gp->EpetraMatrix()->NumMyRows(); ++i)
  {
    int row = sysmat_gp->EpetraMatrix()->GRID(i);

    // check if there is already a value set, otherwise set at least a diagonal entry
    if (dbctoggle.MyGID(row))
    {
      if (diagonalblock)
      {
        double v = 1.0;
#ifdef DEBUG
        int err = sysmat_gp->EpetraMatrix()->InsertGlobalValues(row,1,&v,&row);
        if (err<0) dserror("Epetra_CrsMatrix::InsertGlobalValues returned err=%d",err);
#else
        sysmat_gp->EpetraMatrix()->InsertGlobalValues(row,1,&v,&row);
#endif
      }
    }
  }

  //-------------------------------------------------------------------------------
  // need to export residual_col to systemvector1 (residual_)
  Epetra_Vector res_tmp(residual_gp->Map(),false);
  Epetra_Export exporter(residual_gp_col->Map(),res_tmp.Map());
  int err2 = res_tmp.Export(*residual_gp_col,exporter,Add);
  if (err2) dserror("Export using exporter returned err=%d",err2);
  residual_gp->Update(1.0,res_tmp,1.0);

  //-------------------------------------------------------------------------------
  // finalize the complete matrix
  // REMARK: for EpetraFECrs matrices Complete() calls the GlobalAssemble() routine to gather entries from all processors
  sysmat_gp->Complete();

  return;
}


Teuchos::RCP<Epetra_Vector> FLD::XFluid::StdVelnp()
{
     Teuchos::RCP<Epetra_Vector> initvec = Teuchos::rcp(new Epetra_Vector(*xdiscret_->InitialDofRowMap(),true));
     LINALG::Export(*(state_->velnp_), *initvec);
     return initvec;
}

Teuchos::RCP<Epetra_Vector> FLD::XFluid::StdVeln()
{
     Teuchos::RCP<Epetra_Vector> initvec = Teuchos::rcp(new Epetra_Vector(*xdiscret_->InitialDofRowMap(),true));
     LINALG::Export(*(state_->veln_), *initvec);
     return initvec;
}


/*----------------------------------------------------------------------*
 |  Evaluate errors compared to an analytical solution     schott 09/12 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::vector<double> > FLD::XFluid::EvaluateErrorComparedToAnalyticalSol()
{

  TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::EvaluateErrorComparedToAnalyticalSol" );

  // this functions provides a general implementation for calculating error norms between computed solutions
  // and an analytical solution which is implemented or given by a function in the input file

  // how is the analytical solution available (implemented of via function?)
  INPAR::FLUID::CalcError calcerr = DRT::INPUT::get<INPAR::FLUID::CalcError>(*params_,"calculate error");

  if(calcerr != INPAR::FLUID::no_error_calculation)
  {

    // define the norms that have to be computed

    //-------------------------------------------------------------------------------------------------------------------
    // domain error norms w.r.t incompressible Navier-Stokes/ Oseen equations
    //
    // standard domain errors
    // 1.   || u - u_h ||_L2(Omega)              =   standard L2-norm for velocity
    // 2.   || grad( u - u_h ) ||_L2(Omega)      =   standard H1-seminorm for velocity
    // 3.   || u - u_h ||_H1(Omega)              =   standard H1-norm for velocity
    //                                           =   sqrt( || u - u_h ||^2_L2(Omega) + || grad( u - u_h ) ||^2_L2(Omega) )
    // 4.   || p - p_h ||_L2(Omega)              =   standard L2-norm for pressure
    //
    // viscosity-scaled domain errors
    // 5.   || nu^(+1/2) grad( u - u_h ) ||_L2(Omega)      =   visc-scaled H1-seminorm for velocity
    //                                                     =   nu^(+1/2) * || grad( u - u_h ) ||_L2(Omega) (for homogeneous visc)
    // 6.   || nu^(-1/2) (p - p_h) ||_L2(Omega)            =   visc-scaled L2-norm for pressure
    //                                                     =   nu^(-1/2) * || p - p_h ||_L2(Omega) (for homogeneous visc)
    // 7.   || sigma^(+1/2) ( u - u_h ) ||_L2(Omega)       =   sigma-scaled L2-norm for velocity
    //                                                     =   sigma^(+1/2) * || u - u_h ||_L2(Omega) (for homogeneous sigma)
    // 8.   || Phi^(+1/2) (p - p_h) ||_L2(Omega)           =   Phi-scaled L2-norm for pressure
    //                                                     =   Phi^(+1/2) * || p - p_h ||_L2(Omega) (for homogeneous Phi)
    // with Phi^{-1} = sigma*CP^2 + |beta|*CP + nu + (|beta|*CP/sqrt(sigma*CP^2 + nu))^2, see Massing,Schott,Wall Oseen paper
    //
    // 9. functional G=sin(x)( u,x - u,x exact ) (Sudhakar)
    //
    //
    //-------------------------------------------------------------------------------------------------------------------
    // interface/boundary error norms at the XFEM-interface, boundary
    // w.r.t Nitsche's method to enforce interface/boundary conditions
    //
    // 1.   || nu^(+1/2) (u - u*) ||_H1/2(Gamma)             =  broken H1/2 Sobolev norm for boundary/coupling condition
    // 2.   || nu^(+1/2) grad( u - u_h )*n ||_H-1/2(Gamma)   =  standard H-1/2 Sobolev norm for normal flux (velocity part)
    // 3.   || nu^(-1/2) (p - p_h)*n ||_H-1/2(Gamma)         =  standard H-1/2 Sobolev norm for normal flux (pressure part)
    // 4.   || (u*n)_inflow (u - u*) ||_L2(Gamma)            =  L^2 Sobolev norm for inflow boundary/coupling condition
    // 5.   || (sigma*h+|u|+nu/h)^(+1/2) (u - u*)*n ||_L2(Gamma) =  L^2 Sobolev norm for mass conservation coupling condition
    //
    //-------------------------------------------------------------------------------------------------------------------
    // errors introduced by stabilizations (edge-based fluid stabilizations and ghost-penalty stabilizations)
    //
    // ...
    //-------------------------------------------------------------------------------------------------------------------

    // number of norms that have to be calculated
    const int num_dom_norms    = 10;
    const int num_interf_norms = 8;
    const int num_stab_norms   = 3;
    Teuchos::RCP<Epetra_SerialDenseVector> glob_dom_norms    = Teuchos::rcp(new Epetra_SerialDenseVector(num_dom_norms));
    Teuchos::RCP<Epetra_SerialDenseVector> glob_interf_norms = Teuchos::rcp(new Epetra_SerialDenseVector(num_interf_norms));
    Teuchos::RCP<Epetra_SerialDenseVector> glob_stab_norms   = Teuchos::rcp(new Epetra_SerialDenseVector(num_stab_norms));

    ComputeErrorNorms(glob_dom_norms,glob_interf_norms,glob_stab_norms);

    // standard domain errors
    double dom_err_vel_L2      = 0.0;            //  || u - u_h ||_L2(Omega)              =   standard L2-norm for velocity
    double dom_err_vel_H1_semi = 0.0;            //  || grad( u - u_h ) ||_L2(Omega)      =   standard H1-seminorm for velocity
    double dom_err_vel_H1      = 0.0;            //  || u - u_h ||_H1(Omega)              =   standard H1-norm for velocity
    double dom_err_pre_L2      = 0.0;            //  || p - p_h ||_L2(Omega)              =   standard L2-norm for pressure

    // sigma-,viscosity-scaled domain errors
    double dom_err_vel_H1_semi_nu_scaled = 0.0;  //  || nu^(+1/2) grad( u - u_h ) ||_L2(Omega)      =   visc-scaled H1-seminorm for velocity
    double dom_err_pre_L2_nu_scaled      = 0.0;  //  || nu^(-1/2) (p - p_h) ||_L2(Omega)            =   visc-scaled L2-norm for pressure
    double dom_err_vel_L2_sigma_scaled   = 0.0;  //  || sigma^(+1/2) ( u - u_h ) ||_L2(Omega)       =   sigma-scaled L2-norm for velocity
    double dom_err_pre_L2_Phi_scaled     = 0.0;  //  || Phi^(+1/2) (p - p_h) ||_L2(Omega)           =   Phi-scaled L2-norm for pressure

    // sudhakar functional for testing integration
    double functional = 0.0;

    // interface errors
    double interf_err_Honehalf    = 0.0;         //  || nu^(+1/2) (u - u*) ||_H1/2(Gamma)             =  broken H1/2 Sobolev norm for boundary/coupling condition
    double interf_err_Hmonehalf_u = 0.0;         //  || nu^(+1/2) grad( u - u_h )*n ||_H-1/2(Gamma)   =  broken H-1/2 Sobolev norm for normal flux (velocity part)
    double interf_err_Hmonehalf_p = 0.0;         //  || nu^(-1/2) (p - p_h)*n ||_H-1/2(Gamma)         =  broken H-1/2 Sobolev norm for normal flux (pressure part)
    double interf_err_inflow      = 0.0;         //  || (u*n)_inflow (u - u*) ||_L2(Gamma)            =  L^2 Sobolev norm for inflow boundary/coupling condition
    double interf_err_mass_cons   = 0.0;         //  || (sigma*h+|u|+nu/h)^(+1/2) (u - u*)*n ||_L2(Gamma) =  L^2 Sobolev norm for mass conservation coupling condition


    dom_err_vel_L2             = sqrt((*glob_dom_norms)[0]);
    dom_err_vel_H1_semi        = sqrt((*glob_dom_norms)[1]);
    dom_err_vel_H1             = sqrt((*glob_dom_norms)[2]);
    dom_err_pre_L2             = sqrt((*glob_dom_norms)[3]);

    dom_err_vel_H1_semi_nu_scaled = sqrt((*glob_dom_norms)[4]);
    dom_err_pre_L2_nu_scaled      = sqrt((*glob_dom_norms)[5]);
    dom_err_vel_L2_sigma_scaled   = sqrt((*glob_dom_norms)[6]);
    dom_err_pre_L2_Phi_scaled     = sqrt((*glob_dom_norms)[7]);

    functional                    = (*glob_dom_norms)[8];

    interf_err_Honehalf           = sqrt((*glob_interf_norms)[0]);
    interf_err_Hmonehalf_u        = sqrt((*glob_interf_norms)[1]);
    interf_err_Hmonehalf_p        = sqrt((*glob_interf_norms)[2]);
    interf_err_inflow             = sqrt((*glob_interf_norms)[3]);
    interf_err_mass_cons          = sqrt((*glob_interf_norms)[4]);

    if (myrank_ == 0)
    {
      {
        std::cout.precision(8);
        IO::cout << IO::endl << "---- error norm for analytical solution Nr. "
             <<  DRT::INPUT::get<INPAR::FLUID::CalcError>(*params_,"calculate error")
             <<  " ----------" << IO::endl;
        IO::cout << "-------------- domain error norms -----------------------"       << IO::endl;
        IO::cout << "|| u - u_h ||_L2(Omega)                               =  " << dom_err_vel_L2                     << IO::endl;
        IO::cout << "|| grad( u - u_h ) ||_L2(Omega)                       =  " << dom_err_vel_H1_semi                << IO::endl;
        IO::cout << "|| u - u_h ||_H1(Omega)                               =  " << dom_err_vel_H1                     << IO::endl;
        IO::cout << "|| p - p_h ||_L2(Omega)                               =  " << dom_err_pre_L2                     << IO::endl;
        IO::cout << "---------sigma-,viscosity-scaled domain error norms -----"       << IO::endl;
        IO::cout << "|| nu^(+1/2) grad( u - u_h ) ||_L2(Omega)             =  " << dom_err_vel_H1_semi_nu_scaled      << IO::endl;
        IO::cout << "|| nu^(-1/2) (p - p_h) ||_L2(Omega)                   =  " << dom_err_pre_L2_nu_scaled           << IO::endl;
        IO::cout << "|| sigma^(+1/2) ( u - u_h ) ||_L2(Omega)              =  " << dom_err_vel_L2_sigma_scaled        << IO::endl;
        IO::cout << "|| Phi^(+1/2) (p - p_h) ||_L2(Omega)                  =  " << dom_err_pre_L2_Phi_scaled          << IO::endl;
        IO::cout << "---------------------------------------------------------"       << IO::endl;
        IO::cout << "-------------- interface/boundary error norms -----------"       << IO::endl;
        IO::cout << "|| nu^(+1/2) (u - u*) ||_H1/2(Gamma)                  =  " << interf_err_Honehalf                << IO::endl;
        IO::cout << "|| nu^(+1/2) grad( u - u_h )*n ||_H-1/2(Gamma)        =  " << interf_err_Hmonehalf_u             << IO::endl;
        IO::cout << "|| nu^(-1/2) (p - p_h)*n ||_H-1/2(Gamma)              =  " << interf_err_Hmonehalf_p             << IO::endl;
        IO::cout << "|| (u*n)_inflow (u - u*) ||_L2(Gamma)                 =  " << interf_err_inflow                  << IO::endl;
        IO::cout << "|| (sigma*h+|u|+nu/h)^(+1/2) (u - u*)*n ||_L2(Gamma)  =  " << interf_err_mass_cons               << IO::endl;
        IO::cout << "---------------------------------------------------------"       << IO::endl;
        IO::cout << "-------------- Error on Functionals from solution  ------------"       << IO::endl;
        IO::cout << " | sin(x) ( u,x - u,x exact ) |                       = " << functional                          << IO::endl;
        IO::cout << "---------------------------------------------------------"       << IO::endl;
      }

      // append error of the last time step to the error file
      if ((step_==stepmax_) or (time_==maxtime_))// write results to file
      {
        std::ostringstream temp;
        const std::string simulation = DRT::Problem::Instance()->OutputControlFile()->FileName();
        const std::string fname = simulation+".xfem_abserror";

        std::ofstream f;
        f.open(fname.c_str(),std::fstream::ate | std::fstream::app);
        f << "#| " << simulation << "\n";
        f << "#| Step"
          << " | Time"
          << " | || u - u_h ||_L2(Omega)"
          << " | || grad( u - u_h ) ||_L2(Omega)"
          << " | || u - u_h ||_H1(Omega)"
          << " | || p - p_h ||_L2(Omega)"
          << " | || nu^(+1/2) grad( u - u_h ) ||_L2(Omega)"
          << " | || nu^(-1/2) (p - p_h) ||_L2(Omega)"
          << " | || sigma^(+1/2) ( u - u_h ) ||_L2(Omega)"
          << " | || Phi^(+1/2) (p - p_h) ||_L2(Omega)"
          << " | || nu^(+1/2) (u - u*) ||_H1/2(Gamma)"
          << " | || nu^(+1/2) grad( u - u_h )*n ||_H-1/2(Gamma)"
          << " | || nu^(-1/2) (p - p_h)*n ||_H-1/2(Gamma)"
          << " | || (u*n)_inflow (u - u*) ||_L2(Gamma)"
          << " | || (sigma*h+|u|+nu/h)^(+1/2) (u - u*)*n ||_L2(Gamma)"
          << " |  | sin(x) ( u,x - u,x exact ) | "
          << " |\n";
        f << step_ << " "
          << time_ << " "
          << dom_err_vel_L2 << " "
          << dom_err_vel_H1_semi << " "
          << dom_err_vel_H1 << " "
          << dom_err_pre_L2 << " "
          << dom_err_vel_H1_semi_nu_scaled << " "
          << dom_err_pre_L2_nu_scaled << " "
          << dom_err_vel_L2_sigma_scaled << " "
          << dom_err_pre_L2_Phi_scaled << " "
          << interf_err_Honehalf << " "
          << interf_err_Hmonehalf_u << " "
          << interf_err_Hmonehalf_p << " "
          << interf_err_inflow << " "
          << interf_err_mass_cons << " "
          << functional << " "
          <<"\n";
        f.flush();
        f.close();
      }

      std::ostringstream temp;
      const std::string simulation = DRT::Problem::Instance()->OutputControlFile()->FileName();
      const std::string fname = simulation+"_time.xfem_abserror";

      if(step_==1)
      {
        std::ofstream f;
        f.open(fname.c_str());

        f << "#| Step"
          << " | Time"
          << " | || u - u_h ||_L2(Omega)"
          << " | || grad( u - u_h ) ||_L2(Omega)"
          << " | || u - u_h ||_H1(Omega)"
          << " | || p - p_h ||_L2(Omega)"
          << " | || nu^(+1/2) grad( u - u_h ) ||_L2(Omega)"
          << " | || nu^(-1/2) (p - p_h) ||_L2(Omega)"
          << " | || sigma^(+1/2) ( u - u_h ) ||_L2(Omega)"
          << " | || Phi^(+1/2) (p - p_h) ||_L2(Omega)"
          << " | || nu^(+1/2) (u - u*) ||_H1/2(Gamma)"
          << " | || nu^(+1/2) grad( u - u_h )*n ||_H-1/2(Gamma)"
          << " | || nu^(-1/2) (p - p_h)*n ||_H-1/2(Gamma)"
          << " | || (u*n)_inflow (u - u*) ||_L2(Gamma)"
          << " | || (sigma*h+|u|+nu/h)^(+1/2) (u - u*)*n ||_L2(Gamma)"
          << " |  | sin(x) ( u,x - u,x exact ) | "
          << " |\n";
        f << step_ << " "
          << time_ << " "
          << dom_err_vel_L2 << " "
          << dom_err_vel_H1_semi << " "
          << dom_err_vel_H1 << " "
          << dom_err_pre_L2 << " "
          << dom_err_vel_H1_semi_nu_scaled << " "
          << dom_err_pre_L2_nu_scaled << " "
          << dom_err_vel_L2_sigma_scaled << " "
          << dom_err_pre_L2_Phi_scaled << " "
          << interf_err_Honehalf << " "
          << interf_err_Hmonehalf_u << " "
          << interf_err_Hmonehalf_p << " "
          << interf_err_inflow << " "
          << interf_err_mass_cons << " "
          << functional << " "
          <<"\n";

        f.flush();
        f.close();
      }
      else
      {
        std::ofstream f;
        f.open(fname.c_str(),std::fstream::ate | std::fstream::app);
        f << step_ << " "
          << time_ << " "
          << dom_err_vel_L2 << " "
          << dom_err_vel_H1_semi << " "
          << dom_err_vel_H1 << " "
          << dom_err_pre_L2 << " "
          << dom_err_vel_H1_semi_nu_scaled << " "
          << dom_err_pre_L2_nu_scaled << " "
          << dom_err_vel_L2_sigma_scaled << " "
          << dom_err_pre_L2_Phi_scaled << " "
          << interf_err_Honehalf << " "
          << interf_err_Hmonehalf_u << " "
          << interf_err_Hmonehalf_p << " "
          << interf_err_inflow << " "
          << interf_err_mass_cons << " "
          << functional << " "
          <<"\n";

        f.flush();
        f.close();
      }
    } // myrank = 0
  }

  return Teuchos::null;
}

void FLD::XFluid::ComputeErrorNorms(
  Teuchos::RCP<Epetra_SerialDenseVector> glob_dom_norms,
  Teuchos::RCP<Epetra_SerialDenseVector> glob_interf_norms,
  Teuchos::RCP<Epetra_SerialDenseVector> glob_stab_norms)
{
  // number of norms that have to be calculated
  const int num_dom_norms    = glob_dom_norms->Length();
  const int num_interf_norms = glob_interf_norms->Length();
  const int num_stab_norms   = glob_stab_norms->Length();

  Epetra_SerialDenseVector cpu_dom_norms(num_dom_norms);
  Epetra_SerialDenseVector cpu_interf_norms(num_interf_norms);
  Epetra_SerialDenseVector cpu_stab_norms(num_stab_norms);

  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("u and p at time n+1 (converged)", state_->velnp_);

  condition_manager_->SetState();

  // evaluate domain error norms and interface/boundary error norms at XFEM-interface
  // loop row elements
  const int numrowele = discret_->NumMyRowElements();
  for (int i=0; i<numrowele; ++i)
  {

    // local element-wise squared error norms
    Epetra_SerialDenseVector ele_dom_norms(num_dom_norms);
    Epetra_SerialDenseVector ele_interf_norms(num_interf_norms);


    // pointer to current element
    DRT::Element* actele = discret_->lRowElement(i);

    Teuchos::RCP<MAT::Material> mat = actele->Material();

    DRT::ELEMENTS::Fluid * ele = dynamic_cast<DRT::ELEMENTS::Fluid *>( actele );

    GEO::CUT::ElementHandle * e = state_->Wizard()()->GetElement( actele );

    DRT::Element::LocationArray la( 1 );

    DRT::ELEMENTS::FluidEleInterface * impl = DRT::ELEMENTS::FluidFactory::ProvideImplXFEM( actele->Shape(), "xfem");

    // xfem element
    if ( e!=NULL )
    {

      std::vector< GEO::CUT::plain_volumecell_set > cell_sets;
      std::vector< std::vector<int> > nds_sets;
      std::vector<std::vector< DRT::UTILS::GaussIntegration > >intpoints_sets;

      bool has_xfem_integration_rule =
          e->GetCellSets_DofSets_GaussPoints( cell_sets, nds_sets, intpoints_sets, false); //(include_inner=false)

      if(cell_sets.size() != nds_sets.size()) dserror("number of cell_sets and nds_sets not equal!");

      // loop over volume cells
      for( std::vector< GEO::CUT::plain_volumecell_set>::iterator s=cell_sets.begin();
          s!=cell_sets.end();
          s++)
      {
        GEO::CUT::plain_volumecell_set & cells = *s;
        const int set_counter = s - cell_sets.begin();
        const std::vector<int> & nds = nds_sets[set_counter];

        // get element location vector, dirichlet flags and ownerships
        actele->LocationVector(*discret_,nds,la,false);

        //------------------------------------------------------------
        // Evaluate interface integral errors
        // do cut interface condition

        // maps of sid and corresponding boundary cells ( for quadratic elements: collected via volumecells of subelements)
        std::map<int, std::vector<GEO::CUT::BoundaryCell*> > bcells;
        std::map<int, std::vector<DRT::UTILS::GaussIntegration> > bintpoints;

        for ( GEO::CUT::plain_volumecell_set::iterator i=cells.begin(); i!=cells.end(); ++i )
        {
          GEO::CUT::VolumeCell * vc = *i;
          if ( vc->Position()==GEO::CUT::Point::outside )
          {
              vc->GetBoundaryCells( bcells );
          }

          const int cellcount = std::distance(cells.begin(), i);

          if(!has_xfem_integration_rule) // use standard integration!!!
          {
            // get element location vector, dirichlet flags and ownerships
            actele->LocationVector(*discret_,la,false);

            Epetra_SerialDenseMatrix elemat1;
            Epetra_SerialDenseMatrix elemat2;
            Epetra_SerialDenseVector elevec2;
            Epetra_SerialDenseVector elevec3;
            params_->set<int>("action",FLD::calc_fluid_error);
            impl->EvaluateService(ele,
                *params_,
                mat,
                *discret_,
                la[0].lm_,
                elemat1,
                elemat2,
                ele_dom_norms,
                elevec2,
                elevec3);
          }
          else
          {
            if(cell_sets.size() != intpoints_sets.size()) dserror("number of cell_sets and intpoints_sets not equal!");

            //------------------------------------------------------------
            // Evaluate domain integral errors
            impl->ComputeError(ele,
                *params_,
                mat,
                *discret_,
                la[0].lm_,
                ele_dom_norms,
                intpoints_sets[set_counter][cellcount]
            );
          }
        }

        if ( bcells.size() > 0 )
        {
          // get boundary cell Gaussian points
          e->BoundaryCellGaussPointsLin( bcells, bintpoints);

          if(CouplingMethod() == INPAR::XFEM::Hybrid_LM_Cauchy_stress or
             CouplingMethod() == INPAR::XFEM::Hybrid_LM_viscous_stress or
             CouplingMethod() == INPAR::XFEM::Nitsche)
          {
            impl->ComputeErrorInterface(
                ele,
                *discret_,
                la[0].lm_,
                condition_manager_,
                mat,
                ele_interf_norms,
                bcells,
                bintpoints,
                cells,
                *params_
            );
          }
        } // bcells
      } // end of loop over volume-cell sets
    }
    // standard (no xfem) element
    else
    {
      // get element location vector, dirichlet flags and ownerships
      actele->LocationVector(*discret_,la,false);

      Epetra_SerialDenseMatrix elemat1;
      Epetra_SerialDenseMatrix elemat2;
      Epetra_SerialDenseVector elevec2;
      Epetra_SerialDenseVector elevec3;
      params_->set<int>("action",FLD::calc_fluid_error);
      impl->EvaluateService(ele,
          *params_,
          mat,
          *discret_,
          la[0].lm_,
          elemat1,
          elemat2,
          ele_dom_norms,
          elevec2,
          elevec3);
    }

    // sum up (on each processor)
    cpu_interf_norms += ele_interf_norms;

    // sum up (on each processor)
    cpu_dom_norms += ele_dom_norms;

  }//end loop over fluid elements

  //--------------------------------------------------------
  // reduce and sum over all procs
  for (int i=0; i<num_dom_norms; ++i) (*glob_dom_norms)(i) = 0.0;
  discret_->Comm().SumAll(cpu_dom_norms.Values(), glob_dom_norms->Values(), num_dom_norms);

  for (int i=0; i<num_interf_norms; ++i) (*glob_interf_norms)(i) = 0.0;
  discret_->Comm().SumAll(cpu_interf_norms.Values(), glob_interf_norms->Values(), num_interf_norms);


  //--------------------------------------------------------
  discret_->ClearState();
  condition_manager_->ClearState();

}


/*----------------------------------------------------------------------*
 |  check xfluid input parameters/ safety checks           schott 05/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::CheckXFluidParams() const
{
  // ----------------------------------------------------------------------
  // check XFLUID DYNAMIC/GENERAL parameter list
  // ----------------------------------------------------------------------

  return;
}


/*----------------------------------------------------------------------*
 |  Print fluid stabilization parameters                   schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::PrintStabilizationDetails() const
{
  FluidImplicitTimeInt::PrintStabilizationDetails();
  // output of interface stabilization details
  if (myrank_==0)
  {
    Teuchos::ParameterList *  interfstabparams          =&(params_->sublist("XFLUID DYNAMIC/STABILIZATION"));

    //---------------------------------------------------------------------------------------------

    IO::cout << "+------------------------------------------------------------------------------------+" << IO::endl;
    IO::cout << "                              INTERFACE-STABILIZATION                       \n" << IO::endl;
    IO::cout << "Stabilization type:      " << interfstabparams->get<std::string>("COUPLING_METHOD") << "\n";

    if(coupling_method_ == INPAR::XFEM::Hybrid_LM_Cauchy_stress or coupling_method_ == INPAR::XFEM::Hybrid_LM_viscous_stress)
      IO::cout << "HYBRID_LM_L2_PROJ:       " << interfstabparams->get<std::string>("HYBRID_LM_L2_PROJ") << "\n";

    if(coupling_method_ == INPAR::XFEM::Nitsche)
    {
      IO::cout << "NIT_STAB_FAC:                      " << interfstabparams->get<double>("NIT_STAB_FAC") << "\n";
      IO::cout << "VISC_STAB_TRACE_ESTIMATE:          " << interfstabparams->get<std::string>("VISC_STAB_TRACE_ESTIMATE") << "\n";
      IO::cout << "VISC_STAB_HK:                      " << interfstabparams->get<std::string>("VISC_STAB_HK")  << "\n";
    }

    if (coupling_method_ != INPAR::XFEM::Hybrid_LM_Cauchy_stress)
      IO::cout << "VISC_ADJOINT_SYMMETRY:             " << interfstabparams->get<std::string>("VISC_ADJOINT_SYMMETRY") << "\n";

    IO::cout << "MASS_CONSERVATION_COMBO:           " << interfstabparams->get<std::string>("MASS_CONSERVATION_COMBO") << "\n";
    IO::cout << "MASS_CONSERVATION_SCALING:         " << interfstabparams->get<std::string>("MASS_CONSERVATION_SCALING") << "\n";

    IO::cout << "GHOST_PENALTY_STAB:                " << interfstabparams->get<std::string>("GHOST_PENALTY_STAB") << "\n";
    IO::cout << "GHOST_PENALTY_TRANSIENT_STAB:      " << interfstabparams->get<std::string>("GHOST_PENALTY_TRANSIENT_STAB") << "\n";
    IO::cout << "GHOST_PENALTY_FAC:                 " << interfstabparams->get<double>("GHOST_PENALTY_FAC") << "\n";
    IO::cout << "GHOST_PENALTY_TRANSIENT_FAC:       " << interfstabparams->get<double>("GHOST_PENALTY_TRANSIENT_FAC") << "\n";
    IO::cout << "GHOST_PENALTY_2nd_STAB:            " << interfstabparams->get<std::string>("GHOST_PENALTY_2nd_STAB") << "\n";
    IO::cout << "GHOST_PENALTY_2nd_STAB_NORMAL:     " << interfstabparams->get<std::string>("GHOST_PENALTY_2nd_STAB_NORMAL") << "\n";


    IO::cout << "CONV_STAB_SCALING:                 " << interfstabparams->get<std::string>("CONV_STAB_SCALING") << "\n";

    IO::cout << "IS_PSEUDO_2D:                      " << interfstabparams->get<std::string>("IS_PSEUDO_2D") << "\n";
    IO::cout << "+------------------------------------------------------------------------------------+\n" << IO::endl;

  }

}


/*----------------------------------------------------------------------*
 |  Print information about current time step to screen    schott 02/15 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::PrintTimeStepInfo()
{

  // -------------------------------------------------------------------
  //                       output to screen
  // -------------------------------------------------------------------
  if (myrank_==0)
  {
    switch (timealgo_)
    {
    case INPAR::FLUID::timeint_stationary:
      printf("Stationary Fluid Solver - STEP = %4d/%4d \n",step_,stepmax_);
      break;
    case INPAR::FLUID::timeint_one_step_theta:
      printf("TIME: %11.4E/%11.4E  DT = %11.4E   One-Step-Theta  (theta = %11.2E)  STEP = %4d/%4d \n",
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
    {
      dserror("parameter out of range: IOP\n");
      break;
    }
    } /* end of switch(timealgo) */
  }
}


/*----------------------------------------------------------------------*
 |  Timeloop()                                             schott 02/15 |
 *----------------------------------------------------------------------*/
bool FLD::XFluid::NotFinished()
{
  // -------------------------------------------------------------------
  //                    stop criterium for timeloop
  // -------------------------------------------------------------------

  if(timealgo_ == INPAR::FLUID::timeint_stationary)
    return step_<stepmax_;
  else
    return step_<stepmax_ and time_<maxtime_;
}


/*----------------------------------------------------------------------*
 |  Timeloop()                                             schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::TimeLoop()
{
  if(myrank_ == 0) printf("START TIMELOOP (FLD::XFluid::TimeLoop) -- MAXTIME = %11.4E -- STEPMAX %4d\n\n",maxtime_,stepmax_);

  FluidImplicitTimeInt::TimeLoop();

  // print the results of time measurements
  Teuchos::TimeMonitor::summarize();
}



/*------------------------------------------------------------------------------------------------*
 | prepare a fluid time step                                                         schott 07/11 |
 *------------------------------------------------------------------------------------------------*/
void FLD::XFluid::PrepareTimeStep()
{

  if(myrank_ == 0) IO::cout << "PrepareTimeStep (FLD::XFluid::PrepareTimeStep) " << IO::endl;

  // -------------------------------------------------------------------
  //              reset counters used within timestep
  // -------------------------------------------------------------------
  // reset the state-class iterator for the new time step
  state_it_ = 0;
  itnum_out_= 0;


  // -------------------------------------------------------------------
  //              set time dependent parameters
  // -------------------------------------------------------------------
  IncrementTimeAndStep();

  condition_manager_->IncrementTimeAndStep(dta_);

  // -------------------------------------------------------------------
  // set time parameters dependent on time integration scheme and step
  // -------------------------------------------------------------------
  SetTheta();


  // -------------------------------------------------------------------
  //                     do explicit predictor step
  // -------------------------------------------------------------------
  DoPredictor();


  // -------------------------------------------------------------------
  //               set time parameter for element call
  // -------------------------------------------------------------------
  SetElementTimeParameter();

}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::XFluid::SetTheta()
{
  // Sets theta_ to a specific value for bdf2 and calculates
  // a pseudo-theta for genalpha (the latter in case of startalgo_)
  if (timealgo_ == INPAR::FLUID::timeint_stationary)
  {
    theta_   = 1.0;
    omtheta_ = 0.0;
  }
  else
  {
    // safety
    if (step_ < 1)
      dserror("number of time step is wrong");

    // do a backward Euler step for a user-defined number of starting steps
    if (step_<=numstasteps_)
    {
      if (myrank_==0)
      {
        std::cout << "Starting algorithm for OST active. "
            << "Performing step "<< step_ <<" of "<< numstasteps_
            << " Backward Euler starting steps"<<std::endl;
      }
      theta_   = 1.0;
      omtheta_ = 1.0-theta_;
    }
    else
    {
      // for OST
      if(timealgo_ == INPAR::FLUID::timeint_one_step_theta)
      {
        theta_ = params_->get<double>("theta");
        omtheta_ = 1.0-theta_;
      }

      // for BDF2, theta is set by the time-step sizes, 2/3 for const. dt
      if (timealgo_==INPAR::FLUID::timeint_bdf2)
      {
        theta_ = (dta_+dtp_)/(2.0*dta_ + dtp_);
        omtheta_ = 0.0;
      }
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::XFluid::DoPredictor()
{
  // no predictor in first time step
  if (step_>1)
  {
    if (predictor_ != "TangVel")
    {
      ExplicitPredictor();
    }
    else
    {
      PredictTangVelConsistAcc();
    }
  }
}


/*----------------------------------------------------------------------*
 |  prepare the nonlinear solver                           schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::PrepareXFEMSolve()
{
  // TODO: do we need to call PrepareXFEMSolve for each Newton increment when solving a monolithic system
  // can we shift this to PrepareTimeStep()?
  // -------------------------------------------------------------------
  // set new interface positions and possible values for XFEM Weak Dirichlet and Neumann BCs
  // -------------------------------------------------------------------
  condition_manager_->PrepareSolve();

  output_service_->GmshOutputDiscretization(eval_eos_, step_);
  // -------------------------------------------------------------------
  //  perform CUT, transform vectors from old dofset to new dofset and set state vectors
  // -------------------------------------------------------------------

  CutAndSetStateVectors();


  // -------------------------------------------------------------------
  //                 set old part of righthandside
  // -------------------------------------------------------------------
  SetOldPartOfRighthandside();


  // -------------------------------------------------------------------
  //         evaluate Dirichlet and Neumann boundary conditions
  // -------------------------------------------------------------------
  SetDirichletNeumannBC();

}


/*----------------------------------------------------------------------*
 |  solve the nonlinear problem                            schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::Solve()
{
  PrepareXFEMSolve();

  // ---------------------------------------------- nonlinear iteration
  // ------------------------------- stop nonlinear iteration when both
  //                                 increment-norms are below this bound
  const double  ittol        = params_->get<double>("tolerance for nonlin iter");

  //------------------------------ turn adaptive solver tolerance on/off
  const bool   isadapttol    = params_->get<bool>("ADAPTCONV");
  const double adaptolbetter = params_->get<double>("ADAPTCONV_BETTER",0.01);

  int  itnum = 0;
  bool stopnonliniter = false;

  dtsolve_  = 0.0;
  dtele_    = 0.0;
  dtfilter_ = 0.0;

  if(myrank_ == 0)
    printf("----------------------XFLUID-------  time step %2d ----------------------------------------\n", step_);

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
      // get cpu time
      const double tcpu=Teuchos::Time::wallTime();

      AssembleMatAndRHS( itnum );

      // end time measurement for element
      dtele_=Teuchos::Time::wallTime()-tcpu;
    }

    // blank residual DOFs which are on Dirichlet BC
    // We can do this because the values at the dirichlet positions
    // are not used anyway.
    // We could avoid this though, if velrowmap_ and prerowmap_ would
    // not include the dirichlet values as well. But it is expensive
    // to avoid that.

    output_service_->GmshResidualOutputDebug( "DEBUG_residual_wo_DBC", step_, itnum, state_ );


    // apply Dirichlet conditions to the residual vector by setting zeros into the residual
    state_->DBCMapExtractor()->InsertCondVector(
        state_->DBCMapExtractor()->ExtractCondVector(state_->Zeros()),
        state_->Residual());

    output_service_->GmshResidualOutputDebug( "DEBUG_residual", step_, itnum, state_ );

    if (updateprojection_)
    {
      // even if not ALE, we always need to update projection vectors due to changed cuts
      UpdateKrylovSpaceProjection();
    }

    // remove contributions of pressure mode
    // that would not vanish due to the projection
    if (projector_ != Teuchos::null)
      projector_->ApplyPT(*state_->Residual());

    if (ConvergenceCheck(itnum,itemax_,ittol))
      break;

    //--------- Apply Dirichlet boundary conditions to system of equations
    //          residual displacements are supposed to be zero at
    //          boundary conditions
    state_->IncVel()->PutScalar(0.0);
    LINALG::ApplyDirichlettoSystem(
        state_->SystemMatrix(),state_->IncVel(),state_->Residual(),
        state_->Zeros(),*(state_->DBCMapExtractor()->CondMap()));


//#if 1
//    const double cond_number = LINALG::Condest(static_cast<LINALG::SparseMatrix&>(*state_->sysmat_),Ifpack_Cheap, 1000);
//    // computation of significant digits might be completely bogus, so don't take it serious
//    const double tmp = std::abs(std::log10(cond_number*1.11022e-16));
//    const int sign_digits = (int)floor(tmp);
//    if (!myrank_)
//      IO::cout << " cond est: " << std::scientific << cond_number << ", max.sign.digits: " << sign_digits;
//#endif


    //-------solve for residual displacements to correct incremental displacements
    {
      // get cpu time
      const double tcpusolve=Teuchos::Time::wallTime();

      // do adaptive linear solver tolerance (not in first solve)
      if (isadapttol && itnum>1)
      {
        double currresidual = std::max(vresnorm_,presnorm_);
        currresidual = std::max(currresidual,incvelnorm_L2_/velnorm_L2_);
        currresidual = std::max(currresidual,incprenorm_L2_/prenorm_L2_);
        solver_->AdaptTolerance(ittol,currresidual,adaptolbetter);
      }

#if 0
      // print matrix in matlab format

            // matrix printing options (DEBUGGING!)
      IO::cout << "print matrix in matlab format to sparsematrix.mtl";

            Teuchos::RCP<LINALG::SparseMatrix> A = state_->SystemMatrix();
            if (A != Teuchos::null)
            {
              // print to file in matlab format
              const std::string fname = "sparsematrix.mtl";
              LINALG::PrintMatrixInMatlabFormat(fname,*(A->EpetraMatrix()));
              // print to screen
//              (A->EpetraMatrix())->Print(cout);
              // print sparsity pattern to file
              LINALG::PrintSparsityToPostscript( *(A->EpetraMatrix()) );
            }
            else
            {
              Teuchos::RCP<LINALG::BlockSparseMatrixBase> A = state_->BlockSystemMatrix();
              const std::string fname = "sparsematrix.mtl";
              LINALG::PrintBlockMatrixInMatlabFormat(fname,*(A));
            }

            IO::cout << " ...done" << IO::endl;
            // ScaleLinearSystem();  // still experimental (gjb 04/10)
#endif

     // scale system prior to solver call
     if (fluid_infnormscaling_!= Teuchos::null)
       fluid_infnormscaling_->ScaleSystem(state_->SystemMatrix(), *(state_->Residual()));

      // if Krylov space projection is used, check whether constant pressure
      // is in nullspace of sysmat_
     CheckMatrixNullspace();

     solver_->Solve(state_->SystemMatrix()->EpetraOperator(),state_->IncVel(),state_->Residual(),true,itnum==1, projector_);

     //TODO: here needed because of apply Dirichlet with explicit Dirichlet flag!? CHECK THIS
     solver_->Reset();


      // unscale solution
      if (fluid_infnormscaling_!= Teuchos::null)
        fluid_infnormscaling_->UnscaleSolution(state_->SystemMatrix(), *(state_->IncVel()),*(state_->Residual()));

      solver_->ResetTolerance();

      // end time measurement for solver
      dtsolve_ = Teuchos::Time::wallTime()-tcpusolve;
    }

    output_service_->GmshIncrementOutputDebug( "DEBUG_icnr", step_, itnum, state_ );

    // -------------------------------------------------------------------
    // update velocity and pressure values by increments
    // -------------------------------------------------------------------
    UpdateByIncrement();

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

    std::cout << "MAXNUMENTRIES: " << state_->sysmat_->EpetraMatrix()->MaxNumEntries() << std::endl;

  }

  // Reset the solver and so release the system matrix' pointer (enables to delete the state_->systemmatrix)
  solver_->Reset();


#if(0)
  if(gmsh_debug_out_)
  {
    std::string filename_base = "DEBUG_SOL";

    //--------------------------------------------------------------------
    Teuchos::RCP<const Epetra_Vector> output_col_velnp = DRT::UTILS::GetColVersionOfRowVector(discret_, state_->velnp_ );
    const std::string prefix("SOL");
    output_service_->GmshOutput(filename_base, prefix, step_, state_it_, state_->Wizard(), output_col_velnp );

    //--------------------------------------------------------------------
    condition_manager_->GmshOutput(filename_base, step_, gmsh_step_diff_, gmsh_debug_out_screen_);
  }
#endif

}

bool FLD::XFluid::ConvergenceCheck(
  int          itnum,
  int          itemax,
  const double ittol)
{
  bool stopnonliniter = false;

  incvelnorm_L2_ = 0.0;
  incprenorm_L2_ = 0.0;

  velnorm_L2_ = 0.0;
  prenorm_L2_ = 0.0;

  vresnorm_ = 0.0;
  presnorm_ = 0.0;

  Teuchos::RCP<Epetra_Vector> onlyvel = state_->VelPresSplitter()->ExtractOtherVector(
    state_->Residual());
  onlyvel->Norm2(&vresnorm_);

  state_->VelPresSplitter()->ExtractOtherVector(state_->IncVel(),onlyvel);
  onlyvel->Norm2(&incvelnorm_L2_);

  state_->VelPresSplitter()->ExtractOtherVector(state_->Velnp(),onlyvel);
  onlyvel->Norm2(&velnorm_L2_);

  Teuchos::RCP<Epetra_Vector> onlypre = state_->VelPresSplitter()->ExtractCondVector(
    state_->Residual());
  onlypre->Norm2(&presnorm_);

  state_->VelPresSplitter()->ExtractCondVector(state_->IncVel(),onlypre);
  onlypre->Norm2(&incprenorm_L2_);

  state_->VelPresSplitter()->ExtractCondVector(state_->Velnp(),onlypre);
  onlypre->Norm2(&prenorm_L2_);

  // care for the case that nothing really happens in the velocity
  // or pressure field
  if (velnorm_L2_ < 1e-5) velnorm_L2_ = 1.0;
  if (prenorm_L2_ < 1e-5) prenorm_L2_ = 1.0;

  //-------------------------------------------------- output to screen
  /* special case of very first iteration step:
      - solution increment is not yet available
      - convergence check is not required (we solve at least once!)    */
  if (itnum == 1)
  {
    if (myrank_ == 0)
    {
      printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   |      --      |      --      |",
             itnum,itemax,ittol,vresnorm_,presnorm_);
      printf(" (      --     ,te=%10.3E",dtele_);
      if (turbmodel_==INPAR::FLUID::dynamic_smagorinsky)
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
    if (vresnorm_ <= ittol and presnorm_ <= ittol and
        incvelnorm_L2_/velnorm_L2_ <= ittol and incprenorm_L2_/prenorm_L2_ <= ittol)
    {
      stopnonliniter = true;
      if (myrank_ == 0)
      {
        printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   | %10.3E   | %10.3E   |",
               itnum,itemax,ittol,vresnorm_,presnorm_,
               incvelnorm_L2_/velnorm_L2_,incprenorm_L2_/prenorm_L2_);
        printf(" (ts=%10.3E,te=%10.3E",dtsolve_,dtele_);
        if (turbmodel_==INPAR::FLUID::dynamic_smagorinsky)
        {
          printf(",tf=%10.3E",dtfilter_);
        }
        printf(")\n");
        printf("+------------+-------------------+--------------+--------------+--------------+--------------+\n");

        FILE* errfile = params_->get<FILE*>("err file",NULL);
        if (errfile!=NULL)
        {
          fprintf(errfile,"fluid solve:   %3d/%3d  tol=%10.3E[L_2 ]  vres=%10.3E  pres=%10.3E  vinc=%10.3E  pinc=%10.3E\n",
                  itnum,itemax,ittol,vresnorm_,presnorm_,
                  incvelnorm_L2_/velnorm_L2_,incprenorm_L2_/prenorm_L2_);
        }
      }
    }
    else // if not yet converged
      if (myrank_ == 0)
      {
        printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   | %10.3E   | %10.3E   |",
               itnum,itemax,ittol,vresnorm_,presnorm_,
               incvelnorm_L2_/velnorm_L2_,incprenorm_L2_/prenorm_L2_);
        printf(" (ts=%10.3E,te=%10.3E",dtsolve_,dtele_);
        if (turbmodel_==INPAR::FLUID::dynamic_smagorinsky)
        {
          printf(",tf=%10.3E",dtfilter_);
        }
        printf(")\n");
      }
  }

  // warn if itemax is reached without convergence, but proceed to
  // next timestep...
  if ((itnum == itemax) and (vresnorm_ > ittol or presnorm_ > ittol or
                             incvelnorm_L2_ /velnorm_L2_ > ittol or
                             incprenorm_L2_/prenorm_L2_ > ittol))
  {
    stopnonliniter=true;
    if (myrank_ == 0)
    {
      printf("+---------------------------------------------------------------+\n");
      printf("|            >>>>>> not converged in itemax steps!              |\n");
      printf("+---------------------------------------------------------------+\n");

      FILE* errfile = params_->get<FILE*>("err file",NULL);
      if (errfile!=NULL)
      {
        fprintf(errfile,"fluid unconverged solve:   %3d/%3d  tol=%10.3E[L_2 ]  vres=%10.3E  pres=%10.3E  vinc=%10.3E  pinc=%10.3E\n",
                itnum,itemax,ittol,vresnorm_,presnorm_,
                incvelnorm_L2_/velnorm_L2_,incprenorm_L2_/prenorm_L2_);
      }
    }
  }

  return stopnonliniter;
}

void FLD::XFluid::LinearSolve()
{
  dserror("LinearSolve not implemented for Xfluid");
}


void FLD::XFluid::InitKrylovSpaceProjection()
{
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
    if (myrank_ == 0)
      std::cout << "\nSetup of KrylovSpaceProjection in fluid field\n" << std::endl;
  }
  else if (numfluid == 0)
  {
    updateprojection_ = false;
    projector_ = Teuchos::null;
  }
  else
    dserror("Received more than one KrylovSpaceCondition for fluid field");
  return;
}



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*--------------------------------------------------------------------------*
 | setup Krylov projector including first fill                    nis Feb13 |
 *--------------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::XFluid::SetupKrylovSpaceProjection(DRT::Condition* kspcond)
{
  /*
   * Krylov space projection in the XFEM
   * - generally, the Krylov space projection is possible, if there are no perturbations introduced by
   *   inaccurate integration
   * - the kernel vector c (0,0,0,1; 0,0,0,1; ....), however filled in this way for all dofsets in case of multiple dofsets
   * - if the projection fails, then there is is maybe an inconsistency between the volume and surface integration
   *   on cut elements (either you choose a smaller VOLUME-tolerance in cut_tolerance or choose DirectDivergence instead
   *   of the Tesselation subtetrahedralization, then the surface will be triangulated independent of the integration cells
   * - otherwise there could be further geometric! inconsistencies in the transformation in case of warped volume elements
   */

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

  // set flag for projection update true only if ALE and integral weights
  if (alefluid_ and (*weighttype=="integration"))
    updateprojection_ = true;

  projector_ = Teuchos::rcp(new LINALG::KrylovProjector(activemodeids,weighttype,discret_->DofRowMap()));

  // update the projector
  UpdateKrylovSpaceProjection();

} // XFluid::SetupKrylovSpaceProjection

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*--------------------------------------------------------------------------*
 | update projection vectors w_ and c_ for Krylov projection      nis Feb13 |
 *--------------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::XFluid::UpdateKrylovSpaceProjection()
{
  // get Teuchos::RCP to kernel vector of projector
  Teuchos::RCP<Epetra_MultiVector> c = projector_->GetNonConstKernel();
  Teuchos::RCP<Epetra_Vector> c0 = Teuchos::rcp((*c)(0),false);
  c0->PutScalar(0.0);

  // extract vector of pressure-dofs
  Teuchos::RCP<Epetra_Vector> presmode = state_->velpressplitter_->ExtractCondVector(*c0);

  const std::string* weighttype = projector_->WeightType();

  // compute w_ as defined in dat-file
  if(*weighttype == "pointvalues")
  {
    // Smart xfluid people put dserror here. I guess they had there reasons. KN
    dserror("Pointvalues for weights is not supported for xfluid, choose integration in dat-file");

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
    mode_params.set<int>("action",FLD::integrate_shape);

    if (alefluid_)
    {
      discret_->SetState("dispnp",state_->dispnp_);
    }

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
    IntegrateShapeFunction(mode_params, *discret_, w0);

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

} // XFluid::UpdateKrylovSpaceProjection

/*--------------------------------------------------------------------------*
 | check if constant pressure mode is in kernel of sysmat_     nissen Jan13 |
 *--------------------------------------------------------------------------*/
void FLD::XFluid::CheckMatrixNullspace()
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

    state_->sysmat_->Apply(*c,result);

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



/*--------------------------------------------------------------------------*
 | update the veln-vector with the stepinc to obtain new iteration velnp,   |
 | cut and set new state-vectors, perform time-integration, apply bcs       |
 | evaluate the fluid at the new interface position            schott 08/14 |
 *--------------------------------------------------------------------------*/
void FLD::XFluid::Evaluate(
  Teuchos::RCP<const Epetra_Vector> stepinc ///< solution increment between time step n and n+1, stepinc has to match the current xfluid dofmaps
  )
{
  //--------------------------------------------------------------------------------------------
  // FIRST: update the current velnp vector with the increment from the monolithic solve
  //--------------------------------------------------------------------------------------------

  if (stepinc!=Teuchos::null) // non-first call, when a step increment is already available (also when restarting the global monolithic Newto)
  {
    //-----------------------------
    // update the velnp vector such that the new iteration is stored in velnp
    //-----------------------------
    // set the new solution we just got. Note: the solution we got here
    // is the time step increment which means the sum of all iteration
    // increments of the time step.

    // Take Dirichlet values from last velnp and add stepinc to veln for non-Dirichlet values.
    // * the stepinc should contain the Dirichlet values, however, when using an iterative solver the Dirichlet values
    //   of Newton increment might just be approximately zero. In order to strictly set the Dirichlet values to zero
    //   we set them here again.
    // * for each call of PrepareXFEMSolve (see below) the velnp-vector obtains accurate Dirichlet values
    // * therefore we directly can copy the Dirichlet values from the last iteration
    // * further, in the next PrepareXFEMSolve()-call, after performing time-integration,
    //   the DBCs are set again in velnp

    Teuchos::RCP<Epetra_Vector> velnp_tmp = LINALG::CreateVector(*discret_->DofRowMap(),true);

    state_->incvel_->Update(1.0, *stepinc, -1.0, *state_->velnp_, 0.0);
    state_->incvel_->Update(1.0, *state_->veln_, 1.0);

    // update the current u^(n+1,i+1) = u^n + (u^(n+1,i+1)-u^n) = veln_ + stepinc
    velnp_tmp->Update(1.0, *state_->veln_, 1.0, *stepinc, 0.0);

    // take the Dirichlet values from velnp and insert them in velnp_tmp
    state_->dbcmaps_->InsertCondVector(state_->dbcmaps_->ExtractCondVector(state_->velnp_), velnp_tmp );

    // set the whole vector with u^(n+1,i+1) including the Dirichlet values to velnp_
    state_->velnp_->Update(1.0, *velnp_tmp, 0.0);
  }
  else // the first call in a new time-step
  {
    // for the first call in a new time-step the initialization of velnp_ is not allowed as
    // velnp_ includes a predicted solution (set in PrepareTimeStep).
    // This predicted solution does not include the DBCs yet, however, in the following PrepareXFEMSolve()-call
    // veln_ and the predicted solution velnp_ will be mapped to the new interface position and afterwards
    // DBCs will be set in velnp_.
  }


  //--------------------------------------------------------------------------------------------
  // SECOND:
  // - cut at the new interface position
  // - create new state vectors
  // - did the dofsets change between last Newton iteration and current Newton iteration?
  // - perform time-integration between t^n and t^(n+1) at current interface position (which updates veln) and
  // - transform current iteration velnp_ip to new interface position by a simple copy when dofsets did not change
  //   or via a pseudo-time-integration as a kind of predictor in case that restart of the Newton is necessary
  //   (this includes an update of the permutation map necessary in the monolithic approach for updating the stepinc)
  // TODO: - apply a fluid predictor based on the new interface position
  // - set history values
  // - apply Dirichlet and Neumann boundary conditions
  //--------------------------------------------------------------------------------------------

  output_service_->GmshIncrementOutputDebug( "DEBUG_icnr", step_, itnum_out_, state_ );

  // TODO:maybe we can choose a more intelligent update such that we can reuse graphs of the matrix during the monolithic
  // xfsi solve...
  // currently we use fixed itnum = 1, it is okay as a new graph of the systemmatrix is created in the state-class evaluate routine
  int itnum = 1;
  itnum_out_++;


  PrepareXFEMSolve();


  //--------------------------------------------------------------------------------------------
  // THIRD: evaluate systemmatrix and rhs
  //--------------------------------------------------------------------------------------------


  // -------------------------------------------------------------------
  // call elements to calculate system matrix and RHS
  // -------------------------------------------------------------------
  {
    // get cpu time
    const double tcpu=Teuchos::Time::wallTime();

    AssembleMatAndRHS( itnum );

    // end time measurement for element
    dtele_=Teuchos::Time::wallTime()-tcpu;
  }


  // -------------------------------------------------------------------
  // write gmsh debug output for fluid residual directly after the fluid is evaluated
  // -------------------------------------------------------------------
  output_service_->GmshResidualOutputDebug( "DEBUG_residual_wo_DBC", step_, itnum_out_, state_ );
  output_service_->GmshSolutionOutputDebug( "DEBUG_sol", step_, itnum_out_, state_ );

  return;
}



/*----------------------------------------------------------------------*
 |  time update                                            schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::TimeUpdate()
{
  if(timealgo_ == INPAR::FLUID::timeint_stationary) return;


  if(myrank_ == 0) IO::cout << "FLD::XFluid::TimeUpdate " << IO::endl;

  Teuchos::ParameterList *  stabparams=&(params_->sublist("RESIDUAL-BASED STABILIZATION"));

  if(stabparams->get<std::string>("TDS") == "time_dependent")
  { dserror("check this implementation");
    const double tcpu=Teuchos::Time::wallTime();

    if(myrank_==0)
    {
      std::cout << "time update for subscales";
    }

    // call elements to calculate system matrix and rhs and assemble
    // this is required for the time update of the subgrid scales and
    // makes sure that the current subgrid scales correspond to the
    // current residual
    AssembleMatAndRHS();

    // create the parameters for the discretization
    Teuchos::ParameterList eleparams;

    // update time paramters
    SetGamma(eleparams);


    eleparams.set("dt"     ,dta_    );

    // call loop over elements to update subgrid scales
    discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);

    if(myrank_==0)
    {
      std::cout << "("<<Teuchos::Time::wallTime()-tcpu<<")\n";
    }
  }

  // Compute accelerations
  {
    Teuchos::RCP<Epetra_Vector> onlyaccn  = state_->velpressplitter_->ExtractOtherVector(state_->accn_ );
    Teuchos::RCP<Epetra_Vector> onlyaccnp = state_->velpressplitter_->ExtractOtherVector(state_->accnp_);
    Teuchos::RCP<Epetra_Vector> onlyvelnm = state_->velpressplitter_->ExtractOtherVector(state_->velnm_);
    Teuchos::RCP<Epetra_Vector> onlyveln  = state_->velpressplitter_->ExtractOtherVector(state_->veln_ );
    Teuchos::RCP<Epetra_Vector> onlyvelnp = state_->velpressplitter_->ExtractOtherVector(state_->velnp_);

    CalculateAcceleration(
        onlyvelnp,
        onlyveln ,
        onlyvelnm,
        onlyaccn ,
        onlyaccnp);

    // copy back into global vector
    LINALG::Export(*onlyaccnp,*state_->accnp_);
  }


  // update old acceleration
  state_->accn_->Update(1.0,*state_->accnp_,0.0);

  // velocities/pressures of this step become most recent
  // velocities/pressures of the last step
  state_->velnm_->Update(1.0,*state_->veln_ ,0.0);
  state_->veln_ ->Update(1.0,*state_->velnp_,0.0);

  if (alefluid_)
  {
    // displacements of this step becomes most recent
    // displacements of the last step
    dispnm_->Update(1.0,*dispn_,0.0);
    dispn_->Update(1.0,*dispnp_,0.0);

    // gridvelocities of this step become most recent
    // gridvelocities of the last step
    gridvn_->Update(1.0,*gridvnp_,0.0);
  }

  // update of interface fields (interface velocity and interface displacements)
  condition_manager_->UpdateStateVectors();

} //XFluid::TimeUpdate()


/*----------------------------------------------------------------------*
 |  cut at interface positions, transform vectors, perform              |
 | time integration and set new vectors                    schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::CutAndSetStateVectors()
{
  const bool screen_out = false;

  //------------------------------------------------------------------------------------
  // not required for stationary time integration
  if( timealgo_ == INPAR::FLUID::timeint_stationary ) return;

  //------------------------------------------------------------------------------------
  // not required if neither the background mesh nor the interfaces move

  // get info from condition_manager_ if at least one coupling object has moving interfaces
  const bool has_moving_interface = condition_manager_->HasMovingInterface();
  const bool moving_meshes = ( has_moving_interface or alefluid_ );

  if( !moving_meshes ) return;
  //------------------------------------------------------------------------------------


  if(myrank_==0)
  {
    // counter will be increased when the new state class is created

    IO::cout << "======================================================\n";
    IO::cout << "CutAndSetStateVectors: state-class iterator: " <<  state_it_ +1 << "\n";
    IO::cout << "======================================================\n";
  }


  if(step_ <= 0) return; // do not perform XFEM-time-integration for step 0

  //------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------
  //                             XFEM TIME-INTEGRATION
  //------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------


  // TODO: ADAPT for partitioned fsi

  bool firstcall_in_timestep = false;

  if(state_it_ == 0) firstcall_in_timestep=true;

  //----------------------------------------------------------------
  //---------------- STORE OLD STATE DATA --------------------------
  //----------------------------------------------------------------

  // save state data from the last time-step before the first iteration in a new time step is done and
  // save state data from the last (Newton, partitioned) iteration-step
  XTimint_StoreOldStateData(firstcall_in_timestep);

  //----------------------------------------------------------------
  //------------  NEW STATE CLASS including CUT  -------------------
  //----------------------------------------------------------------

  // create new state class object
  // state_it_ has been increased by one now
  // performs cut at current interface position and creates new vectors and a new system-matrix
  CreateState();


  //----------------------------------------------------------------
  //-------- TRANSFER veln_Int_n -> veln_Int_n+1_i+1  --------------
  //----------------------------------------------------------------

  // Transfer vectors from old time-step t^n w.r.t dofset and interface position from t^n
  // to vectors w.r.t current dofset and interface position
  XTimint_DoTimeStepTransfer(screen_out);


  //----------------------------------------------------------------
  //-------- TRANSFER velnp_Int_n+1_i -> velnp_Int_n+1_i+1  --------
  //----------------------------------------------------------------

  // Transfer vectors within the same time-step t^n+1 w.r.t dofset and interface position from last iteration
  // to vectors w.r.t current dofset and interface position
  //
  // NOTE:
  // fluid predictor has been called in PrepareTimeStep, therefore veln_ != velnp_, so we have to map both vectors,
  // also in the first call of a new time-step.
  // When SL is necessary to map velnp_, it might worsen the quality of the predicted solution:
  // * for partitioned FSI:
  //   it is possible to start the Fluid-Newton from veln_ (use a steady-state predictor afterwards),
  //   this usually yields more iterations however it does not influence the Convergence-behaviour of the staggered scheme
  // * for monolithic FSI:
  //   remark that in case that SL has to be used for mapping velnp_ it is NOT reasonable to restart the Newton from veln_
  //   since then we loose the whole information of the fluid-increments and convergence is not guaranteed at all!
  //TODO: what to do then?

  bool increment_tranfer_success = XTimint_DoIncrementStepTransfer(screen_out, firstcall_in_timestep);


#if(1) // just possible for partitioned FSI, the usage for pure fluids overwrites the fluid-predictor
  //------------------------------------------------------------------------------------
  //      set initial start vectors for new time step (steady-state predictor)
  //------------------------------------------------------------------------------------

  if(!increment_tranfer_success)
  {
    //velocity as start value for first Newton step
    state_->velnp_->Update(1.0,*state_->veln_,0.0);  // use old velocity as start value
    state_->accnp_->Update(1.0,*state_->accn_,0.0);  // use old velocity as start value
  }
#endif

  //---------------------------------- GMSH SOLUTION OUTPUT (reference/predicted solution fields for pressure, velocity, acc) ------------------------

  // write gmsh-output for reference solution fields
  // reference solution output

  //-------------
  // output for the reference solution veln
  output_service_->GmshSolutionOutputPrevious( "TIMINT_N_", step_, state_, state_it_ );

  //-------------
  // output for the predicted iteration velnp
  output_service_->GmshSolutionOutput( "TIMINT_NP_", step_, state_, state_it_ );


  if(myrank_==0 and screen_out) std::cout << "finished CutAndSetStateVectors()" << std::endl;


  return;
}



/*----------------------------------------------------------------------*
 |  store state data from old time-step t^n                schott 04/14 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::XTimint_StoreOldStateData(const bool firstcall_in_timestep)
{

  if(firstcall_in_timestep)
  {
    // store the solution of the old time step t^n w.r.t the old interface position
    veln_Intn_ = Teuchos::rcp(new Epetra_Vector(*discret_->DofRowMap()));
    *veln_Intn_  = *(state_->veln_);
    accn_Intn_ = Teuchos::rcp(new Epetra_Vector(*discret_->DofRowMap()));
    *accn_Intn_  = *(state_->accn_);

    // for BDF2
    velnm_Intn_ = Teuchos::rcp(new Epetra_Vector(*discret_->DofRowMap()));
    *velnm_Intn_  = *(state_->velnm_);

    // safe the old wizard and dofset w.r.t the interface position of the last time-step
    wizard_Intn_ = state_->Wizard();
    dofset_Intn_ = state_->DofSet();

    // safe the old dofmap
    dofcolmap_Intn_ = Teuchos::rcp(new Epetra_Map(*discret_->DofColMap()));
  }

  //------------------------------------------
  // store the last velocity solution w.r.t the last interface position (last XFSI iteration or last time-step solution for first-call)
  // to get mapped as fluid predictor for next XFSI iteration
  velnp_Intnpi_ = Teuchos::rcp(new Epetra_Vector(*discret_->DofRowMap()));
  *velnp_Intnpi_ = *state_->velnp_;

  // get the wizard w.r.t the last interface position (last XFSI iteration)
  wizard_Intnpi_ = state_->Wizard();
  dofset_Intnpi_ = state_->DofSet();

  return;
}



/*----------------------------------------------------------------------*
 |  is a restart of the global monolithic system necessary?             |
 |                                                         schott 08/14 |
 *----------------------------------------------------------------------*/
bool FLD::XFluid::XTimint_CheckForMonolithicNewtonRestart(
    const bool                        timint_ghost_penalty,    ///< dofs have to be reconstructed via ghost penalty reconstruction techniques
    const bool                        timint_semi_lagrangean,  ///< dofs have to be reconstructed via semi-Lagrangean reconstruction techniques
    Teuchos::RCP<DRT::Discretization> dis,                     ///< discretization
    Teuchos::RCP<XFEM::XFEMDofSet>    dofset_i,                ///< dofset last iteration
    Teuchos::RCP<XFEM::XFEMDofSet>    dofset_ip,               ///< dofset current iteration
    const bool                        screen_out               ///< screen output?
)
{
  discret_->Comm().Barrier();
  TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XTimint_CheckForMonolithicNewtonRestart" );

  // is a Newton restart necessary? initialize
  bool restart_necessary = false;


  // Restart the global monolithic system in the case that for at least one node the number of dofsets has changed
  // or for at least one node Semi-Lagrangean (SL) or Ghost-Penalty (GP) techniques have to be used
  // to transfer data between the current and last Newton iteration
  // Remark
  // * that pure copying is also possible when the global system changes (e.g. copy 1 ghost set -to-> 2 ghost sets)
  // * that SL or GP usually changes the increment/residual very much, such that the convergence seems to
  //   stagnate or diverge. Therefore we perform a restart to indicate the larger manipulation of the system

  //---------------
  // check if the dofsets changed
  const bool dofsets_changed = XTimint_ChangedDofsets( dis, dofset_i, dofset_ip);

  if(myrank_ == 0 and screen_out)
  {
    if(dofsets_changed) IO::cout << " CHANGING DOFSETS in the last two iterations " << IO::endl;
    else                IO::cout << " NON-CHANGING DOFSETS in the last two iterations " << IO::endl;
  }

  //---------------
  // restart of global monolithic Newton necessary?
  const bool pure_copying_possible = (!timint_ghost_penalty and !timint_semi_lagrangean);

  if( !pure_copying_possible or dofsets_changed )
  {
    restart_necessary = true;
  }
  else
  {
    restart_necessary = false;
  }

  if(myrank_ == 0 and screen_out)
  {
    if(restart_necessary) IO::cout << " RESTART of NEWTON necessary if not the first run after restarting/starting a timestep "     << IO::endl;
    else                  IO::cout << " RESTART of NEWTON not necessary " << IO::endl;
  }

  return restart_necessary;
}



/*----------------------------------------------------------------------*
 |  did the dofsets change?                                schott 08/14 |
 *----------------------------------------------------------------------*/
bool FLD::XFluid::XTimint_ChangedDofsets(
    Teuchos::RCP<DRT::Discretization> dis,                       ///< discretization
    Teuchos::RCP<XFEM::XFEMDofSet>    dofset,                    ///< first dofset
    Teuchos::RCP<XFEM::XFEMDofSet>    dofset_other               ///< other dofset
)
{
  //---------------
  // changed dofsets on this proc?
  // Use overloaded == operator for XFEM::XFEMDofset, comparison based on number of dofsets per node
  int changed_dofsets_proc_count = (int)(*dofset != *dofset_other);

  // assume changed dofsets
  int changed_dofsets_glob_max = 0;

  // check if at least one proc has changed dofsets? (maximum or sum of counts > 0)
  dis->Comm().MaxAll(&changed_dofsets_proc_count, &changed_dofsets_glob_max, 1);
  const bool changed_dofsets_glob = (changed_dofsets_glob_max > 0);

  return changed_dofsets_glob;
}



/*----------------------------------------------------------------------*
 | Transfer vectors from old time-step t^n w.r.t dofset and             |
 | interface position from t^n to vectors w.r.t current dofset and      |
 | interface position                                      schott 08/14 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::XTimint_DoTimeStepTransfer(const bool screen_out)
{
  //---------------------------------------------------------------
  if(myrank_==0 and screen_out) IO::cout << "XFEM::TIMEINTEGRATION: ..." << IO::endl;

  //---------------------------------------------------------------
  if(timealgo_ !=  INPAR::FLUID::timeint_one_step_theta) dserror("check which vectors have to be reconstructed for non-OST scheme");

  //---------------------------------------------------------------
  const Epetra_Map* newdofrowmap = discret_->DofRowMap();

  // all vectors that have to be transferred from old dofset at t^n to new dofset at t^(n+1=
  std::vector<Teuchos::RCP<const Epetra_Vector> > oldRowStateVectors;
  std::vector<Teuchos::RCP<Epetra_Vector> >       newRowStateVectors;

  // reconstruction map for nodes and its dofsets - how do we have to reconstruct the sinlge dofs
  std::map<int, std::vector<INPAR::XFEM::XFluidTimeInt> >             node_to_reconstr_method;
  std::map<INPAR::XFEM::XFluidTimeInt, std::map<int,std::set<int> > > reconstr_method_to_node;
  // vector of DOF-IDs which are Dirichlet BCs for ghost penalty reconstruction method
  Teuchos::RCP<std::set<int> > dbcgids = Teuchos::rcp(new std::set<int>());

  //------------------------------------------------------------------------------------
  // set interface state vectors for mesh coupling objects
  //------------------------------------------------------------------------------------
  condition_manager_->SetStateDisplacement();  // set idispnp, idispn and idispnpi vectors

  //------------------------------------------------------------------------------------
  // STEP 1: CopyDofsToNewMap and determine RECONSTRUCTION METHOD for missing values
  //------------------------------------------------------------------------------------
  //
  // REMARK:
  // * do this for row nodes only
  // * the cut information around the node should be available, since the cut is performed for col elements
  // * after transferring data from old interface position to new interface position the col vectors have to get
  //   exported from row vectors
  //------------------------------------------------------------------------------------

  //-----------------------------time integration----------------------

  // create time integration class just locally not to keep pointers to dofset and wizard...
  Teuchos::RCP<XFEM::XFluidTimeInt> xfluid_timeint = Teuchos::rcp(new
      XFEM::XFluidTimeInt(
          false, // is_newton_increment_transfer?
          discret_,
          condition_manager_,
          wizard_Intn_,
          state_->Wizard(),
          dofset_Intn_,
          state_->DofSet(),
          xfluid_timintapproach_,            // use the chosen approach as defined in the input file
          node_to_reconstr_method,
          reconstr_method_to_node,
          step_,
          xfluid_timint_check_interfacetips_)
  );

  {
    if(myrank_==0 and screen_out) IO::cout << "\t ...TransferVectorsToNewMap - TimeStepTransfer...";

    // --------------------------------------------
    // transfer of vectors from the old time step at the old interface position/dofset from t_n
    // to the current interface position/dofset at t_(n+1,i+1)
    //
    // vec_n(Gamma_n) -> vec_n(Gamma_n+1,i+1)

    //---------------------------------------------------------------
    // set old row state vectors at time step t^n that have to be updated to new interface position

    oldRowStateVectors.clear();
    newRowStateVectors.clear();

    oldRowStateVectors.push_back(veln_Intn_);
    newRowStateVectors.push_back(state_->veln_);

    if(timealgo_ == INPAR::FLUID::timeint_one_step_theta)
    {
      oldRowStateVectors.push_back(accn_Intn_);
      newRowStateVectors.push_back(state_->accn_);
    }
    else if(timealgo_ == INPAR::FLUID::timeint_bdf2)
    {
      oldRowStateVectors.push_back(velnm_Intn_);
      newRowStateVectors.push_back(state_->velnm_);
      oldRowStateVectors.push_back(accn_Intn_);
      newRowStateVectors.push_back(state_->accn_);
    }
    else
      dserror("check which vectors have to be reconstructed for non-OST and non-BDF2-scheme");

    XTimint_TransferVectorsBetweenSteps(
        xfluid_timeint,
        oldRowStateVectors,
        newRowStateVectors,
        dbcgids,
        false,
        screen_out);

  } // TransferDofsToNewMap

  if (xfluid_timintapproach_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_or_Proj_AND_GHOST_by_Proj_or_Copy_or_GP)
  {
    // project from another mesh, if possible (only for multimesh fluid)
    bool projection_success = XTimint_ProjectFromEmbeddedDiscretization(
        xfluid_timeint,
        newRowStateVectors,
        Teuchos::null,
        screen_out);

    if (!projection_success)
    {
      if (myrank_ == 0 and screen_out)
        IO::cout << "Reassigment of single-dof time integration approach after projection FAILED in some cases." << IO::endl;

      // we have nodes for which projection failed --> Correct the labels for those!
      XTimint_CorrectiveTransferVectorsBetweenSteps(
        xfluid_timeint,
        xfluid_timintapproach_,
        oldRowStateVectors,
        newRowStateVectors,
        dbcgids,
        screen_out);

      if (!xfluid_timeint->Get_NodeToDofMap_For_Reconstr(INPAR::XFEM::Xf_TimeInt_by_PROJ_from_DIS).empty())
        dserror("Even though projection failed, some nodes still demand projection. No alternatives found for e.g. %d",
            xfluid_timeint->Get_NodeToDofMap_For_Reconstr(INPAR::XFEM::Xf_TimeInt_by_PROJ_from_DIS).begin()->first);
    }
  }

  //------------------------------------------------------------------------------------
  //    GHOST PENALTY RECONSTRUCTION and/or SEMILAGRANGE RECONSTRUCTION necessary?
  //------------------------------------------------------------------------------------
  // decide if semi-Lagrangean back-tracking or ghost-penalty reconstruction has to be performed on any processor

  bool timint_ghost_penalty   = false;
  bool timint_semi_lagrangean = false;

  XTimint_GetReconstructStatus(xfluid_timeint, timint_ghost_penalty, timint_semi_lagrangean);

  //------------------------------------------------------------------------------------
  // STEP 2:               SEMILAGRANGE RECONSTRUCTION of std values
  //------------------------------------------------------------------------------------
  if(timint_semi_lagrangean)
  {
    Teuchos::RCP<Epetra_Vector> dispnpcol = Teuchos::null;
    Teuchos::RCP<Epetra_Vector> dispncol = Teuchos::null;

    if (alefluid_)
    {
      Teuchos::RCP<Epetra_Vector> dispnpcol = Teuchos::rcp(new Epetra_Vector(*DiscretisationXFEM()->InitialDofColMap()));
      Teuchos::RCP<Epetra_Vector> dispncol = Teuchos::rcp(new Epetra_Vector(*DiscretisationXFEM()->InitialDofColMap()));

      LINALG::Export(*dispnp_,*dispnpcol); //dispnp row->col
      LINALG::Export(*dispn_,*dispncol);   //dispn row->col
    }

    XTimint_SemiLagrangean(
        newRowStateVectors,             ///< vectors to be reconstructed
        newdofrowmap,                   ///< dofrowmap at current interface position
        oldRowStateVectors,             ///< vectors from which we reconstruct values (same order of vectors as in newRowStateVectors)
        dispnpcol,                      ///< displacement col - vector timestep n
        dispncol,                       ///< displacement row - vector timestep n+1
        &*dofcolmap_Intn_,              ///< dofcolmap at time and interface position t^n
        node_to_reconstr_method,        ///< reconstruction map for nodes and its dofsets
        screen_out                      ///< screen output?
        );

  } //SEMILAGRANGE RECONSTRUCTION of std values



  //------------------------------------------------------------------------------------
  // STEP 3:            GHOST PENALTY RECONSTRUCTION of ghost values
  //------------------------------------------------------------------------------------
  if(timint_ghost_penalty)
  {

    XTimint_GhostPenalty(
        newRowStateVectors,         ///< vectors to be reconstructed
        newdofrowmap,               ///< dofrowmap
        dbcgids,                    ///< dbc global ids
        screen_out                  ///< screen output?
        );

  }

  condition_manager_->ClearState();

  return;
}



/*----------------------------------------------------------------------*
 | Transfer vectors at current time-step t^(n+1) w.r.t dofset and       |
 | interface position from last iteration i to vectors w.r.t            |
 | current dofset and interface position (i+1)                          |
 | return, if increment step tranfer was successful!       schott 08/14 |
 *----------------------------------------------------------------------*/
bool FLD::XFluid::XTimint_DoIncrementStepTransfer(
    const bool screen_out,
    const bool firstcall_in_timestep
)
{

  const bool check_for_newton_restart =true;

  //------ CHANGING DOFSETS COMPARED TO LAST ITERATION? -----------

  // check for changing dofsets.
  // This is just required for new Newton increments to decide if a restart of the Newton has to be performed,
  // however, not for the first solve where the new interface position is given by the structural predictor and
  // at least one monolithic solve has to be performed before we can
  // decide if the Newton has to be restarted


  // MONOLITHIC XFSI
  // check if the dofmaps between last monolithic Newton iteration i and new Newton iteration i+1 changed in the fluid
  // dofmaps did not change when:
  //        1. the number of nodal dofsets for each node is the same for both iterations
  //        2. the time-integration identified respective nodal dofsets between Newton iterations,
  //           such that values of the nodal dofsets could be simply copied between the two iterations
  //           (note: between two Newton iterations with non-changing dofsets the ordering of respective ghost-dofsets can change
  //                  (as the cut cannot guarantee for the same order of ghost sets for slightly different interface positions).
  //                  Further a copy between a std dofset at one iteration and ghost dofsets at the other iteration can be reasonable,
  //                  in that case the dofsets did not change their meaning, however PERMUTATIONS of dofsets of single nodes
  //                  have to be taken into account, see PERMUTATIONS in fsi_xfem_monolithic)

  //---------------------------------------------------------------



  //---------------------------------------------------------------
  const Epetra_Map* newdofrowmap = discret_->DofRowMap();

  // all vectors that have to be transferred from old dofset to new dofset
  // vec_n+1(Gamma_n+1,i) -> vec_n+1(Gamma_n+1,i+1)
  std::vector<Teuchos::RCP<const Epetra_Vector> > rowStateVectors_npi;
  std::vector<Teuchos::RCP<Epetra_Vector> >       rowStateVectors_npip;

  // reconstruction map for nodes and its dofsets - how do we have to reconstruct the sinlge dofs
  std::map<int, std::vector<INPAR::XFEM::XFluidTimeInt> >             node_to_reconstr_method;
  std::map<INPAR::XFEM::XFluidTimeInt, std::map<int,std::set<int> > > reconstr_method_to_node;

  // vector of DOF-IDs which are Dirichlet BCs for ghost penalty reconstruction method
  Teuchos::RCP<std::set<int> > dbcgids = Teuchos::rcp(new std::set<int>());

  //------------------------------------------------------------------------------------
  // set interface state vectors for mesh coupling objects
  //------------------------------------------------------------------------------------
  condition_manager_->SetStateDisplacement();  // set idispnp, idispn and idispnpi vectors

  //------------------------------------------------------------------------------------
  // STEP 1: CopyDofsToNewMap and determine RECONSTRUCTION METHOD for missing values
  //------------------------------------------------------------------------------------
  //
  // REMARK:
  // * do this for row nodes only
  // * the cut information around the node should be available, since the cut is performed for col elements
  // * after transferring data from old interface position to new interface position the col vectors have to get
  //   exported from row vectors
  //------------------------------------------------------------------------------------

  INPAR::XFEM::XFluidTimeIntScheme timint_method;

  if(firstcall_in_timestep) // for the first iteration we allow the standard reconstruction method as we again reconstruct w.r.t t^n
    timint_method = xfluid_timintapproach_;
  else // for further iterations we just allow for simple copying and ghost-penalty reconstruction
    // for monolithic fsi and also for partitioned fsi it is the best not to allow semi-lagrangean
    timint_method = INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_AND_GHOST_by_Copy_or_GP;

  //-----------------------------time integration----------------------

  // create time integration class just locally not to keep pointers to dofset and wizard...
  Teuchos::RCP<XFEM::XFluidTimeInt> xfluid_timeint = Teuchos::rcp(new
      XFEM::XFluidTimeInt(
          true, // is_newton_increment_transfer?
          discret_,
          condition_manager_,
          wizard_Intnpi_,
          state_->Wizard(),
          dofset_Intnpi_,
          state_->DofSet(),
          timint_method,
          node_to_reconstr_method,
          reconstr_method_to_node,
          step_,
          xfluid_timint_check_interfacetips_)
  );

  {
    if(myrank_==0 and screen_out) IO::cout << "\t ...TransferVectorsToNewMap - IncrementStepTransfer...";

    // --------------------------------------------
    // transfer for the current iteration solution between last interface position of iteration i
    // and the current interface position at iteration i+1

    rowStateVectors_npi.clear();
    rowStateVectors_npip.clear();

    // transform the last Newton iteration
    rowStateVectors_npi.push_back(velnp_Intnpi_);
    rowStateVectors_npip.push_back(state_->velnp_);

    //Note: for reconstruction w.r.t last increment, do not use any semi-lagrangean approach
    XTimint_TransferVectorsBetweenSteps(
        xfluid_timeint,
        rowStateVectors_npi,
        rowStateVectors_npip,
        dbcgids,
        true, // fill the permutation map
        screen_out
        );

  }

  if (xfluid_timintapproach_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_or_Proj_AND_GHOST_by_Proj_or_Copy_or_GP)
  {
    // project from another mesh, if possible (only for multimesh fluid)
    bool projection_success = XTimint_ProjectFromEmbeddedDiscretization(
        xfluid_timeint,
        rowStateVectors_npip,
        Teuchos::null,
        screen_out);

    if (!projection_success)
    {
      if (myrank_ == 0 and screen_out)
        IO::cout << "Reassigment of single-dof time integration approach after projection FAILED in some cases." << IO::endl;

      // we have nodes for which projection failed --> correct the labels for those!
      XTimint_CorrectiveTransferVectorsBetweenSteps(
        xfluid_timeint,
        xfluid_timintapproach_,
        rowStateVectors_npi,
        rowStateVectors_npip,
        dbcgids,
        screen_out);

      if (!xfluid_timeint->Get_NodeToDofMap_For_Reconstr(INPAR::XFEM::Xf_TimeInt_by_PROJ_from_DIS).empty())
        dserror("Even though projection failed, some nodes still hold a projection label. No alternatives found for e.g. %d",
            xfluid_timeint->Get_NodeToDofMap_For_Reconstr(INPAR::XFEM::Xf_TimeInt_by_PROJ_from_DIS).begin()->first);
    }
  }

  //------------------------------------------------------------------------------------
  //    GHOST PENALTY RECONSTRUCTION and/or SEMILAGRANGE RECONSTRUCTION necessary?
  //------------------------------------------------------------------------------------
  // decide if semi-Lagrangean back-tracking or ghost-penalty reconstruction has to be performed on any processor

  bool timint_ghost_penalty   = false;
  bool timint_semi_lagrangean = false;

  XTimint_GetReconstructStatus(xfluid_timeint, timint_ghost_penalty, timint_semi_lagrangean);

  if( timint_semi_lagrangean )
  {
    if( firstcall_in_timestep ) // allow for semi-lagrangean in the first iteration
    {
      Teuchos::RCP<Epetra_Vector> dispnpcol = Teuchos::null;
      Teuchos::RCP<Epetra_Vector> dispncol = Teuchos::null;

      if (alefluid_)
      {
        Teuchos::RCP<Epetra_Vector> dispnpcol = Teuchos::rcp(new Epetra_Vector(*DiscretisationXFEM()->InitialDofColMap()));
        Teuchos::RCP<Epetra_Vector> dispncol = Teuchos::rcp(new Epetra_Vector(*DiscretisationXFEM()->InitialDofColMap()));

        LINALG::Export(*dispnp_,*dispnpcol); //dispnp row->col
        LINALG::Export(*dispn_,*dispncol);   //dispn row->col
      }

      XTimint_SemiLagrangean(
          rowStateVectors_npip,           ///< vectors to be reconstructed
          newdofrowmap,                   ///< dofrowmap at current interface position
          rowStateVectors_npi,            ///< vectors from which we reconstruct values (same order of vectors as in newRowStateVectors)
          dispnpcol,                      ///< displacement col - vector timestep n
          dispncol,                       ///< displacement row - vector timestep n+1
          &*dofcolmap_Intn_,              ///< dofcolmap at time and interface position t^n
          node_to_reconstr_method,        ///< reconstruction map for nodes and its dofsets
          screen_out                      ///< screen output?
          );

    }
    else
    {
      // How to perform a good prediction as startvalue when restarting the monolithic Newton is required
      // and simple copying is not possible???

      IO::cout << "check, how we can get the best predicted velnpip when simple copying + ghost penalty is not sufficient! " << IO::endl;

      // in this case SEMILAGRANGE is probably not reasonable as it is a mapping within the same timestep
      // reconstruct the missing values purely via Ghost-Penalty? GP-Faces sufficient? -> maybe use more faces
      dserror("using a Semi-lagrangean technique for reconstructing w.r.t last increment not reasonable, as the last increment is already an approximation to the actual solution at the same timestep!");

      return false; // apply the steady-state predictor in first time-step again instead, then we loose the information of the actual fluid predictor
    }
  }

  //------------------------------------------------------------------------------------
  // STEP 3:            GHOST PENALTY RECONSTRUCTION of ghost values
  //------------------------------------------------------------------------------------
  if(timint_ghost_penalty)
  {

    XTimint_GhostPenalty(
        rowStateVectors_npip,       ///< vectors to be reconstructed
        newdofrowmap,               ///< dofrowmap
        dbcgids,                    ///< dbc global ids
        screen_out                  ///< screen output?
    );

  }

  //------------------------------------------------------------------------------------
  // decide if the monolithic Newton has to be restarted, in case of the first iteration after a restart this information is
  // not used in the Newton loop
  //------------------------------------------------------------------------------------


  newton_restart_monolithic_ = false;

  if(check_for_newton_restart)
  {
    newton_restart_monolithic_ = XTimint_CheckForMonolithicNewtonRestart(
        timint_ghost_penalty,    ///< dofs have to be reconstructed via ghost-penalty reconstruction techniques
        timint_semi_lagrangean,  ///< dofs have to be reconstructed via semi-Lagrangean reconstruction techniques
        discret_,                ///< discretization
        dofset_Intnpi_,          ///< dofset last iteration
        state_->DofSet(),        ///< dofset current iteration
        screen_out               ///< screen output?
    );
  }

  condition_manager_->ClearState();

  return true;
}



/*----------------------------------------------------------------------*
 |  transfer vectors between two time-steps or Newton steps             |
 |                                                         schott 04/14 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::XTimint_TransferVectorsBetweenSteps(
    const Teuchos::RCP<XFEM::XFluidTimeInt> &                              xfluid_timeint,           ///< xfluid time integration class
    std::vector<Teuchos::RCP<const Epetra_Vector> >&                       oldRowStateVectors,       /// row map based vectors w.r.t old interface position
    std::vector<Teuchos::RCP<Epetra_Vector> >&                             newRowStateVectors,       /// row map based vectors w.r.t new interface position
    Teuchos::RCP<std::set<int> >                                           dbcgids,                  /// set of dof gids that must not be changed by ghost penalty reconstruction
    bool                                                                   fill_permutation_map,
    bool                                                                   screen_out
)
{
  const bool reconstruct_method_output = false;

  xfluid_timeint->TransferDofsToNewMap(oldRowStateVectors, newRowStateVectors, dbcgids);

  if(fill_permutation_map) permutation_map_ = xfluid_timeint->GetPermutationMap();

  if(myrank_==0 and screen_out) std::cout << " done\n" << std::flush;

  xfluid_timeint->SetAndPrintStatus(screen_out);

  if(reconstruct_method_output) xfluid_timeint->Output();
}

/*----------------------------------------------------------------------*
 |  transfer vectors between two time-steps or Newton steps             |
 |  (second run in case of failure in first attempt)       kruse 04/15  |
 *----------------------------------------------------------------------*/
void FLD::XFluid::XTimint_CorrectiveTransferVectorsBetweenSteps(
    const Teuchos::RCP<XFEM::XFluidTimeInt> &                              xfluid_timeint,           ///< xfluid time integration class
    const INPAR::XFEM::XFluidTimeIntScheme                                 xfluid_timintapproach,    /// xfluid_timintapproch
    std::vector<Teuchos::RCP<const Epetra_Vector> >&                       oldRowStateVectors,       ///< row map based vectors w.r.t old interface position
    std::vector<Teuchos::RCP<Epetra_Vector> >&                             newRowStateVectors,       ///< row map based vectors w.r.t new interface position
    Teuchos::RCP<std::set<int> >                                           dbcgids,                  ///< set of dof gids that must not be changed by ghost penalty reconstruction
    bool                                                                   screen_out                ///< output to screen
)
{
  std::map<int,std::set<int> >& reconstr_map = xfluid_timeint->Get_NodeToDofMap_For_Reconstr(INPAR::XFEM::Xf_TimeInt_by_PROJ_from_DIS);

  std::vector<int> failed_nodevec;
  failed_nodevec.reserve(reconstr_map.size());
  for (std::map<int,std::set<int> >::const_iterator in = reconstr_map.begin(); in != reconstr_map.end(); ++ in)
  {
    failed_nodevec.push_back(in->first);
  }

  xfluid_timeint->TransferDofsToNewMap(oldRowStateVectors, newRowStateVectors, dbcgids, failed_nodevec);

  xfluid_timeint->SetAndPrintStatus(screen_out);
}

/*----------------------------------------------------------------------*
 | decide if semi-Lagrangean back-tracking or ghost-penalty            |
 | reconstruction has to be performed on any processor    schott 08/14 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::XTimint_GetReconstructStatus(
    const Teuchos::RCP<XFEM::XFluidTimeInt> & xfluid_timeint, ///< xfluid time integration class
    bool & timint_ghost_penalty,         ///< do we have to perform ghost penalty reconstruction of ghost values?
    bool & timint_semi_lagrangean        ///< do we have to perform semi-Lagrangean reconstruction of standard values?
    )
{
  //------------------------------------------------------------------------------------
  // decide if semi-lagrangean back-tracking or ghost-penalty reconstruction has to be performed on any processor
  // if at least one proc has to do any reconstruction all procs has to call the routine

  int proc_timint_ghost_penalty   = 0;
  int proc_timint_semi_lagrangean = 0;

  if(xfluid_timeint == Teuchos::null) dserror("xfluid_timint_ - class not available here!");

  std::map<INPAR::XFEM::XFluidTimeInt, int>& reconstr_count =  xfluid_timeint->Get_Reconstr_Counts();

  std::map<INPAR::XFEM::XFluidTimeInt, int>::iterator it;

  if((it = reconstr_count.find(INPAR::XFEM::Xf_TimeInt_GHOST_by_GP)) != reconstr_count.end())
    proc_timint_ghost_penalty = it->second;
  if((it = reconstr_count.find(INPAR::XFEM::Xf_TimeInt_STD_by_SL)) != reconstr_count.end())
    proc_timint_semi_lagrangean = it->second;

  // parallel communication if at least one node has to do a semilagrangean backtracking or ghost penalty reconstruction
  int glob_timint_ghost_penalty   = 0;
  int glob_timint_semi_lagrangean = 0;

  discret_->Comm().SumAll(&proc_timint_ghost_penalty, &glob_timint_ghost_penalty, 1);
  discret_->Comm().SumAll(&proc_timint_semi_lagrangean, &glob_timint_semi_lagrangean, 1);


  //------------------------------------------------------------------------------------

  timint_ghost_penalty = (glob_timint_ghost_penalty>0);
  timint_semi_lagrangean = (glob_timint_semi_lagrangean>0);

  //------------------------------------------------------------------------------------
  return;
}



/*----------------------------------------------------------------------*
 | create DBC and free map and return their common extractor            |
 |                                                         schott 08/14 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::MapExtractor> FLD::XFluid::CreateDBCMapExtractor(
    const Teuchos::RCP< const std::set<int> >  dbcgids,                    ///< dbc global dof ids
    const Epetra_Map*                          dofrowmap                   ///< dofrowmap
)
{
  // create DBC and free map and build their common extractor

  // build map of Dirichlet DOFs
  int nummyelements = 0;
  int* myglobalelements = NULL;
  std::vector<int> dbcgidsv;
  if (dbcgids->size() > 0)
  {
    dbcgidsv.reserve(dbcgids->size());
    dbcgidsv.assign(dbcgids->begin(),dbcgids->end());
    nummyelements = dbcgidsv.size();
    myglobalelements = &(dbcgidsv[0]);
  }
  Teuchos::RCP<Epetra_Map> dbcmap
  = Teuchos::rcp(new Epetra_Map(-1, nummyelements, myglobalelements, dofrowmap->IndexBase(), dofrowmap->Comm()));

  // build the map extractor of Dirichlet-conditioned and free DOFs
  return Teuchos::rcp(new LINALG::MapExtractor(*dofrowmap, dbcmap));
}



/*----------------------------------------------------------------------*
 | create new dbc maps for ghost penalty reconstruction and             |
 | reconstruct value which are not fixed by DBCs           schott 08/14 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::XTimint_GhostPenalty(
    std::vector<Teuchos::RCP<Epetra_Vector> >& rowVectors,                 ///< vectors to be reconstructed
    const Epetra_Map*                          dofrowmap,                  ///< dofrowmap
    const Teuchos::RCP<const std::set<int> >   dbcgids,                    ///< dbc global ids
    const bool                                 screen_out                  ///< screen output?
)
{
  if(myrank_==0 and screen_out) std::cout << "\t ...Ghost Penalty Reconstruction..." << std::endl;

  //----------------------------------------
  // object holds maps/subsets for DOFs subjected to Dirichlet BCs
  // which will not be modified by the ghost-penalty reconstruction
  Teuchos::RCP<LINALG::MapExtractor> ghost_penaly_dbcmaps = CreateDBCMapExtractor(dbcgids, dofrowmap);

  //----------------------------------------
  // perform ghost-penalty reconstruction for all vectors
  for(std::vector<Teuchos::RCP<Epetra_Vector> >::iterator vecs_it = rowVectors.begin();
      vecs_it != rowVectors.end();
      vecs_it++)
  {
    // reconstruct values using ghost penalty approach
    XTimint_ReconstructGhostValues(*vecs_it, ghost_penaly_dbcmaps, screen_out);
  }


  if(myrank_==0 and screen_out) std::cout << " done\n" << std::flush;

  return;
}

/*----------------------------------------------------------------------*
 |  reconstruct ghost values via ghost penalties           schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::XTimint_ReconstructGhostValues(
    Teuchos::RCP<Epetra_Vector>                vec,                        ///< vector to be reconstructed
    Teuchos::RCP<LINALG::MapExtractor>         ghost_penaly_dbcmaps,       ///< which dofs are fixed during the ghost-penalty reconstruction?
    const bool                                 screen_out                  ///< screen output?
)
{
  discret_->Comm().Barrier();

  TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XTimint_ReconstructGhostValues" );

  // ---------------------------------------------- setup solver

  Teuchos::RCP<Teuchos::ParameterList> solverparams = Teuchos::rcp(new Teuchos::ParameterList);

#if(0) // use direct solver
  solverparams->set("solver","umfpack");
#else  // use iterative solver
  Teuchos::ParameterList& azlist = solverparams->sublist("Aztec Parameters");
  azlist.set<int>("reuse", 0);
  azlist.set<double>("AZ_tol", 1.0e-12);
  //azlist.set("AZ_kspace",inparams.get<int>("AZSUB")); //AZSUB
  azlist.set("AZ_overlap",0);

#if 0 // suppress output
  azlist.set("AZ_output",0);
#else
  azlist.set("AZ_output",50);
#endif
  solverparams->sublist("IFPACK Parameters");

  solverparams->set("solver","aztec");
#endif

  Teuchos::RCP<LINALG::Solver> solver_gp =
      Teuchos::rcp(new LINALG::Solver( solverparams, discret_->Comm(), DRT::Problem::Instance()->ErrorFile()->Handle() ));

  // ---------------------------------------------- new matrix and vectors

  // TODO: use the matrix more than once when this step becomes expensive!

  // get a good estimate for the non-zeros!
  // create a map (Dirichlet values get ones, non-Dirichlet values get the 162)

  int numentries_dbc_row = 1;
  int numentries_ghost_penalty_row = 162;

  std::vector<int> numentries( state_->xfluiddofrowmap_->NumMyElements() );

  const Epetra_Map & rowmap = *state_->xfluiddofrowmap_;
  const Epetra_Map & condmap = *(ghost_penaly_dbcmaps->CondMap());

  for ( unsigned i=0; i<numentries.size(); ++i )
  {
    int gid = rowmap.GID(i);
    int dbclid = condmap.LID(gid);
    if (dbclid < 0) // non-dbc-row
      numentries[i] = numentries_ghost_penalty_row;
    else // dbc-row
      numentries[i] = numentries_dbc_row;
  }

  // note: we use explicitdirichlet =  false, as we don't want to create a new sysmat when applying Dirichlet bcs
  // note: savegraph = true as we assemble the matrix more than once
  Teuchos::RCP<LINALG::SparseMatrix> sysmat_gp =
      Teuchos::rcp(new LINALG::SparseMatrix(*state_->xfluiddofrowmap_,numentries,false,true,LINALG::SparseMatrix::FE_MATRIX));


  Teuchos::RCP<Epetra_Vector> zeros_gp        = LINALG::CreateVector(*state_->xfluiddofrowmap_,true);
  Teuchos::RCP<Epetra_Vector> residual_gp     = LINALG::CreateVector(*state_->xfluiddofrowmap_,true);
  Teuchos::RCP<Epetra_Vector> incvel_gp       = LINALG::CreateVector(*state_->xfluiddofrowmap_,true);

  dtsolve_  = 0.0;
  dtele_    = 0.0;
  dtfilter_ = 0.0;

  if (myrank_ == 0 and screen_out)
  {
    printf("\n+++++++++++++++++++++ Gradient Penalty Ghost value reconstruction ++++++++++++++++++++++++++++\n");
  }

  // do only one solve (as the system is linear!)
  {
    discret_->Comm().Barrier();

    // get cpu time
    const double tcpu=Teuchos::Time::wallTime();

    // evaluate routine
    AssembleMatAndRHS_GradientPenalty(ghost_penaly_dbcmaps, sysmat_gp, residual_gp, vec);

    // end time measurement for element
    dtele_=Teuchos::Time::wallTime()-tcpu;
  }

  // blank residual DOFs which are on Dirichlet BC
  // We can do this because the values at the dirichlet positions
  // are not used anyway.
  // We could avoid this though, if velrowmap_ and prerowmap_ would
  // not include the dirichlet values as well. But it is expensive
  // to avoid that.

  {
    discret_->Comm().Barrier();

    TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XTimint_ReconstructGhostValues::ghost_penaly_dbcmaps->InsertCondVector" );

    ghost_penaly_dbcmaps->InsertCondVector(ghost_penaly_dbcmaps->ExtractCondVector(zeros_gp), residual_gp);
  }

  //--------- Apply Dirichlet boundary conditions to system of equations
  //          residual displacements are supposed to be zero at
  //          boundary conditions
  incvel_gp->PutScalar(0.0);

  {
    discret_->Comm().Barrier();
    TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XTimint_ReconstructGhostValues::ApplyDirichlettoSystem" );

    Epetra_Time tf(discret_->Comm());

    LINALG::ApplyDirichlettoSystem(sysmat_gp,incvel_gp,residual_gp,zeros_gp,*(ghost_penaly_dbcmaps->CondMap()));

    //IO::cout << "fluid time : ApplyDirichlettoSystem " << tf.ElapsedTime() << IO::endl;
  }

  //-------solve for residual displacements to correct incremental displacements
  {
    discret_->Comm().Barrier();

    TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XTimint_ReconstructGhostValues::Solve" );

    // get cpu time
    const double tcpusolve=Teuchos::Time::wallTime();

    solver_gp->Solve(sysmat_gp->EpetraOperator(),incvel_gp,residual_gp,true, true);

    // end time measurement for solver
    dtsolve_ = Teuchos::Time::wallTime()-tcpusolve;
  }

  // -------------------------------------------------------------------
  // update velocity and pressure values by increments
  // -------------------------------------------------------------------
  vec->Update(1.0,*incvel_gp,1.0);

  return;
} // ReconstructGhostValues


/*----------------------------------------------------------------------*
 |  reconstruct standard values via semi-Lagrangean method schott 08/14 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::XTimint_SemiLagrangean(
    std::vector<Teuchos::RCP<Epetra_Vector> >&               newRowStateVectors,      ///< vectors to be reconstructed
    const Epetra_Map*                                        newdofrowmap,            ///< dofrowmap at current interface position
    std::vector<Teuchos::RCP<const Epetra_Vector> >&         oldRowStateVectors,      ///< vectors from which we reconstruct values (same order of vectors as in newRowStateVectors)
    Teuchos::RCP<Epetra_Vector>                              dispn,                   ///< displacement initial col - vector timestep n //set to Teuchos::null if no ale displacments
    Teuchos::RCP<Epetra_Vector>                              dispnp,                  ///< displacement initial col - vector timestep n+1 //if Teuchos::null ... --> no ale displacments
    const Epetra_Map*                                        olddofcolmap,            ///< dofcolmap at time and interface position t^n
    std::map<int, std::vector<INPAR::XFEM::XFluidTimeInt> >& node_to_reconstr_method, ///< reconstruction map for nodes and its dofsets
    const bool                                               screen_out               ///< screen output?
)
{

  if(myrank_==0 and screen_out) std::cout << "\t ...SemiLagrangean Reconstruction...";

  Teuchos::RCP<XFEM::MeshCoupling> mc_coupl = condition_manager_->GetMeshCoupling(mc_idx_);
  Teuchos::RCP<DRT::Discretization> bounddis = mc_coupl->GetCutterDis();

  condition_manager_->SetStateDisplacement();

  //--------------------------------------------------------
  // export veln row vector from t^n to a col vector

  Teuchos::RCP<Epetra_Vector> veln_col = Teuchos::rcp(new Epetra_Vector(*olddofcolmap,true));
  LINALG::Export(*veln_Intn_,*veln_col);

  //--------------------------------------------------------
  // export row vectors from t^n to col vectors
  // Important: export the vectors used for Semi-Lagrangean method after transfer between interface processors above
  std::vector<Teuchos::RCP<Epetra_Vector> > oldColStateVectorsn;

  for(std::vector<Teuchos::RCP<const Epetra_Vector> >::iterator vec_it = oldRowStateVectors.begin();
      vec_it != oldRowStateVectors.end();
      vec_it++)
  {
    Teuchos::RCP<Epetra_Vector> vec_col = Teuchos::rcp(new Epetra_Vector(*olddofcolmap,true));
    LINALG::Export(**vec_it,*vec_col);
    oldColStateVectorsn.push_back(vec_col);
  }


  // TODO: set this param
  int totalitnumFRS_ = 0;
  int itemaxFRS_ = 5;
  Teuchos::RCP<XFEM::XFLUID_STD> timeIntStd_ = Teuchos::null;

  INPAR::XFEM::XFluidTimeInt xfemtimeint_ = INPAR::XFEM::Xf_TimeInt_STD_by_SL;

  if (totalitnumFRS_==0) // construct time int classes once every time step
  {
    // basic time integration data
    Teuchos::RCP<XFEM::XFLUID_TIMEINT_BASE> timeIntData = Teuchos::null;

    timeIntData = Teuchos::rcp(new XFEM::XFLUID_TIMEINT_BASE(
        discret_,
        bounddis,
        wizard_Intn_,
        state_->Wizard(),
        dofset_Intn_,
        state_->DofSet(),
        oldColStateVectorsn,
        dispn,
        dispnp,
        *dofcolmap_Intn_,
        *newdofrowmap,
        Teuchos::null));

    //Safty check (both displacements have to exist or not --> based on that ale fluid is activated)
    if ((dispn != Teuchos::null and dispnp == Teuchos::null) or (dispn == Teuchos::null and dispnp != Teuchos::null))
      dserror("FLD::XFluid::XTimint_SemiLagrangean: dispn or dispnp indicate ale fluid!");

    switch (xfemtimeint_)
    {
    case INPAR::XFEM::Xf_TimeInt_STD_by_SL:
    {
      // time integration data for standard dofs, semi-lagrangean approach
      timeIntStd_ = Teuchos::rcp(new XFEM::XFLUID_SemiLagrange(
          *timeIntData,
          node_to_reconstr_method,
          xfemtimeint_,
          veln_col,
          dta_,
          theta_,
          true));
      break;
    }
    default:
    {
      dserror("unknown recomputation approach in XFEM time integration not implemented");
      break;
    }
    }

    totalitnumFRS_++;

    timeIntStd_->type(totalitnumFRS_,itemaxFRS_); // update algorithm handling
    timeIntStd_->compute(newRowStateVectors);     // call computation

  } //totalit

  condition_manager_->ClearState();

  if(myrank_==0) std::cout << " done\n" << std::flush;

  return;
}

/*----------------------------------------------------------------------*
 | calculate lift&drag forces                              schott 01/15 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::LiftDrag() const
{
  // initially check whether computation of lift and drag values is required
  if (params_->get<bool>("LIFTDRAG"))
  {
    condition_manager_->LiftDrag(step_,time_);
  }
}


/// return time integration factor
double FLD::XFluid::TimIntParam() const
{
  double retval = 0.0;
  switch (TimIntScheme())
  {
  case INPAR::FLUID::timeint_afgenalpha:
  case INPAR::FLUID::timeint_npgenalpha:
    // this is the interpolation weight for quantities from last time step
    retval = 1.0 - alphaF_;
  break;
  case INPAR::FLUID::timeint_one_step_theta:
    // this is the interpolation weight for quantities from last time step
    retval = 0.0;
  break;
  case INPAR::FLUID::timeint_bdf2:
    // this is the interpolation weight for quantities from last time step
    retval = 0.0;
  break;
  case INPAR::FLUID::timeint_stationary:
    // this is the interpolation weight for quantities from last time step
    retval = 0.0;
  break;
  default:
    dserror("Unknown time integration scheme");
  break;
  }
  return retval;
}

/*----------------------------------------------------------------------*
 |  write solution output                                  schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::Output()
{
  const bool write_restart_data = step_!=0 and uprestart_ != 0 and step_%uprestart_ == 0;

  //---------------------------------- GMSH SOLUTION OUTPUT (solution fields for pressure, velocity) ------------------------

  // write gmsh-output for solution fields
  // solution output
  output_service_->GmshSolutionOutput("SOL", step_, state_);

  //---------------------------------- GMSH DISCRET OUTPUT (extended output for EOS) ------------------------
  output_service_->GmshOutputEOS(step_,edgestab_);

  //---------------------------------- PARAVIEW SOLUTION OUTPUT (solution fields for pressure, velocity) ------------------------

  if (step_%upres_ == 0)
  {
    output_service_->Output(
      step_,
      time_,
      write_restart_data,
      state_,
      dispnp_,
      gridvnp_);
  }


  return;
}


/*----------------------------------------------------------------------*
 |  set an initial flow field                              schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::SetInitialFlowField(
  const INPAR::FLUID::InitialField initfield,
  const int startfuncno
  )
{

  // no set initial flow field for restart
  if(step_ != 0) return;

  if(myrank_ == 0) std::cout << "SetInitialFlowField " << std::endl;

  // initial field by (undisturbed) function (init==2)
  // or disturbed function (init==3)
  if (initfield == INPAR::FLUID::initfield_field_by_function/* or
      initfield == INPAR::FLUID::initfield_disturbed_field_from_function*/)
  {
    if(myrank_ == 0) std::cout << "SetInitialFlowField with function number " << startfuncno << std::endl;

    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = discret_->lRowNode(lnodeid);
      // the set of degrees of freedom associated with the node
      const std::vector<int> nodedofset = discret_->Dof(0,lnode);

      if (nodedofset.size()!=0)
      {
          for(int dof=0;dof<(int)nodedofset.size();++dof)
          {
            int gid = nodedofset[dof];

            double initialval=DRT::Problem::Instance()->Funct(startfuncno-1).Evaluate(dof%4,lnode->X(),time_,NULL);
            state_->velnp_->ReplaceGlobalValues(1,&initialval,&gid);
          }
      }
    }

    // initialize veln_ as well.
    state_->veln_->Update(1.0,*state_->velnp_ ,0.0);
    state_->velnm_->Update(1.0,*state_->velnp_ ,0.0);

    state_->accnp_->PutScalar(0.0);
    state_->accn_->PutScalar(0.0);

  }
  // special initial function: Beltrami flow (3-D)
  else if (initfield == INPAR::FLUID::initfield_beltrami_flow)
  {
    const Epetra_Map* dofrowmap = discret_->DofRowMap();

    int err = 0;

    const int npredof = numdim_;

    double         p;
    std::vector<double> u  (numdim_);
    std::vector<double> xyz(numdim_);

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
      std::vector<int> nodedofset = discret_->Dof(0,lnode);

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
      p = -a*a/2.0 * dens *
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
        err += state_->velnp_->ReplaceMyValues(1,&(u[nveldof]),&lid);
        err += state_->veln_ ->ReplaceMyValues(1,&(u[nveldof]),&lid);
        err += state_->velnm_->ReplaceMyValues(1,&(u[nveldof]),&lid);
      }

      // set initial pressure
      const int gid = nodedofset[npredof];
      int lid = dofrowmap->LID(gid);
      err += state_->velnp_->ReplaceMyValues(1,&p,&lid);
      err += state_->veln_ ->ReplaceMyValues(1,&p,&lid);
      err += state_->velnm_->ReplaceMyValues(1,&p,&lid);
    } // end loop nodes lnodeid

    if (err!=0) dserror("dof not on proc");
  }
  //----------------------------------------------------------------------------------------------
  // flame-vortex interaction problem: two counter-rotating vortices (2-D) moving the flame front
  //----------------------------------------------------------------------------------------------
  else if(initfield == INPAR::FLUID::initfield_flame_vortex_interaction)
  {
    //TODO: shift this function to the condition-manager!

    //Only supported for 1 levelset so far.
    if(condition_manager_->NumLevelSetCoupling() != 1)
      dserror("There is either no LevelSetCoupling or more than 1. Exactly 1 is expected and supported at this point!");

    Teuchos::RCP<XFEM::LevelSetCoupling>           levelset_condition   = condition_manager_->GetLevelSetCoupling("XFEMLevelsetCombustion");
    Teuchos::RCP<XFEM::LevelSetCouplingCombustion> combustion_condition = Teuchos::rcp_dynamic_cast<XFEM::LevelSetCouplingCombustion>(levelset_condition, true);


    // vector of DOF-IDs which are Dirichlet BCs for ghost penalty reconstruction method
    Teuchos::RCP<std::set<int> > dbcgids = Teuchos::rcp(new std::set<int>());

    const Teuchos::RCP<GEO::CutWizard>   & wizard = state_->Wizard();
    const Teuchos::RCP<XFEM::XFEMDofSet> & dofset = state_->DofSet();
    const Epetra_Map * dofrowmap = dofset->DofRowMap();

    //------------------------
    // get material parameters
    //------------------------
    // arbitrarily take first node on this proc
    DRT::Node* lnode = discret_->lRowNode(0);
    // get list of adjacent elements of the first node
    DRT::Element** elelist = lnode->Elements();
    DRT::Element* ele = elelist[0]; // (arbitrary!) first element
    // get material from first (arbitrary!) element adjacent to this node
    const Teuchos::RCP<MAT::Material> material = ele->Material();
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
    const double dens_u = mat0->Density(); // outside, master, (i for i<j convention)
    if (dens_u != 1.161) dserror("unburnt density should be 1.161 for the 'flame-vortex-interaction' case");
    const double dens_b = mat1->Density(); // inside, slave, (j for i<j convention)
    if (dens_b != 0.157) dserror("burnt density should be 0.157 for the 'flame-vortex-interaction' case");

    double dens = dens_u;

    // get interface specific parameter from condition manager and underlying conditions
    const XFEM::EleCoupCond & coup_cond = combustion_condition->GetCouplingCondition(ele->Id());
    double flamespeed = coup_cond.second->GetDouble("laminar_flamespeed");

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
    if (flamespeed != 1.0) dserror("flame speed should be 1.0 for the 'flame-vortex-interaction' case");
    // vortex strength C (scaled by laminar burning velocity)
    const double C = 70.0*flamespeed;
    // (squared) vortex radius R
    const double R_squared = 16.0;

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

      // get the node from the cut wizard
      const int gid=lnode->Id();
      GEO::CUT::Node * cut_node = wizard->GetNode(gid);

      // ask for the number of dofsets
      const int numDofSets = cut_node->NumDofSets();

      const std::vector<Teuchos::RCP<GEO::CUT::NodalDofSet> > & nodaldofsets = cut_node->NodalDofSets();

      // set values just for the standard dofset, all ghost sets are determined by a ghost-penalty time integration solve!
      for(int i=0; i<numDofSets; ++i)
      {
        //-------------------------------------------
        // STOP FOR GHOSTSETS
        //-------------------------------------------
        if(!nodaldofsets[i]->Is_Standard_DofSet())
          continue; // do nothing for ghost dofsets!

        //-------------------------------------------
        // just FOR STD SETS
        //-------------------------------------------
        GEO::CUT::Point::PointPosition pos = nodaldofsets[i]->Position();

        //-------------------------------------------
        // get values for flamevortex interaction
        //-------------------------------------------

        // compute preliminary values for both vortices
        const double r_squared_left  = ((xyz(0)-xyz0_left(0))*(xyz(0)-xyz0_left(0))
            +(xyz(1)-xyz0_left(1))*(xyz(1)-xyz0_left(1)))/R_squared;
        const double r_squared_right = ((xyz(0)-xyz0_right(0))*(xyz(0)-xyz0_right(0))
            +(xyz(1)-xyz0_right(1))*(xyz(1)-xyz0_right(1)))/R_squared;

        //----------------------------------------
        // set density with respect to flame front
        //----------------------------------------
        const double tmp_vortex = -0.5*(C*C/R_squared)*(exp(-r_squared_left) + exp(-r_squared_right));

        if (pos == GEO::CUT::Point::inside) // plus/burnt domain -> burnt material ( GEO::CUT::Position is inside ) / slave side
        {
          dens = dens_b;
          pres = 0.0; // matching the zero pressure condition at outflow
        }
        else if( pos == GEO::CUT::Point::outside )// minus/unburnt domain -> unburnt material / master side
        {
          dens = dens_u;
          double jump_dens_inverse = 1.0/dens_u - 1.0/dens_b;
          pres = -flamespeed*flamespeed * dens_u*dens_u * jump_dens_inverse;
        }
        else
        {
          dserror("what to do now?");
        }
        pres += tmp_vortex;

        //----------------------------------------------
        // compute components of initial velocity vector
        //----------------------------------------------
#if(1)
        vel(0) = (C/R_squared)*(-(xyz(1)-xyz0_left(1))*exp(-r_squared_left/2.0)
            +(xyz(1)-xyz0_right(1))*exp(-r_squared_right/2.0));
        vel(1) = (C/R_squared)*( (xyz(0)-xyz0_left(0))*exp(-r_squared_left/2.0)
            -(xyz(0)-xyz0_right(0))*exp(-r_squared_right/2.0))
            + flamespeed*dens_u/dens;
        // 2D problem -> vel_z = 0.0
        vel(2) = 0.0;
#else // without vortices (just for testing)
        vel(0) = 0.0;
        vel(1) = flamespeed*dens_u/dens;
        // 2D problem -> vel_z = 0.0
        vel(2) = 0.0;
#endif
        // access standard FEM dofset (3 x vel + 1 x pressure) to get std-dof IDs for this node
        std::vector<int> std_dofs;
        dofset->Dof(std_dofs, lnode, i );

        //-----------------------------------------
        // set components of initial velocity field
        //-----------------------------------------
        for(int idim=0;idim<nsd+1;idim++)
        {
          const int gid = std_dofs[idim];
          //local node id
          int lid = dofrowmap->LID(gid);
          if(idim==3)
          { // pressure dof
            err += state_->velnp_->ReplaceMyValues(1,&pres,&lid);
          }
          else
          { // velocity dof
            err += state_->velnp_->ReplaceMyValues(1,&vel(idim),&lid);
          }

          // set Dirichlet BC for ghost penalty reconstruction
          if (dbcgids != Teuchos::null)
            (*dbcgids).insert(gid);
          if (err!=0) dserror("dof not on proc");

        }
      } // loop nodal dofsets
    } // end loop nodes lnodeid


    // reconstruct ghost values / use the ghost penalty reconstruction technique as used within the XFEM time integration
    std::vector<Teuchos::RCP<Epetra_Vector> > rowStateVectors_npip;
    rowStateVectors_npip.push_back(state_->velnp_);

    XTimint_GhostPenalty(
        rowStateVectors_npip,       ///< vectors to be reconstructed
        dofrowmap,                  ///< dofrowmap
        dbcgids,                    ///< dbc global ids
        true                        ///< screen output?
    );

    // set also veln and velnm
    // initialize veln_ and velnm_ as well.
    state_->veln_ ->Update(1.0,*state_->velnp_ ,0.0);
    state_->velnm_->Update(1.0,*state_->velnp_ ,0.0);
  }
  else
  {
    dserror("Only initial fields auch as a zero field, initial fields by (un-)disturbed functions, flamevortes and Beltrami flow!");
  }

  //---------------------------------- GMSH START OUTPUT (reference solution fields for pressure, velocity) ------------------------

  // write gmsh-output for start fields
  output_service_->GmshSolutionOutputPrevious( "START", step_, state_);

  return;
} // end SetInitialFlowField


void FLD::XFluid::WriteAccess_GeometricQuantities(
    Teuchos::RCP<Epetra_Vector> &      scalaraf,
    Teuchos::RCP<Epetra_MultiVector> & smoothed_gradphiaf,
    Teuchos::RCP<Epetra_Vector> &      curvatureaf
)
{
  condition_manager_->WriteAccess_GeometricQuantities(scalaraf, smoothed_gradphiaf, curvatureaf);
}

void FLD::XFluid::ExportGeometricQuantities()
{
  condition_manager_->ExportGeometricQuantities();
}

const Teuchos::RCP<const Epetra_Vector> FLD::XFluid::GetTransportVelocity()
{
  //Only supported for 1 levelset so far.
  if(condition_manager_->NumLevelSetCoupling() != 1)
    dserror("There is either no LevelSetCoupling or more than 1. Exactly 1 is expected and supported at this point!");

  Teuchos::RCP<XFEM::LevelSetCouplingTwoPhase> levelsetcoupling_twophase =
      Teuchos::rcp_dynamic_cast<XFEM::LevelSetCouplingTwoPhase>(condition_manager_->GetLevelSetCoupling(0), true);

  const Teuchos::RCP<const Epetra_Vector> & transport_velocity =
      levelsetcoupling_twophase->ComputeTransportVelocity(
          state_->Wizard(),
          state_->velnp_);

  return transport_velocity;
}

// -------------------------------------------------------------------
// set general fluid parameter (AE 01/2011)
// -------------------------------------------------------------------
void FLD::XFluid::SetDirichletNeumannBC()
{

  Teuchos::ParameterList eleparams;

  // other parameters needed by the elements
  eleparams.set("total time",time_);

  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("velaf",state_->velnp_);
  // predicted dirichlet values
  // velnp then also holds prescribed new dirichlet values
  discret_->EvaluateDirichlet(eleparams,state_->velnp_,Teuchos::null,Teuchos::null,Teuchos::null,state_->dbcmaps_);

  discret_->ClearState();

  if (alefluid_)
  {
    discret_->SetState("dispnp",state_->dispnp_);
  }

  // set thermodynamic pressure
  eleparams.set("thermodynamic pressure",thermpressaf_);

  state_->neumann_loads_->PutScalar(0.0);
  discret_->SetState("scaaf",state_->scaaf_);

  XFEM::EvaluateNeumann( eleparams, discret_, state_->neumann_loads_ );

  discret_->ClearState();

}


void FLD::XFluid::AssembleMatAndRHS()
{

}



/*----------------------------------------------------------------------*
 * Explicit predictor                                   rasthofer 12/13 |
 *----------------------------------------------------------------------*/
void FLD::XFluid::ExplicitPredictor()
{

  if(discret_->Comm().MyPID()==0)
  {
    printf("fluid: using explicit predictor %s",predictor_.c_str());
  }

  if (predictor_=="steady_state")
  {
    // steady state predictor
    //
    //       n+1    n
    //      u    = u
    //       (0)
    //
    //  and
    //
    //       n+1    n
    //      p    = p
    //       (0)

    // this has already been done in TimeUpdate()
  }
  else if(predictor_=="zero_acceleration")
  {
    // zero acceleration predictor
    //
    //       n+1    n                   n
    //      u    = u  + (1-gamma)*dt*acc
    //       (0)
    //
    //  and
    //
    //       n+1    n
    //      p    = p
    //       (0)
    //
    state_->velnp_->Update(1.0,*state_->veln_,0.0);

    // split between acceleration and pressure
    Teuchos::RCP<Epetra_Vector> inc = state_->velpressplitter_->ExtractOtherVector(state_->accn_);
    inc->Scale((1.0-theta_)*dta_);

    state_->velpressplitter_->AddOtherVector(inc,state_->velnp_);
  }
  else if(predictor_=="constant_acceleration")
  {
    // constant acceleration predictor
    //
    //       n+1    n         n
    //      u    = u  + dt*acc
    //       (0)
    //
    //  and
    //
    //       n+1    n
    //      p    = p
    //       (0)
    //
    state_->velnp_->Update(1.0,*state_->veln_,0.0);

    Teuchos::RCP<Epetra_Vector> inc = state_->velpressplitter_->ExtractOtherVector(state_->accn_);
    inc->Scale(dta_);

    state_->velpressplitter_->AddOtherVector(inc,state_->velnp_);
  }
  else if(predictor_=="constant_increment")
  {
    dserror("not supported for XFEM as we need to transform also velnm? Maybe it is possible! Check this!");

    // constant increment predictor
    //
    //       n+1      n    n-1
    //      u    = 2*u  - u
    //       (0)
    //
    //  and
    //
    //       n+1    n
    //      p    = p
    //       (0)
    //
    state_->velnp_->Update(1.0,*state_->veln_,0.0);

    Teuchos::RCP<Epetra_Vector> un  = state_->velpressplitter_->ExtractOtherVector(state_->veln_ );
    Teuchos::RCP<Epetra_Vector> unm = state_->velpressplitter_->ExtractOtherVector(state_->velnm_);
    unm->Scale(-1.0);

    state_->velpressplitter_->AddOtherVector(un ,state_->velnp_);
    state_->velpressplitter_->AddOtherVector(unm,state_->velnp_);
  }
  else if(predictor_=="explicit_second_order_midpoint")
  {
    // the conventional explicit second order predictor (assuming constant dt)
    // also known as leapfrog integration
    /*
    //                        /          n    n-1 \
    //       n+1    n        |      n   u  - u     |
    //      u    = u  + dt * | 2*acc  - ---------  |
    //       (0)             |             dt      |
    //                        \                   /
    // respectively
    //
    //       n+1    n-1               n
    //      u    = u    + 2 * dt * acc
    //       (0)
    //
    //  and
    //
    //       n+1    n
    //      p    = p
    //       (0)
    */
    state_->velnp_->Update(1.0,*state_->veln_,0.0);

    // split between acceleration and pressure
    Teuchos::RCP<Epetra_Vector> unm = state_->velpressplitter_->ExtractOtherVector(state_->velnm_);
    Teuchos::RCP<Epetra_Vector> an  = state_->velpressplitter_->ExtractOtherVector(state_->accn_ );

    unm->Update(2.0*dta_,*an,1.0);

    state_->velpressplitter_->InsertOtherVector(unm,state_->velnp_);
  }
  else
    dserror("Unknown fluid predictor %s", predictor_.c_str());

  if(discret_->Comm().MyPID()==0)
  {
    printf("\n");
  }

  return;
}


/*------------------------------------------------------------------------------------------------*
 |
 *------------------------------------------------------------------------------------------------*/
void FLD::XFluid::PredictTangVelConsistAcc()
{
  // message to screen
  if(discret_->Comm().MyPID()==0)
  {
    std::cout << "fluid: doing TangVel predictor" << std::endl;
  }

  // total time required for evaluation of Dirichlet conditions
  Teuchos::ParameterList eleparams;
  eleparams.set("total time",time_);

  // initialize
  state_->velnp_->Update(1.0, *state_->veln_, 0.0);
  state_->accnp_->Update(1.0, *state_->accn_, 0.0);
  state_->incvel_->PutScalar(0.0);

  // for solution increments on Dirichlet boundary
  Teuchos::RCP<Epetra_Vector> dbcinc
    = LINALG::CreateVector(*(discret_->DofRowMap()), true);

  // copy last converged solution
  dbcinc->Update(1.0, *state_->veln_, 0.0);

  // get Dirichlet values at t_{n+1}
  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("velnp",state_->velnp_);

  // predicted Dirichlet values
  // velnp_ then also holds prescribed new dirichlet values
  discret_->EvaluateDirichlet(eleparams,state_->velnp_,Teuchos::null,Teuchos::null,Teuchos::null);

  // subtract the displacements of the last converged step
  // DBC-DOFs hold increments of current step
  // free-DOFs hold zeros
  dbcinc->Update(-1.0, *state_->veln_, 1.0);

  // -------------------------------------------------------------------
  // compute residual forces residual_ and stiffness sysmat_
  // at velnp_, etc which are unchanged

  // -------------------------------------------------------------------
  // set old part of righthandside
  SetOldPartOfRighthandside();

  // -------------------------------------------------------------------
  // evaluate Dirichlet and Neumann boundary conditions
  SetDirichletNeumannBC();

  // -------------------------------------------------------------------
  // assemble matrix and rhs based on the last interface position (note, this is done before a new state class is created after performing the predictor!)
  AssembleMatAndRHS( 1 );


  // add linear reaction forces to residual
  // linear reactions
  Teuchos::RCP<Epetra_Vector> freact
    = LINALG::CreateVector(*(discret_->DofRowMap()), true);
  state_->sysmat_->Multiply(false, *dbcinc, *freact);

  // add linear reaction forces due to prescribed Dirichlet BCs
  state_->residual_->Update(1.0, *freact, 1.0);

  // extract reaction forces
  freact->Update(1.0, *state_->residual_, 0.0);
  state_->dbcmaps_->InsertOtherVector(state_->dbcmaps_->ExtractOtherVector(state_->zeros_), freact);

  // blank residual at DOFs on Dirichlet BC
  state_->dbcmaps_->InsertCondVector(state_->dbcmaps_->ExtractCondVector(state_->zeros_), state_->residual_);

  // apply Dirichlet BCs to system of equations
  state_->incvel_->PutScalar(0.0);
  state_->sysmat_->Complete();
  LINALG::ApplyDirichlettoSystem(state_->sysmat_, state_->incvel_, state_->residual_,
                                 Teuchos::null, state_->zeros_, *(state_->dbcmaps_->CondMap()));

  // solve for incvel_
  solver_->Solve(state_->sysmat_->EpetraOperator(), state_->incvel_, state_->residual_, true, true);

  // set Dirichlet increments in solution increments
  state_->incvel_->Update(1.0, *dbcinc, 1.0);

  // update end-point velocities and pressure
  UpdateIterIncrementally(state_->incvel_);

  // keep pressure values from previous time step
  state_->velpressplitter_->InsertCondVector(state_->velpressplitter_->ExtractCondVector(state_->veln_),state_->velnp_);

  // Note: accelerations on Dirichlet DOFs are not set.

  // reset to zero
  state_->incvel_->PutScalar(0.0);

  // free the system matrix to get the matrix deleted
  solver_->Reset();

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
//Overloaded in TimIntPoro and TimIntRedModels bk 12/13
void FLD::XFluid::UpdateIterIncrementally(
  Teuchos::RCP<const Epetra_Vector> vel
  )
{
  // set the new solution we just got
  if (vel != Teuchos::null)
  {
    // Take Dirichlet values from velnp and add vel to veln for non-Dirichlet
    // values.
    Teuchos::RCP<Epetra_Vector> aux = LINALG::CreateVector(
        *(discret_->DofRowMap(0)), true);
    aux->Update(1.0, *state_->velnp_, 1.0, *vel, 0.0);
    //    dbcmaps_->InsertOtherVector(dbcmaps_->ExtractOtherVector(aux), velnp_);
    state_->dbcmaps_->InsertCondVector(state_->dbcmaps_->ExtractCondVector(state_->velnp_), aux);

    *state_->velnp_ = *aux;
  }

  return;
}

// -------------------------------------------------------------------
// Read Restart data
// -------------------------------------------------------------------
void FLD::XFluid::ReadRestart(int step)
{

  //-------- fluid discretization
  IO::DiscretizationReader reader(discret_,step);
  time_ = reader.ReadDouble("time");
  step_ = reader.ReadInt("step");

  if(myrank_ == 0) IO::cout << "ReadRestart for fluid dis (time="<<time_<<" ; step="<<step_<<")" << IO::endl;

  if(myrank_ == 0) {IO::cout << RED << "Warning: For Restart we Cut the configuration of the last time step with the final (in best case converged)" <<
      " solution, without restart the configuration used would be one newton step ealier! --> This might lead to problems if the solution is no " <<
      " converged an therefore the dofset coming from restart and during simulation differ!" << END_COLOR << IO::endl;}

  if (alefluid_)
  {
    reader.ReadVector(dispnp_,"full_dispnp_res");
    reader.ReadVector(dispn_,"full_dispnp_res"); //as update() was called anyway before output...
    reader.ReadVector(gridvnp_,"full_gridvnp_res");
    reader.ReadVector(gridvn_,"full_gridvnp_res"); //as update() was called anyway before output...
    //state-vectors in state will be set in the creation of a new state
    CreateInitialState(); //Create an State with the deformed Fluid Mesh (otherwise state vectors wouldn't fit)
  }
  reader.ReadVector(state_->velnp_,"velnp_res");
  reader.ReadVector(state_->velnm_,"velnm_res");
  reader.ReadVector(state_->veln_, "veln_res" );
  reader.ReadVector(state_->accnp_,"accnp_res");
  reader.ReadVector(state_->accn_ ,"accn_res" );

#if(0)
  std::cout << "velnp_ " << *(state_->velnp_) << std::endl;
  std::cout << "veln_ "  << *(state_->veln_)  << std::endl;
  std::cout << "accnp_ " << *(state_->accnp_) << std::endl;
  std::cout << "accn_ "  << *(state_->accn_)  << std::endl;
#endif

  // set element time parameter after restart:
  // Here it is already needed by AVM3 and impedance boundary condition!!
  SetElementTimeParameter();

  // ensure that the overall dof numbering is identical to the one
  // that was used when the restart data was written. Especially
  // in case of multiphysics problems & periodic boundary conditions
  // it is better to check the consistency of the maps here:
  if (not (discret_->DofRowMap())->SameAs(state_->velnp_->Map()))
    dserror("Global dof numbering in maps does not match");
  if (not (discret_->DofRowMap())->SameAs(state_->veln_->Map()))
    dserror("Global dof numbering in maps does not match");
  if (not (discret_->DofRowMap())->SameAs(state_->accn_->Map()))
    dserror("Global dof numbering in maps does not match");


  // write gmsh-output for start fields
  // reference solution output
  output_service_->GmshSolutionOutputPrevious( "RESTART", step_, state_ );

  // set the new time and step also to the coupling objects
  condition_manager_->SetTimeAndStep(time_,step_);

}

// -------------------------------------------------------------------
// -------------------------------------------------------------------
Teuchos::RCP<XFEM::MeshCoupling> FLD::XFluid::GetMeshCoupling(const std::string & condname)
{
  return condition_manager_->GetMeshCoupling(condname);
}

// -------------------------------------------------------------------
// -------------------------------------------------------------------
Teuchos::RCP<LINALG::SparseMatrix> FLD::XFluid::C_sx_Matrix(const std::string & cond_name)
{
  const int coup_idx = condition_manager_->GetCouplingIndex(cond_name);
  return state_->coup_state_[coup_idx]->C_sx_;
}

// -------------------------------------------------------------------
// -------------------------------------------------------------------
Teuchos::RCP<LINALG::SparseMatrix> FLD::XFluid::C_xs_Matrix(const std::string & cond_name)
{
  const int coup_idx = condition_manager_->GetCouplingIndex(cond_name);
  return state_->coup_state_[coup_idx]->C_xs_;
}

// -------------------------------------------------------------------
// -------------------------------------------------------------------
Teuchos::RCP<LINALG::SparseMatrix> FLD::XFluid::C_ss_Matrix(const std::string & cond_name)
{
  const int coup_idx = condition_manager_->GetCouplingIndex(cond_name);
  return state_->coup_state_[coup_idx]->C_ss_;
}

// -------------------------------------------------------------------
// -------------------------------------------------------------------
Teuchos::RCP<Epetra_Vector> FLD::XFluid::RHS_s_Vec(const std::string & cond_name)
{
  const int coup_idx = condition_manager_->GetCouplingIndex(cond_name);
  return state_->coup_state_[coup_idx]->rhC_s_;
}

/*------------------------------------------------------------------------------------------------*
 | create field test
 *------------------------------------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> FLD::XFluid::CreateFieldTest()
{
  return Teuchos::rcp(new FLD::XFluidResultTest(*this));
}


void FLD::XFluid::GenAlphaIntermediateValues()
{
  //       n+alphaM                n+1                      n
  //    acc         = alpha_M * acc     + (1-alpha_M) *  acc
  //       (i)                     (i)
  {
    // extract the degrees of freedom associated with velocities
    // only these are allowed to be updated, otherwise you will
    // run into trouble in loma, where the 'pressure' component
    // is used to store the acceleration of the temperature
    Teuchos::RCP<Epetra_Vector> onlyaccn  = state_->velpressplitter_->ExtractOtherVector(state_->accn_ );
    Teuchos::RCP<Epetra_Vector> onlyaccnp = state_->velpressplitter_->ExtractOtherVector(state_->accnp_);

    Teuchos::RCP<Epetra_Vector> onlyaccam = Teuchos::rcp(new Epetra_Vector(onlyaccnp->Map()));

    onlyaccam->Update((alphaM_),*onlyaccnp,(1.0-alphaM_),*onlyaccn,0.0);

    // copy back into global vector
    LINALG::Export(*onlyaccam,*state_->accam_);
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
  state_->velaf_->Update((alphaF_),*state_->velnp_,(1.0-alphaF_),*state_->veln_,0.0);
}

void FLD::XFluid::GenAlphaUpdateAcceleration()
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
  Teuchos::RCP<Epetra_Vector> onlyaccn  = state_->velpressplitter_->ExtractOtherVector(state_->accn_ );
  Teuchos::RCP<Epetra_Vector> onlyveln  = state_->velpressplitter_->ExtractOtherVector(state_->veln_ );
  Teuchos::RCP<Epetra_Vector> onlyvelnp = state_->velpressplitter_->ExtractOtherVector(state_->velnp_);

  Teuchos::RCP<Epetra_Vector> onlyaccnp = Teuchos::rcp(new Epetra_Vector(onlyaccn->Map()));

  const double fact1 = 1.0/(gamma_*dta_);
  const double fact2 = 1.0 - (1.0/gamma_);
  onlyaccnp->Update(fact2,*onlyaccn,0.0);
  onlyaccnp->Update(fact1,*onlyvelnp,-fact1,*onlyveln,1.0);

  // copy back into global vector
  LINALG::Export(*onlyaccnp,*state_->accnp_);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::XFluid::UpdateGridv()
{
  // get order of accuracy of grid velocity determination
  // from input file data
  const Teuchos::ParameterList& fluiddynparams =  DRT::Problem::Instance()->FluidDynamicParams();
  const int order = DRT::INPUT::IntegralValue<INPAR::FLUID::Gridvel>(fluiddynparams, "GRIDVEL");

  Teuchos::RCP<Epetra_Vector> gridv = Teuchos::rcp(new Epetra_Vector(dispnp_->Map(),true));

  switch (order)
  {
    case INPAR::FLUID::BE:
      /* get gridvelocity from BE time discretisation of mesh motion:
           -> cheap
           -> easy
           -> limits FSI algorithm to first order accuracy in time

                  x^n+1 - x^n
             uG = -----------
                    Delta t                        */
      gridvnp_->Update(1/dta_, *dispnp_, -1/dta_, *dispn_, 0.0);
    break;
    case INPAR::FLUID::BDF2:
      /* get gridvelocity from BDF2 time discretisation of mesh motion:
           -> requires one more previous mesh position or displacement
           -> somewhat more complicated
           -> allows second order accuracy for the overall flow solution  */
      gridvnp_->Update(1.5/dta_, *dispnp_, -2.0/dta_, *dispn_, 0.0);
      gridvnp_->Update(0.5/dta_, *dispnm_, 1.0);
    break;
    case INPAR::FLUID::OST:
    {
      /* get gridvelocity from OST time discretisation of mesh motion:
         -> needed to allow consistent linearization of FPSI problem  */
      const double theta = fluiddynparams.get<double>("THETA");
      gridvnp_->Update(1/(theta*dta_), *dispnp_, -1/(theta*dta_), *dispn_, 0.0);
      gridvnp_->Update(-((1.0/theta)-1.0),*gridvn_,1.0);
    }
    break;
    default:
      dserror("Unknown or invalid type of grid velocity determination. Fix GRIDVEL section of your input file.");
    break;
  }

}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::XFluid::UpdateByIncrement()
{
  state_->Velnp()->Update(1.0,*state_->IncVel(),1.0);
  double f_norm = 0;
  state_->Velnp()->Norm2(&f_norm);
//  std::cout << std::setprecision(14) << f_norm << std::endl;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::XFluid::SetOldPartOfRighthandside()
{
  SetOldPartOfRighthandside(
      state_->veln_,
      state_->velnm_,
      state_->accn_,
      timealgo_,
      dta_,
      theta_,
      state_->hist_);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::XFluid::SetOldPartOfRighthandside(
    const Teuchos::RCP<Epetra_Vector>&   veln,
    const Teuchos::RCP<Epetra_Vector>&   velnm,
    const Teuchos::RCP<Epetra_Vector>&   accn,
    const INPAR::FLUID::TimeIntegrationScheme timealgo,
    const double                         dta,
    const double                         theta,
    Teuchos::RCP<Epetra_Vector>&         hist
    )
{
  /*!
    \brief Set the part of the righthandside belonging to the last
           timestep for incompressible or low-Mach-number flow

       for low-Mach-number flow: distinguish momentum and continuity part
       (continuity part only meaningful for low-Mach-number flow)

       Stationary/af-generalized-alpha:

                     mom: hist_ = 0.0
                    (con: hist_ = 0.0)

       One-step-Theta:

                     mom: hist_ = veln_  + dt*(1-Theta)*accn_
                    (con: hist_ = densn_ + dt*(1-Theta)*densdtn_)

       BDF2: for constant time step:

                     mom: hist_ = 4/3 veln_  - 1/3 velnm_
                    (con: hist_ = 4/3 densn_ - 1/3 densnm_)


   */
  switch (timealgo)
  {
  case INPAR::FLUID::timeint_stationary: /* Stationary algorithm */
  case INPAR::FLUID::timeint_afgenalpha: /* Af-generalized-alpha time integration */
  case INPAR::FLUID::timeint_npgenalpha:
    hist->PutScalar(0.0);
    break;

  case INPAR::FLUID::timeint_one_step_theta: /* One step Theta time integration */
    hist->Update(1.0, *veln, dta*(1.0-theta), *accn, 0.0);
    break;

  case INPAR::FLUID::timeint_bdf2:    /* 2nd order backward differencing BDF2 */
    hist->Update(4./3., *veln, -1./3., *velnm, 0.0);
    break;

  default:
  {
    dserror("Time integration scheme unknown!");
    break;
  }
  }
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::XFluid::SetGamma(Teuchos::ParameterList& eleparams)
{
  if (timealgo_==INPAR::FLUID::timeint_afgenalpha)
  {
    eleparams.set("gamma"  ,gamma_);
  }
  else if (timealgo_==INPAR::FLUID::timeint_one_step_theta)
  {
    eleparams.set("gamma"  ,theta_);
  }
  else if (timealgo_==INPAR::FLUID::timeint_bdf2)
  {
    eleparams.set("gamma"  ,1.0);
  }
  else dserror("unknown timealgo_");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::XFluid::SetStateTimInt()
{
  // set scheme-specific element parameters and vector values
  if (timealgo_==INPAR::FLUID::timeint_afgenalpha)
    discret_->SetState("velaf",state_->velaf_);
  else
    discret_->SetState("velaf",state_->velnp_);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::XFluid::CalculateAcceleration(
    const Teuchos::RCP<const Epetra_Vector>    velnp,
    const Teuchos::RCP<const Epetra_Vector>    veln,
    const Teuchos::RCP<const Epetra_Vector>    velnm,
    const Teuchos::RCP<const Epetra_Vector>    accn,
    const Teuchos::RCP<Epetra_Vector>          accnp
)
{
    /*

    Following formulations are for n+1; acceleration values, however, are
    directly stored in vectors at time n (velocity has not yet been updated).

    One-step-Theta:

     acc(n+1) = (vel(n+1)-vel(n)) / (Theta * dt(n)) - (1/Theta -1) * acc(n)


    BDF2:

                   2*dt(n)+dt(n-1)                  dt(n)+dt(n-1)
     acc(n+1) = --------------------- vel(n+1) - --------------- vel(n)
                 dt(n)*[dt(n)+dt(n-1)]              dt(n)*dt(n-1)

                         dt(n)
               + ----------------------- vel(n-1)
                 dt(n-1)*[dt(n)+dt(n-1)]

    */

    switch (timealgo_)
    {
      case INPAR::FLUID::timeint_stationary: /* no accelerations for stationary problems*/
      {
        accnp->PutScalar(0.0);
        break;
      }
      case INPAR::FLUID::timeint_one_step_theta: /* One-step-theta time integration */
      {
        const double fact1 = 1.0/(theta_*dta_);
        const double fact2 =-1.0/theta_ +1.0;   /* = -1/Theta + 1 */

        accnp->Update( fact1,*velnp,0.0);
        accnp->Update(-fact1,*veln ,1.0);
        accnp->Update( fact2,*accn,1.0);
        break;
      }
      case INPAR::FLUID::timeint_bdf2:    /* 2nd order backward differencing BDF2 */
      {
        //TODO: computed, even though not really used afterwards! CHECK!!!
        if (dta_*dtp_ < EPS15) dserror("Zero time step size!!!!!");
        const double sum = dta_ + dtp_;

        accnp->Update((2.0*dta_+dtp_)/(dta_*sum),*velnp, -sum/(dta_*dtp_),*veln ,0.0);
        accnp->Update(dta_/(dtp_*sum),*velnm,1.0);
        break;
      }
      case INPAR::FLUID::timeint_afgenalpha: /* Af-generalized-alpha time integration */
      case INPAR::FLUID::timeint_npgenalpha:
      {
        // do nothing: new acceleration is calculated at beginning of next time step
        break;
      }
      default:
      {
        dserror("Time integration scheme unknown!");
        break;
      }
    }

  return;
}
