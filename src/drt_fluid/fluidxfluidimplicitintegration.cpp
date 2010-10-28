/*----------------------------------------------------------------------*/
/*!
\file fluidxfluidimplicitintegration.cpp


<pre>
Maintainer: Shadan Shahmiri
            shahmiri@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET


#undef WRITEOUTSTATISTICS

#include "fluidxfluidimplicitintegration.H"
#include "time_integration_scheme.H"
#include "../linalg/linalg_ana.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_mapextractor.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_condition_utils.H"
#include "fluid_utils.H"
#include "fluid_utils_mapextractor.H"

#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_xfem/dof_management.H"
#include "../drt_xfem/dof_distribution_switcher.H"
#include "../drt_xfem/enrichment_utils.H"
#include "../drt_xfem/element_ansatz.H"
#include "../drt_geometry/position_array.H"
#include "../drt_f3/xfluid3_interpolation.H"
#include "../drt_xdiff3/xdiff3_interpolation.H"
#include "../drt_xdiff3/xdiff3_interpolation.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_io/io_control.H"
#include "../drt_inpar/inpar_xfem.H"
#include "../drt_inpar/inpar_fluid.H"
#include "../drt_lib/drt_dofset_transparent.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "../drt_inpar/inpar_xfem.H"
#include "../drt_lib/drt_elementtype.H"
#include "../drt_xdiff3/xdiff3.H"
FLD::FluidXFluidImplicitTimeInt::FluidXFluidImplicitTimeInt(RefCountPtr<DRT::Discretization> fluiddis,
                                                            RefCountPtr<DRT::Discretization> xfluiddis,
                                                            LINALG::Solver&                  solver,
                                                            ParameterList&                   params,
                                                            IO::DiscretizationWriter&        output,
                                                            bool                             alefluid)
  :
  // call constructor for "nontrivial" objects
  fluiddis_(fluiddis),
  xfluiddis_(xfluiddis),
  solver_ (solver),
  params_ (params),
  xparams_(params.sublist("XFEM")),
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
  inrelaxation_(false)
{

  // -------------------------------------------------------------------
  // get the processor ID from the communicator
  // -------------------------------------------------------------------
  myrank_  = fluiddis_->Comm().MyPID();

  // time measurement: initialization
  TEUCHOS_FUNC_TIME_MONITOR(" + initialization");

  cout << "FLD::FluidXFluidImplicitTimeInt::FluidXFluidImplicitTimeInt()" << endl;
  
  // -------------------------------------------------------------------
  // get the basic parameters first
  // -------------------------------------------------------------------

  // physical type of fluid flow (incompressible, varying density, loma, Boussinesq approximation)
  physicaltype_ = params_.get<INPAR::FLUID::PhysicalType>("Physical Type");
  // type of time-integration
  timealgo_ = params_.get<INPAR::FLUID::TimeIntegrationScheme>("time int algo");
  // time-step size
  dtp_ = dta_ = params_.get<double>("time step size");
  // maximum number of timesteps
  stepmax_  = params_.get<int>   ("max number timesteps");
  // maximum simulation time
  maxtime_  = params_.get<double>("total time");
  // parameter theta for time-integration schemes
  theta_    = params_.get<double>("theta");
  
  if (timealgo_ == INPAR::FLUID::timeint_stationary and params_.get<double>("time step size") != 1.0)
    dserror("Timestep size (delta t) has to be 1.0 for stationary computations!");
  
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
  newton_ = params_.get<INPAR::FLUID::LinearisationAction>("Linearisation");

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
  pbcmapmastertoslave_ = params_.get<RCP<map<int,vector<int> > > >("periodic bc");
  fluiddis_->ComputeNullSpaceIfNecessary(solver_.Params(),true);

  // ensure that degrees of freedom in the discretization have been set
  if (!fluiddis_->Filled() || !fluiddis_->HaveDofs()) fluiddis_->FillComplete();

  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  // -------------------------------------------------------------------
  const Epetra_Map* fluiddofrowmap = fluiddis_->DofRowMap();

  // -------------------------------------------------------------------
  // get a vector layout from the discretization for a vector which only
  // contains the velocity dofs and for one vector which only contains
  // pressure degrees of freedom.
  // -------------------------------------------------------------------
  numdim_ = params_.get<int>("number of velocity degrees of freedom");

  FLD::UTILS::SetupFluidSplit(*fluiddis_,numdim_,velpressplitter_);

  // -------------------------------------------------------------------
  // create empty system matrix --- stiffness and mass are assembled in
  // one system matrix!
  // -------------------------------------------------------------------

  // This is a first estimate for the number of non zeros in a row of
  // the matrix. Assuming a structured 3d-fluid mesh we have 27 adjacent
  // nodes with 4 dofs each. (27*4=108)
  // We do not need the exact number here, just for performance reasons
  // a 'good' estimate

  if (not params_.get<int>("Simple Preconditioner",0) && not params_.get<int>("AMG BS Preconditioner",0))
  {
    // initialize standard (stabilized) system matrix
    fluidsysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*fluiddofrowmap,108,false,true));
  }
  else
  {
    Teuchos::RCP<LINALG::BlockSparseMatrix<FLD::UTILS::VelPressSplitStrategy> > blocksysmat =
    Teuchos::rcp(new LINALG::BlockSparseMatrix<FLD::UTILS::VelPressSplitStrategy>(velpressplitter_,velpressplitter_,108,false,true));
    blocksysmat->SetNumdim(numdim_);
    fluidsysmat_ = blocksysmat;
  }

  // sysmat might be singular (if we have a purely Dirichlet constrained
  // problem, the pressure mode is defined only up to a constant)
  // in this case, we need a basis vector for the nullspace/kernel
  vector<DRT::Condition*> KSPcond;
  fluiddis_->GetCondition("KrylovSpaceProjection",KSPcond);
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
    w_       = LINALG::CreateVector(*fluiddofrowmap,true);
    c_       = LINALG::CreateVector(*fluiddofrowmap,true);
    kspsplitter_.Setup(*fluiddis_);
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

  // velocity/pressure at time n+1, n and n-1
  fluidstate_.velnp_ = LINALG::CreateVector(*fluiddofrowmap,true);
  fluidstate_.veln_  = LINALG::CreateVector(*fluiddofrowmap,true);
  fluidstate_.velnm_ = LINALG::CreateVector(*fluiddofrowmap,true);

  // acceleration/(scalar time derivative) at time n+1 and n
  fluidstate_.accnp_ = LINALG::CreateVector(*fluiddofrowmap,true);
  fluidstate_.accn_  = LINALG::CreateVector(*fluiddofrowmap,true);

  // velocity/pressure at time n+alpha_F
  fluidstate_.velaf_ = LINALG::CreateVector(*fluiddofrowmap,true);

  // acceleration/(scalar time derivative) at time n+alpha_M/(n+alpha_M/n)
  fluidstate_.accam_ = LINALG::CreateVector(*fluiddofrowmap,true);

  // scalar at time n+alpha_F/n+1 and n+alpha_M/n
  // (only required for low-Mach-number case)
  fluidstate_.scaaf_ = LINALG::CreateVector(*fluiddofrowmap,true);
  fluidstate_.scaam_ = LINALG::CreateVector(*fluiddofrowmap,true);

  // history vector
  hist_ = LINALG::CreateVector(*fluiddofrowmap,true);

  if (alefluid_)
  {
    fluidstate_.dispnp_ = LINALG::CreateVector(*fluiddofrowmap,true);
    fluidstate_.dispn_  = LINALG::CreateVector(*fluiddofrowmap,true);
    fluidstate_.dispnm_ = LINALG::CreateVector(*fluiddofrowmap,true);
    gridv_  = LINALG::CreateVector(*fluiddofrowmap,true);
  }

  // Vectors associated to boundary conditions
  // -----------------------------------------

  // a vector of zeros to be used to enforce zero dirichlet boundary conditions
  zeros_   = LINALG::CreateVector(*fluiddofrowmap,true);

  // object holds maps/subsets for DOFs subjected to Dirichlet BCs and otherwise
  dbcmaps_ = Teuchos::rcp(new LINALG::MapExtractor());
  {
    ParameterList eleparams;
    // other parameters needed by the elements
    eleparams.set("total time",time_);
    fluiddis_->EvaluateDirichlet(eleparams, zeros_, Teuchos::null, Teuchos::null,
                                Teuchos::null, dbcmaps_);

    zeros_->PutScalar(0.0); // just in case of change
  }

  // the vector containing body and surface forces
  neumann_loads_= LINALG::CreateVector(*fluiddofrowmap,true);

  // Vectors used for solution process
  // ---------------------------------

  // rhs: standard (stabilized) residual vector (rhs for the incremental form)
  //residual_     = LINALG::CreateVector(*dofrowmap,true);
  fluidresidual_     = LINALG::CreateVector(*fluiddofrowmap,true);
  trueresidual_ = LINALG::CreateVector(*fluiddofrowmap,true);

  // Nonlinear iteration increment vector
  incvel_ = LINALG::CreateVector(*fluiddofrowmap,true);

  // get constant density variable for incompressible flow
  ParameterList eleparams;
  eleparams.set("action","get_density");
  eleparams.set("Physical Type", physicaltype_);
  fluiddis_->Evaluate(eleparams,null,null,null,null,null);
  density_ = eleparams.get("density", 1.0);
  if (density_ < EPS15) dserror("received zero or negative density value");


  // initialize all thermodynamic pressure values and its time derivative
  // to one or zero, respectively
  // -> they are kept this way for incompressible flow
  thermpressaf_   = 1.0;
  thermpressam_   = 1.0;
  thermpressdtam_ = 0.0;
  
  //XFEM stuff
  {
    ParameterList eleparams;
    eleparams.set("action","set_output_mode");
    eleparams.set("output_mode",true);
    xfluiddis_->Evaluate(eleparams);
  }

  xfluiddis_->FillComplete();
  
  //output_.WriteMesh(0,0.0);
  
  if (xfluiddis_->Comm().MyPID()==0)
  {
    if (not xparams_.get<bool>("DLM_condensation"))
    {
      std::cout << RED_LIGHT << "DLM_condensation turned off!" << END_COLOR << endl;
    }
  }
  
  // store a dofset with the complete fluid unknowns
  xfluiddofset_out_.Reset();
  xfluiddofset_out_.AssignDegreesOfFreedom(*xfluiddis_,0,0);
  // split based on complete fluid field
  FLD::UTILS::SetupFluidSplit(*xfluiddis_,xfluiddofset_out_,3,xfluidvelpressplitterForOutput_);
  
  xfluidstate_.nodalDofDistributionMap_.clear();
  xfluidstate_.elementalDofDistributionMap_.clear();
  
  {
    ParameterList eleparams;
    eleparams.set("action","set_output_mode");
    eleparams.set("output_mode",false);
    xfluiddis_->Evaluate(eleparams);
  }
    
  // print information about elements
  // (to double-check and log that correct input has been read)
  std::set<DRT::Element::DiscretizationType> distypeset;
  std::set<DRT::ElementType*> etypeset;
  for (int i=0; i<xfluiddis_->NumMyColElements(); ++i)
  {
    distypeset.insert(xfluiddis_->lColElement(i)->Shape());
    etypeset.insert(&xfluiddis_->lColElement(i)->ElementType());
  }
  
  physprob_.fieldset_.insert(XFEM::PHYSICS::Velx);
  physprob_.fieldset_.insert(XFEM::PHYSICS::Vely);
  physprob_.fieldset_.insert(XFEM::PHYSICS::Velz);
  physprob_.fieldset_.insert(XFEM::PHYSICS::Pres);
  switch (xparams_.get<INPAR::XFEM::BoundaryIntegralType>("EMBEDDED_BOUNDARY"))
  {
  case INPAR::XFEM::BoundaryTypeSigma:
    physprob_.elementAnsatz_  = rcp<XFLUID::FluidElementAnsatz>(new XFLUID::FluidElementAnsatz());
    break;
  case INPAR::XFEM::BoundaryTypeTauPressure:
    physprob_.elementAnsatz_  = rcp<XFLUID::FluidElementAnsatzWithExtraElementPressure>(new XFLUID::FluidElementAnsatzWithExtraElementPressure());
    break;
  default:
    dserror("unknown boundary type");
  }
  
  {
    cout << "Element shapes in xfluid discretization: ";
    xfluiddis_->Comm().Barrier();
    bool moreThanOne = false;
    for (std::set<DRT::Element::DiscretizationType>::const_iterator iter = distypeset.begin(); iter != distypeset.end(); ++iter)
    {
      if (moreThanOne)  cout << ", ";
      cout << DRT::DistypeToString(*iter);
      moreThanOne = true;
    }
    if (xfluiddis_->Comm().MyPID()==0)
    {
      cout << endl;
    }
    xfluiddis_->Comm().Barrier();
  }
  
  {
    cout << "Element types in xfluid discretization: ";
    xfluiddis_->Comm().Barrier();
    bool moreThanOnee = false;
    for (std::set<DRT::ElementType*>::const_iterator iter = etypeset.begin(); iter != etypeset.end(); ++iter)
    {
      if (moreThanOnee)  cout << ", ";
      cout << ( *iter )->Name();
      moreThanOnee = true;
    }
    if (xfluiddis_->Comm().MyPID()==0)
    {
      cout << endl;
    }
    xfluiddis_->Comm().Barrier();
  }  
} // FluidXFluidImplicitTimeInt::FluidXFluidImplicitTimeInt


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
void FLD::FluidXFluidImplicitTimeInt::Integrate(
    const Teuchos::RCP<DRT::Discretization>    fluidxfluidboundarydis
    )
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
  if (timealgo_==INPAR::FLUID::timeint_stationary) SolveStationaryProblem(fluidxfluidboundarydis);
  else TimeLoop(fluidxfluidboundarydis);

  // print the results of time measurements
  TimeMonitor::summarize();

  return;
} // FluidXFluidImplicitTimeInt::Integrate


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | contains the time loop                                    gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidXFluidImplicitTimeInt::TimeLoop(
    const Teuchos::RCP<DRT::Discretization> fluidxfluidboundarydis
        )
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

  const Epetra_Map* fluidxfluid_dofcolmap = fluidxfluidboundarydis->DofColMap();

  fluidxfluidboundarydis->SetState("idispcolnp", LINALG::CreateVector(*fluidxfluid_dofcolmap,true));
  fluidxfluidboundarydis->SetState("ivelcolnp", LINALG::CreateVector(*fluidxfluid_dofcolmap,true));
  fluidxfluidboundarydis->SetState("idispcoln", LINALG::CreateVector(*fluidxfluid_dofcolmap,true));
  fluidxfluidboundarydis->SetState("ivelcoln", LINALG::CreateVector(*fluidxfluid_dofcolmap,true));
  fluidxfluidboundarydis->SetState("ivelcolnm", LINALG::CreateVector(*fluidxfluid_dofcolmap,true));
  fluidxfluidboundarydis->SetState("iacccoln", LINALG::CreateVector(*fluidxfluid_dofcolmap,true));


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
        printf("TIME: %11.4E/%11.4E  DT = %11.4E   One-Step-Theta    STEP = %4d/%4d \n",
              time_,maxtime_,dta_,step_,stepmax_);
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
    NonlinearSolve(fluidxfluidboundarydis);

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
void FLD::FluidXFluidImplicitTimeInt::PrepareTimeStep()
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
  TIMEINT_THETA_BDF2::SetOldPartOfRighthandside(fluidstate_.veln_,fluidstate_.velnm_, fluidstate_.accn_,
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
        fluidstate_.veln_,
        fluidstate_.velnm_,
        fluidstate_.accn_,
        velpressplitter_,
        timealgo_,
        theta_,
        dta_,
        dtp_,
        fluidstate_.velnp_,
        fluiddis_->Comm());
    }
  }

  // -------------------------------------------------------------------
  //  evaluate Dirichlet and Neumann boundary conditions
  // -------------------------------------------------------------------
  {
    ParameterList eleparams;

    // total time required for Dirichlet conditions
    eleparams.set("total time",time_);

    // set vector values needed by elements
    fluiddis_->ClearState();
    fluiddis_->SetState("velnp",fluidstate_.velnp_);

    // predicted dirichlet values
    // velnp then also holds prescribed new dirichlet values
    fluiddis_->EvaluateDirichlet(eleparams,fluidstate_.velnp_,null,null,null);

    fluiddis_->ClearState();

    // set all parameters and states required for Neumann conditions
 //   eleparams.set("Physical Type",physicaltype_);
    eleparams.set("thermpress at n+1",thermpressaf_);
    if (timealgo_==INPAR::FLUID::timeint_afgenalpha)
    {
      eleparams.set("total time",time_-(1-alphaF_)*dta_);
      eleparams.set("thsl",1.0);
    }
    else
    {
      eleparams.set("total time",time_);
      eleparams.set("thsl",theta_*dta_);
    }

    // evaluate Neumann conditions
    neumann_loads_->PutScalar(0.0);
    fluiddis_->SetState("scanp",fluidstate_.scaaf_);
    fluiddis_->EvaluateNeumann(eleparams,*neumann_loads_);
    fluiddis_->ClearState();
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
    GenAlphaUpdateAcceleration();

    // ----------------------------------------------------------------
    // compute values at intermediate time steps
    // ----------------------------------------------------------------
    GenAlphaIntermediateValues();
  }

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
void FLD::FluidXFluidImplicitTimeInt::NonlinearSolve(
    const Teuchos::RCP<DRT::Discretization> fluidxfluidboundarydis
    )
{
  cout << "FLD::FluidXFluidImplicitTimeInt::NonlinearSolve()" << endl;
  
  inrelaxation_ = false;
  dirichletlines_ = Teuchos::null;
  // Do not remove meshmatrix_ here as we want to reuse its graph.
  // (We pay for the memory anyway if we use it, we might as well keep it.)
  //meshmatrix_ = Teuchos::null;

  // time measurement: nonlinear iteration
  TEUCHOS_FUNC_TIME_MONITOR("   + nonlin. iteration/lin. solve");
  
  {
    ParameterList eleparams;
    eleparams.set("action","reset");
    xfluiddis_->Evaluate(eleparams);
  }

  ComputeInterfaceAndSetDOFs(fluidxfluidboundarydis);

  // ---------------------------------------------- nonlinear iteration
  // ------------------------------- stop nonlinear iteration when both
  //                                 increment-norms are below this bound
  const double  ittol     =params_.get<double>("tolerance for nonlin iter");

  //------------------------------ turn adaptive solver tolerance on/off
  const bool   isadapttol    = params_.get<bool>("ADAPTCONV",true);
  const double adaptolbetter = params_.get<double>("ADAPTCONV_BETTER",0.01);


  int  itnum = 0;
  int  itemax = 0;
  bool stopnonliniter = false;

  // currently default for turbulent channel flow: only one iteration before sampling
  itemax  = params_.get<int>   ("max nonlin iter steps");

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

      fluidsysmat_->Zero();

      // create the parameters for the discretization
      ParameterList eleparams;

      // add Neumann loads
      residual_->Update(1.0,*neumann_loads_,0.0);

      // set general element parameters
      eleparams.set("dt",dta_);
      eleparams.set("theta",theta_);
      eleparams.set("omtheta",omtheta_);
      eleparams.set("form of convective term",convform_);
      eleparams.set("fs subgrid viscosity",fssgv_);
      eleparams.set("Linearisation",newton_);
      eleparams.set("Physical Type", physicaltype_);

      // parameters for stabilization
      eleparams.sublist("STABILIZATION") = params_.sublist("STABILIZATION");

      // parameters for stabilization
      eleparams.sublist("TURBULENCE MODEL") = params_.sublist("TURBULENCE MODEL");

      eleparams.set("thermpress at n+alpha_F/n+1",thermpressaf_);
      eleparams.set("thermpress at n+alpha_M/n",thermpressam_);
      eleparams.set("thermpressderiv at n+alpha_M/n+1",thermpressdtam_);

      // set general vector values needed by elements
      fluiddis_->ClearState();
      fluiddis_->SetState("hist" ,hist_ );
      fluiddis_->SetState("accam",fluidstate_.accam_);
      fluiddis_->SetState("scaaf",fluidstate_.scaaf_);
      fluiddis_->SetState("scaam",fluidstate_.scaam_);
      if (alefluid_)
      {
        fluiddis_->SetState("dispnp", fluidstate_.dispnp_);
        fluiddis_->SetState("gridv", gridv_);
      }

      // set scheme-specific element parameters and vector values
      if (timealgo_==INPAR::FLUID::timeint_stationary)
      {
        eleparams.set("action","calc_fluid_stationary_systemmat_and_residual");
        eleparams.set("using generalized-alpha time integration",false);
        eleparams.set("total time",time_);
        eleparams.set("is stationary", true);

        fluiddis_->SetState("velaf",fluidstate_.velnp_);
      }
      else if (timealgo_==INPAR::FLUID::timeint_afgenalpha)
      {
        eleparams.set("action","calc_fluid_afgenalpha_systemmat_and_residual");
        eleparams.set("using generalized-alpha time integration",true);
        eleparams.set("total time",time_-(1-alphaF_)*dta_);
        eleparams.set("is stationary", false);
        eleparams.set("alphaF",alphaF_);
        eleparams.set("alphaM",alphaM_);
        eleparams.set("gamma",gamma_);

        fluiddis_->SetState("velaf",fluidstate_.velaf_);
      }
      else
      {
        eleparams.set("action","calc_fluid_systemmat_and_residual");
        eleparams.set("using generalized-alpha time integration",false);
        eleparams.set("total time",time_);
        eleparams.set("is stationary", false);

        fluiddis_->SetState("velaf",fluidstate_.velnp_);
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
        fluiddis_->Evaluate(eleparams,fluidsysmat_,null,residual_,null,null);

        fluiddis_->ClearState();

        // account for potential Neumann inflow terms
        if (neumanninflow_)
        {
          // create parameter list
          ParameterList condparams;

          // action for elements
          condparams.set("action","calc_Neumann_inflow");
          condparams.set("thsl",theta_*dta_);
      //    condparams.set("Physical Type",physicaltype_);
          condparams.set("thermpress at n+alpha_F/n+1",thermpressaf_);

          // set vector values needed by elements
          fluiddis_->ClearState();
          fluiddis_->SetState("scaaf",fluidstate_.scaaf_);
          if (timealgo_==INPAR::FLUID::timeint_afgenalpha)
          {
            condparams.set("using generalized-alpha time integration",true);
            fluiddis_->SetState("velaf",fluidstate_.velaf_);
          }
          else
          {
            condparams.set("using generalized-alpha time integration",false);
            fluiddis_->SetState("velaf",fluidstate_.velnp_);
          }

          std::string condstring("FluidNeumannInflow");
          fluiddis_->EvaluateCondition(condparams,fluidsysmat_,Teuchos::null,residual_,Teuchos::null,Teuchos::null,condstring);
          fluiddis_->ClearState();
        }

        if (timealgo_==INPAR::FLUID::timeint_afgenalpha)
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
        fluidsysmat_->Complete();
      }

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

    velpressplitter_.ExtractOtherVector(fluidstate_.velnp_,onlyvel);
    onlyvel->Norm2(&velnorm_L2);

    Teuchos::RCP<Epetra_Vector> onlypre = velpressplitter_.ExtractCondVector(residual_);
    onlypre->Norm2(&presnorm);

    velpressplitter_.ExtractCondVector(incvel_,onlypre);
    onlypre->Norm2(&incprenorm_L2);

    velpressplitter_.ExtractCondVector(fluidstate_.velnp_,onlypre);
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
      LINALG::ApplyDirichlettoSystem(fluidsysmat_,incvel_,residual_,zeros_,*(dbcmaps_->CondMap()));
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

      solver_.Solve(fluidsysmat_->EpetraOperator(),incvel_,residual_,true,itnum==1, w_, c_, project_);
      solver_.ResetTolerance();

      // end time measurement for solver
      dtsolve_ = Teuchos::Time::wallTime()-tcpusolve;
    }

    // -------------------------------------------------------------------
    // update velocity and pressure values by increments
    // -------------------------------------------------------------------
    fluidstate_.velnp_->Update(1.0,*incvel_,1.0);

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
} // FluidImplicitTimeInt::NonlinearSolve

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | compute values at intermediate time steps for gen.-alpha    vg 02/09 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidXFluidImplicitTimeInt::GenAlphaIntermediateValues()
{
  //       n+alphaM                n+1                      n
  //    acc         = alpha_M * acc     + (1-alpha_M) *  acc
  //       (i)                     (i)
  {
    // extract the degrees of freedom associated with velocities
    // only these are allowed to be updated, otherwise you will
    // run into trouble in loma, where the 'pressure' component
    // is used to store the acceleration of the temperature
    Teuchos::RCP<Epetra_Vector> onlyaccn  = velpressplitter_.ExtractOtherVector(fluidstate_.accn_ );
    Teuchos::RCP<Epetra_Vector> onlyaccnp = velpressplitter_.ExtractOtherVector(fluidstate_.accnp_);

    Teuchos::RCP<Epetra_Vector> onlyaccam = rcp(new Epetra_Vector(onlyaccnp->Map()));

    onlyaccam->Update((alphaM_),*onlyaccnp,(1.0-alphaM_),*onlyaccn,0.0);

    // copy back into global vector
    LINALG::Export(*onlyaccam,*fluidstate_.accam_);
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
  fluidstate_.velaf_->Update((alphaF_),*fluidstate_.velnp_,(1.0-alphaF_),*fluidstate_.veln_,0.0);

} // FluidXFluidImplicitTimeInt::GenAlphaIntermediateValues


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | call elements to calculate system matrix/rhs and assemble   vg 02/09 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidXFluidImplicitTimeInt::AssembleMatAndRHS()
{
  dtele_    = 0.0;
  dtfilter_ = 0.0;

  // time measurement: element
  TEUCHOS_FUNC_TIME_MONITOR("      + element calls");

  // get cpu time
  const double tcpu=Teuchos::Time::wallTime();

  fluidsysmat_->Zero();

  // create the parameters for the discretization
  ParameterList eleparams;

  // add Neumann loads
  residual_->Update(1.0,*neumann_loads_,0.0);
  
  // set general element parameters
  eleparams.set("dt",dta_);
  eleparams.set("theta",theta_);
  eleparams.set("omtheta",omtheta_);
  eleparams.set("form of convective term",convform_);
  eleparams.set("fs subgrid viscosity",fssgv_);
  eleparams.set("Linearisation",newton_);
 // eleparams.set("Physical Type", physicaltype_);

  // parameters for stabilization
  eleparams.sublist("STABILIZATION") = params_.sublist("STABILIZATION");

  // parameters for turbulence model
  eleparams.sublist("TURBULENCE MODEL") = params_.sublist("TURBULENCE MODEL");

  eleparams.set("thermpress at n+alpha_F/n+1",thermpressaf_);
  eleparams.set("thermpress at n+alpha_M/n",thermpressam_);
  eleparams.set("thermpressderiv at n+alpha_M/n+1",thermpressdtam_);

  // set general vector values needed by elements
  fluiddis_->ClearState();
  fluiddis_->SetState("hist" ,hist_ );
  fluiddis_->SetState("accam",fluidstate_.accam_);
  fluiddis_->SetState("scaaf",fluidstate_.scaaf_);
  fluiddis_->SetState("scaam",fluidstate_.scaam_);
  if (alefluid_)
  {
    fluiddis_->SetState("dispnp", fluidstate_.dispnp_);
    fluiddis_->SetState("gridv", gridv_);
  }

  // set scheme-specific element parameters and vector values
  if (timealgo_==INPAR::FLUID::timeint_stationary)
  {
    eleparams.set("action","calc_fluid_stationary_systemmat_and_residual");
    eleparams.set("using generalized-alpha time integration",false);
    eleparams.set("total time",time_);
    eleparams.set("is stationary", true);

    fluiddis_->SetState("velaf",fluidstate_.velnp_);
  }
  else if (timealgo_==INPAR::FLUID::timeint_afgenalpha)
  {
    eleparams.set("action","calc_fluid_afgenalpha_systemmat_and_residual");
    eleparams.set("using generalized-alpha time integration",true);
    eleparams.set("total time",time_-(1-alphaF_)*dta_);
    eleparams.set("is stationary", false);
    eleparams.set("alphaF",alphaF_);
    eleparams.set("alphaM",alphaM_);
    eleparams.set("gamma",gamma_);

    fluiddis_->SetState("velaf",fluidstate_.velaf_);
  }
  else
  {
    eleparams.set("action","calc_fluid_systemmat_and_residual");
    eleparams.set("using generalized-alpha time integration",false);
    eleparams.set("total time",time_);
    eleparams.set("is stationary", false);

    fluiddis_->SetState("velaf",fluidstate_.velnp_);
  }

  // call standard loop over elements
  fluiddis_->Evaluate(eleparams,fluidsysmat_,null,residual_,null,null);
  fluiddis_->ClearState();

  // account for potential Neumann inflow terms
  if (neumanninflow_)
  {
    // create parameter list
    ParameterList condparams;

    // action for elements
    condparams.set("action","calc_Neumann_inflow");
    condparams.set("thsl",theta_*dta_);
    condparams.set("Physical Type", physicaltype_);
    condparams.set("thermpress at n+alpha_F/n+1",thermpressaf_);

    // set vector values needed by elements
    fluiddis_->ClearState();
    fluiddis_->SetState("scaaf",fluidstate_.scaaf_);
    if (timealgo_==INPAR::FLUID::timeint_afgenalpha)
    {
      condparams.set("using generalized-alpha time integration",true);
      fluiddis_->SetState("velaf",fluidstate_.velaf_);
    }
    else
    {
      condparams.set("using generalized-alpha time integration",false);
      fluiddis_->SetState("velaf",fluidstate_.velnp_);
    }

    std::string condstring("FluidNeumannInflow");
    fluiddis_->EvaluateCondition(condparams,fluidsysmat_,Teuchos::null,residual_,Teuchos::null,Teuchos::null,condstring);
    fluiddis_->ClearState();
  }

  if (timealgo_==INPAR::FLUID::timeint_afgenalpha)
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
  fluidsysmat_->Complete();

  // end time measurement for element
  dtele_=Teuchos::Time::wallTime()-tcpu;

} // FluidXFluidImplicitTimeInt::AssembleMatAndRHS


////<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
////<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
////<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
///*----------------------------------------------------------------------*
// | update acceleration for generalized-alpha time integration  vg 02/09 |
// *----------------------------------------------------------------------*/
////<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
////<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
////<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidXFluidImplicitTimeInt::GenAlphaUpdateAcceleration()
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
  Teuchos::RCP<Epetra_Vector> onlyaccn  = velpressplitter_.ExtractOtherVector(fluidstate_.accn_ );
  Teuchos::RCP<Epetra_Vector> onlyveln  = velpressplitter_.ExtractOtherVector(fluidstate_.veln_ );
  Teuchos::RCP<Epetra_Vector> onlyvelnp = velpressplitter_.ExtractOtherVector(fluidstate_.velnp_);

  Teuchos::RCP<Epetra_Vector> onlyaccnp = rcp(new Epetra_Vector(onlyaccn->Map()));

  const double fact1 = 1.0/(gamma_*dta_);
  const double fact2 = 1.0 - (1.0/gamma_);
  onlyaccnp->Update(fact2,*onlyaccn,0.0);
  onlyaccnp->Update(fact1,*onlyvelnp,-fact1,*onlyveln,1.0);

  // copy back into global vector
  LINALG::Export(*onlyaccnp,*fluidstate_.accnp_);

} // FluidXFluidImplicitTimeInt::GenAlphaUpdateAcceleration

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
void FLD::FluidXFluidImplicitTimeInt::TimeUpdate()
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

    eleparams.set("dt"  ,dta_);

    // call loop over elements to update subgrid scales
    fluiddis_->Evaluate(eleparams,null,null,null,null,null);

    if(myrank_==0)
    {
      cout << "("<<Teuchos::Time::wallTime()-tcpu<<")\n";
    }
  }

  // compute accelerations
  {
    Teuchos::RCP<Epetra_Vector> onlyaccn  = velpressplitter_.ExtractOtherVector(fluidstate_.accn_ );
    Teuchos::RCP<Epetra_Vector> onlyaccnp = velpressplitter_.ExtractOtherVector(fluidstate_.accnp_);
    Teuchos::RCP<Epetra_Vector> onlyvelnm = velpressplitter_.ExtractOtherVector(fluidstate_.velnm_);
    Teuchos::RCP<Epetra_Vector> onlyveln  = velpressplitter_.ExtractOtherVector(fluidstate_.veln_ );
    Teuchos::RCP<Epetra_Vector> onlyvelnp = velpressplitter_.ExtractOtherVector(fluidstate_.velnp_);

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
    LINALG::Export(*onlyaccnp,*fluidstate_.accnp_);
  }

  // update old acceleration
  fluidstate_.accn_->Update(1.0,*fluidstate_.accnp_,0.0);

  // velocities/pressures of this step become most recent
  // velocities/pressures of the last step
  fluidstate_.velnm_->Update(1.0,*fluidstate_.veln_ ,0.0);
  fluidstate_.veln_ ->Update(1.0,*fluidstate_.velnp_,0.0);

  if (alefluid_)
  {
    fluidstate_.dispnm_->Update(1.0,*fluidstate_.dispn_,0.0);
    fluidstate_.dispn_ ->Update(1.0,*fluidstate_.dispnp_,0.0);
  }

  fluiddis_->ClearState();
  fluiddis_->SetState("velnp",fluidstate_.velnp_);
  fluiddis_->SetState("hist",hist_);

  return;
}// FluidXFluidImplicitTimeInt::TimeUpdate

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
void FLD::FluidXFluidImplicitTimeInt::StatisticsAndOutput()
{
  // time measurement: output and statistics
  TEUCHOS_FUNC_TIME_MONITOR("      + output and statistics");

  // -------------------------------------------------------------------
  //          calculate lift'n'drag forces from the residual
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

  return;
} // FluidXFluidImplicitTimeInt::StatisticsAndOutput


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | output of solution vector to binio                        gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidXFluidImplicitTimeInt::Output()
{

  //  ART_exp_timeInt_->Output();
  // output of solution
  if (step_%upres_ == 0)
  {
    // step number and time
    output_.NewStep(step_,time_);

    // velocity/pressure vector
    output_.WriteVector("velnp",fluidstate_.velnp_);

    // (hydrodynamic) pressure
    Teuchos::RCP<Epetra_Vector> pressure = velpressplitter_.ExtractCondVector(fluidstate_.velnp_);
    pressure->Scale(density_);
    output_.WriteVector("pressure", pressure);

    //output_.WriteVector("residual", trueresidual_);
    if (alefluid_) output_.WriteVector("dispnp", fluidstate_.dispnp_);

    //only perform stress calculation when output is needed
    if (writestresses_)
    {
      RCP<Epetra_Vector> traction = CalcStresses();
      output_.WriteVector("traction",traction);

      cout<<"Writing stresses"<<endl;
      //only perform wall shear stress calculation when output is needed
      if (write_wall_shear_stresses_)
      {
        RCP<Epetra_Vector> wss = CalcWallShearStresses();
        output_.WriteVector("wss",wss);
      }
    }


    // write domain decomposition for visualization (only once!)
    if (step_==upres_) output_.WriteElementData();

    if (uprestart_ != 0 && step_%uprestart_ == 0) //add restart data
    {
      // acceleration vector at time n+1 and n, velocity/pressure vector at time n and n-1
      output_.WriteVector("accnp",fluidstate_.accnp_);
      output_.WriteVector("accn", fluidstate_.accn_);
      output_.WriteVector("veln", fluidstate_.veln_);
      output_.WriteVector("velnm",fluidstate_.velnm_);

      if (alefluid_)
      {
        output_.WriteVector("dispn", fluidstate_.dispn_);
        output_.WriteVector("dispnm",fluidstate_.dispnm_);
      }
    }

  }
  // write restart also when uprestart_ is not a integer multiple of upres_
  else if (uprestart_ != 0 && step_%uprestart_ == 0)
  {
    // step number and time
    output_.NewStep(step_,time_);

    // velocity/pressure vector
    output_.WriteVector("velnp",fluidstate_.velnp_);

    //output_.WriteVector("residual", trueresidual_);
    if (alefluid_)
    {
      output_.WriteVector("dispnp", fluidstate_.dispnp_);
      output_.WriteVector("dispn", fluidstate_.dispn_);
      output_.WriteVector("dispnm",fluidstate_.dispnm_);
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
    output_.WriteVector("accnp",fluidstate_.accnp_);
    output_.WriteVector("accn", fluidstate_.accn_);
    output_.WriteVector("veln", fluidstate_.veln_);
    output_.WriteVector("velnm",fluidstate_.velnm_);
  }

  return;
} // FluidXFluidImplicitTimeInt::Output


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |                                                             kue 04/07|
 -----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidXFluidImplicitTimeInt::ReadRestart(int step)
{
  //  ART_exp_timeInt_->ReadRestart(step);
  IO::DiscretizationReader reader(fluiddis_,step);
  time_ = reader.ReadDouble("time");
  step_ = reader.ReadInt("step");

  reader.ReadVector(fluidstate_.velnp_,"velnp");
  reader.ReadVector(fluidstate_.veln_, "veln");
  reader.ReadVector(fluidstate_.velnm_,"velnm");
  reader.ReadVector(fluidstate_.accnp_,"accnp");
  reader.ReadVector(fluidstate_.accn_ ,"accn");

  if (alefluid_)
  {
    reader.ReadVector(fluidstate_.dispnp_,"dispnp");
    reader.ReadVector(fluidstate_.dispn_ , "dispn");
    reader.ReadVector(fluidstate_.dispnm_,"dispnm");
  }
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
void FLD::FluidXFluidImplicitTimeInt::UpdateGridv()
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
      gridv_->Update(1/dta_, *fluidstate_.dispnp_, -1/dta_, *fluidstate_.dispn_, 0.0);
    break;
    case 2:
      /* get gridvelocity from BDF2 time discretisation of mesh motion:
           -> requires one more previous mesh position or displacemnt
           -> somewhat more complicated
           -> allows second order accuracy for the overall flow solution  */
      gridv_->Update(1.5/dta_, *fluidstate_.dispnp_, -2.0/dta_, *fluidstate_.dispn_, 0.0);
      gridv_->Update(0.5/dta_, *fluidstate_.dispnm_, 1.0);
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
void FLD::FluidXFluidImplicitTimeInt::SetInitialFlowField(
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
    for(int lnodeid=0;lnodeid<fluiddis_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = fluiddis_->lRowNode(lnodeid);
      // the set of degrees of freedom associated with the node
      const vector<int> nodedofset = fluiddis_->Dof(lnode);

      for(int index=0;index<numdim_+1;++index)
      {
        int gid = nodedofset[index];

        double initialval=DRT::Problem::Instance()->Funct(startfuncno-1).Evaluate(index,lnode->X(),0.0,NULL);

        fluidstate_.velnp_->ReplaceGlobalValues(1,&initialval,&gid);
        fluidstate_.veln_ ->ReplaceGlobalValues(1,&initialval,&gid);
      }
    }

    // add random perturbation of certain percentage to function
    if (initfield == INPAR::FLUID::initfield_disturbed_field_from_function)
    {
      const Epetra_Map* dofrowmap = fluiddis_->DofRowMap();

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
      for(int lnodeid=0;lnodeid<fluiddis_->NumMyRowNodes();++lnodeid)
      {
        // get the processor local node
        DRT::Node*  lnode      = fluiddis_->lRowNode(lnodeid);
        // the set of degrees of freedom associated with the node
        vector<int> nodedofset = fluiddis_->Dof(lnode);

        for(int index=0;index<numdim_;++index)
        {
          int gid = nodedofset[index];
          int lid = dofrowmap->LID(gid);

          thisvel=(*fluidstate_.velnp_)[lid];
          if (mybmvel*mybmvel < thisvel*thisvel) mybmvel=thisvel;
        }
      }

      // the noise is proportional to the bulk mean velocity of the
      // undisturbed initial field (=2/3*maximum velocity)
      mybmvel=2*mybmvel/3;
      fluiddis_->Comm().MaxAll(&mybmvel,&bmvel,1);

      // loop all nodes on the processor
      for(int lnodeid=0;lnodeid<fluiddis_->NumMyRowNodes();++lnodeid)
      {
        // get the processor local node
        DRT::Node*  lnode      = fluiddis_->lRowNode(lnodeid);
        // the set of degrees of freedom associated with the node
        vector<int> nodedofset = fluiddis_->Dof(lnode);

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

          err += fluidstate_.velnp_->SumIntoGlobalValues(1,&noise,&gid);
          err += fluidstate_.veln_ ->SumIntoGlobalValues(1,&noise,&gid);
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
    const Epetra_Map* dofrowmap = fluiddis_->DofRowMap();

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
    for(int lnodeid=0;lnodeid<fluiddis_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = fluiddis_->lRowNode(lnodeid);

      // the set of degrees of freedom associated with the node
      vector<int> nodedofset = fluiddis_->Dof(lnode);

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
        err += fluidstate_.velnp_->ReplaceMyValues(1,&(u[nveldof]),&lid);
        err += fluidstate_.veln_ ->ReplaceMyValues(1,&(u[nveldof]),&lid);
        err += fluidstate_.velnm_->ReplaceMyValues(1,&(u[nveldof]),&lid);
      }
    } // end loop nodes lnodeid

    if (err!=0) dserror("dof not on proc");
  }
  // special initial function: Beltrami flow (3-D)
  else if (initfield == INPAR::FLUID::initfield_beltrami_flow)
  {
    const Epetra_Map* dofrowmap = fluiddis_->DofRowMap();

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
    for(int lnodeid=0;lnodeid<fluiddis_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = fluiddis_->lRowNode(lnodeid);

      // the set of degrees of freedom associated with the node
      vector<int> nodedofset = fluiddis_->Dof(lnode);

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
        err += fluidstate_.velnp_->ReplaceMyValues(1,&(u[nveldof]),&lid);
        err += fluidstate_.veln_ ->ReplaceMyValues(1,&(u[nveldof]),&lid);
        err += fluidstate_.velnm_->ReplaceMyValues(1,&(u[nveldof]),&lid);
      }

      // set initial pressure
      const int gid = nodedofset[npredof];
      int lid = dofrowmap->LID(gid);
      err += fluidstate_.velnp_->ReplaceMyValues(1,&p,&lid);
      err += fluidstate_.veln_ ->ReplaceMyValues(1,&p,&lid);
      err += fluidstate_.velnm_->ReplaceMyValues(1,&p,&lid);
    } // end loop nodes lnodeid

    if (err!=0) dserror("dof not on proc");
  }
  // special initial function: test case due to Bochev et al. (2007) (2-D)
  else if (initfield == INPAR::FLUID::initfield_bochev_test)
  {
    const Epetra_Map* dofrowmap = fluiddis_->DofRowMap();

    int err = 0;

    // check whether present flow is indeed two-dimensional
    if (numdim_!=2) dserror("Bochev test case is a two-dimensional flow!");

    // define vectors for velocity and pressure field as well as node coordinates
    vector<double> up(numdim_+1);
    vector<double> xy(numdim_);

    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<fluiddis_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = fluiddis_->lRowNode(lnodeid);

      // the set of degrees of freedom associated with the node
      vector<int> nodedofset = fluiddis_->Dof(lnode);

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
        err += fluidstate_.velnp_->ReplaceMyValues(1,&(up[ndof]),&lid);
        err += fluidstate_.veln_ ->ReplaceMyValues(1,&(up[ndof]),&lid);
        err += fluidstate_.velnm_->ReplaceMyValues(1,&(up[ndof]),&lid);
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

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | solve stationary fluid problem                              gjb 10/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidXFluidImplicitTimeInt::SolveStationaryProblem(
    const Teuchos::RCP<DRT::Discretization> fluidxfluidboundarydis
    )
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

    // -------------------------------------------------------------------
    //         evaluate dirichlet and neumann boundary conditions
    // -------------------------------------------------------------------
    {
      ParameterList eleparams;

      // other parameters needed by the elements
      eleparams.set("total time",time_);
      eleparams.set("delta time",dta_);
      eleparams.set("thsl",1.0); // no timefac in stationary case
      eleparams.set("form of convective term",convform_);
      eleparams.set("fs subgrid viscosity",fssgv_);
 //     eleparams.set("Physical Type", physicaltype_);

      eleparams.set("thermodynamic pressure",thermpressaf_);

      // set vector values needed by elements
      fluiddis_->ClearState();
      fluiddis_->SetState("velaf",fluidstate_.velnp_);
      // predicted dirichlet values
      // velnp then also holds prescribed new dirichlet values
      fluiddis_->EvaluateDirichlet(eleparams,fluidstate_.velnp_,null,null,null);

      fluiddis_->ClearState();

      // evaluate Neumann b.c.
      //eleparams.set("inc_density",density_);

      neumann_loads_->PutScalar(0.0);
      fluiddis_->SetState("scanp",fluidstate_.scaaf_);
      fluiddis_->EvaluateNeumann(eleparams,*neumann_loads_);
      fluiddis_->ClearState();
    }

    // -------------------------------------------------------------------
    //                     solve nonlinear equation system
    // -------------------------------------------------------------------
    NonlinearSolve(fluidxfluidboundarydis);

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

} // FluidXFluidImplicitTimeInt::SolveStationaryProblem


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  calculate traction vector at (Dirichlet) boundary (public) gjb 07/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
Teuchos::RCP<Epetra_Vector> FLD::FluidXFluidImplicitTimeInt::CalcStresses()
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

} // FluidXFluidImplicitTimeInt::CalcStresses()


/*----------------------------------------------------------------------*
 | Destructor dtor (public)                                  gammi 04/07|
 *----------------------------------------------------------------------*/
FLD::FluidXFluidImplicitTimeInt::~FluidXFluidImplicitTimeInt()
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
void FLD::FluidXFluidImplicitTimeInt::LiftDrag() const
{
  // in this map, the results of the lift drag calculation are stored
  RCP<map<int,vector<double> > > liftdragvals;

  FLD::UTILS::LiftDrag(*fluiddis_,*trueresidual_,params_,liftdragvals);

  if (liftdragvals!=Teuchos::null and fluiddis_->Comm().MyPID() == 0)
    FLD::UTILS::WriteLiftDragToFile(time_, step_, *liftdragvals);

  return;
}


/*----------------------------------------------------------------------*
 | compute flow rates through desired boundary parts        u.may 01/10 |
 *----------------------------------------------------------------------*/
void FLD::FluidXFluidImplicitTimeInt::ComputeFlowRates() const
{
  vector<DRT::Condition*> flowratecond;
  string condstring;

  if(numdim_ == 2)
  {
    condstring = "LineFlowRate";
    fluiddis_->GetCondition("LineFlowRate", flowratecond);
    // if no flowrate condition is present we do not compute anything
    if((int) flowratecond.size()== 0)
      return;
  }
  else if (numdim_ == 3)
  {
    condstring = "SurfFlowRate";
    fluiddis_->GetCondition("SurfFlowRate", flowratecond);
    // if no flowrate condition is present we do not compute anything
    if((int) flowratecond.size()== 0)
      return;
  }
  else
    dserror("flow rate computation is not implemented for the 1D case");

  const std::map<int,double> flowrates = FLD::UTILS::ComputeFlowRates(*fluiddis_, fluidstate_.velnp_, condstring);

  // write to file
  if(fluiddis_->Comm().MyPID() == 0)
    FLD::UTILS::WriteFlowRatesToFile(time_, step_, flowrates );

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FLD::FluidXFluidImplicitTimeInt::IntegrateInterfaceShape(std::string condname)
{
  ParameterList eleparams;
  // set action for elements
  eleparams.set("action","integrate_Shapefunction");

  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  const Epetra_Map* dofrowmap = fluiddis_->DofRowMap();

  // create vector (+ initialization with zeros)
  Teuchos::RCP<Epetra_Vector> integratedshapefunc = LINALG::CreateVector(*dofrowmap,true);

  // call loop over elements
  fluiddis_->ClearState();
  if (alefluid_)
  {
    fluiddis_->SetState("dispnp", fluidstate_.dispnp_);
  }
  fluiddis_->EvaluateCondition(eleparams,integratedshapefunc,condname);
  fluiddis_->ClearState();

  return integratedshapefunc;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::FluidXFluidImplicitTimeInt::AddDirichCond(const Teuchos::RCP<const Epetra_Map> maptoadd)
{
  std::vector<Teuchos::RCP<const Epetra_Map> > condmaps;
  condmaps.push_back(maptoadd);
  condmaps.push_back(dbcmaps_->CondMap());
  Teuchos::RCP<Epetra_Map> condmerged = LINALG::MultiMapExtractor::MergeMaps(condmaps);
  *dbcmaps_ = LINALG::MapExtractor(*(fluiddis_->DofRowMap()), condmerged);
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> FLD::FluidXFluidImplicitTimeInt::VelocityRowMap()
{ return velpressplitter_.OtherMap(); }

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> FLD::FluidXFluidImplicitTimeInt::PressureRowMap()
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
Teuchos::RCP<Epetra_Vector> FLD::FluidXFluidImplicitTimeInt::CalcWallShearStresses()
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
  const Epetra_Map* dofrowmap = fluiddis_->DofRowMap();
  
  //vector ndnorm0 with pressure-entries is needed for EvaluateCondition
  Teuchos::RCP<Epetra_Vector> ndnorm0 = LINALG::CreateVector(*dofrowmap,true);
  
  //call loop over elements, note: normal vectors do not yet have length = 1.0
  fluiddis_->ClearState();
  if (alefluid_)
  {
    fluiddis_->SetState("dispnp", fluidstate_.dispnp_);
  }
  // evaluate the normals of the surface
  fluiddis_->EvaluateCondition(eleparams,ndnorm0,"FluidStressCalc");
  fluiddis_->ClearState();

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
    
    // normalise the normal vector
    for (int j = 0; j < numdim_; j++)
    {
      (*ndnorm0)[i+j] /=  L;
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

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | setup the variables to do a new time loop                 u.kue 06/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidXFluidImplicitTimeInt::PrepareTimeLoop(const Teuchos::RCP<DRT::Discretization>  cutterdiscret)
{

  cout << "FLD::FluidXFluidImplicitTimeInt::PrepareTimeLoop()" << endl;

  {
    ParameterList eleparams;
    eleparams.set("action","reset");
    xfluiddis_->Evaluate(eleparams);
  }
  
  ComputeInterfaceAndSetDOFs(cutterdiscret);
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
Teuchos::RCP<XFEM::InterfaceHandleXFSI> FLD::FluidXFluidImplicitTimeInt::ComputeInterfaceAndSetDOFs(
    const Teuchos::RCP<DRT::Discretization>  cutterdiscret
    )
{
  // dump old matrix to save memory while we construct a new matrix
  xfluidsysmat_ = Teuchos::null;
  dofmanager_np_ = Teuchos::null;

  // within this routine, no parallel re-distribution is allowed to take place
  // before and after this function, it's ok to do that

  const Teuchos::ParameterList& xfemparams = DRT::Problem::Instance()->XFEMGeneralParams();
  const bool gmshdebugout = getIntegralValue<int>(xfemparams,"GMSH_DEBUG_OUT")==1;
  const bool screen_out = getIntegralValue<int>(xfemparams,"GMSH_DEBUG_OUT_SCREEN")==1;

  if (gmshdebugout)
  {
    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("Fluid_Fluid_Coupling", 1, 0, screen_out, xfluiddis_->Comm().MyPID());
    std::ofstream gmshfilecontent(filename.c_str());
    IO::GMSH::disToStream("XFluidFluidboundary", 0.0, cutterdiscret,gmshfilecontent);
    IO::GMSH::disToStream("XFluid", 0.0, xfluiddis_, gmshfilecontent);
    IO::GMSH::disToStream("Fluid", 0.0, fluiddis_,gmshfilecontent);
    gmshfilecontent.close();
  }

  //dummy vectors / here not needed
  const vector<int> MovingFluideleGIDs;
  const vector<int> MovingFluidNodeGIDs;
  
   // compute Intersection
  ih_np_ = rcp(new XFEM::InterfaceHandleXFSI(xfluiddis_, cutterdiscret, MovingFluideleGIDs));

  ih_np_->toGmsh(step_);

  // apply enrichments
  dofmanager_np_ = rcp(new XFEM::DofManager(ih_np_, physprob_.fieldset_, *physprob_.elementAnsatz_, xparams_, MovingFluidNodeGIDs));

  // tell elements about the dofs and the integration
  TransferDofInformationToElements(ih_np_, dofmanager_np_);

  // print global and element dofmanager to Gmsh
  dofmanager_np_->toGmsh(step_);


  // get old dofmap, compute new one and get the new one, too
  const Epetra_Map olddofrowmap = *xfluiddis_->DofRowMap();
  xfluiddis_->FillComplete(true,false,true);
  const Epetra_Map& newdofrowmap = *xfluiddis_->DofRowMap();

  {
    const std::map<XFEM::DofKey<XFEM::onNode>, XFEM::DofGID> oldNodalDofDistributionMap(xfluidstate_.nodalDofDistributionMap_);
    const std::map<XFEM::DofKey<XFEM::onElem>, XFEM::DofGID> oldElementalDofDistributionMap(xfluidstate_.elementalDofDistributionMap_);
    dofmanager_np_->fillDofRowDistributionMaps(xfluidstate_.nodalDofDistributionMap_,xfluidstate_.elementalDofDistributionMap_);

    // create switcher
    const XFEM::DofDistributionSwitcher dofswitch(
            ih_np_, dofmanager_np_,
            olddofrowmap, newdofrowmap,
            oldNodalDofDistributionMap,xfluidstate_.nodalDofDistributionMap_,
            oldElementalDofDistributionMap,xfluidstate_.elementalDofDistributionMap_
            );

    // --------------------------------------------
    // switch state vectors to new dof distribution
    // --------------------------------------------

   // cout0_ << " ->  Initialize system vectors..." << endl;
    // accelerations
    dofswitch.mapVectorToNewDofDistribution(xfluidstate_.accnp_);
    dofswitch.mapVectorToNewDofDistribution(xfluidstate_.accn_);

    // velocities and pressures
    dofswitch.mapVectorToNewDofDistribution(xfluidstate_.velnp_);
    dofswitch.mapVectorToNewDofDistribution(xfluidstate_.veln_);
    dofswitch.mapVectorToNewDofDistribution(xfluidstate_.velnm_);

    // if dofs appear, extrapolate from the interface to the newly created dofs
    if (Step() > 1)
    {
      dofswitch.extrapolateOldTimeStepValues(ih_np_->cutterdis(), *ih_np_->cutterposn(), ih_np_->cutterdis()->GetState("ivelcoln") , xfluidstate_.veln_ );
      dofswitch.extrapolateOldTimeStepValues(ih_np_->cutterdis(), *ih_np_->cutterposn(), ih_np_->cutterdis()->GetState("ivelcolnm"), xfluidstate_.velnm_);
      dofswitch.extrapolateOldTimeStepValues(ih_np_->cutterdis(), *ih_np_->cutterposn(), ih_np_->cutterdis()->GetState("iacccoln") , xfluidstate_.accn_ );
    }
    PlotVectorFieldToGmsh(DRT::UTILS::GetColVersionOfRowVector(xfluiddis_, xfluidstate_.veln_) , "sol_field_veln_extrapolated" ,"Velocity Solution (Physical) n"    ,false, Step(), Time());
    PlotVectorFieldToGmsh(DRT::UTILS::GetColVersionOfRowVector(xfluiddis_, xfluidstate_.velnm_), "sol_field_velnm_extrapolated","Velocity Solution (Physical) n-1"  ,false, Step(), Time());
    PlotVectorFieldToGmsh(DRT::UTILS::GetColVersionOfRowVector(xfluiddis_, xfluidstate_.accn_) , "sol_field_accn_extrapolated" ,"Acceleration Solution (Physical) n",false, Step(), Time());
  }

  // --------------------------------------------------
  // create remaining vectors with new dof distribution
  // --------------------------------------------------
  //hist_         = LINALG::CreateVector(newdofrowmap,true);

  // object holds maps/subsets for DOFs subjected to Dirichlet BCs and otherwise
  dbcmaps_ = Teuchos::rcp(new LINALG::MapExtractor());
  {
    ParameterList eleparams;
    // other parameters needed by the elements
    eleparams.set("total time",time_);
    Teuchos::RCP<Epetra_Vector> tmp = LINALG::CreateVector(newdofrowmap,true);
    xfluiddis_->EvaluateDirichletXFEM(eleparams, tmp, Teuchos::null, Teuchos::null,
                                Teuchos::null, dbcmaps_);
  }

  // -------------------------------------------------------------------
  // get a vector layout from the discretization for a vector which only
  // contains the velocity dofs and for one vector which only contains
  // pressure degrees of freedom.
  // -------------------------------------------------------------------
  FLD::UTILS::SetupXFluidSplit(*xfluiddis_,dofmanager_np_,velpressplitter_);

//  // project old interpolated velocity vector onto divergence free space
//  if (xparams_.get<bool>("INCOMP_PROJECTION"))
//    if (timealgo_ != INPAR::FLUID::timeint_stationary)
//      ProjectOldTimeStepValues(velpressplitter_);

  neumann_loads_= LINALG::CreateVector(newdofrowmap,true);

  // ---------------------------------
  // Vectors used for solution process
  // ---------------------------------
  residual_     = LINALG::CreateVector(newdofrowmap,true);
  trueresidual_ = LINALG::CreateVector(newdofrowmap,true);
  incvel_       = LINALG::CreateVector(newdofrowmap,true);

  // -------------------------------------------------------------------
  // create empty system matrix --- stiffness and mass are assembled in
  // one system matrix!
  // -------------------------------------------------------------------
  //cout0_ << " ->  Initialize system matrix..." << endl;

  // initialize system matrix
  xfluidsysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(newdofrowmap,0,false,true));

  // print information about dofs
  if (xfluiddis_->Comm().NumProc() == 1)
  {
    const int numdof = newdofrowmap.NumGlobalElements();
    const int numnodaldof = dofmanager_np_->NumNodalDof();
    cout << "Xfluid DOF report: numdof (nodal and elemental) = " << numdof << ", numstressdof = "<< (numdof - numnodaldof) << endl;
  }
  else
  {
    if(xfluiddis_->Comm().MyPID() == 0)
    {
      cout << "XFluid DOF report: numdof = " << newdofrowmap.NumGlobalElements() << endl;
    }
  }

 // cout0_ << "Setup phase done!" << endl;

  return ih_np_;
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidXFluidImplicitTimeInt::TransferDofInformationToElements(
    const Teuchos::RCP<XFEM::InterfaceHandleXFSI>  ih,
    const Teuchos::RCP<XFEM::DofManager> dofmanager
    )
{
  ParameterList eleparams;
  eleparams.set("action","store_xfem_info");
  eleparams.set("dofmanager",dofmanager);
  eleparams.set("DLM_condensation",xparams_.get<bool>("DLM_condensation"));
  eleparams.set("boundaryRatioLimit",xparams_.get<double>("boundaryRatioLimit"));
  eleparams.set("volumeRatioLimit",xparams_.get<double>("volumeRatioLimit"));
  eleparams.set("EMBEDDED_BOUNDARY",xparams_.get<INPAR::XFEM::BoundaryIntegralType>("EMBEDDED_BOUNDARY"));
  eleparams.set("interfacehandle",ih);
  xfluiddis_->Evaluate(eleparams);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FLD::FluidXFluidImplicitTimeInt::PlotVectorFieldToGmsh(
    const Teuchos::RCP<const Epetra_Vector>   vectorfield,
    const std::string filestr,
    const std::string name_in_gmsh,
    const bool plot_to_gnuplot,
    const int step,
    const double time
    ) const
{

  const Teuchos::ParameterList& xfemparams = DRT::Problem::Instance()->XFEMGeneralParams();
  const bool gmshdebugout = getIntegralValue<int>(xfemparams,"GMSH_DEBUG_OUT")==1;
  const bool screen_out = getIntegralValue<int>(xfemparams,"GMSH_DEBUG_OUT_SCREEN")==1;

  if (gmshdebugout)
  {
    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles(filestr, step, 5, screen_out, xfluiddis_->Comm().MyPID());
    std::ofstream gmshfilecontent(filename.c_str());

    {
      gmshfilecontent << "View \" " << name_in_gmsh << "\" {\n";
      for (int i=0; i<xfluiddis_->NumMyColElements(); ++i)
      {
        const DRT::Element* actele = xfluiddis_->lColElement(i);

        // create local copy of information about dofs
        const XFEM::ElementDofManager eledofman(*actele, physprob_.elementAnsatz_->getElementAnsatz(actele->Shape()), *dofmanager_np_);

        vector<int> lm;
        vector<int> lmowner;
        actele->LocationVector(*xfluiddis_, lm, lmowner);

        // extract local values from the global vector
        vector<double> myvelnp(lm.size());
        DRT::UTILS::ExtractMyValues(*vectorfield, myvelnp, lm);

        const vector<int>& dofposvelx =
          eledofman.LocalDofPosPerField(XFEM::PHYSICS::Velx);
        const vector<int>& dofposvely =
          eledofman.LocalDofPosPerField(XFEM::PHYSICS::Vely);
        const vector<int>& dofposvelz =
          eledofman.LocalDofPosPerField(XFEM::PHYSICS::Velz);

        const int numparam = eledofman.NumDofPerField(XFEM::PHYSICS::Velx);
        LINALG::SerialDenseMatrix elementvalues(3, numparam);
        for (int iparam=0; iparam<numparam; ++iparam)
        {
          elementvalues(0, iparam) = myvelnp[dofposvelx[iparam]];
          elementvalues(1, iparam) = myvelnp[dofposvely[iparam]];
          elementvalues(2, iparam) = myvelnp[dofposvelz[iparam]];
        }

          const GEO::DomainIntCells& domainintcells =
            dofmanager_np_->getInterfaceHandle()->GetDomainIntCells(actele);
          for (GEO::DomainIntCells::const_iterator cell =
            domainintcells.begin(); cell != domainintcells.end(); ++cell)
          {
            LINALG::SerialDenseMatrix cellvalues(3, DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
            XFEM::computeVectorCellNodeValues(*actele, &*ih_np_, eledofman,
                *cell, XFEM::PHYSICS::Velx, elementvalues, cellvalues);
            IO::GMSH::cellWithVectorFieldToStream(
                cell->Shape(), cellvalues, cell->CellNodalPosXYZ(), gmshfilecontent);
          }
//          const GEO::BoundaryIntCells& boundaryintcells =
//              ih_np_->GetBoundaryIntCells(actele->Id());
//          // draw boundary integration cells with values
//          for (GEO::BoundaryIntCells::const_iterator cell =
//            boundaryintcells.begin(); cell != boundaryintcells.end(); ++cell)
//          {
//            LINALG::SerialDenseMatrix cellvalues(3, DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
//
//            const DRT::Element* boundaryele = ih_np_->GetBoundaryEle(cell->GetSurfaceEleGid());
//            const int label = ih_np_->GetLabelPerBoundaryElementId(boundaryele->Id());
//
//            XFEM::computeVectorCellNodeValues(*actele, ih_np_, eledofman,
//                *cell, XFEM::PHYSICS::Velx, label, elementvalues, cellvalues);
//            IO::GMSH::cellWithVectorFieldToStream(
//                cell->Shape(), cellvalues, cell->CellNodalPosXYZ(), gmshfilecontent);
//          }

          // draw uncut element
          {
            LINALG::SerialDenseMatrix elevalues(3, DRT::UTILS::getNumberOfElementNodes(actele->Shape()),true);
            static LINALG::Matrix<3,27> xyze_ele;
            GEO::fillInitialPositionArray(actele, xyze_ele);
            IO::GMSH::cellWithVectorFieldToStream(
                            actele->Shape(), elevalues, xyze_ele, gmshfilecontent);
          }

//        }
        //if (dofmanager_np_->getInterfaceHandle()->ElementIntersected(elegid) and not ele_to_textfile and ele_to_textfile2)
//        if (elegid == 14 and elementvalues.N() > 0 and plot_to_gnuplot)
        if (actele->Id() == 1 and elementvalues.N() > 0 and plot_to_gnuplot)
        {
          std::ofstream f;
          const std::string fname = DRT::Problem::Instance()->OutputControlFile()->FileName() + "."+ filestr + ".outflow_gnuplot.txt";
          if (step <= 1)
            f.open(fname.c_str(),std::fstream::trunc);
          else
            f.open(fname.c_str(),std::fstream::ate | std::fstream::app);
          f << time << "  " << elementvalues(0,0) << "\n";
          f.close();
        }
      }
      gmshfilecontent << "};\n";
    }
    gmshfilecontent.close();
    if (screen_out) std::cout << " done" << endl;
  }
}
#endif /* CCADISCRET       */
