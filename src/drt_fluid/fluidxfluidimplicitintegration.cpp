/*----------------------------------------------------------------------*/
/*!
\file fluidxfluidimplicitintegration.cpp


<pre>
maintainer: shadan shahmiri
            shahmiri@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef  CCADISCRET


#undef   WRITEOUTSTATISTICS

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
#include "../drt_fluid/fluid_utils_mapextractor.H"

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
  fluidsolverparams_(DRT::Problem::Instance()->FluidSolverParams()),
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

  xfluidoutput_ = rcp(new IO::DiscretizationWriter(xfluiddis_));

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

  if (timealgo_== INPAR::FLUID::timeint_afgenalpha or timealgo_ == INPAR::FLUID::timeint_gen_alpha)
    dserror("Generalized alpha integration scheme is not available in FluidXFluid coubling!!");

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

  FLD::UTILS::SetupFluidSplit(*fluiddis_,numdim_,fluidvelpressplitter_);

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
    Teuchos::rcp(new LINALG::BlockSparseMatrix<FLD::UTILS::VelPressSplitStrategy>(fluidvelpressplitter_,fluidvelpressplitter_,108,false,true));
    blocksysmat->SetNumdim(numdim_);
    fluidsysmat_ = blocksysmat;
  }

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
    fluidstate_.dispnmm_ = LINALG::CreateVector(*fluiddofrowmap,true);

    gridv_  = LINALG::CreateVector(*fluiddofrowmap,true);
  }

  // Vectors associated to boundary conditions
  // -----------------------------------------

  // a vector of zeros to be used to enforce zero dirichlet boundary conditions
  fluidzeros_   = LINALG::CreateVector(*fluiddofrowmap,true);

  // object holds maps/subsets for DOFs subjected to Dirichlet BCs and otherwise
  fluiddbcmaps_ = Teuchos::rcp(new LINALG::MapExtractor());
  {
    ParameterList eleparams;
    // other parameters needed by the elements
    eleparams.set("total time",time_);
    fluiddis_->EvaluateDirichlet(eleparams, fluidzeros_, Teuchos::null, Teuchos::null,
                                Teuchos::null, fluiddbcmaps_);

    fluidzeros_->PutScalar(0.0); // just in case of change
  }

  // the vector containing body and surface forces
  fluid_neumann_loads_= LINALG::CreateVector(*fluiddofrowmap,true);

  // Vectors used for solution process
  // ---------------------------------
  fluidresidual_     = LINALG::CreateVector(*fluiddofrowmap,true);
  fluidtrueresidual_ = LINALG::CreateVector(*fluiddofrowmap,true);

  // Nonlinear iteration increment vector
  fluidincvel_ = LINALG::CreateVector(*fluiddofrowmap,true);

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

  //xfluid vectors
  {
    ParameterList eleparams;
    eleparams.set("action","set_output_mode");
    eleparams.set("output_mode",true);
    xfluiddis_->Evaluate(eleparams);
  }

  xfluiddis_->FillComplete();

  xfluidoutput_->WriteMesh(0,0.0);

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

  // ---------------------------------------------------------------------
  // set general fluid parameter defined before
  // ---------------------------------------------------------------------
  SetElementGeneralFluidParameter();

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
  cout << "FluidXFluidImplicitTimeInt::TimeLoop()" << endl;

  // time measurement: time loop
  TEUCHOS_FUNC_TIME_MONITOR(" + time loop");

  const Epetra_Map* fxfboundary_dofcolmap = fluidxfluidboundarydis->DofColMap();
  fluidxfluidboundarydis->SetState("idispcolnp", LINALG::CreateVector(*fxfboundary_dofcolmap,true));
  fluidxfluidboundarydis->SetState("ivelcolnp", LINALG::CreateVector(*fxfboundary_dofcolmap,true));
  fluidxfluidboundarydis->SetState("idispcoln", LINALG::CreateVector(*fxfboundary_dofcolmap,true));
  fluidxfluidboundarydis->SetState("ivelcoln", LINALG::CreateVector(*fxfboundary_dofcolmap,true));
  fluidxfluidboundarydis->SetState("ivelcolnm", LINALG::CreateVector(*fxfboundary_dofcolmap,true));
  fluidxfluidboundarydis->SetState("iacccoln", LINALG::CreateVector(*fxfboundary_dofcolmap,true));

  while (step_<stepmax_ or time_<maxtime_)
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
  cout << "FLD::FluidXFluidImplicitTimeInt::PrepareTimeStep()" << endl;

  // update interface handle
  ih_n_ = ih_np_;

  xfluidstate_.nodelabeln_ = xfluidstate_.nodelabelnp_;

  const Teuchos::ParameterList& xfemparams = DRT::Problem::Instance()->XFEMGeneralParams();
  const bool gmshdebugout = getIntegralValue<int>(xfemparams,"GMSH_DEBUG_OUT")==1;

  // Gmsh_Output for Node positions in old time step
  if (gmshdebugout)
  {
    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("node_positions",step_, 5, false, xfluiddis_->Comm().MyPID());
    std::ofstream gmshfilecontent(filename.c_str());
    gmshfilecontent << "View \" " << "Node->Position \" {\n";
    for (int i=0; i<xfluiddis_->NumMyColNodes(); ++i)
    {
      const DRT::Node* actnode = xfluiddis_->lColNode(i);
      const LINALG::Matrix<3,1> pos(actnode->X());
      map<int,int>::const_iterator iter = xfluidstate_.nodelabeln_.find(actnode->Id());
      int label = iter->second;
      IO::GMSH::cellWithScalarToStream(DRT::Element::point1, label, pos, gmshfilecontent);
    }
    gmshfilecontent << "};\n";
  }

  // -------------------------------------------------------------------
  //              set time dependent parameters
  // -------------------------------------------------------------------
  step_ += 1;
  time_ += dta_;

  if (params_.get<INPAR::FLUID::TimeIntegrationScheme>("time int algo") == INPAR::FLUID::timeint_stationary)
  {
    timealgo_ = INPAR::FLUID::timeint_stationary;
    theta_ = 1.0;
  }
  else
  {
    // do a backward Euler step for the first timestep
    if (step_==1)
    {
      timealgo_ = INPAR::FLUID::timeint_one_step_theta;
      theta_ = params_.get<double>("start theta");
    }
    else
    {
      timealgo_ = params_.get<INPAR::FLUID::TimeIntegrationScheme>("time int algo");
      theta_ = params_.get<double>("theta");

      // for BDF2, theta is set by the time-step sizes, 2/3 for const. dt
      if (timealgo_==INPAR::FLUID::timeint_bdf2) theta_ = (dta_+dtp_)/(2.0*dta_ + dtp_);
    }
  }

  // -------------------------------------------------------------------
  //  Set time parameter for element call
  // -------------------------------------------------------------------
  SetElementTimeParameter();

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
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidXFluidImplicitTimeInt::PrepareNonlinearSolve()
{
  cout << "FLD::FluidXFluidImplicitTimeInt::PrepareNonlinearSolve()" << endl;

  // -------------------------------------------------------------------
  //         evaluate Dirichlet and Neumann boundary conditions for xfluid
  // -------------------------------------------------------------------
  {
    ParameterList eleparams;

    // total time required for Dirichlet conditions
    eleparams.set("total time",time_);
    eleparams.set("delta time",dta_);
    eleparams.set("thsl",theta_*dta_);

    // set vector values needed by elements
    xfluiddis_->ClearState();
    xfluiddis_->SetState("velnp",xfluidstate_.velnp_);
    // predicted dirichlet values
    // velnp then also holds prescribed new dirichlet values
    xfluiddis_->EvaluateDirichletXFEM(eleparams,xfluidstate_.velnp_,null,null,null,xfluiddbcmaps_);
    xfluiddis_->ClearState();

    // set all parameters and states required for Neumann conditions
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

    xfluid_neumann_loads_->PutScalar(0.0);

      // evaluate Neumann conditions
    xfluiddis_->EvaluateNeumann(eleparams,*xfluid_neumann_loads_);
    xfluiddis_->ClearState();
  }

  // -------------------------------------------------------------------
  //  evaluate Dirichlet and Neumann boundary conditions
  // -------------------------------------------------------------------
  {
    ParameterList eleparams;

    // total time required for Dirichlet conditions
    eleparams.set("total time",time_);

    // set vector values needed by elements
    fluiddis_->SetState("velnp",fluidstate_.velnp_);
    // predicted dirichlet values
    // velnp then also holds prescribed new dirichlet values
    fluiddis_->EvaluateDirichlet(eleparams,fluidstate_.velnp_,null,null,null);

    // set all parameters and states required for Neumann conditions
 //   eleparams.set("Physical Type",physicaltype_); // brauche ich das?!!
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
    fluid_neumann_loads_->PutScalar(0.0);
    fluiddis_->SetState("scanp",fluidstate_.scaaf_); // brauche ich das?!!
    fluiddis_->EvaluateNeumann(eleparams,*fluid_neumann_loads_);
    fluiddis_->ClearState();
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

  PrepareFluidXFluidBoundaryDis(fluidxfluidboundarydis);
  ComputeFluidXFluidInterfaceAccelerationsAndVelocities();
  ComputeInterfaceAndSetDOFs(fluidxfluidboundarydis);

  // merge the fluid and xfluid maps
  RCP<Epetra_Map> fluiddofrowmap = rcp(new Epetra_Map(*fluiddis_->DofRowMap()));
  RCP<Epetra_Map> xfluiddofrowmap = rcp(new Epetra_Map(*xfluiddis_->DofRowMap()));
  std::vector<Teuchos::RCP<const Epetra_Map> > maps;
  maps.push_back(fluiddofrowmap);
  maps.push_back(xfluiddofrowmap);
  fluidxfluidrowmap_ = LINALG::MultiMapExtractor::MergeMaps(maps);
  fluidxfluidsplitter_.Setup(*fluidxfluidrowmap_,fluiddofrowmap,xfluiddofrowmap);

  fluidresidual_  = LINALG::CreateVector(*fluiddofrowmap,true);
  fluidsysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*fluiddofrowmap,108,false,true));

  // split velocity and pressure
  FLD::UTILS::SetupFluidXFluidVelPresSplit(*fluiddis_,numdim_,*xfluiddis_,dofmanager_np_,fluidxfluidvelpressplitter_,fluidxfluidrowmap_);

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
    w_       = LINALG::CreateVector(*fluidxfluidrowmap_,true);
    c_       = LINALG::CreateVector(*fluidxfluidrowmap_,true);
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

  Teuchos::RCP<LINALG::Solver> fluidxfluidsolver =
        rcp(new LINALG::Solver(fluidsolverparams_,
            xfluiddis_->Comm(),
            DRT::Problem::Instance()->ErrorFile()->Handle()));

  // sanity check
  if (fluidxfluidsolver->Params().get<string>("solver") == "aztec"
      and fluidxfluidsolver->Params().isSublist("ML Parameters")
      and not xparams_.get<bool>("DLM_condensation")
      )
  {
      dserror("for MLFLUID2, you need to set \"DLM_CONDENSATION  yes\" ");
  }

  PrepareNonlinearSolve();

  fluidxfluidsysmat_   = Teuchos::rcp(new LINALG::SparseMatrix(*fluidxfluidrowmap_,0,false,true));
  fluidxfluidresidual_ = LINALG::CreateVector(*fluidxfluidrowmap_,true);
  fluidxfluidincvel_   = LINALG::CreateVector(*fluidxfluidrowmap_,true);

  fluidxfluidstate_.velnp_ = LINALG::CreateVector(*fluidxfluidrowmap_,true);

  // ---------------------------------------------- nonlinear iteration
  // ------------------------------- stop nonlinear iteration when both
  //                                 increment-norms are below this bound
  const double  ittol     =params_.get<double>("tolerance for nonlin iter");

  //------------------------------ turn adaptive solver tolerance on/off
  const bool   isadapttol    = params_.get<bool>("ADAPTCONV",true);
  const double adaptolbetter = params_.get<double>("ADAPTCONV_BETTER",0.01);


  int  itnum = 0;
  int  itemax = params_.get<int>   ("max nonlin iter steps");
  bool stopnonliniter = false;

  dtsolve_  = 0.0;
  dtele_    = 0.0;
 // dtfilter_ = 0.0;

  if (myrank_ == 0)
  {
    std::cout << "+------------+------------------+---------+---------+---------+---------+---------+---------+\n";
    std::cout << "|  step/max  |  tol     [norm]  | vel-res | pre-res | fullres | vel-inc | pre-inc | fullinc |" << std::endl;
  }

  const Teuchos::RCP<Epetra_Vector> iforcecolnp = LINALG::CreateVector(*fluidxfluidboundarydis->DofColMap(),true);

  fluidincvel_->PutScalar(0.0);
  xfluidresidual_->PutScalar(0.0);
  fluidresidual_->PutScalar(0.0);
  fluidxfluidresidual_->PutScalar(0.0);

  // increment of the old iteration step - used for update of condensed element stresses
  Teuchos::RCP<Epetra_Vector> oldinc = LINALG::CreateVector(*xfluiddis_->DofRowMap(),true);

  while (stopnonliniter==false)
  {
    itnum++;

    // Insert fluid and xfluid vectors to fluidxfluid
    fluidxfluidsplitter_.InsertXFluidVector(xfluidstate_.velnp_,fluidxfluidstate_.velnp_);
    fluidxfluidsplitter_.InsertFluidVector(fluidstate_.velnp_,fluidxfluidstate_.velnp_);

    // -------------------------------------------------------------------
    // call elements to calculate system matrix and RHS
    // -------------------------------------------------------------------
    {
      // time measurement: element
      TEUCHOS_FUNC_TIME_MONITOR("      + element calls");

      // get cpu time
      const double tcpu=Teuchos::Time::wallTime();

      fluidsysmat_->Zero();
      xfluidsysmat_->Zero();
      fluidxfluidsysmat_->Zero();

      // create the parameters for the discretization
      ParameterList eleparams;

      // add Neumann loads
      fluidresidual_->Update(1.0,*fluid_neumann_loads_,0.0);
      xfluidresidual_->Update(1.0,*xfluid_neumann_loads_,0.0);

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
      fluiddis_->SetState("velnp",fluidstate_.velnp_);
      fluiddis_->SetState("veln" ,fluidstate_.veln_);
      fluiddis_->SetState("velnm",fluidstate_.velnm_);
      fluiddis_->SetState("accn" ,fluidstate_.accn_);

      eleparams.set("timealgo",timealgo_);
      if(params_.get<INPAR::FLUID::LinearisationAction>("Linearisation") == INPAR::FLUID::Newton)
        eleparams.set("include reactive terms for linearisation",true);
      else if (params_.get<INPAR::FLUID::LinearisationAction>("Linearisation") == INPAR::FLUID::minimal)
        dserror("LinearisationAction minimal is not defined in the XFEM formulation");
      else
        eleparams.set("include reactive terms for linearisation",false);

      xfluiddis_->ClearState();
      xfluiddis_->SetState("velnp",xfluidstate_.velnp_);
      xfluiddis_->SetState("veln" ,xfluidstate_.veln_);
      xfluiddis_->SetState("velnm",xfluidstate_.velnm_);
      xfluiddis_->SetState("accn" ,xfluidstate_.accn_);

      fluidxfluidboundarydis->SetState("iveln"  ,fxfiveln_);
      fluidxfluidboundarydis->SetState("ivelnm" ,fxfivelnm_);
      fluidxfluidboundarydis->SetState("iaccn"  ,fxfiaccn_);

      if (alefluid_)
      {
        fluiddis_->SetState("dispnp", fluidstate_.dispnp_);
        fluiddis_->SetState("gridv", gridv_);
      }

      xfluiddis_->SetState("velpres nodal iterinc",oldinc);

      // reset interface force and let the elements fill it
      iforcecolnp->PutScalar(0.0);
      eleparams.set("interface force",iforcecolnp);

      // set scheme-specific element parameters and vector values
      if (timealgo_==INPAR::FLUID::timeint_stationary)
      {
        eleparams.set("action","calc_fluid_systemmat_and_residual");
        eleparams.set("using generalized-alpha time integration",false);
        eleparams.set("total time",time_);
        eleparams.set("is stationary", true);

        fluiddis_->SetState("velaf",fluidstate_.velnp_);
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
        // call standard loop over fluid-elements
        fluiddis_->Evaluate(eleparams,fluidsysmat_,null,fluidresidual_,null,null);

        double L2 = 0.0;
        eleparams.set("L2",L2);
        eleparams.set("DLM_condensation",xparams_.get<bool>("DLM_condensation"));
        eleparams.set("monolithic_FSI",false);
        eleparams.set("EMBEDDED_BOUNDARY",xparams_.get<INPAR::XFEM::BoundaryIntegralType>("EMBEDDED_BOUNDARY"));
        fluidxfluidboundarydis->SetState("interface nodal iterinc", fxfiincvel_);

        eleparams.set("action","fluidxfluidCoupling");

        Cud_ = Teuchos::rcp(new LINALG::SparseMatrix(*xfluiddofrowmap,0,false,false));
        Cdu_ = Teuchos::rcp(new LINALG::SparseMatrix(*fluidxfluidboundarydis->DofRowMap(),0,false,false));
        Cdd_ = Teuchos::rcp(new LINALG::SparseMatrix(*fluidxfluidboundarydis->DofRowMap(),0,false,false));
        rhsd_ = LINALG::CreateVector(*fluidxfluidboundarydis->DofRowMap(),true);

        // call standard loop over xfluid-elements
        MonolithicMultiDisEvaluate(ih_np_, xfluiddis_, fluidxfluidboundarydis, eleparams,
                    Cud_, Cdu_, Cdd_, rhsd_, xfluidsysmat_, xfluidresidual_, true, true);

        Cud_->Complete(*fluidxfluidboundarydis->DofRowMap(),*xfluiddofrowmap);
        Cdu_->Complete(*xfluiddofrowmap,*fluidxfluidboundarydis->DofRowMap());
        Cdd_->Complete(*fluidxfluidboundarydis->DofRowMap(),*fluidxfluidboundarydis->DofRowMap());

        // adding rhsd_ to fluidresidual_
        // in *vector[LID] immer mit LID!!! Hier rhsd_->Map().LID(rhsdgid) muss wieder gleich iter sein!
        for (int iter=0; iter<rhsd_->MyLength();++iter)
        {
          int rhsdgid = rhsd_->Map().GID(iter);
          if (rhsd_->Map().MyGID(rhsdgid) == false) dserror("rhsd_ should be on all prossesors");
          if (fluidresidual_->Map().MyGID(rhsdgid))
            (*fluidresidual_)[fluidresidual_->Map().LID(rhsdgid)]=(*fluidresidual_)[fluidresidual_->Map().LID(rhsdgid)] +
            (*rhsd_)[rhsd_->Map().LID(rhsdgid)];
        }

        // account for potential Neumann inflow terms
        if (neumanninflow_)
        {
          // create parameter list
          ParameterList condparams;

          // action for elements
          condparams.set("action","calc_Neumann_inflow");
          condparams.set("thsl",theta_*dta_);
          condparams.set("thermpress at n+alpha_F/n+1",thermpressaf_);

          // set vector values needed by elements
          fluiddis_->ClearState();
          xfluiddis_->ClearState();

          condparams.set("using generalized-alpha time integration",false);
          fluiddis_->SetState("velaf",fluidstate_.velnp_);
          fluiddis_->SetState("velaf",xfluidstate_.velnp_);

          std::string condstring("FluidNeumannInflow");
          fluiddis_->EvaluateCondition(condparams,fluidsysmat_,Teuchos::null,fluidresidual_,Teuchos::null,Teuchos::null,condstring);
          fluiddis_->ClearState();

          xfluiddis_->EvaluateCondition(condparams,xfluidsysmat_,Teuchos::null,xfluidresidual_,Teuchos::null,Teuchos::null,condstring);
          xfluiddis_->ClearState();
        }

        // scaling to get true residual vector for all other schemes
        fluidtrueresidual_->Update(ResidualScaling(),*fluidresidual_,0.0);
        xfluidtrueresidual_->Update(ResidualScaling(),*xfluidresidual_,0.0);

        xfluiddis_->ClearState();
        fluiddis_->ClearState();

        // finalize the complete matrix
        fluidsysmat_->Complete();
        xfluidsysmat_->Complete();
        fluidxfluidsysmat_->Zero();
      }

      // end time measurement for element
      xfluiddis_->Comm().Barrier();
      dtele_=Teuchos::Time::wallTime()-tcpu;
    }

    // blank residual DOFs which are on Dirichlet BC
    // We can do this because the values at the dirichlet positions
    // are not used anyway.
    // We could avoid this though, if velrowmap_ and prerowmap_ would
    // not include the dirichlet values as well. But it is expensive
    // to avoid that.
    fluiddbcmaps_->InsertCondVector(fluiddbcmaps_->ExtractCondVector(fluidzeros_), fluidresidual_);
    const Teuchos::RCP<const Epetra_Vector> xfluidzeros = LINALG::CreateVector(*xfluiddis_->DofRowMap(),true);
    xfluiddbcmaps_->InsertCondVector(xfluiddbcmaps_->ExtractCondVector(xfluidzeros), xfluidresidual_);

    // insert fluid and xfluid residuals to fluidxfluidresidual
    fluidxfluidsplitter_.InsertXFluidVector(xfluidresidual_,fluidxfluidresidual_);
    fluidxfluidsplitter_.InsertFluidVector(fluidresidual_,fluidxfluidresidual_);

    double vresnorm;
    double incvelnorm_L2;
    double velnorm_L2;

    Teuchos::RCP<Epetra_Vector> fluidxfluidonlyvel = fluidxfluidvelpressplitter_.ExtractOtherVector(fluidxfluidresidual_);
    fluidxfluidonlyvel->Norm2(&vresnorm);

    fluidxfluidvelpressplitter_.ExtractOtherVector(fluidxfluidincvel_,fluidxfluidonlyvel);
    fluidxfluidonlyvel->Norm2(&incvelnorm_L2);

    fluidxfluidvelpressplitter_.ExtractOtherVector(fluidxfluidstate_.velnp_,fluidxfluidonlyvel);
    fluidxfluidonlyvel->Norm2(&velnorm_L2);


    double incprenorm_L2;
    double presnorm;
    double prenorm_L2;

    Teuchos::RCP<Epetra_Vector> fluidxfluidonlypre = fluidxfluidvelpressplitter_.ExtractCondVector(fluidxfluidresidual_);
    fluidxfluidonlypre->Norm2(&presnorm);

    fluidxfluidvelpressplitter_.ExtractCondVector(fluidxfluidincvel_,fluidxfluidonlypre);
    fluidxfluidonlypre->Norm2(&incprenorm_L2);

    fluidxfluidvelpressplitter_.ExtractCondVector(fluidxfluidstate_.velnp_,fluidxfluidonlypre);
    fluidxfluidonlypre->Norm2(&prenorm_L2);


    double incfullnorm_L2;
    double fullnorm_L2;
    double fullresnorm;

    Epetra_Vector full(*fluidxfluidrowmap_);
    Epetra_Import importer(*fluidxfluidrowmap_,fluidxfluidresidual_->Map());

    int err = full.Import(*fluidxfluidresidual_,importer,Insert);
    if (err) dserror("Import using importer returned err=%d",err);
    full.Norm2(&fullresnorm);

    err = full.Import(*fluidxfluidincvel_,importer,Insert);
    if (err) dserror("Import using importer returned err=%d",err);
    full.Norm2(&incfullnorm_L2);

    err = full.Import(*fluidxfluidstate_.velnp_,importer,Insert);
    if (err) dserror("Import using importer returned err=%d",err);
    full.Norm2(&fullnorm_L2);

    // care for the case that nothing really happens in the velocity
    // or pressure field
    if (velnorm_L2 < 1e-5)  velnorm_L2 = 1.0;
    if (prenorm_L2 < 1e-5)  prenorm_L2 = 1.0;
    if (fullnorm_L2 < 1e-5) fullnorm_L2 = 1.0;

    //-------------------------------------------------- output to screen
    /* special case of very first iteration step:
        - solution increment is not yet available
        - convergence check is not required (we solve at least once!)    */
    if (itnum == 1)
    {
      if (myrank_ == 0)
      {
        printf("|  %3d/%3d   | %9.2E[L_2 ]  | %7.1e | %7.1e | %7.1e |    --   |    --   |    --   |",
                       itnum,itemax,ittol,vresnorm,presnorm,fullresnorm);
        printf(" (      --     ,te=%9.2E",dtele_);
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
      if (vresnorm <= ittol and presnorm <= ittol and fullresnorm <= ittol and
          incvelnorm_L2/velnorm_L2 <= ittol and incprenorm_L2/prenorm_L2 <= ittol and
          incfullnorm_L2/fullnorm_L2 <= ittol)
      {
        stopnonliniter=true;
        if (myrank_ == 0)
        {
          printf("|  %3d/%3d   | %9.2e[L_2 ]  | %7.1e | %7.1e | %7.1e | %7.1e | %7.1e | %7.1e |",
                  itnum,itemax,ittol,vresnorm,presnorm,fullresnorm,
                  incvelnorm_L2/velnorm_L2,incprenorm_L2/prenorm_L2,incfullnorm_L2/fullnorm_L2);
          printf(" (ts=%10.3E,te=%10.3E)\n",dtsolve_,dtele_);
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
          printf("|  %3d/%3d   | %9.2E[L_2 ]  | %7.1e | %7.1e | %7.1e | %7.1e | %7.1e | %7.1e |",
                  itnum,itemax,ittol,vresnorm,presnorm,fullresnorm,
                  incvelnorm_L2/velnorm_L2,incprenorm_L2/prenorm_L2,incfullnorm_L2/fullnorm_L2);
          printf(" (ts=%10.3E,te=%10.3E)",dtsolve_,dtele_);
        }
    }

    // warn if itemax is reached without convergence, but proceed to
    // next timestep...
    if ((itnum > itemax) and stopnonliniter == false)
    {
      stopnonliniter=true;
      if (myrank_ == 0)
      {
        printf("\n");
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
    fluidxfluidincvel_->PutScalar(0.0);

    // Add the fluid & xfluid & couple-matrices to fluidxfluidsysmat
    fluidxfluidsysmat_->Add(*fluidsysmat_,false,1.0,0.0);
    fluidxfluidsysmat_->Add(*xfluidsysmat_,false,1.0,1.0);
    fluidxfluidsysmat_->Add(*Cud_,false,1.0,1.0);
    fluidxfluidsysmat_->Add(*Cdu_,false,1.0,1.0);
    fluidxfluidsysmat_->Add(*Cdd_,false,1.0,1.0);
    fluidxfluidsysmat_->Complete();

    {
      // time measurement: application of dbc
      TEUCHOS_FUNC_TIME_MONITOR("      + apply DBC");
      //fluid-xfluid dbc
      Teuchos::RCP<Epetra_Vector> fluidxfluiddbc = LINALG::CreateVector(*fluidxfluidrowmap_,true);

      //build a merged map from fluid and xfluid dbc-maps
      std::vector<Teuchos::RCP<const Epetra_Map> > maps;
      maps.push_back(xfluiddbcmaps_->CondMap());
      maps.push_back(fluiddbcmaps_->CondMap());
      Teuchos::RCP<Epetra_Map> fluidxfluiddbcmaps = LINALG::MultiMapExtractor::MergeMaps(maps);

      const Teuchos::RCP<const Epetra_Vector> fluidxfluidzeros = LINALG::CreateVector(*fluidxfluidrowmap_,true);
      LINALG::ApplyDirichlettoSystem(fluidxfluidsysmat_,fluidxfluidincvel_,fluidxfluidresidual_,fluidxfluidzeros,*fluidxfluiddbcmaps);
    }
    if (myrank_ == 0)
    {
      cout << endl;
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

      fluidxfluidsolver->Solve(fluidxfluidsysmat_->EpetraOperator(),fluidxfluidincvel_,fluidxfluidresidual_,true,itnum==1, w_, c_, project_);
      fluidxfluidsolver->ResetTolerance();

      // end time measurement for solver
      dtsolve_ = Teuchos::Time::wallTime()-tcpusolve;
    }

    // -------------------------------------------------------------------
    // update velocity and pressure values by increments
    // -------------------------------------------------------------------
    fluidxfluidstate_.velnp_->Update(1.0,*fluidxfluidincvel_,1.0);
    xfluidstate_.velnp_  = fluidxfluidsplitter_.ExtractXFluidVector(fluidxfluidstate_.velnp_);
    fluidstate_.velnp_   = fluidxfluidsplitter_.ExtractFluidVector(fluidxfluidstate_.velnp_);

    // extract residual
    xfluidresidual_ = fluidxfluidsplitter_.ExtractXFluidVector(fluidxfluidresidual_);
    fluidresidual_  = fluidxfluidsplitter_.ExtractFluidVector(fluidxfluidresidual_);

    //------------------- store nodal increment for element stress update
    xfluidincvel_        = fluidxfluidsplitter_.ExtractXFluidVector(fluidxfluidincvel_);
    fluidincvel_         = fluidxfluidsplitter_.ExtractFluidVector(fluidxfluidincvel_);
    oldinc->Update(1.0,*xfluidincvel_,0.0);

    // Update the fluid material velocity along the interface (fxfivelnp_), source (in): fluidstate_.velnp_
    LINALG::Export(*fluidstate_.velnp_,*fxfivelnp_);
    LINALG::Export(*fluidincvel_,*fxfiincvel_);
    fluidxfluidboundarydis->SetState("ivelcolnp",fxfivelnp_);
  }

  CalculateAcceleration();

  //Gmsh Output
  OutputToGmsh(step_, time_);
  MovingFluidOutput();
  FluidXFluidBoundaryOutput(fluidxfluidboundarydis);

  xfluidresidual_ = Teuchos::null;
  fluidresidual_ = Teuchos::null;
  w_ = Teuchos::null;
  c_ = Teuchos::null;
  fluidxfluidsysmat_ = Teuchos::null;
  xfluidsysmat_ = Teuchos::null;
  fluidsysmat_ = Teuchos::null;

} // FluidImplicitTimeInt::NonlinearSolve

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::FluidXFluidImplicitTimeInt::CalculateAcceleration()
{

 // compute acceleration at n+1 for xfluid
  TIMEINT_THETA_BDF2::CalculateAcceleration(
    xfluidstate_.velnp_, xfluidstate_.veln_, xfluidstate_.velnm_, xfluidstate_.accn_,
    timealgo_, step_, theta_, dta_, dtp_,
    xfluidstate_.accnp_);

  // compute acceleration at n+1 for fluid
  Teuchos::RCP<Epetra_Vector> onlyaccn  = fluidvelpressplitter_.ExtractOtherVector(fluidstate_.accn_ );
  Teuchos::RCP<Epetra_Vector> onlyaccnp = fluidvelpressplitter_.ExtractOtherVector(fluidstate_.accnp_);
  Teuchos::RCP<Epetra_Vector> onlyvelnm = fluidvelpressplitter_.ExtractOtherVector(fluidstate_.velnm_);
  Teuchos::RCP<Epetra_Vector> onlyveln  = fluidvelpressplitter_.ExtractOtherVector(fluidstate_.veln_ );
  Teuchos::RCP<Epetra_Vector> onlyvelnp = fluidvelpressplitter_.ExtractOtherVector(fluidstate_.velnp_);

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

  return;
}//FLD::FluidXFluidImplicitTimeInt::CalculateAcceleratio
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
  cout << "FLD::FluidXFluidImplicitTimeInt::TimeUpdate()" << endl;
  ParameterList *  stabparams=&(params_.sublist("STABILIZATION"));

  if(stabparams->get<string>("TDS") == "time_dependent")
  {
    const double tcpu=Teuchos::Time::wallTime();

    if(myrank_==0)
    {
      cout << "time update for subscales";
    }

    // create the parameters for the discretization
    ParameterList eleparams;
    // action for elements
    eleparams.set("action","time update for subscales");

    // update time paramters
    if (timealgo_==INPAR::FLUID::timeint_one_step_theta)
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

    if(myrank_==0)
    {
      cout << "("<<Teuchos::Time::wallTime()-tcpu<<")\n";
    }
  }

  // update old acceleration
  fluidstate_.accn_->Update(1.0,*fluidstate_.accnp_,0.0);

  // velocities/pressures of this step become most recent
  // velocities/pressures of the last step
  fluidstate_.velnm_->Update(1.0,*fluidstate_.veln_ ,0.0);
  fluidstate_.veln_ ->Update(1.0,*fluidstate_.velnp_,0.0);

  if (alefluid_)
  {
    fluidstate_.dispnmm_->Update(1.0,*fluidstate_.dispnm_,0.0);
    fluidstate_.dispnm_->Update(1.0,*fluidstate_.dispn_,0.0);
    fluidstate_.dispn_ ->Update(1.0,*fluidstate_.dispnp_,0.0);
  }

  fluiddis_->ClearState();
  fluiddis_->SetState("velnp",fluidstate_.velnp_);
  fluiddis_->SetState("hist",hist_);

  // export the fluid mesh displacement to FluidXFluidBoundaryDis displacement
  if (alefluid_)
  {
    fxfidispn_->Update(1.0,*fxfidispnp_,0.0);
    LINALG::Export(*fluidstate_.dispnp_,*fxfidispnp_);
  }

  // Update the fluid material velocity along the interface at time step n and
  // n-1, source (in): fluidstate_.veln_ and fluidstate_.velnm_
  LINALG::Export(*fluidstate_.veln_,*fxfiveln_);
  LINALG::Export(*fluidstate_.velnm_,*fxfivelnm_);

  // xfluid
  // update old acceleration
  xfluidstate_.accn_->Update(1.0,*xfluidstate_.accnp_,0.0);

  // velocities/pressures of this step become most recent
  // velocities/pressures of the last step
  xfluidstate_.velnm_->Update(1.0,*xfluidstate_.veln_ ,0.0);
  xfluidstate_.veln_ ->Update(1.0,*xfluidstate_.velnp_,0.0);
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
  //LiftDrag();

  // -------------------------------------------------------------------
  //                        compute flow rates
  // -------------------------------------------------------------------
  //ComputeFlowRates();

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
  const bool write_visualization_data = step_%upres_ == 0;
  const bool write_restart_data       = uprestart_ != 0 and step_%uprestart_ == 0;

  // output of solution
  if (write_visualization_data)
  {
    // step number and time
    output_.NewStep(step_,time_);

    // velocity/pressure vector
    output_.WriteVector("velnp",fluidstate_.velnp_);

    // (hydrodynamic) pressure
    Teuchos::RCP<Epetra_Vector> pressure = fluidvelpressplitter_.ExtractCondVector(fluidstate_.velnp_);
    pressure->Scale(density_);
    output_.WriteVector("pressure", pressure);

    //output_.WriteVector("residual", trueresidual_);
    if (alefluid_) output_.WriteVector("dispnp", fluidstate_.dispnm_);

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

    if (write_restart_data) //add restart data
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

    //xfluid output
    {
      xfluidoutput_->NewStep(step_,time_);

      Teuchos::RCP<Epetra_Vector> xfluidvelnp_out = dofmanager_np_->fillPhysicalOutputVector(
             *xfluidstate_.velnp_, xfluiddofset_out_, xfluidstate_.nodalDofDistributionMap_, physprob_.fieldset_);

      xfluidoutput_->WriteVector("velocity_smoothed",xfluidvelnp_out);

      Teuchos::RCP<Epetra_Vector> xfluidpressure = xfluidvelpressplitterForOutput_.ExtractCondVector(xfluidvelnp_out);
      xfluidpressure->Scale(density_);
      xfluidoutput_->WriteVector("pressure_smoothed", xfluidpressure);

      if (step_==upres_) xfluidoutput_->WriteElementData();

      if (write_restart_data) //add restart data
      {
        xfluidoutput_->WriteVector("velnp", xfluidstate_.velnp_);
        xfluidoutput_->WriteVector("veln" , xfluidstate_.veln_);
        xfluidoutput_->WriteVector("velnm", xfluidstate_.velnm_);
        xfluidoutput_->WriteVector("accnp", xfluidstate_.accnp_);
        xfluidoutput_->WriteVector("accn" , xfluidstate_.accn_);
      }
    }

  }
  // write restart also when uprestart_ is not a integer multiple of upres_
  else if (write_restart_data)
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
  const Teuchos::RCP<DRT::Discretization> cutterdiscret,
  const INPAR::FLUID::InitialField initfield,
  const int startfuncno
  )
{
  // initial field by (undisturbed) function (init==2)
  // or disturbed function (init==3)
  // create zero displacement vector to use initial position of interface
  {
    const Epetra_Map* fluidxfluid_dofcolmap = cutterdiscret->DofColMap();

    cutterdiscret->SetState("idispcolnp", LINALG::CreateVector(*fluidxfluid_dofcolmap,true));
    cutterdiscret->SetState("ivelcolnp", LINALG::CreateVector(*fluidxfluid_dofcolmap,true));
    cutterdiscret->SetState("idispcoln", LINALG::CreateVector(*fluidxfluid_dofcolmap,true));
    cutterdiscret->SetState("ivelcoln", LINALG::CreateVector(*fluidxfluid_dofcolmap,true));
    cutterdiscret->SetState("ivelcolnm", LINALG::CreateVector(*fluidxfluid_dofcolmap,true));
    cutterdiscret->SetState("iacccoln", LINALG::CreateVector(*fluidxfluid_dofcolmap,true));

    ComputeInterfaceAndSetDOFs(cutterdiscret);
    cutterdiscret->ClearState();
  }
  if (initfield == INPAR::FLUID::initfield_field_by_function or
      initfield == INPAR::FLUID::initfield_disturbed_field_from_function)
  {
    const int numdim = params_.get<int>("number of velocity degrees of freedom");

    // loop all xfluid nodes on the processor
    for(int lnodeid=0;lnodeid<xfluiddis_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = xfluiddis_->lRowNode(lnodeid);
      // the set of degrees of freedom associated with the node
      const vector<int> nodedofset = xfluiddis_->Dof(lnode);
      if (nodedofset.size() > 0)
      {
        for(int index=0;index<numdim+1;++index)
        {
          int gid = nodedofset[index];
          double initialval=DRT::Problem::Instance()->Funct(startfuncno-1).Evaluate(index,lnode->X(),0.0,NULL);
          xfluidstate_.velnp_->ReplaceGlobalValues(1,&initialval,&gid);
          xfluidstate_.veln_ ->ReplaceGlobalValues(1,&initialval,&gid);
        }
      }
    }

    // loop all fluid nodes on the processor
    for(int flnodeid=0;flnodeid<fluiddis_->NumMyRowNodes();flnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = fluiddis_->lRowNode(flnodeid);
      // the set of degrees of freedom associated with the node
      const vector<int> nodedofset = fluiddis_->Dof(lnode);

      for(int index=0;index<numdim+1;++index)
      {
        int gid = nodedofset[index];
        double initialval=DRT::Problem::Instance()->Funct(startfuncno-1).Evaluate(index,lnode->X(),0.0,NULL);
        fluidstate_.velnp_->ReplaceGlobalValues(1,&initialval,&gid);
        fluidstate_.veln_ ->ReplaceGlobalValues(1,&initialval,&gid);
      }
    }
  }

  // add random perturbation of certain percentage to function
  else if (initfield == INPAR::FLUID::initfield_disturbed_field_from_function)
  {
    dserror ("initfield_disturbed_field_from_function not adapted for fluidxfluid");
    const Epetra_Map* dofrowmap = fluiddis_->DofRowMap();

    int err = 0;

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

      // check whether a periodic boundary condition is active on this nodedofset
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

  // special initial function: two counter-rotating vortices (2-D) and flame front
  // for flame-vortex interaction problem
  else if (initfield == INPAR::FLUID::initfield_flame_vortex_interaction)
  {
    dserror ("INPAR::FLUID::initfield_flame_vortex_interaction not adapted for fluidxfluid");
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
    dserror ("INPAR::FLUID::initfield_beltrami_flow not adapted for fluidxfluid");
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
    dserror ("INPAR::FLUID::initfield_bochev_test not adapted for fluidxfluid");
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
  const Epetra_Map* fxfboundary_dofcolmap = fluidxfluidboundarydis->DofColMap();
  fluidxfluidboundarydis->SetState("idispcolnp", LINALG::CreateVector(*fxfboundary_dofcolmap,true));
  fluidxfluidboundarydis->SetState("ivelcolnp", LINALG::CreateVector(*fxfboundary_dofcolmap,true));
  fluidxfluidboundarydis->SetState("idispcoln", LINALG::CreateVector(*fxfboundary_dofcolmap,true));
  fluidxfluidboundarydis->SetState("ivelcoln", LINALG::CreateVector(*fxfboundary_dofcolmap,true));
  fluidxfluidboundarydis->SetState("ivelcolnm", LINALG::CreateVector(*fxfboundary_dofcolmap,true));
  fluidxfluidboundarydis->SetState("iacccoln", LINALG::CreateVector(*fxfboundary_dofcolmap,true));

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
   //  Set time parameter for element call
   // -------------------------------------------------------------------
   SetElementTimeParameter();


    // -------------------------------------------------------------------
    //                     solve nonlinear equation system
    // -------------------------------------------------------------------
    NonlinearSolve(fluidxfluidboundarydis);

    // -------------------------------------------------------------------
    //         calculate lift'n'drag forces from the residual
    // -------------------------------------------------------------------
    //LiftDrag();

    // -------------------------------------------------------------------
    //                        compute flow rates
    // -------------------------------------------------------------------
    //ComputeFlowRates();

    // -------------------------------------------------------------------
    //                         output of solution
    // -------------------------------------------------------------------
    Output();

    {
      // compute acceleration at n+1 for xfluid
      TIMEINT_THETA_BDF2::CalculateAcceleration(
          xfluidstate_.velnp_, xfluidstate_.veln_, xfluidstate_.velnm_, xfluidstate_.accn_,
          timealgo_, step_, theta_, dta_, dtp_,
          xfluidstate_.accnp_);

      // compute acceleration at n+1 for fluid
      Teuchos::RCP<Epetra_Vector> onlyaccn  = fluidvelpressplitter_.ExtractOtherVector(fluidstate_.accn_ );
      Teuchos::RCP<Epetra_Vector> onlyaccnp = fluidvelpressplitter_.ExtractOtherVector(fluidstate_.accnp_);
      Teuchos::RCP<Epetra_Vector> onlyvelnm = fluidvelpressplitter_.ExtractOtherVector(fluidstate_.velnm_);
      Teuchos::RCP<Epetra_Vector> onlyveln  = fluidvelpressplitter_.ExtractOtherVector(fluidstate_.veln_ );
      Teuchos::RCP<Epetra_Vector> onlyvelnp = fluidvelpressplitter_.ExtractOtherVector(fluidstate_.velnp_);

      TIMEINT_THETA_BDF2::CalculateAcceleration(onlyvelnp, onlyveln, onlyvelnm, onlyaccn ,
                                                  timealgo_, step_ , theta_, dta_,
                                                  dtp_, onlyaccnp);

      // copy back into global vector
      LINALG::Export(*onlyaccnp,*fluidstate_.accnp_);
    }

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
      (*integratedshapefunc)[i] = (*fluidtrueresidual_)[i]/(*integratedshapefunc)[i];
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

  FLD::UTILS::LiftDrag(*fluiddis_,*fluidtrueresidual_,params_,liftdragvals);

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
  condmaps.push_back(fluiddbcmaps_->CondMap());
  Teuchos::RCP<Epetra_Map> condmerged = LINALG::MultiMapExtractor::MergeMaps(condmaps);
  *fluiddbcmaps_ = LINALG::MapExtractor(*(fluiddis_->DofRowMap()), condmerged);
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> FLD::FluidXFluidImplicitTimeInt::VelocityRowMap()
{ return fluidvelpressplitter_.OtherMap(); }

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> FLD::FluidXFluidImplicitTimeInt::PressureRowMap()
{ return fluidvelpressplitter_.CondMap(); }


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
  cout << "FLD::FluidXFluidImplicitTimeInt::ComputeInterfaceAndSetDOFs" << endl;
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

  //dummy vectors  here not needed
  const vector<int> MovingFluideleGIDs;
  const vector<int> MovingFluidNodeGIDs;

  // map of element id of the fluid-patch element to it's patchbox
  std::map<int,GEO::CUT::BoundingBox> patchboxes;
  if (Step() > 1 and alefluid_)
    CreatePatchBoxes(patchboxes);

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

    // save the label of each fluid-node in xfluidstate_.nodelabelnp
    // (label=0 is fluid, label=1 is structure)
    for (int node=0; node<xfluiddis_->NumMyColNodes(); ++node)
    {
      const DRT::Node* actnode = xfluiddis_->lColNode(node);
      const LINALG::Matrix<3,1> pos(actnode->X());
      xfluidstate_.nodelabelnp_[actnode->Id()] = ih_np_->PositionWithinConditionNP(pos);
    }

    // create switcher
    const XFEM::DofDistributionSwitcher dofswitch(
            ih_np_, dofmanager_np_,
            olddofrowmap, newdofrowmap,
            oldNodalDofDistributionMap,xfluidstate_.nodalDofDistributionMap_,
            oldElementalDofDistributionMap,xfluidstate_.elementalDofDistributionMap_,
            xfluidstate_.nodelabeln_);

    // --------------------------------------------
    // switch state vectors to new dof distribution
    // --------------------------------------------

   // cout0_ << " ->  Initialize system vectors..." << endl;
    // accelerations
    dofswitch.mapVectorToNewDofDistributionFluidXFluid(xfluidstate_.accnp_,xfluidstate_.nodelabeln_);
    dofswitch.mapVectorToNewDofDistributionFluidXFluid(xfluidstate_.accn_,xfluidstate_.nodelabeln_);

    // velocities and pressures
    dofswitch.mapVectorToNewDofDistributionFluidXFluid(xfluidstate_.velnp_,xfluidstate_.nodelabeln_);
    dofswitch.mapVectorToNewDofDistributionFluidXFluid(xfluidstate_.veln_,xfluidstate_.nodelabeln_);
    dofswitch.mapVectorToNewDofDistributionFluidXFluid(xfluidstate_.velnm_,xfluidstate_.nodelabeln_);

    // if unknown dofs appear, project the value of fluid patch to background fluid
    if (Step() > 1 and alefluid_)
    {
      //Gmsh output before Interpolation
      PlotVectorFieldToGmsh(DRT::UTILS::GetColVersionOfRowVector(xfluiddis_, xfluidstate_.veln_) , "sol_xfluid_veln_beforeProject" ,"Velocity Solution (Physical) n" ,false, Step(), Time());
      PlotVectorFieldToGmsh(DRT::UTILS::GetColVersionOfRowVector(xfluiddis_, xfluidstate_.velnm_) , "sol_xfluid_velnm_beforeProject" ,"Velocity Solution (Physical) n-1" ,false, Step(), Time());
      PlotVectorFieldToGmsh(DRT::UTILS::GetColVersionOfRowVector(xfluiddis_, xfluidstate_.accn_) , "sol_xfluid_accn_beforeProject" ,"Acc Solution (Physical) n-1" ,false, Step(), Time());

      //Project the values of the last tme step from embedded fluid to background fluid
      dofswitch.projectEmbeddedtoBackgroundfluid(fluiddis_, patchboxes, fluidstate_.veln_,xfluidstate_.veln_, fluidstate_.dispnm_, step_);
      dofswitch.projectEmbeddedtoBackgroundfluid(fluiddis_, patchboxes, fluidstate_.velnm_,xfluidstate_.velnm_, fluidstate_.dispnm_, step_);
      dofswitch.projectEmbeddedtoBackgroundfluid(fluiddis_, patchboxes, fluidstate_.accn_,xfluidstate_.accn_,fluidstate_.dispnm_, step_);

      //Gmsh output after Interpolation
      PlotVectorFieldToGmsh(DRT::UTILS::GetColVersionOfRowVector(xfluiddis_, xfluidstate_.veln_) , "sol_xfluid_veln_afterProject" ,"Velocity Solution (Physical) n"    ,false, Step(), Time());
      PlotVectorFieldToGmsh(DRT::UTILS::GetColVersionOfRowVector(xfluiddis_, xfluidstate_.velnm_) , "sol_xfluid_velnm_afterProject" ,"Velocity Solution (Physical) n-1"    ,false, Step(), Time());
      PlotVectorFieldToGmsh(DRT::UTILS::GetColVersionOfRowVector(xfluiddis_, xfluidstate_.accn_) , "sol_xfluid_accn_afterProject" ,"Acc Solution (Physical) n-1"    ,false, Step(), Time());
    }
  }

  // --------------------------------------------------
  // create remaining vectors with new dof distribution
  // --------------------------------------------------
  //hist_         = LINALG::CreateVector(newdofrowmap,true);

  // object holds maps/subsets for DOFs subjected to Dirichlet BCs and otherwise
  xfluiddbcmaps_ = Teuchos::rcp(new LINALG::MapExtractor());
  {
    ParameterList eleparams;
    // other parameters needed by the elements
    eleparams.set("total time",time_);
    Teuchos::RCP<Epetra_Vector> tmp = LINALG::CreateVector(newdofrowmap,true);
    xfluiddis_->EvaluateDirichletXFEM(eleparams, tmp, Teuchos::null, Teuchos::null,
                                Teuchos::null, xfluiddbcmaps_);
  }

  // -------------------------------------------------------------------
  // get a vector layout from the discretization for a vector which only
  // contains the velocity dofs and for one vector which only contains
  // pressure degrees of freedom.
  // -------------------------------------------------------------------
  FLD::UTILS::SetupXFluidSplit(*xfluiddis_,dofmanager_np_,xfluidvelpressplitter_);

//  // project old interpolated velocity vector onto divergence free space
//  if (xparams_.get<bool>("INCOMP_PROJECTION"))
//    if (timealgo_ != INPAR::FLUID::timeint_stationary)
//      ProjectOldTimeStepValues(velpressplitter_);

  xfluid_neumann_loads_= LINALG::CreateVector(newdofrowmap,true);

  // ---------------------------------
  // Vectors used for solution process
  // ---------------------------------
  xfluidresidual_     = LINALG::CreateVector(newdofrowmap,true);
  xfluidtrueresidual_ = LINALG::CreateVector(newdofrowmap,true);
  xfluidincvel_       = LINALG::CreateVector(newdofrowmap,true);

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
        vector<int> lmstride;
        actele->LocationVector(*xfluiddis_, lm, lmowner, lmstride);

        // extract local values from the global vector
        vector<double> myvelnp(lm.size());
        DRT::UTILS::ExtractMyValues(*vectorfield, myvelnp, lm);

        const int numparam = eledofman.NumDofPerField(XFEM::PHYSICS::Velx);
        LINALG::SerialDenseMatrix elementvalues(3, numparam);

        if ( numparam > 0 )
        {
          const vector<int>& dofposvelx =
            eledofman.LocalDofPosPerField(XFEM::PHYSICS::Velx);
          const vector<int>& dofposvely =
            eledofman.LocalDofPosPerField(XFEM::PHYSICS::Vely);
          const vector<int>& dofposvelz =
            eledofman.LocalDofPosPerField(XFEM::PHYSICS::Velz);

          for (int iparam=0; iparam<numparam; ++iparam)
          {
            elementvalues(0, iparam) = myvelnp[dofposvelx[iparam]];
            elementvalues(1, iparam) = myvelnp[dofposvely[iparam]];
            elementvalues(2, iparam) = myvelnp[dofposvelz[iparam]];
          }
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
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FLD::FluidXFluidImplicitTimeInt::PrepareFluidXFluidBoundaryDis(const Teuchos::RCP<DRT::Discretization>  fluidxfluidboundarydis)
{
  // put vectors into boundary discretization (SetState generates col vector automatically)
  fluidxfluidboundarydis->SetState("idispcolnp",fxfidispnp_);
  fluidxfluidboundarydis->SetState("idispcoln" ,fxfidispn_);
  fluidxfluidboundarydis->SetState("ivelcolnp" ,fxfivelnp_);
  fluidxfluidboundarydis->SetState("ivelcoln"  ,fxfiveln_);
  fluidxfluidboundarydis->SetState("ivelcolnm" ,fxfivelnm_);
  fluidxfluidboundarydis->SetState("iacccoln"  ,fxfiaccn_);
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FLD::FluidXFluidImplicitTimeInt::PrepareFluidXFluidBoundaryDofset(const Teuchos::RCP<DRT::Discretization>  fluidxfluidboundarydis)
{
  ComputeInterfaceAndSetDOFs(fluidxfluidboundarydis);

  RCP<DRT::DofSet> newdofset = rcp(new DRT::TransparentDofSet(fluiddis_));
  fluidxfluidboundarydis->ReplaceDofSet(newdofset);
  const int error = fluidxfluidboundarydis->FillComplete();
  if (error) dserror("fluidxfluidboundarydis->FillComplete() returned error=%d",error);

  DRT::UTILS::PrintParallelDistribution(*fluidxfluidboundarydis);

  // mesh output needed for post processing in paraview
  //FluidFluidboundaryoutput_ = rcp(new IO::DiscretizationWriter(FluidFluidboundarydis_));

  // create fluid-fluid interface DOF vectors
  fxfidispnp_    = LINALG::CreateVector(*fluidxfluidboundarydis->DofRowMap(),true);
  fxfivelnp_     = LINALG::CreateVector(*fluidxfluidboundarydis->DofRowMap(),true);
  fxfitrueresnp_ = LINALG::CreateVector(*fluidxfluidboundarydis->DofRowMap(),true);
  fxfidispn_     = LINALG::CreateVector(*fluidxfluidboundarydis->DofRowMap(),true);
  fxfiveln_      = LINALG::CreateVector(*fluidxfluidboundarydis->DofRowMap(),true);
  fxfivelnm_     = LINALG::CreateVector(*fluidxfluidboundarydis->DofRowMap(),true);
  fxfiaccnp_     = LINALG::CreateVector(*fluidxfluidboundarydis->DofRowMap(),true);
  fxfiaccn_      = LINALG::CreateVector(*fluidxfluidboundarydis->DofRowMap(),true);
  fxfiincvel_    = LINALG::CreateVector(*fluidxfluidboundarydis->DofRowMap(),true);
}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FLD::FluidXFluidImplicitTimeInt::ComputeFluidXFluidInterfaceAccelerationsAndVelocities()
{
  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();
  const double dt = fsidyn.get<double>("TIMESTEP");

  // use only one of these switches
  //#define USE_FSI_VEL_AND_COMPUTE_ONLY_ACCNP
  #define COMPUTE_VEL_AND_ACCNP_VIA_BETA_NEWMARK
  //#define COMPUTE_VEL_AND_ACCNP_VIA_ONE_STEP_THETA

  #ifdef USE_FSI_VEL_AND_COMPUTE_ONLY_ACCNP
    double theta = 0.66;
    if (Step() == 1)
      theta = 1.0;
    iaccnp_->Update(-(1.0-theta)/(theta),*iaccn_,0.0);
    iaccnp_->Update(1.0/(theta*dt),*ivelnp_,-1.0/(theta*dt),*iveln_,1.0);
  #endif
  #ifdef COMPUTE_VEL_AND_ACCNP_VIA_BETA_NEWMARK
    double beta;
    double gamma;
    if (DRT::INPUT::IntegralValue<int>(fsidyn,"SECONDORDER") == 1)
    {
      if (Step() == 1)
      {
        gamma = 1.0;
        beta = gamma/2.0;
      }
      else
      {
        gamma = 0.66;
        beta = gamma/2.0;
      }
    }
    else
    {
      gamma = 1.0;
      beta = gamma/2.0;
    }

    // compute acceleration at timestep n+1
    fxfiaccnp_->Update(-(1.0-(2.0*beta))/(2.0*beta),*fxfiaccn_,0.0);
    fxfiaccnp_->Update(-1.0/(beta*dt),*fxfiveln_,1.0);
    fxfiaccnp_->Update(1.0/(beta*dt*dt),*fxfidispnp_,-1.0/(beta*dt*dt),*fxfidispn_,1.0);

    //NOTE: For fluid-fluid Coupling the interface velocity is updated at the end of while-loop in Nonlinearsolve.

  #endif
  #ifdef COMPUTE_VEL_AND_ACCNP_VIA_ONE_STEP_THETA

    double theta_vel = 0.5;
    double theta_acc = 0.66;
    if (DRT::INPUT::IntegralValue<int>(fsidyn,"SECONDORDER") == 1)
    {
      if (Step() == 1)
      {
        theta_vel = 1.0;
        theta_acc = 1.0;
      }
      else
      {
        theta_vel = 0.66;
        theta_acc = 0.66;
      }
    }
    else
    {
      theta_vel = 1.0;
      theta_acc = 1.0;
    }

    // For fluid-fluid Coupling the interface velocity is updated at the end of while-loop in Nonlinearsolve.
    //fxfivelnp_->Update(-(1.0-theta_vel)/theta_vel,*fluidfluidstate_.fiveln_,0.0);
    //fxfivelnp_->Update(1.0/(theta_vel*dt),*fluidfluidstate_.fidispnp_,-1.0/(theta_vel*dt),*fxfidispn_,1.0);

    // compute acceleration at timestep n+1
    fxfiaccnp_->Update(-(1.0-theta_acc)/theta_acc,*fxfiaccn_,0.0);
    fxfiaccnp_->Update(1.0/(theta_acc*dt),*fxfivelnp_,-1.0/(theta_acc*dt),*fxfiveln_,1.0);
  #endif
}
/*-----------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------*/
static void AdaptElementMatrix(
    const int&                dim1,
    const int&                dim2,
    Epetra_SerialDenseMatrix& elematrix
    )
{
  if (elematrix.M()!=dim1 or elematrix.N()!=dim2)
  {
    elematrix.Shape(dim1,dim2);
  }
  else
    memset(elematrix.A(),0,dim1*dim2*sizeof(double));
}

/*-----------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------*/
void FLD::FluidXFluidImplicitTimeInt::MonolithicMultiDisEvaluate(
    Teuchos::RCP<XFEM::InterfaceHandleXFSI> ih,
    const Teuchos::RCP<DRT::Discretization> xfluiddis,
    const Teuchos::RCP<DRT::Discretization> ifacedis,
    Teuchos::ParameterList&                 params,
    Teuchos::RCP<LINALG::SparseOperator>    Cud,
    Teuchos::RCP<LINALG::SparseOperator>    Cdu,
    Teuchos::RCP<LINALG::SparseOperator>    Cdd,
    Teuchos::RCP<Epetra_Vector>             rhsd,
    Teuchos::RCP<LINALG::SparseOperator>    systemmatrix1,
    Teuchos::RCP<Epetra_Vector>             systemvector1,
    const bool assemble_coupling_matrices,
    const bool fluidfluidcoupling)
{
  TEUCHOS_FUNC_TIME_MONITOR("FLD::XFluidImplicitTimeInt::MonolithicMultiDisEvaluate");

  if (!xfluiddis->Filled()) dserror("xfluiddis->FillComplete() was not called");
  if (!xfluiddis->HaveDofs()) dserror("xfluiddis->AssignDegreesOfFreedom() was not called");

  if (!xfluiddis->Filled()) dserror("xfluiddis->FillComplete() was not called");
  if (!xfluiddis->HaveDofs()) dserror("xfluiddis->AssignDegreesOfFreedom() was not called");

  systemmatrix1->Zero();

  // define element matrices and vectors
  Epetra_SerialDenseMatrix elematrix1;
  Epetra_SerialDenseMatrix elematrix2;
  Epetra_SerialDenseVector elevector1;
  Epetra_SerialDenseVector elevector2;
  Epetra_SerialDenseVector elevector3;

  // loop over column elements
  const int numcolele = xfluiddis->NumMyColElements();
  for (int i=0; i<numcolele; ++i)
  {
    DRT::Element* xele = xfluiddis->lColElement(i);
    const bool intersected = ih->ElementIntersected(xele->Id());

    vector<int> fluidlm;
    vector<int> fluidlmowner;
    vector<int> fluidlmstride;
    xele->LocationVector(*xfluiddis, fluidlm, fluidlmowner, fluidlmstride);

    // get dimension of element matrices and vectors
    // Reshape element matrices and vectors and init to zero
    const std::size_t fluideledim = fluidlm.size();
    AdaptElementMatrix(fluideledim, fluideledim, elematrix1);
    if (elevector1.Length()!=(int)fluideledim)
      elevector1.Size(fluideledim);
    else
      memset(elevector1.Values(),0,fluideledim*sizeof(double));


    const RCP<vector<int> > ifacepatchlm = rcp(new vector<int>());
    const RCP<vector<int> > ifacepatchlmowner = rcp(new vector<int>());
    if (intersected)
    {
       ih->GetInterfacepatchLocationVectors(*xele, ifacepatchlm, ifacepatchlmowner);
    }
    params.set("ifacepatchlm",ifacepatchlm);
    params.set("ifacepatchlmowner",ifacepatchlmowner);

    {
      TEUCHOS_FUNC_TIME_MONITOR("FLD::XFluidImplicitTimeInt::MonolithicMultiDisEvaluate elements");

      // call the element evaluate method
      const int err = xele->Evaluate(params,*xfluiddis,fluidlm,elematrix1,elematrix2,
                      elevector1,elevector2,elevector3);
      if (err) dserror("Proc %d: Element %d returned err=%d",xfluiddis->Comm().MyPID(),xele->Id(),err);
    }

    {
      TEUCHOS_FUNC_TIME_MONITOR("FLD::XFluidImplicitTimeInt::MonolithicMultiDisEvaluate assemble");

      if (intersected and assemble_coupling_matrices)
      {
        // create a dummy stride vector that is correct
        vector<int> stride(1); stride[0] = (int)ifacepatchlm->size();
        Cud->Assemble(-1, stride, *params.get<RCP<Epetra_SerialDenseMatrix> >("Cud"), fluidlm      , fluidlmowner      , *ifacepatchlm);
        stride[0] = (int)fluidlm.size();
        Cdu->Assemble(-1, stride, *params.get<RCP<Epetra_SerialDenseMatrix> >("Cdu"), *ifacepatchlm, *ifacepatchlmowner, fluidlm);
        stride[0] = (int)ifacepatchlm->size();
        Cdd->Assemble(-1, stride, *params.get<RCP<Epetra_SerialDenseMatrix> >("Cdd"), *ifacepatchlm, *ifacepatchlmowner, *ifacepatchlm);
        LINALG::Assemble(*rhsd, *params.get<RCP<Epetra_SerialDenseVector> >("rhsd"), *ifacepatchlm, *ifacepatchlmowner);
      }

      // create a dummy stride vector that is correct
      vector<int> stride(1); stride[0] = (int)fluidlm.size();
      systemmatrix1->Assemble(-1,stride,elematrix1,fluidlm,fluidlmowner);
      LINALG::Assemble(*systemvector1,elevector1,fluidlm,fluidlmowner);
    }

  } // element loop
  return;
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
void FLD::FluidXFluidImplicitTimeInt::OutputToGmsh(
    const int step,
    const double time
    ) const
{
  const Teuchos::ParameterList& xfemparams = DRT::Problem::Instance()->XFEMGeneralParams();
  const bool gmshdebugout = getIntegralValue<int>(xfemparams,"GMSH_DEBUG_OUT")==1;
  const bool screen_out = getIntegralValue<int>(xfemparams,"GMSH_DEBUG_OUT_SCREEN")==1;

  // get a copy on column parallel distribution
  Teuchos::RCP<const Epetra_Vector> output_col_velnp = DRT::UTILS::GetColVersionOfRowVector(xfluiddis_, xfluidstate_.velnp_);

  if (gmshdebugout and (this->physprob_.fieldset_.find(XFEM::PHYSICS::Pres) != this->physprob_.fieldset_.end()))
  {
    cout << "FluidXFluidImplicitTimeInt::OutputToGmsh()" << endl;

    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("sol_field_xfluid_pres", step, 5, screen_out, xfluiddis_->Comm().MyPID());
    std::ofstream gmshfilecontent(filename.c_str());

    const XFEM::PHYSICS::Field field = XFEM::PHYSICS::Pres;

    {
      gmshfilecontent << "View \" " << "Pressure Solution (Physical) \" {\n";
      for (int i=0; i<xfluiddis_->NumMyColElements(); ++i)
      {
        const DRT::Element* actele = xfluiddis_->lColElement(i);

        // create local copy of information about dofs
        const XFEM::ElementDofManager eledofman(*actele,physprob_.elementAnsatz_->getElementAnsatz(actele->Shape()),*dofmanager_np_);

        vector<int> lm;
        vector<int> lmowner;
        vector<int> lmstride;
        actele->LocationVector(*(xfluiddis_), lm, lmowner, lmstride);

        // extract local values from the global vector
        vector<double> myvelnp(lm.size());
        DRT::UTILS::ExtractMyValues(*output_col_velnp, myvelnp, lm);

        const int numparam = eledofman.NumDofPerField(field);
        if (numparam > 0)
        {
          const vector<int>& dofpos = eledofman.LocalDofPosPerField(field);

          LINALG::SerialDenseVector elementvalues(numparam);
          for (int iparam=0; iparam<numparam; ++iparam)
            elementvalues(iparam) = myvelnp[dofpos[iparam]];

          const GEO::DomainIntCells& domainintcells =
              dofmanager_np_->getInterfaceHandle()->GetDomainIntCells(actele);
          for (GEO::DomainIntCells::const_iterator cell =
              domainintcells.begin(); cell != domainintcells.end(); ++cell)
          {
            LINALG::SerialDenseVector cellvalues(DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
            XFEM::computeScalarCellNodeValuesFromNodalUnknowns(*actele, &*dofmanager_np_->getInterfaceHandle(), eledofman,
                *cell, field, elementvalues, cellvalues);
            IO::GMSH::cellWithScalarFieldToStream(
                cell->Shape(), cellvalues, cell->CellNodalPosXYZ(), gmshfilecontent);
          }
        }
      }
      gmshfilecontent << "};\n";
    }
    gmshfilecontent.close();
    if (screen_out) std::cout << " done" << endl;
  }

#if 1
  if (gmshdebugout)
  {
    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("solution_field_sigma_disc", step, 5, screen_out, xfluiddis_->Comm().MyPID());
    if (screen_out) std::cout << endl;
    const std::string filenamexx = IO::GMSH::GetNewFileNameAndDeleteOldFiles("solution_field_sigmaxx_disc", step, 5, screen_out, xfluiddis_->Comm().MyPID());
    if (screen_out) std::cout << endl;
    const std::string filenameyy = IO::GMSH::GetNewFileNameAndDeleteOldFiles("solution_field_sigmayy_disc", step, 5, screen_out, xfluiddis_->Comm().MyPID());
    if (screen_out) std::cout << endl;
    const std::string filenamezz = IO::GMSH::GetNewFileNameAndDeleteOldFiles("solution_field_sigmazz_disc", step, 5, screen_out, xfluiddis_->Comm().MyPID());
    if (screen_out) std::cout << endl;
    const std::string filenamexy = IO::GMSH::GetNewFileNameAndDeleteOldFiles("solution_field_sigmaxy_disc", step, 5, screen_out, xfluiddis_->Comm().MyPID());
    if (screen_out) std::cout << endl;
    const std::string filenamexz = IO::GMSH::GetNewFileNameAndDeleteOldFiles("solution_field_sigmaxz_disc", step, 5, screen_out, xfluiddis_->Comm().MyPID());
    if (screen_out) std::cout << endl;
    const std::string filenameyz = IO::GMSH::GetNewFileNameAndDeleteOldFiles("solution_field_sigmayz_disc", step, 5, screen_out, xfluiddis_->Comm().MyPID());
    std::ofstream gmshfilecontent(  filename.c_str());
    std::ofstream gmshfilecontentxx(filenamexx.c_str());
    std::ofstream gmshfilecontentyy(filenameyy.c_str());
    std::ofstream gmshfilecontentzz(filenamezz.c_str());
    std::ofstream gmshfilecontentxy(filenamexy.c_str());
    std::ofstream gmshfilecontentxz(filenamexz.c_str());
    std::ofstream gmshfilecontentyz(filenameyz.c_str());

    XFEM::PHYSICS::Field fieldxx = XFEM::PHYSICS::Sigmaxx;
    XFEM::PHYSICS::Field fieldyy = XFEM::PHYSICS::Sigmayy;
    XFEM::PHYSICS::Field fieldzz = XFEM::PHYSICS::Sigmazz;
    XFEM::PHYSICS::Field fieldxy = XFEM::PHYSICS::Sigmaxy;
    XFEM::PHYSICS::Field fieldxz = XFEM::PHYSICS::Sigmaxz;
    XFEM::PHYSICS::Field fieldyz = XFEM::PHYSICS::Sigmayz;
    switch (xparams_.get<INPAR::XFEM::BoundaryIntegralType>("EMBEDDED_BOUNDARY"))
    {
    case INPAR::XFEM::BoundaryTypeSigma:
      fieldxx = XFEM::PHYSICS::Sigmaxx;
      fieldyy = XFEM::PHYSICS::Sigmayy;
      fieldzz = XFEM::PHYSICS::Sigmazz;
      fieldxy = XFEM::PHYSICS::Sigmaxy;
      fieldxz = XFEM::PHYSICS::Sigmaxz;
      fieldyz = XFEM::PHYSICS::Sigmayz;
      break;
    case INPAR::XFEM::BoundaryTypeTauPressure:
      fieldxx = XFEM::PHYSICS::Tauxx;
      fieldyy = XFEM::PHYSICS::Tauyy;
      fieldzz = XFEM::PHYSICS::Tauzz;
      fieldxy = XFEM::PHYSICS::Tauxy;
      fieldxz = XFEM::PHYSICS::Tauxz;
      fieldyz = XFEM::PHYSICS::Tauyz;
      break;
    default:
      dserror("unknown boundary type");
    }

    std::size_t numplot = 0;
    {
      gmshfilecontent   << "View \" " << "Discontinous Stress Solution (Physical) \" {\n";
      gmshfilecontentxx << "View \" " << "Discontinous Stress (xx) Solution (Physical) \" {\n";
      gmshfilecontentyy << "View \" " << "Discontinous Stress (yy) Solution (Physical) \" {\n";
      gmshfilecontentzz << "View \" " << "Discontinous Stress (zz) Solution (Physical) \" {\n";
      gmshfilecontentxy << "View \" " << "Discontinous Stress (xy) Solution (Physical) \" {\n";
      gmshfilecontentxz << "View \" " << "Discontinous Stress (xz) Solution (Physical) \" {\n";
      gmshfilecontentyz << "View \" " << "Discontinous Stress (yz) Solution (Physical) \" {\n";
      for (int i=0; i<xfluiddis_->NumMyColElements(); ++i)
      {
        const DRT::Element* actele = xfluiddis_->lColElement(i);

        // create local copy of information about dofs
        const XFEM::ElementDofManager eledofman(*actele,physprob_.elementAnsatz_->getElementAnsatz(actele->Shape()),*dofmanager_np_);

        vector<int> lm;
        vector<int> lmowner;
        vector<int> lmstride;
        actele->LocationVector(*(xfluiddis_), lm, lmowner, lmstride);

        // extract local values from the global vector
        vector<double> myvelnp(lm.size());
        DRT::UTILS::ExtractMyValues(*output_col_velnp, myvelnp, lm);

        const int numparam = eledofman.NumDofPerField(fieldxx);
        if (numparam == 0)
          continue;

        numplot++;

        const vector<int>& dofposxx = eledofman.LocalDofPosPerField(fieldxx);
        const vector<int>& dofposyy = eledofman.LocalDofPosPerField(fieldyy);
        const vector<int>& dofposzz = eledofman.LocalDofPosPerField(fieldzz);
        const vector<int>& dofposxy = eledofman.LocalDofPosPerField(fieldxy);
        const vector<int>& dofposxz = eledofman.LocalDofPosPerField(fieldxz);
        const vector<int>& dofposyz = eledofman.LocalDofPosPerField(fieldyz);

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
          dofmanager_np_->getInterfaceHandle()->GetDomainIntCells(actele);
        for (GEO::DomainIntCells::const_iterator cell =
          domainintcells.begin(); cell != domainintcells.end(); ++cell)
        {
          //LINALG::SerialDenseMatrix xyze_cell(3, cell->NumNode());
          //cell->NodalPosXYZ(*actele, xyze_cell);
          // TODO remove
          const LINALG::SerialDenseMatrix& xyze_cell = cell->CellNodalPosXYZ();

          {
          LINALG::SerialDenseMatrix cellvalues(9,DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
          XFEM::computeTensorCellNodeValuesFromElementUnknowns(*actele, &*dofmanager_np_->getInterfaceHandle(), eledofman,
              *cell, fieldxx, elementvalues, cellvalues);
           IO::GMSH::cellWithTensorFieldToStream(cell->Shape(), cellvalues, xyze_cell, gmshfilecontent);
          }

          {
          LINALG::SerialDenseVector cellvaluexx(DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
          XFEM::computeScalarCellNodeValuesFromElementUnknowns(*actele, &*dofmanager_np_->getInterfaceHandle(), eledofman, *cell, fieldxx, elementvaluexx, cellvaluexx);
          IO::GMSH::cellWithScalarFieldToStream(cell->Shape(), cellvaluexx, xyze_cell, gmshfilecontentxx);
          }
          {
          LINALG::SerialDenseVector cellvalueyy(DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
          XFEM::computeScalarCellNodeValuesFromElementUnknowns(*actele, &*dofmanager_np_->getInterfaceHandle(), eledofman, *cell, fieldyy, elementvalueyy, cellvalueyy);
          IO::GMSH::cellWithScalarFieldToStream(cell->Shape(), cellvalueyy, xyze_cell, gmshfilecontentyy);
          }
          {
          LINALG::SerialDenseVector cellvaluezz(DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
          XFEM::computeScalarCellNodeValuesFromElementUnknowns(*actele, &*dofmanager_np_->getInterfaceHandle(), eledofman, *cell, fieldzz, elementvaluezz, cellvaluezz);
          IO::GMSH::cellWithScalarFieldToStream(cell->Shape(), cellvaluezz, xyze_cell, gmshfilecontentzz);
          }
          {
          LINALG::SerialDenseVector cellvaluexy(DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
          XFEM::computeScalarCellNodeValuesFromElementUnknowns(*actele, &*dofmanager_np_->getInterfaceHandle(), eledofman, *cell, fieldxy, elementvaluexy, cellvaluexy);
          IO::GMSH::cellWithScalarFieldToStream(cell->Shape(), cellvaluexy, xyze_cell, gmshfilecontentxy);
          }
          {
          LINALG::SerialDenseVector cellvaluexz(DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
          XFEM::computeScalarCellNodeValuesFromElementUnknowns(*actele, &*dofmanager_np_->getInterfaceHandle(), eledofman, *cell, fieldxz, elementvaluexz, cellvaluexz);
          IO::GMSH::cellWithScalarFieldToStream(cell->Shape(), cellvaluexz, xyze_cell, gmshfilecontentxz);
          }
          {
          LINALG::SerialDenseVector cellvalueyz(DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
          XFEM::computeScalarCellNodeValuesFromElementUnknowns(*actele, &*dofmanager_np_->getInterfaceHandle(), eledofman, *cell, fieldyz, elementvalueyz, cellvalueyz);
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
    gmshfilecontent.close();
    gmshfilecontentxx.close();
    gmshfilecontentyy.close();
    gmshfilecontentzz.close();
    gmshfilecontentxy.close();
    gmshfilecontentxz.close();
    gmshfilecontentyz.close();

    if (numplot == 0)
    {
      std::remove(filename.c_str());
      std::remove(filenamexx.c_str());
      std::remove(filenameyy.c_str());
      std::remove(filenamezz.c_str());
      std::remove(filenamexy.c_str());
      std::remove(filenamexz.c_str());
      std::remove(filenameyz.c_str());
    }
    if (screen_out) std::cout << " done" << endl;
  }
#endif

  if (this->physprob_.fieldset_.find(XFEM::PHYSICS::Velx) != this->physprob_.fieldset_.end())
  {
    PlotVectorFieldToGmsh(DRT::UTILS::GetColVersionOfRowVector(xfluiddis_, xfluidstate_.velnp_), "sol_field_xfluid_vel_np","Velocity Solution (Physical) n+1",true, step, time);
    PlotVectorFieldToGmsh(DRT::UTILS::GetColVersionOfRowVector(xfluiddis_, xfluidstate_.veln_),  "sol_field_xfluid_vel_n","Velocity Solution (Physical) n",true, step, time);
    PlotVectorFieldToGmsh(DRT::UTILS::GetColVersionOfRowVector(xfluiddis_, xfluidstate_.velnm_), "sol_field_xfluid_vel_nm","Velocity Solution (Physical) n-1",false, step, time);
    PlotVectorFieldToGmsh(DRT::UTILS::GetColVersionOfRowVector(xfluiddis_, xfluidstate_.accnp_), "sol_field_xfluid_acc_np","Acceleration Solution (Physical) n+1",true, step, time);
    PlotVectorFieldToGmsh(DRT::UTILS::GetColVersionOfRowVector(xfluiddis_, xfluidstate_.accn_),  "sol_field_xfluid_acc_n","Acceleration Solution (Physical) n",true, step, time);
  }
}
/*----------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void FLD::FluidXFluidImplicitTimeInt::MovingFluidOutput()
{
  cout << "FLD::FluidXFluidImplicitTimeInt::MovingFluidOutput()" << endl;
  const Teuchos::ParameterList& xfemparams = DRT::Problem::Instance()->XFEMGeneralParams();
  const bool gmshdebugout = getIntegralValue<int>(xfemparams,"GMSH_DEBUG_OUT")==1;
  const bool screen_out = getIntegralValue<int>(xfemparams,"GMSH_DEBUG_OUT_SCREEN")==1;

  // get column version of fluidstate_.velnp_
  Teuchos::RCP<const Epetra_Vector> col_fluidvelnp = DRT::UTILS::GetColVersionOfRowVector(fluiddis_, fluidstate_.velnp_);

  // get column version of fluidstate_.dispn_
  Teuchos::RCP<const Epetra_Vector> col_fluiddispn;
  if (alefluid_)
  {
    col_fluiddispn = DRT::UTILS::GetColVersionOfRowVector(fluiddis_, fluidstate_.dispn_);
  }

  if (gmshdebugout)
  {
    // fluid velocity (n+1) Gmsh-output
    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("sol_field_fluid_vel_np", step_, 5, screen_out, fluiddis_->Comm().MyPID());
    std::ofstream gmshfilecontent(filename.c_str());
    {
      gmshfilecontent << "View \" " << "Velocity Solution n+1" << "\" {\n";
      for (int i=0; i<fluiddis_->NumMyColElements(); ++i)
      {
        const DRT::Element* actele = fluiddis_->lColElement(i);

        vector<int> lm;
        vector<int> lmowner;
        vector<int> lmstride;
        actele->LocationVector(*fluiddis_, lm, lmowner, lmstride);

        // extract local values from the global vector
        vector<double> myvelnp(lm.size());
        DRT::UTILS::ExtractMyValues(*col_fluidvelnp, myvelnp, lm);

        vector<double> mydispn(lm.size());
        if (alefluid_)
        {
          DRT::UTILS::ExtractMyValues(*col_fluiddispn, mydispn, lm);
        }

        int numnode = actele->NumNode();
        const int nsd = 3;
        LINALG::SerialDenseMatrix elementvalues(nsd, numnode);
        LINALG::SerialDenseMatrix elementpositions(nsd,numnode);

        for (int inode=0; inode<numnode; ++inode)
        {
          elementvalues(0, inode) = myvelnp[0+(inode*4)];
          elementvalues(1, inode) = myvelnp[1+(inode*4)];
          elementvalues(2, inode) = myvelnp[2+(inode*4)];

          elementpositions(0,inode) = actele->Nodes()[inode]->X()[0];
          elementpositions(1,inode) = actele->Nodes()[inode]->X()[1];
          elementpositions(2,inode) = actele->Nodes()[inode]->X()[2];

          if (alefluid_)
          {
            elementpositions(0,inode) += mydispn[0+(inode*4)];
            elementpositions(1,inode) += mydispn[1+(inode*4)];
            elementpositions(2,inode) += mydispn[2+(inode*4)];
          }
        }

        IO::GMSH::cellWithVectorFieldToStream(
            actele->Shape(), elementvalues, elementpositions, gmshfilecontent);
      }
      gmshfilecontent << "};\n";
    }
    gmshfilecontent.close();

    // fluid velocity (n) Gmsh-output
    const std::string filename_veln = IO::GMSH::GetNewFileNameAndDeleteOldFiles("sol_field_fluid_vel_n", step_, 5, screen_out, fluiddis_->Comm().MyPID());
    std::ofstream gmshfilecontent_veln(filename_veln.c_str());

    // get column version of fluidstate_.veln_
    Teuchos::RCP<const Epetra_Vector> col_fluidveln = DRT::UTILS::GetColVersionOfRowVector(fluiddis_, fluidstate_.veln_);

    // get column version of fluidstate_.dispn_
    Teuchos::RCP<const Epetra_Vector> col_fluiddispnm;
    if (alefluid_)
    {
      col_fluiddispnm = DRT::UTILS::GetColVersionOfRowVector(fluiddis_, fluidstate_.dispnm_);
    }

    {
      gmshfilecontent_veln << "View \" " << "Velocity Solution n" << "\" {\n";

      for (int i=0; i<fluiddis_->NumMyColElements(); ++i)
      {
        const DRT::Element* actele = fluiddis_->lColElement(i);

        vector<int> lm;
        vector<int> lmowner;
        vector<int> lmstride;
        actele->LocationVector(*fluiddis_, lm, lmowner, lmstride);

        // extract local values from the global vector
        vector<double> myveln(lm.size());
        DRT::UTILS::ExtractMyValues(*col_fluidveln, myveln, lm);

        vector<double> mydispnm(lm.size());
        if (alefluid_)
        {
          DRT::UTILS::ExtractMyValues(*col_fluiddispnm, mydispnm, lm);
        }

        int numparam = actele->NumNode();
        LINALG::SerialDenseMatrix xyze(3,numparam);
        LINALG::SerialDenseMatrix elementvalues(3, numparam);
        for (int inode=0; inode<numparam; ++inode)
        {
          elementvalues(0, inode) = myveln[0+(inode*4)];
          elementvalues(1, inode) = myveln[1+(inode*4)];
          elementvalues(2, inode) = myveln[2+(inode*4)];

          xyze(0,inode) = actele->Nodes()[inode]->X()[0];
          xyze(1,inode) = actele->Nodes()[inode]->X()[1];
          xyze(2,inode) = actele->Nodes()[inode]->X()[2];

          if (alefluid_ and step_>1)
          {
            xyze(0,inode) += mydispnm[0+(inode*4)];
            xyze(1,inode) += mydispnm[1+(inode*4)];
            xyze(2,inode) += mydispnm[2+(inode*4)];
          }
        }

        IO::GMSH::cellWithVectorFieldToStream(
            actele->Shape(), elementvalues, xyze, gmshfilecontent_veln);
      }
      gmshfilecontent_veln << "};\n";
    }
    gmshfilecontent_veln.close();

    // fluid velocity (n-1) Gmsh-output
    const std::string filename_velnm = IO::GMSH::GetNewFileNameAndDeleteOldFiles("sol_field_fluid_vel_nm", step_, 5, screen_out, fluiddis_->Comm().MyPID());
    std::ofstream gmshfilecontent_velnm(filename_velnm.c_str());

    // get column version of fluidstate_.velnm_
    Teuchos::RCP<const Epetra_Vector> col_fluidvelnm = DRT::UTILS::GetColVersionOfRowVector(fluiddis_, fluidstate_.velnm_);

    // get column version of fluidstate_.dispnmm_
    Teuchos::RCP<const Epetra_Vector> col_fluiddispnmm;
    if (alefluid_)
    {
      col_fluiddispnmm = DRT::UTILS::GetColVersionOfRowVector(fluiddis_, fluidstate_.dispnmm_);
    }

    {
      gmshfilecontent_velnm << "View \" " << "Velocity Solution n" << "\" {\n";

      for (int i=0; i<fluiddis_->NumMyColElements(); ++i)
      {
        const DRT::Element* actele = fluiddis_->lColElement(i);

        vector<int> lm;
        vector<int> lmowner;
        vector<int> lmstride;
        actele->LocationVector(*fluiddis_, lm, lmowner, lmstride);

        // extract local values from the global vector
        vector<double> myvelnm(lm.size());
        DRT::UTILS::ExtractMyValues(*col_fluidvelnm, myvelnm, lm);

        vector<double> mydispnmm(lm.size());
        if (alefluid_)
        {
          DRT::UTILS::ExtractMyValues(*col_fluiddispnmm, mydispnmm, lm);
        }

        int numparam = actele->NumNode();
        LINALG::SerialDenseMatrix xyze(3,numparam);
        LINALG::SerialDenseMatrix elementvalues(3, numparam);
        for (int inode=0; inode<numparam; ++inode)
        {
          elementvalues(0, inode) = myvelnm[0+(inode*4)];
          elementvalues(1, inode) = myvelnm[1+(inode*4)];
          elementvalues(2, inode) = myvelnm[2+(inode*4)];

          xyze(0,inode) = actele->Nodes()[inode]->X()[0];
          xyze(1,inode) = actele->Nodes()[inode]->X()[1];
          xyze(2,inode) = actele->Nodes()[inode]->X()[2];

          if (alefluid_ and step_>1)
          {
            xyze(0,inode) += mydispnmm[0+(inode*4)];
            xyze(1,inode) += mydispnmm[1+(inode*4)];
            xyze(2,inode) += mydispnmm[2+(inode*4)];
          }
        }

        IO::GMSH::cellWithVectorFieldToStream(
            actele->Shape(), elementvalues, xyze, gmshfilecontent_velnm);
      }
      gmshfilecontent_velnm << "};\n";
    }
    gmshfilecontent_velnm.close();

  // fluid pressure Gmsh-output
    const std::string filename_pres = IO::GMSH::GetNewFileNameAndDeleteOldFiles("sol_field_fluid_pres_np", step_, 5, screen_out, fluiddis_->Comm().MyPID());
    std::ofstream gmshfilecontent_pres(filename_pres.c_str());

    {
      gmshfilecontent_pres << "View \" " << "Pressure Solution n+1" << "\" {\n";
      for (int i=0; i<fluiddis_->NumMyColElements(); ++i)
      {
        const DRT::Element* actele = fluiddis_->lColElement(i);

        vector<int> lm;
        vector<int> lmowner;
        vector<int> lmstride;
        actele->LocationVector(*fluiddis_, lm, lmowner, lmstride);

        //extract local values from the global vector
        vector<double> myvelnp(lm.size());
        DRT::UTILS::ExtractMyValues(*col_fluidvelnp, myvelnp, lm);

        vector<double> mydispn(lm.size());
        if (alefluid_)
        {
          DRT::UTILS::ExtractMyValues(*col_fluiddispn, mydispn, lm);
        }

        int numparam = actele->NumNode();
        LINALG::SerialDenseMatrix xyze(3,numparam);
        LINALG::SerialDenseVector elementvalues(numparam);

        for (int inode=0; inode<numparam; ++inode)
        {
          elementvalues(inode) = myvelnp[3+(inode*4)];

          xyze(0,inode) = actele->Nodes()[inode]->X()[0];
          xyze(1,inode) = actele->Nodes()[inode]->X()[1];
          xyze(2,inode) = actele->Nodes()[inode]->X()[2];

          if (alefluid_ and step_>1)
          {
            xyze(0,inode) += mydispn[0+(inode*4)];
            xyze(1,inode) += mydispn[1+(inode*4)];
            xyze(2,inode) += mydispn[2+(inode*4)];
          }
        }

        IO::GMSH::cellWithScalarFieldToStream(actele->Shape(), elementvalues, xyze, gmshfilecontent_pres);
      }
     gmshfilecontent_pres << "};\n";
    }

    gmshfilecontent_pres.close();
    if (screen_out) std::cout << " done" << endl;
  }
}
/*-----------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------*/
void FLD::FluidXFluidImplicitTimeInt::FluidXFluidBoundaryOutput(const Teuchos::RCP<DRT::Discretization> fluidxfluidboundarydis)
{
  // the fluid-fluid-boundary output

  // create interface DOF vectors using the fluid parallel distribution
  Teuchos::RCP<Epetra_Vector> fivelnpcol   = LINALG::CreateVector(*fluidxfluidboundarydis->DofColMap(),true);
  Teuchos::RCP<Epetra_Vector> fivelncol    = LINALG::CreateVector(*fluidxfluidboundarydis->DofColMap(),true);
  Teuchos::RCP<Epetra_Vector> fiaccnpcol   = LINALG::CreateVector(*fluidxfluidboundarydis->DofColMap(),true);
  Teuchos::RCP<Epetra_Vector> fiaccncol    = LINALG::CreateVector(*fluidxfluidboundarydis->DofColMap(),true);
  Teuchos::RCP<Epetra_Vector> fidispnpcol  = LINALG::CreateVector(*fluidxfluidboundarydis->DofColMap(),true);
  Teuchos::RCP<Epetra_Vector> fitruerescol = LINALG::CreateVector(*fluidxfluidboundarydis->DofColMap(),true);

  // map to fluid parallel distribution
  LINALG::Export(*fxfivelnp_,*fivelnpcol);
  LINALG::Export(*fxfiveln_ ,*fivelncol);
  LINALG::Export(*fxfiaccnp_,*fiaccnpcol);
  LINALG::Export(*fxfiaccn_,*fiaccncol);
  LINALG::Export(*fxfidispnp_,*fidispnpcol);

  // print redundant arrays on proc 0 (Gmsh)
  if (fluidxfluidboundarydis->Comm().MyPID() == 0)
  {
    PrintFluidFluidBoundaryVectorField(fluidxfluidboundarydis,fidispnpcol, fivelnpcol  , "sol_fluidfluidiface_velnp", "fluidfluidboundary velocity n+1");
    PrintFluidFluidBoundaryVectorField(fluidxfluidboundarydis,fidispnpcol, fivelncol   , "sol_fluidfluidiface_veln" , "fluidfluidboundary velocity n");
    PrintFluidFluidBoundaryVectorField(fluidxfluidboundarydis,fidispnpcol, fiaccnpcol  , "sol_fluidfluidiface_accnp", "fluidfluidboundary acceleration n+1");
    PrintFluidFluidBoundaryVectorField(fluidxfluidboundarydis,fidispnpcol, fiaccncol   , "sol_fluidfluidiface_accn" , "fluidfluidboundary acceleration n");
    PrintFluidFluidBoundaryVectorField(fluidxfluidboundarydis,fidispnpcol, fidispnpcol , "sol_fluidfluidiface_disp" , "fluidfluidboundary displacement n+1");
    PrintFluidFluidBoundaryVectorField(fluidxfluidboundarydis,fidispnpcol, fitruerescol, "sol_fluidfluidiface_force", "fluidfluidboundary force");
   }

  if (fluidxfluidboundarydis->Comm().MyPID() == 0 && fidispnpcol->MyLength() >= 3)
  {
    std::ofstream f;
    const std::string fname = DRT::Problem::Instance()->OutputControlFile()->FileName()
                                     + ".outfluidfluiddispnp.txt";
    if (Step() <= 1)
      f.open(fname.c_str(),std::fstream::trunc);
    else
      f.open(fname.c_str(),std::fstream::ate | std::fstream::app);

    f << Time() << " " << (*fidispnpcol)[0] << "  " << "\n";

    f.close();
  }

  if (fluidxfluidboundarydis->Comm().MyPID() == 0 && fivelnpcol->MyLength() >= 3)
  {
    std::ofstream f;
    const std::string fname = DRT::Problem::Instance()->OutputControlFile()->FileName()
                                     + ".outfluidfluidvelnp.txt";
    if (Step() <= 1)
      f.open(fname.c_str(),std::fstream::trunc);
    else
      f.open(fname.c_str(),std::fstream::ate | std::fstream::app);

    f << Time() << " " << (*fivelnpcol)[0] << "  " << "\n";

    f.close();
  }

  if (fluidxfluidboundarydis->Comm().MyPID() == 0 && fivelncol->MyLength() >= 3)
  {
    std::ofstream f;
    const std::string fname = DRT::Problem::Instance()->OutputControlFile()->FileName()
                                     + ".outfluidfluidveln.txt";
    if (Step() <= 1)
      f.open(fname.c_str(),std::fstream::trunc);
    else
      f.open(fname.c_str(),std::fstream::ate | std::fstream::app);

    f << Time() << " " << (*fivelncol)[0] << "  " << "\n";

    f.close();
  }

  if (fluidxfluidboundarydis->Comm().MyPID() == 0 && fiaccnpcol->MyLength() >= 3)
  {
    std::ofstream f;
    const std::string fname = DRT::Problem::Instance()->OutputControlFile()->FileName()
                                     + ".outfluidfluidaccnp.txt";
    if (Step() <= 1)
      f.open(fname.c_str(),std::fstream::trunc);
    else
      f.open(fname.c_str(),std::fstream::ate | std::fstream::app);

    f << Time() << " " << (*fiaccnpcol)[0] << "  " << "\n";

    f.close();
  }

  if (fluidxfluidboundarydis->Comm().MyPID() == 0 && fiaccncol->MyLength() >= 3)
  {
    std::ofstream f;
    const std::string fname = DRT::Problem::Instance()->OutputControlFile()->FileName()
                                 + ".outfluidfluidaccn.txt";
    if (Step() <= 1)
      f.open(fname.c_str(),std::fstream::trunc);
    else
      f.open(fname.c_str(),std::fstream::ate | std::fstream::app);

    f << Time() << " " << (*fiaccncol)[0] << "  " << "\n";

    f.close();
  }
}
/*----------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------*/
void FLD::FluidXFluidImplicitTimeInt::PrintFluidFluidBoundaryVectorField(
    const Teuchos::RCP<DRT::Discretization> fluidxfluidboundarydis,
    const Teuchos::RCP<Epetra_Vector>   displacementfield,
    const Teuchos::RCP<Epetra_Vector>   vectorfield,
    const std::string filestr,
    const std::string name_in_gmsh
    ) const
{
  const Teuchos::ParameterList& xfemparams = DRT::Problem::Instance()->XFEMGeneralParams();
  const bool gmshdebugout = getIntegralValue<int>(xfemparams,"GMSH_DEBUG_OUT")==1;
  const bool screen_out = getIntegralValue<int>(xfemparams,"GMSH_DEBUG_OUT_SCREEN")==1;

  if (gmshdebugout)
  {
    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles(filestr, Step(), 5, screen_out, fluidxfluidboundarydis->Comm().MyPID());
    std::ofstream gmshfilecontent(filename.c_str());

    {
      gmshfilecontent << "View \" " << name_in_gmsh << " \" {\n";
      for (int i=0; i<fluidxfluidboundarydis->NumMyColElements(); ++i)
      {
        const DRT::Element* bele = fluidxfluidboundarydis->lColElement(i);
        vector<int> lm;
        vector<int> lmowner;
        vector<int> lmstride;
        bele->LocationVector(*fluidxfluidboundarydis, lm, lmowner, lmstride);

        // extract local values from the global vector
        vector<double> myvelnp(lm.size());
        DRT::UTILS::ExtractMyValues(*vectorfield, myvelnp, lm);

        vector<double> mydisp(lm.size());
        DRT::UTILS::ExtractMyValues(*displacementfield, mydisp, lm);

        const int nsd = 3;
        const int numnode = bele->NumNode();
        LINALG::SerialDenseMatrix elementvalues(nsd,numnode);
        LINALG::SerialDenseMatrix elementpositions(nsd,numnode);
        int counter = 0;
        for (int iparam=0; iparam<numnode; ++iparam)
        {
          const DRT::Node* node = bele->Nodes()[iparam];
          const double* pos = node->X();
          for (int isd = 0; isd < nsd; ++isd)
          {
            elementvalues(isd,iparam) = myvelnp[counter];
            elementpositions(isd,iparam) = pos[isd] + mydisp[counter];
            counter++;
          }
        }

        IO::GMSH::cellWithVectorFieldToStream(
            bele->Shape(), elementvalues, elementpositions, gmshfilecontent);
      }
      gmshfilecontent << "};\n";
    }
    gmshfilecontent.close();
    if (screen_out) std::cout << " done" << endl;
  }
}

// -------------------------------------------------------------------
// set general fluid parameter (AE 01/2011)
// -------------------------------------------------------------------

void FLD::FluidXFluidImplicitTimeInt::SetElementGeneralFluidParameter()
{
  ParameterList eleparams;

  eleparams.set("action","set_general_fluid_parameter");

  // set general element parameters
  eleparams.set("form of convective term",convform_);
  eleparams.set("fs subgrid viscosity",fssgv_);
  eleparams.set("Linearisation",newton_);
  eleparams.set("Physical Type", physicaltype_);

  // parameter for stabilization
  eleparams.sublist("STABILIZATION") = params_.sublist("STABILIZATION");

  // parameter for turbulent flow
  eleparams.sublist("TURBULENCE MODEL") = params_.sublist("TURBULENCE MODEL");

  //set time integration scheme
  eleparams.set("TimeIntegrationScheme", timealgo_);

  // call standard loop over elements
  fluiddis_->Evaluate(eleparams,null,null,null,null,null);
  return;
}


// -------------------------------------------------------------------
// set general time parameter (AE 01/2011)
// -------------------------------------------------------------------

void FLD::FluidXFluidImplicitTimeInt::SetElementTimeParameter()
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
  // is not available in FluidXFluid coupling
  /*
  else if (timealgo_==INPAR::FLUID::timeint_afgenalpha)
  {
    eleparams.set("total time",time_-(1-alphaF_)*dta_);
    eleparams.set("alphaF",alphaF_);
    eleparams.set("alphaM",alphaM_);
    eleparams.set("gamma",gamma_);
  }*/
  else
  {
    eleparams.set("total time",time_);
  }

  // call standard loop over elements
  fluiddis_->Evaluate(eleparams,null,null,null,null,null);
  return;
}
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
void FLD::FluidXFluidImplicitTimeInt::CreatePatchBoxes(
  std::map<int,GEO::CUT::BoundingBox>& patchboxes
  ) const
{
  cout << "CreatePatchBoxes" << endl;
  // get column version of the displacement vector
  Teuchos::RCP<const Epetra_Vector> col_fluiddispnm = DRT::UTILS::GetColVersionOfRowVector(fluiddis_, fluidstate_.dispnm_);

  // Map of all boxes of patch fluid discretization
  for (int pele=0; pele<fluiddis_->NumMyColElements(); ++pele)
  {
    const DRT::Element* actpele = fluiddis_->lColElement(pele);
    const DRT::Node*const* pelenodes = actpele->Nodes();
    //cout << " pele id " << actpele->Id() << endl;

    vector<int> lm;
    vector<int> lmowner;
    vector<int> lmstride;
    actpele->LocationVector(*fluiddis_, lm, lmowner, lmstride);

    vector<double> mydispnm(lm.size());
    DRT::UTILS::ExtractMyValues(*col_fluiddispnm, mydispnm, lm);

    GEO::CUT::BoundingBox patchbox;
    // patchboxes
    for (int pnode = 0; pnode < actpele->NumNode(); ++pnode)
    {
      // the coordinates of the actuall node
      LINALG::Matrix<3,1> pnodepos(true);
      pnodepos(0,0) = pelenodes[pnode]->X()[0] + mydispnm[0+(pnode*4)];
      pnodepos(1,0) = pelenodes[pnode]->X()[1] + mydispnm[1+(pnode*4)];
      pnodepos(2,0) = pelenodes[pnode]->X()[2] + mydispnm[2+(pnode*4)];;

      //cout <<  pnodepos(0,0) << " " << pnodepos(1,0) << " " << pnodepos(2,0)<< endl;

      // fill the patchbox
      patchbox.AddPoint(pnodepos);
    }
    //patchbox.Print();
    patchboxes[actpele->Id()] = patchbox;
  }
}

#endif /* CCADISCRET       */
