/*!----------------------------------------------------------------------
\file combust_fluidimplicitintegration.cpp
\brief class holding implicit time integration schemes for combustion problems

This class is a merger of the standard fluid time integration and the XFSI time integration classes.
Thus, a mayor part of the code is a duplicate, but the class also contains some new and modified
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
#ifdef CCADISCRET

#include <stdio.h>

#include "combust_defines.H"
#include "combust_fluidimplicitintegration.H"
//#include "combust_defines.H"
#include "combust3_interpolation.H"
#include "../drt_fluid/time_integration_scheme.H"
#include "../drt_fluid/fluid_utils.H"
#include "../drt_geometry/position_array.H"
#include "../drt_xfem/dof_management.H"
#include "../drt_xfem/dof_distribution_switcher.H"
#include "../drt_xfem/startvalues.H"
#include "../drt_xfem/enrichment_utils.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/matlist.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_function.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../linalg/linalg_ana.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>


/*------------------------------------------------------------------------------------------------*
 | constructor                                                                        henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
FLD::CombustFluidImplicitTimeInt::CombustFluidImplicitTimeInt(
    Teuchos::RCP<DRT::Discretization> actdis,
    LINALG::Solver&                   solver,
    ParameterList&                    params,
    IO::DiscretizationWriter&         output) :
  // call constructor for "nontrivial" objects
  discret_(actdis),
  solver_ (solver),
  params_ (params),
  xparams_(params.sublist("XFEM")),
//  output_ (output),
  output_ (rcp(new IO::DiscretizationWriter(actdis))), // so ist es bei Axel
  myrank_(discret_->Comm().MyPID()),
  cout0_(discret_->Comm(), std::cout),
  combusttype_(Teuchos::getIntegralValue<INPAR::COMBUST::CombustionType>(params_.sublist("COMBUSTION FLUID"),"COMBUSTTYPE")),
  veljumptype_(Teuchos::getIntegralValue<INPAR::COMBUST::VelocityJumpType>(params_.sublist("COMBUSTION FLUID"),"VELOCITY_JUMP_TYPE")),
  normaltensionjumptype_(Teuchos::getIntegralValue<INPAR::COMBUST::NormalTensionJumpType>(params_.sublist("COMBUSTION FLUID"),"NORMAL_TENSION_JUMP_TYPE")),
  flamespeed_(params_.sublist("COMBUSTION FLUID").get<double>("LAMINAR_FLAMESPEED")),
  nitschevel_(params_.sublist("COMBUSTION FLUID").get<double>("NITSCHE_VELOCITY")),
  nitschepres_(params_.sublist("COMBUSTION FLUID").get<double>("NITSCHE_PRESSURE")),
  condensation_(xparams_.get<bool>("DLM_condensation")),
  surftensapprox_(Teuchos::getIntegralValue<INPAR::COMBUST::SurfaceTensionApprox>(params_.sublist("COMBUSTION FLUID"),"SURFTENSAPPROX")),
  surftenscoeff_(params_.sublist("COMBUSTION FLUID").get<double>("SURFTENSCOEFF")),
  connected_interface_(Teuchos::getIntegralValue<int>(params_.sublist("COMBUSTION FLUID"),"CONNECTED_INTERFACE")),
  step_(0),
  time_(0.0),
  stepmax_ (params_.get<int>   ("max number timesteps")),
  maxtime_ (params_.get<double>("total time")),
  dta_     (params_.get<double> ("time step size")),
  dtp_     (params_.get<double> ("time step size")),
  timealgo_(params_.get<INPAR::FLUID::TimeIntegrationScheme>("time int algo")),
  //startalgo_(params_.get<INPAR::FLUID::TimeIntegrationScheme>("start time int algo")),
  theta_   (params_.get<double>("theta")),
  initstatsol_(Teuchos::getIntegralValue<int>(params_.sublist("COMBUSTION FLUID"),"INITSTATSOL")),
  itemax_(params_.get<int>("max nonlin iter steps")),
  extrapolationpredictor_(params.get("do explicit predictor",false)),
  uprestart_(params.get("write restart every", -1)),
  upres_(params.get("write solution every", -1)),
  writestresses_(params.get<int>("write stresses", 0)),
  flamefront_(Teuchos::null)
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
  //------------------------------------------------------------------------------------------------
  // future development: connect degrees of freedom for periodic boundary conditions     henke 01/09
  //------------------------------------------------------------------------------------------------

  //------------------------------------------------------------------------------------------------
  // prepare XFEM (initial degree of freedom management)
  //------------------------------------------------------------------------------------------------
  physprob_.xfemfieldset_.clear();
  // declare physical fields to be enriched (discontinuous XFEM fields)
  physprob_.xfemfieldset_.insert(XFEM::PHYSICS::Velx);
  physprob_.xfemfieldset_.insert(XFEM::PHYSICS::Vely);
  physprob_.xfemfieldset_.insert(XFEM::PHYSICS::Velz);
  physprob_.xfemfieldset_.insert(XFEM::PHYSICS::Pres);
  // define approach for extra stress field (stress-based Lagrange Multiplier approach)
  physprob_.elementAnsatz_ = rcp(new COMBUST::TauPressureAnsatz());

  // create dummy instance of interfacehandle holding no flamefront and hence no integration cells
  Teuchos::RCP<COMBUST::InterfaceHandleCombust> ihdummy = rcp(new COMBUST::InterfaceHandleCombust(discret_,Teuchos::null,Teuchos::null));
  // create dummy instance of dof manager assigning standard enrichments to all nodes
  const Teuchos::RCP<XFEM::DofManager> dofmanagerdummy = rcp(new XFEM::DofManager(ihdummy,physprob_.xfemfieldset_,*physprob_.elementAnsatz_,xparams_));

  // pass dof information to elements (no enrichments yet, standard FEM!)
  TransferDofInformationToElements(ihdummy, dofmanagerdummy);

  // save dofmanager to be able to plot Gmsh stuff in Output()
  dofmanagerForOutput_ = dofmanagerdummy;

//  {
//    ParameterList eleparams;
//    eleparams.set("action","set_output_mode");
//    eleparams.set("output_mode",true);
//    discret_->Evaluate(eleparams);
//  }

  // ensure that degrees of freedom in the discretization have been set
  discret_->FillComplete();

  output_->WriteMesh(0,0.0);

  // store a dofset with the complete fluid unknowns
  standarddofset_ = Teuchos::rcp(new DRT::DofSet());
  standarddofset_->Reset();
  standarddofset_->AssignDegreesOfFreedom(*discret_,0,0);
  // split based on complete fluid field
  FLD::UTILS::SetupFluidSplit(*discret_,*standarddofset_,3,velpressplitterForOutput_);

//  {
//    ParameterList eleparams;
//    eleparams.set("action","set_output_mode");
//    eleparams.set("output_mode",false);
//    discret_->Evaluate(eleparams);
//  }

  //------------------------------------------------------------------------------------------------
  // get dof layout from the discretization to construct vectors and matrices
  //------------------------------------------------------------------------------------------------

  // parallel dof distribution contained in dofrowmap: local (LID) <-> global (GID) dof numbering
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // get layout of velocity and pressure dofs in a vector
  const int numdim = params_.get<int>("number of velocity degrees of freedom");
  FLD::UTILS::SetupFluidSplit(*discret_,numdim,velpressplitter_);

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
  state_.accnp_ = LINALG::CreateVector(*dofrowmap,true);
  state_.accn_  = LINALG::CreateVector(*dofrowmap,true);

  state_.nodalDofDistributionMap_.clear();
  state_.elementalDofDistributionMap_.clear();

  dofmanagerForOutput_->fillDofRowDistributionMaps(
      state_.nodalDofDistributionMap_,
      state_.elementalDofDistributionMap_);

  // history vector: scheint in state-Konstrukt nicht nötig zu sein!
//  hist_ = LINALG::CreateVector(*dofrowmap,true);

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
  // get density from elements
  //------------------------------------------------------------------------------------------------
//  {
//    ParameterList eleparams;
//    eleparams.set("action","get_density");
//    std::cout << "Warning: two-phase flows have different densities, evaluate(get_density) returns 1.0" << std::endl;
//    discret_->Evaluate(eleparams,null,null,null,null,null);
//    density_ = eleparams.get<double>("density");
//    if (density_ <= 0.0) dserror("received negative or zero density value from elements");
//  }
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
 | prepare a fluid time step                                                          henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::PrepareTimeStep()
{
  // update interface handle
  //ih_n_ = ih_np_;

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
  dta_      = params_.get<double> ("time step size");
  dtp_      = params_.get<double> ("time step size");
  theta_    = params_.get<double>("theta");
  timealgo_ = params_.get<INPAR::FLUID::TimeIntegrationScheme>("time int algo");
  itemax_   = params_.get<int>("max nonlin iter steps");

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
      itemax_ = params_.get<int>("max nonlin iter steps init stat sol");
    }
    // compute first (instationary) time step differently
    // remark: usually backward Euler (theta=1.0) to damp inital pertubations
    else if (step_==1)
    {
      // get starting 'theta' for first time step
      theta_ = params_.get<double>("start theta");
      cout0_ << "/!\\ warning: first timestep computed with theta =  " << theta_ << endl;
    }
    // regular time step
    else if (step_ > 1)
    {
      theta_ = params_.get<double>("theta");
    }
    else
      dserror("number of time step is wrong");
    break;
  }
  case INPAR::FLUID::timeint_bdf2:
  {
    // do a backward Euler step for the first time step
    if (step_==1)
    {
      timealgo_ = INPAR::FLUID::timeint_one_step_theta;
      theta_ = params_.get<double>("start theta");
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
  if (extrapolationpredictor_)
  {
    if (step_>1)
    {
      double timealgo_constant=theta_;

      TIMEINT_THETA_BDF2::ExplicitPredictor(
        "default",
        state_.veln_, 
        state_.velnm_, 
        state_.accn_,
        velpressplitter_,
        timealgo_, 
        timealgo_constant,
        dta_, 
        dtp_,
        state_.velnp_,
        discret_->Comm());
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
    discret_->SetState("velnp",state_.velnp_);
    // predicted dirichlet values
    // velnp then also holds prescribed new dirichlet values
    discret_->EvaluateDirichletXFEM(eleparams,state_.velnp_,null,null,null,dbcmaps_);
    discret_->ClearState();

    // evaluate Neumann conditions
    neumann_loads_->PutScalar(0.0);
    discret_->EvaluateNeumann(eleparams,*neumann_loads_);
    discret_->ClearState();
  }

}

/*------------------------------------------------------------------------------------------------*
 | hand over information about (XFEM) degrees of freedom to elements                     ag 04/09 |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::TransferDofInformationToElements(
    const Teuchos::RCP<COMBUST::InterfaceHandleCombust> interfacehandle,
    const Teuchos::RCP<XFEM::DofManager> dofmanager
    )
{
  ParameterList eleparams;
  eleparams.set("action","store_xfem_info");
  eleparams.set("dofmanager",dofmanager);
  eleparams.set("DLM_condensation",xparams_.get<bool>("DLM_condensation"));
  eleparams.set("boundaryRatioLimit",xparams_.get<double>("boundaryRatioLimit"));
  eleparams.set("interfacehandle",interfacehandle);
  discret_->Evaluate(eleparams);
}

/*------------------------------------------------------------------------------------------------*
 | import geometrical information about the interface (integration cells) from the combustion     |
 | algorithm and incorporate it into the fluid field                                  henke 03/09 |
 |
 | remark: Within this routine, no parallel re-distribution is allowed to take place. Before and  |
 | after this function, it's ok to do so.                                                    axel |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::IncorporateInterface(
       const Teuchos::RCP<COMBUST::InterfaceHandleCombust>& interfacehandleNP,
       const Teuchos::RCP<COMBUST::InterfaceHandleCombust>& interfacehandleN)
{
  // information about interface is imported via ADAPTER::FluidCombust::ImportInterface()
  if (false) // flamefront_ == Teuchos::null)
  {
    dserror("combustion time integration scheme cannot see flame front");
  }

  //cout << *(state_.veln_)<< endl;
  /* momentan gebe ich den ElementAnsatz (Lagrange Multiplier Zeug) nicht an den DofManager weiter,
   * vielleicht brauche ich es aber. Hängt das hier eigentlich nicht direkt von InputParametern ab? */
  //const COMBUST::CombustElementAnsatz elementAnsatz;

  // build instance of DofManager with information about the interface from the interfacehandle
  // remark: DofManager is rebuilt in every inter-field iteration step, because number and position
  // of enriched degrees of freedom change constantly.
  const Teuchos::RCP<XFEM::DofManager> dofmanager = rcp(new XFEM::DofManager(interfacehandleNP,physprob_.xfemfieldset_,*physprob_.elementAnsatz_,xparams_));

  //if (true) // INPAR::XFEM::timeintegration_semilagrangian
  //{
  // temporarely save old dofmanager
  const RCP<XFEM::DofManager> olddofmanager = dofmanagerForOutput_;
  //}

  // save dofmanager to be able to plot Gmsh stuff in Output()
  dofmanagerForOutput_ = dofmanager;

  // tell elements about the dofs and the integration
  TransferDofInformationToElements(interfacehandleNP, dofmanager);

  // print global and element dofmanager to Gmsh
  std::cout<< "Dofmanager and InterfaceHandle to Gmsh" << std::endl;
  dofmanager->toGmsh(step_);
  interfacehandleNP->toGmsh(step_);

  // get old dofmaps, compute a new one and get the new one, too
  const Epetra_Map olddofrowmap = *discret_->DofRowMap();
  //if (true) // INPAR::XFEM::timeintegration_semilagrangian
  //{
  const Epetra_Map olddofcolmap = *discret_->DofColMap();
  map<XFEM::DofKey<XFEM::onNode>,XFEM::DofGID> oldNodalDofColDistrib;
  olddofmanager->fillNodalDofColDistributionMap(oldNodalDofColDistrib);
  //}
  // assign degrees of freedom
  // remark: - assign degrees of freedom (first slot)
  //         - build geometry for (Neumann) boundary conditions (third slot);
  //           without Neumann boundary conditions Fillcomplete(true,false,false) will also work
  discret_->FillComplete(true,false,true);
  const Epetra_Map& newdofrowmap = *discret_->DofRowMap();

  // remark: 'true' is needed to prevent iterative solver from crashing
  discret_->ComputeNullSpaceIfNecessary(solver_.Params(),true);

  { // area for dofswitcher and startvalues
    //const std::map<XFEM::DofKey<XFEM::onNode>, XFEM::DofGID> oldNodalDofDistributionMap(state_.nodalDofDistributionMap_);
    std::map<XFEM::DofKey<XFEM::onNode>, XFEM::DofGID> oldNodalDofDistributionMap(state_.nodalDofDistributionMap_);
    const std::map<XFEM::DofKey<XFEM::onElem>, XFEM::DofGID> oldElementalDofDistributionMap(state_.elementalDofDistributionMap_);
    dofmanager->fillDofRowDistributionMaps(
        state_.nodalDofDistributionMap_,
        state_.elementalDofDistributionMap_);

    //    cout << "elemental dof set with size " <<
    //        state_.elementalDofDistributionMap_.size() <<
    //        state_.nodalDofDistributionMap_.size() << endl;
    //    for (std::map<XFEM::DofKey<XFEM::onElem>, XFEM::DofGID>::const_iterator i=oldElementalDofDistributionMap.begin();
    //        i!=oldElementalDofDistributionMap.end();
    //        i++)
    //      cout << i->first << ", " << i->second << endl;
    //    cout << "elemental dofset done " << endl;
    //    discret_->Comm().Barrier();

    // create switcher
    const XFEM::DofDistributionSwitcher dofswitch(
        interfacehandleNP, dofmanager,
        olddofrowmap, newdofrowmap,
        oldNodalDofDistributionMap, state_.nodalDofDistributionMap_,
        oldElementalDofDistributionMap, state_.elementalDofDistributionMap_
    );

    // ----------------------------------------------------------
    // extract enrichment dofkeys and values before they are lost
    // ----------------------------------------------------------
    vector<RCP<Epetra_Vector> > oldColStateVectors;

    RCP<Epetra_Vector> accnp = rcp(new Epetra_Vector(olddofcolmap,true));
    LINALG::Export(*state_.accnp_,*accnp);
    oldColStateVectors.push_back(accnp);

    RCP<Epetra_Vector> accn = rcp(new Epetra_Vector(olddofcolmap,true));
    LINALG::Export(*state_.accn_,*accn);
    oldColStateVectors.push_back(accn);

    RCP<Epetra_Vector> velnp = rcp(new Epetra_Vector(olddofcolmap,true));
    LINALG::Export(*state_.velnp_,*velnp);
    oldColStateVectors.push_back(velnp);

    RCP<Epetra_Vector> veln = rcp(new Epetra_Vector(olddofcolmap,true));
    LINALG::Export(*state_.veln_,*veln);
    oldColStateVectors.push_back(veln);

    //    discret_->Comm().Barrier();
    //    cout << "acceleration is " << *state_.accn_;
    //    cout << "velocity is " << *state_.veln_;
    //
    //    map<int,set<XFEM::DofKey<XFEM::onNode> > > oldEnrNodalDofSet;
    //    map<XFEM::DofKey<XFEM::onNode>,XFEM::DofGID> oldEnrNodalDofDistrib;
    //    vector<map<int,double> > oldEnrValues(stateVector.size());
    //
    //    {
    ////      cout << "here1" << endl;
    ////      for (map<XFEM::DofKey<XFEM::onNode>,XFEM::DofGID>::const_iterator i=oldNodalDofColDistrib.begin();
    ////          i != oldNodalDofColDistrib.end();i++)
    ////        cout << i->first << ", " << i->second << endl;
    //      dofswitch.extractEnrMapCombust(
    //          stateVector,
    //          olddofcolmap,
    //          oldNodalDofColDistrib,
    //          oldEnrNodalDofSet,
    //          oldEnrNodalDofDistrib,
    //          oldEnrValues
    //      );
    //    }

    //    cout << "here2" << endl;
    //    cout << "enrsizes are: ";
    //    for (size_t i=0;i<oldEnrValues.size();i++)
    //      cout << oldEnrValues[i].size() << ", ";
    //    for (map<XFEM::DofKey<XFEM::onNode>,XFEM::DofGID>::const_iterator i=oldEnrNodalDofDistrib.begin();
    //        i != oldEnrNodalDofDistrib.end();
    //        i++)
    //      cout << i->first << ", " << i->second << endl;
    //    for (map<int,double>::const_iterator i=oldEnrValues[3].begin();
    //        i!=oldEnrValues[3].end();
    //        i++)
    //      cout << i->first << ", " << i->second << endl;

    // --------------------------------------------
    // switch state vectors to new dof distribution
    // --------------------------------------------
    //cout << "old solution was " << *state_.velnp_;
    cout0_ << " Initialize system vectors..." << endl;
    // accelerations at time n and n-1

    dofswitch.mapVectorToNewDofDistributionCombust(state_.accnp_);
    dofswitch.mapVectorToNewDofDistributionCombust(state_.accn_);

    // velocities and pressures at time n+1, n and n-1
    dofswitch.mapVectorToNewDofDistributionCombust(state_.velnp_); // use old velocity as start value
    dofswitch.mapVectorToNewDofDistributionCombust(state_.veln_);
    dofswitch.mapVectorToNewDofDistributionCombust(state_.velnm_);

if (false) // INPAR::XFEM::timeintegration_semilagrangian
{
    // here the interface is known for the first time and so
    // all initial values can be set including enrichment values
    if (step_ == 1)
      SetEnrichmentField(dofmanager,newdofrowmap);

    if (step_ > 1)
    {
      const size_t oldsize = oldColStateVectors.size();
      vector<RCP<Epetra_Vector> > newRowStateVectors;
      newRowStateVectors.push_back(state_.accnp_);
      newRowStateVectors.push_back(state_.accn_);
      newRowStateVectors.push_back(state_.velnp_);
      newRowStateVectors.push_back(state_.veln_);
      const size_t newsize = newRowStateVectors.size();
      if (oldsize != newsize)
        dserror("stateVector sizes are different! Fix this!");

      //        cout << "old phi is " << *flamefront_->Phin();
      //        cout << "new phi is " << *flamefront_->Phinp();

      //        cout << "vel before setting new startvalues " << *stateVector[2];

      //        cout << "in algorithm old phi is " << *flamefront_->Phin();
      //        cout << "in algorithm new phi is " << *flamefront_->Phinp();
//TODO activate new input parameter
//      bool start_val_semilagrange = Teuchos::getIntegralValue<int>(params_.sublist("COMBUSTION FLUID"),"START_VAL_SEMILAGRANGE");
//      bool start_val_enrichment   = Teuchos::getIntegralValue<int>(params_.sublist("COMBUSTION FLUID"),"START_VAL_ENRICHMENT");

      //        if( (start_val_semilagrange == true) || (start_val_enrichment == true) )
      //        {
      XFEM::Startvalues startval(
          discret_,
          olddofmanager,
          dofmanager,
          oldColStateVectors,
          newRowStateVectors,
          veln,
          flamefront_,
          interfacehandleN,
          olddofcolmap,
          oldNodalDofColDistrib,
          newdofrowmap,
          state_.nodalDofDistributionMap_);
      //        }

      // initialize system vectors for nodes which changed interface side
//      if(start_val_semilagrange == true)
//        startval.semiLagrangeExtrapolation(flamefront_,dta_);

      //        cout << "here before setting enrichment values" << endl;
//      if(start_val_enrichment == true)
        startval.setEnrichmentValues();

      //        cout << "vel after setting new startvalues " << *stateVector[2];

      OutputToGmsh("start_field_pres","start_field_vel",Step(), Time());
    }
}
  } // end area for dofswitcher and startvalues

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
    discret_->EvaluateDirichletXFEM(eleparams, zeros_, Teuchos::null, Teuchos::null,
        Teuchos::null, dbcmaps_);
    zeros_->PutScalar(0.0); // just in case of change
  }

  neumann_loads_= LINALG::CreateVector(newdofrowmap,true);

  // ---------------------------------
  // Vectors used for solution process
  // ---------------------------------
  residual_     = LINALG::CreateVector(newdofrowmap,true);
  trueresidual_ = Teuchos::null;
  incvel_       = LINALG::CreateVector(newdofrowmap,true);


  // -------------------------------------------------------------------
  // get a vector layout from the discretization for a vector which only
  // contains the velocity dofs and for one vector which only contains
  // pressure degrees of freedom.
  // -------------------------------------------------------------------
  FLD::UTILS::SetupXFluidSplit(*discret_,dofmanager,velpressplitter_);

  // -------------------------------------------------------------------
  // create empty system matrix --- stiffness and mass are assembled in
  // one system matrix!
  // -------------------------------------------------------------------
  cout0_ << " Initialize system matrix..." << endl;

  // initialize system matrix
  sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(newdofrowmap,0,false,true));

}


/*------------------------------------------------------------------------------------------------*
 | import geometrical information about the interface (integration cells) from the combustion     |
 | algorithm and incorporate it into the fluid field                                  henke 03/09 |
 |
 | remark: Within this routine, no parallel re-distribution is allowed to take place. Before and  |
 | after this function, it's ok to do so.                                                    axel |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::StoreFlameFront(const Teuchos::RCP<COMBUST::FlameFront>& flamefront)
{
  flamefront_ = flamefront;
  return;
}


/*------------------------------------------------------------------------------------------------*
 | get convection velocity vector for transfer to scalar transport field              henke 07/09 |
 *------------------------------------------------------------------------------------------------*/
const Teuchos::RCP<Epetra_Vector> FLD::CombustFluidImplicitTimeInt::ConVelnp()
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
  Teuchos::RCP<Epetra_Vector> convel = dofmanagerForOutput_->transformXFEMtoStandardVector(
                                         *state_.velnp_, *standarddofset_,
                                         state_.nodalDofDistributionMap_, outputfields);
  return convel;
}

/*------------------------------------------------------------------------------------------------*
 | get history vector for transfer to scalar transport field                      rasthofer 01/10 |
 | needed for subgrid-velocity                                                                    |
 *------------------------------------------------------------------------------------------------*/
const Teuchos::RCP<Epetra_Vector> FLD::CombustFluidImplicitTimeInt::Hist()
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
  if (timealgo_==INPAR::FLUID::timeint_one_step_theta)
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
  PrepareNonlinearSolve();

  TEUCHOS_FUNC_TIME_MONITOR("   + nonlin. iteration/lin. solve");

  // ---------------------------------------------- nonlinear iteration
  // ------------------------------- stop nonlinear iteration when both
  //                                 increment-norms are below this bound
  const double  ittol     = params_.get<double>("tolerance for nonlin iter");

  //------------------------------ turn adaptive solver tolerance on/off
  const bool   isadapttol    = params_.get<bool>("ADAPTCONV");
  const double adaptolbetter = params_.get<double>("ADAPTCONV_BETTER");

  //const bool fluidrobin = params_.get<bool>("fluidrobin", false);

  int               itnum = 0;
  bool              stopnonliniter = false;

  double dtsolve = 0.0;
  double dtele   = 0.0;

  // out to screen
  PrintTimeStepInfo();
  
  // action for elements
  if (timealgo_!=INPAR::FLUID::timeint_stationary and theta_ < 1.0)
  {
    cout0_ << "* Warning! Works reliably only for Backward Euler time discretization! *" << endl;
  }

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

    printf("+------------+-------------------+--------------+--------------+--------------+--------------+--------------+--------------+\n");
    printf("|- step/max -|- tol      [norm] -|-- vel-res ---|-- pre-res ---|-- fullres ---|-- vel-inc ---|-- pre-inc ---|-- fullinc ---|\n");
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
      if (timealgo_==INPAR::FLUID::timeint_stationary)
        eleparams.set("action","calc_fluid_stationary_systemmat_and_residual");
      else
        eleparams.set("action","calc_fluid_systemmat_and_residual");

      // flag for type of combustion problem
      eleparams.set("combusttype",combusttype_);
      eleparams.set("veljumptype",veljumptype_);
      eleparams.set("normaltensionjumptype",normaltensionjumptype_);
      eleparams.set("flamespeed",flamespeed_);
      eleparams.set("nitschevel",nitschevel_);
      eleparams.set("nitschepres",nitschepres_);
      eleparams.set("DLM_condensation",condensation_);

      // parameters for two-phase flow problems with surface tension
      eleparams.set("surftensapprox",surftensapprox_);
      eleparams.set("surftenscoeff",surftenscoeff_);
      eleparams.set("connected_interface",connected_interface_);

      // other parameters that might be needed by the elements
      //eleparams.set("total time",time_);
      //eleparams.set("thsl",theta_*dta_);
      eleparams.set("timealgo",timealgo_);
      eleparams.set("dt",dta_);
      eleparams.set("theta",theta_);

#ifdef SUGRVEL_OUTPUT
      //eleparams.set("step",step_);
#endif

      //eleparams.set("include reactive terms for linearisation",params_.get<bool>("Use reaction terms for linearisation"));
      //type of linearisation: include reactive terms for linearisation
      if(params_.get<INPAR::FLUID::LinearisationAction>("Linearisation") == INPAR::FLUID::Newton)
        eleparams.set("include reactive terms for linearisation",true);
      else if (params_.get<INPAR::FLUID::LinearisationAction>("Linearisation") == INPAR::FLUID::minimal)
        dserror("LinearisationAction minimal is not defined in the combustion formulation");
      else
        eleparams.set("include reactive terms for linearisation",false);

      // parameters for stabilization
      eleparams.sublist("STABILIZATION") = params_.sublist("STABILIZATION");

      // parameters for stabilization
      eleparams.sublist("TURBULENCE MODEL") = params_.sublist("TURBULENCE MODEL");

      // set vector values needed by elements
      discret_->ClearState();
      discret_->SetState("velnp",state_.velnp_);
      discret_->SetState("veln" ,state_.veln_);
      discret_->SetState("velnm",state_.velnm_);
      discret_->SetState("accn" ,state_.accn_);

      discret_->SetState("velpres nodal iterinc",oldinc_);

      // convergence check at itemax is skipped for speedup if
      // CONVCHECK is set to L_2_norm_without_residual_at_itemax
      if ((itnum != itemax_)
          ||
          (params_.get<string>("CONVCHECK")
           !=
           "L_2_norm_without_residual_at_itemax"))
      {
        // call standard loop over elements
        discret_->Evaluate(eleparams,sysmat_,residual_);

        discret_->ClearState();

        // scaling to get true residual vector for all other schemes
//        trueresidual_->Update(ResidualScaling(),*residual_,0.0);

        // finalize the complete matrix
        sysmat_->Complete();
      }

      // end time measurement for element
      dtele=Teuchos::Time::wallTime()-tcpu;

    } // end of element call

    // blank residual DOFs which are on Dirichlet BC
    // We can do this because the values at the dirichlet positions
    // are not used anyway.
    // We could avoid this though, if velrowmap_ and prerowmap_ would
    // not include the dirichlet values as well. But it is expensive
    // to avoid that.
    dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), residual_);

    double incvelnorm_L2 = 0.0;
    double velnorm_L2 = 0.0;
    double vresnorm = 0.0;

    Teuchos::RCP<Epetra_Vector> onlyvel = velpressplitter_.ExtractOtherVector(residual_);
    onlyvel->Norm2(&vresnorm);

    velpressplitter_.ExtractOtherVector(incvel_,onlyvel);
    onlyvel->Norm2(&incvelnorm_L2);

    velpressplitter_.ExtractOtherVector(state_.velnp_,onlyvel);
    onlyvel->Norm2(&velnorm_L2);

    double incprenorm_L2 = 0.0;
    double prenorm_L2 = 0.0;
    double presnorm = 0.0;

    Teuchos::RCP<Epetra_Vector> onlypre = velpressplitter_.ExtractCondVector(residual_);
    onlypre->Norm2(&presnorm);

    velpressplitter_.ExtractCondVector(incvel_,onlypre);
    onlypre->Norm2(&incprenorm_L2);

    velpressplitter_.ExtractCondVector(state_.velnp_,onlypre);
    onlypre->Norm2(&prenorm_L2);

    double incfullnorm_L2 = 0.0;
    double fullnorm_L2 = 0.0;
    double fullresnorm = 0.0;

    const Epetra_Map* dofrowmap       = discret_->DofRowMap();
    Epetra_Vector full(*dofrowmap);
    Epetra_Import importer(*dofrowmap,residual_->Map());

    int err = full.Import(*residual_,importer,Insert);
    if (err) dserror("Import using importer returned err=%d",err);
    full.Norm2(&fullresnorm);

    err = full.Import(*incvel_,importer,Insert);
    if (err) dserror("Import using importer returned err=%d",err);
    full.Norm2(&incfullnorm_L2);

    err = full.Import(*state_.velnp_,importer,Insert);
    if (err) dserror("Import using importer returned err=%d",err);
    full.Norm2(&fullnorm_L2);

    // care for the case that nothing really happens in the velocity
    // or pressure field
    if (velnorm_L2 < 1e-5) velnorm_L2 = 1.0;
    if (prenorm_L2 < 1e-5) prenorm_L2 = 1.0;
    if (fullnorm_L2 < 1e-5) fullnorm_L2 = 1.0;

    //-------------------------------------------------- output to screen
    /* special case of very first iteration step:
        - solution increment is not yet available
        - convergence check is not required (we solve at least once!)    */
    if (itnum == 1)
    {
      if (myrank_ == 0)
      {
        printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   | %10.3E   |      --      |      --      |      --      |",
               itnum,itemax_,ittol,vresnorm,presnorm,fullresnorm);
        printf(" (      --     ,te=%10.3E",dtele);
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
      if (vresnorm <= ittol and
          presnorm <= ittol and
          fullresnorm <= ittol and
          incvelnorm_L2/velnorm_L2 <= ittol and
          incprenorm_L2/prenorm_L2 <= ittol and
          incfullnorm_L2/fullnorm_L2 <= ittol)
      {
        stopnonliniter=true;
        if (myrank_ == 0)
        {
          printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   | %10.3E   | %10.3E   | %10.3E   | %10.3E   |",
                 itnum,itemax_,ittol,vresnorm,presnorm,fullresnorm,
                 incvelnorm_L2/velnorm_L2,incprenorm_L2/prenorm_L2,incfullnorm_L2/fullnorm_L2);
          printf(" (ts=%10.3E,te=%10.3E",dtsolve,dtele);
          printf(")\n");
          printf("+------------+-------------------+--------------+--------------+--------------+--------------+--------------+--------------+\n");

          FILE* errfile = params_.get<FILE*>("err file");
          if (errfile!=NULL)
          {
            fprintf(errfile,"fluid solve:   %3d/%3d  tol=%10.3E[L_2 ]  vres=%10.3E  pres=%10.3E  vinc=%10.3E  pinc=%10.3E\n",
                    itnum,itemax_,ittol,vresnorm,presnorm,
                    incvelnorm_L2/velnorm_L2,incprenorm_L2/prenorm_L2);
          }
        }
        break;
      }
      else // if not yet converged
        if (myrank_ == 0)
        {
          printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   | %10.3E   | %10.3E   | %10.3E   | %10.3E   |",
                 itnum,itemax_,ittol,vresnorm,presnorm,fullresnorm,
                 incvelnorm_L2/velnorm_L2,incprenorm_L2/prenorm_L2,incfullnorm_L2/fullnorm_L2);
          printf(" (ts=%10.3E,te=%10.3E",dtsolve,dtele);
          printf(")\n");
        }
    }

    // warn if itemax is reached without convergence, but proceed to
    // next timestep...
    if ((itnum == itemax_) and (vresnorm > ittol or presnorm > ittol or
                             fullresnorm > ittol or
                             incvelnorm_L2/velnorm_L2 > ittol or
                             incprenorm_L2/prenorm_L2 > ittol or
                             incfullnorm_L2/fullnorm_L2 > ittol))
    {
      stopnonliniter=true;
      if (myrank_ == 0)
      {
        printf("+---------------------------------------------------------------+\n");
        printf("|            >>>>>> not converged in itemax steps!              |\n");
        printf("+---------------------------------------------------------------+\n");

        FILE* errfile = params_.get<FILE*>("err file");
        if (errfile!=NULL)
        {
          fprintf(errfile,"fluid unconverged solve:   %3d/%3d  tol=%10.3E[L_2 ]  vres=%10.3E  pres=%10.3E  vinc=%10.3E  pinc=%10.3E\n",
                  itnum,itemax_,ittol,vresnorm,presnorm,
                  incvelnorm_L2/velnorm_L2,incprenorm_L2/prenorm_L2);
        }
      }
      break;
    }

    // stop if NaNs occur
    if (std::isnan(vresnorm) or
        std::isnan(presnorm) or
        std::isnan(fullresnorm) or
        std::isnan(incvelnorm_L2/velnorm_L2) or
        std::isnan(incprenorm_L2/prenorm_L2) or
        std::isnan(incfullnorm_L2/fullnorm_L2))
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
        solver_.AdaptTolerance(ittol,currresidual,adaptolbetter);
      }
      solver_.Solve(sysmat_->EpetraOperator(),incvel_,residual_,true,itnum==1);
      solver_.ResetTolerance();

      // end time measurement for solver
      // remark: to get realistic times, this barrier results in the longest time being printed
      discret_->Comm().Barrier();
      dtsolve = Teuchos::Time::wallTime()-tcpusolve;
    }

    //------------------------------------------------ update (u,p) trial
    state_.velnp_->Update(1.0,*incvel_,1.0);
    //------------------- store nodal increment for element stress update
    oldinc_->Update(1.0,*incvel_,0.0);
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
  // compute error norms // schott Aug 6, 2010
  //--------------------
  INPAR::COMBUST::NitscheError errortype = Teuchos::getIntegralValue<INPAR::COMBUST::NitscheError>(params_.sublist("COMBUSTION FLUID"),"NITSCHE_ERROR");
  if(errortype != INPAR::COMBUST::nitsche_error_none)
    FLD::CombustFluidImplicitTimeInt::EvaluateErrorComparedToAnalyticalSol_Nitsche(errortype);

  oldinc_ = Teuchos::null;
  incvel_ = Teuchos::null;
  residual_ = Teuchos::null;
  zeros_ = Teuchos::null;
  sysmat_ = Teuchos::null;

} // CombustImplicitTimeInt::NonlinearSolve


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

  // TODO @Florian: copied from Axel, needed here?
  oldinc_= Teuchos::null;

  return;
}// FluidImplicitTimeInt::TimeUpdate

/*------------------------------------------------------------------------------------------------*
 | lift'n'drag forces, statistics time sample and output of solution and statistics   gammi 11/08 |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::StatisticsAndOutput()
{
  // time measurement: output and statistics
  TEUCHOS_FUNC_TIME_MONITOR("      + output and statistics");

  // -------------------------------------------------------------------
  //          calculate lift'n'drag forces from the residual
  // -------------------------------------------------------------------
  LiftDrag();

  // -------------------------------------------------------------------
  //          calculate flow through surfaces
  // -------------------------------------------------------------------
  ComputeSurfaceFlowrates();

  // -------------------------------------------------------------------
  //   add calculated velocity to mean value calculation (statistics)
  // -------------------------------------------------------------------
//  statisticsmanager_->DoTimeSample(step_,time_);

  // -------------------------------------------------------------------
  //                         output of solution
  // -------------------------------------------------------------------
  Output();

  // -------------------------------------------------------------------
  //          dumping of turbulence statistics if required
  // -------------------------------------------------------------------
//  statisticsmanager_->DoOutput(step_);

  return;
} // CombustFluidImplicitTimeInt::StatisticsAndOutput

/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::Output()
{


  const bool write_visualization_data = step_%upres_ == 0;
  const bool write_restart_data       = uprestart_ != 0 and step_%uprestart_ == 0;

  //-------------------------------------------- output of solution

  if (write_visualization_data or write_restart_data)
  {
    output_->NewStep(step_,time_);
  }

  if (write_visualization_data)  //write solution for visualization
  {
    // transform velocity XFEM vector to (standard FEM) output velocity vector
    Teuchos::RCP<Epetra_Vector> velnp_out = dofmanagerForOutput_->transformXFEMtoStandardVector(
        *state_.velnp_, *standarddofset_, state_.nodalDofDistributionMap_, physprob_.xfemfieldset_);

    // write physical fields on full domain including voids etc.
    if (physprob_.xfemfieldset_.find(XFEM::PHYSICS::Velx) != physprob_.xfemfieldset_.end())
    {
      // output velocity field for visualization
      output_->WriteVector("velocity_smoothed", velnp_out);

      // output (hydrodynamic) pressure for visualization
      Teuchos::RCP<Epetra_Vector> pressure = velpressplitterForOutput_.ExtractCondVector(velnp_out);
//      pressure->Scale(density_);
      output_->WriteVector("pressure_smoothed", pressure);

      //output_->WriteVector("residual", trueresidual_);

      //only perform stress calculation when output is needed
      if (writestresses_)
      {
        Teuchos::RCP<Epetra_Vector> traction = CalcStresses();
        output_->WriteVector("traction",traction);
      }
    }
    else if (physprob_.xfemfieldset_.find(XFEM::PHYSICS::Temp) != physprob_.xfemfieldset_.end())
    {
      output_->WriteVector("temperature_smoothed", velnp_out);
    }

    // write domain decomposition for visualization
    output_->WriteElementData();
  }

  // write restart
  if (write_restart_data)
  {
   std::cout << "Write restart" << std::endl;
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
  }

//  if (discret_->Comm().NumProc() == 1)
  if (step_ % 1 == 0 or step_== 1) //write every 5th time step only
  {
    OutputToGmsh("solution_field_pressure","solution_field_velocity",step_, time_);
  }

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
//    // remark: pressure scaling was removed in COMBUST;
//    //         we always compute the real (hydrodynamic) pressure [N/m^2] (not p/density!)
//    // pressure->Scale(density_);
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

  std::cout << "Read restart" << std::endl;

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

}

/*------------------------------------------------------------------------------------------------*
 | write output to Gmsh postprocessing files                                          henke 10/09 |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::OutputToGmsh(
    char* presName,
    char* velName,
    const int step,
    const double time
    ) const
{
  const Teuchos::ParameterList& xfemparams = DRT::Problem::Instance()->XFEMGeneralParams();
  const bool gmshdebugout = (bool)getIntegralValue<int>(xfemparams,"GMSH_DEBUG_OUT");

  const bool screen_out = true;

  // get a copy on columnmn parallel distribution
  Teuchos::RCP<const Epetra_Vector> output_col_velnp = DRT::UTILS::GetColVersionOfRowVector(discret_, state_.velnp_);


  if (gmshdebugout and (this->physprob_.xfemfieldset_.find(XFEM::PHYSICS::Pres) != this->physprob_.xfemfieldset_.end()))
  {
    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles(presName, step, 500, screen_out, discret_->Comm().MyPID());
    std::ofstream gmshfilecontent(filename.c_str());

    const XFEM::PHYSICS::Field field = XFEM::PHYSICS::Pres;
    {
      gmshfilecontent << "View \" " << "Pressure Solution (Physical) \" {\n";
      for (int i=0; i<discret_->NumMyRowElements(); ++i)
      {
        const DRT::Element* actele = discret_->lRowElement(i);

        // create local copy of information about dofs
        const XFEM::ElementDofManager eledofman(*actele,physprob_.elementAnsatz_->getElementAnsatz(actele->Shape()),*dofmanagerForOutput_);

        vector<int> lm;
        vector<int> lmowner;
        actele->LocationVector(*(discret_), lm, lmowner);

        // extract local values from the global vector
        vector<double> myvelnp(lm.size());
        DRT::UTILS::ExtractMyValues(*output_col_velnp, myvelnp, lm);

        const int numparam = eledofman.NumDofPerField(field);
        const vector<int>& dofpos = eledofman.LocalDofPosPerField(field);

        LINALG::SerialDenseMatrix elementvalues(1,numparam);
        for (int iparam=0; iparam<numparam; ++iparam)
          elementvalues(0,iparam) = myvelnp[dofpos[iparam]];

        //------------------------------------------------------------------------------------------
        // extract local level-set (G-function) values from global vector
        //------------------------------------------------------------------------------------------
        // get pointer to vector holding G-function values at the fluid nodes
        const Teuchos::RCP<Epetra_Vector> phinp = flamefront_->Phinp();
#ifdef DEBUG
        // get map of this vector
        const Epetra_BlockMap& phimap = phinp->Map();
        // check, whether this map is still identical with the current node map in the discretization
        if (not phimap.SameAs(*discret_->NodeColMap())) dserror("node column map has changed!");
#endif

        size_t numnode = actele->NumNode();
        vector<double> myphinp(numnode);
        // extract G-function values to element level
        DRT::UTILS::ExtractMyNodeBasedValues(actele, myphinp, *phinp);
        const GEO::DomainIntCells& domainintcells =
          dofmanagerForOutput_->getInterfaceHandle()->GetDomainIntCells(actele);

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
            XFEM::InterpolateCellValuesFromElementValuesLevelSet(*actele, eledofman, *cell, myphinp, field,
              elementvalues, cellvalues);
            break;
          }
          case INPAR::COMBUST::combusttype_twophaseflow:
          {
            //
            XFEM::InterpolateCellValuesFromElementValuesLevelSetKink(*actele, eledofman, *cell, myphinp, field,
              elementvalues, cellvalues);
            break;
          }
          case INPAR::COMBUST::combusttype_twophaseflow_surf:
          {
            // schott May 17, 2010
            // there is interpolation function for combined enrichments (jumps in pressure and kinks in velocity field)
            XFEM::InterpolateCellValuesFromElementValuesLevelSetKinkJump(*actele, eledofman, *cell, myphinp, field,
              elementvalues, cellvalues);
            break;
          }
          default:
            dserror("unknown type of combustion problem!");
          }
          //
          const size_t numnodes = DRT::UTILS::getNumberOfElementNodes(cell->Shape());
          LINALG::SerialDenseVector cellvaluesvec(numnodes);
          for (size_t inode=0; inode<numnodes; ++inode)
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
  if (gmshdebugout and (this->physprob_.xfemfieldset_.find(XFEM::PHYSICS::Temp) != this->physprob_.xfemfieldset_.end()) )
  {
    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("solution_field_temperature", step, 5, screen_out, discret_->Comm().MyPID());
    std::ofstream gmshfilecontent(filename.c_str());

    {
      const XFEM::PHYSICS::Field field = XFEM::PHYSICS::Temp;
      gmshfilecontent << "View \" " << "Temperature Solution (Physical) \" {\n";
      for (int i=0; i<discret_->NumMyRowElements(); ++i)
      {
        const DRT::Element* actele = discret_->lRowElement(i);

        // create local copy of information about dofs
        const XFEM::ElementDofManager eledofman(*actele,physprob_.elementAnsatz_->getElementAnsatz(actele->Shape()),*dofmanagerForOutput_);

        vector<int> lm;
        vector<int> lmowner;
        actele->LocationVector(*(discret_), lm, lmowner);

        // extract local values from the global vector
        vector<double> myvelnp(lm.size());
        DRT::UTILS::ExtractMyValues(*output_col_velnp, myvelnp, lm);

        const int numparam = eledofman.NumDofPerField(field);
        const vector<int>& dofpos = eledofman.LocalDofPosPerField(field);

        LINALG::SerialDenseMatrix elementvalues(1,numparam);
        for (int iparam=0; iparam<numparam; ++iparam)
          elementvalues(0,iparam) = myvelnp[dofpos[iparam]];

        //------------------------------------------------------------------------------------------
        // extract local level-set (G-function) values from global vector
        //------------------------------------------------------------------------------------------
        // get pointer to vector holding G-function values at the fluid nodes
        const Teuchos::RCP<Epetra_Vector> phinp = flamefront_->Phinp();
#ifdef DEBUG
        // get map of this vector
        const Epetra_BlockMap& phimap = phinp->Map();
        // check, whether this map is still identical with the current node map in the discretization
        if (not phimap.SameAs(*discret_->NodeColMap())) dserror("node column map has changed!");
#endif

        size_t numnode = actele->NumNode();
        vector<double> myphinp(numnode);
        // extract G-function values to element level
        DRT::UTILS::ExtractMyNodeBasedValues(actele, myphinp, *phinp);

        const GEO::DomainIntCells& domainintcells =
          dofmanagerForOutput_->getInterfaceHandle()->GetDomainIntCells(actele);
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
            XFEM::InterpolateCellValuesFromElementValuesLevelSet(*actele, eledofman, *cell, myphinp, field,
              elementvalues, cellvalues);
            break;
          }
          case INPAR::COMBUST::combusttype_twophaseflow:
          {
            //
            XFEM::InterpolateCellValuesFromElementValuesLevelSetKink(*actele, eledofman, *cell, myphinp, field,
              elementvalues, cellvalues);
            break;
          }
          case INPAR::COMBUST::combusttype_twophaseflow_surf:
          {
            // schott May 17, 2010
            // there is interpolation function for combined enrichments (jumps in pressure and kinks in velocity field)
            XFEM::InterpolateCellValuesFromElementValuesLevelSetKinkJump(*actele, eledofman, *cell, myphinp, field,
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
#if 0
  if (gmshdebugout)
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
        const DRT::Element* actele = discret_->lRowElement(i);

        static LINALG::Matrix<3,27> xyze_xfemElement;
        GEO::fillInitialPositionArray(actele,xyze_xfemElement);

        const map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> element_ansatz(COMBUST::getElementAnsatz(actele->Shape()));

        // create local copy of information about dofs
        const XFEM::ElementDofManager eledofman(*actele,element_ansatz,*dofmanagerForOutput_);

        vector<int> lm;
        vector<int> lmowner;
        actele->LocationVector(*(discret_), lm, lmowner);

        // extract local values from the global vector
        vector<double> myvelnp(lm.size());
        DRT::UTILS::ExtractMyValues(*output_col_velnp, myvelnp, lm);

        const int numparam = eledofman.NumDofPerField(field);
        const vector<int>& dofpos = eledofman.LocalDofPosPerField(field);

        LINALG::SerialDenseVector elementvalues(numparam);
        for (int iparam=0; iparam<numparam; ++iparam)
          elementvalues(iparam) = myvelnp[dofpos[iparam]];

        const GEO::DomainIntCells& domainintcells =
          dofmanagerForOutput_->getInterfaceHandle()->GetDomainIntCells(actele);
        for (GEO::DomainIntCells::const_iterator cell =
          domainintcells.begin(); cell != domainintcells.end(); ++cell)
        {
          LINALG::SerialDenseVector cellvalues(DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
          XFEM::computeScalarCellNodeValuesFromElementUnknowns(*actele, dofmanagerForOutput_->getInterfaceHandle(), eledofman,
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
  if (gmshdebugout)
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
        const DRT::Element* actele = discret_->lRowElement(i);

        // create local copy of information about dofs
        const XFEM::ElementDofManager eledofman(*actele,physprob_.elementAnsatz_->getElementAnsatz(actele->Shape()),*dofmanagerForOutput_);

        vector<int> lm;
        vector<int> lmowner;
        actele->LocationVector(*(discret_), lm, lmowner);

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
          dofmanagerForOutput_->getInterfaceHandle()->GetDomainIntCells(actele);
        for (GEO::DomainIntCells::const_iterator cell =
          domainintcells.begin(); cell != domainintcells.end(); ++cell)
        {
          //LINALG::SerialDenseMatrix xyze_cell(3, cell->NumNode());
          //cell->NodalPosXYZ(*actele, xyze_cell);
          // TODO remove
          const LINALG::SerialDenseMatrix& xyze_cell = cell->CellNodalPosXYZ();

          {
          LINALG::SerialDenseMatrix cellvalues(9,DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
          XFEM::computeTensorCellNodeValuesFromElementUnknowns(*actele, dofmanagerForOutput_->getInterfaceHandle(), eledofman,
              *cell, field, elementvalues, cellvalues);
           IO::GMSH::cellWithTensorFieldToStream(cell->Shape(), cellvalues, xyze_cell, gmshfilecontent);
          }

          {
          LINALG::SerialDenseVector cellvaluexx(DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
          XFEM::computeScalarCellNodeValuesFromElementUnknowns(*actele, dofmanagerForOutput_->getInterfaceHandle(), eledofman, *cell, field, elementvaluexx, cellvaluexx);
          IO::GMSH::cellWithScalarFieldToStream(cell->Shape(), cellvaluexx, xyze_cell, gmshfilecontentxx);
          }
          {
          LINALG::SerialDenseVector cellvalueyy(DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
          XFEM::computeScalarCellNodeValuesFromElementUnknowns(*actele, dofmanagerForOutput_->getInterfaceHandle(), eledofman, *cell, field, elementvalueyy, cellvalueyy);
          IO::GMSH::cellWithScalarFieldToStream(cell->Shape(), cellvalueyy, xyze_cell, gmshfilecontentyy);
          }
          {
          LINALG::SerialDenseVector cellvaluezz(DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
          XFEM::computeScalarCellNodeValuesFromElementUnknowns(*actele, dofmanagerForOutput_->getInterfaceHandle(), eledofman, *cell, field, elementvaluezz, cellvaluezz);
          IO::GMSH::cellWithScalarFieldToStream(cell->Shape(), cellvaluezz, xyze_cell, gmshfilecontentzz);
          }
          {
          LINALG::SerialDenseVector cellvaluexy(DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
          XFEM::computeScalarCellNodeValuesFromElementUnknowns(*actele, dofmanagerForOutput_->getInterfaceHandle(), eledofman, *cell, field, elementvaluexy, cellvaluexy);
          IO::GMSH::cellWithScalarFieldToStream(cell->Shape(), cellvaluexy, xyze_cell, gmshfilecontentxy);
          }
          {
          LINALG::SerialDenseVector cellvaluexz(DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
          XFEM::computeScalarCellNodeValuesFromElementUnknowns(*actele, dofmanagerForOutput_->getInterfaceHandle(), eledofman, *cell, field, elementvaluexz, cellvaluexz);
          IO::GMSH::cellWithScalarFieldToStream(cell->Shape(), cellvaluexz, xyze_cell, gmshfilecontentxz);
          }
          {
          LINALG::SerialDenseVector cellvalueyz(DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
          XFEM::computeScalarCellNodeValuesFromElementUnknowns(*actele, dofmanagerForOutput_->getInterfaceHandle(), eledofman, *cell, field, elementvalueyz, cellvalueyz);
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

  if (this->physprob_.xfemfieldset_.find(XFEM::PHYSICS::Velx) != this->physprob_.xfemfieldset_.end())
  {
    PlotVectorFieldToGmsh(DRT::UTILS::GetColVersionOfRowVector(discret_, state_.velnp_), velName,"Velocity Solution (Physical) n+1",true, step, time);
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
  const Teuchos::ParameterList& xfemparams = DRT::Problem::Instance()->XFEMGeneralParams();
  const bool gmshdebugout = (bool)getIntegralValue<int>(xfemparams,"GMSH_DEBUG_OUT");

  const bool screen_out = true;

  if (gmshdebugout)
  {
    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles(filestr, step, 500, screen_out, discret_->Comm().MyPID());
    std::ofstream gmshfilecontent(filename.c_str());

    {
      gmshfilecontent << "View \" " << name_in_gmsh << "\" {\n";
      for (int i=0; i<discret_->NumMyRowElements(); ++i)
      {
        const DRT::Element* actele = discret_->lRowElement(i);

        // create local copy of information about dofs
        const COMBUST::TauPressureAnsatz elementAnsatz;
        const XFEM::ElementDofManager eledofman(*actele,elementAnsatz.getElementAnsatz(actele->Shape()),*dofmanagerForOutput_);

        vector<int> lm;
        vector<int> lmowner;
        actele->LocationVector(*discret_, lm, lmowner);

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

        //------------------------------------------------------------------------------------------
        // extract local level-set (G-function) values from global vector
        //------------------------------------------------------------------------------------------
        // get pointer to vector holding G-function values at the fluid nodes
        const Teuchos::RCP<Epetra_Vector> phinp = flamefront_->Phinp();
#ifdef DEBUG
        // get map of this vector
        const Epetra_BlockMap& phimap = phinp->Map();
        // check, whether this map is still identical with the current node map in the discretization
        if (not phimap.SameAs(*discret_->NodeColMap())) dserror("node column map has changed!");
#endif

        size_t numnode = actele->NumNode();
        vector<double> myphinp(numnode);
        // extract G-function values to element level
        DRT::UTILS::ExtractMyNodeBasedValues(actele, myphinp, *phinp);

        const GEO::DomainIntCells& domainintcells =
          dofmanagerForOutput_->getInterfaceHandle()->GetDomainIntCells(actele);
        for (GEO::DomainIntCells::const_iterator cell =
          domainintcells.begin(); cell != domainintcells.end(); ++cell)
        {
          LINALG::SerialDenseMatrix cellvalues(3, DRT::UTILS::getNumberOfElementNodes(cell->Shape()));

          const XFEM::PHYSICS::Field field = XFEM::PHYSICS::Velx;
          switch(combusttype_)
          {
          case INPAR::COMBUST::combusttype_premixedcombustion:
          case INPAR::COMBUST::combusttype_twophaseflowjump:
          {
            //
            XFEM::InterpolateCellValuesFromElementValuesLevelSet(*actele, eledofman, *cell, myphinp, field,
              elementvalues, cellvalues);
            break;
          }
          case INPAR::COMBUST::combusttype_twophaseflow:
          {
            //
            XFEM::InterpolateCellValuesFromElementValuesLevelSetKink(*actele, eledofman, *cell, myphinp, field,
              elementvalues, cellvalues);
            break;
          }
          case INPAR::COMBUST::combusttype_twophaseflow_surf:
          {
            // plots for Velx, Vely, Velz ... (Kink enrichment)
            XFEM::InterpolateCellValuesFromElementValuesLevelSetKinkJump(*actele, eledofman, *cell, myphinp, field,
              elementvalues, cellvalues);
            break;
          }
          default:
            dserror("unknown type of combustion problem!");
          }

          IO::GMSH::cellWithVectorFieldToStream(
              cell->Shape(), cellvalues, cell->CellNodalPosXYZ(), gmshfilecontent);
        }
//          const GEO::BoundaryIntCells& boundaryintcells =
//              dofmanagerForOutput_->getInterfaceHandle()->GetBoundaryIntCells(actele->Id());
//          // draw boundary integration cells with values
//          for (GEO::BoundaryIntCells::const_iterator cell =
//            boundaryintcells.begin(); cell != boundaryintcells.end(); ++cell)
//          {
//            LINALG::SerialDenseMatrix cellvalues(3, DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
//
//            // the 'label' is set to -1 for combustion problems   henke 10/09
//            XFEM::computeVectorCellNodeValues(*actele, ih_np_, eledofman,
//                *cell, XFEM::PHYSICS::Velx, -1, elementvalues, cellvalues);
//            IO::GMSH::cellWithVectorFieldToStream(
//                cell->Shape(), cellvalues, cell->CellNodalPosXYZ(), gmshfilecontent);
//          }

        // draw uncut element
//        {
//          LINALG::SerialDenseMatrix elevalues(3, DRT::UTILS::getNumberOfElementNodes(actele->Shape()),true);
//          static LINALG::Matrix<3,27> xyze_ele;
//          GEO::fillInitialPositionArray(actele, xyze_ele);
//          IO::GMSH::cellWithVectorFieldToStream(
//                          actele->Shape(), elevalues, xyze_ele, gmshfilecontent);
//        }

//        }
        //if (dofmanagerForOutput_->getInterfaceHandle()->ElementIntersected(elegid) and not ele_to_textfile and ele_to_textfile2)
//        if (elegid == 14 and elementvalues.N() > 0 and plot_to_gnuplot)
        if (actele->Id() == 1 and elementvalues.N() > 0 and plot_to_gnuplot)
        {
          //std::cout << elementvalues << std::endl;
          std::ofstream f;
          const std::string fname = DRT::Problem::Instance()->OutputControlFile()->FileName()
                                  + ".outflowvel.txt";
          if (step <= 1)
            f.open(fname.c_str(),std::fstream::trunc);
          else
            f.open(fname.c_str(),std::fstream::ate | std::fstream::app);

          //f << time_ << " " << (-1.5*std::sin(0.1*2.0*time_* M_PI) * M_PI*0.1) << "  " << elementvalues(0,0) << endl;
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
        err += state_.velnp_->ReplaceMyValues(1,&(u[nveldof]),&lid);
        err += state_.veln_ ->ReplaceMyValues(1,&(u[nveldof]),&lid);
        err += state_.velnm_->ReplaceMyValues(1,&(u[nveldof]),&lid);
     }

      // initial pressure
      const int gid = nodedofset[npredof];
      int lid = dofrowmap->LID(gid);
      err += state_.velnp_->ReplaceMyValues(1,&p,&lid);
      err += state_.veln_ ->ReplaceMyValues(1,&p,&lid);
      err += state_.velnm_->ReplaceMyValues(1,&p,&lid);

    } // end loop nodes lnodeid
    if(err!=0)
    {
      dserror("dof not on proc");
    }
    break;
  }
  case INPAR::COMBUST::initfield_field_by_function:
  case INPAR::COMBUST::initfield_disturbed_field_by_function:
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

        double initialval=DRT::Problem::Instance()->Funct(initfuncno-1).Evaluate(index,lnode->X(),0.0,NULL);

        state_.velnp_->ReplaceGlobalValues(1,&initialval,&gid);
        state_.veln_ ->ReplaceGlobalValues(1,&initialval,&gid);
      }
    }

    // add random perturbation
    if(initfield == INPAR::COMBUST::initfield_disturbed_field_by_function)
    {
      const int numdim = params_.get<int>("number of velocity degrees of freedom");

      const Epetra_Map* dofrowmap = discret_->DofRowMap();

      int err =0;

      // random noise is perc percent of the initial profile

      double perc = params_.sublist("TURBULENCE MODEL").get<double>("CHAN_AMPL_INIT_DIST");

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
    xyz0_left(0)  = 37.5; // x-coordinate left vortex
    xyz0_left(1)  = 75.0; // y-coordinate left vortex
    xyz0_left(2)  = 0.0;  // z-coordinate is 0 (2D problem)
    xyz0_right(0) = 62.5; // x-coordinate right vortex
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

    // get G-function value vector on fluid NodeColMap
    const Teuchos::RCP<Epetra_Vector> phinp = flamefront_->Phinp();

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
      const int lid = phinp->Map().LID(lnode->Id());
      const double gfuncval = (*phinp)[lid];

//cout << "G-function value " << gfuncval << endl;

      //----------------------------------------
      // set density with respect to flame front
      //----------------------------------------
      if (gfuncval > 0.0) // plus/burnt domain -> burnt material
      {
        dens = dens_b;
        pres = 0.0;
      }
      else // minus/unburnt domain -> unburnt material
      {
        dens = dens_u;
        // TODO @Florian check
        pres = -flamespeed_*flamespeed_*dens_u*dens_u*(1.0/dens_b - 1.0/dens_u);
      }
      //----------------------------------------------
      // compute components of initial velocity vector
      //----------------------------------------------
      // compute preliminary values for both vortices
      double r_squared_left  = ((xyz(0)-xyz0_left(0))*(xyz(0)-xyz0_left(0))
                               +(xyz(1)-xyz0_left(1))*(xyz(1)-xyz0_left(1)))/R_squared;
      double r_squared_right = ((xyz(0)-xyz0_right(0))*(xyz(0)-xyz0_right(0))
                               +(xyz(1)-xyz0_right(1))*(xyz(1)-xyz0_right(1)))/R_squared;

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
        err += state_.velnp_->ReplaceMyValues(1,&vel(idim),&lid);
        err += state_.veln_ ->ReplaceMyValues(1,&vel(idim),&lid);
        err += state_.velnm_->ReplaceMyValues(1,&vel(idim),&lid);
        if(idim==3)
        {
          err += state_.velnp_->ReplaceMyValues(1,&pres,&lid);
          err += state_.veln_ ->ReplaceMyValues(1,&pres,&lid);
          err += state_.velnm_->ReplaceMyValues(1,&pres,&lid);
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
    const Epetra_Map newdofrowmap)
{
  // get G-function value vector on fluid NodeColMap
  const Teuchos::RCP<Epetra_Vector> phinp = flamefront_->Phinp();
  
  
  // initial field modification for flame_vortex_interaction
//  for (int nodeid=0;nodeid<discret_->NumMyRowNodes();nodeid++) // loop over element nodes
//  {
//    DRT::Node* lnode = discret_->lRowNode(nodeid);
//    const std::set<XFEM::FieldEnr>& fieldenrset(dofmanager->getNodeDofSet(lnode->Id()));
//    for (set<XFEM::FieldEnr>::const_iterator fieldenr = fieldenrset.begin();
//        fieldenr != fieldenrset.end();++fieldenr)
//    {
//      const XFEM::DofKey<XFEM::onNode> newdofkey(lnode->Id(), *fieldenr);
//      const int newdofpos = state_.nodalDofDistributionMap_.find(newdofkey)->second;
//      if (fieldenr->getEnrichment().Type() == XFEM::Enrichment::typeJump)
//      {
//        if (fieldenr->getField() == XFEM::PHYSICS::Velx)
//          /* nothing to do in flame_vortex example*/;
//        else if (fieldenr->getField() == XFEM::PHYSICS::Vely)
//        {
//          (*state_.veln_)[newdofrowmap.LID(newdofpos)] = 3.19745;
//          (*state_.velnp_)[newdofrowmap.LID(newdofpos)] = 3.19745;
//        }
//        else if (fieldenr->getField() == XFEM::PHYSICS::Velz)
//          /* nothing to do in flame_vortex example*/;
//        else if (fieldenr->getField() == XFEM::PHYSICS::Pres)
//        {
//          (*state_.veln_)[newdofrowmap.LID(newdofpos)] = 3.712242;
//          (*state_.velnp_)[newdofrowmap.LID(newdofpos)] = 3.712242;
//        }
//      } // end if jump enrichment
//    } // end loop over fieldenr
//  } // end loop over element nodes
  
  
  // initial field modification for collapse_flame
  const int nsd = 3;
  for (int nodeid=0;nodeid<discret_->NumMyRowNodes();nodeid++) // loop over element nodes
  {
    DRT::Node* lnode = discret_->lRowNode(nodeid);
    
    LINALG::Matrix<nsd,1> coords(lnode->X());
//    coords(1) -= 0.5;
    double coordsnorm = sqrt(coords(0)*coords(0)+coords(1)*coords(1));
    
    const int lid = phinp->Map().LID(lnode->Id());
    const double gfuncval = (*phinp)[lid];
    
    const std::set<XFEM::FieldEnr>& fieldenrset(dofmanager->getNodeDofSet(lnode->Id()));
    for (set<XFEM::FieldEnr>::const_iterator fieldenr = fieldenrset.begin();
        fieldenr != fieldenrset.end();++fieldenr)
    {
      const XFEM::DofKey<XFEM::onNode> newdofkey(lnode->Id(), *fieldenr);
      const int newdofpos = state_.nodalDofDistributionMap_.find(newdofkey)->second;
      if (fieldenr->getEnrichment().Type() == XFEM::Enrichment::typeJump)
      {
        if (fieldenr->getField() == XFEM::PHYSICS::Velx)
        {
          // halbe Sprunghoehe von 1.0
          (*state_.veln_)[newdofrowmap.LID(newdofpos)] = coords(0)/(2*coordsnorm);
          (*state_.velnp_)[newdofrowmap.LID(newdofpos)] = coords(0)/(2*coordsnorm);
        }
        else if (fieldenr->getField() == XFEM::PHYSICS::Vely)
        {
          (*state_.veln_)[newdofrowmap.LID(newdofpos)] = coords(1)/(2*coordsnorm);
          (*state_.velnp_)[newdofrowmap.LID(newdofpos)] = coords(1)/(2*coordsnorm);
        }
        else if (fieldenr->getField() == XFEM::PHYSICS::Velz)
        {
          (*state_.veln_)[newdofrowmap.LID(newdofpos)] = 0;//coords(2)/(2*coordsnorm);
          (*state_.velnp_)[newdofrowmap.LID(newdofpos)] = 0;//coords(2)/(2*coordsnorm);
        }
        else if (fieldenr->getField() == XFEM::PHYSICS::Pres)
        {
          (*state_.veln_)[newdofrowmap.LID(newdofpos)] = -2.25;
          (*state_.velnp_)[newdofrowmap.LID(newdofpos)] = -2.25;
        }
      } // end if jump enrichment
      else if (fieldenr->getEnrichment().Type() == XFEM::Enrichment::typeStandard)
      {
        if (gfuncval>=0)
        {
          if (fieldenr->getField() == XFEM::PHYSICS::Velx)
          {
            (*state_.veln_)[newdofrowmap.LID(newdofpos)] = sqrt(0.1)*coords(0)/(coordsnorm*coordsnorm);
            (*state_.velnp_)[newdofrowmap.LID(newdofpos)] = sqrt(0.1)*coords(0)/(coordsnorm*coordsnorm);
          }
          else if (fieldenr->getField() == XFEM::PHYSICS::Vely)
          {
            (*state_.veln_)[newdofrowmap.LID(newdofpos)] = sqrt(0.1)*coords(1)/(coordsnorm*coordsnorm);
            (*state_.velnp_)[newdofrowmap.LID(newdofpos)] = sqrt(0.1)*coords(1)/(coordsnorm*coordsnorm);
          }
          else if (fieldenr->getField() == XFEM::PHYSICS::Velz)
          {
            (*state_.veln_)[newdofrowmap.LID(newdofpos)] = 0;//sqrt(0.1)*coords(2)/(coordsnorm*coordsnorm);
            (*state_.velnp_)[newdofrowmap.LID(newdofpos)] = 0;//sqrt(0.1)*coords(2)/(coordsnorm*coordsnorm);
          }
          else if (fieldenr->getField() == XFEM::PHYSICS::Pres)
          {
            (*state_.veln_)[newdofrowmap.LID(newdofpos)] = -2.7*sqrt(0.1)*sqrt(0.1)/(coordsnorm*coordsnorm);
            (*state_.velnp_)[newdofrowmap.LID(newdofpos)] = -2.7*sqrt(0.1)*sqrt(0.1)/(coordsnorm*coordsnorm);
          }
        }
        else
        {
          if (fieldenr->getField() == XFEM::PHYSICS::Velx)
          {
            (*state_.veln_)[newdofrowmap.LID(newdofpos)] = 0.0;
            (*state_.velnp_)[newdofrowmap.LID(newdofpos)] = 0.0;
          }
          if (fieldenr->getField() == XFEM::PHYSICS::Vely)
          {
            (*state_.veln_)[newdofrowmap.LID(newdofpos)] = 0.0;
            (*state_.velnp_)[newdofrowmap.LID(newdofpos)] = 0.0;
          }
          else if (fieldenr->getField() == XFEM::PHYSICS::Velz)
          {
            (*state_.veln_)[newdofrowmap.LID(newdofpos)] = 0.0;
            (*state_.velnp_)[newdofrowmap.LID(newdofpos)] = 0.0;
          }
          else if (fieldenr->getField() == XFEM::PHYSICS::Pres)
          {
            (*state_.veln_)[newdofrowmap.LID(newdofpos)] = 1.8;
            (*state_.velnp_)[newdofrowmap.LID(newdofpos)] = 1.8;
          }
        }
      }
    } // end loop over fieldenr
  } // end loop over element nodes
  OutputToGmsh("mod_start_field_pres","mod_start_field_vel",Step(), Time());
  
  
}


/*------------------------------------------------------------------------------------------------*
 | compute mesh dependent error norms for Nitsche's method                           schott 05/10 |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::EvaluateErrorComparedToAnalyticalSol_Nitsche(INPAR::COMBUST::NitscheError& NitscheErrorType)
{
  // schott Jun 9, 2010
  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("velnp",state_.velnp_);
  discret_->SetState("veln" ,state_.veln_);
  discret_->SetState("velnm",state_.velnm_);
  discret_->SetState("accn" ,state_.accn_);
  discret_->SetState("nodal increment",oldinc_);

  // create parameters for discretization
  ParameterList eleparams;

  // schott May 26, 2010
  // a new ActionType in combust3.H:	"calc_nitsche_error"
  eleparams.set("action", "calc_nitsche_error");
  eleparams.set("Nitsche_Compare_Analyt", NitscheErrorType);
  // switch different test cases -> set "flowproblem" for elements


  // set parameters for parts of the whole Nitsche-error (mesh-dependent norms), here norms not square rooted
  eleparams.set<double>("L2 integrated velocity domain error", 0.0);
  eleparams.set<double>("L2 integrated grad_velocity domain error", 0.0);
  eleparams.set<double>("H-1/2 integrated viscosity interface error", 0.0);
  eleparams.set<double>("H1/2 integrated velocity jump interface error", 0.0);
  eleparams.set<double>("L2 integrated pressure domain error", 0.0);
  eleparams.set<double>("L2 integrated grad_pressure domain error", 0.0);
  eleparams.set<double>("L2 integrated weighted pressure domain error", 0.0);
  eleparams.set<double>("Nitsche integrated error", 0.0); // the whole Nitsche error

  // TODO @Benedikt: check wether H1 pressure norm is needed ????

  // TODO @Benedikt: relative error norms!!!

  // call loop over elements (but do not assemble anything)
  // discret_->Evaluate calls combust3_evaluate for each element
  discret_->Evaluate(eleparams,null,null,null,null,null); // only actele->Evaluate is called (see combust3_evaluate.cpp)

  discret_->ClearState();

  double locVelDomErr 		= eleparams.get<double>("L2 integrated velocity domain error");
  double locGradVelDomErr 	= eleparams.get<double>("L2 integrated grad_velocity domain error");
  double locViscInterfErr 	= eleparams.get<double>("H-1/2 integrated viscosity interface error");
  double locVelJumpInterfErr 	= eleparams.get<double>("H1/2 integrated velocity jump interface error");
  double locPresDomErr		= eleparams.get<double>("L2 integrated pressure domain error");
  double locGradPresDomErr    = eleparams.get<double>("L2 integrated grad_pressure domain error");
  double locWeightPresDomErr	= eleparams.get<double>("L2 integrated weighted pressure domain error");
  double locNitscheErr		= eleparams.get<double>("Nitsche integrated error");

  // initialize global errors
  double VelDomErr 			= 0.0;
  double GradVelDomErr 		= 0.0;
  double ViscInterfErr	 	= 0.0;
  double VelJumpInterfErr 	= 0.0;
  double PresDomErr			= 0.0;
  double GradPresDomErr		= 0.0;
  double WeightPresDomErr		= 0.0;
  double NitscheErr			= 0.0;

  // TODO @Benedikt: check wether Comm() works also for interface integrals ?
  discret_->Comm().SumAll(&locVelDomErr,&VelDomErr,1);			// sum over processors, each list (list of a processor) has length 1
  discret_->Comm().SumAll(&locGradVelDomErr,&GradVelDomErr,1);
  discret_->Comm().SumAll(&locViscInterfErr,&ViscInterfErr,1);
  discret_->Comm().SumAll(&locVelJumpInterfErr,&VelJumpInterfErr,1);
  discret_->Comm().SumAll(&locPresDomErr,&PresDomErr,1);
  discret_->Comm().SumAll(&locGradPresDomErr, &GradPresDomErr,1);
  discret_->Comm().SumAll(&locWeightPresDomErr,&WeightPresDomErr,1);
  discret_->Comm().SumAll(&locNitscheErr,&NitscheErr,1);


  // for the norms, we need the square roots
  VelDomErr 			= sqrt(VelDomErr);
  GradVelDomErr 		= sqrt(GradVelDomErr);
  ViscInterfErr 		= sqrt(ViscInterfErr);
  VelJumpInterfErr 	= sqrt(VelJumpInterfErr);
  PresDomErr 			= sqrt(PresDomErr);
  GradPresDomErr		= sqrt(GradPresDomErr);
  WeightPresDomErr 	= sqrt(WeightPresDomErr);
  NitscheErr 			= sqrt(NitscheErr);

  //TODO: check if this output is after the right processor
  if (myrank_ == 0)
  {
    printf("\n======================================================================="
        "\n======================= absolute Nitsche errors ======================="
        "\n======= compare analytical solution with approximated solution=========");
    printf("\n  || u-u_h ||_L2(Omega)\t\t\t\t\t%15.8e", 						VelDomErr);
    printf("\n  || sqrt(mu)grad(u-u_h) ||_L2(Omega1 U Omega2) \t%15.8e", 		GradVelDomErr);
    printf("\n  || {2mu* E(u-u_h)*n} ||_H-1/2(Interface) \t\t%15.8e",			ViscInterfErr);
    printf("\n  || |[u-u_h]| ||_H1/2(Interface)\t\t\t%15.8e", 				VelJumpInterfErr);
    printf("\n  || p-p_h ||_L2(Omega)\t\t\t\t\t%15.8e",						PresDomErr);
    printf("\n  ||grad(p-p_h) ||_L2(Omega1 U Omega2)\t\t\t%15.8e",			GradPresDomErr);
    printf("\n  || 1/sqrt(mu_max) * (p-p_h) ||_L2(Omega)\t\t%15.8e",			WeightPresDomErr);
    printf("\n ||| (u-u_h, p-p_h) |||_Nitsche(Omega)\t\t\t%15.8e",			NitscheErr);
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
  INPAR::FLUID::InitialField calcerr = params_.get<INPAR::FLUID::InitialField>("eval err for analyt sol");

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
      printf("\n  L2_err for beltrami flow:  velocity %15.8e  pressure %15.8e\n\n",
             velerr,preerr);
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

/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FLD::CombustFluidImplicitTimeInt::CalcStresses()
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

} // CombustFluidImplicitTimeInt::CalcStresses()

/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::LiftDrag() const
{
  // in this map, the results of the lift drag calculation are stored
  RCP<map<int,vector<double> > > liftdragvals;

  FLD::UTILS::LiftDrag(*discret_,*trueresidual_,params_,liftdragvals);

  return;
} // CombustFluidImplicitTimeInt::LiftDrag

/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::ComputeSurfaceFlowrates() const
{

//  const map<int,double> volumeflowratepersurface = FLD::UTILS::ComputeSurfaceFlowrates(*discret_, state_.velnp_);
//
//  if (not volumeflowratepersurface.empty())
//  {
//    cout << "Number of flow rate conditions... " << volumeflowratepersurface.size() << endl;
//  }
//
//  double overall_flowrate = 0.0;
//  std::map<int,double>::const_iterator entry;
//  for(entry = volumeflowratepersurface.begin(); entry != volumeflowratepersurface.end(); ++entry )
//  {
//    const int condID = entry->first;
//    const double value = entry->second;
//    overall_flowrate += value;
//    if (myrank_ == 0)
//    {
//      cout << " - flowrate for label " << condID << ":  " <<  scientific << value << endl;
//    }
//  }
//  if (not volumeflowratepersurface.empty())
//  {
//    cout << " - flowrate over all boundaries: " << overall_flowrate << endl;
//  }

  return;
}


/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FLD::CombustFluidImplicitTimeInt::IntegrateInterfaceShape(std::string condname)
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
  discret_->EvaluateCondition(eleparams,integratedshapefunc,condname);
  discret_->ClearState();

  return integratedshapefunc;
}

/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::UseBlockMatrix(
    Teuchos::RCP<std::set<int> > condelements,
    const LINALG::MultiMapExtractor& domainmaps,
    const LINALG::MultiMapExtractor& rangemaps,
    bool splitmatrix
    )
{
  dserror("not tested");
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
    //shapederivatives_ = mat;
  }
}

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
  case INPAR::FLUID::timeint_gen_alpha :
    return "  Gen. Alpha  ";
    break;
  default :
    dserror("Cannot cope with name %d", term);
    return "";
    break;
  }
}


#endif // CCADISCRET
