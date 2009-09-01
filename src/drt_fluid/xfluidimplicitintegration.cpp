/*!----------------------------------------------------------------------
\file xfluidimplicitintegration.cpp
\brief Control routine for fluid (in)stationary solvers,

     including instationary solvers based on

     o one-step-theta time-integration scheme

     o two-step BDF2 time-integration scheme
       (with potential one-step-theta start algorithm)

     and stationary solver.

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <stdio.h>


#include <Ifpack_IC.h>

#include "xfluidimplicitintegration.H"
#include "time_integration_scheme.H"

#include "../drt_io/io_control.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/linalg_ana.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_function.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_xfem/dof_management.H"
#include "../drt_xfem/dof_distribution_switcher.H"
#include "../drt_xfem/enrichment_utils.H"
#include "../drt_xfem/element_ansatz.H"
//#include "../drt_xfem/load_balancing.H"
#include "../drt_geometry/position_array.H"
#include "fluid_utils.H"
#include "../drt_f3/xfluid3_interpolation.H"
#include "../drt_xdiff3/xdiff3_interpolation.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "../drt_io/io_gmsh.H"
#include <Teuchos_TimeMonitor.hpp>


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FLD::XFluidImplicitTimeInt::XFluidImplicitTimeInt(
    const Teuchos::RCP<DRT::Discretization> actdis,
    const ParameterList&                    params
//    IO::DiscretizationWriter&               output
    ) :
  // call constructor for "nontrivial" objects
  discret_(actdis),
  fluidsolverparams_(DRT::Problem::Instance()->FluidSolverParams()),
  solver_ (rcp(new LINALG::Solver(DRT::Problem::Instance()->FluidSolverParams(),
                                  actdis->Comm(),
                                  DRT::Problem::Instance()->ErrorFile()->Handle()))),
  params_ (params),
  xparams_(params.sublist("XFEM")),
  output_ (rcp(new IO::DiscretizationWriter(actdis))),
  myrank_(actdis->Comm().MyPID()),
  cout0_(actdis->Comm(), std::cout),
  alefluid_(false),
  time_(0.0),
  step_(0),
  stepmax_(params.get<int>   ("max number timesteps")),
  maxtime_(params.get<double>("total time")),
  timealgo_(params.get<FLUID_TIMEINTTYPE>("time int algo")),
  itemax_(params.get<int>("max nonlin iter steps")),
  extrapolationpredictor_(params.get<bool>("do explicit predictor")),
  uprestart_(params.get<int>("write restart every")),
  upres_(params.get<int>("write solution every")),
  writestresses_(params.get<int>("write stresses")),
  surfacesplitter_(NULL)
{

  // time measurement: initialization
  TEUCHOS_FUNC_TIME_MONITOR(" + initialization");

  cout << "FLD::XFluidImplicitTimeInt::XFluidImplicitTimeInt()" << endl;

  // -------------------------------------------------------------------
  // get the basic parameters first
  // -------------------------------------------------------------------
  // time-step size
  dtp_ = dta_ = params_.get<double>("time step size");
  if (timealgo_ == timeint_stationary and params_.get<double>("time step size") != 1.0)
    dserror("Timestep size (delta t) has to be 1.0 for stationary computations!");

  // parameter theta for time-integration schemes
  theta_    = params_.get<double>("theta");
  if (timealgo_ == timeint_stationary)
    theta_ = 1.0;
  // af-generalized-alpha parameters: gamma_ = 0.5 + alphaM_ - alphaF_
  // (may be reset below when starting algorithm is used)
  alphaM_   = params_.get<double>("alpha_M");
  alphaF_   = params_.get<double>("alpha_F");
  gamma_    = params_.get<double>("gamma");

  // create empty cutter discretization
#if 0
  vector<string> conditions;
  Teuchos::RCP<DRT::Discretization> emptyboundarydis = DRT::UTILS::CreateDiscretizationFromCondition(
      actdis, "dummydiscretizationame", "empty boundary", "BELE3", conditions);
  Teuchos::RCP<Epetra_Vector> tmpdisp1 = LINALG::CreateVector(*emptyboundarydis->DofRowMap(),true);
  Teuchos::RCP<Epetra_Vector> tmpdisp2 = LINALG::CreateVector(*emptyboundarydis->DofRowMap(),true);
  emptyboundarydis->SetState("idispcolnp",tmpdisp1);
  emptyboundarydis->SetState("idispcoln",tmpdisp2);
  // intersection with empty cutter will result in a complete fluid domain with no holes or intersections
  Teuchos::RCP<XFEM::InterfaceHandleXFSI> ih =
    rcp(new XFEM::InterfaceHandleXFSI(discret_,emptyboundarydis));

  // apply enrichments
  std::set<XFEM::PHYSICS::Field> fieldset;
  fieldset.insert(XFEM::PHYSICS::Velx);
  fieldset.insert(XFEM::PHYSICS::Vely);
  fieldset.insert(XFEM::PHYSICS::Velz);
  fieldset.insert(XFEM::PHYSICS::Pres);
  const XFLUID::FluidElementAnsatz elementAnsatz;
  Teuchos::RCP<XFEM::DofManager> dofmanager =
    rcp(new XFEM::DofManager(ih, fieldset, elementAnsatz, xparams_));

  // tell elements about the dofs and the integration
  TransferDofInformationToElements(ih, dofmanager);
#else
  {
    ParameterList eleparams;
    eleparams.set("action","set_output_mode");
    eleparams.set("output_mode",true);
    discret_->Evaluate(eleparams);
  }
#endif

  discret_->FillComplete();


  output_->WriteMesh(0,0.0);


  // sanity check
  if (  solver_->Params().get<string>("solver") == "aztec"
    and solver_->Params().isSublist("ML Parameters")
    and not xparams_.get<bool>("DLM_condensation")
    )
  {
    dserror("for MLFLUID2, you need to set \"DLM_CONDENSATION  yes\" ");
  }

  if (actdis->Comm().MyPID()==0)
  {
    if (not xparams_.get<bool>("DLM_condensation"))
    {
      std::cout << RED_LIGHT << "DLM_condensation turned off!" << END_COLOR << endl;
    }
  }

  // store a dofset with the complete fluid unknowns
  dofset_out_.Reset();
  dofset_out_.AssignDegreesOfFreedom(*discret_,0);
  // split based on complete fluid field
  FLD::UTILS::SetupFluidSplit(*discret_,dofset_out_,3,velpressplitterForOutput_);

  state_.nodalDofDistributionMap_.clear();
  state_.elementalDofDistributionMap_.clear();

  {
    ParameterList eleparams;
    eleparams.set("action","set_output_mode");
    eleparams.set("output_mode",false);
    discret_->Evaluate(eleparams);
  }

  // get constant density variable for incompressible flow
  {
    ParameterList eleparams;
    eleparams.set("action","get_density");
    discret_->Evaluate(eleparams);
    density_ = eleparams.get<double>("density");
    if (density_ <= 0.0) dserror("received negative or zero density value from elements");
  }

  // print information about elements
  // (to double-check and log that correct input has been read)
  std::set<DRT::Element::DiscretizationType> distypeset;
  std::set<DRT::Element::ElementType> etypeset;
  for (int i=0; i<discret_->NumMyColElements(); ++i)
  {
    distypeset.insert(discret_->lColElement(i)->Shape());
    etypeset.insert(discret_->lColElement(i)->Type());
  }

  if (etypeset.count(DRT::Element::element_xdiff3) > 0)
  {
    diffusion_problem_ = true;
  }
  else
  {
    diffusion_problem_ = false;
  }


  physprob_.fieldset_.clear();
  if (diffusion_problem_)
  {
    cout0_ << endl<<endl << "Computing a diffusion problem!" << endl << endl;
    physprob_.fieldset_.insert(XFEM::PHYSICS::Temp);
    physprob_.elementAnsatzp_ = rcp(new XDIFF::Diff3ElementAnsatz());
  }
  else
  {
    physprob_.fieldset_.insert(XFEM::PHYSICS::Velx);
    physprob_.fieldset_.insert(XFEM::PHYSICS::Vely);
    physprob_.fieldset_.insert(XFEM::PHYSICS::Velz);
    physprob_.fieldset_.insert(XFEM::PHYSICS::Pres);
    physprob_.elementAnsatzp_ = rcp(new XFLUID::FluidElementAnsatz());
  }

  {
    cout0_ << "Element shapes in xfluid discretization: ";
    discret_->Comm().Barrier();
    bool moreThanOne = false;
    for (std::set<DRT::Element::DiscretizationType>::const_iterator iter = distypeset.begin(); iter != distypeset.end(); ++iter)
    {
      if (moreThanOne)  cout << ", ";
      cout << DRT::DistypeToString(*iter);
      moreThanOne = true;
    }
    if (actdis->Comm().MyPID()==0)
    {
      cout << endl;
    }
    discret_->Comm().Barrier();
  }

  {
    cout0_ << "Element types in xfluid discretization: ";
    discret_->Comm().Barrier();
    bool moreThanOnee = false;
    for (std::set<DRT::Element::ElementType>::const_iterator iter = etypeset.begin(); iter != etypeset.end(); ++iter)
    {
      if (moreThanOnee)  cout << ", ";
      cout << (*iter);
      moreThanOnee = true;
    }
    if (actdis->Comm().MyPID()==0)
    {
      cout << endl;
    }
    discret_->Comm().Barrier();
  }

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
void FLD::XFluidImplicitTimeInt::Integrate(
    const Teuchos::RCP<DRT::Discretization> cutterdiscret
        )
{
  // bound for the number of startsteps
  //const int    numstasteps         =params_.get<int>   ("number of start steps");

  // output of stabilization details
  if (myrank_==0)
  {
    const ParameterList *  stabparams=&(params_.sublist("STABILIZATION"));

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
  if (timealgo_==timeint_stationary) SolveStationaryProblem(cutterdiscret);
  else TimeLoop(cutterdiscret);

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
void FLD::XFluidImplicitTimeInt::TimeLoop(
    const Teuchos::RCP<DRT::Discretization> cutterdiscret
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
    /* additionally it remains to mention that for the linearised
       fluid the stbilisation is hard coded to be SUPG/PSPG */
  }

  const Epetra_Map* fluidsurface_dofcolmap = cutterdiscret->DofColMap();
  Teuchos::RCP<Epetra_Vector> idispcolnp  = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
  Teuchos::RCP<Epetra_Vector> ivelcolnp   = LINALG::CreateVector(*fluidsurface_dofcolmap,true);

  Teuchos::RCP<Epetra_Vector> idispcoln   = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
  Teuchos::RCP<Epetra_Vector> ivelcoln    = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
  Teuchos::RCP<Epetra_Vector> ivelcolnm   = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
  Teuchos::RCP<Epetra_Vector> iacccoln    = LINALG::CreateVector(*fluidsurface_dofcolmap,true);

  cutterdiscret->SetState("idispcolnp",idispcolnp);
  cutterdiscret->SetState("ivelcolnp",ivelcolnp);

  cutterdiscret->SetState("idispcoln",idispcoln);
  cutterdiscret->SetState("ivelcoln",ivelcoln);
  cutterdiscret->SetState("ivelcolnm",ivelcolnm);
  cutterdiscret->SetState("iacccoln",iacccoln);

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
      NonlinearSolve(cutterdiscret);
      break;
    case 1:
      // -----------------------------------------------------------------
      //                     solve linearised equation
      // -----------------------------------------------------------------
      //LinearSolve(cutterdiscret);
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
void FLD::XFluidImplicitTimeInt::PrepareTimeStep()
{

  // update interface handle
  ih_n_ = ih_np_;

  // update old acceleration
  if (state_.accn_ != Teuchos::null)
    state_.accn_->Update(1.0,*state_.accnp_,0.0);

  // velocities/pressures of this step become most recent
  // velocities/pressures of the last step
  if (state_.velnm_ != Teuchos::null)
    state_.velnm_->Update(1.0,*state_.veln_ ,0.0);
  if (state_.veln_ != Teuchos::null)
    state_.veln_ ->Update(1.0,*state_.velnp_,0.0);

  // -------------------------------------------------------------------
  //              set time dependent parameters
  // -------------------------------------------------------------------
  step_ += 1;
  time_ += dta_;

  if (params_.get<FLUID_TIMEINTTYPE>("time int algo") == timeint_stationary)
  {
    timealgo_ = timeint_stationary;
    theta_ = 1.0;
  }
  else
  {
    // do a backward Euler step for the first timestep
    if (step_==1)
    {
      timealgo_ = timeint_one_step_theta;
      theta_ = params_.get<double>("start theta");
    }
    else
    {
      timealgo_ = params_.get<FLUID_TIMEINTTYPE>("time int algo");
      theta_ = params_.get<double>("theta");

      // for BDF2, theta is set by the time-step sizes, 2/3 for const. dt
      if (timealgo_==timeint_bdf2) theta_ = (dta_+dtp_)/(2.0*dta_ + dtp_);
    }
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
          state_.veln_, state_.velnm_, state_.accn_,
          timealgo_, dta_, dtp_,
          state_.velnp_);
    }
  }
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::XFluidImplicitTimeInt::PrepareNonlinearSolve()
{

  // -------------------------------------------------------------------
  // set part of the rhs vector belonging to the old timestep
  // -------------------------------------------------------------------
//  TIMEINT_THETA_BDF2::SetOldPartOfRighthandside(
//      state_.veln_, state_.velnm_, state_.accn_,
//          timealgo_, dta_, theta_,
//          hist_);

  // -------------------------------------------------------------------
  //         evaluate Dirichlet and Neumann boundary conditions
  // -------------------------------------------------------------------
  {
    ParameterList eleparams;

    // total time required for Dirichlet conditions
    eleparams.set("total time",time_);
    eleparams.set("delta time",dta_);
    eleparams.set("thsl",theta_*dta_);

    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("velnp",state_.velnp_);
    // predicted dirichlet values
    // velnp then also holds prescribed new dirichlet values
    discret_->EvaluateDirichlet(eleparams,state_.velnp_,null,null,null,dbcmaps_);
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

    neumann_loads_->PutScalar(0.0);

    // evaluate Neumann conditions
    discret_->EvaluateNeumann(eleparams,*neumann_loads_);
    discret_->ClearState();
  }

}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::XFluidImplicitTimeInt::TransferDofInformationToElements(
    const Teuchos::RCP<XFEM::InterfaceHandleXFSI>  ih,
    const Teuchos::RCP<XFEM::DofManager> dofmanager
    )
{
  ParameterList eleparams;
  eleparams.set("action","store_xfem_info");
  eleparams.set("dofmanager",dofmanager);
  eleparams.set("DLM_condensation",xparams_.get<bool>("DLM_condensation"));
  eleparams.set("boundaryRatioLimit",xparams_.get<double>("boundaryRatioLimit"));
  eleparams.set("interfacehandle",ih);
  discret_->Evaluate(eleparams);
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
Teuchos::RCP<XFEM::InterfaceHandleXFSI> FLD::XFluidImplicitTimeInt::ComputeInterfaceAndSetDOFs(
    const Teuchos::RCP<DRT::Discretization>  cutterdiscret
    )
{
  // dump old matrix to save memory while we construct a new matrix
  sysmat_ = Teuchos::null;

//  {
//    ParameterList eleparams;
//    eleparams.set("action","reset");
//    discret_->Evaluate(eleparams);
//  }
  dofmanagerForOutput_ = Teuchos::null;

  // within this routine, no parallel re-distribution is allowed to take place
  // before and after this function, it's ok to do that

  // compute Intersection
  ih_np_ = rcp(new XFEM::InterfaceHandleXFSI(discret_, cutterdiscret));
  ih_np_->toGmsh(step_);

  // apply enrichments
  const Teuchos::RCP<XFEM::DofManager> dofmanager =
      rcp(new XFEM::DofManager(ih_np_, physprob_.fieldset_, *physprob_.elementAnsatzp_, xparams_));

  // save dofmanager to be able to plot Gmsh stuff in Output()
  dofmanagerForOutput_ = dofmanager;

  // tell elements about the dofs and the integration
  TransferDofInformationToElements(ih_np_, dofmanager);

  // print global and element dofmanager to Gmsh
  dofmanager->toGmsh(step_);


  // get old dofmap, compute new one and get the new one, too
  const Epetra_Map olddofrowmap = *discret_->DofRowMap();
  discret_->FillComplete(true,false,false);
  const Epetra_Map& newdofrowmap = *discret_->DofRowMap();

  {
    const std::map<XFEM::DofKey<XFEM::onNode>, XFEM::DofGID> oldNodalDofDistributionMap(state_.nodalDofDistributionMap_);
    const std::map<XFEM::DofKey<XFEM::onElem>, XFEM::DofGID> oldElementalDofDistributionMap(state_.elementalDofDistributionMap_);
    dofmanager->fillDofDistributionMaps(
        state_.nodalDofDistributionMap_,
        state_.elementalDofDistributionMap_);

    // create switcher
    const XFEM::DofDistributionSwitcher dofswitch(
            ih_np_, dofmanager,
            olddofrowmap, newdofrowmap,
            oldNodalDofDistributionMap, state_.nodalDofDistributionMap_,
            oldElementalDofDistributionMap, state_.elementalDofDistributionMap_
            );

    // --------------------------------------------
    // switch state vectors to new dof distribution
    // --------------------------------------------

    cout0_ << " ->  Initialize system vectors..." << endl;
    // accelerations
    dofswitch.mapVectorToNewDofDistribution(state_.accnp_);
    dofswitch.mapVectorToNewDofDistribution(state_.accn_);

    // velocities and pressures
    dofswitch.mapVectorToNewDofDistribution(state_.velnp_);
    dofswitch.mapVectorToNewDofDistribution(state_.veln_);
    dofswitch.mapVectorToNewDofDistribution(state_.velnm_);
  }

//  // fix ghost fluid values
//  const int standard_label = 0;
//  const XFEM::Enrichment enr_std(standard_label, XFEM::Enrichment::typeStandard);
//  for (int i=0; i<discret_->NumMyColNodes(); ++i)
//  {
//    const DRT::Node* node = discret_->lColNode(i);
//    const LINALG::Matrix<3,1> nodalpos(node->X());
//    const bool is_in_fluid = (ih->PositionWithinConditionNP(nodalpos) == 0);
//    const bool was_in_fluid = (ih->PositionWithinConditionN(nodalpos) == 0);
//
//    const bool voidenrichment_in_set = EnrichmentInNodalDofSet(node->Id(), XFEM::Enrichment::typeVoid, nodalDofSet);
//
//    if (not voidenrichment_in_set)
//    {
//      const int label = ih.PositionWithinConditionNP(nodalpos);
//      const bool in_fluid = (0 == label);
//
//      if (in_fluid)
//      {
//        for (std::set<XFEM::PHYSICS::Field>::const_iterator field = fieldset.begin();field != fieldset.end();++field)
//        {
//          nodalDofSet[node->Id()].insert(XFEM::FieldEnr(*field, enr_std));
//        }
//      }
//    }
//  };

  // --------------------------------------------------
  // create remaining vectors with new dof distribution
  // --------------------------------------------------
  //hist_         = LINALG::CreateVector(newdofrowmap,true);

  zeros_        = LINALG::CreateVector(newdofrowmap,true);

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
  FLD::UTILS::SetupXFluidSplit(*discret_,dofmanager,velpressplitter_);

  // -------------------------------------------------------------------
  // create empty system matrix --- stiffness and mass are assembled in
  // one system matrix!
  // -------------------------------------------------------------------
  cout0_ << " ->  Initialize system matrix..." << endl;

  // initialize system matrix
  sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(newdofrowmap,0,false,true));

#if 0

  // rotation of matrix to improve matrix condition number

  // find nodes connected to enrichments
  set<int> enr_node_gids;
  for (int iele=0; iele < discret_->NumMyRowElements(); iele++)
  {
    const DRT::Element* ele = discret_->lRowElement(iele);
    if (ih->ElementIntersected(ele->Id()))
    {
      for (int inode=0;inode<ele->NumNode();inode++)
      {
        enr_node_gids.insert(ele->NodeIds()[inode]);
      }
    }
  }

  set<int> ext_enr_ele_gids;
  for (set<int>::const_iterator nodeid = enr_node_gids.begin(); nodeid != enr_node_gids.end();nodeid++)
  {
    const DRT::Node* node = discret_->gNode(*nodeid);
    for (int iele=0;iele<node->NumElement();iele++)
    {
      ext_enr_ele_gids.insert(node->Elements()[iele]->Id());
    }
  }

  set<int> ext_enr_node_gids;
  for (set<int>::const_iterator eleid = ext_enr_ele_gids.begin(); eleid != ext_enr_ele_gids.end();eleid++)
  {
    const DRT::Element* ele = discret_->lRowElement(*eleid);
    for (int inode=0;inode<ele->NumNode();inode++)
    {
      ext_enr_node_gids.insert(ele->NodeIds()[inode]);
    }
  }

  FLD::UTILS::SetupEnrichmentSplit(*discret_,dofmanager,ext_enr_node_gids, normalenrichedsplitter_);

  const Teuchos::RCP<const Epetra_Map>& enrichmap = normalenrichedsplitter_.CondMap();
  const Teuchos::RCP<const Epetra_Map>& normalmap = normalenrichedsplitter_.OtherMap();


//  velpressplitter_.ExtractOtherVector(incvel_,onlyvel);

  // create small dofmap
//  Epetra_Map newdofrowmap;


  // create smaller matrix A
  enrichsysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*enrichmap,108,false,true));
#endif

  // print information about dofs
  if (discret_->Comm().NumProc() == 1)
  {
    const int numdof = newdofrowmap.NumGlobalElements();
    const int numnodaldof = dofmanager->NumNodalDof();
    cout << " DOF report: numdof = " << numdof << ", numstressdof = "<< (numdof - numnodaldof) << endl;
  }
  else
  {
    if(discret_->Comm().MyPID() == 0)
    {
      cout << " DOF report: numdof = " << newdofrowmap.NumGlobalElements() << endl;
    }
  }

//  cout0_ << *state_.accnp_ << endl;
//  cout0_ << *state_.accn_ << endl;
//
//  cout0_ << *state_.velnp_ << endl;
//  cout0_ << *state_.veln_ << endl;
//  cout0_ << *state_.velnm_ << endl;

  cout0_ << "Setup phase done!" << endl;

  return ih_np_;
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | contains Krylov Projection                                gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::XFluidImplicitTimeInt::KrylovSpaceProjection(const int & numdim)
{
  if (project_)
  {
    DRT::Condition* KSPcond=discret_->GetCondition("KrylovSpaceProjection");

    // in this case, we want to project out some zero pressure modes
    const string* definition = KSPcond->Get<string>("weight vector definition");

    // zero w and c
    w_->PutScalar(0.0);
    c_->PutScalar(0.0);

    if(*definition == "pointvalues")
    {
      // get pressure
      const vector<double>* mode = KSPcond->Get<vector<double> >("mode");

      for(int rr=0;rr<numdim;++rr)
      {
        if(abs((*mode)[rr]>1e-14))
        {
          dserror("expecting only an undetermined pressure");
        }
      }

      int predof = numdim;

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
      LINALG::Export(*presmode,*w_);

      // export to vector of ones
      presmode->PutScalar(1.0);
      LINALG::Export(*presmode,*c_);
    }
    else if(*definition == "integration")
    {
      ParameterList mode_params;

      // set action for elements
      mode_params.set("action","integrate_shape");

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

      for(int rr=0;rr<numdim;++rr)
      {
        if(abs((*mode)[rr]>1e-14))
        {
          dserror("expecting only an undetermined pressure");
        }
      }

      Teuchos::RCP<Epetra_Vector> presmode = velpressplitter_.ExtractCondVector(*w_);

      // export to vector of ones
      presmode->PutScalar(1.0);
      LINALG::Export(*presmode,*c_);
    }
    else
    {
      dserror("unknown definition of weight vector w for restriction of Krylov space");
    }
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
void FLD::XFluidImplicitTimeInt::NonlinearSolve(
    const Teuchos::RCP<DRT::Discretization> cutterdiscret
    )
{
  cout << "FLD::XFluidImplicitTimeInt::NonlinearSolve()" << endl;
  // time measurement: nonlinear iteration
  TEUCHOS_FUNC_TIME_MONITOR("   + nonlin. iteration/lin. solve");

  {
    ParameterList eleparams;
    eleparams.set("action","reset");
    discret_->Evaluate(eleparams);
  }

  // since number of dofs can change, we don't allow reuse of preconditioners over timestep borders
  solver_ = Teuchos::null;

  ComputeInterfaceAndSetDOFs(cutterdiscret);

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
    cout0_ << "Krylov projection active..." << endl;
    project_ = true;
    w_       = LINALG::CreateVector(*discret_->DofRowMap(),true);
    c_       = LINALG::CreateVector(*discret_->DofRowMap(),true);
  }
  else if (numfluid == 0)
  {
    project_ = false;
    w_       = Teuchos::null;
    c_       = Teuchos::null;
  }
  else
    dserror("Received more than one KrylovSpaceCondition for fluid field");



  solver_ = rcp(new LINALG::Solver(fluidsolverparams_,
                                   discret_->Comm(),
                                   DRT::Problem::Instance()->ErrorFile()->Handle()));
  discret_->ComputeNullSpaceIfNecessary(solver_->Params());

//  DRT::PAR::LoadBalancer balancer(discret_);
//  balancer.Partition();
//
//  ComputeInterfaceAndSetDOFs(cutterdiscret);
//
  PrepareNonlinearSolve();

  // ---------------------------------------------- nonlinear iteration
  // ------------------------------- stop nonlinear iteration when both
  //                                 increment-norms are below this bound
  const double  ittol     =params_.get<double>("tolerance for nonlin iter");

  //------------------------------ turn adaptive solver tolerance on/off
  const bool   isadapttol    = params_.get<bool>("ADAPTCONV");
  const double adaptolbetter = params_.get<double>("ADAPTCONV_BETTER");

  //const bool fluidrobin = params_.get<bool>("fluidrobin", false);

  int               itnum = 0;
  const int         itemax = params_.get<int>("max nonlin iter steps");
  bool              stopnonliniter = false;


  double dtsolve = 0.0;
  double dtele   = 0.0;

  // action for elements
  if (timealgo_!=timeint_stationary and theta_ < 1.0)
  {
    cout0_ << "* Warning! Works reliable only for Backward Euler time discretization! *" << endl;
  }

#if 0
  LINALG::Matrix<3,1> x_test;
  x_test(0) = 0.5;
  x_test(1) = 0.5;
  x_test(2) = 0.01;

  double min_dist = 1.0e12;
  int min_node = -1;

  // brute force node finding for test cases
  for (int i=0; i<discret_->NumMyColNodes(); ++i)
  {
    const DRT::Node* node = discret_->lColNode(i);
    const LINALG::Matrix<3,1> x_node(node->X());

    LINALG::Matrix<3,1> distance_vector;
    // vector pointing away from the node towards physCoord
    distance_vector.Update(1.0, x_test, -1.0, x_node);

    // absolute distance between point and node
    const double distance = distance_vector.Norm2();

    if (distance < min_dist)
    {
      min_dist = distance;
      min_node = node->Id();
    }
  }
#endif

  if (myrank_ == 0)
  {
    printf("+------------+------------------+---------+---------+---------+---------+---------+---------+\n");
    printf("|  step/max  |  tol     [norm]  | vel-res | pre-res | fullres | vel-inc | pre-inc | fullinc |\n");
  }

  const Teuchos::RCP<Epetra_Vector> iforcecolnp = LINALG::CreateVector(*cutterdiscret->DofColMap(),true);

  incvel_->PutScalar(0.0);
  residual_->PutScalar(0.0);

  // increment of the old iteration step - used for update of condensed element stresses
  oldinc_ = LINALG::CreateVector(*discret_->DofRowMap(),true);


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
      discret_->Comm().Barrier();
      const double tcpu=ds_cputime();

      // add Neumann loads
      residual_->Update(1.0,*neumann_loads_,0.0);

      // create the parameters for the discretization
      ParameterList eleparams;

      // action for elements
      if (timealgo_==timeint_stationary)
        eleparams.set("action","calc_fluid_stationary_systemmat_and_residual");
      else
        eleparams.set("action","calc_fluid_systemmat_and_residual");

      // set general element parameters
      //eleparams.set("total time",time_);
      //eleparams.set("thsl",theta_*dta_);
      eleparams.set("timealgo",timealgo_);
      eleparams.set("dt",dta_);
      eleparams.set("theta",theta_);
      eleparams.set("include reactive terms for linearisation",params_.get<bool>("Use reaction terms for linearisation"));

      // parameters for stabilization
      eleparams.sublist("STABILIZATION") = params_.sublist("STABILIZATION");

      // parameters for stabilization
      eleparams.sublist("TURBULENCE MODEL") = params_.sublist("TURBULENCE MODEL");

      // set general vector values needed by elements
      discret_->ClearState();
      discret_->SetState("velnp",state_.velnp_);
      discret_->SetState("veln" ,state_.veln_);
      discret_->SetState("velnm",state_.velnm_);
      discret_->SetState("accn" ,state_.accn_);

      discret_->SetState("nodal increment",oldinc_);

      eleparams.set("interface velocity",cutterdiscret->GetState("ivelcolnp"));

      // reset interface force and let the elements fill it
      iforcecolnp->PutScalar(0.0);
      eleparams.set("interface force",iforcecolnp);

      double L2 = 0.0;
      eleparams.set("L2",L2);

      eleparams.set("DLM_condensation",xparams_.get<bool>("DLM_condensation"));
      eleparams.set("FAST_INTEGRATION",xparams_.get<bool>("FAST_INTEGRATION"));
      eleparams.set("boundaryRatioLimit",xparams_.get<double>("boundaryRatioLimit"));
      eleparams.set("monolithic_FSI",false);

      // convergence check at itemax is skipped for speedup if
      // CONVCHECK is set to L_2_norm_without_residual_at_itemax
      if ((itnum != itemax)
          ||
          (params_.get<string>("CONVCHECK")
           !=
           "L_2_norm_without_residual_at_itemax"))
      {
        sysmat_->Zero();
        // call standard loop over elements
        discret_->Evaluate(eleparams,sysmat_,residual_);
        discret_->ClearState();

        // account for potential Neumann inflow terms
        if (params_.get<string>("Neumann inflow") == "yes")
        {
          // create parameter list
          ParameterList condparams;

          // action for elements
          condparams.set("action","calc_Neumann_inflow");
          condparams.set("thsl",theta_*dta_);

          // set vector values needed by elements
          discret_->ClearState();
          condparams.set("using generalized-alpha time integration",false);
          discret_->SetState("velnp",state_.velnp_);

          std::string condstring("FluidNeumannInflow");
          discret_->EvaluateCondition(condparams,sysmat_,Teuchos::null,residual_,Teuchos::null,Teuchos::null,condstring);
          discret_->ClearState();
        }

        // finalize the complete matrix
        sysmat_->Complete();
      }

      {
        const double L2_result = sqrt(eleparams.get<double>("L2"));
        //cout << "L2 norm = " << scientific << L2_result << endl;

        const std::string fname("L2.txt");
        std::ofstream f(fname.c_str(),std::fstream::trunc);
        f.setf(ios::scientific,ios::floatfield);
        f.precision(12);
        f << L2_result;
        f.close();
      }

#if 0
      std::vector<double> myvel(1);
      std::vector<int> gdofs(1);
      discret_->Dof(discret_->gNode(min_node),0,gdofs);
      DRT::UTILS::ExtractMyValues(*state_.velnp_,myvel,gdofs);
      const double temp = myvel[0];
      cout << "Temp("<< x_test(0) <<", "<< x_test(1) <<", "<< x_test(02) <<") = " << scientific << temp << endl;

      {
        std::ofstream f;
        const std::string fname = "Temp.txt";
        f.open(fname.c_str(),std::fstream::trunc);
        f.setf(ios::scientific,ios::floatfield);
        f.precision(12);
        f << temp;
        f.close();
      }
#endif

      // end time measurement for element
      // to get realistic element times, this barrier results in the longest
      // time being printed (processors with no intersected elements are too fast...)
      discret_->Comm().Barrier();
      dtele=ds_cputime()-tcpu;
    }

    // blank residual DOFs which are on Dirichlet BC
    dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), residual_);

    double incvelnorm_L2;
    double velnorm_L2;
    double vresnorm;

    Teuchos::RCP<Epetra_Vector> onlyvel = velpressplitter_.ExtractOtherVector(residual_);
    onlyvel->Norm2(&vresnorm);

    velpressplitter_.ExtractOtherVector(incvel_,onlyvel);
    onlyvel->Norm2(&incvelnorm_L2);

    velpressplitter_.ExtractOtherVector(state_.velnp_,onlyvel);
    onlyvel->Norm2(&velnorm_L2);

    double incprenorm_L2;
    double prenorm_L2;
    double presnorm;

    Teuchos::RCP<Epetra_Vector> onlypre = velpressplitter_.ExtractCondVector(residual_);
    onlypre->Norm2(&presnorm);

    velpressplitter_.ExtractCondVector(incvel_,onlypre);
    onlypre->Norm2(&incprenorm_L2);

    velpressplitter_.ExtractCondVector(state_.velnp_,onlypre);
    onlypre->Norm2(&prenorm_L2);

    double incfullnorm_L2;
    double fullnorm_L2;
    double fullresnorm;

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
        printf(" (      --    , te=%9.2E)",dtele);
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
          printf("|  %3d/%3d   | %9.2e[L_2 ]  | %7.1e | %7.1e | %7.1e | %7.1e | %7.1e | %7.1e |",
                 itnum,itemax,ittol,vresnorm,presnorm,fullresnorm,
                 incvelnorm_L2/velnorm_L2,incprenorm_L2/prenorm_L2,incfullnorm_L2/fullnorm_L2);
          printf(" (ts=%9.2e, te=%9.2e)\n",dtsolve,dtele);
          printf("+------------+------------------+---------+---------+---------+---------+---------+---------+\n");

          FILE* errfile = params_.get<FILE*>("err file");
          if (errfile!=NULL)
          {
            fprintf(errfile,"fluid solve:   %3d/%3d  tol=%9.2E[L_2 ]  vres=%9.2E  pres=%9.2E  vinc=%9.2E  pinc=%9.2E\n",
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
          printf(" (ts=%9.2e, te=%9.2e)",dtsolve,dtele);
        }
    }

    // warn if itemax is reached without convergence, but proceed to
    // next timestep...
    if ((itnum == itemax) and (vresnorm > ittol or presnorm > ittol or
                             fullresnorm > ittol or
                             incvelnorm_L2/velnorm_L2 > ittol or
                             incprenorm_L2/prenorm_L2 > ittol or
                             incfullnorm_L2/fullnorm_L2 > ittol))
    {
      stopnonliniter=true;
      if (myrank_ == 0)
      {
        printf("\n");
        printf("+---------------------------------------------------------------+\n");
        printf("|            >>>>>> not converged in itemax steps!              |\n");
        printf("+---------------------------------------------------------------+\n");

        FILE* errfile = params_.get<FILE*>("err file");
        if (errfile!=NULL)
        {
          fprintf(errfile,"fluid unconverged solve:   %3d/%3d  tol=%9.2E[L_2 ]  vres=%9.2E  pres=%9.2E  vinc=%9.2E  pinc=%9.2E\n",
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

    // estimate condition number of the tangent stiffness matrix
    if (xparams_.get<bool>("CONDEST"))
    {
      const double cond_number = LINALG::Condest(dynamic_cast<LINALG::SparseMatrix&>(*sysmat_),Ifpack_GMRES, 100);
      // computation of significant digits might be completely bogus, so don't take it serious
      const double tmp = std::abs(std::log10(cond_number*1.11022e-16));
      const int sign_digits = (int)std::floor(tmp);
      cout0_ << " cond est: " << scientific << cond_number << ", max.sign.digits: " << sign_digits;
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
      discret_->Comm().Barrier();
      const double tcpusolve=ds_cputime();

      // do adaptive linear solver tolerance (not in first solve)
      if (isadapttol and itnum>1)
      {
        double currresidual = max(vresnorm,presnorm);
        currresidual = max(currresidual,incvelnorm_L2/velnorm_L2);
        currresidual = max(currresidual,incprenorm_L2/prenorm_L2);
        solver_->AdaptTolerance(ittol,currresidual,adaptolbetter);
      }

      const int numdim = 3;
      KrylovSpaceProjection(numdim);
      solver_->Solve(sysmat_->EpetraOperator(), incvel_, residual_, true, itnum == 1, w_, c_, project_);
      solver_->ResetTolerance();

      // end time measurement for solver
      // to get realistic times, this barrier results in the longest
      // time being printed
      discret_->Comm().Barrier();
      dtsolve = ds_cputime()-tcpusolve;
    }

    // -------------------------------------------------------------------
    // update velocity and pressure values by increments
    // -------------------------------------------------------------------
    state_.velnp_->Update(1.0,*incvel_,1.0);

    //------------------- store nodal increment for element stress update
    oldinc_->Update(1.0,*incvel_,0.0);
  }

  // macht der FSI algorithmus
  iforcecolnp->Scale(-ResidualScaling());

  cutterdiscret->SetState("iforcenp", iforcecolnp);

  oldinc_ = Teuchos::null;
  incvel_ = Teuchos::null;
  sysmat_ = Teuchos::null;

} // FluidImplicitTimeInt::NonlinearSolve



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | build linear system matrix and rhs                        u.kue 11/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::XFluidImplicitTimeInt::Evaluate(
    const Teuchos::RCP<DRT::Discretization> cutterdiscret,
    Teuchos::RCP<const Epetra_Vector> fluidvel)
{
  cout << "FLD::XFluidImplicitTimeInt::Evaluate()" << endl;

  const bool firstcall = state_.velnp_ == Teuchos::null;

  if (state_.velnp_ != Teuchos::null)
    state_.velnp_->Update(1.0,*fluidvel,1.0);

  if (firstcall)
  {
    cout << "Resetting!!!!!!!!!!!!!!!!!!!" << endl;
    ParameterList eleparams;
    eleparams.set("action","reset");
    discret_->Evaluate(eleparams);
  }

  Teuchos::RCP<XFEM::InterfaceHandleXFSI> ih = ComputeInterfaceAndSetDOFs(cutterdiscret);

  PrepareNonlinearSolve();

  const Epetra_Map& fluiddofrowmap = *discret_->DofRowMap();
  const Epetra_Map& ifacedofrowmap = *cutterdiscret->DofRowMap();

  Cuu_ = Teuchos::rcp(new LINALG::SparseMatrix(fluiddofrowmap,0,false,false));
  Mud_ = Teuchos::rcp(new LINALG::SparseMatrix(fluiddofrowmap,0,false,false));
  Mdu_ = Teuchos::rcp(new LINALG::SparseMatrix(ifacedofrowmap,0,false,false));
  Cdd_ = Teuchos::rcp(new LINALG::SparseMatrix(ifacedofrowmap,0,false,false));
  RHSd_ = LINALG::CreateVector(ifacedofrowmap,true);

  // action for elements
  if (timealgo_!=timeint_stationary and theta_ < 1.0)
  {
    cout0_ << "* Warning! Works reliable only for Backward Euler time discretization! *" << endl;
  }

  const Teuchos::RCP<Epetra_Vector> iforcecolnp = LINALG::CreateVector(*cutterdiscret->DofColMap(),true);

  incvel_->PutScalar(0.0);

  // increment of the old iteration step - used for update of condensed element stresses
  if (oldinc_ == Teuchos::null)
    oldinc_ = LINALG::CreateVector(*discret_->DofRowMap(),true);
  else
  {
    if (oldinc_->MyLength() != incvel_->MyLength())
      oldinc_ = LINALG::CreateVector(*discret_->DofRowMap(),true);
  }


  // add Neumann loads
  residual_->Update(1.0,*neumann_loads_,0.0);

  // create the parameters for the discretization
  ParameterList eleparams;

  // action for elements
  if (timealgo_==timeint_stationary)
    eleparams.set("action","calc_fluid_stationary_systemmat_and_residual");
  else
    eleparams.set("action","calc_fluid_systemmat_and_residual");

  // set general element parameters
  //eleparams.set("total time",time_);
  //eleparams.set("thsl",theta_*dta_);
  eleparams.set("timealgo",timealgo_);
  eleparams.set("dt",dta_);
  eleparams.set("theta",theta_);
  eleparams.set("include reactive terms for linearisation",params_.get<bool>("Use reaction terms for linearisation"));

  // parameters for stabilization
  eleparams.sublist("STABILIZATION") = params_.sublist("STABILIZATION");

  // parameters for stabilization
  eleparams.sublist("TURBULENCE MODEL") = params_.sublist("TURBULENCE MODEL");

  // set general vector values needed by elements
  discret_->ClearState();
  discret_->SetState("velnp",state_.velnp_);
  discret_->SetState("veln" ,state_.veln_);
  discret_->SetState("velnm",state_.velnm_);
  discret_->SetState("accn" ,state_.accn_);

  discret_->SetState("nodal increment",oldinc_);

  // reset interface force and let the elements fill it
  eleparams.set("interface force",iforcecolnp);

  double L2 = 0.0;
  eleparams.set("L2",L2);

  eleparams.set("DLM_condensation",xparams_.get<bool>("DLM_condensation"));
  eleparams.set("FAST_INTEGRATION",xparams_.get<bool>("FAST_INTEGRATION"));
  eleparams.set("boundaryRatioLimit",xparams_.get<double>("boundaryRatioLimit"));
  eleparams.set<bool>("monolithic_FSI",true);

  // call loop over elements
  //discret_->Evaluate(eleparams,sysmat_,residual_);
  MonolithicMultiDisEvaluate(ih, discret_,cutterdiscret, eleparams,
      Cuu_, Mud_, Mdu_, Cdd_, RHSd_, sysmat_, residual_);
  discret_->ClearState();

  // finalize the system matrix
  sysmat_->Complete();
  cout << "doing Cuu_" << endl;
  Cuu_->Complete(fluiddofrowmap, fluiddofrowmap);
  cout << "doing Mud_" << endl;
  Mud_->Complete(ifacedofrowmap, fluiddofrowmap);
  cout << "doing Mdu_" << endl;
  Mdu_->Complete(fluiddofrowmap, ifacedofrowmap);
  cout << "doing Cdd_" << endl;
  Cdd_->Complete(ifacedofrowmap, ifacedofrowmap);

  // blank residual DOFs which are on Dirichlet BC
  dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), residual_);
  trueresidual_->Update(ResidualScaling(),*residual_,0.0);

  // Apply dirichlet boundary conditions to system of equations
  // residual displacements are supposed to be zero at boundary
  // conditions
  incvel_->PutScalar(0.0);
  LINALG::ApplyDirichlettoSystem(sysmat_,incvel_,residual_,zeros_,*(dbcmaps_->CondMap()));

  //------------------- store nodal increment for element stress update
  oldinc_->Update(1.0,*incvel_,0.0);

  // macht der FSI algorithmus
  iforcecolnp->Scale(-ResidualScaling());

  cutterdiscret->SetState("iforcenp", iforcecolnp);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::map<std::string, Teuchos::RCP<LINALG::SparseMatrix> > FLD::XFluidImplicitTimeInt::CouplingMatrices()
{
  std::map<std::string, Teuchos::RCP<LINALG::SparseMatrix> > matrices;
  matrices.insert(make_pair("Cuu",Cuu_));
  matrices.insert(make_pair("Mud",Mud_));
  matrices.insert(make_pair("Mdu",Mdu_));
  matrices.insert(make_pair("Cdd",Cdd_));
  return matrices;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::map<std::string, Teuchos::RCP<Epetra_Vector> > FLD::XFluidImplicitTimeInt::CouplingVectors()
{
  std::map<std::string, Teuchos::RCP<Epetra_Vector> > vectors;
  vectors.insert(make_pair("rhsd",RHSd_));
  return vectors;
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
void FLD::XFluidImplicitTimeInt::TimeUpdate()
{

  // compute acceleration
  TIMEINT_THETA_BDF2::CalculateAcceleration(
      state_.velnp_, state_.veln_, state_.velnm_, state_.accn_,
          timealgo_, step_, theta_, dta_, dtp_,
          state_.accnp_);

  oldinc_= Teuchos::null;

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
void FLD::XFluidImplicitTimeInt::StatisticsAndOutput()
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
  ComputeSurfaceFlowRates();

  // -------------------------------------------------------------------
  //          calculate impuls rate through surfaces
  // -------------------------------------------------------------------
  ComputeSurfaceImpulsRates();

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
//  statisticsmanager_->DoOutput(output_,step_);

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
void FLD::XFluidImplicitTimeInt::Output()
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
    Teuchos::RCP<Epetra_Vector> velnp_out = dofmanagerForOutput_->fillPhysicalOutputVector(
        *state_.velnp_, dofset_out_, state_.nodalDofDistributionMap_, physprob_.fieldset_);

    // write physical fields on full domain including voids etc.
    if (physprob_.fieldset_.find(XFEM::PHYSICS::Velx) != physprob_.fieldset_.end())
    {
      // output velocity field for visualization
      output_->WriteVector("velocity", velnp_out);

      // output (hydrodynamic) pressure for visualization
      Teuchos::RCP<Epetra_Vector> pressure = velpressplitterForOutput_.ExtractCondVector(velnp_out);
      pressure->Scale(density_);
      output_->WriteVector("pressure", pressure);

      //output_->WriteVector("residual", trueresidual_);

      //only perform stress calculation when output is needed
      if (writestresses_)
      {
        Teuchos::RCP<Epetra_Vector> traction = CalcStresses();
        output_->WriteVector("traction",traction);
      }
    }
    else if (physprob_.fieldset_.find(XFEM::PHYSICS::Temp) != physprob_.fieldset_.end())
    {
      output_->WriteVector("temperature", velnp_out);
    }

    // write domain decomposition for visualization
    output_->WriteElementData();
  }

  // write restart
  if (write_restart_data)
  {
    output_->WriteVector("velnp", state_.velnp_);
    output_->WriteVector("veln" , state_.veln_);
    output_->WriteVector("velnm", state_.velnm_);
    output_->WriteVector("accnp", state_.accnp_);
    output_->WriteVector("accn" , state_.accn_);
  }

  if (discret_->Comm().NumProc() == 1)
  {
    OutputToGmsh(step_, time_);
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
void FLD::XFluidImplicitTimeInt::ReadRestart(
    int step,
    const Teuchos::RCP<DRT::Discretization> cutterdiscret
    )
{

  const int output_test_step = 999999;

  const Teuchos::RCP<XFEM::InterfaceHandleXFSI> ih =
    rcp(new XFEM::InterfaceHandleXFSI(discret_, cutterdiscret));
  ih->toGmsh(output_test_step);

  // apply enrichments
  std::set<XFEM::PHYSICS::Field> fieldset;
  fieldset.insert(XFEM::PHYSICS::Velx);
  fieldset.insert(XFEM::PHYSICS::Vely);
  fieldset.insert(XFEM::PHYSICS::Velz);
  fieldset.insert(XFEM::PHYSICS::Pres);
  const XFLUID::FluidElementAnsatz elementAnsatz;
  const Teuchos::RCP<XFEM::DofManager> dofmanager =
    rcp(new XFEM::DofManager(ih, fieldset, elementAnsatz, xparams_));

  // save dofmanager to be able to plot Gmsh stuff in Output()
  dofmanagerForOutput_ = dofmanager;

  // tell elements about the dofs and the integration
  TransferDofInformationToElements(ih, dofmanager);

  // print global and element dofmanager to Gmsh
  dofmanager->toGmsh(output_test_step);


  // get old dofmap, compute new one and get the new one, too
  discret_->FillComplete();

  dofmanager->fillDofDistributionMaps(
      state_.nodalDofDistributionMap_,
      state_.elementalDofDistributionMap_);



  IO::DiscretizationReader reader(discret_,step);
  time_ = reader.ReadDouble("time");
  step_ = reader.ReadInt("step");

  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // Vectors passed to the element
  // -----------------------------
  // velocity/pressure at time n+1, n and n-1
  state_.velnp_        = LINALG::CreateVector(*dofrowmap,true);
  state_.veln_         = LINALG::CreateVector(*dofrowmap,true);
  state_.velnm_        = LINALG::CreateVector(*dofrowmap,true);

  // acceleration at time n+1 and n
  state_.accnp_        = LINALG::CreateVector(*dofrowmap,true);
  state_.accn_         = LINALG::CreateVector(*dofrowmap,true);

  reader.ReadVector(state_.velnp_,"velnp");
  reader.ReadVector(state_.veln_, "veln");
  reader.ReadVector(state_.velnm_,"velnm");
  reader.ReadVector(state_.accnp_ ,"accnp");
  reader.ReadVector(state_.accn_ ,"accn");

  if (discret_->Comm().NumProc() == 1)
  {
    OutputToGmsh(output_test_step, time_);
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
void FLD::XFluidImplicitTimeInt::OutputToGmsh(
    const int step,
    const double time
    ) const
{
  const Teuchos::ParameterList& xfemparams = DRT::Problem::Instance()->XFEMGeneralParams();
  const bool gmshdebugout = getIntegralValue<int>(xfemparams,"GMSH_DEBUG_OUT")==1;

  const bool screen_out = true;

  if (gmshdebugout and (this->physprob_.fieldset_.find(XFEM::PHYSICS::Pres) != this->physprob_.fieldset_.end()))
  {
    cout << "XFluidImplicitTimeInt::OutputToGmsh()" << endl;

    std::stringstream filename;
    std::stringstream filenamedel;
    filename    << DRT::Problem::Instance()->OutputControlFile()->FileName() << ".solution_field_pressure_" << std::setw(5) << setfill('0') << step   << ".pos";
    filenamedel << DRT::Problem::Instance()->OutputControlFile()->FileName() << ".solution_field_pressure_" << std::setw(5) << setfill('0') << step-5 << ".pos";
    std::remove(filenamedel.str().c_str());
    if (screen_out) std::cout << "writing " << left << std::setw(50) <<filename.str()<<"...";
    std::ofstream f_system(filename.str().c_str());

    const XFEM::PHYSICS::Field field = XFEM::PHYSICS::Pres;

    {
      stringstream gmshfilecontent;
      gmshfilecontent << "View \" " << "Pressure Solution (Physical) \" {\n";
      for (int i=0; i<discret_->NumMyColElements(); ++i)
      {

        const DRT::Element* actele = discret_->lColElement(i);

        static LINALG::Matrix<3,27> xyze_xfemElement;
        GEO::fillInitialPositionArray(actele,xyze_xfemElement);

        // create local copy of information about dofs
        const XFEM::ElementDofManager eledofman(*actele,physprob_.elementAnsatzp_->getElementAnsatz(actele->Shape()),*dofmanagerForOutput_);

        vector<int> lm;
        vector<int> lmowner;
        actele->LocationVector(*(discret_), lm, lmowner);

        // extract local values from the global vector
        vector<double> myvelnp(lm.size());
        DRT::UTILS::ExtractMyValues(*state_.velnp_, myvelnp, lm);

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
          XFEM::computeScalarCellNodeValuesFromNodalUnknowns(*actele, dofmanagerForOutput_->getInterfaceHandle(), eledofman,
              *cell, field, elementvalues, cellvalues);
          //LINALG::SerialDenseMatrix xyze_cell(3, cell->NumNode());
          //cell->NodalPosXYZ(*actele, xyze_cell);
          // TODO remove
          const LINALG::SerialDenseMatrix& xyze_cell = cell->CellNodalPosXYZ();
          gmshfilecontent << IO::GMSH::cellWithScalarFieldToString(
              cell->Shape(), cellvalues, xyze_cell) << "\n";
        }
//        if (actele->Id() == 1 and elementvalues.Length() > 0)
//        {
//          //std::cout << elementvalues << "\n";
//          std::ofstream f;
//          const std::string fname = DRT::Problem::Instance()->OutputControlFile()->FileName()
//                                  + ".outflowpres.txt";
//          if (step <= 1)
//            f.open(fname.c_str(),std::fstream::trunc);
//          else
//            f.open(fname.c_str(),std::fstream::ate | std::fstream::app);
//
//          //f << time_ << " " << (-1.5*std::sin(0.1*2.0*time_* PI) * PI*0.1) << "  " << elementvalues(0,0) << endl;
//          f << time << "  " << elementvalues(0) << "\n";
//
//          f.close();
//        }
      }
      gmshfilecontent << "};\n";
      f_system << gmshfilecontent.str();
    }
    f_system.close();
    if (screen_out) std::cout << " done" << endl;
  }
  if (gmshdebugout and (this->physprob_.fieldset_.find(XFEM::PHYSICS::Temp) != this->physprob_.fieldset_.end()) )
  {
    std::stringstream filename;
    std::stringstream filenamedel;
    filename    << DRT::Problem::Instance()->OutputControlFile()->FileName() << ".solution_field_temperature_" << std::setw(5) << setfill('0') << step   << ".pos";
    filenamedel << DRT::Problem::Instance()->OutputControlFile()->FileName() << ".solution_field_temperature_" << std::setw(5) << setfill('0') << step-5 << ".pos";
    std::remove(filenamedel.str().c_str());
    if (screen_out) std::cout << "writing " << left << std::setw(50) <<filename.str()<<"...";
    std::ofstream f_system(filename.str().c_str());

    {
      const XFEM::PHYSICS::Field field = XFEM::PHYSICS::Temp;
      stringstream gmshfilecontent;
      gmshfilecontent << "View \" " << "Temperature Solution (Physical) \" {\n";
      for (int i=0; i<discret_->NumMyColElements(); ++i)
      {

        const DRT::Element* actele = discret_->lColElement(i);

        static LINALG::Matrix<3,27> xyze_xfemElement;
        GEO::fillInitialPositionArray(actele,xyze_xfemElement);

        // create local copy of information about dofs
        const XFEM::ElementDofManager eledofman(*actele,physprob_.elementAnsatzp_->getElementAnsatz(actele->Shape()),*dofmanagerForOutput_);

        vector<int> lm;
        vector<int> lmowner;
        actele->LocationVector(*(discret_), lm, lmowner);

        // extract local values from the global vector
        vector<double> myvelnp(lm.size());
        DRT::UTILS::ExtractMyValues(*state_.velnp_, myvelnp, lm);

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
          XFEM::computeScalarCellNodeValuesFromNodalUnknowns(*actele, dofmanagerForOutput_->getInterfaceHandle(), eledofman,
              *cell, field, elementvalues, cellvalues);
          const LINALG::SerialDenseMatrix& xyze_cell = cell->CellNodalPosXYZ();

          gmshfilecontent << IO::GMSH::cellWithScalarFieldToString(
              cell->Shape(), cellvalues, xyze_cell) << "\n";
        }
      }
      gmshfilecontent << "};\n";
      f_system << gmshfilecontent.str();
    }
    f_system.close();
    if (screen_out) std::cout << " done" << endl;
  }

#if 0
  if (gmshdebugout)
  {
    std::stringstream filename;
    std::stringstream filenamedel;
    const std::string filebase = DRT::Problem::Instance()->OutputControlFile()->FileName();
    filename    << filebase << ".solution_field_pressure_disc_" << std::setw(5) << setfill('0') << step   << ".pos";
    filenamedel << filebase << ".solution_field_pressure_disc_" << std::setw(5) << setfill('0') << step-5 << ".pos";
    std::remove(filenamedel.str().c_str());
    if (screen_out) std::cout << "writing " << std::left << std::setw(50) <<filename.str()<<"...";
    std::ofstream f_system(filename.str().c_str());

    const XFEM::PHYSICS::Field field = XFEM::PHYSICS::DiscPres;

    {
      stringstream gmshfilecontent;
      gmshfilecontent << "View \" " << "Discontinous Pressure Solution (Physical) \" {"
      << "\n";
      for (int i=0; i<discret_->NumMyColElements(); ++i)
      {
        const DRT::Element* actele = discret_->lColElement(i);

        static LINALG::Matrix<3,27> xyze_xfemElement;
        GEO::fillInitialPositionArray(actele,xyze_xfemElement);

        const map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> element_ansatz(XFLUID::getElementAnsatz(actele->Shape()));

        // create local copy of information about dofs
        const XFEM::ElementDofManager eledofman(*actele,element_ansatz,*dofmanagerForOutput_);

        vector<int> lm;
        vector<int> lmowner;
        actele->LocationVector(*(discret_), lm, lmowner);

        // extract local values from the global vector
        vector<double> myvelnp(lm.size());
        DRT::UTILS::ExtractMyValues(*state_.velnp_, myvelnp, lm);

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
          // LINALG::SerialDenseMatrix xyze_cell(3, cell->NumNode());
          //cell->NodalPosXYZ(*actele, xyze_cell);
          // TODO remove
          const LINALG::SerialDenseMatrix& xyze_cell = cell->CellNodalPosXYZ();
          gmshfilecontent << IO::GMSH::cellWithScalarFieldToString(
              cell->Shape(), cellvalues, xyze_cell) << "\n";
        }
      }
      gmshfilecontent << "};\n";
      f_system << gmshfilecontent.str();
    }
    f_system.close();
    if (screen_out) std::cout << " done" << endl;
  }
#endif
#if 0
  if (gmshdebugout)
  {
    std::stringstream filename;
    std::stringstream filenamexx;
    std::stringstream filenameyy;
    std::stringstream filenamezz;
    std::stringstream filenamexy;
    std::stringstream filenamexz;
    std::stringstream filenameyz;
    std::stringstream filenamedel;
    std::stringstream filenamexxdel;
    std::stringstream filenameyydel;
    std::stringstream filenamezzdel;
    std::stringstream filenamexydel;
    std::stringstream filenamexzdel;
    std::stringstream filenameyzdel;
    //filename   << "solution_tau_disc_"   << std::setw(5) << setfill('0') << step_ << ".pos";
    const std::string filebase = DRT::Problem::Instance()->OutputControlFile()->FileName();
    filename   << filebase << ".solution_field_sigma_disc_"   << std::setw(5) << setfill('0') << step << ".pos";
    filenamexx << filebase << ".solution_field_sigmaxx_disc_" << std::setw(5) << setfill('0') << step << ".pos";
    filenameyy << filebase << ".solution_field_sigmayy_disc_" << std::setw(5) << setfill('0') << step << ".pos";
    filenamezz << filebase << ".solution_field_sigmazz_disc_" << std::setw(5) << setfill('0') << step << ".pos";
    filenamexy << filebase << ".solution_field_sigmaxy_disc_" << std::setw(5) << setfill('0') << step << ".pos";
    filenamexz << filebase << ".solution_field_sigmaxz_disc_" << std::setw(5) << setfill('0') << step << ".pos";
    filenameyz << filebase << ".solution_field_sigmayz_disc_" << std::setw(5) << setfill('0') << step << ".pos";
    filenamedel   << filebase << ".solution_field_sigma_disc_"   << std::setw(5) << setfill('0') << step-5 << ".pos";
    filenamexxdel << filebase << ".solution_field_sigmaxx_disc_" << std::setw(5) << setfill('0') << step-5 << ".pos";
    filenameyydel << filebase << ".solution_field_sigmayy_disc_" << std::setw(5) << setfill('0') << step-5 << ".pos";
    filenamezzdel << filebase << ".solution_field_sigmazz_disc_" << std::setw(5) << setfill('0') << step-5 << ".pos";
    filenamexydel << filebase << ".solution_field_sigmaxy_disc_" << std::setw(5) << setfill('0') << step-5 << ".pos";
    filenamexzdel << filebase << ".solution_field_sigmaxz_disc_" << std::setw(5) << setfill('0') << step-5 << ".pos";
    filenameyzdel << filebase << ".solution_field_sigmayz_disc_" << std::setw(5) << setfill('0') << step-5 << ".pos";
    std::remove(filenamedel.str().c_str());
    std::remove(filenamexxdel.str().c_str());
    std::remove(filenameyydel.str().c_str());
    std::remove(filenamezzdel.str().c_str());
    std::remove(filenamexydel.str().c_str());
    std::remove(filenamexzdel.str().c_str());
    std::remove(filenameyzdel.str().c_str());
    if (screen_out) std::cout << "writing " << std::left << std::setw(50) <<"stresses"<<"..."<<flush;

    std::ofstream f_system(  filename.str().c_str());
    std::ofstream f_systemxx(filenamexx.str().c_str());
    std::ofstream f_systemyy(filenameyy.str().c_str());
    std::ofstream f_systemzz(filenamezz.str().c_str());
    std::ofstream f_systemxy(filenamexy.str().c_str());
    std::ofstream f_systemxz(filenamexz.str().c_str());
    std::ofstream f_systemyz(filenameyz.str().c_str());

    const XFEM::PHYSICS::Field field = XFEM::PHYSICS::Sigmaxx;

    {
      stringstream gmshfilecontent;
      stringstream gmshfilecontentxx;
      stringstream gmshfilecontentyy;
      stringstream gmshfilecontentzz;
      stringstream gmshfilecontentxy;
      stringstream gmshfilecontentxz;
      stringstream gmshfilecontentyz;
      gmshfilecontent   << "View \" " << "Discontinous Stress Solution (Physical) \" {" << endl;
      gmshfilecontentxx << "View \" " << "Discontinous Stress (xx) Solution (Physical) \" {\n";
      gmshfilecontentyy << "View \" " << "Discontinous Stress (yy) Solution (Physical) \" {\n";
      gmshfilecontentzz << "View \" " << "Discontinous Stress (zz) Solution (Physical) \" {\n";
      gmshfilecontentxy << "View \" " << "Discontinous Stress (xy) Solution (Physical) \" {\n";
      gmshfilecontentxz << "View \" " << "Discontinous Stress (xz) Solution (Physical) \" {\n";
      gmshfilecontentyz << "View \" " << "Discontinous Stress (yz) Solution (Physical) \" {\n";
      for (int i=0; i<discret_->NumMyColElements(); ++i)
      {
        const DRT::Element* actele = discret_->lColElement(i);

        static LINALG::Matrix<3,27> xyze_xfemElement;
        GEO::fillInitialPositionArray(actele,xyze_xfemElement);

        const map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> element_ansatz(XFLUID::getElementAnsatz(actele->Shape()));

        // create local copy of information about dofs
        const XFEM::ElementDofManager eledofman(*actele,element_ansatz,*dofmanagerForOutput_);

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
          gmshfilecontent << IO::GMSH::cellWithTensorFieldToString(cell->Shape(), cellvalues, xyze_cell) << endl;
          }

          {
          LINALG::SerialDenseVector cellvaluexx(DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
          XFEM::computeScalarCellNodeValuesFromElementUnknowns(*actele, dofmanagerForOutput_->getInterfaceHandle(), eledofman, *cell, field, elementvaluexx, cellvaluexx);
          gmshfilecontentxx << IO::GMSH::cellWithScalarFieldToString(cell->Shape(), cellvaluexx, xyze_cell) << "\n";
          }
          {
          LINALG::SerialDenseVector cellvalueyy(DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
          XFEM::computeScalarCellNodeValuesFromElementUnknowns(*actele, dofmanagerForOutput_->getInterfaceHandle(), eledofman, *cell, field, elementvalueyy, cellvalueyy);
          gmshfilecontentyy << IO::GMSH::cellWithScalarFieldToString(cell->Shape(), cellvalueyy, xyze_cell) << "\n";
          }
          {
          LINALG::SerialDenseVector cellvaluezz(DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
          XFEM::computeScalarCellNodeValuesFromElementUnknowns(*actele, dofmanagerForOutput_->getInterfaceHandle(), eledofman, *cell, field, elementvaluezz, cellvaluezz);
          gmshfilecontentzz << IO::GMSH::cellWithScalarFieldToString(cell->Shape(), cellvaluezz, xyze_cell) << "\n";
          }
          {
          LINALG::SerialDenseVector cellvaluexy(DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
          XFEM::computeScalarCellNodeValuesFromElementUnknowns(*actele, dofmanagerForOutput_->getInterfaceHandle(), eledofman, *cell, field, elementvaluexy, cellvaluexy);
          gmshfilecontentxy << IO::GMSH::cellWithScalarFieldToString(cell->Shape(), cellvaluexy, xyze_cell) << "\n";
          }
          {
          LINALG::SerialDenseVector cellvaluexz(DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
          XFEM::computeScalarCellNodeValuesFromElementUnknowns(*actele, dofmanagerForOutput_->getInterfaceHandle(), eledofman, *cell, field, elementvaluexz, cellvaluexz);
          gmshfilecontentxz << IO::GMSH::cellWithScalarFieldToString(cell->Shape(), cellvaluexz, xyze_cell) << "\n";
          }
          {
          LINALG::SerialDenseVector cellvalueyz(DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
          XFEM::computeScalarCellNodeValuesFromElementUnknowns(*actele, dofmanagerForOutput_->getInterfaceHandle(), eledofman, *cell, field, elementvalueyz, cellvalueyz);
          gmshfilecontentyz << IO::GMSH::cellWithScalarFieldToString(cell->Shape(), cellvalueyz, xyze_cell) << "\n";
          }
        }
      }
      gmshfilecontent   << "};" << endl;
      gmshfilecontentxx << "};\n";
      gmshfilecontentyy << "};\n";
      gmshfilecontentzz << "};\n";
      gmshfilecontentxy << "};\n";
      gmshfilecontentxz << "};\n";
      gmshfilecontentyz << "};\n";
      f_system   << gmshfilecontent.str();
      f_systemxx << gmshfilecontentxx.str();
      f_systemyy << gmshfilecontentyy.str();
      f_systemzz << gmshfilecontentzz.str();
      f_systemxy << gmshfilecontentxy.str();
      f_systemxz << gmshfilecontentxz.str();
      f_systemyz << gmshfilecontentyz.str();
    }
    if (screen_out) std::cout << " done" << endl;
  }
#endif

  if (this->physprob_.fieldset_.find(XFEM::PHYSICS::Velx) != this->physprob_.fieldset_.end())
  {
    PlotVectorFieldToGmsh(state_.velnp_, ".solution_field_velocity_np_","Velocity Solution (Physical) n+1",true, step, time);
    PlotVectorFieldToGmsh(state_.veln_,  ".solution_field_velocity_n_","Velocity Solution (Physical) n",false, step, time);
//    PlotVectorFieldToGmsh(state_.velnm_, ".solution_field_velocity_nm_","Velocity Solution (Physical) n-1",false, step, time);
    PlotVectorFieldToGmsh(state_.accnp_, ".solution_field_acceleration_np_","Acceleration Solution (Physical) n+1",false, step, time);
    PlotVectorFieldToGmsh(state_.accn_,  ".solution_field_acceleration_n_","Acceleration Solution (Physical) n",false, step, time);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FLD::XFluidImplicitTimeInt::PlotVectorFieldToGmsh(
    const Teuchos::RCP<Epetra_Vector>   vectorfield,
    const std::string filestr,
    const std::string name_in_gmsh,
    const bool plot_to_gnuplot,
    const int step,
    const double time
    ) const
{

  const Teuchos::ParameterList& xfemparams = DRT::Problem::Instance()->XFEMGeneralParams();
  const bool gmshdebugout = (bool)getIntegralValue<int>(xfemparams,"GMSH_DEBUG_OUT");

  const bool screen_out = false;

  if (gmshdebugout)
  {

    bool ele_to_textfile = false;
    std::stringstream filename;
    std::stringstream filenamedel;
    const std::string filebase = DRT::Problem::Instance()->OutputControlFile()->FileName();
    filename    << filebase << filestr << std::setw(5) << std::setfill('0') << step << ".pos";
    filenamedel << filebase << filestr << std::setw(5) << std::setfill('0') << step-5 << ".pos";
    std::remove(filenamedel.str().c_str());
    if (screen_out) std::cout << "writing " << std::left << std::setw(50) <<filename.str()<<"...";
    std::ofstream f_system(filename.str().c_str());

    {
      stringstream gmshfilecontent;
      gmshfilecontent << "View \" " << name_in_gmsh << "\" {\n";
      for (int i=0; i<discret_->NumMyColElements(); ++i)
      {
        const DRT::Element* actele = discret_->lColElement(i);
        const int elegid = actele->Id();

        static LINALG::Matrix<3,27> xyze_xfemElement;
        GEO::fillInitialPositionArray(actele,xyze_xfemElement);

        // create local copy of information about dofs
        const XFLUID::FluidElementAnsatz elementAnsatz;
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

//        if (!dofmanagerForOutput_->getInterfaceHandle()->ElementIntersected(elegid))
//        {
          const GEO::DomainIntCells& domainintcells =
            dofmanagerForOutput_->getInterfaceHandle()->GetDomainIntCells(actele);
          for (GEO::DomainIntCells::const_iterator cell =
            domainintcells.begin(); cell != domainintcells.end(); ++cell)
          {
            LINALG::SerialDenseMatrix cellvalues(3, DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
            //std::cout << cellvalues << endl;
            XFEM::computeVectorCellNodeValues(*actele, dofmanagerForOutput_->getInterfaceHandle(), eledofman,
                *cell, XFEM::PHYSICS::Velx, elementvalues, cellvalues);
            //LINALG::SerialDenseMatrix xyze_cell(3, cell->NumNode());
            //cell->NodalPosXYZ(*actele, xyze_cell);
            // TODO remove
            const LINALG::SerialDenseMatrix& xyze_cell = cell->CellNodalPosXYZ();
            gmshfilecontent << IO::GMSH::cellWithVectorFieldToString(
                cell->Shape(), cellvalues, xyze_cell) << "\n";
          }
//        }
//        else
//        {
          const GEO::BoundaryIntCells& boundaryintcells =
            dofmanagerForOutput_->getInterfaceHandle()->GetBoundaryIntCells(elegid);
          // draw boundary integration cells with values
          for (GEO::BoundaryIntCells::const_iterator cell =
            boundaryintcells.begin(); cell != boundaryintcells.end(); ++cell)
          {
            LINALG::SerialDenseMatrix cellvalues(3, DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
            //std::cout << cellvalues << endl;
            XFEM::computeVectorCellNodeValues(*actele, dofmanagerForOutput_->getInterfaceHandle(), eledofman,
                *cell, XFEM::PHYSICS::Velx, elementvalues, cellvalues);
            //LINALG::SerialDenseMatrix xyze_cell(3, cell->NumNode());
            //cell->NodalPosXYZ(*actele, xyze_cell);
            // TODO remove
            const LINALG::SerialDenseMatrix& xyze_cell = cell->CellNodalPosXYZ();
            gmshfilecontent << IO::GMSH::cellWithVectorFieldToString(
                cell->Shape(), cellvalues, xyze_cell) << "\n";
          }

          // draw uncutted element
          {
            LINALG::SerialDenseMatrix elevalues(3, DRT::UTILS::getNumberOfElementNodes(actele->Shape()),true);
            static LINALG::Matrix<3,27> xyze_ele;
            GEO::fillInitialPositionArray(actele, xyze_ele);
            const GEO::DomainIntCell cell(actele->Shape(),GEO::InitialPositionArray(actele));
            gmshfilecontent << IO::GMSH::cellWithVectorFieldToString(
                actele->Shape(), elevalues, xyze_ele) << "\n";
          }

//        }
        //if (dofmanagerForOutput_->getInterfaceHandle()->ElementIntersected(elegid) and not ele_to_textfile and ele_to_textfile2)
//        if (elegid == 14 and elementvalues.N() > 0 and plot_to_gnuplot)
        if (elegid == 1 and elementvalues.N() > 0 and plot_to_gnuplot)
        {
          ele_to_textfile = true;
          //std::cout << elementvalues << std::endl;
          std::ofstream f;
          const std::string fname = DRT::Problem::Instance()->OutputControlFile()->FileName()
                                  + ".outflowvel.txt";
          if (step <= 1)
            f.open(fname.c_str(),std::fstream::trunc);
          else
            f.open(fname.c_str(),std::fstream::ate | std::fstream::app);

          //f << time_ << " " << (-1.5*std::sin(0.1*2.0*time_* PI) * PI*0.1) << "  " << elementvalues(0,0) << endl;
          f << time << "  " << elementvalues(0,0) << "\n";

          f.close();
        }

      }
      gmshfilecontent << "};\n";
      f_system << gmshfilecontent.str();
    }
    f_system.close();
    if (screen_out) std::cout << " done" << endl;
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
void FLD::XFluidImplicitTimeInt::SetInitialFlowField(
    const Teuchos::RCP<DRT::Discretization> cutterdiscret,
    int whichinitialfield,
    int startfuncno
    )
{
  cout << "SetInitialFlowField" << endl;
  // create zero displacement vector to use initial position of interface
  {
    const Epetra_Map* fluidsurface_dofcolmap = cutterdiscret->DofColMap();
    Teuchos::RCP<Epetra_Vector> idispcolnp  = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
    Teuchos::RCP<Epetra_Vector> ivelcolnp   = LINALG::CreateVector(*fluidsurface_dofcolmap,true);

    Teuchos::RCP<Epetra_Vector> idispcoln   = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
    Teuchos::RCP<Epetra_Vector> ivelcoln    = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
    Teuchos::RCP<Epetra_Vector> ivelcolnm   = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
    Teuchos::RCP<Epetra_Vector> iacccoln    = LINALG::CreateVector(*fluidsurface_dofcolmap,true);

    cutterdiscret->SetState("idispcolnp",idispcolnp);
    cutterdiscret->SetState("ivelcolnp",ivelcolnp);

    cutterdiscret->SetState("idispcoln",idispcoln);
    cutterdiscret->SetState("ivelcoln",ivelcoln);
    cutterdiscret->SetState("ivelcolnm",ivelcolnm);
    cutterdiscret->SetState("iacccoln",iacccoln);

    ComputeInterfaceAndSetDOFs(cutterdiscret);
    cutterdiscret->ClearState();
  }

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

        double initialval=DRT::Problem::Instance()->Funct(startfuncno-1).Evaluate(index,lnode->X(),0.0,NULL);

        state_.velnp_->ReplaceGlobalValues(1,&initialval,&gid);
        state_.veln_ ->ReplaceGlobalValues(1,&initialval,&gid);
      }
    }

    // add random perturbation
    if(whichinitialfield==3)
    {
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
  }
  else
  {
    dserror("no other initial fields than zero, function and beltrami are available up to now");
  }

  return;
} // end SetInitialFlowField


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | evaluate error for test cases with analytical solutions   gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::XFluidImplicitTimeInt::EvaluateErrorComparedToAnalyticalSol()
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
void FLD::XFluidImplicitTimeInt::SolveStationaryProblem(
    const Teuchos::RCP<DRT::Discretization> cutterdiscret
    )
{

  const Epetra_Map* fluidsurface_dofcolmap = cutterdiscret->DofColMap();
  const Teuchos::RCP<Epetra_Vector> idispcolnp  = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
  const Teuchos::RCP<Epetra_Vector> idispcoln   = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
  const Teuchos::RCP<Epetra_Vector> ivelcolnp   = LINALG::CreateVector(*fluidsurface_dofcolmap,true); // one could give a velocity here to have stationary flow over the interface
  const Teuchos::RCP<Epetra_Vector> ivelcoln    = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
  const Teuchos::RCP<Epetra_Vector> ivelcolnm   = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
  const Teuchos::RCP<Epetra_Vector> iacccoln    = LINALG::CreateVector(*fluidsurface_dofcolmap,true);

  cutterdiscret->SetState("idispcolnp", idispcolnp);
  cutterdiscret->SetState("idispcoln", idispcoln);
  cutterdiscret->SetState("ivelcolnp",ivelcolnp);
  cutterdiscret->SetState("ivelcoln", ivelcoln);
  cutterdiscret->SetState("ivelcolnm",ivelcolnm);
  cutterdiscret->SetState("iacccoln", iacccoln);

  ComputeInterfaceAndSetDOFs(cutterdiscret);

  PrepareNonlinearSolve();

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

      // set vector values needed by elements
      discret_->ClearState();
      discret_->SetState("velnp",state_.velnp_);
      // predicted dirichlet values
      // velnp then also holds prescribed new dirichlet values
      discret_->EvaluateDirichlet(eleparams,state_.velnp_,null,null,null,dbcmaps_);
      discret_->ClearState();

      // evaluate Neumann b.c.
      neumann_loads_->PutScalar(0.0);
      discret_->EvaluateNeumann(eleparams,*neumann_loads_);
      discret_->ClearState();
    }

    // -------------------------------------------------------------------
    //                     solve nonlinear equation system
    // -------------------------------------------------------------------
    NonlinearSolve(cutterdiscret);

    // -------------------------------------------------------------------
    //                         output of solution
    // -------------------------------------------------------------------
    StatisticsAndOutput();

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
Teuchos::RCP<Epetra_Vector> FLD::XFluidImplicitTimeInt::CalcStresses()
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
FLD::XFluidImplicitTimeInt::~XFluidImplicitTimeInt()
{
  return;
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FLD::XFluidImplicitTimeInt::LiftDrag() const
{
  // in this map, the results of the lift drag calculation are stored
  RCP<map<int,vector<double> > > liftdragvals;

  FLD::UTILS::LiftDrag(*discret_,*trueresidual_,params_,liftdragvals);

  if (liftdragvals!=Teuchos::null and discret_->Comm().MyPID() == 0)
    FLD::UTILS::WriteLiftDragToFile(time_, step_, *liftdragvals);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FLD::XFluidImplicitTimeInt::ComputeSurfaceFlowRates() const
{

  const map<int,double> volumeflowratepersurface = FLD::UTILS::ComputeSurfaceFlowRates(*discret_, state_.velnp_);

  if (not volumeflowratepersurface.empty())
  {
    cout << "Number of flow rate conditions... " << volumeflowratepersurface.size() << endl;
  }

  double overall_flowrate = 0.0;
  std::map<int,double>::const_iterator entry;
  for(entry = volumeflowratepersurface.begin(); entry != volumeflowratepersurface.end(); ++entry )
  {
    const int condID = entry->first;
    const double value = entry->second;
    overall_flowrate += value;
    if (myrank_ == 0)
    {
      cout << " - flowrate for label " << condID << ":  " <<  scientific << value << endl;
    }
  }
  if (not volumeflowratepersurface.empty())
  {
    cout << " - flowrate over all boundaries: " << overall_flowrate << endl;
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FLD::XFluidImplicitTimeInt::ComputeSurfaceImpulsRates() const
{

  const map<int,LINALG::Matrix<3,1> > impulsratepersurface = FLD::UTILS::ComputeSurfaceImpulsRates(*discret_, state_.velnp_);

  if (not impulsratepersurface.empty())
  {
    cout << "Number of impuls rate conditions... " << impulsratepersurface.size() << endl;
  }

  LINALG::Matrix<3,1> overall_flowrate(true);
  for(std::map<int,LINALG::Matrix<3,1> >::const_iterator entry = impulsratepersurface.begin();
      entry != impulsratepersurface.end();
      ++entry )
  {
    const int condID = entry->first;
    const LINALG::Matrix<3,1> value = entry->second;
    overall_flowrate += value;
    if (myrank_ == 0)
    {
      cout << " - impulsrate for label " << condID << ":  " <<  value << endl;
    }
  }
  if (not impulsratepersurface.empty())
  {
    cout << " - impulsrate over all boundaries: " << overall_flowrate << endl;
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FLD::XFluidImplicitTimeInt::IntegrateInterfaceShape(std::string condname)
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


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::XFluidImplicitTimeInt::UseBlockMatrix(Teuchos::RCP<std::set<int> > condelements,
                                                const LINALG::MultiMapExtractor& domainmaps,
                                                const LINALG::MultiMapExtractor& rangemaps,
                                                bool splitmatrix)
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

void FLD::XFluidImplicitTimeInt::MonolithicMultiDisEvaluate(
    Teuchos::RCP<XFEM::InterfaceHandleXFSI> ih,
    const Teuchos::RCP<DRT::Discretization> fluiddis,
    const Teuchos::RCP<DRT::Discretization> ifacedis,
    Teuchos::ParameterList&                 params,
    Teuchos::RCP<LINALG::SparseOperator>    Cuu,
    Teuchos::RCP<LINALG::SparseOperator>    Mud,
    Teuchos::RCP<LINALG::SparseOperator>    Mdu,
    Teuchos::RCP<LINALG::SparseOperator>    Cdd,
    Teuchos::RCP<Epetra_Vector>             RHSd,
    Teuchos::RCP<LINALG::SparseOperator>    systemmatrix1,
    Teuchos::RCP<Epetra_Vector>             systemvector1)
{
  TEUCHOS_FUNC_TIME_MONITOR("FLD::XFluidImplicitTimeInt::MonolithicMultiDisEvaluate");

  if (!fluiddis->Filled()) dserror("fluiddis->FillComplete() was not called");
  if (!fluiddis->HaveDofs()) dserror("fluiddis->AssignDegreesOfFreedom() was not called");

  if (!fluiddis->Filled()) dserror("fluiddis->FillComplete() was not called");
  if (!fluiddis->HaveDofs()) dserror("fluiddis->AssignDegreesOfFreedom() was not called");

  // define element matrices and vectors
  RCP<Epetra_SerialDenseMatrix> elematrixCuu = rcp(new Epetra_SerialDenseMatrix());
  RCP<Epetra_SerialDenseMatrix> elematrixMud = rcp(new Epetra_SerialDenseMatrix());
  RCP<Epetra_SerialDenseMatrix> elematrixMdu = rcp(new Epetra_SerialDenseMatrix());
  RCP<Epetra_SerialDenseMatrix> elematrixCdd = rcp(new Epetra_SerialDenseMatrix());
  RCP<Epetra_SerialDenseVector> elevectorRHSd = rcp(new Epetra_SerialDenseVector());
  Epetra_SerialDenseMatrix elematrix1;
  Epetra_SerialDenseMatrix elematrix2;
  Epetra_SerialDenseVector elevector1;
  Epetra_SerialDenseVector elevector2;
  Epetra_SerialDenseVector elevector3;



  // loop over column elements
  const int numcolele = fluiddis->NumMyColElements();
  for (int i=0; i<numcolele; ++i)
  {
    DRT::Element* actele = fluiddis->lColElement(i);

    const bool intersected = ih->ElementIntersected(actele->Id());

    vector<int> fluidlm;
    vector<int> fluidlmowner;
    actele->LocationVector(*fluiddis, fluidlm, fluidlmowner);

    vector<int> ifacepatchlm;
    vector<int> ifacepatchlmowner;
    std::set<int> begids;

    if (intersected)
    {
      begids = ih->GetIntersectingBoundaryElementsGID(actele->Id());
      for (std::set<int>::const_iterator begid = begids.begin(); begid != begids.end(); ++begid)
      {
        DRT::Element* bele = ifacedis->gElement(*begid);

        vector<int> ifacelm;
        vector<int> ifacelmowner;
        bele->LocationVector(*ifacedis,ifacelm,ifacelmowner);

        ifacepatchlm.reserve( ifacepatchlm.size() + ifacelm.size());
        ifacepatchlm.insert( ifacepatchlm.end(), ifacelm.begin(), ifacelm.end());

        ifacepatchlmowner.reserve( ifacepatchlmowner.size() + ifacelmowner.size());
        ifacepatchlmowner.insert( ifacepatchlmowner.end(), ifacelmowner.begin(), ifacelmowner.end());
      }
    }

    // get dimension of element matrices and vectors
    // Reshape element matrices and vectors and init to zero
    const std::size_t fluideledim = fluidlm.size();
    const std::size_t ifacepatchdim = ifacepatchlm.size();

    if (intersected)
    {
      AdaptElementMatrix(fluideledim,   fluideledim,   *elematrixCuu);
      AdaptElementMatrix(fluideledim,   ifacepatchdim, *elematrixMud);
      AdaptElementMatrix(ifacepatchdim, fluideledim,   *elematrixMdu);
      AdaptElementMatrix(ifacepatchdim, ifacepatchdim, *elematrixCdd);
      if (elevectorRHSd->Length()!=(int)ifacepatchdim)
        elevectorRHSd->Size(ifacepatchdim);
      else
        memset(elevectorRHSd->Values(),0,ifacepatchdim*sizeof(double));
    }

    AdaptElementMatrix(fluideledim,   fluideledim,   elematrix1);

    if (elevector1.Length()!=(int)fluideledim)
      elevector1.Size(fluideledim);
    else
      memset(elevector1.Values(),0,fluideledim*sizeof(double));



    {
      TEUCHOS_FUNC_TIME_MONITOR("DRT::Discretization::Evaluate elements");

      if (intersected)
      {
        params.set("Cuu",elematrixCuu);
        params.set("Mud",elematrixMud);
        params.set("Mdu",elematrixMdu);
        params.set("Cdd",elematrixCdd);
        params.set("rhsd",elevectorRHSd);
        params.set("ifacepatchlm",&ifacepatchlm);
      }

      // call the element evaluate method
      const int err = actele->Evaluate(params,*fluiddis,fluidlm,elematrix1,elematrix2,
                                       elevector1,elevector2,elevector3);
      if (err) dserror("Proc %d: Element %d returned err=%d",fluiddis->Comm().MyPID(),actele->Id(),err);
    }

    {
      TEUCHOS_FUNC_TIME_MONITOR("DRT::Discretization::Evaluate assemble");

      if (intersected)
      {
        Cdd->Assemble(-1, *elematrixCdd, ifacepatchlm, ifacepatchlmowner, ifacepatchlm);
        Mdu->Assemble(-1, *elematrixMdu, ifacepatchlm, ifacepatchlmowner, fluidlm);
        Mud->Assemble(-1, *elematrixMud, fluidlm     , fluidlmowner     , ifacepatchlm);
        Cuu->Assemble(-1, *elematrixCuu, fluidlm     , fluidlmowner     , fluidlm);
        LINALG::Assemble(*RHSd,*elevectorRHSd,ifacepatchlm,ifacepatchlmowner);
      }
      systemmatrix1->Assemble(-1,elematrix1,fluidlm,fluidlmowner);
      LINALG::Assemble(*systemvector1,elevector1,fluidlm,fluidlmowner);
    }

  } // for (int i=0; i<numcolele; ++i)

  return;
}




#endif /* CCADISCRET       */
