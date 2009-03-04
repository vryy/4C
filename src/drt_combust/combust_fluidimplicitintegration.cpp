/*!----------------------------------------------------------------------
\file combust_fluidimplicitintegration.cpp
\brief not documented yet!

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <stdio.h>

#include "combust_fluidimplicitintegration.H"

#include "../drt_fluid/time_integration_scheme.H"
#include "../drt_io/io_control.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/linalg_ana.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_function.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_geometry/position_array.H"
#include "../drt_xfem/dof_management.H"
#include "../drt_xfem/dof_distribution_switcher.H"
#include "../drt_xfem/enrichment_utils.H"
#include "../drt_fluid/fluid_utils.H"
#include "combust3_interpolation.H"
#include "combust_interface.H" // nötig?
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "../drt_io/io_gmsh.H"
#include <Teuchos_TimeMonitor.hpp>


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
  output_ (output),
  myrank_(discret_->Comm().MyPID()),
  cout0_(discret_->Comm(), std::cout),
  time_(0.0),
  step_(0),
  stepmax_(params_.get<int>   ("max number timesteps")),
  maxtime_(params_.get<double>("total time")),
  timealgo_(params_.get<FLUID_TIMEINTTYPE>("time int algo")),
  extrapolationpredictor_(params.get("do explicit predictor",true)),
  uprestart_(params.get("write restart every", -1)),
  upres_(params.get("write solution every", -1)),
  writestresses_(params.get<int>("write stresses", 0))
{

  // time measurement: initialization
  TEUCHOS_FUNC_TIME_MONITOR(" + initialization");
  /*----------------------------------------------------------------------------------------------*
   * comment missing! Axels comment: get the basic parameters first
   *----------------------------------------------------------------------------------------------*/
  timealgo_ = params_.get<FLUID_TIMEINTTYPE>("time int algo");
  dtp_ = dta_ = params_.get<double>("time step size");
  theta_    = params_.get<double>("theta");
  // af-generalized-alpha parameters: gamma_ = 0.5 + alphaM_ - alphaF_
  // (may be reset below when starting algorithm is used)
  alphaM_   = params_.get<double>("alpha_M");
  alphaF_   = params_.get<double>("alpha_F");
  gamma_    = params_.get<double>("gamma");

  // not needed for COMBUST: create empty cutter discretization
//  Teuchos::RCP<DRT::Discretization> emptyboundarydis_ = DRT::UTILS::CreateDiscretizationFromCondition(
//      actdis, "schnackelzappel", "DummyBoundary", "BELE3", vector<string>(0));
//  Teuchos::RCP<Epetra_Vector> tmpdisp = LINALG::CreateVector(*emptyboundarydis_->DofRowMap(),true);
//  emptyboundarydis_->SetState("idispcolnp",tmpdisp);
//  emptyboundarydis_->SetState("idispcoln",tmpdisp);

  // soll der DofManager hier bleiben, oder soll er in den Combustion Algorithmus wandern? henke 10/08
  // apply enrichments
  std::set<XFEM::PHYSICS::Field> fieldset;
  fieldset.insert(XFEM::PHYSICS::Velx);
  fieldset.insert(XFEM::PHYSICS::Vely);
  fieldset.insert(XFEM::PHYSICS::Velz);
  fieldset.insert(XFEM::PHYSICS::Pres);
  const COMBUST::CombustElementAnsatz elementAnsatz;
  Teuchos::RCP<XFEM::DofManager> dofmanager = rcp(new XFEM::DofManager(null,fieldset,elementAnsatz,params_));
  /*----------------------------------------------------------------------------------------------*
   * comment missing! Axels comment: tell elements about the dofs and the integration  
   *----------------------------------------------------------------------------------------------*/
  {
    ParameterList eleparams;
    eleparams.set("action","store_xfem_info");
    eleparams.set("dofmanager",dofmanager);
    eleparams.set("interfacehandle",null);  // klären, wie das interfacehandle hier rein kommt!
    discret_->Evaluate(eleparams,null,null,null,null,null);
  }

  /*----------------------------------------------------------------------------------------------*
   * comment missing!
   *----------------------------------------------------------------------------------------------*/
  discret_->FillComplete();

  //output_.WriteMesh(0,0.0);

  // store a dofset with the complete fluid unknowns
  dofset_out_.Reset();
  dofset_out_.AssignDegreesOfFreedom(*discret_,0);
  // split based on complete fluid field
  FLD::UTILS::SetupFluidSplit(*discret_,dofset_out_,3,velpressplitterForOutput_);

  state_.nodalDofDistributionMap_.clear();
  state_.elementalDofDistributionMap_.clear();

  /*----------------------------------------------------------------------------------------------*
   * get density from elements
   *----------------------------------------------------------------------------------------------*/
  {
    ParameterList eleparams;
    eleparams.set("action","get_density");
    discret_->Evaluate(eleparams,null,null,null,null,null);
    density_ = eleparams.get<double>("density");
    if (density_ <= 0.0) dserror("received negative or zero density value from elements");
  }

} // FluidImplicitTimeInt::FluidImplicitTimeInt

/*------------------------------------------------------------------------------------------------*
 | destructor: Steht in xfluidimplicitintegration.cpp weiter unten!                   henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
FLD::CombustFluidImplicitTimeInt::~CombustFluidImplicitTimeInt()
{
  return;
}

/*------------------------------------------------------------------------------------------------*
 | Don't use! Switch over transient/stationary is in combust_dyn!                     henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::Integrate()
{
  dserror("Don't use this function! Its switch over transient/stationary scheme is in dynamic routine!");
/*
  // bound for the number of startsteps
  const int    numstasteps         =params_.get<int>   ("number of start steps");

  // output of stabilization details
  ParameterList *  stabparams=&(params_.sublist("STABILIZATION"));

  cout0_ << "Stabilization type         : " << stabparams->get<string>("STABTYPE") << "\n";
  cout0_ << "                             " << stabparams->get<string>("TDS")<< "\n";
  cout0_ << "\n";

  if(stabparams->get<string>("TDS") == "quasistatic")
  {
    if(stabparams->get<string>("TRANSIENT")=="yes_transient")
    {
      dserror("The quasistatic version of the residual-based stabilization currently does not support the incorporation of the transient term.");
    }
  }
  cout0_ <<  "                             " << "TRANSIENT       = " << stabparams->get<string>("TRANSIENT")      <<"\n";
  cout0_ <<  "                             " << "SUPG            = " << stabparams->get<string>("SUPG")           <<"\n";
  cout0_ <<  "                             " << "PSPG            = " << stabparams->get<string>("PSPG")           <<"\n";
  cout0_ <<  "                             " << "VSTAB           = " << stabparams->get<string>("VSTAB")          <<"\n";
  cout0_ <<  "                             " << "CSTAB           = " << stabparams->get<string>("CSTAB")          <<"\n";
  cout0_ <<  "                             " << "CROSS-STRESS    = " << stabparams->get<string>("CROSS-STRESS")   <<"\n";
  cout0_ <<  "                             " << "REYNOLDS-STRESS = " << stabparams->get<string>("REYNOLDS-STRESS")<<"\n";
  cout0_ << "\n";

  if (timealgo_==timeint_stationary)
    // stationary case
    SolveStationaryProblem(cutterdiscret);

  else  // instationary case
  {
    // start procedure
    if (step_<numstasteps)
    {
      if (numstasteps>stepmax_)
      {
        dserror("more startsteps than steps");
      }

      dserror("no starting steps supported");
    }

    // continue with the final time integration
    TimeLoop(cutterdiscret);
  }

  // print the results of time measurements
  TimeMonitor::summarize();

  return;
*/
} // FluidImplicitTimeInt::Integrate

/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::TimeLoop()
{

/*
 * This function is commented out because it should not be used by any combustion simulation.
 * The time loop to be called is the one in the class COMBUST::Algorithm. It accesses directly e.g.
 * NonlinearSolve in this class CombustFluidImplicitTimeInt.
 */

  dserror("Don't use this function! Use COMBUST::Algorithm::TimeLoop()");
/*
  // time measurement: time loop
  TEUCHOS_FUNC_TIME_MONITOR(" + time loop");

  // how do we want to solve or fluid equations?
  const int dyntype    =params_.get<int>   ("type of nonlinear solve");

  if (dyntype==1)
  {
//    if (alefluid_)
//      dserror("no ALE possible with linearised fluid");
//      additionally it remains to mention that for the linearised
//       fluid the stbilisation is hard coded to be SUPG/PSPG 
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
  cutterdiscret->SetState("iacccoln",iacccoln);

  while (step_<stepmax_ and time_<maxtime_)
  {
    PrepareTimeStep();
    // -------------------------------------------------------------------
    //                         out to screen
    // -------------------------------------------------------------------
    if (myrank_==0)
    {
      switch (timealgo_)
      {
      case timeint_one_step_theta:
        printf("TIME: %11.4E/%11.4E  DT = %11.4E  One-Step-Theta  STEP = %4d/%4d \n",
              time_,maxtime_,dta_,step_,stepmax_);
        break;
      case timeint_bdf2:
        printf("TIME: %11.4E/%11.4E  DT = %11.4E     BDF2         STEP = %4d/%4d \n",
               time_,maxtime_,dta_,step_,stepmax_);
        break;
      default:
        dserror("parameter out of range: IOP\n");
      } // end of switch(timealgo) 
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
      //LinearSolve(cutterdiscret,idispcol);
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
*/
} // FluidImplicitTimeInt::TimeLoop

/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::PrepareTimeStep()
{
  // -------------------------------------------------------------------
  //              set time dependent parameters
  // -------------------------------------------------------------------
  step_ += 1;
  time_ += dta_;

  // for bdf2 theta is set by the timestepsizes, 2/3 for const. dt
  if (timealgo_==timeint_bdf2)
  {
    // for bdf2 theta is set  by the timestepsizes, 2/3 for const. dt
    if (params_.get<FLUID_TIMEINTTYPE>("time int algo")==timeint_bdf2)
    {
      theta_ = (dta_+dtp_)/(2.0*dta_ + dtp_);
    }
    
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
}

/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::PrepareNonlinearSolve()
{

  // -------------------------------------------------------------------
  // set part of the rhs vector belonging to the old timestep
  // -------------------------------------------------------------------
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
      TIMEINT_THETA_BDF2::ExplicitPredictor(
          state_.veln_, state_.velnm_, state_.accn_,
              timealgo_, dta_, dtp_,
              state_.velnp_);
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
    discret_->EvaluateDirichlet(eleparams,state_.velnp_,null,null,null,dbcmaps_);
    discret_->ClearState();

    // evaluate Neumann conditions
    eleparams.set("total time",time_);
    eleparams.set("thsl",theta_*dta_);

    neumann_loads_->PutScalar(0.0);
    discret_->EvaluateNeumann(eleparams,*neumann_loads_);
    discret_->ClearState();
  }

}

/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 |
 | Within this routine, no parallel re-distribution is allowed to take place. Before and after this 
 | function, it's ok to do that.
 | Calling this function multiple times always results in the same solution vectors (axels comment)
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::IncorporateInterface(Teuchos::RCP<COMBUST::InterfaceHandleCombust> interfacehandle)
{
  // import information about interface from Adapter class ADAPTER::FluidCombust
  interfacehandle_ = interfacehandle;

//  std::cout << "tree after interfaceconstructor" << endl;
//  interfacehandle_->PrintTreeInformation(step_);
//  interfacehandle_->toGmsh(step_);

  // apply enrichments
  std::set<XFEM::PHYSICS::Field> fieldset;
  fieldset.insert(XFEM::PHYSICS::Velx);
  fieldset.insert(XFEM::PHYSICS::Vely);
  fieldset.insert(XFEM::PHYSICS::Velz);
  fieldset.insert(XFEM::PHYSICS::Pres);
  COMBUST::CombustElementAnsatz elementAnsatz;
  Teuchos::RCP<XFEM::DofManager> dofmanager = rcp(new XFEM::DofManager(interfacehandle_,fieldset,elementAnsatz,params_));
  
  // save dofmanager to be able to plot Gmsh stuff in Output()
  dofmanagerForOutput_ = dofmanager;

  // tell elements about the dofs and the integration
  {
      ParameterList eleparams;
      eleparams.set("action","store_xfem_info");
      eleparams.set("dofmanager",dofmanager);
      eleparams.set("interfacehandle",interfacehandle_);
      discret_->Evaluate(eleparams,null,null,null,null,null);
  }

  // print global and element dofmanager to Gmsh
//  dofmanager->toGmsh(step_);


  // store old (proc-overlapping) dofmap, compute new one and return it
  Epetra_Map olddofrowmap = *discret_->DofRowMap();
  discret_->FillComplete();
  Epetra_Map newdofrowmap = *discret_->DofRowMap();

  // print information about dofs
  const int numdof = newdofrowmap.NumGlobalElements();
  const int numnodaldof = dofmanager->NumNodalDof();
  cout0_ << "numdof = " << numdof << ", numstressdof = "<< (numdof - numnodaldof) << endl; 

  discret_->ComputeNullSpaceIfNecessary(solver_.Params());

  XFEM::NodalDofPosMap oldNodalDofDistributionMap(state_.nodalDofDistributionMap_);
  XFEM::ElementalDofPosMap oldElementalDofDistributionMap(state_.elementalDofDistributionMap_);
  dofmanager->fillDofDistributionMaps(
      state_.nodalDofDistributionMap_,
      state_.elementalDofDistributionMap_);

  // create switcher
  const XFEM::DofDistributionSwitcher dofswitch(
          interfacehandle_, dofmanager,
          olddofrowmap, newdofrowmap,
          oldNodalDofDistributionMap, state_.nodalDofDistributionMap_,
          oldElementalDofDistributionMap, state_.elementalDofDistributionMap_
          );

  // --------------------------------------------
  // switch state vectors to new dof distribution
  // --------------------------------------------

  // accelerations at time n and n-1
  dofswitch.mapVectorToNewDofDistribution(state_.accnp_);
  dofswitch.mapVectorToNewDofDistribution(state_.accn_);

  // velocities and pressures at time n+1, n and n-1
  dofswitch.mapVectorToNewDofDistribution(state_.velnp_); // use old velocity as start value
  dofswitch.mapVectorToNewDofDistribution(state_.veln_);
  dofswitch.mapVectorToNewDofDistribution(state_.velnm_);

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
    discret_->EvaluateDirichlet(eleparams, zeros_, Teuchos::null, Teuchos::null, 
                                Teuchos::null, dbcmaps_);
    zeros_->PutScalar(0.0); // just in case of change
  }

  neumann_loads_= LINALG::CreateVector(newdofrowmap,true);

  // ---------------------------------
  // Vectors used for solution process
  // ---------------------------------
  residual_     = LINALG::CreateVector(newdofrowmap,true);
  //trueresidual_ = LINALG::CreateVector(newdofrowmap,true);
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

  // This is a first estimate for the number of non zeros in a row of
  // the matrix. Assuming a structured 3d-fluid mesh we have 27 adjacent
  // nodes with 4 dofs each. (27*4=108)
  // We do not need the exact number here, just for performance reasons
  // a 'good' estimate

//  if (not params_.get<int>("Simple Preconditioner",0))
//  {
  // initialize standard (stabilized) system matrix
  sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(newdofrowmap,108,false,true));
//  }
//  else
//  {
//    const int numdim = params_.get<int>("number of velocity degrees of freedom");
//    Teuchos::RCP<LINALG::BlockSparseMatrix<VelPressSplitStrategy> > blocksysmat =
//      Teuchos::rcp(new LINALG::BlockSparseMatrix<VelPressSplitStrategy>(velpressplitter_,velpressplitter_,108,false,true));
//    blocksysmat->SetNumdim(numdim);
//    sysmat_ = blocksysmat;
//  }

}

/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::NonlinearSolve()
{
  //IncorporateInterface();

  PrepareNonlinearSolve();

  TEUCHOS_FUNC_TIME_MONITOR("   + nonlin. iteration/lin. solve");

  // ---------------------------------------------- nonlinear iteration
  // ------------------------------- stop nonlinear iteration when both
  //                                 increment-norms are below this bound
  const double  ittol     =params_.get<double>("tolerance for nonlin iter");

  //------------------------------ turn adaptive solver tolerance on/off
  const bool   isadapttol    = params_.get<bool>("ADAPTCONV",true);
  const double adaptolbetter = params_.get<double>("ADAPTCONV_BETTER",0.01);

  //const bool fluidrobin = params_.get<bool>("fluidrobin", false);

  int               itnum = 0;
  const int         itemax = params_.get<int>("max nonlin iter steps");
  bool              stopnonliniter = false;

  double dtsolve = 0.0;
  double dtele   = 0.0;

  // get new interface velocity
// const Teuchos::RCP<const Epetra_Vector> ivelcolnp = cutterdiscret->GetState("ivelcolnp");
// const Teuchos::RCP<const Epetra_Vector> ivelcoln  = cutterdiscret->GetState("ivelcoln");

  // rausgenommen, damit es kompiliert henke 14.11.08
  //if (myrank_ == 0 && ivelcolnp->MyLength() >= 3)
  {
//    std::cout << "applying interface velocity ivelcolnp[0] = " << (*ivelcolnp)[0] << std::endl;
//    std::cout << "applying interface velocity ivelcolnp[1] = " << (*ivelcolnp)[1] << std::endl;
//    std::cout << "applying interface velocity ivelcolnp[2] = " << (*ivelcolnp)[2] << std::endl;
    std::ofstream f;
    const std::string fname = DRT::Problem::Instance()->OutputControlFile()->FileName()
                            + ".outifacevelnp.txt";
    if (step_ <= 1)
      f.open(fname.c_str(),std::fstream::trunc);
    else
      f.open(fname.c_str(),std::fstream::ate | std::fstream::app);

//    f << time_ << " " << (*ivelcolnp)[0] << "  " << "\n";

    f.close();
  }

  // rausgenommen, damit es kompiliert henke 14.11.08
  //if (myrank_ == 0 && ivelcoln->MyLength() >= 3)
  {
//    std::cout << "applying interface velocity ivelcoln[0] = " << (*ivelcoln)[0] << std::endl;
//    std::cout << "applying interface velocity ivelcoln[1] = " << (*ivelcoln)[1] << std::endl;
//    std::cout << "applying interface velocity ivelcoln[2] = " << (*ivelcoln)[2] << std::endl;
    std::ofstream f;
    const std::string fname = DRT::Problem::Instance()->OutputControlFile()->FileName()
                            + ".outifaceveln.txt";
    if (step_ <= 1)
      f.open(fname.c_str(),std::fstream::trunc);
    else
      f.open(fname.c_str(),std::fstream::ate | std::fstream::app);

//    f << time_ << " " << (*ivelcoln)[0] << "  " << "\n";

    f.close();
  }

  // rausgenommen, damit es kompiliert henke 14.11.08
  //if (myrank_ == 0 && ivelcolnp->MyLength() >= 3)
  {
    std::ofstream f;
    const std::string fname = DRT::Problem::Instance()->OutputControlFile()->FileName()
                            + ".outifaceanalytischvel.txt";
    if (step_ <= 1)
      f.open(fname.c_str(),std::fstream::trunc);
    else
      f.open(fname.c_str(),std::fstream::ate | std::fstream::app);

    const double periodendauer = 10.0;
    f << time_ << " " << (-1.5*std::sin(2.0*time_* PI/periodendauer) * PI/periodendauer) << "\n";

    f.close();
  }

  // action for elements
  if (timealgo_!=timeint_stationary)
  {
    cout0_ << "******************************************************" << endl;
    cout0_ << "* Warning! Does not work for moving boundaries, yet! *" << endl;
    cout0_ << "******************************************************" << endl;
  }

  if (myrank_ == 0)
  {

    printf("+------------+-------------------+--------------+--------------+--------------+--------------+--------------+--------------+\n");
    printf("|- step/max -|- tol      [norm] -|-- vel-res ---|-- pre-res ---|-- fullres ---|-- vel-inc ---|-- pre-inc ---|-- fullinc ---|\n");
  }
  
  // this is a hack to make the code compile! There should be no cutterdis in here! henke 10/08
  Teuchos::RCP<DRT::Discretization>      cutterdiscret;
  cutterdiscret = null;
  const Epetra_Map* fluidsurface_dofcolmap = cutterdiscret->DofColMap();
  const Teuchos::RCP<Epetra_Vector> iforcecolnp = LINALG::CreateVector(*fluidsurface_dofcolmap,true);

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

      // create the parameters for the discretization
      ParameterList eleparams;

      // action for elements
      if (timealgo_==timeint_stationary)
        eleparams.set("action","calc_fluid_stationary_systemmat_and_residual");
      else
        eleparams.set("action","calc_fluid_systemmat_and_residual");

      // other parameters that might be needed by the elements
      //eleparams.set("total time",time_);
      //eleparams.set("thsl",theta_*dta_);
      eleparams.set("timealgo",timealgo_);
      eleparams.set("dt",dta_);
      eleparams.set("theta",theta_);
      eleparams.set("include reactive terms for linearisation",params_.get<bool>("Use reaction terms for linearisation",false));

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

      // give interface velocity to elements
//      eleparams.set("interface velocity",ivelcolnp);

      // reset interface force and let the elements fill it
      iforcecolnp->PutScalar(0.0);
      eleparams.set("interface force",iforcecolnp);

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

        // get physical surface force
        iforcecolnp->Scale(ResidualScaling());

        // finalize the complete matrix
        sysmat_->Complete();
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
    if (velnorm_L2 < 1e-5)
    {
      velnorm_L2 = 1.0;
    }
    if (prenorm_L2 < 1e-5)
    {
      prenorm_L2 = 1.0;
    }
    if (fullnorm_L2 < 1e-5)
    {
      fullnorm_L2 = 1.0;
    }

    //-------------------------------------------------- output to screen
    /* special case of very first iteration step:
        - solution increment is not yet available
        - convergence check is not required (we solve at least once!)    */
    if (itnum == 1)
    {
      if (myrank_ == 0)
      {
        printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   | %10.3E   |      --      |      --      |      --      |",
               itnum,itemax,ittol,vresnorm,presnorm,fullresnorm);
        printf(" (      --     ,te=%10.3E",dtele);
        printf(")");
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
          printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   | %10.3E   | %10.3E   | %10.3E   | %10.3E   |",
                 itnum,itemax,ittol,vresnorm,presnorm,fullresnorm,
                 incvelnorm_L2/velnorm_L2,incprenorm_L2/prenorm_L2,incfullnorm_L2/fullnorm_L2);
          printf(" (ts=%10.3E,te=%10.3E",dtsolve,dtele);
          printf(")\n");
          printf("+------------+-------------------+--------------+--------------+--------------+--------------+--------------+--------------+\n");

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
          printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   | %10.3E   | %10.3E   | %10.3E   | %10.3E   |",
                 itnum,itemax,ittol,vresnorm,presnorm,fullresnorm,
                 incvelnorm_L2/velnorm_L2,incprenorm_L2/prenorm_L2,incfullnorm_L2/fullnorm_L2);
          printf(" (ts=%10.3E,te=%10.3E",dtsolve,dtele);
          printf(")");
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

    //------------------------------------------------ update (u,p) trial
    state_.velnp_->Update(1.0,*incvel_,1.0);

    //cout << "*iforcecol" << endl;
    //cout << *iforcecol << endl;
  }

  // macht der FSI algorithmus
  iforcecolnp->Scale(-1.0);

//  cutterdiscret->SetState("iforcenp", iforcecolnp);

} // FluidImplicitTimeInt::NonlinearSolve

/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::Evaluate(Teuchos::RCP<const Epetra_Vector> vel)
{
  sysmat_->Zero();
  dserror("no monolithic FSI tested, check first!");
/*
  // set the new solution we just got
  if (vel!=Teuchos::null)
  {
    // Take Dirichlet values from velnp and add vel to veln for non-Dirichlet
    // values.
    Teuchos::RCP<Epetra_Vector> aux = LINALG::CreateVector(*(discret_->DofRowMap()),true);
    aux->Update(1.0, *state_.veln_, 1.0, *vel, 0.0);
    dbcmaps_->InsertOtherVector(dbcmaps_->ExtractOtherVector(aux), state_.velnp_);
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

  // other parameters that might be needed by the elements
  eleparams.set("total time",time_);
  eleparams.set("thsl",theta_*dta_);
  eleparams.set("dt",dta_);

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

  // call loop over elements
  discret_->Evaluate(eleparams,sysmat_,residual_);
  discret_->ClearState();

  // finalize the system matrix
  sysmat_->Complete();

  trueresidual_->Update(ResidualScaling(),*residual_,0.0);

  // Apply dirichlet boundary conditions to system of equations
  // residual displacements are supposed to be zero at boundary
  // conditions
  incvel_->PutScalar(0.0);
  LINALG::ApplyDirichlettoSystem(sysmat_,incvel_,residual_,zeros_,*(dbcmaps_->CondMap()));
  */
}

/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::TimeUpdate()
{

  // compute acceleration
  TIMEINT_THETA_BDF2::CalculateAcceleration(
      state_.velnp_, state_.veln_, state_.velnm_, state_.accn_,
          timealgo_, step_, theta_, dta_, dtp_,
          state_.accnp_);

  // update old acceleration
  state_.accn_->Update(1.0,*state_.accnp_,0.0);
  
  // velocities/pressures of this step become most recent 
  // velocities/pressures of the last step
  state_.velnm_->Update(1.0,*state_.veln_ ,0.0);
  state_.veln_ ->Update(1.0,*state_.velnp_,0.0);

  return;
}// FluidImplicitTimeInt::TimeUpdate

/*------------------------------------------------------------------------------------------------*
 | gammi 11/08 |
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

  //-------------------------------------------- output of solution

  if (step_%upres_ == 0)  //write solution
  {
    output_.NewStep    (step_,time_);
    //output_.WriteVector("velnp", state_.velnp_);
    std::set<XFEM::PHYSICS::Field> fields_out;
    fields_out.insert(XFEM::PHYSICS::Velx);
    fields_out.insert(XFEM::PHYSICS::Vely);
    fields_out.insert(XFEM::PHYSICS::Velz);
    fields_out.insert(XFEM::PHYSICS::Pres);
    Teuchos::RCP<Epetra_Vector> velnp_out = dofmanagerForOutput_->fillPhysicalOutputVector(
        *state_.velnp_, dofset_out_, state_.nodalDofDistributionMap_, fields_out);
    output_.WriteVector("velnp", velnp_out);
//    Teuchos::RCP<Epetra_Vector> accn_out = dofmanagerForOutput_->fillPhysicalOutputVector(
//        *state_.accn_, dofset_out_, nodalDofDistributionMap_, fields_out);
//    output_.WriteVector("accn", accn_out);

    // output (hydrodynamic) pressure
    Teuchos::RCP<Epetra_Vector> pressure = velpressplitterForOutput_.ExtractCondVector(velnp_out);
    pressure->Scale(density_);
    output_.WriteVector("pressure", pressure);

    //output_.WriteVector("residual", trueresidual_);

    //only perform stress calculation when output is needed
    if (writestresses_)
    {
      dserror("not supported, yet");
      Teuchos::RCP<Epetra_Vector> traction = CalcStresses();
      output_.WriteVector("traction",traction);
    }

    // write domain decomposition for visualization (only once!)
    if (step_==upres_)
     output_.WriteElementData();

    if (uprestart_ != 0 && step_%uprestart_ == 0) //add restart data
    {
      output_.WriteVector("accn", state_.accn_);
      output_.WriteVector("veln", state_.veln_);
      output_.WriteVector("velnm", state_.velnm_);
    }
  }

  // write restart also when uprestart_ is not a integer multiple of upres_
  else if (uprestart_ != 0 && step_%uprestart_ == 0)
  {
    output_.NewStep    (step_,time_);
    output_.WriteVector("velnp", state_.velnp_);
    //output_.WriteVector("residual", trueresidual_);

    //only perform stress calculation when output is needed
    if (writestresses_)
    {
      Teuchos::RCP<Epetra_Vector> traction = CalcStresses();
      output_.WriteVector("traction",traction);
    }

    output_.WriteVector("accn", state_.accn_);
    output_.WriteVector("veln", state_.veln_);
    output_.WriteVector("velnm", state_.velnm_);
  }

  if (discret_->Comm().NumProc() == 1)
  {
    OutputToGmsh();
  }
  return;
} // FluidImplicitTimeInt::Output

/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::ReadRestart(int step)
{
  IO::DiscretizationReader reader(discret_,step);
  time_ = reader.ReadDouble("time");
  step_ = reader.ReadInt("step");

  reader.ReadVector(state_.velnp_,"velnp");
  reader.ReadVector(state_.veln_, "veln");
  reader.ReadVector(state_.velnm_,"velnm");
  reader.ReadVector(state_.accn_ ,"accn");

}

/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::OutputToGmsh()
{
  const Teuchos::ParameterList& xfemparams = DRT::Problem::Instance()->XFEMGeneralParams();
  const bool gmshdebugout = (bool)getIntegralValue<int>(xfemparams,"GMSH_DEBUG_OUT");

  const bool screen_out = false;

  if (gmshdebugout)
  {
    cout << "CombustFluidImplicitTimeInt::OutputToGmsh()" << endl;

    std::stringstream filename;
    std::stringstream filenamedel;
    filename << DRT::Problem::Instance()->OutputControlFile()->FileName() << "_solution_pressure_" << std::setw(5) << setfill('0') << step_ << ".pos";
    filenamedel << DRT::Problem::Instance()->OutputControlFile()->FileName() << "_solution_pressure_" << std::setw(5) << setfill('0') << step_-5 << ".pos";
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
        const int elegid = actele->Id();

        LINALG::SerialDenseMatrix xyze_xfemElement(GEO::InitialPositionArray(actele));

        // create local copy of information about dofs
        const COMBUST::CombustElementAnsatz elementAnsatz;
        const XFEM::ElementDofManager eledofman(*actele,elementAnsatz.getElementAnsatz(actele->Shape()),*dofmanagerForOutput_);

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
          //LINALG::SerialDenseMatrix xyze_cell(3, cell->NumNode());
          //cell->NodalPosXYZ(*actele, xyze_cell);
          // TODO remove
          gmshfilecontent << IO::GMSH::cellWithScalarFieldToString(
              cell->Shape(), cellvalues, xyze_cell) << "\n";
        }
        if (elegid == 1 and elementvalues.Length() > 0)
        {
          //std::cout << elementvalues << "\n";
          std::ofstream f;
          const std::string fname = DRT::Problem::Instance()->OutputControlFile()->FileName()
                                  + ".outflowpres.txt";
          if (step_ <= 1)
            f.open(fname.c_str(),std::fstream::trunc);
          else
            f.open(fname.c_str(),std::fstream::ate | std::fstream::app);

          //f << time_ << " " << (-1.5*std::sin(0.1*2.0*time_* PI) * PI*0.1) << "  " << elementvalues(0,0) << endl;
          f << time_ << "  " << elementvalues(0) << "\n";

          f.close();
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
    filename << DRT::Problem::Instance()->OutputControlFile()->FileName() << "_solution_pressure_disc_" << std::setw(5) << setfill('0') << step_
    << ".pos";
    filenamedel << DRT::Problem::Instance()->OutputControlFile()->FileName() << "_solution_pressure_disc_" << std::setw(5) << setfill('0') << step_-5 << ".pos";
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
        const int elegid = actele->Id();

        LINALG::SerialDenseMatrix xyze_xfemElement(GEO::InitialPositionArray(actele));

        const map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> element_ansatz(COMBUST::getElementAnsatz(actele->Shape()));

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
        if(elementvalues.Length() != 0)
        {
          //cout << "eleval DiscPres" << endl;
          //cout << elementvalues << endl;
        }
        const GEO::DomainIntCells& domainintcells =
          dofmanagerForOutput_->getInterfaceHandle()->GetDomainIntCells(actele);
        for (GEO::DomainIntCells::const_iterator cell =
          domainintcells.begin(); cell != domainintcells.end(); ++cell)
        {
          LINALG::SerialDenseVector cellvalues(DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
          XFEM::computeScalarCellNodeValuesFromElementUnknowns(*actele, dofmanagerForOutput_->getInterfaceHandle(), eledofman,
              *cell, field, elementvalues, cellvalues);
          
          const LINALG::SerialDenseMatrix& xyze_cell = cell->CellNodalPosXYZ();
          //LINALG::SerialDenseMatrix xyze_cell(3, cell->NumNode());
          //cell->NodalPosXYZ(*actele, xyze_cell);
          //TODO remove
          if(elementvalues.Length() != 0)
          {
            //cout << "cellvalues DiscPres" << endl;
            //cout << cellvalues << endl;
          }
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
    //std::stringstream filename;
    std::stringstream filenamexx;
    std::stringstream filenameyy;
    std::stringstream filenamezz;
    std::stringstream filenamexy;
    std::stringstream filenamexz;
    std::stringstream filenameyz;
    std::stringstream filenamexxdel;
    std::stringstream filenameyydel;
    std::stringstream filenamezzdel;
    std::stringstream filenamexydel;
    std::stringstream filenamexzdel;
    std::stringstream filenameyzdel;
    //filename   << "solution_tau_disc_"   << std::setw(5) << setfill('0') << step_ << ".pos";
    const std::string filebase = DRT::Problem::Instance()->OutputControlFile()->FileName();
    filenamexx << filebase << "_solution_tauxx_disc_" << std::setw(5) << setfill('0') << step_ << ".pos";
    filenameyy << filebase << "_solution_tauyy_disc_" << std::setw(5) << setfill('0') << step_ << ".pos";
    filenamezz << filebase << "_solution_tauzz_disc_" << std::setw(5) << setfill('0') << step_ << ".pos";
    filenamexy << filebase << "_solution_tauxy_disc_" << std::setw(5) << setfill('0') << step_ << ".pos";
    filenamexz << filebase << "_solution_tauxz_disc_" << std::setw(5) << setfill('0') << step_ << ".pos";
    filenameyz << filebase << "_solution_tauyz_disc_" << std::setw(5) << setfill('0') << step_ << ".pos";
    filenamexxdel << filebase << "_solution_tauxx_disc_" << std::setw(5) << setfill('0') << step_-5 << ".pos";
    filenameyydel << filebase << "_solution_tauyy_disc_" << std::setw(5) << setfill('0') << step_-5 << ".pos";
    filenamezzdel << filebase << "_solution_tauzz_disc_" << std::setw(5) << setfill('0') << step_-5 << ".pos";
    filenamexydel << filebase << "_solution_tauxy_disc_" << std::setw(5) << setfill('0') << step_-5 << ".pos";
    filenamexzdel << filebase << "_solution_tauxz_disc_" << std::setw(5) << setfill('0') << step_-5 << ".pos";
    filenameyzdel << filebase << "_solution_tauyz_disc_" << std::setw(5) << setfill('0') << step_-5 << ".pos";
    std::remove(filenamexxdel.str().c_str());
    std::remove(filenameyydel.str().c_str());
    std::remove(filenamezzdel.str().c_str());
    std::remove(filenamexydel.str().c_str());
    std::remove(filenamexzdel.str().c_str());
    std::remove(filenameyzdel.str().c_str());
    if (screen_out) std::cout << "writing " << std::left << std::setw(50) <<"stresses"<<"..."<<flush;

    //std::ofstream f_system(  filename.str().c_str());
    std::ofstream f_systemxx(filenamexx.str().c_str());
    std::ofstream f_systemyy(filenameyy.str().c_str());
    std::ofstream f_systemzz(filenamezz.str().c_str());
    std::ofstream f_systemxy(filenamexy.str().c_str());
    std::ofstream f_systemxz(filenamexz.str().c_str());
    std::ofstream f_systemyz(filenameyz.str().c_str());

    const XFEM::PHYSICS::Field field = XFEM::PHYSICS::Sigmaxx;

    {
      //stringstream gmshfilecontent;
      stringstream gmshfilecontentxx;
      stringstream gmshfilecontentyy;
      stringstream gmshfilecontentzz;
      stringstream gmshfilecontentxy;
      stringstream gmshfilecontentxz;
      stringstream gmshfilecontentyz;
      //gmshfilecontent << "View \" " << "Discontinous Viscous Stress Solution (Physical) \" {" << endl;
      gmshfilecontentxx << "View \" " << "Discontinous Viscous Stress (xx) Solution (Physical) \" {\n";
      gmshfilecontentyy << "View \" " << "Discontinous Viscous Stress (yy) Solution (Physical) \" {\n";
      gmshfilecontentzz << "View \" " << "Discontinous Viscous Stress (zz) Solution (Physical) \" {\n";
      gmshfilecontentxy << "View \" " << "Discontinous Viscous Stress (xy) Solution (Physical) \" {\n";
      gmshfilecontentxz << "View \" " << "Discontinous Viscous Stress (xz) Solution (Physical) \" {\n";
      gmshfilecontentyz << "View \" " << "Discontinous Viscous Stress (yz) Solution (Physical) \" {\n";
      for (int i=0; i<discret_->NumMyColElements(); ++i)
      {

        const DRT::Element* actele = discret_->lColElement(i);
        const int elegid = actele->Id();

        LINALG::SerialDenseMatrix xyze_xfemElement(GEO::InitialPositionArray(actele));

        const map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> element_ansatz(COMBUST::getElementAnsatz(actele->Shape()));

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
//          {
//          LINALG::SerialDenseMatrix cellvalues(9,DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
//          XFEM::computeTensorCellNodeValuesFromElementUnknowns(*actele, dofmanagerForOutput_->getInterfaceHandle(), eledofman,
//              *cell, field, elementvalues, cellvalues);
//          gmshfilecontent << IO::GMSH::cellWithTensorFieldToString(cell->Shape(), cellvalues, xyze_cell) << endl;
//          }

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
      //gmshfilecontent   << "};" << endl;
      gmshfilecontentxx << "};\n";
      gmshfilecontentyy << "};\n";
      gmshfilecontentzz << "};\n";
      gmshfilecontentxy << "};\n";
      gmshfilecontentxz << "};\n";
      gmshfilecontentyz << "};\n";
      //f_system   << gmshfilecontent.str();
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


  PlotVectorFieldToGmsh(state_.velnp_, "_solution_velocity_","Velocity Solution (Physical) n+1",true);
//  PlotVectorFieldToGmsh(state_.veln_,  "_solution_velocity_old_step_","Velocity Solution (Physical) n",false);
//  PlotVectorFieldToGmsh(state_.velnm_, "_solution_velocity_old2_step_","Velocity Solution (Physical) n-1",false);
//  PlotVectorFieldToGmsh(state_.accn_,  "_solution_acceleration_old_step_","Acceleration Solution (Physical) n",false);
}

/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::PlotVectorFieldToGmsh(
    const Teuchos::RCP<Epetra_Vector>   vectorfield,
    const std::string filestr,
    const std::string name_in_gmsh,
    const bool plot_to_gnuplot
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
    filename << filebase << filestr << std::setw(5) << std::setfill('0') << step_ << ".pos";
    filenamedel << filebase << filestr << std::setw(5) << std::setfill('0') << step_-5 << ".pos";
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

        LINALG::SerialDenseMatrix xyze_xfemElement = GEO::InitialPositionArray(actele);
                      
        // create local copy of information about dofs
        const COMBUST::CombustElementAnsatz elementAnsatz;
        const XFEM::ElementDofManager eledofman(*actele,elementAnsatz.getElementAnsatz(actele->Shape()),*dofmanagerForOutput_);

        vector<int> lm;
        vector<int> lmowner;
        actele->LocationVector(*(discret_), lm, lmowner);

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

        if (!dofmanagerForOutput_->getInterfaceHandle()->ElementIntersected(elegid))
        {
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
        }
        else
        {
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
            const LINALG::SerialDenseMatrix xyze_ele(GEO::InitialPositionArray(actele));
            const GEO::DomainIntCell cell(actele->Shape(), xyze_ele);
            gmshfilecontent << IO::GMSH::cellWithVectorFieldToString(
                actele->Shape(), elevalues, xyze_ele) << "\n";
          }

        }
        //if (dofmanagerForOutput_->getInterfaceHandle()->ElementIntersected(elegid) and not ele_to_textfile and ele_to_textfile2)
        if (elegid == 1 and elementvalues.N() > 0 and plot_to_gnuplot)
        {
          ele_to_textfile = true;
          //std::cout << elementvalues << std::endl;
          std::ofstream f;
          const std::string fname = DRT::Problem::Instance()->OutputControlFile()->FileName()
                                  + ".outflowvel.txt";
          if (step_ <= 1)
            f.open(fname.c_str(),std::fstream::trunc);
          else
            f.open(fname.c_str(),std::fstream::ate | std::fstream::app);

          //f << time_ << " " << (-1.5*std::sin(0.1*2.0*time_* PI) * PI*0.1) << "  " << elementvalues(0,0) << endl;
          f << time_ << "  " << elementvalues(0,0) << "\n";

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

/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::SetInitialFlowField(
    Teuchos::RCP<DRT::Discretization> cutterdiscret,
    int whichinitialfield,
    int startfuncno)
{
  /* Die member Funktion in CombustFluid ist nicht implementiert und erzeugt einen dserror!
   * Es kann auf diese Funktion eigentlich nur direkt zugegriffen werden (ist das so? warum?)*/

  // create zero displacement vector to use initial position of interface
  {
    const Epetra_Map* fluidsurface_dofcolmap = cutterdiscret->DofColMap();
    Teuchos::RCP<Epetra_Vector> idispcolnp  = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
    Teuchos::RCP<Epetra_Vector> ivelcolnp   = LINALG::CreateVector(*fluidsurface_dofcolmap,true);

    Teuchos::RCP<Epetra_Vector> idispcoln   = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
    Teuchos::RCP<Epetra_Vector> ivelcoln    = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
    Teuchos::RCP<Epetra_Vector> ivelcolnm   = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
    Teuchos::RCP<Epetra_Vector> iacccoln    = LINALG::CreateVector(*fluidsurface_dofcolmap,true);

    //IncorporateInterface();
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

        double initialval=DRT::UTILS::FunctionManager::Instance().Funct(startfuncno-1).Evaluate(index,lnode->X());

        state_.velnp_->ReplaceGlobalValues(1,&initialval,&gid);
        state_.veln_ ->ReplaceGlobalValues(1,&initialval,&gid);
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

/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void FLD::CombustFluidImplicitTimeInt::EvaluateErrorComparedToAnalyticalSol()
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
/*
  const Epetra_Map* fluidsurface_dofcolmap = cutterdiscret->DofColMap();
  const Teuchos::RCP<Epetra_Vector> idispcolnp  = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
  const Teuchos::RCP<Epetra_Vector> idispcoln   = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
  const Teuchos::RCP<Epetra_Vector> ivelcolnp   = LINALG::CreateVector(*fluidsurface_dofcolmap,true); // one could give a velocity here to have stationary flow over the interface
  const Teuchos::RCP<Epetra_Vector> ivelcoln    = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
  const Teuchos::RCP<Epetra_Vector> ivelcolnm   = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
  const Teuchos::RCP<Epetra_Vector> iacccoln    = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
  const Teuchos::RCP<Epetra_Vector> itruerescol = LINALG::CreateVector(*fluidsurface_dofcolmap,true);

  IncorporateInterface();

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
    NonlinearSolve();

    // -------------------------------------------------------------------
    //                         output of solution
    // -------------------------------------------------------------------
    StatisticsAndOutput();

  } // end of pseudo time loop
*/
} // CombustFluidImplicitTimeInt::SolveStationaryProblem

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
    //meshmovematrix_ = mat;
  }
}

#endif // CCADISCRET
