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

#include "xfluidimplicitintegration.H"
#include "time_integration_scheme.H"

#include "../drt_lib/linalg_ana.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_function.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_xfem/interface.H"
#include "../drt_xfem/dof_management.H"
#include "../drt_xfem/dof_distribution_switcher.H"
#include "../drt_xfem/field_enriched.H"
#include "../drt_xfem/enrichment_utils.H"
#include "fluid_utils.H"
#include "../drt_f3/xfluid3_interpolation.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "../drt_io/io_gmsh.H"
//#include <ctime>
#include <Teuchos_TimeMonitor.hpp>


extern "C" /* stuff which is c and is accessed from c++ */
{
#include "../headers/standardtypes.h"
}
extern struct _FILES  allfiles;

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Constructor (public)                                     gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
FLD::XFluidImplicitTimeInt::XFluidImplicitTimeInt(
    Teuchos::RCP<DRT::Discretization> actdis,
    LINALG::Solver&                   solver,
    ParameterList&                    params,
    IO::DiscretizationWriter&         output,
    const bool                        alefluid
    ) :
  // call constructor for "nontrivial" objects
  discret_(actdis),
  solver_ (solver),
  params_ (params),
  output_ (output),
  myrank_(discret_->Comm().MyPID()),
  alefluid_(alefluid),
  time_(0.0),
  step_(0),
  stepmax_(params_.get<int>   ("max number timesteps")),
  maxtime_(params_.get<double>("total time")),
  extrapolationpredictor_(params.get("do explicit predictor",true)),
  uprestart_(params.get("write restart every", -1)),
  upres_(params.get("write solution every", -1)),
  writestresses_(params.get<int>("write stresses", 0))
{

  // time measurement: initialization
  TEUCHOS_FUNC_TIME_MONITOR(" + initialization");

  // -------------------------------------------------------------------
  // get the basic parameters first
  // -------------------------------------------------------------------
  timealgo_ = params_.get<FLUID_TIMEINTTYPE>("time int algo");
  dtp_ = dta_ = params_.get<double>("time step size");
  theta_    = params_.get<double>("theta");

  // create empty cutter discretization
  Teuchos::RCP<DRT::Discretization> emptyboundarydis_ = DRT::UTILS::CreateDiscretizationFromCondition(
      actdis, "schnackelzappel", "DummyBoundary", "BELE3", vector<string>(0));
  Teuchos::RCP<Epetra_Vector> tmpdisp = LINALG::CreateVector(*emptyboundarydis_->DofRowMap(),true);
  emptyboundarydis_->SetState("idispcolnp",tmpdisp);
  emptyboundarydis_->SetState("idispcoln",tmpdisp);
  // intersection with empty cutter will result in a complete fluid domain with no holes or intersections
  Teuchos::RCP<XFEM::InterfaceHandle> ih = rcp(new XFEM::InterfaceHandle(discret_,emptyboundarydis_));
  // apply enrichments
  Teuchos::RCP<XFEM::DofManager> dofmanager = rcp(new XFEM::DofManager(ih));
  // tell elements about the dofs and the integration
  {
    ParameterList eleparams;
    eleparams.set("action","store_xfem_info");
    eleparams.set("dofmanager",dofmanager);
    eleparams.set("interfacehandle",ih);
    eleparams.set("assemble matrix 1",false);
    eleparams.set("assemble matrix 2",false);
    eleparams.set("assemble vector 1",false);
    eleparams.set("assemble vector 2",false);
    eleparams.set("assemble vector 3",false);
    discret_->Evaluate(eleparams,null,null,null,null,null);
  }
  discret_->FillComplete();

  //output_.WriteMesh(0,0.0);

  // store a dofset with the complete fluid unknowns
  dofset_out_.Reset();
  dofset_out_.AssignDegreesOfFreedom(*discret_,0);
  // split based on complete fluid field
  FLD::UTILS::SetupFluidSplit(*discret_,dofset_out_,3,velpressplitterForOutput_);

  state_.nodalDofDistributionMap_.clear();
  state_.elementalDofDistributionMap_.clear();
  
  // get density from elements
  {
    ParameterList eleparams;
    eleparams.set("action","get_density");
    discret_->Evaluate(eleparams,null,null,null,null,null);
    density_ = eleparams.get("density", 0.0);
    if (density_ <= 0.0) dserror("received negative density value from elements");
  }

} // FluidImplicitTimeInt::FluidImplicitTimeInt


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::XFluidImplicitTimeInt::Integrate(
    Teuchos::RCP<DRT::Discretization> cutterdiscret
        )
{
  // bound for the number of startsteps
  const int    numstasteps         =params_.get<int>   ("number of start steps");

  // output of stabilization details
  if (myrank_==0)
  {
    ParameterList *  stabparams=&(params_.sublist("STABILIZATION"));

    std::cout << "Stabilization type         : " << stabparams->get<string>("STABTYPE") << "\n";
    std::cout << "                             " << stabparams->get<string>("TDS")<< "\n";
    std::cout << "\n";

    if(stabparams->get<string>("TDS") == "quasistatic")
    {
      if(stabparams->get<string>("TRANSIENT")=="yes_transient")
      {
        dserror("The quasistatic version of the residual-based stabilization currently does not support the incorporation of the transient term.");
      }
    }
    std::cout <<  "                             " << "TRANSIENT       = " << stabparams->get<string>("TRANSIENT")      <<"\n";
    std::cout <<  "                             " << "SUPG            = " << stabparams->get<string>("SUPG")           <<"\n";
    std::cout <<  "                             " << "PSPG            = " << stabparams->get<string>("PSPG")           <<"\n";
    std::cout <<  "                             " << "VSTAB           = " << stabparams->get<string>("VSTAB")          <<"\n";
    std::cout <<  "                             " << "CSTAB           = " << stabparams->get<string>("CSTAB")          <<"\n";
    std::cout <<  "                             " << "CROSS-STRESS    = " << stabparams->get<string>("CROSS-STRESS")   <<"\n";
    std::cout <<  "                             " << "REYNOLDS-STRESS = " << stabparams->get<string>("REYNOLDS-STRESS")<<"\n";
    std::cout << "\n";
  }

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
} // FluidImplicitTimeInt::Integrate



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::XFluidImplicitTimeInt::TimeLoop(
    Teuchos::RCP<DRT::Discretization> cutterdiscret
        )
{
  // time measurement: time loop
  TEUCHOS_FUNC_TIME_MONITOR(" + time loop");

  // how do we want to solve or fluid equations?
  const int dyntype    =params_.get<int>   ("type of nonlinear solve");

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

    // time measurement: output and statistics
    TEUCHOS_FUNC_TIME_MONITOR("      + output and statistics");

    // -------------------------------------------------------------------
    // evaluate error for test flows with analytical solutions
    // -------------------------------------------------------------------
    EvaluateErrorComparedToAnalyticalSol();

    // -------------------------------------------------------------------
    //                         output of solution
    // -------------------------------------------------------------------
    Output();

    // -------------------------------------------------------------------
    //                    calculate lift'n'drag forces
    // -------------------------------------------------------------------
    const int liftdrag = params_.get<int>("liftdrag");

    if (liftdrag == 0); // do nothing, we don't want lift & drag
    if (liftdrag == 1)
      dserror("how did you manage to get here???");
    if (liftdrag == 2)
      LiftDrag();

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
void FLD::XFluidImplicitTimeInt::PrepareTimeStep()
{
  // -------------------------------------------------------------------
  //              set time dependent parameters
  // -------------------------------------------------------------------
  step_ += 1;
  time_ += dta_;

  // for bdf2 theta is set  by the timestepsizes, 2/3 for const. dt
  if (timealgo_==timeint_bdf2)
  {
    theta_ = (dta_+dtp_)/(2.0*dta_ + dtp_);
  }
  
  // do a backward Euler step for the first timestep
  if (step_==1)
  {
    theta_ = 1.0;
  }
  else
  {
    theta_ = params_.get<double>("theta");
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
              dta_, dtp_,
              state_.velnp_);
    }
  }

  // -------------------------------------------------------------------
  //         evaluate dirichlet and neumann boundary conditions
  // -------------------------------------------------------------------
  {
    ParameterList eleparams;

    // choose what to assemble
    eleparams.set("assemble matrix 1",false);
    eleparams.set("assemble matrix 2",false);
    eleparams.set("assemble vector 1",true);
    eleparams.set("assemble vector 2",false);
    eleparams.set("assemble vector 3",false);
    // other parameters needed by the elements
    eleparams.set("total time",time_);
    eleparams.set("delta time",dta_);
    eleparams.set("thsl",theta_*dta_);

    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("velnp",state_.velnp_);
    // predicted dirichlet values
    // velnp then also holds prescribed new dirichlet values
    // dirichtoggle is 1 for dirichlet dofs, 0 elsewhere
    discret_->EvaluateDirichlet(eleparams,state_.velnp_,null,null,dirichtoggle_);
    discret_->ClearState();

    // evaluate Neumann conditions
    eleparams.set("total time",time_);
    eleparams.set("thsl",theta_*dta_);

    neumann_loads_->PutScalar(0.0);
    discret_->EvaluateNeumann(eleparams,*neumann_loads_);
    discret_->ClearState();
  }

  //----------------------- compute an inverse of the dirichtoggle vector
  invtoggle_->PutScalar(1.0);
  invtoggle_->Update(-1.0,*dirichtoggle_,1.0);
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::XFluidImplicitTimeInt::ComputeInterfaceAndSetDOFs(
    Teuchos::RCP<DRT::Discretization>  cutterdiscret
    )
{
  // within this routine, no parallel re-distribution is allowed to take place
  // before and after this function, it's ok to do that

  // calling this function multiple times always results in the same solution vectors

  // compute Intersection
  Teuchos::RCP<XFEM::InterfaceHandle> ih = rcp(new XFEM::InterfaceHandle(discret_, cutterdiscret));
//  cout << "tree after interfaceconstructor" << endl;
//  ih->PrintTreeInformation(step_);
  ih->toGmsh(step_);

  // apply enrichments
  Teuchos::RCP<XFEM::DofManager> dofmanager = rcp(new XFEM::DofManager(ih));

  // save to be able to plot Gmsh stuff in Output()
  dofmanagerForOutput_ = dofmanager;

  // tell elements about the dofs and the integration
  {
      ParameterList eleparams;
      eleparams.set("action","store_xfem_info");
      eleparams.set("dofmanager",dofmanager);
      eleparams.set("interfacehandle",ih);
      eleparams.set("assemble matrix 1",false);
      eleparams.set("assemble matrix 2",false);
      eleparams.set("assemble vector 1",false);
      eleparams.set("assemble vector 2",false);
      eleparams.set("assemble vector 3",false);
      discret_->Evaluate(eleparams,null,null,null,null,null);
  }

  // print global and element dofmanager to Gmsh
  dofmanager->toGmsh(ih, step_);


  // store old (proc-overlapping) dofmap, compute new one and return it
  Epetra_Map olddofrowmap = *discret_->DofRowMap();
  discret_->FillComplete();
  Epetra_Map newdofrowmap = *discret_->DofRowMap();

  discret_->ComputeNullSpaceIfNecessary(solver_.Params());

  XFEM::NodalDofPosMap oldNodalDofDistributionMap(state_.nodalDofDistributionMap_);
  XFEM::ElementalDofPosMap oldElementalDofDistributionMap(state_.elementalDofDistributionMap_);
  dofmanager->fillDofDistributionMaps(
      state_.nodalDofDistributionMap_,
      state_.elementalDofDistributionMap_);

  std::cout << "switching " << endl;

  // create switcher
  const XFEM::DofDistributionSwitcher dofswitch(
          ih, dofmanager,
          olddofrowmap, newdofrowmap,
          oldNodalDofDistributionMap, state_.nodalDofDistributionMap_,
          oldElementalDofDistributionMap, state_.elementalDofDistributionMap_
          );

  // --------------------------------------------
  // switch state vectors to new dof distribution
  // --------------------------------------------

  // rigid body hack - assume structure is rigid and has uniform acceleration and velocity
  const Epetra_Vector& ivelcoln = *cutterdiscret->GetState("ivelcoln");
  const Epetra_Vector& iacccoln = *cutterdiscret->GetState("iacccoln");
  
  BlitzVec3 rigidveln;
  BlitzVec3 rigidaccn;

  if (ivelcoln.GlobalLength() > 3)
  {
    rigidveln(0) = (ivelcoln)[0];
    rigidveln(1) = (ivelcoln)[1];
    rigidveln(2) = (ivelcoln)[2];
    
    // falsch!!!
    rigidaccn(0) = (iacccoln)[0];
    rigidaccn(1) = (iacccoln)[1];
    rigidaccn(2) = (iacccoln)[2];
    
  }
  else
  {
    std::cout << "Could not compute rigid body velocity/acceleration. Set them to zero..." << std::endl;
    rigidveln = 0.0;
    rigidaccn = 0.0;
  }
  cout << "rigidveln " << rigidveln << endl;
  cout << "rigidaccn " << rigidaccn << endl;


  // accelerations at time n and n-1
  dofswitch.mapVectorToNewDofDistribution(state_.accn_, rigidaccn);

  // velocities and pressures at time n+1, n and n-1
  dofswitch.mapVectorToNewDofDistribution(state_.velnp_, rigidveln); // use old velocity as start value
  dofswitch.mapVectorToNewDofDistribution(state_.veln_, rigidveln);
  dofswitch.mapVectorToNewDofDistribution(state_.velnm_);

//  for (int i=0; i<ih->xfemdis()->NumMyColNodes(); ++i)
//  {
//    const DRT::Node* node = ih->xfemdis()->lColNode(i);
//    const BlitzVec3 nodalpos(toBlitzArray(node->X()));
//
//    bool is_in_fluid;
//    if (ih->PositionWithinConditionNP(nodalpos) == 0)
//    {
//      is_in_fluid = true;
//    }
//    else
//    {
//      is_in_fluid = false;
//    }
//    
//    bool was_in_fluid = false;
//    if (ih->PositionWithinConditionN(nodalpos) == 0)
//    {
//      was_in_fluid = true;
//    }
//    else
//    {
//      was_in_fluid = false;
//    }
//    
//    if (not is_in_fluid)
//    {
//      std::set<XFEM::FieldEnr> nodalDofSet = dofmanager->getNodeDofSet(node->Id());
//      if (not nodalDofSet.empty())
//      {
//        if (nodalDofSet.size() == 4)
//        {
//          for (std::set<XFEM::FieldEnr>::const_iterator iter=nodalDofSet.begin();iter!=nodalDofSet.end();iter++)
//          {
//            XFEM::Enrichment enr = iter->getEnrichment();
//            XFEM::PHYSICS::Field field = iter->getField();
//            XFEM::DofKey<XFEM::onNode> newdofkey(node->Id(),*iter);
//            
//            std::map<XFEM::DofKey<XFEM::onNode>, XFEM::DofPos>::const_iterator newdof = state_.nodalDofDistributionMap_.find(newdofkey);
//            if (newdof == state_.nodalDofDistributionMap_.end())
//            {
//              dserror("should definitely be there!");
//            }
//            const int newdofpos = newdof->second;
//            
//            if (field == XFEM::PHYSICS::Velx)
//            {
//              (*state_.veln_)[newdofrowmap.LID(newdofpos)] = rigidveln(0);
//              (*state_.accn_)[newdofrowmap.LID(newdofpos)] = rigidaccn(0);
//            }
//            else if (field == XFEM::PHYSICS::Vely)
//            {
//              (*state_.veln_)[newdofrowmap.LID(newdofpos)] = rigidveln(1);
//              (*state_.accn_)[newdofrowmap.LID(newdofpos)] = rigidaccn(1);
//            }
//            else if (field == XFEM::PHYSICS::Velz)
//            {
//              (*state_.veln_)[newdofrowmap.LID(newdofpos)] = rigidveln(2);
//              (*state_.accn_)[newdofrowmap.LID(newdofpos)] = rigidaccn(2);
//            }
//            else
//            {
//              (*state_.veln_)[newdofrowmap.LID(newdofpos)] = 0.0;
//            }
//  
//          }
//        }
//        else
//        {
//          dserror("works only for 4 dofs per node");
//        }
//      }
//      
//    }
//
//    
//  };
  
//  if (alefluid_)
//  {
//      dofswitch.mapVectorToNewDofDistribution(state_.dispnp_);
//      dofswitch.mapVectorToNewDofDistribution(state_.dispn_);
//      dofswitch.mapVectorToNewDofDistribution(state_.dispnm_);
//      dofswitch.mapVectorToNewDofDistribution(gridv_);
//  }

  // --------------------------------------------------
  // create remaining vectors with new dof distribution
  // --------------------------------------------------
  //hist_         = LINALG::CreateVector(newdofrowmap,true);

//  gridv_        = LINALG::CreateVector(newdofrowmap,true);

  dirichtoggle_ = LINALG::CreateVector(newdofrowmap,true);
  invtoggle_    = LINALG::CreateVector(newdofrowmap,true);

  zeros_        = LINALG::CreateVector(newdofrowmap,true);

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


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::XFluidImplicitTimeInt::NonlinearSolve(
    Teuchos::RCP<DRT::Discretization> cutterdiscret
    )
{

  ComputeInterfaceAndSetDOFs(cutterdiscret);

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
  bool              stopnonliniter = false;

  const int itemax  = params_.get<int>   ("max nonlin iter steps");

  double dtsolve = 0.0;
  double dtele   = 0.0;

  // get new interface velocity
  const Teuchos::RCP<const Epetra_Vector> ivelcolnp = cutterdiscret->GetState("ivelcolnp");
  const Teuchos::RCP<const Epetra_Vector> ivelcoln  = cutterdiscret->GetState("ivelcoln");
  
  if (myrank_ == 0 && ivelcolnp->MyLength() >= 3)
  {
    std::cout << "applying interface velocity ivelcolnp[0] = " << (*ivelcolnp)[0] << std::endl;
    std::cout << "applying interface velocity ivelcolnp[1] = " << (*ivelcolnp)[1] << std::endl;
    std::cout << "applying interface velocity ivelcolnp[2] = " << (*ivelcolnp)[2] << std::endl;
    std::ofstream f;
    if (step_ <= 1)
      f.open("outifacevelnp.txt",std::fstream::trunc);
    else
      f.open("outifacevelnp.txt",std::fstream::ate | std::fstream::app);

    f << time_ << " " << (*ivelcolnp)[0] << "  " << endl;

    f.close();
  }
  
  if (myrank_ == 0 && ivelcoln->MyLength() >= 3)
  {
    std::cout << "applying interface velocity ivelcoln[0] = " << (*ivelcoln)[0] << std::endl;
    std::cout << "applying interface velocity ivelcoln[1] = " << (*ivelcoln)[1] << std::endl;
    std::cout << "applying interface velocity ivelcoln[2] = " << (*ivelcoln)[2] << std::endl;
    std::ofstream f;
    if (step_ <= 1)
      f.open("outifaceveln.txt",std::fstream::trunc);
    else
      f.open("outifaceveln.txt",std::fstream::ate | std::fstream::app);

    f << time_ << " " << (*ivelcoln)[0] << "  " << endl;

    f.close();
  }

  if (myrank_ == 0 && ivelcolnp->MyLength() >= 3)
  {
    std::ofstream f;
    if (step_ <= 1)
      f.open("outifaceanalytischvel.txt",std::fstream::trunc);
    else
      f.open("outifaceanalytischvel.txt",std::fstream::ate | std::fstream::app);

    const double periodendauer = 10.0;
    f << time_ << " " << (-1.5*std::sin(2.0*time_* PI/periodendauer) * PI/periodendauer) << endl;

    f.close();
  }



  if (myrank_ == 0)
  {
    // action for elements
    if (timealgo_==timeint_stationary)
    {

    }
    else
    {
      cout << "******************************************************" << endl;
      cout << "* Warning! Does not work for moving boundaries, yet! *" << endl;
      cout << "******************************************************" << endl;
    }

    printf("+------------+-------------------+--------------+--------------+--------------+--------------+--------------+--------------+\n");
    printf("|- step/max -|- tol      [norm] -|-- vel-res ---|-- pre-res ---|-- fullres ---|-- vel-inc ---|-- pre-inc ---|-- fullinc ---|\n");
  }

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
      {
        eleparams.set("action","calc_fluid_systemmat_and_residual");
      }

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
      eleparams.set("interface velocity",ivelcolnp);
      //cout << "interface velocity" << endl;
      //cout << *ivelcol << endl;

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

        // How to extract the density from the fluid material?
        //trueresidual_->Update(density_/dta_/theta_,*residual_,0.0);
        iforcecolnp->Scale(density_/dta_/theta_);

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
    {
      Epetra_Vector residual(*residual_);
      residual_->Multiply(1.0,*invtoggle_,residual,0.0);
    }

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
          printf(")\n");
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

    //--------- Apply dirichlet boundary conditions to system of equations
    //          residual discplacements are supposed to be zero at
    //          boundary conditions
    incvel_->PutScalar(0.0);
    {
      // time measurement: application of dbc
      TEUCHOS_FUNC_TIME_MONITOR("      + apply DBC");
      LINALG::ApplyDirichlettoSystem(sysmat_,incvel_,residual_,zeros_,dirichtoggle_);
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
  
  cutterdiscret->SetState("iforcenp", iforcecolnp);


  const int nsd = 3;
  const Epetra_Map* dofcolmap = cutterdiscret->DofColMap();
  BlitzVec3 c;
  c = 0.0;
  for (int inode = 0; inode < cutterdiscret->NumMyColNodes(); ++inode)
  {
    const DRT::Node* node = cutterdiscret->lColNode(inode);
    const std::vector<int> dof = cutterdiscret->Dof(node);
    for (int isd = 0; isd < nsd; ++isd)
    {
      const double val = (*iforcecolnp)[dofcolmap->LID(dof[isd])];
      c(isd) -= val; // minus to get correct sign of lift and drag (force acting on the body)
    }

  }
  
  {
    std::stringstream s;
    std::stringstream header;
    
    header << left  << std::setw(10) << "Time" 
           << right << std::setw(16) << "F_x"
           << right << std::setw(16) << "F_y"
           << right << std::setw(16) << "F_z";
    s << left  << std::setw(10) << scientific << time_ 
      << right << std::setw(16) << scientific << c(0)
      << right << std::setw(16) << scientific << c(1)
      << right << std::setw(16) << scientific << c(2);
    
    std::ofstream f;
    if (step_ <= 1)
    {
      f.open("liftdrag.txt",std::fstream::trunc);
      //f << header.str() << endl;
    }
    else
    {
      f.open("liftdrag.txt",std::fstream::ate | std::fstream::app);
    }
    f << s.str() << endl;
    f.close();
    
    //cout << header.str() << endl << s.str() << endl;
  }
  
  if (myrank_ == 0 && iforcecolnp->MyLength() >= 3)
  {
    std::ofstream f;
    if (step_ <= 1)
      f.open("outifaceforce.txt",std::fstream::trunc);
    else
      f.open("outifaceforce.txt",std::fstream::ate | std::fstream::app);

    f << time_ << " " << (*iforcecolnp)[0] << "  " << endl;

    f.close();
  }
  

} // FluidImplicitTimeInt::NonlinearSolve



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::XFluidImplicitTimeInt::Evaluate(Teuchos::RCP<const Epetra_Vector> vel)
{
  sysmat_->Zero();
  dserror("no monolithic FSI tested, check first!");

  // set the new solution we just got
  if (vel!=Teuchos::null)
  {
    const int len = vel->MyLength();

    // Take Dirichlet values from velnp and add vel to veln for non-Dirichlet
    // values.
    //
    // There is no epetra operation for this! Maybe we could have such a beast
    // in ANA?

    double* veln  = &(*state_.veln_)[0];
    double* velnp = &(*state_.velnp_)[0];
    double* dt    = &(*dirichtoggle_)[0];
    double* idv   = &(*invtoggle_)[0];
    const double* incvel = &(*vel)[0];

    //------------------------------------------------ update (u,p) trial
    for (int i=0; i<len; ++i)
    {
      velnp[i] = velnp[i]*dt[i] + (veln[i] + incvel[i])*idv[i];
    }
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
  eleparams.set("include reactive terms for linearisation",params_.get<bool>("Use reaction terms for linearisation",false));

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

  trueresidual_->Update(density_/dta_/theta_,*residual_,0.0);

  // Apply dirichlet boundary conditions to system of equations
  // residual displacements are supposed to be zero at boundary
  // conditions
  incvel_->PutScalar(0.0);
  LINALG::ApplyDirichlettoSystem(sysmat_,incvel_,residual_,zeros_,dirichtoggle_);
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::XFluidImplicitTimeInt::TimeUpdate()
{

  // prev. acceleration becomes (n-1)-accel. of next time step
  const Teuchos::RCP<Epetra_Vector> accn_tmp = rcp(new Epetra_Vector(*state_.accn_));
  
  // compute acceleration 
  // note a(n+1) is directly stored in a(n),
  // hence we use a(n-1) as a(n) to save a vector copy (see line above)
  TIMEINT_THETA_BDF2::CalculateAcceleration(
      state_.velnp_, state_.veln_, state_.velnm_, accn_tmp,
          timealgo_, step_, theta_, dta_, dtp_,
          state_.accn_);

  // solution of this step becomes most recent solution of the last step
  state_.velnm_->Update(1.0,*state_.veln_ ,0.0);
  state_.veln_ ->Update(1.0,*state_.velnp_,0.0);

  return;
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::XFluidImplicitTimeInt::Output()
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

    // output real pressure
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

    if (step_%uprestart_ == 0) //add restart data
    {
      output_.WriteVector("accn", state_.accn_);
      output_.WriteVector("veln", state_.veln_);
      output_.WriteVector("velnm", state_.velnm_);
    }
  }

  // write restart also when uprestart_ is not a integer multiple of upres_
  else if (step_%uprestart_ == 0)
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


  OutputToGmsh();

  return;
} // FluidImplicitTimeInt::Output


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::XFluidImplicitTimeInt::ReadRestart(int step)
{
  dserror("check wich data was written. one might need 2 discretization writers: \n \
           one for the output and one for the restart with changing vectors.\n \
           Problem is, the numdofs are written during WriteMesh(). is that used for restart ore not?");

  IO::DiscretizationReader reader(discret_,step);
  time_ = reader.ReadDouble("time");
  step_ = reader.ReadInt("step");

  reader.ReadVector(state_.velnp_,"velnp");
  reader.ReadVector(state_.veln_, "veln");
  reader.ReadVector(state_.velnm_,"velnm");
  reader.ReadVector(state_.accn_ ,"accn");

}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::XFluidImplicitTimeInt::OutputToGmsh()
{
  const Teuchos::ParameterList& xfemparams = DRT::Problem::Instance()->XFEMGeneralParams();
  const bool gmshdebugout = (bool)getIntegralValue<int>(xfemparams,"GMSH_DEBUG_OUT");

  if (gmshdebugout)
  {
    cout << "XFluidImplicitTimeInt::OutputToGmsh()" << endl;

    std::stringstream filename;
    std::stringstream filenamedel;
    filename << allfiles.outputfile_kenner << "_solution_pressure_" << std::setw(5) << setfill('0') << step_ << ".pos";
    filenamedel << allfiles.outputfile_kenner << "_solution_pressure_" << std::setw(5) << setfill('0') << step_-5 << ".pos";
    std::remove(filenamedel.str().c_str());
    std::cout << "writing " << left << std::setw(50) <<filename.str()<<"...";
    std::ofstream f_system(filename.str().c_str());

    const XFEM::PHYSICS::Field field = XFEM::PHYSICS::Pres;

    {
      stringstream gmshfilecontent;
      gmshfilecontent << "View \" " << "Pressure Solution (Physical) \" {"
      << endl;
      for (int i=0; i<discret_->NumMyColElements(); ++i)
      {

        const DRT::Element* actele = discret_->lColElement(i);
        const int elegid = actele->Id();

        BlitzMat xyze_xfemElement(DRT::UTILS::InitialPositionArrayBlitz(actele));

        const map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> element_ansatz(XFLUID::getElementAnsatz(actele->Shape()));

        // create local copy of information about dofs
        const XFEM::ElementDofManager eledofman =
          dofmanagerForOutput_->constructElementDofManager(*actele,
              element_ansatz);

        vector<int> lm;
        vector<int> lmowner;
        actele->LocationVector(*(discret_), lm, lmowner);

        // extract local values from the global vector
        vector<double> myvelnp(lm.size());
        DRT::UTILS::ExtractMyValues(*state_.velnp_, myvelnp, lm);

        const int numparam = eledofman.NumDofPerField(field);
        const vector<int>& dofpos = eledofman.LocalDofPosPerField(field);

        BlitzVec elementvalues(numparam);
        for (int iparam=0; iparam<numparam; ++iparam)
          elementvalues(iparam) = myvelnp[dofpos[iparam]];

        const GEO::DomainIntCells& domainintcells =
          dofmanagerForOutput_->getInterfaceHandle()->GetDomainIntCells(elegid, actele->Shape());
        for (GEO::DomainIntCells::const_iterator cell =
          domainintcells.begin(); cell != domainintcells.end(); ++cell)
        {
          BlitzVec cellvalues(DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
          XFEM::computeScalarCellNodeValuesFromNodalUnknowns(*actele, dofmanagerForOutput_->getInterfaceHandle(), eledofman,
              *cell, field, elementvalues, cellvalues);
          BlitzMat xyze_cell(3, cell->NumNode());
          cell->NodalPosXYZ(*actele, xyze_cell);
          gmshfilecontent << IO::GMSH::cellWithScalarFieldToString(
              cell->Shape(), cellvalues, xyze_cell) << endl;
        }
        if (elegid == 1 and elementvalues.size() > 0)
        {
          //std::cout << elementvalues << std::endl;
          std::ofstream f;
          if (step_ <= 1)
            f.open("outflowpres.txt",std::fstream::trunc);
          else
            f.open("outflowpres.txt",std::fstream::ate | std::fstream::app);

          //f << time_ << " " << (-1.5*std::sin(0.1*2.0*time_* PI) * PI*0.1) << "  " << elementvalues(0,0) << endl;
          f << time_ << "  " << elementvalues(0) << endl;

          f.close();
        }
      }
      gmshfilecontent << "};" << endl;
      f_system << gmshfilecontent.str();
    }
    f_system.close();
    std::cout << " done" << endl;
  }
#if 0
  if (gmshdebugout)
  {
    std::stringstream filename;
    filename << allfiles.outputfile_kenner << "_solution_pressure_disc_" << std::setw(5) << setfill('0') << step_
    << ".pos";
    std::cout << "writing " << std::left << std::setw(50) <<filename.str()<<"...";
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

        BlitzMat xyze_xfemElement(DRT::UTILS::InitialPositionArrayBlitz(actele));

        const map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> element_ansatz(XFLUID::getElementAnsatz(actele->Shape()));

        // create local copy of information about dofs
        const XFEM::ElementDofManager eledofman =
          dofmanagerForOutput_->constructElementDofManager(*actele,
              element_ansatz);

        vector<int> lm;
        vector<int> lmowner;
        actele->LocationVector(*(discret_), lm, lmowner);

        // extract local values from the global vector
        vector<double> myvelnp(lm.size());
        DRT::UTILS::ExtractMyValues(*state_.velnp_, myvelnp, lm);

        const int numparam = eledofman.NumDofPerField(field);
        const vector<int>& dofpos = eledofman.LocalDofPosPerField(field);

        BlitzVec elementvalues(numparam);
        for (int iparam=0; iparam<numparam; ++iparam)
          elementvalues(iparam) = myvelnp[dofpos[iparam]];
        if(elementvalues.size() != 0)
        {
          //cout << "eleval DiscPres" << endl;
          //cout << elementvalues << endl;
        }
        const XFEM::DomainIntCells& domainintcells =
          dofmanagerForOutput_->getInterfaceHandle()->GetDomainIntCells(elegid, actele->Shape());
        for (XFEM::DomainIntCells::const_iterator cell =
          domainintcells.begin(); cell != domainintcells.end(); ++cell)
        {
          BlitzVec cellvalues(DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
          XFEM::computeScalarCellNodeValuesFromElementUnknowns(*actele, dofmanagerForOutput_->getInterfaceHandle(), eledofman,
              *cell, field, elementvalues, cellvalues);
          BlitzMat xyze_cell(3, cell->NumNode());
          cell->NodalPosXYZ(*actele, xyze_cell);
          if(elementvalues.size() != 0)
          {
            //cout << "cellvalues DiscPres" << endl;
            //cout << cellvalues << endl;
          }
          gmshfilecontent << IO::GMSH::cellWithScalarFieldToString(
              cell->Shape(), cellvalues, xyze_cell) << endl;
        }
      }
      gmshfilecontent << "};" << endl;
      f_system << gmshfilecontent.str();
    }
    f_system.close();
    std::cout << " done" << endl;
  }
#endif
#if 1
  if (gmshdebugout)
  {
    //std::stringstream filename;
    std::stringstream filenamexx;
    std::stringstream filenameyy;
    std::stringstream filenamezz;
    std::stringstream filenamexy;
    std::stringstream filenamexz;
    std::stringstream filenameyz;
    //filename   << "solution_tau_disc_"   << std::setw(5) << setfill('0') << step_ << ".pos";
    filenamexx << allfiles.outputfile_kenner << "_solution_tauxx_disc_" << std::setw(5) << setfill('0') << step_ << ".pos";
    filenameyy << allfiles.outputfile_kenner << "_solution_tauyy_disc_" << std::setw(5) << setfill('0') << step_ << ".pos";
    filenamezz << allfiles.outputfile_kenner << "_solution_tauzz_disc_" << std::setw(5) << setfill('0') << step_ << ".pos";
    filenamexy << allfiles.outputfile_kenner << "_solution_tauxy_disc_" << std::setw(5) << setfill('0') << step_ << ".pos";
    filenamexz << allfiles.outputfile_kenner << "_solution_tauxz_disc_" << std::setw(5) << setfill('0') << step_ << ".pos";
    filenameyz << allfiles.outputfile_kenner << "_solution_tauyz_disc_" << std::setw(5) << setfill('0') << step_ << ".pos";
    std::cout << "writing " << std::left << std::setw(50) <<"stresses"<<"...";
    flush(cout);
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
      gmshfilecontentxx << "View \" " << "Discontinous Viscous Stress (xx) Solution (Physical) \" {" << endl;
      gmshfilecontentyy << "View \" " << "Discontinous Viscous Stress (yy) Solution (Physical) \" {" << endl;
      gmshfilecontentzz << "View \" " << "Discontinous Viscous Stress (zz) Solution (Physical) \" {" << endl;
      gmshfilecontentxy << "View \" " << "Discontinous Viscous Stress (xy) Solution (Physical) \" {" << endl;
      gmshfilecontentxz << "View \" " << "Discontinous Viscous Stress (xz) Solution (Physical) \" {" << endl;
      gmshfilecontentyz << "View \" " << "Discontinous Viscous Stress (yz) Solution (Physical) \" {" << endl;
      for (int i=0; i<discret_->NumMyColElements(); ++i)
      {

        const DRT::Element* actele = discret_->lColElement(i);
        const int elegid = actele->Id();

        BlitzMat xyze_xfemElement(DRT::UTILS::InitialPositionArrayBlitz(actele));

        const map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> element_ansatz(XFLUID::getElementAnsatz(actele->Shape()));

        // create local copy of information about dofs
        const XFEM::ElementDofManager eledofman =
          dofmanagerForOutput_->constructElementDofManager(*actele,
              element_ansatz);

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

        BlitzMat elementvalues(9,numparam);
        for (int iparam=0; iparam<numparam; ++iparam) elementvalues(0,iparam) = myvelnp[dofposxx[iparam]];
        for (int iparam=0; iparam<numparam; ++iparam) elementvalues(1,iparam) = myvelnp[dofposxy[iparam]];
        for (int iparam=0; iparam<numparam; ++iparam) elementvalues(2,iparam) = myvelnp[dofposxz[iparam]];
        for (int iparam=0; iparam<numparam; ++iparam) elementvalues(3,iparam) = myvelnp[dofposxy[iparam]];
        for (int iparam=0; iparam<numparam; ++iparam) elementvalues(4,iparam) = myvelnp[dofposyy[iparam]];
        for (int iparam=0; iparam<numparam; ++iparam) elementvalues(5,iparam) = myvelnp[dofposyz[iparam]];
        for (int iparam=0; iparam<numparam; ++iparam) elementvalues(6,iparam) = myvelnp[dofposxz[iparam]];
        for (int iparam=0; iparam<numparam; ++iparam) elementvalues(7,iparam) = myvelnp[dofposyz[iparam]];
        for (int iparam=0; iparam<numparam; ++iparam) elementvalues(8,iparam) = myvelnp[dofposzz[iparam]];

        BlitzVec elementvaluexx(numparam); for (int iparam=0; iparam<numparam; ++iparam) elementvaluexx(iparam) = myvelnp[dofposxx[iparam]];
        BlitzVec elementvalueyy(numparam); for (int iparam=0; iparam<numparam; ++iparam) elementvalueyy(iparam) = myvelnp[dofposyy[iparam]];
        BlitzVec elementvaluezz(numparam); for (int iparam=0; iparam<numparam; ++iparam) elementvaluezz(iparam) = myvelnp[dofposzz[iparam]];
        BlitzVec elementvaluexy(numparam); for (int iparam=0; iparam<numparam; ++iparam) elementvaluexy(iparam) = myvelnp[dofposxy[iparam]];
        BlitzVec elementvaluexz(numparam); for (int iparam=0; iparam<numparam; ++iparam) elementvaluexz(iparam) = myvelnp[dofposxz[iparam]];
        BlitzVec elementvalueyz(numparam); for (int iparam=0; iparam<numparam; ++iparam) elementvalueyz(iparam) = myvelnp[dofposyz[iparam]];


        const GEO::DomainIntCells& domainintcells =
          dofmanagerForOutput_->getInterfaceHandle()->GetDomainIntCells(elegid, actele->Shape());
        for (GEO::DomainIntCells::const_iterator cell =
          domainintcells.begin(); cell != domainintcells.end(); ++cell)
        {
          BlitzMat xyze_cell(3, cell->NumNode());
          cell->NodalPosXYZ(*actele, xyze_cell);

//          {
//          BlitzMat cellvalues(9,DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
//          XFEM::computeTensorCellNodeValuesFromElementUnknowns(*actele, dofmanagerForOutput_->getInterfaceHandle(), eledofman,
//              *cell, field, elementvalues, cellvalues);
//          gmshfilecontent << IO::GMSH::cellWithTensorFieldToString(cell->Shape(), cellvalues, xyze_cell) << endl;
//          }

          {
          BlitzVec cellvaluexx(DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
          XFEM::computeScalarCellNodeValuesFromElementUnknowns(*actele, dofmanagerForOutput_->getInterfaceHandle(), eledofman, *cell, field, elementvaluexx, cellvaluexx);
          gmshfilecontentxx << IO::GMSH::cellWithScalarFieldToString(cell->Shape(), cellvaluexx, xyze_cell) << endl;
          }
          {
          BlitzVec cellvalueyy(DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
          XFEM::computeScalarCellNodeValuesFromElementUnknowns(*actele, dofmanagerForOutput_->getInterfaceHandle(), eledofman, *cell, field, elementvalueyy, cellvalueyy);
          gmshfilecontentyy << IO::GMSH::cellWithScalarFieldToString(cell->Shape(), cellvalueyy, xyze_cell) << endl;
          }
          {
          BlitzVec cellvaluezz(DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
          XFEM::computeScalarCellNodeValuesFromElementUnknowns(*actele, dofmanagerForOutput_->getInterfaceHandle(), eledofman, *cell, field, elementvaluezz, cellvaluezz);
          gmshfilecontentzz << IO::GMSH::cellWithScalarFieldToString(cell->Shape(), cellvaluezz, xyze_cell) << endl;
          }
          {
          BlitzVec cellvaluexy(DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
          XFEM::computeScalarCellNodeValuesFromElementUnknowns(*actele, dofmanagerForOutput_->getInterfaceHandle(), eledofman, *cell, field, elementvaluexy, cellvaluexy);
          gmshfilecontentxy << IO::GMSH::cellWithScalarFieldToString(cell->Shape(), cellvaluexy, xyze_cell) << endl;
          }
          {
          BlitzVec cellvaluexz(DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
          XFEM::computeScalarCellNodeValuesFromElementUnknowns(*actele, dofmanagerForOutput_->getInterfaceHandle(), eledofman, *cell, field, elementvaluexz, cellvaluexz);
          gmshfilecontentxz << IO::GMSH::cellWithScalarFieldToString(cell->Shape(), cellvaluexz, xyze_cell) << endl;
          }
          {
          BlitzVec cellvalueyz(DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
          XFEM::computeScalarCellNodeValuesFromElementUnknowns(*actele, dofmanagerForOutput_->getInterfaceHandle(), eledofman, *cell, field, elementvalueyz, cellvalueyz);
          gmshfilecontentyz << IO::GMSH::cellWithScalarFieldToString(cell->Shape(), cellvalueyz, xyze_cell) << endl;
          }
        }
      }
      //gmshfilecontent   << "};" << endl;
      gmshfilecontentxx << "};" << endl;
      gmshfilecontentyy << "};" << endl;
      gmshfilecontentzz << "};" << endl;
      gmshfilecontentxy << "};" << endl;
      gmshfilecontentxz << "};" << endl;
      gmshfilecontentyz << "};" << endl;
      //f_system   << gmshfilecontent.str();
      f_systemxx << gmshfilecontentxx.str();
      f_systemyy << gmshfilecontentyy.str();
      f_systemzz << gmshfilecontentzz.str();
      f_systemxy << gmshfilecontentxy.str();
      f_systemxz << gmshfilecontentxz.str();
      f_systemyz << gmshfilecontentyz.str();
    }
    std::cout << " done" << endl;
  }
#endif


  PlotVectorFieldToGmsh(state_.velnp_, "_solution_velocity_","Velocity Solution (Physical) n+1",true);
  PlotVectorFieldToGmsh(state_.veln_,  "_solution_velocity_old_step_","Velocity Solution (Physical) n",false);
  PlotVectorFieldToGmsh(state_.velnm_, "_solution_velocity_old2_step_","Velocity Solution (Physical) n-1",false);
  PlotVectorFieldToGmsh(state_.accn_,  "_solution_acceleration_old_step_","Acceleration Solution (Physical) n",false);
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::XFluidImplicitTimeInt::PlotVectorFieldToGmsh(
    const Teuchos::RCP<Epetra_Vector>   vectorfield,
    const std::string filestr,
    const std::string name_in_gmsh,
    const bool plot_to_gnuplot
    ) const
{

  const Teuchos::ParameterList& xfemparams = DRT::Problem::Instance()->XFEMGeneralParams();
  const bool gmshdebugout = (bool)getIntegralValue<int>(xfemparams,"GMSH_DEBUG_OUT");

  if (gmshdebugout)
  {

    bool ele_to_textfile = false;
    std::stringstream filename;
    std::stringstream filenamedel;
    filename << allfiles.outputfile_kenner << filestr << std::setw(5) << std::setfill('0') << step_ << ".pos";
    filenamedel << allfiles.outputfile_kenner << filestr << std::setw(5) << std::setfill('0') << step_-5 << ".pos";
    std::remove(filenamedel.str().c_str());
    std::cout << "writing " << std::left << std::setw(50) <<filename.str()<<"...";
    std::ofstream f_system(filename.str().c_str());

    {
      stringstream gmshfilecontent;
      gmshfilecontent << "View \" " << name_in_gmsh << "\" {" << endl;
      for (int i=0; i<discret_->NumMyColElements(); ++i)
      {
        const DRT::Element* actele = discret_->lColElement(i);
        const int elegid = actele->Id();

        BlitzMat xyze_xfemElement(DRT::UTILS::InitialPositionArrayBlitz(actele));

        const map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> element_ansatz(XFLUID::getElementAnsatz(actele->Shape()));

        // create local copy of information about dofs
        const XFEM::ElementDofManager eledofman =
          dofmanagerForOutput_->constructElementDofManager(*actele,
              element_ansatz);

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
        BlitzMat elementvalues(3, numparam);
        for (int iparam=0; iparam<numparam; ++iparam)
        {
          elementvalues(0, iparam) = myvelnp[dofposvelx[iparam]];
          elementvalues(1, iparam) = myvelnp[dofposvely[iparam]];
          elementvalues(2, iparam) = myvelnp[dofposvelz[iparam]];
        }

        if (!dofmanagerForOutput_->getInterfaceHandle()->ElementIntersected(elegid))
        {
          const GEO::DomainIntCells& domainintcells =
            dofmanagerForOutput_->getInterfaceHandle()->GetDomainIntCells(elegid, actele->Shape());
          for (GEO::DomainIntCells::const_iterator cell =
            domainintcells.begin(); cell != domainintcells.end(); ++cell)
          {
            BlitzMat cellvalues(3, DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
            //std::cout << cellvalues << endl;
            XFEM::computeVectorCellNodeValues(*actele, dofmanagerForOutput_->getInterfaceHandle(), eledofman,
                *cell, XFEM::PHYSICS::Velx, elementvalues, cellvalues);
            BlitzMat xyze_cell(3, cell->NumNode());
            cell->NodalPosXYZ(*actele, xyze_cell);
            gmshfilecontent << IO::GMSH::cellWithVectorFieldToString(
                cell->Shape(), cellvalues, xyze_cell) << endl;
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
            {
              BlitzMat cellvalues(3, DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
              //std::cout << cellvalues << endl;
              XFEM::computeVectorCellNodeValues(*actele, dofmanagerForOutput_->getInterfaceHandle(), eledofman,
                  *cell, XFEM::PHYSICS::Velx, elementvalues, cellvalues);
              BlitzMat xyze_cell(3, cell->NumNode());
              cell->NodalPosXYZ(*actele, xyze_cell);
              gmshfilecontent << IO::GMSH::cellWithVectorFieldToString(
                  cell->Shape(), cellvalues, xyze_cell) << endl;
            }
          }

          // draw uncutted element
          {
            BlitzMat elevalues(3, DRT::UTILS::getNumberOfElementNodes(actele->Shape()));
            const GEO::DomainIntCell cell(actele->Shape());
            elevalues = 0.0;
//            XFEM::computeVectorCellNodeValues(*actele, dofmanagerForOutput_->getInterfaceHandle(), eledofman,
//                              cell, XFEM::PHYSICS::Velx, elementvalues, elevalues);

            const BlitzMat xyze_ele(DRT::UTILS::InitialPositionArrayBlitz(actele));
            gmshfilecontent << IO::GMSH::cellWithVectorFieldToString(
                actele->Shape(), elevalues, xyze_ele) << endl;
          }

        }
        //if (dofmanagerForOutput_->getInterfaceHandle()->ElementIntersected(elegid) and not ele_to_textfile and ele_to_textfile2)
        if (elegid == 1 and elementvalues.size() > 0 and plot_to_gnuplot)
        {
          ele_to_textfile = true;
          //std::cout << elementvalues << std::endl;
          std::ofstream f;
          if (step_ <= 1)
            f.open("outflowvel.txt",std::fstream::trunc);
          else
            f.open("outflowvel.txt",std::fstream::ate | std::fstream::app);

          //f << time_ << " " << (-1.5*std::sin(0.1*2.0*time_* PI) * PI*0.1) << "  " << elementvalues(0,0) << endl;
          f << time_ << "  " << elementvalues(0,0) << endl;

          f.close();
        }

      }
      gmshfilecontent << "};" << endl;
      f_system << gmshfilecontent.str();
    }
    f_system.close();
    std::cout << " done" << endl;
  }
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::XFluidImplicitTimeInt::SetInitialFlowField(
    Teuchos::RCP<DRT::Discretization> cutterdiscret,
    int whichinitialfield,
    int startfuncno
    )
{
  // create zero displacement vector to use initial position of interface
  {
    const Epetra_Map* fluidsurface_dofcolmap = cutterdiscret->DofColMap();
    Teuchos::RCP<Epetra_Vector> idispcolnp  = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
    Teuchos::RCP<Epetra_Vector> ivelcolnp   = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
    
    Teuchos::RCP<Epetra_Vector> idispcoln   = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
    Teuchos::RCP<Epetra_Vector> ivelcoln    = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
    Teuchos::RCP<Epetra_Vector> iacccoln    = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
    
    cutterdiscret->SetState("idispcolnp",idispcolnp);
    cutterdiscret->SetState("ivelcolnp",ivelcolnp);
    
    cutterdiscret->SetState("idispcoln",idispcoln);
    cutterdiscret->SetState("ivelcoln",ivelcoln);
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

        double initialval=DRT::UTILS::FunctionManager::Instance().Funct(startfuncno-1).Evaluate(index,lnode->X());

        state_.velnp_->ReplaceGlobalValues(1,&initialval,&gid);
        state_.veln_ ->ReplaceGlobalValues(1,&initialval,&gid);
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
void FLD::XFluidImplicitTimeInt::SolveStationaryProblem(
    Teuchos::RCP<DRT::Discretization> cutterdiscret
    )
{

  const Epetra_Map* fluidsurface_dofcolmap = cutterdiscret->DofColMap();
  const Teuchos::RCP<Epetra_Vector> idispcolnp  = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
  const Teuchos::RCP<Epetra_Vector> idispcoln   = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
  const Teuchos::RCP<Epetra_Vector> ivelcolnp   = LINALG::CreateVector(*fluidsurface_dofcolmap,true); // one could give a velocity here to have stationary flow over the interface
  const Teuchos::RCP<Epetra_Vector> ivelcoln    = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
  const Teuchos::RCP<Epetra_Vector> iacccoln    = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
  const Teuchos::RCP<Epetra_Vector> itruerescol = LINALG::CreateVector(*fluidsurface_dofcolmap,true);

  cutterdiscret->SetState("idispcolnp", idispcolnp);
  cutterdiscret->SetState("idispcoln", idispcoln);
  cutterdiscret->SetState("ivelcolnp",ivelcolnp);
  cutterdiscret->SetState("ivelcoln", ivelcoln);
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
     // dirichtoggle is 1 for dirichlet dofs, 0 elsewhere
     discret_->EvaluateDirichlet(eleparams,state_.velnp_,null,null,dirichtoggle_);
     discret_->ClearState();

     // evaluate Neumann b.c.
     neumann_loads_->PutScalar(0.0);
     discret_->EvaluateNeumann(eleparams,*neumann_loads_);
     discret_->ClearState();
   }

   // compute an inverse of the dirichtoggle vector
   invtoggle_->PutScalar(1.0);
   invtoggle_->Update(-1.0,*dirichtoggle_,1.0);


    // -------------------------------------------------------------------
    //                     solve nonlinear equation system
    // -------------------------------------------------------------------
    NonlinearSolve(cutterdiscret);

    // -------------------------------------------------------------------
    //                    calculate lift'n'drag forces
    // -------------------------------------------------------------------
    int liftdrag = params_.get<int>("liftdrag");

    if(liftdrag == 0); // do nothing, we don't want lift & drag
    if(liftdrag == 1)
      dserror("how did you manage to get here???");
    if(liftdrag == 2)
      LiftDrag();

    // -------------------------------------------------------------------
    //                         output of solution
    // -------------------------------------------------------------------
    Output();

  } // end of time loop

} // FluidImplicitTimeInt::SolveStationaryProblem


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



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*!
\brief calculate lift&drag forces and angular moments

Lift and drag forces are based upon the right hand side true-residual entities
of the corresponding nodes. The contribution of the end node of a line is entirely
added to a present L&D force.

Idea of this routine:

create

map< label, std::set<DRT::Node*> >

which is a set of nodes to each L&D Id
nodal forces of all the nodes within one set are added to one L&D force

Notice: Angular moments obtained from lift&drag forces currently refer to the
        initial configuration, i.e. are built with the coordinates X of a particular
        node irrespective of its current position.
*/
void FLD::XFluidImplicitTimeInt::LiftDrag() const
{
  std::map< const int, std::set<DRT::Node* > > ldnodemap;
  std::map< const int, const std::vector<double>* > ldcoordmap;

  // allocate and initialise LiftDrag conditions
  std::vector<DRT::Condition*> ldconds;
  discret_->GetCondition("LIFTDRAG",ldconds);

  // space dimension of the problem
  const int ndim = params_.get<int>("number of velocity degrees of freedom");

  // there is an L&D condition if it has a size
  if( ldconds.size() )
  {

    // prepare output
    if (myrank_==0)
    {
      std::cout << "Lift and drag calculation:" << "\n";
      if (ndim == 2)
      {
        std::cout << "lift'n'drag Id      F_x             F_y             M_z :" << "\n";
      }
      if (ndim == 3)
      {
        std::cout << "lift'n'drag Id      F_x             F_y             F_z           ";
        std::cout << "M_x             M_y             M_z :" << "\n";
      }
    }

    // sort data
    for( unsigned i=0; i<ldconds.size(); ++i) // loop L&D conditions (i.e. lines in .dat file)
    {
      /* get label of present LiftDrag condition  */
      const unsigned int label = ldconds[i]->Getint("label");
      /* get new nodeset for new label OR:
         return pointer to nodeset for known label ... */
      std::set<DRT::Node*>& nodes = ldnodemap[label];

      // centre coordinates to present label
      ldcoordmap[label] = ldconds[i]->Get<vector<double> >("centerCoord");

      /* get pointer to its nodal Ids*/
      const vector<int>* ids = ldconds[i]->Get<vector<int> >("Node Ids");

      /* put all nodes belonging to the L&D line or surface into 'nodes' which are
         associated with the present label */
      for (unsigned j=0; j<ids->size(); ++j)
      {
        // give me present node Id
        const int node_id = (*ids)[j];
        // put it into nodeset of actual label if node is new and mine
        if( discret_->HaveGlobalNode(node_id) && discret_->gNode(node_id)->Owner()==myrank_ )
	  nodes.insert(discret_->gNode(node_id));
      }
    } // end loop over conditions


    // now step the label map
    for( std::map< const int, std::set<DRT::Node*> >::const_iterator labelit = ldnodemap.begin();
         labelit != ldnodemap.end(); ++labelit )
    {
      const std::set<DRT::Node*>& nodes = labelit->second; // pointer to nodeset of present label
      const int label = labelit->first;                    // the present label
      std::vector<double> values(6,0.0);             // vector with lift&drag forces
      std::vector<double> resultvec(6,0.0);          // vector with lift&drag forces after communication

      // get also pointer to centre coordinates
      const std::vector<double>* centerCoord = ldcoordmap[label];

      // loop all nodes within my set
      for( std::set<DRT::Node*>::const_iterator actnode = nodes.begin(); actnode != nodes.end(); ++actnode)
      {
        const double* x = (*actnode)->X(); // pointer to nodal coordinates
        const Epetra_BlockMap& rowdofmap = trueresidual_->Map();
        const std::vector<int> dof = discret_->Dof(*actnode);

        std::vector<double> distances (3);
        for (unsigned j=0; j<3; ++j)
        {
          distances[j]= x[j]-(*centerCoord)[j];
        }
        // get nodal forces
        const double fx = (*trueresidual_)[rowdofmap.LID(dof[0])];
        const double fy = (*trueresidual_)[rowdofmap.LID(dof[1])];
        const double fz = (*trueresidual_)[rowdofmap.LID(dof[2])];
        values[0] += fx;
        values[1] += fy;
        values[2] += fz;

        // calculate nodal angular momenta
        values[3] += distances[1]*fz-distances[2]*fy;
        values[4] += distances[2]*fx-distances[0]*fz;
        values[5] += distances[0]*fy-distances[1]*fx;
      } // end: loop over nodes

      // care for the fact that we are (most likely) parallel
      trueresidual_->Comm().SumAll (&(values[0]), &(resultvec[0]), 6);

      // do the output
      if (myrank_==0)
      {
        if (ndim == 2)
	{
          std::cout << "     " << label << "         ";
          std::cout << std::scientific << resultvec[0] << "    ";
          std::cout << std::scientific << resultvec[1] << "    ";
          std::cout << std::scientific << resultvec[5];
          std::cout << "\n";
        }
        if (ndim == 3)
	{
          std::cout << "     " << label << "         ";
          std::cout << std::scientific << resultvec[0] << "    ";
          std::cout << std::scientific << resultvec[1] << "    ";
          std::cout << std::scientific << resultvec[2] << "    ";
          std::cout << std::scientific << resultvec[3] << "    ";
          std::cout << std::scientific << resultvec[4] << "    ";
          std::cout << std::scientific << resultvec[5];
          std::cout << "\n";
	}
      }

    } // end: loop over L&D labels
    if (myrank_== 0)
    {
      std::cout << "\n";
    }
  }
}//FluidImplicitTimeInt::LiftDrag



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
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


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::XFluidImplicitTimeInt::UseBlockMatrix(
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
  if (params_.get<bool>("mesh movement linearization"))
  {
    // allocate special mesh moving matrix
    mat = Teuchos::rcp(new LINALG::BlockSparseMatrix<FLD::UTILS::InterfaceSplitStrategy>(domainmaps,rangemaps,108,false,true));
    mat->SetCondElements(condelements);
    //meshmovematrix_ = mat;
  }
}

#endif /* CCADISCRET       */
