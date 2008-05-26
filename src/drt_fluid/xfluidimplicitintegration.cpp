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
#include "../drt_xfem/interface.H"
#include "../drt_xfem/dof_management.H"
#include "../drt_xfem/dof_distribution_switcher.H"
#include "../drt_xfem/field_enriched.H"
#include "../drt_xfem/enrichment_utils.H"
#include "fluid_utils.H"
#include "../drt_f3/xfluid3_interpolation.H"


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Constructor (public)                                     gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
XFluidImplicitTimeInt::XFluidImplicitTimeInt(
		RCP<DRT::Discretization> actdis,
		LINALG::Solver&       solver,
		ParameterList&        params,
		IO::DiscretizationWriter& output,
		const bool alefluid) :
  // call constructor for "nontrivial" objects
  discret_(actdis),
  solver_ (solver),
  params_ (params),
  output_ (output),
  alefluid_(alefluid),
  time_(0.0),
  step_(0),
  extrapolationpredictor_(params.get("do explicit predictor",true)),
  restartstep_(0),
  uprestart_(params.get("write restart every", -1)),
  writestep_(0),
  upres_(params.get("write solution every", -1)),
  writestresses_(params.get<int>("write stresses", 0))
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
  timealgo_ = params_.get<FLUID_TIMEINTTYPE>("time int algo");
  dtp_ = dta_ = params_.get<double>("time step size");
  stepmax_  = params_.get<int>   ("max number timesteps");
  maxtime_  = params_.get<double>("total time");
  theta_    = params_.get<double>("theta");
  
  // create empty cutter discretization
  Teuchos::RCP<DRT::Discretization> emptyboundarydis_ = DRT::UTILS::CreateDiscretizationFromCondition(
      actdis, "schnackelzappel", "DummyBoundary", "BELE3", vector<string>(0));
  Teuchos::RCP<Epetra_Vector> tmpdisp = LINALG::CreateVector(*emptyboundarydis_->DofRowMap(),true);
  // intersection with empty cutter will result in a complete fluid domain with no holes or intersections
  Teuchos::RCP<XFEM::InterfaceHandle> ih = rcp(new XFEM::InterfaceHandle(discret_,emptyboundarydis_,*tmpdisp));
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
  FLUID_UTILS::SetupFluidSplit(*discret_,dofset_out_,3,velpressplitterForOutput_);
  
  nodalDofDistributionMap_.clear();
  elementalDofDistributionMap_.clear();

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
void XFluidImplicitTimeInt::Integrate(
        RCP<DRT::Discretization> cutterdiscret
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
void XFluidImplicitTimeInt::TimeLoop(
        RCP<DRT::Discretization> cutterdiscret
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
  Teuchos::RCP<Epetra_Vector> ivelcol     = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
  Teuchos::RCP<Epetra_Vector> idispcol    = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
  Teuchos::RCP<Epetra_Vector> itruerescol = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
  
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
      NonlinearSolve(cutterdiscret,ivelcol,idispcol,itruerescol);
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
    //
    // One-step-Theta: (step>1)
    //
    //  accn_  = (velnp_-veln_) / (Theta * dt) - (1/Theta -1) * accn_
    //  "(n+1)"
    //
    //  velnm_ =veln_
    //  veln_  =velnp_
    //
    // BDF2:           (step>1)
    //
    //               2*dt(n)+dt(n-1)		  dt(n)+dt(n-1)
    //  accn_   = --------------------- velnp_ - --------------- veln_
    //             dt(n)*[dt(n)+dt(n-1)]	  dt(n)*dt(n-1)
    //
    //                     dt(n)
    //           + ----------------------- velnm_
    //             dt(n-1)*[dt(n)+dt(n-1)]
    //
    //
    //  velnm_ =veln_
    //  veln_  =velnp_
    //
    // BDF2 and  One-step-Theta: (step==1)
    //
    // The given formulas are only valid from the second timestep. In the
    // first step, the acceleration is calculated simply by
    //
    //  accn_  = (velnp_-veln_) / (dt)
    //
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
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | setup the variables to do a new time step                 u.kue 06/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void XFluidImplicitTimeInt::PrepareTimeStep()
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
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | setup the variables to do a new nonlinear iteration      a.ger 064/08|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void XFluidImplicitTimeInt::PrepareNonlinearSolve()
{

  // -------------------------------------------------------------------
  // set part of the rhs vector belonging to the old timestep
  //
  //
  //         One-step-Theta:
  //
  //                 hist_ = veln_ + dta*(1-Theta)*accn_
  //
  //
  //         BDF2: for constant time step:
  //
  //                   hist_ = 4/3 veln_ - 1/3 velnm_
  //
  // -------------------------------------------------------------------
  TIMEINT_THETA_BDF2::SetOldPartOfRighthandside(
      state_.veln_, state_.velnm_, state_.accn_,
          timealgo_, dta_, theta_,
          hist_);

  // -------------------------------------------------------------------
  //                     do explicit predictor step
  //
  //                     +-                                      -+
  //                     | /     dta \          dta  veln_-velnm_ |
  // velnp_ =veln_ + dta | | 1 + --- | accn_ - ----- ------------ |
  //                     | \     dtp /          dtp     dtp       |
  //                     +-                                      -+
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
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | contains the nonlinear iteration loop                     gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void XFluidImplicitTimeInt::ComputeInterfaceAndSetDOFs(
        RCP<DRT::Discretization>  cutterdiscret,
        const Epetra_Vector&      idispcol
        )
{
  // within this routine, no parallel re-distribution is allowed to take place
  // before and after this function, it's ok to do that

  // calling this function multiple times always results in the same solution vectors

  // compute Intersection
  RCP<XFEM::InterfaceHandle> ih = rcp(new XFEM::InterfaceHandle(discret_,cutterdiscret,idispcol));
  ih->toGmsh(step_);
  ihForOutput_ = ih;

  // apply enrichments
  RCP<XFEM::DofManager> dofmanager = rcp(new XFEM::DofManager(ih));
  dofmanager->toGmsh(ih, step_);
  dofmanagerForOutput_ = dofmanager; // save to be able to plot Gmsh stuff in Output()

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

  // debug info: print ele dofmanager information
  std::stringstream filename;
  filename << "eledofman_check_" << std::setw(5) << setfill('0') << step_ << ".pos";
  std::cout << "writing '"<<filename.str()<<"'...";
  {
    std::ofstream f_system(filename.str().c_str());
      //f_system << IO::GMSH::disToString("Fluid", 0.0, ih->xfemdis(), ih->elementalDomainIntCells());
//      f_system << IO::GMSH::disToString("Solid", 1.0, ih->cutterdis());
//      {
//          // draw elements with associated gid
//          stringstream gmshfilecontent;
//          gmshfilecontent << "View \" " << "Element->Id() \" {" << endl;
//          for (int i=0; i<ih->xfemdis()->NumMyColElements(); ++i)
//          {
//              DRT::Element* actele = ih->xfemdis()->lColElement(i);
//              gmshfilecontent << IO::GMSH::elementToString(double(actele->Id()), actele) << endl;
//          };
//          gmshfilecontent << "};" << endl;
//          f_system << gmshfilecontent.str();
//      }
      {
          stringstream gmshfilecontent;
          gmshfilecontent << "View \" " << " NumDofPerElement() in element \" {" << endl;
          for (int i=0; i<ih->xfemdis()->NumMyColElements(); ++i)
          {
              DRT::Element* actele = ih->xfemdis()->lColElement(i);
              //const int ele_gid = actele->Id();
              //double val = 0.0;
              //std::map<int, const std::set<XFEM::FieldEnr> >::const_iterator blub = elementalDofs_.find(ele_gid);
              const double val = actele->NumDofPerElement();
              gmshfilecontent << IO::GMSH::elementToString(val, actele) << endl;

          };
          gmshfilecontent << "};" << endl;
          f_system << gmshfilecontent.str();
      }
//      {
//          stringstream gmshfilecontent;
//          gmshfilecontent << "View \" " << " eleDofManager->label in element \" {" << endl;
//          for (int i=0; i<ih->xfemdis()->NumMyColElements(); ++i)
//          {
//              const DRT::Element* actele = ihForOutput_->xfemdis()->lColElement(i);
//              const int numvirtualnodes = DRT::UTILS::getNumberOfElementNodes(actele->Shape());
//              
//              // create local copy of information about dofs
//              const XFEM::ElementDofManager eleDofManager = dofmanager->constructElementDofManager(*actele, numvirtualnodes);
//              
//              //const int ele_gid = actele->Id();
//              double val = 0.0;
//              
//              if (actele->NumDofPerElement() > 0)
//              {
//                const std::set<XFEM::FieldEnr> enrfieldset = eleDofManager.FieldEnrSetPerVirtualElementNode(0);
//                XFEM::FieldEnr firstenrfield = *(enrfieldset.begin());
//                const int label = firstenrfield.getEnrichment().XFEMConditionLabel();
//                val = label;
//              }
//              gmshfilecontent << IO::GMSH::elementToString(val, actele);
//
//          };
//          gmshfilecontent << "};" << endl;
//          f_system << gmshfilecontent.str();
//      }
      //f_system << IO::GMSH::getConfigString(2);
      f_system.close();
  }
  std::cout << "done" << endl;
  
  
  // store old (proc-overlapping) dofmap, compute new one and return it
  Epetra_Map olddofrowmap = *discret_->DofRowMap();
  discret_->FillComplete();
  Epetra_Map newdofrowmap = *discret_->DofRowMap();

  discret_->ComputeNullSpaceIfNecessary(solver_.Params());

  XFEM::NodalDofPosMap oldNodalDofDistributionMap(nodalDofDistributionMap_);
  XFEM::ElementalDofPosMap oldElementalDofDistributionMap(elementalDofDistributionMap_);
  dofmanager->fillDofDistributionMaps(nodalDofDistributionMap_,elementalDofDistributionMap_);

  std::cout << "switching " << endl;

  // create switcher
  const XFEM::DofDistributionSwitcher dofswitch(
          ih, dofmanager,
          olddofrowmap, newdofrowmap,
          oldNodalDofDistributionMap, nodalDofDistributionMap_,
          oldElementalDofDistributionMap, elementalDofDistributionMap_
          );

  // --------------------------------------------
  // switch state vectors to new dof distribution
  // --------------------------------------------

  // accelerations at time n and n-1
  dofswitch.mapVectorToNewDofDistribution(state_.accn_);
  dofswitch.mapVectorToNewDofDistribution(state_.accnm_);

  // velocities and pressures at time n+1, n and n-1
  dofswitch.mapVectorToNewDofDistribution(state_.velnp_);
  dofswitch.mapVectorToNewDofDistribution(state_.veln_);
  dofswitch.mapVectorToNewDofDistribution(state_.velnm_);

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
  hist_         = LINALG::CreateVector(newdofrowmap,true);

//  gridv_        = LINALG::CreateVector(newdofrowmap,true);
  
  dirichtoggle_ = LINALG::CreateVector(newdofrowmap,true);
  invtoggle_    = LINALG::CreateVector(newdofrowmap,true);

  zeros_        = LINALG::CreateVector(newdofrowmap,true);

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
  FLUID_UTILS::SetupXFluidSplit(*discret_,dofmanager,velpressplitter_);

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
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | contains the nonlinear iteration loop                     gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void XFluidImplicitTimeInt::NonlinearSolve(
        RCP<DRT::Discretization> cutterdiscret,
        RCP<Epetra_Vector>       idispcol,
        RCP<Epetra_Vector>       ivelcol,
        RCP<Epetra_Vector>       iforcecol
        )
{

  //idispcol->PutScalar(-0.3);
  //ivelcol->PutScalar(-0.5);
  
  ComputeInterfaceAndSetDOFs(cutterdiscret,*idispcol);

  PrepareNonlinearSolve();

  // time measurement: nonlinear iteration
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

  if (myrank_ == 0 && ivelcol->MyLength() >= 3)
  {
    std::cout << "applying interface velocity ivelcol[0] = " << (*ivelcol)[0] << std::endl;
    std::cout << "applying interface velocity ivelcol[1] = " << (*ivelcol)[1] << std::endl;
    std::cout << "applying interface velocity ivelcol[2] = " << (*ivelcol)[2] << std::endl;
    std::ofstream f;
    if (step_ <= 1)
      f.open("outifacevel.txt",std::fstream::trunc);
    else
      f.open("outifacevel.txt",std::fstream::ate | std::fstream::app);
    
    f << step_ << " " << (*ivelcol)[0] << "  " << endl;
    
    f.close();
    
    printf("+------------+-------------------+--------------+--------------+--------------+--------------+--------------+--------------+\n");
    printf("|- step/max -|- tol      [norm] -|-- vel-res ---|-- pre-res ---|-- fullres ---|-- vel-inc ---|-- pre-inc ---|-- fullinc ---|\n");
  }

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
      eleparams.set("total time",time_);
      eleparams.set("thsl",theta_*dta_);
      eleparams.set("dt",dta_);
      eleparams.set("include reactive terms for linearisation",params_.get<bool>("Use reaction terms for linearisation",false));

      // parameters for stabilization
      eleparams.sublist("STABILIZATION") = params_.sublist("STABILIZATION");

      // parameters for stabilization
      eleparams.sublist("TURBULENCE MODEL") = params_.sublist("TURBULENCE MODEL");

      // set vector values needed by elements
      discret_->ClearState();
      discret_->SetState("velnp",state_.velnp_);

      discret_->SetState("hist"  ,hist_ );
//      if (alefluid_)
//      {
//        discret_->SetState("dispnp", state_.dispnp_);
//        discret_->SetState("gridv", gridv_);
//      }
      // give interface velocity to elements
      eleparams.set("interface velocity",ivelcol);
      //cout << "interface velocity" << endl;
      //cout << *ivelcol << endl;
      // give interface force to elements
      iforcecol->PutScalar(0.0);
      eleparams.set("interface force",iforcecol);

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

        density_ = eleparams.get("density", 0.0);
        if (density_ <= 0.0) dserror("recieved illegal density value");

        // How to extract the density from the fluid material?
        trueresidual_->Update(density_/dta_/theta_,*residual_,0.0);
        iforcecol->Scale(density_/dta_/theta_);

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

    // free surface update
    if (alefluid_ and freesurface_->Relevant())
    {
      //using namespace LINALG::ANA;

      Teuchos::RefCountPtr<Epetra_Vector> fsvelnp = freesurface_->ExtractCondVector(state_.velnp_);
//      Teuchos::RefCountPtr<Epetra_Vector> fsdisp = freesurface_->ExtractCondVector(state_.dispn_);
      Teuchos::RefCountPtr<Epetra_Vector> fsdispnp = Teuchos::rcp(new Epetra_Vector(*freesurface_->CondMap()));

      // this is first order
      //*fsdispnp = *fsdisp + dta_*(*fsvelnp);
//      fsdispnp->Update(1.0,*fsdisp,dta_,*fsvelnp,0.0);

//      freesurface_->InsertCondVector(fsdispnp,state_.dispnp_);
//      freesurface_->InsertCondVector(fsvelnp,gridv_);
    }
    //cout << "*iforcecol" << endl;
    //cout << *iforcecol << endl;
  }
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
void XFluidImplicitTimeInt::Evaluate(Teuchos::RCP<const Epetra_Vector> vel)
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
  discret_->SetState("hist" ,hist_ );
//  if (alefluid_)
//  {
//    discret_->SetState("dispnp", state_.dispnp_);
//    discret_->SetState("gridv", gridv_);
//  }

  // call loop over elements
  discret_->Evaluate(eleparams,sysmat_,residual_);
  discret_->ClearState();

  density_ = eleparams.get("density", 0.0);

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
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                           gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void XFluidImplicitTimeInt::TimeUpdate()
{

  // update acceleration
  TIMEINT_THETA_BDF2::CalculateAcceleration(
      state_.velnp_, state_.veln_, state_.velnm_,
          timealgo_, step_, theta_, dta_, dtp_,
          state_.accn_, state_.accnm_);

  // solution of this step becomes most recent solution of the last step
  state_.velnm_->Update(1.0,*state_.veln_ ,0.0);
  state_.veln_ ->Update(1.0,*state_.velnp_,0.0);

//  if (alefluid_)
//  {
//    state_.dispnm_->Update(1.0,*state_.dispn_,0.0);
//    state_.dispn_ ->Update(1.0,*state_.dispnp_,0.0);
//  }

  return;
}// FluidImplicitTimeInt::TimeUpdate

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | output of solution vector to binio                        gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void XFluidImplicitTimeInt::Output()
{

  //-------------------------------------------- output of solution

  //increase counters
  restartstep_ += 1;
  writestep_ += 1;

  if (writestep_ == upres_)  //write solution
  {
    writestep_= 0;

    output_.NewStep    (step_,time_);
    //output_.WriteVector("velnp", state_.velnp_);
    std::set<XFEM::PHYSICS::Field> fields_out;
    fields_out.insert(XFEM::PHYSICS::Velx);
    fields_out.insert(XFEM::PHYSICS::Vely);
    fields_out.insert(XFEM::PHYSICS::Velz);
    fields_out.insert(XFEM::PHYSICS::Pres);
    Teuchos::RCP<Epetra_Vector> velnp_out = dofmanagerForOutput_->fillPhysicalOutputVector(
        *ihForOutput_,*state_.velnp_, dofset_out_, nodalDofDistributionMap_, fields_out);
    //cout << *velnp_out << endl;
    output_.WriteVector("velnp", velnp_out);
    
    // output real pressure
    Teuchos::RCP<Epetra_Vector> pressure = velpressplitterForOutput_.ExtractCondVector(velnp_out);
    pressure->Scale(density_);
    output_.WriteVector("pressure", pressure);

    //output_.WriteVector("residual", trueresidual_);
//    if (alefluid_)
//      output_.WriteVector("dispnp", state_.dispnp_);

    //only perform stress calculation when output is needed
    if (writestresses_)
    {
      dserror("not supported, yet");
      RCP<Epetra_Vector> traction = CalcStresses();
      output_.WriteVector("traction",traction);
    }

    // write domain decomposition for visualization (only once!)
    if (step_==upres_)
     output_.WriteElementData();

    if (restartstep_ == uprestart_) //add restart data
    {
      restartstep_ = 0;

      output_.WriteVector("accn", state_.accn_);
      output_.WriteVector("veln", state_.veln_);
      output_.WriteVector("velnm", state_.velnm_);

//      if (alefluid_)
//      {
//        output_.WriteVector("dispn", state_.dispn_);
//        output_.WriteVector("dispnm",state_.dispnm_);
//      }
    }
  }
    
   
  // write restart also when uprestart_ is not a integer multiple of upres_
  if ((restartstep_ == uprestart_) && (writestep_ > 0))
  {
    restartstep_ = 0;

    output_.NewStep    (step_,time_);
    output_.WriteVector("velnp", state_.velnp_);
    //output_.WriteVector("residual", trueresidual_);
//    if (alefluid_)
//    {
//      output_.WriteVector("dispnp", state_.dispnp_);
//      output_.WriteVector("dispn", state_.dispn_);
//      output_.WriteVector("dispnm",state_.dispnm_);
//    }

    //only perform stress calculation when output is needed
    if (writestresses_)
    {
      RCP<Epetra_Vector> traction = CalcStresses();
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
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | output of solution vector to gmsh                          acki 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void XFluidImplicitTimeInt::OutputToGmsh()
{
    {
      std::stringstream filename;
      filename << "solution_pressure_" << std::setw(5) << setfill('0') << step_
          << ".pos";
      std::cout << "writing '"<<filename.str()<<"'...";
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
          //const int numnodeele = actele->Shape();

          const DRT::Element::DiscretizationType stressdistype =
              XFLUID::getStressInterpolationType3D(actele->Shape());
          const int numeleparam =
              DRT::UTILS::getNumberOfElementNodes(stressdistype);
          BlitzMat xyze_xfemElement(DRT::UTILS::PositionArrayBlitz(actele));
          
          // create local copy of information about dofs
          const XFEM::ElementDofManager eledofman =
              dofmanagerForOutput_->constructElementDofManager(*actele,
                  numeleparam);
          
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
          
          const XFEM::DomainIntCells& domainintcells =
              ihForOutput_->GetDomainIntCells(elegid, actele->Shape());
          for (XFEM::DomainIntCells::const_iterator cell =
              domainintcells.begin(); cell != domainintcells.end(); ++cell)
          {
            BlitzVec cellvalues(DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
            XFEM::computeScalarCellNodeValues(*actele, ihForOutput_, eledofman,
                *cell, field, elementvalues, cellvalues);
            BlitzMat xyze_cell(3, cell->NumNode());
            cell->NodalPosXYZ(*actele, xyze_cell);
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
    
    {
      std::stringstream filename;
      filename << "solution_velocity_" << std::setw(5) << setfill('0') << step_
          << ".pos";
      std::cout << "writing '"<<filename.str()<<"'...";
      std::ofstream f_system(filename.str().c_str());
      
      {
        stringstream gmshfilecontent;
        gmshfilecontent << "View \" " << "Velocity Solution (Physical) \" {"
            << endl;
        for (int i=0; i<discret_->NumMyColElements(); ++i)
        {
          const DRT::Element* actele = discret_->lColElement(i);
          const int elegid = actele->Id();

          const DRT::Element::DiscretizationType stressdistype =
              XFLUID::getStressInterpolationType3D(actele->Shape());
          const int numeleparam =
              DRT::UTILS::getNumberOfElementNodes(stressdistype);
          BlitzMat xyze_xfemElement(DRT::UTILS::PositionArrayBlitz(actele));
          
          // create local copy of information about dofs
          const XFEM::ElementDofManager eledofman =
              dofmanagerForOutput_->constructElementDofManager(*actele,
                  numeleparam);
          
          vector<int> lm;
          vector<int> lmowner;
          actele->LocationVector(*(discret_), lm, lmowner);
          
          // extract local values from the global vector
          vector<double> myvelnp(lm.size());
          DRT::UTILS::ExtractMyValues(*state_.velnp_, myvelnp, lm);
          
          
          const vector<int>& dofposvelx =
              eledofman.LocalDofPosPerField(XFEM::PHYSICS::Velx);
          const vector<int>& dofposvely =
              eledofman.LocalDofPosPerField(XFEM::PHYSICS::Vely);
          const vector<int>& dofposvelz =
              eledofman.LocalDofPosPerField(XFEM::PHYSICS::Velz);
          
          const int numparam = eledofman.NumDofPerField(XFEM::PHYSICS::Velx);
          blitz::Array<double,2> elementvalues(3, numparam);
          for (int iparam=0; iparam<numparam; ++iparam)
          {
            elementvalues(0, iparam) = myvelnp[dofposvelx[iparam]];
            elementvalues(1, iparam) = myvelnp[dofposvely[iparam]];
            elementvalues(2, iparam) = myvelnp[dofposvelz[iparam]];
          }
          
          if (!ihForOutput_->ElementIntersected(elegid))
          {
            const XFEM::DomainIntCells& domainintcells =
                ihForOutput_->GetDomainIntCells(elegid, actele->Shape());
            for (XFEM::DomainIntCells::const_iterator cell =
                domainintcells.begin(); cell != domainintcells.end(); ++cell)
            {
              BlitzMat cellvalues(3, DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
              //std::cout << cellvalues << endl;
              XFEM::computeVectorCellNodeValues(*actele, ihForOutput_, eledofman,
                  *cell, XFEM::PHYSICS::Velx, elementvalues, cellvalues);
              BlitzMat xyze_cell(3, cell->NumNode());
              cell->NodalPosXYZ(*actele, xyze_cell);
              gmshfilecontent << IO::GMSH::cellWithVectorFieldToString(
                  cell->Shape(), cellvalues, xyze_cell) << endl;
            }
          }
          else
          {
            const XFEM::BoundaryIntCells& boundaryintcells =
                ihForOutput_->GetBoundaryIntCells(elegid);
            // draw boundary integration cells with values
            for (XFEM::BoundaryIntCells::const_iterator cell =
              boundaryintcells.begin(); cell != boundaryintcells.end(); ++cell)
            {
              {
                BlitzMat cellvalues(3, DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
                //std::cout << cellvalues << endl;
                XFEM::computeVectorCellNodeValues(*actele, ihForOutput_, eledofman,
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
              const XFEM::DomainIntCell cell(actele->Shape());
              XFEM::computeVectorCellNodeValues(*actele, ihForOutput_, eledofman,
                  cell, XFEM::PHYSICS::Velx, elementvalues, elevalues);
              
              const BlitzMat xyze_ele(DRT::UTILS::PositionArrayBlitz(actele));
              gmshfilecontent << IO::GMSH::cellWithVectorFieldToString(
                  actele->Shape(), elevalues, xyze_ele) << endl;
            }
            
          }
          
          if (elegid == 1)
          {
            //std::cout << elementvalues << std::endl;
            std::ofstream f;
            if (step_ <= 1)
              f.open("outflowvel.txt",std::fstream::trunc);
            else
              f.open("outflowvel.txt",std::fstream::ate | std::fstream::app);
            
            //f << step_ << " " << (-1.5*std::sin(0.1*2.0*time_* PI) * PI*0.1) << "  " << elementvalues(0,0) << endl;
            f << step_ << "  " << elementvalues(0,0) << endl;
            
            f.close();
          }
          
        }
        gmshfilecontent << "};" << endl;
        f_system << gmshfilecontent.str();
      }
      f_system.close();
      std::cout << " done" << endl;
      //exit(1);
    }
    

}
  
  

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |                                                             kue 04/07|
 -----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void XFluidImplicitTimeInt::ReadRestart(int step)
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

//  if (alefluid_)
//  {
//    reader.ReadVector(state_.dispnp_,"dispnp");
//    reader.ReadVector(state_.dispn_ , "dispn");
//    reader.ReadVector(state_.dispnm_,"dispnm");
//  }
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |                                                           chfoe 01/08|
 -----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//void XFluidImplicitTimeInt::UpdateGridv()
//{
//  // get order of accuracy of grid velocity determination
//  // from input file data
//  const int order  = params_.get<int>("order gridvel");
//
//  TIMEINT_THETA_BDF2::ComputeGridVelocity(
//      state_.dispnp_, state_.dispn_, state_.dispnm_,
//          order, step_, theta_, dta_, dtp_,
//          gridv_);
//}




//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  set initial flow field for test cases                    gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void XFluidImplicitTimeInt::SetInitialFlowField(
        RCP<DRT::Discretization> cutterdiscret,
  int whichinitialfield,
  int startfuncno
  )
{
  // create zero displacement vector to use initial position of interface
  const Epetra_Map* fluidsurface_dofcolmap = cutterdiscret->DofColMap();
  Teuchos::RCP<Epetra_Vector> idispcol     = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
    
  ComputeInterfaceAndSetDOFs(cutterdiscret,*idispcol);

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
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | evaluate error for test cases with analytical solutions   gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void XFluidImplicitTimeInt::EvaluateErrorComparedToAnalyticalSol()
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
    // choose what to assemble --- nothing
    eleparams.set("assemble matrix 1",false);
    eleparams.set("assemble matrix 2",false);
    eleparams.set("assemble vector 1",false);
    eleparams.set("assemble vector 2",false);
    eleparams.set("assemble vector 3",false);
    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("u and p at time n+1 (converged)",state_.velnp_);

    // call loop over elements
    discret_->Evaluate(eleparams,sysmat_,null,residual_,null,null);
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
void XFluidImplicitTimeInt::SolveStationaryProblem(
        RCP<DRT::Discretization> cutterdiscret
        )
{

  const Epetra_Map* fluidsurface_dofcolmap = cutterdiscret->DofColMap();
  Teuchos::RCP<Epetra_Vector> ivelcol     = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
  Teuchos::RCP<Epetra_Vector> idispcol    = LINALG::CreateVector(*fluidsurface_dofcolmap,true); // one could give a velocity here to have stationary flow over the interface
  Teuchos::RCP<Epetra_Vector> itruerescol = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
  
  ComputeInterfaceAndSetDOFs(cutterdiscret,*idispcol);

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

     // choose what to assemble
     eleparams.set("assemble matrix 1",false);
     eleparams.set("assemble matrix 2",false);
     eleparams.set("assemble vector 1",true);
     eleparams.set("assemble vector 2",false);
     eleparams.set("assemble vector 3",false);
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
    NonlinearSolve(cutterdiscret,ivelcol,idispcol,itruerescol);

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
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  calculate traction vector at (Dirichlet) boundary (public) gjb 07/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
Teuchos::RCP<Epetra_Vector> XFluidImplicitTimeInt::CalcStresses()
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
XFluidImplicitTimeInt::~XFluidImplicitTimeInt()
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
/*!
\brief calculate lift&drag forces and angular momenta

Lift and drag forces are based upon the right hand side true-residual entities
of the corresponding nodes. The contribution of the end node of a line is entirely
added to a present L&D force.

Idea of this routine:

create

map< label, set<DRT::Node*> >

which is a set of nodes to each L&D Id
nodal forces of all the nodes within one set are added to one L&D force

Notice: Angular moments obtained from lift&drag forces currently refere to the
        initial configuration, i.e. are built with the coordinates X of a particular
        node irrespective of its current position.
*/
void XFluidImplicitTimeInt::LiftDrag() const
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



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> XFluidImplicitTimeInt::IntegrateInterfaceShape(std::string condname)
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
//  if (alefluid_)
//  {
//    discret_->SetState("dispnp", state_.dispnp_);
//  }
  discret_->EvaluateCondition(eleparams,integratedshapefunc,condname);
  discret_->ClearState();

  return integratedshapefunc;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void XFluidImplicitTimeInt::UseBlockMatrix(Teuchos::RCP<std::set<int> > condelements,
                                          const LINALG::MultiMapExtractor& domainmaps,
                                          const LINALG::MultiMapExtractor& rangemaps)
{
#if 0
  sysmat_ =
    Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(domainmaps,rangemaps,108,false,true));
#else
  //const int numdim = params_.get<int>("number of velocity degrees of freedom");

  Teuchos::RCP<LINALG::BlockSparseMatrix<FLUID_UTILS::InterfaceSplitStrategy> > mat =
    Teuchos::rcp(new LINALG::BlockSparseMatrix<FLUID_UTILS::InterfaceSplitStrategy>(domainmaps,rangemaps,108,false,true));
  mat->SetCondElements(condelements);
  sysmat_ = mat;
#endif
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void XFluidImplicitTimeInt::LinearRelaxationSolve(Teuchos::RCP<Epetra_Vector> relax)
{
  //
  // This method is really stupid, but simple. We calculate the fluid
  // elements twice here. First because we need the global matrix to
  // do a linear solve. We do not have any RHS other than the one from
  // the Dirichlet condition at the FSI interface.
  //
  // After that we recalculate the matrix so that we can get the
  // reaction forces at the FSI interface easily.
  //
  // We do more work that required, but we do not need any special
  // element code to perform the steepest descent calculation. This is
  // quite a benefit as the special code in the old discretization was
  // a real nightmare.
  //

  //relax_->PutScalar(0.0);
  //interface_.InsertCondVector(ivel,relax_);

  const Epetra_Map* dofrowmap = discret_->DofRowMap();
  Teuchos::RCP<Epetra_Vector> griddisp = LINALG::CreateVector(*dofrowmap,false);

  // set the grid displacement independent of the trial value at the
  // interface
  //griddisp->Update(1., *state_.dispnp_, -1., *state_.dispn_, 0.);

  // dirichtoggle_ has already been set up

  // zero out the stiffness matrix
  sysmat_->Zero();

  // zero out residual, no neumann bc
  residual_->PutScalar(0.0);

  ParameterList eleparams;
  eleparams.set("action","calc_fluid_systemmat_and_residual");
  eleparams.set("total time",time_);
  eleparams.set("thsl",theta_*dta_);
  eleparams.set("using stationary formulation",false);
  eleparams.set("include reactive terms for linearisation",params_.get<bool>("Use reaction terms for linearisation",false));

  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("velnp",state_.velnp_);
  discret_->SetState("hist"  ,zeros_);
  discret_->SetState("dispnp", griddisp);
//  discret_->SetState("gridv", gridv_);

  // call loop over elements
  discret_->Evaluate(eleparams,sysmat_,residual_);
  discret_->ClearState();

  // finalize the system matrix
  sysmat_->Complete();

  // No, we do not want to have any rhs. There cannot be any.
  residual_->PutScalar(0.0);

  //--------- Apply dirichlet boundary conditions to system of equations
  //          residual discplacements are supposed to be zero at
  //          boundary conditions
  incvel_->PutScalar(0.0);
  LINALG::ApplyDirichlettoSystem(sysmat_,incvel_,residual_,relax,dirichtoggle_);

  //-------solve for residual displacements to correct incremental displacements
  solver_.Solve(sysmat_->EpetraOperator(),incvel_,residual_,true,true);

  // and now we need the reaction forces

  // zero out the stiffness matrix
  sysmat_->Zero();

  // zero out residual, no neumann bc
  residual_->PutScalar(0.0);

  eleparams.set("action","calc_fluid_systemmat_and_residual");
  eleparams.set("thsl",theta_*dta_);
  eleparams.set("using stationary formulation",false);

  // set vector values needed by elements
  discret_->ClearState();
  //discret_->SetState("velnp",incvel_);
  discret_->SetState("velnp",state_.velnp_);
  discret_->SetState("hist"  ,zeros_);
  discret_->SetState("dispnp", griddisp);
//  discret_->SetState("gridv", gridv_);

  // call loop over elements
  discret_->Evaluate(eleparams,sysmat_,residual_);
  discret_->ClearState();

  density_ = eleparams.get("density", 0.0);

  // finalize the system matrix
  sysmat_->Complete();

  if (sysmat_->Apply(*incvel_, *trueresidual_)!=0)
    dserror("sysmat_->Apply() failed");
  trueresidual_->Scale(-density_/dta_/theta_);
}


#endif /* CCADISCRET       */
