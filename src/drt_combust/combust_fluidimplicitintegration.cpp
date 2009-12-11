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
  xparams_(params.sublist("XFEM")),
//  output_ (output),
  output_ (rcp(new IO::DiscretizationWriter(actdis))), // so ist es bei Axel
  myrank_(discret_->Comm().MyPID()),
  cout0_(discret_->Comm(), std::cout),
  combusttype_(Teuchos::getIntegralValue<INPAR::COMBUST::CombustionType>(params_.sublist("COMBUSTION FLUID"),"COMBUSTTYPE")),
  flamespeed_(params_.sublist("COMBUSTION FLUID").get<double>("LAMINAR_FLAMESPEED")),
  nitschevel_(params_.sublist("COMBUSTION FLUID").get<double>("NITSCHE_VELOCITY")),
  nitschepres_(params_.sublist("COMBUSTION FLUID").get<double>("NITSCHE_PRESSURE")),
  step_(0),
  time_(0.0),
  stepmax_ (params_.get<int>   ("max number timesteps")),
  maxtime_ (params_.get<double>("total time")),
  dta_     (params_.get<double> ("time step size")),
  dtp_     (params_.get<double> ("time step size")),
  timealgo_(params_.get<FLUID_TIMEINTTYPE>("time int algo")),
  theta_   (params_.get<double>("theta")),
  extrapolationpredictor_(params.get("do explicit predictor",true)),
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
  // set time integration parameters for stationary simulation
  //------------------------------------------------------------------------------------------------
  if (timealgo_ == timeint_stationary)
  {
    dta_ = 1.0;
    dtp_ = 1.0;
    theta_ = 1.0;
    cout0_ << "parameters 'theta' and 'time step size' have been set to 1.0 for stationary problem " << endl;
  }
  //------------------------------------------------------------------------------------------------
  // future: connect degrees of freedom for periodic boundary conditions                 henke 01/09
  //------------------------------------------------------------------------------------------------

  //------------------------------------------------------------------------------------------------
  // prepare XFEM (initial degree of freedom management)
  //------------------------------------------------------------------------------------------------
  physprob_.xfemfieldset_.clear();
  // declare physical fields to be enriched (XFEM fields)
  physprob_.xfemfieldset_.insert(XFEM::PHYSICS::Velx);
  physprob_.xfemfieldset_.insert(XFEM::PHYSICS::Vely);
  physprob_.xfemfieldset_.insert(XFEM::PHYSICS::Velz);
  physprob_.xfemfieldset_.insert(XFEM::PHYSICS::Pres);
  physprob_.elementAnsatzp_ = rcp(new COMBUST::CombustElementAnsatz());

  // create dummy instance of interfacehandle holding no flamefront and hence no integration cells
  Teuchos::RCP<COMBUST::InterfaceHandleCombust> ihdummy = rcp(new COMBUST::InterfaceHandleCombust(discret_,Teuchos::null,Teuchos::null));
  // create dummy instance of dof manager assigning standard enrichments to all nodes
  const Teuchos::RCP<XFEM::DofManager> dofmanagerdummy = rcp(new XFEM::DofManager(ihdummy,physprob_.xfemfieldset_,xparams_));

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
//  std::cout << "dofrowmap: " << *dofrowmap << "\n" << endl;

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

} // CombustFluidImplicitTimeInt::CombustFluidImplicitTimeInt

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
  step_ += 1;
  time_ += dta_;

  if (params_.get<FLUID_TIMEINTTYPE>("time int algo") == timeint_stationary)
  {
    // for stationary problems PrepareTimeStep() should only be called for output related reasons
    cout << "Warning: 'time' and 'time step' are set to 1.0 and 1 for output control file" << endl;
    timealgo_ = timeint_stationary;
    step_ = 1;
    time_ = 1.0;
    theta_ = 1.0;
  }
  else if (params_.get<FLUID_TIMEINTTYPE>("time int algo") == timeint_bdf2)
  {
    // do a backward Euler step for the first time step
    if (step_==1)
    {
      timealgo_ = timeint_one_step_theta;
      theta_ = params_.get<double>("start theta");
    }
    else
    {
      timealgo_ = timeint_bdf2;
      // for BDF2, theta is set by the time-step sizes, 2/3 for const. dt
      theta_ = (dta_+dtp_)/(2.0*dta_ + dtp_);
    }
  }
  else
  {
    timealgo_ = params_.get<FLUID_TIMEINTTYPE>("time int algo");
    theta_ = params_.get<double>("theta");
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
       const Teuchos::RCP<COMBUST::InterfaceHandleCombust>& interfacehandle)
{
  // information about interface is imported via ADAPTER::FluidCombust::ImportInterface()

  /* momentan gebe ich den ElementAnsatz (Lagrange Multiplier Zeug) nicht an den DofManager weiter,
   * vielleicht brauche ich es aber. Hängt das hier eigentlich nicht direkt von InputParametern ab? */
//  const COMBUST::CombustElementAnsatz elementAnsatz;

  // build instance of DofManager with information about the interface from the interfacehandle
  // remark: DofManager is rebuilt in every inter-field iteration step, because number and position
  // of enriched degrees of freedom change constantly.
  const Teuchos::RCP<XFEM::DofManager> dofmanager = rcp(new XFEM::DofManager(interfacehandle,physprob_.xfemfieldset_,xparams_));

  // save dofmanager to be able to plot Gmsh stuff in Output()
  dofmanagerForOutput_ = dofmanager;

  // tell elements about the dofs and the integration
  TransferDofInformationToElements(interfacehandle, dofmanager);

  // print global and element dofmanager to Gmsh
  std::cout<< "Dofmanager and InterfaceHandle to Gmsh" << std::endl;
  dofmanager->toGmsh(step_);
  interfacehandle->toGmsh(step_);


  // get old dofmap, compute new one and get the new one, too
  const Epetra_Map olddofrowmap = *discret_->DofRowMap();
  // assign degrees of freedom
  // remark: - assign degrees of freedom (first slot)
  //         - build geometry for (Neumann) boundary conditions (third slot);
  //           without Neumann boundary conditions Fillcomplete(true,false,false) will also work
  discret_->FillComplete(true,false,true);
  const Epetra_Map& newdofrowmap = *discret_->DofRowMap();

  discret_->ComputeNullSpaceIfNecessary(solver_.Params());

  {
  const std::map<XFEM::DofKey<XFEM::onNode>, XFEM::DofGID> oldNodalDofDistributionMap(state_.nodalDofDistributionMap_);
  const std::map<XFEM::DofKey<XFEM::onElem>, XFEM::DofGID> oldElementalDofDistributionMap(state_.elementalDofDistributionMap_);
  dofmanager->fillDofDistributionMaps(
      state_.nodalDofDistributionMap_,
      state_.elementalDofDistributionMap_);

  // create switcher
  const XFEM::DofDistributionSwitcher dofswitch(
          interfacehandle, dofmanager,
          olddofrowmap, newdofrowmap,
          oldNodalDofDistributionMap, state_.nodalDofDistributionMap_,
          oldElementalDofDistributionMap, state_.elementalDofDistributionMap_
          );

  // --------------------------------------------
  // switch state vectors to new dof distribution
  // --------------------------------------------

    cout0_ << " Initialize system vectors..." << endl;
  // accelerations at time n and n-1
  dofswitch.mapVectorToNewDofDistributionCombust(state_.accnp_);
  dofswitch.mapVectorToNewDofDistributionCombust(state_.accn_);

  // velocities and pressures at time n+1, n and n-1
  dofswitch.mapVectorToNewDofDistributionCombust(state_.velnp_); // use old velocity as start value
  dofswitch.mapVectorToNewDofDistributionCombust(state_.veln_);
  dofswitch.mapVectorToNewDofDistributionCombust(state_.velnm_);
  }
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
  const int         itemax = params_.get<int>("max nonlin iter steps");
  bool              stopnonliniter = false;

  double dtsolve = 0.0;
  double dtele   = 0.0;

  // action for elements
  if (timealgo_!=timeint_stationary and theta_ < 1.0)
  {
    cout0_ << "* Warning! Works reliable only for Backward Euler time discretization! *" << endl;
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

      // flag for type of combustion problem
      eleparams.set("combusttype",combusttype_);
      eleparams.set("flamespeed",flamespeed_);
      eleparams.set("nitschevel",nitschevel_);
      eleparams.set("nitschepres",nitschepres_);

      // other parameters that might be needed by the elements
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

      // set vector values needed by elements
      discret_->ClearState();
      discret_->SetState("velnp",state_.velnp_);
      discret_->SetState("veln" ,state_.veln_);
      discret_->SetState("velnm",state_.velnm_);
      discret_->SetState("accn" ,state_.accn_);

      // convergence check at itemax is skipped for speedup if
      // CONVCHECK is set to L_2_norm_without_residual_at_itemax
      if ((itnum != itemax)
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
      dtele=ds_cputime()-tcpu;

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
                 itnum,itemax,ittol,vresnorm,presnorm,fullresnorm,
                 incvelnorm_L2/velnorm_L2,incprenorm_L2/prenorm_L2,incfullnorm_L2/fullnorm_L2);
          printf(" (ts=%10.3E,te=%10.3E",dtsolve,dtele);
          printf(")\n");
          printf("+------------+-------------------+--------------+--------------+--------------+--------------+--------------+--------------+\n");

          FILE* errfile = params_.get<FILE*>("err file");
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

        FILE* errfile = params_.get<FILE*>("err file");
        if (errfile!=NULL)
        {
          fprintf(errfile,"fluid unconverged solve:   %3d/%3d  tol=%10.3E[L_2 ]  vres=%10.3E  pres=%10.3E  vinc=%10.3E  pinc=%10.3E\n",
                  itnum,itemax,ittol,vresnorm,presnorm,
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
      const double tcpusolve=ds_cputime();

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
      dtsolve = ds_cputime()-tcpusolve;
    }

    //------------------------------------------------ update (u,p) trial
    state_.velnp_->Update(1.0,*incvel_,1.0);
  }
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
//    Teuchos::RCP<Epetra_Vector> velnp_out = dofmanagerForOutput_->fillPhysicalOutputVector(
//        *state_.velnp_, dofset_out_, state_.nodalDofDistributionMap_, physprob_.fieldset_);

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
    output_->WriteVector("velnp", state_.velnp_);
    output_->WriteVector("veln" , state_.veln_);
    output_->WriteVector("velnm", state_.velnm_);
    output_->WriteVector("accnp", state_.accnp_);
    output_->WriteVector("accn" , state_.accn_);
  }

//  if (discret_->Comm().NumProc() == 1)
  {
    OutputToGmsh(step_, time_);
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
void FLD::CombustFluidImplicitTimeInt::OutputToGmsh(
    const int step,
    const double time
    ) const
{
  const Teuchos::ParameterList& xfemparams = DRT::Problem::Instance()->XFEMGeneralParams();
  const bool gmshdebugout = getIntegralValue<int>(xfemparams,"GMSH_DEBUG_OUT")==1;

  const bool screen_out = true;

  // get a copy on columnmn parallel distribution
  Teuchos::RCP<const Epetra_Vector> output_col_velnp = DRT::UTILS::GetColVersionOfRowVector(discret_, state_.velnp_);


  if (gmshdebugout and (this->physprob_.xfemfieldset_.find(XFEM::PHYSICS::Pres) != this->physprob_.xfemfieldset_.end()))
  {
    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("solution_field_pressure", step, 50, screen_out, discret_->Comm().MyPID());
    std::ofstream gmshfilecontent(filename.c_str());

    const XFEM::PHYSICS::Field field = XFEM::PHYSICS::Pres;
    {
      gmshfilecontent << "View \" " << "Pressure Solution (Physical) \" {\n";
      for (int i=0; i<discret_->NumMyRowElements(); ++i)
      {
        const DRT::Element* actele = discret_->lRowElement(i);

        // create local copy of information about dofs
        const XFEM::ElementDofManager eledofman(*actele,physprob_.elementAnsatzp_->getElementAnsatz(actele->Shape()),*dofmanagerForOutput_);

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
        const XFEM::ElementDofManager eledofman(*actele,physprob_.elementAnsatzp_->getElementAnsatz(actele->Shape()),*dofmanagerForOutput_);

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
        const XFEM::ElementDofManager eledofman(*actele,physprob_.elementAnsatzp_->getElementAnsatz(actele->Shape()),*dofmanagerForOutput_);

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
    PlotVectorFieldToGmsh(DRT::UTILS::GetColVersionOfRowVector(discret_, state_.velnp_), "solution_field_velocity_np","Velocity Solution (Physical) n+1",true, step, time);
    if (timealgo_ != timeint_stationary)
    {
//      PlotVectorFieldToGmsh(DRT::UTILS::GetColVersionOfRowVector(discret_, state_.velnm_), "solution_field_velocity_nm","Velocity Solution (Physical) n-1",false, step, time);
//      PlotVectorFieldToGmsh(DRT::UTILS::GetColVersionOfRowVector(discret_, state_.veln_), "solution_field_velocity_n","Velocity Solution (Physical) n",false, step, time);
//      PlotVectorFieldToGmsh(DRT::UTILS::GetColVersionOfRowVector(discret_, state_.accn_), "solution_field_acceleration_n","Acceleration Solution (Physical) n",false, step, time);
      PlotVectorFieldToGmsh(DRT::UTILS::GetColVersionOfRowVector(discret_, state_.accnp_), "solution_field_acceleration_np","Acceleration Solution (Physical) n+1",false, step, time);
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
  cout << "PlotVectorFieldToGmsh" << endl;

  const Teuchos::ParameterList& xfemparams = DRT::Problem::Instance()->XFEMGeneralParams();
  const bool gmshdebugout = getIntegralValue<int>(xfemparams,"GMSH_DEBUG_OUT")==1;

  const bool screen_out = true;

  if (gmshdebugout)
  {
    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles(filestr, step, 50, screen_out, discret_->Comm().MyPID());
    std::ofstream gmshfilecontent(filename.c_str());

    {
      gmshfilecontent << "View \" " << name_in_gmsh << "\" {\n";
      for (int i=0; i<discret_->NumMyRowElements(); ++i)
      {
        const DRT::Element* actele = discret_->lRowElement(i);

        // create local copy of information about dofs
        const COMBUST::CombustElementAnsatz elementAnsatz;
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
    int whichinitialfield,
    int startfuncno)
{
  std::cout << "SetInitialFlowField() wird ausgeführt!" << endl;

  // create zero displacement vector to use initial position of interface
  {
    //IncorporateInterface();
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
  }
  else
  {
    dserror("no other initial fields than zero, function and beltrami are available up to now");
  }

  return;
} // CombustFluidImplicitTimeInt::SetInitialFlowField()

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

#endif // CCADISCRET
