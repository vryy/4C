/*----------------------------------------------------------------------*/
/*!
\file scatra_timint_ost.cpp
\brief One-Step-Theta time-integration scheme

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#include "scatra_timint_ost.H"
#include "../drt_scatra_ele/scatra_ele_action.H"
#include "turbulence_hit_scalar_forcing.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include "../drt_inpar/inpar_elch.H"
#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"
#include "../linalg/linalg_solver.H"
#include "../drt_fluid_turbulence/dyn_smag.H"
#include "../drt_fluid_turbulence/dyn_vreman.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"


/*----------------------------------------------------------------------*
 |  Constructor (public)                                      gjb 08/08 |
 *----------------------------------------------------------------------*/
SCATRA::TimIntOneStepTheta::TimIntOneStepTheta(
  Teuchos::RCP<DRT::Discretization>      actdis,
  Teuchos::RCP<LINALG::Solver>           solver,
  Teuchos::RCP<Teuchos::ParameterList>   params,
  Teuchos::RCP<Teuchos::ParameterList>   extraparams,
  Teuchos::RCP<IO::DiscretizationWriter> output)
: ScaTraTimIntImpl(actdis,solver,params,extraparams,output),
  theta_(params_->get<double>("THETA")),
  fsphinp_(Teuchos::null)
{
  // DO NOT DEFINE ANY STATE VECTORS HERE (i.e., vectors based on row or column maps)
  // this is important since we have problems which require an extended ghosting
  // this has to be done before all state vectors are initialized
  return;
}


/*----------------------------------------------------------------------*
 |  initialize time integration                         rasthofer 09/13 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::Init()
{
  // initialize base class
  ScaTraTimIntImpl::Init();

  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  // -------------------------------------------------------------------
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // fine-scale vector at time n+1
  if (fssgd_ != INPAR::SCATRA::fssugrdiff_no or turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales)
    fsphinp_ = LINALG::CreateVector(*dofrowmap,true);

  // -------------------------------------------------------------------
  // set element parameters
  // -------------------------------------------------------------------
  // note: - this has to be done before element routines are called
  //       - order is important here: for savety checks in SetElementGeneralScaTraParameter(),
  //         we have to konw the time-integration parameters
  SetElementTimeParameter();
  SetElementGeneralScaTraParameter();
  SetElementTurbulenceParameter();

  //TODO: SCATRA_ELE_CLEANING
  // initialize time-dependent electrode kinetics variables (galvanostatic mode or double layer contribution)
  //if (IsElch(scatratype_))
  //  ComputeTimeDerivPot0(true);

  // Important: this adds the required ConditionID's to the single conditions.
  // It is necessary to do this BEFORE ReadRestart() is called!
  // Output to screen and file is suppressed
  //OutputElectrodeInfo(false,false);

  // setup krylov
  PrepareKrylovProjection();

  // -------------------------------------------------------------------
  // initialize forcing for homogeneous isotropic turbulence
  // -------------------------------------------------------------------
  // note: this constructor has to be called after the forcing_ vector has
  //       been initialized; this is done in ScaTraTimIntImpl::Init() called before
  if (special_flow_ == "scatra_forced_homogeneous_isotropic_turbulence")
  {
    if (extraparams_->sublist("TURBULENCE MODEL").get<std::string>("SCALAR_FORCING")=="isotropic")
    {
      homisoturb_forcing_ = Teuchos::rcp(new SCATRA::HomIsoTurbScalarForcing(this));
      // initialize forcing algorithm
      homisoturb_forcing_->SetInitialSpectrum(DRT::INPUT::IntegralValue<INPAR::SCATRA::InitialField>(*params_,"INITIALFIELD"));
    }
  }

  return;
}


/*----------------------------------------------------------------------*
| Destructor dtor (public)                                    gjb 08/08 |
*-----------------------------------------------------------------------*/
SCATRA::TimIntOneStepTheta::~TimIntOneStepTheta()
{
  return;
}


/*----------------------------------------------------------------------*
 |  set time parameter for element evaluation (usual call)   ehrl 11/13 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::SetElementTimeParameter()
{
  Teuchos::ParameterList eleparams;

  eleparams.set<int>("action",SCATRA::set_time_parameter);
  // set type of scalar transport problem (after preevaluate evaluate, which need scatratype is called)
  eleparams.set<int>("scatratype",scatratype_);

  eleparams.set<bool>("using generalized-alpha time integration",false);
  eleparams.set<bool>("using stationary formulation",false);
  eleparams.set<bool>("incremental solver",incremental_);

  eleparams.set<double>("time-step length",dta_);
  eleparams.set<double>("total time",time_);
  eleparams.set<double>("time factor",theta_*dta_);
  eleparams.set<double>("alpha_F",1.0);

  // call standard loop over elements
  discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);

  return;
}


/*----------------------------------------------------------------------*
 |  set time parameter for element evaluation                ehrl 11/13 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::SetElementTimeParameterForForcedIncrementalSolve()
{
  Teuchos::ParameterList eleparams;

  eleparams.set<int>("action",SCATRA::set_time_parameter);
  // set type of scalar transport problem (after preevaluate evaluate, which need scatratype is called)
  eleparams.set<int>("scatratype",scatratype_);

  eleparams.set<bool>("using generalized-alpha time integration",false);
  eleparams.set<bool>("using stationary formulation",false);
  // this is important to have here and the only difference compared to SetElementTimeParameter()
  eleparams.set<bool>("incremental solver",true);

  eleparams.set<double>("time-step length",dta_);
  eleparams.set<double>("total time",time_);
  eleparams.set<double>("time factor",theta_*dta_);
  eleparams.set<double>("alpha_F",1.0);

  // call standard loop over elements
  discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);

  return;
}


/*----------------------------------------------------------------------*
 |  set time parameter for element evaluation                ehrl 11/13 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::SetElementTimeParameterInitial()
{
  Teuchos::ParameterList eleparams;

  eleparams.set<int>("action",SCATRA::set_time_parameter);
  // set type of scalar transport problem (after preevaluate evaluate, which need scatratype is called)
  eleparams.set<int>("scatratype",scatratype_);

  eleparams.set<bool>("using generalized-alpha time integration",false);
  eleparams.set<bool>("using stationary formulation",false);
  eleparams.set<bool>("incremental solver",incremental_);

  eleparams.set<double>("time-step length",dta_);
  eleparams.set<double>("total time",0.0); // only different!
  eleparams.set<double>("time factor",theta_*dta_);
  eleparams.set<double>("alpha_F",1.0);

  // call standard loop over elements
  discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);

  return;
}


/*----------------------------------------------------------------------*
| Print information about current time step to screen                   |
*-----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::PrintTimeStepInfo()
{
  if (myrank_==0)
  {
    IO::cout << "TIME: "
             << std::setw(11) << std::setprecision(4) << std::scientific << time_ << "/"
             << std::setw(11) << std::setprecision(4) << std::scientific << maxtime_ << "  DT = "
             << std::setw(11) << std::setprecision(4) << std::scientific << dta_ << "  "
             << MethodTitle() << " (theta = "
             << std::setw(3)  << std::setprecision(2) << theta_ << ") STEP = "
             << std::setw(4) << step_ << "/" << std::setw(4) << stepmax_ << IO::endl;
  }
  return;
}


/*----------------------------------------------------------------------*
 | set part of the residual vector belonging to the old timestep        |
 |                                                            gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::SetOldPartOfRighthandside()
{
  // hist_ = phin_ + dt*(1-Theta)*phidtn_
  hist_->Update(1.0, *phin_, dta_*(1.0-theta_), *phidtn_, 0.0);

  // for electrochemical applications
  ElectrodeKineticsSetOldPartOfRHS();

  return;
}


/*----------------------------------------------------------------------*
 | perform an explicit predictor step                         gjb 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::ExplicitPredictor()
{
  phinp_->Update(dta_, *phidtn_,1.0);

  // for the electric potential we just use the 'old' value of last time step
  Teuchos::RCP<Epetra_Vector> onlypot = splitter_->ExtractCondVector(phin_);
  splitter_->InsertCondVector(onlypot, phinp_);

  return;
}


/*----------------------------------------------------------------------*
 | add actual Neumann loads                                             |
 | scaled with a factor resulting from time discretization     vg 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::AddNeumannToResidual()
{
  residual_->Update(theta_*dta_,*neumann_loads_,1.0);
  return;
}


/*----------------------------------------------------------------------*
 | AVM3-based scale separation                                 vg 03/09 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::AVM3Separation()
{
  // time measurement: avm3
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:            + avm3");

  // AVM3 separation
  Sep_->Multiply(false,*phinp_,*fsphinp_);

  // set fine-scale vector
  discret_->SetState("fsphinp",fsphinp_);

  return;
}


/*----------------------------------------------------------------------*
 | dynamic Smagorinsky model                           rasthofer  08/12 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::DynamicComputationOfCs()
{
  if (turbmodel_==INPAR::FLUID::dynamic_smagorinsky)
  {
    // perform filtering and computation of Prt
    // compute averaged values for LkMk and MkMk
    const Teuchos::RCP<const Epetra_Vector> dirichtoggle = DirichletToggle();
    DynSmag_->ApplyFilterForDynamicComputationOfPrt(convel_,phinp_,0.0,dirichtoggle,*extraparams_);
  }

  return;
}


/*----------------------------------------------------------------------*
 | dynamic Smagorinsky model                           krank  09/13     |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::DynamicComputationOfCv()
{
  if (turbmodel_==INPAR::FLUID::dynamic_vreman)
  {
    const Teuchos::RCP<const Epetra_Vector> dirichtoggle = DirichletToggle();
    Vrem_->ApplyFilterForDynamicComputationOfDt(convel_,phinp_,0.0,dirichtoggle,*extraparams_);
  }

  return;
}


/*----------------------------------------------------------------------*
 | add parameters specific for time-integration scheme         vg 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::AddTimeIntegrationSpecificVectors()
{
  discret_->SetState("hist",hist_);
  discret_->SetState("phinp",phinp_);

  return;
}


/*----------------------------------------------------------------------*
 | compute time derivative                                     vg 09/09 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::ComputeTimeDerivative()
{
  // time derivative of phi:
  // phidt(n+1) = (phi(n+1)-phi(n)) / (theta*dt) + (1-(1/theta))*phidt(n)
  const double fact1 = 1.0/(theta_*dta_);
  const double fact2 = 1.0 - (1.0/theta_);
  phidtnp_->Update(fact2,*phidtn_,0.0);
  phidtnp_->Update(fact1,*phinp_,-fact1,*phin_,1.0);

  // We know the first time derivative on Dirichlet boundaries
  // so we do not need an approximation of these values!
  // However, we do not want to break the linear relationship
  // as stated above. We do not want to set Dirichlet values for
  // dependent values like phidtnp_. This turned out to be inconsistent.
  // ApplyDirichletBC(time_,Teuchos::null,phidtnp_);

  return;
}


/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                            gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::Update(const int num)
{
  // compute time derivative at time n+1
  ComputeTimeDerivative();

  // compute flux vector field for later output BEFORE time shift of results
  // is performed below !!
  if (writeflux_!=INPAR::SCATRA::flux_no)
  {
    if (DoOutput() or DoBoundaryFluxStatistics())
      flux_ = CalcFlux(true, num);
  }

  // after the next command (time shift of solutions) do NOT call
  // ComputeTimeDerivative() anymore within the current time step!!!

  // solution of this step becomes most recent solution of the last step
  phin_ ->Update(1.0,*phinp_,0.0);

  // time deriv. of this step becomes most recent time derivative of
  // last step
  phidtn_->Update(1.0,*phidtnp_,0.0);

  // call time update of forcing routine
  if (homisoturb_forcing_ != Teuchos::null)
    homisoturb_forcing_->TimeUpdateForcing();

  // TODO: SCATRA_ELE_CLEANING: muss weg
  // potential time update of time-dependent materials
//  ElementMaterialTimeUpdate();

  return;
}


/*----------------------------------------------------------------------*
 | write additional data required for restart                 gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::OutputRestart()
{
  // additional state vectors that are needed for One-Step-Theta restart
  output_->WriteVector("phidtn", phidtn_);
  output_->WriteVector("phin", phin_);

  // for elch problems with moving boundary
  if (isale_)
    output_->WriteVector("trueresidual", trueresidual_);

  //TODO: SCATRA_ELE_CLEANING: muss weg
  if (scatratype_ == INPAR::SCATRA::scatratype_cardio_monodomain)
  {
   output_->WriteMesh(step_,time_); // add info to control file for reading all variables in restart
  }

  return;
}


/*----------------------------------------------------------------------*
 |                                                            gjb 08/08 |
 -----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::ReadRestart(int step)
{
  IO::DiscretizationReader reader(discret_,step);
  time_ = reader.ReadDouble("time");
  step_ = reader.ReadInt("step");

  if (myrank_==0)
    std::cout<<"Reading ScaTra restart data (time="<<time_<<" ; step="<<step_<<")"<<std::endl;

  // read state vectors that are needed for One-Step-Theta restart
  reader.ReadVector(phinp_, "phinp");
  reader.ReadVector(phin_,  "phin");
  reader.ReadVector(phidtn_,"phidtn");

  // for elch problems with moving boundary
  if(isale_)
    reader.ReadVector(trueresidual_, "trueresidual");

  if (fssgd_ != INPAR::SCATRA::fssugrdiff_no or
      turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales)
    AVM3Preparation();

  if (scatratype_ == INPAR::SCATRA::scatratype_cardio_monodomain)
  {
    reader.ReadVector(activation_time_np_, "activation_time_np");
    reader.ReadMesh(step); // Read all saved data in nodes and elements und call nodal and element Unpacking each global variable has to be read
  }

  return;
}


/*----------------------------------------------------------------------*
 | Initialization procedure before the first time step        gjb 08/08 |
 -----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::PrepareFirstTimeStep()
{
  // evaluate Dirichlet boundary conditions at time t=0
  // the values should match your initial field at the boundary!
  //ApplyDirichletBC(time_,phin_,phidtn_);
  ApplyDirichletBC(time_,phin_,Teuchos::null);

  // compute initial field for electric potential (ELCH)
  CalcInitialPotentialField();

  // for calculation of initial time derivative, we have to switch off all stabilization and
  // turbulence modeling terms
  // therefore, we have another PerEvaluate call here
  Teuchos::ParameterList eleparams;

  eleparams.set<int>("action",SCATRA::set_general_scatra_parameter);

  // set type of scalar transport problem
  eleparams.set<int>("scatratype",scatratype_);

  eleparams.set<int>("form of convective term",convform_);
  eleparams.set("isale",isale_);

  // set flag for writing the flux vector fields
  eleparams.set<int>("writeflux",writeflux_);
  //! set vector containing ids of scalars for which flux vectors are calculated
  eleparams.set<Teuchos::RCP<std::vector<int> > >("writefluxids",writefluxids_);

  // parameters for stabilization
  eleparams.sublist("STABILIZATION") = params_->sublist("STABILIZATION");
  Teuchos::setStringToIntegralParameter<int>("STABTYPE",
      "no_stabilization",
      "type of stabilization (if any)",
      Teuchos::tuple<std::string>("no_stabilization"),
      Teuchos::tuple<std::string>("Do not use any stabilization"),
      Teuchos::tuple<int>(
          INPAR::SCATRA::stabtype_no_stabilization),
          &eleparams.sublist("STABILIZATION"));
  DRT::INPUT::BoolParameter("SUGRVEL","no","potential incorporation of subgrid-scale velocity",&eleparams.sublist("STABILIZATION"));
  DRT::INPUT::BoolParameter("ASSUGRDIFF","no",
        "potential incorporation of all-scale subgrid diffusivity (a.k.a. discontinuity-capturing) term",&eleparams.sublist("STABILIZATION"));

  // call standard loop over elements
  discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);

  // we also have to modify the time-parameter list
  SetElementTimeParameterForForcedIncrementalSolve();

  // add turbulence list here and set model to no model
  Teuchos::ParameterList eleturbparams;

  eleturbparams.set<int>("action",SCATRA::set_turbulence_scatra_parameter);

  // set type of scalar transport problem
  eleturbparams.set<int>("scatratype",scatratype_);

  eleturbparams.sublist("TURBULENCE MODEL") = extraparams_->sublist("TURBULENCE MODEL");
  Teuchos::setStringToIntegralParameter<int>(
      "PHYSICAL_MODEL",
      "no_model",
      "Classical LES approaches require an additional model for\nthe turbulent viscosity.",
      Teuchos::tuple<std::string>("no_model"),
      Teuchos::tuple<std::string>("If classical LES is our turbulence approach, this is a contradiction and should cause a dserror."),
      Teuchos::tuple<int>(INPAR::FLUID::no_model),
      &eleturbparams.sublist("TURBULENCE MODEL"));

  // set model-dependent parameters
  eleturbparams.sublist("SUBGRID VISCOSITY") = extraparams_->sublist("SUBGRID VISCOSITY");
  // and set parameters for multifractal subgrid-scale modeling
  eleturbparams.sublist("MULTIFRACTAL SUBGRID SCALES") = extraparams_->sublist("MULTIFRACTAL SUBGRID SCALES");

  eleturbparams.set<bool>("turbulent inflow",turbinflow_);

  eleturbparams.set<int>("fs subgrid diffusivity",INPAR::SCATRA::fssugrdiff_no);

  // call standard loop over elements
  discret_->Evaluate(eleturbparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);

  // compute time derivative of phi at time t=0
  CalcInitialPhidt();

  // and finally undo our temporary settings
  SetElementGeneralScaTraParameter();
  SetElementTimeParameter();
  SetElementTurbulenceParameter();

  return;
}




