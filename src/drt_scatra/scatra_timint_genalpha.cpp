/*----------------------------------------------------------------------*/
/*!
\file scatra_timint_genalpha.cpp
\brief Generalized-alpha time-integration scheme

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>
*/
/*----------------------------------------------------------------------*/

#include "scatra_timint_genalpha.H"
#include "scatra_ele_action.H"
#include "scatra_utils.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include "../drt_io/io.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"
#include "../drt_fluid/dyn_smag.H"


/*----------------------------------------------------------------------*
 |  Constructor (public)                                       vg 11/08 |
 *----------------------------------------------------------------------*/
SCATRA::TimIntGenAlpha::TimIntGenAlpha(
    Teuchos::RCP<DRT::Discretization>      actdis,
  Teuchos::RCP<LINALG::Solver>           solver,
  Teuchos::RCP<Teuchos::ParameterList>   params,
  Teuchos::RCP<Teuchos::ParameterList>   extraparams,
  Teuchos::RCP<IO::DiscretizationWriter> output)
: ScaTraTimIntImpl(actdis,solver,params,extraparams,output),
  alphaM_(params_->get<double>("ALPHA_M")),
  alphaF_(params_->get<double>("ALPHA_F")),
  gamma_ (params_->get<double>("GAMMA"))
{
  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  // -------------------------------------------------------------------
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // Vectors passed to the element
  // -----------------------------

  // scalar at times n+alpha_F and n+alpha_M
  phiaf_ = LINALG::CreateVector(*dofrowmap,true);
  phiam_ = LINALG::CreateVector(*dofrowmap,true);

  // temporal derivative of scalar at times n+1, n and n+alpha_M
  phidtnp_ = LINALG::CreateVector(*dofrowmap,true);
  phidtn_  = LINALG::CreateVector(*dofrowmap,true);
  phidtam_ = LINALG::CreateVector(*dofrowmap,true);

  // compute specific time factor for generalized-alpha time integration:
  // genalphatimefac = gamma*alpha_F/alpha_M
  if (alphaM_ < EPS12) dserror("factor alpha_M lower than or equal zero");
  genalphafac_ = gamma_/alphaM_;

  // fine-scale vector at time n+alpha_F
  if (fssgd_ != INPAR::SCATRA::fssugrdiff_no or turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales)
    fsphiaf_ = LINALG::CreateVector(*dofrowmap,true);

  // initialize time-dependent electrode kinetics variables (galvanostatic mode)
  ElectrodeKineticsTimeUpdate(true);

  // Important: this adds the required ConditionID's to the single conditions.
  // It is necessary to do this BEFORE ReadRestart() is called!
  // Output to screen and file is suppressed
  OutputElectrodeInfo(false,false);

  // for initializing phiaf_, phiam based on the initial field that was
  // set for phinp_, phin_ in the TimInt base class constructor
  ComputeIntermediateValues();

  return;
}


/*----------------------------------------------------------------------*
| Destructor dtor (public)                                     vg 11/08 |
*-----------------------------------------------------------------------*/
SCATRA::TimIntGenAlpha::~TimIntGenAlpha()
{
  return;
}


/*----------------------------------------------------------------------*
 | set part of the residual vector belonging to the old timestep        |
 |                                                             vg 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::SetOldPartOfRighthandside()
{
  // calculation of history vector only for non-incremental formulation:
  // (History vector is used in both cases, but in incremental case, it
  // contains time derivatives of scalar, see below.)
  // hist_ = phin_ + dt*(1-(gamma/alpha_M))*phidtn_
  if (not incremental_)
    hist_->Update(1.0, *phin_, dta_*(1.0-genalphafac_), *phidtn_, 0.0);

  return;
}


/*----------------------------------------------------------------------*
 | perform an explicit predictor step                          vg 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::ExplicitPredictor()
{
  // constant predictor
  phinp_->Update(1.0,*phin_,0.0);
  return;
}


/*----------------------------------------------------------------------*
 | predict thermodynamic pressure and time derivative          vg 12/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::PredictThermPressure()
{
  // same-thermodynamic-pressure predictor (not required to be performed,
  // since we just updated the thermodynamic pressure, and thus,
  // thermpressnp_ = thermpressn_)
  // prediction of time derivative:
  double fact = (gamma_-1.0)/gamma_;
  thermpressdtnp_ = fact*thermpressdtn_;

  // same-thermodynamic-pressure-derivative predictor (currrently not used)
  //thermpressnp_ += dta_*thermpressdtn_;
  // prediction of time derivative not required (would also not be required
  // to be performed, since we just updated the time derivatives of density,
  // and thus, thermpressdtnp_ = thermpressdtn_)

  return;
}


/*----------------------------------------------------------------------*
 | compute values at intermediate time steps                   vg 09/09 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::ComputeIntermediateValues()
{
  // compute phi at n+alpha_F and n+alpha_M
  phiaf_->Update(alphaF_,*phinp_,(1.0-alphaF_),*phin_,0.0);
  phiam_->Update(alphaM_,*phinp_,(1.0-alphaM_),*phin_,0.0);

  // accelerations are not independent but rather have to be computed
  // from phinp_, phin_ and phidtn_
  ComputeTimeDerivative();

  // compute time derivative of phi at n+alpha_M
  phidtam_->Update(alphaM_,*phidtnp_,(1.0-alphaM_),*phidtn_,0.0);

  return;
}


/*----------------------------------------------------------------------*
 | compute values of therm. pressure at interm. time steps     vg 09/09 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::ComputeThermPressureIntermediateValues()
{
  // thermodynamic pressure at n+alpha_F and n+alpha_M for low-Mach-number case
  // -> required for evaluation of equation of state
  thermpressaf_ = alphaF_*thermpressnp_ + (1.0-alphaF_)*thermpressn_;
  thermpressam_ = alphaM_*thermpressnp_ + (1.0-alphaM_)*thermpressn_;

  // time derivative of thermodyn. press. at n+alpha_F for low-Mach-number case
  // -> required as right-hand-side contribution to temperature equation,
  // hence, evaluated at n+alpha_F
  thermpressdtaf_ = alphaF_*thermpressdtnp_ + (1.0-alphaF_)*thermpressdtn_;

  // time derivative of thermodyn. press. at n+alpha_M for low-Mach-number case
  // -> required for transfer to flow solver and use in continuity equation
  thermpressdtam_ = alphaM_*thermpressdtnp_ + (1.0-alphaM_)*thermpressdtn_;

  return;
}


/*----------------------------------------------------------------------*
 | set time for evaluation of Neumann boundary conditions      vg 12/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::SetTimeForNeumannEvaluation(
  ParameterList& params)
{
  params.set("total time",time_-(1-alphaF_)*dta_);
  return;
}


/*----------------------------------------------------------------------*
 | add actual Neumann loads                                             |
 | scaled with a factor resulting from time discretization     vg 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::AddNeumannToResidual()
{
  residual_->Update(genalphafac_*dta_,*neumann_loads_,1.0);
  return;
}


/*----------------------------------------------------------------------*
 | AVM3-based scale separation                                 vg 03/09 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::AVM3Separation()
{
  // time measurement: avm3
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:            + avm3");

  // AVM3 separation
  Sep_->Multiply(false,*phiaf_,*fsphiaf_);

  // set fine-scale vector
  discret_->SetState("fsphinp",fsphiaf_);

  return;
}


/*----------------------------------------------------------------------*
 | dynamic Smagorinsky model                           rasthofer  08/12 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::DynamicComputationOfCs()
{
  if (turbmodel_==INPAR::FLUID::dynamic_smagorinsky)
  {
    // perform filtering and computation of Prt
    // compute averaged values for LkMk and MkMk
    const Teuchos::RCP<const Epetra_Vector> dirichtoggle = DirichletToggle();
    DynSmag_->ApplyFilterForDynamicComputationOfPrt(convel_,phiaf_,thermpressaf_,dirichtoggle,*extraparams_);
  }

  return;
}


/*----------------------------------------------------------------------*
 | add parameters specific for time-integration scheme         vg 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::AddSpecificTimeIntegrationParameters(
  ParameterList& params)
{
  params.set("using stationary formulation",false);
  params.set("using generalized-alpha time integration",true);
  params.set("total time",time_-(1-alphaF_)*dta_);
  params.set("time factor",genalphafac_*dta_);
  params.set("alpha_F",alphaF_);

  if (scatratype_==INPAR::SCATRA::scatratype_loma)
  {
    params.set("thermodynamic pressure",thermpressaf_);
    params.set("thermodynamic pressure at n+alpha_M",thermpressam_);
    params.set("time derivative of thermodynamic pressure",thermpressdtaf_);
    discret_->SetState("phiam",phiam_);
  }

  discret_->SetState("phinp",phiaf_);
  if (not incremental_)
  {
    discret_->SetState("hist",hist_);
    discret_->SetState("phin",phin_);
  }
  else discret_->SetState("hist",phidtam_);

  return;
}


/*----------------------------------------------------------------------*
 | compute thermodynamic pressure for low-Mach-number flow     vg 12/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::ComputeThermPressure()
{
  // compute temperature at n+alpha_F
  phiaf_->Update(alphaF_,*phinp_,(1.0-alphaF_),*phin_,0.0);

  // define element parameter list
  ParameterList eleparams;

  // DO THIS BEFORE PHINP IS SET (ClearState() is called internally!!!!)
  // compute flux approximation and add it to the parameter list
  AddFluxApproxToParameterList(eleparams,INPAR::SCATRA::flux_diffusive_domain);

  // set scalar values needed by elements
  discret_->ClearState();
  discret_->SetState("phinp",phiaf_);

  // provide velocity field (export to column map necessary for parallel evaluation)
  AddMultiVectorToParameterList(eleparams,"convective velocity field",convel_);
  AddMultiVectorToParameterList(eleparams,"velocity field",vel_);

  // provide displacement field in case of ALE
  eleparams.set("isale",isale_);
  if (isale_) AddMultiVectorToParameterList(eleparams,"dispnp",dispnp_);

  // set action for elements
  eleparams.set<int>("action",SCATRA::calc_domain_and_bodyforce);
  eleparams.set<int>("scatratype",scatratype_);
  eleparams.set("total time",time_-(1-alphaF_)*dta_);

  // variables for integrals of domain and bodyforce
  Teuchos::RCP<Epetra_SerialDenseVector> scalars
    = Teuchos::rcp(new Epetra_SerialDenseVector(2));

  // evaluate domain and bodyforce integral
  discret_->EvaluateScalars(eleparams, scalars);

  // get global integral values
  double pardomint  = (*scalars)[0];
  double parbofint  = (*scalars)[1];

  // set action for elements
  eleparams.set<int>("action",SCATRA::bd_calc_loma_therm_press);

  // variables for integrals of normal velocity and diffusive flux
  double normvelint      = 0.0;
  double normdifffluxint = 0.0;
  eleparams.set("normal velocity integral",normvelint);
  eleparams.set("normal diffusive flux integral",normdifffluxint);

  // evaluate velocity-divergence and diffusive (minus sign!) flux on boundaries
  // We may use the flux-calculation condition for calculation of fluxes for
  // thermodynamic pressure, since it is usually at the same boundary.
  vector<std::string> condnames;
  condnames.push_back("ScaTraFluxCalc");
  for (unsigned int i=0; i < condnames.size(); i++)
  {
    discret_->EvaluateCondition(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,condnames[i]);
  }

  // get integral values on this proc
  normvelint      = eleparams.get<double>("normal velocity integral");
  normdifffluxint = eleparams.get<double>("normal diffusive flux integral");

  // get integral values in parallel case
  double parnormvelint      = 0.0;
  double parnormdifffluxint = 0.0;
  discret_->Comm().SumAll(&normvelint,&parnormvelint,1);
  discret_->Comm().SumAll(&normdifffluxint,&parnormdifffluxint,1);

  // clean up
  discret_->ClearState();

  // compute thermodynamic pressure (with specific heat ratio fixed to be 1.4)
  const double shr  = 1.4;
  const double divt = shr*parnormvelint/pardomint;
  const double lhs  = alphaF_*genalphafac_*dta_*divt;
  const double rhs  = genalphafac_*dta_*(shr-1.0)*(-parnormdifffluxint+parbofint)/pardomint;
  const double hist = thermpressn_
                     - (1.0 - alphaF_)*genalphafac_*dta_*divt*thermpressn_
                     + (1.0 - genalphafac_)*dta_*thermpressdtn_;
  thermpressnp_ = (rhs + hist)/(1.0 + lhs);

  // print out thermodynamic pressure
  if (myrank_ == 0)
  {
    cout << endl;
    cout << "+--------------------------------------------------------------------------------------------+" << endl;
    cout << "Data output for instationary thermodynamic pressure:" << endl;
    cout << "Velocity in-/outflow at indicated boundary: " << parnormvelint << endl;
    cout << "Diffusive flux at indicated boundary: "       << parnormdifffluxint << endl;
    cout << "Thermodynamic pressure: "                     << thermpressnp_ << endl;
    cout << "+--------------------------------------------------------------------------------------------+" << endl;
  }

  // compute time derivative of thermodynamic pressure at time step n+1
  ComputeThermPressureTimeDerivative();

  // compute values at intermediate time steps
  ComputeThermPressureIntermediateValues();

  return;
}


/*----------------------------------------------------------------------*
 | compute time derivative                                     vg 09/09 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::ComputeTimeDerivative()
{
  // time derivative of phi:
  // phidt(n+1) = (phi(n+1)-phi(n)) / (gamma*dt) + (1-(1/gamma))*phidt(n)
  const double fact1 = 1.0/(gamma_*dta_);
  const double fact2 = 1.0 - (1.0/gamma_);
  phidtnp_->Update(fact2,*phidtn_,0.0);
  phidtnp_->Update(fact1,*phinp_,-fact1,*phin_,1.0);

  // We know the first time derivative on Dirichlet boundaries
  // so we do not need an approximation of these values!
  // However, we do not want to break the linear relationship
  // as stated above. We do not want to set Dirichlet values for
  // dependent values like phidtnp_. This turned out to be inconsistent.
  // Such an inconsistency can cause different results for
  // our different Gen. Alpha formulations (linear_full <-> linear_incremental).
  // We don't want this to happen.
  // ApplyDirichletBC(time_,Teuchos::null,phidtnp_);

  return;
}


/*----------------------------------------------------------------------*
 | compute time derivative of thermodynamic pressure           vg 09/09 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::ComputeThermPressureTimeDerivative()
{
  // time derivative of thermodynamic pressure:
  // tpdt(n+1) = (tp(n+1)-tp(n)) / (gamma*dt) + (1-(1/gamma))*tpdt(n)
  const double fact1 = 1.0/(gamma_*dta_);
  const double fact2 = 1.0 - (1.0/gamma_);
  thermpressdtnp_ = fact1*(thermpressnp_-thermpressn_) + fact2*thermpressdtn_;

  return;
}


/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                             vg 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::Update()
{
  // set history variable to zero for not spoiling flux calculation
  //if (not incremental_) hist_->PutScalar(0.0);

  // compute flux vector field for later output BEFORE time shift of results
  // is performed below !!
  if (writeflux_!=INPAR::SCATRA::flux_no)
  {
    if (DoOutput() or DoBoundaryFluxStatistics())
      flux_ = CalcFlux(true);
  }

  // compute time derivative at time n+1
  ComputeTimeDerivative();

  // solution of this step becomes most recent solution of last step
  phin_->Update(1.0,*phinp_,0.0);

  // time deriv. of this step becomes most recent time derivative of
  // last step
  phidtn_->Update(1.0,*phidtnp_,0.0);

  // perform update of time-dependent electrode variables
  ElectrodeKineticsTimeUpdate();

  // potential time update of time-dependent materials
  ElementMaterialTimeUpdate();

  return;
}


/*----------------------------------------------------------------------*
 | update level set after reinitialization              rasthofer 02/10 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::UpdateReinit()
{
  // solution of this step becomes most recent solution of the last step
  phin_ ->Update(1.0,*phinp_,0.0);

  // time deriv. of this step becomes most recent time derivative of
  // last step
  //phidtn_->Update(1.0,*phidtnp_,0.0);

  // compute time derivative at time n (and n+1)
  CalcInitialPhidt();
}


/*----------------------------------------------------------------------*
 | update thermodynamic pressure at n for low-Mach-number flow vg 12/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::UpdateThermPressure()
{
  thermpressn_   = thermpressnp_;
  thermpressdtn_ = thermpressdtnp_;

  return;
}


/*----------------------------------------------------------------------*
 | write additional data required for restart                  vg 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::OutputRestart()
{
  // additional state vectors that are needed for generalized-alpha restart
  output_->WriteVector("phidtnp",phidtnp_);
  output_->WriteVector("phidtn", phidtn_);
  output_->WriteVector("phin",   phin_);

  // for elch problems with moving boundary
  if (isale_)
    output_->WriteVector("trueresidual", trueresidual_);

  // write additional restart data for galvanostatic applications
  if (IsElch(scatratype_))
  {
    if (DRT::INPUT::IntegralValue<int>(extraparams_->sublist("ELCH CONTROL"),"GALVANOSTATIC"))
    {
      // define a vector with all electrode kinetics BCs
      vector<DRT::Condition*> cond;
      discret_->GetCondition("ElectrodeKinetics",cond);

      int condid_cathode = extraparams_->sublist("ELCH CONTROL").get<int>("GSTATCONDID_CATHODE");

      vector<DRT::Condition*>::iterator fool;
      // loop through conditions and find the cathode
      for (fool=cond.begin(); fool!=cond.end(); ++fool)
      {
        DRT::Condition* mycond = (*(fool));
        const int condid = mycond->GetInt("ConditionID");
        if (condid_cathode==condid)
        {
          // electrode potential of the adjusted electrode kinetics BC at time n+1
          double pot = mycond->GetDouble("pot");
          output_->WriteDouble("pot",pot);

          // electrode potential of the adjusted electrode kinetics BC at time n
          double potn = mycond->GetDouble("potn");
          output_->WriteDouble("potn",potn);
        }
      }
    }
  }

  // write additional restart data for loma
  // required for restart of closed systems
  if (scatratype_ == INPAR::SCATRA::scatratype_loma)
  {
    // thermodynamic pressure at time n+1
    output_->WriteDouble("thermpressnp",thermpressnp_);
    // thermodynamic pressure at time n
    output_->WriteDouble("thermpressn",thermpressn_);
    // thermodynamic pressure at time n+alpha_f
    output_->WriteDouble("thermpressaf",thermpressaf_);
    // thermodynamic pressure at time n+alpha_m
    output_->WriteDouble("thermpressam",thermpressam_);
    // time derivative of thermodynamic pressure at time n+1
    output_->WriteDouble("thermpressdtnp",thermpressdtnp_);
    // time derivative of thermodynamic pressure at time n
    output_->WriteDouble("thermpressdtn",thermpressdtn_);
    // time derivative of thermodynamic pressure at time n+alpha_f
    output_->WriteDouble("thermpressdtaf",thermpressdtaf_);
    // time derivative of thermodynamic pressure at time n+alpha_m
    output_->WriteDouble("thermpressdtam",thermpressdtam_);
    // as well as initial mass
    output_->WriteDouble("initialmass",initialmass_);
  }

  return;
}


/*----------------------------------------------------------------------*
 |                                                             vg 11/08 |
 -----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::ReadRestart(int step)
{
  IO::DiscretizationReader reader(discret_,step);
  time_ = reader.ReadDouble("time");
  step_ = reader.ReadInt("step");

  if (myrank_==0)
    cout<<"Reading ScaTra restart data (time="<<time_<<" ; step="<<step_<<")"<<endl;

  // read state vectors that are needed for generalized-alpha restart
  reader.ReadVector(phinp_,  "phinp");
  reader.ReadVector(phin_,   "phin");
  reader.ReadVector(phidtnp_,"phidtnp");
  reader.ReadVector(phidtn_, "phidtn");

  // for elch problems with moving boundary
  if(isale_)
    reader.ReadVector(trueresidual_, "trueresidual");

  // restart for galvanostatic applications
  if (IsElch(scatratype_))
  {
    if (DRT::INPUT::IntegralValue<int>(extraparams_->sublist("ELCH CONTROL"),"GALVANOSTATIC"))
    {
      // define a vector with all electrode kinetics BCs
      vector<DRT::Condition*> cond;
      discret_->GetCondition("ElectrodeKinetics",cond);

      int condid_cathode = extraparams_->sublist("ELCH CONTROL").get<int>("GSTATCONDID_CATHODE");
      vector<DRT::Condition*>::iterator fool;
      bool read_pot=false;

      // read desired values from the .control file and add/set the value to
      // the electrode kinetics boundary condition representing the cathode
      for (fool=cond.begin(); fool!=cond.end(); ++fool)
      {
        DRT::Condition* mycond = (*(fool));
        const int condid = mycond->GetInt("ConditionID");
        if (condid_cathode==condid)
        {
          double pot = reader.ReadDouble("pot");
          mycond->Add("pot",pot);

          double potn = reader.ReadDouble("potn");
          mycond->Add("potn",potn);

          read_pot=true;
          if (myrank_==0)
            cout<<"Successfully read restart data for galvanostatic mode (condid "<<condid<<")"<<endl;
        }
      }
      if (!read_pot)
        dserror("Reading of electrode potential for restart not successful.");
    }
  }

  // restart data of loma problems
  // required for restart of closed systems
  if (scatratype_ == INPAR::SCATRA::scatratype_loma)
  {
    // thermodynamic pressure at time n+1
    thermpressnp_ = reader.ReadDouble("thermpressnp");
    // thermodynamic pressure at time n
    thermpressn_ = reader.ReadDouble("thermpressn");
    // thermodynamic pressure at time n+alpha_f
    thermpressaf_ = reader.ReadDouble("thermpressaf");
    // thermodynamic pressure at time n+alpha_m
    thermpressam_ = reader.ReadDouble("thermpressam");
    // time derivative of thermodynamic pressure at time n+1
    thermpressdtnp_ = reader.ReadDouble("thermpressdtnp");
    // time derivative of thermodynamic pressure at time n
    thermpressdtn_ = reader.ReadDouble("thermpressdtn");
    // time derivative of thermodynamic pressure at time n+alpha_f
    thermpressdtaf_ = reader.ReadDouble("thermpressdtaf");
    // time derivative of thermodynamic pressure at time n+alpha_m
    thermpressdtam_ = reader.ReadDouble("thermpressdtam");
    // as well as initial mass
    initialmass_ = reader.ReadDouble("initialmass");
  }

  if (fssgd_ != INPAR::SCATRA::fssugrdiff_no or
      turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales)
    AVM3Preparation();

  return;
}


/*--------------------------------------------------------------------------------------------*
 | Redistribute the scatra discretization and vectors according to nodegraph   wichmann 02/12 |
 *--------------------------------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::Redistribute(const Teuchos::RCP<Epetra_CrsGraph> nodegraph)
{
  // let the base class do the basic redistribution and transfer of the base class members
  ScaTraTimIntImpl::Redistribute(nodegraph);

  // now do all the ost specfic steps
  const Epetra_Map* newdofrowmap = discret_->DofRowMap();
  Teuchos::RCP<Epetra_Vector> old;

  if (phiaf_ != Teuchos::null)
  {
    old = phiaf_;
    phiaf_ = LINALG::CreateVector(*newdofrowmap,true);
    LINALG::Export(*old, *phiaf_);
  }

  if (phiam_ != Teuchos::null)
  {
    old = phiam_;
    phiam_ = LINALG::CreateVector(*newdofrowmap,true);
    LINALG::Export(*old, *phiam_);
  }

  if (phidtam_ != Teuchos::null)
  {
    old = phidtam_;
    phidtam_ = LINALG::CreateVector(*newdofrowmap,true);
    LINALG::Export(*old, *phidtam_);
  }

  if (fsphiaf_ != Teuchos::null)
  {
    old = fsphiaf_;
    fsphiaf_ = LINALG::CreateVector(*newdofrowmap,true);
    LINALG::Export(*old, *fsphiaf_);
  }

  return;
}


/*----------------------------------------------------------------------*
 | Initialization procedure before the first time step         vg 11/08 |
 -----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::PrepareFirstTimeStep()
{
  // evaluate Dirichlet boundary conditions at time t=0
  // the values should match your initial field at the boundary!
  // ApplyDirichletBC(time_,phin_,phidtn_);
  ApplyDirichletBC(time_,phin_,Teuchos::null);

  // compute initial field for electric potential (ELCH)
  CalcInitialPotentialField();

  // compute time derivative of phi at time t=0
  CalcInitialPhidt();

  return;
}


/*----------------------------------------------------------------------*
 | update of time-dependent variables for electrode kinetics  gjb 11/09 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::ElectrodeKineticsTimeUpdate(const bool init)
{
  if (IsElch(scatratype_))
  {
    if (DRT::INPUT::IntegralValue<int>(extraparams_->sublist("ELCH CONTROL"),"GALVANOSTATIC"))
    {
      vector<DRT::Condition*> cond;
       discret_->GetCondition("ElectrodeKinetics",cond);
       for (size_t i=0; i < cond.size(); i++) // we update simply every condition!
       {
         double potnp = cond[i]->GetDouble("pot");
         if (init) // create and initialize additional b.c. entries if desired
         {
           cond[i]->Add("potn",potnp);
         }
         //double potn = cond[i]->GetDouble("potn");
         // shift status variables
         cond[i]->Add("potn",potnp);

         const double dlcapacitance = cond[i]->GetDouble("dlcap");
         if (dlcapacitance > EPS12)
           dserror("Galvanostatic mode for GenAlpha does not support double-layer capacitance yet.");
       }
    }
  }
  return;
}


