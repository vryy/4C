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

#ifdef CCADISCRET

#include "scatra_timint_genalpha.H"
#include <Teuchos_TimeMonitor.hpp>


/*----------------------------------------------------------------------*
 |  Constructor (public)                                       vg 11/08 |
 *----------------------------------------------------------------------*/
SCATRA::TimIntGenAlpha::TimIntGenAlpha(
  RCP<DRT::Discretization>      actdis,
  RCP<LINALG::Solver>           solver,
  RCP<ParameterList>            params,
  RCP<IO::DiscretizationWriter> output)
: ScaTraTimIntImpl(actdis,solver,params,output),
  alphaM_(params_->get<double>("alpha_M")),
  alphaF_(params_->get<double>("alpha_F")),
  gamma_ (params_->get<double>("gamma"))
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

  // temporal solution derivative at times n+1 and n
  phidtnp_ = LINALG::CreateVector(*dofrowmap,true);
  phidtn_  = LINALG::CreateVector(*dofrowmap,true);

  // only required for low-Mach-number flow
  if (prbtype_ == "loma")
  {
    // time derivative of thermodynamic pressure at n+1, n+alpha_F and n
    // (computed if not constant, otherwise remaining zero)
    thermpressdtnp_ = 0.0;
    thermpressdtaf_ = 0.0;
    thermpressdtn_  = 0.0;
  }

  // compute specific time factor for generalized-alpha time integration:
  // genalphatimefac = gamma*alpha_F/alpha_M
  if (alphaM_ < EPS12) dserror("factor alpha_M lower than or equal zero");
  genalphafac_ = gamma_/alphaM_;

  // fine-scale vector at time n+alpha_F
  if (fssgd_ != "No") fsphiaf_ = LINALG::CreateVector(*dofrowmap,true);

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
 | compute values at intermediate time steps                   vg 02/09 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::ComputeIntermediateValues()
{
  // computations for incremental case: phi at time n+alpha_F and time deriv.
  if (incremental_)
  {
    // computation of time derivative of phi:
    // phidt(n+1) = (phi(n+1)-phi(n)) / (gamma*dt) - (1/gamma -1)*phidt(n)
    const double fact1 = 1.0/(gamma_*dta_);
    const double fact2 = (-1.0/gamma_) +1.0;
    phidtnp_->Update(fact2,*phidtn_,0.0);
    phidtnp_->Update(fact1,*phinp_,-fact1,*phin_,1.0);

    // calculation of time derivative of phi at n+alpha_M, stored on
    // history vector for comfortable later transport to element routine
    hist_->Update(alphaM_,*phidtnp_,(1.0-alphaM_),*phidtn_,0.0);
  }

  // compute phi at n+alpha_F
  phiaf_->Update(alphaF_,*phinp_,(1.0-alphaF_),*phin_,0.0);

  if (prbtype_ == "loma")
  {
    // thermodynamic pressure at n+alpha_F and n+alpha_M for low-Mach-number case
    // -> required for evaluation of equation of state
    thermpressaf_ = alphaF_*thermpressnp_ + (1.0-alphaF_)*thermpressn_;
    thermpressam_ = alphaM_*thermpressnp_ + (1.0-alphaM_)*thermpressn_;

    // time derivative of thermodyn. press. at n+alpha_F for low-Mach-number case
    // -> required as right-hand-side contribution to temperature equation,
    // hence, evaluated at n+alpha_F
    thermpressdtaf_ = alphaF_*thermpressdtnp_ + (1.0-alphaF_)*thermpressdtn_;

    // compute phi at n+alpha_M
    phiam_->Update(alphaM_,*phinp_,(1.0-alphaM_),*phin_,0.0);
  }

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

  if (prbtype_ == "loma")
  {
    params.set("thermodynamic pressure",thermpressaf_);
    params.set("thermodynamic pressure at n+alpha_M",thermpressam_);
    params.set("time derivative of thermodynamic pressure",thermpressdtaf_);
    discret_->SetState("phiam",phiam_);
  }

  discret_->SetState("phinp",phiaf_);
  if (not incremental_) discret_->SetState("phin",phin_);

  return;
}


/*----------------------------------------------------------------------*
 | compute thermodynamic pressure for low-Mach-number flow     vg 12/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::ComputeThermPressure()
{
  // compute temperature at n+alpha_F
  phiaf_->Update(alphaF_,*phinp_,(1.0-alphaF_),*phin_,0.0);

  // set scalar values needed by elements
  discret_->ClearState();
  discret_->SetState("phinp",phiaf_);

  // define element parameter list
  ParameterList eleparams;

  // provide velocity field (export to column map necessary for parallel evaluation)
  AddMultiVectorToParameterList(eleparams,"velocity field",convel_);

  //provide displacement field in case of ALE
  eleparams.set("isale",isale_);
  if (isale_)
    AddMultiVectorToParameterList(eleparams,"dispnp",dispnp_);

  // set action for elements
  eleparams.set("action","calc_domain_and_bodyforce");
  eleparams.set("total time",time_-(1-alphaF_)*dta_);

  // variables for integrals of domain and bodyforce
  Teuchos::RCP<Epetra_SerialDenseVector> scalars
    = Teuchos::rcp(new Epetra_SerialDenseVector(2));

  discret_->EvaluateScalars(eleparams, scalars);

  // get global integral values
  double pardomint  = (*scalars)[0];
  double parbofint  = (*scalars)[1];

  // evaluate domain integral
  // set action for elements
  eleparams.set("action","calc_therm_press");

  // variables for integrals of velocity-divergence and diffusive flux
  double divuint = 0.0;
  double diffint = 0.0;
  eleparams.set("velocity-divergence integral",divuint);
  eleparams.set("diffusive-flux integral",     diffint);

  // evaluate velocity-divergence and rhs on boundaries
  // We may use the flux-calculation condition for calculation of fluxes for
  // thermodynamic pressure, since it is usually at the same boundary.
  vector<std::string> condnames;
  condnames.push_back("FluxCalculation");
  for (unsigned int i=0; i < condnames.size(); i++)
  {
    discret_->EvaluateCondition(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,condnames[i]);
  }

  // get integral values on this proc
  divuint = eleparams.get<double>("velocity-divergence integral");
  diffint = eleparams.get<double>("diffusive-flux integral");

  // get integral values in parallel case
  double pardivuint = 0.0;
  double pardiffint = 0.0;
  discret_->Comm().SumAll(&divuint,&pardivuint,1);
  discret_->Comm().SumAll(&diffint,&pardiffint,1);

  // clean up
  discret_->ClearState();

  // compute thermodynamic pressure (with specific heat ratio fixed to be 1.4)
  const double shr  = 1.4;
  const double divt = shr*pardivuint/pardomint;
  const double lhs  = alphaF_*genalphafac_*dta_*divt;
  const double rhs  = genalphafac_*dta_*(shr-1.0)*(pardiffint+parbofint)/pardomint;
  const double hist = thermpressn_
                     - (1.0 - alphaF_)*genalphafac_*dta_*divt*thermpressn_
                     + (1.0 - genalphafac_)*dta_*thermpressdtn_;
  thermpressnp_ = (rhs + hist)/(1.0 + lhs);

  // compute time derivative of thermodynamic pressure
  // tpdt(n+1) = (tp(n+1)-tp(n))/(gamma*dt)+((gamma-1)/gamma)*tpdt(n)
  const double fact1 = 1.0/(gamma_*dta_);
  const double fact2 = (gamma_-1.0)/gamma_;
  thermpressdtnp_ = fact1*(thermpressnp_-thermpressn_) + fact2*thermpressdtn_;

  // generalized-alpha version for time derivative of thermodynamic pressure
  /*const double lhs  = alphaF_*gamma_*dta_*divt;
  const double rhs  = (shr-1.0)*(pardiffint+parbofint)/pardomint;
  const double hist = (alphaF_-1.0)*divt*thermpressn_
                     + (alphaM_-1.0+alphaF_*(gamma_-1.0)*dta_*divt)*thermpressdtn_;
  thermpressdtnp_ = (rhs + hist)/(alphaM_ + lhs);
  const double fact1 = gamma_*dta_;
  constdouble fact2 = (1.0-gamma_)*dta_;
  thermpressnp_ = thermpressn_ + fact1*thermpressdtnp_ + fact2*thermpressdtn_;*/

  // backward-Euler version for thermodynamic pressure
  /*const double lhs  = dta_*divt;
  const double rhs  = dta_*(shr-1.0)*(pardiffint+parbofint)/pardomint;
  const double hist = thermpressn_;
  thermpressnp_ = (rhs + hist)/(1.0 + lhs);
  thermpressdtnp_ = (thermpressnp_-thermpressn_)/dta_;*/

  // print out thermodynamic pressure
  if (myrank_ == 0)
  {
    cout << endl;
    cout << "+--------------------------------------------------------------------------------------------+" << endl;
    cout << "Data output for instationary thermodynamic pressure:" << endl;
    cout << "Velocity in-/outflow at indicated boundary: " << pardivuint << endl;
    cout << "Diffusive flux at indicated boundary: "       << pardiffint << endl;
    cout << "Thermodynamic pressure: "                     << thermpressnp_ << endl;
    cout << "+--------------------------------------------------------------------------------------------+" << endl;
  }

  return;
}


/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 | and computation of acceleration for next time step          vg 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::Update()
{
  // compute time derivative of phi for non-incremental case:
  // phidt(n) = (phi(n)-phi(n-1)) / (gamma*dt(n)) - (1/gamma -1)*phidt(n-1)
  if (not incremental_)
  {
    const double fact1 = 1.0/(gamma_*dta_);
    const double fact2 = (-1.0/gamma_) +1.0;
    phidtn_->Update(fact1,*phinp_,-fact1,*phin_,fact2);

    // set history variable to zero for not spoiling flux calculation
    hist_->PutScalar(0.0);
  }
  // time deriv. of this step becomes most recent time derivative of
  // last step for incremental solver
  else phidtn_->Update(1.0,*phidtnp_,0.0);

  // solution of this step becomes most recent solution of last step
  phin_->Update(1.0,*phinp_,0.0);

  return;
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

  // read state vectors that are needed for generalized-alpha restart
  reader.ReadVector(phinp_,  "phinp");
  reader.ReadVector(phin_,   "phin");
  reader.ReadVector(phidtnp_,"phidtnp");
  reader.ReadVector(phidtn_, "phidtn");

  return;
}

/*----------------------------------------------------------------------*
 | Initialization procedure before the first time step         vg 11/08 |
 -----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::PrepareFirstTimeStep()
{
  // evaluate Dirichlet boundary conditions at time t=0
  ApplyDirichletBC(time_,phin_,phidtn_);

  // evaluate Neumann boundary conditions at time t=0
  neumann_loads_->PutScalar(0.0);
  ParameterList p;
  p.set("total time",time_);
  discret_->ClearState();
  discret_->EvaluateNeumann(p,*neumann_loads_);
  discret_->ClearState();

  // compute time derivative of phi at time t=0
  CalcInitialPhidt();

  return;
}

#endif /* CCADISCRET */
