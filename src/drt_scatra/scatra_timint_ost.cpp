/*----------------------------------------------------------------------*/
/*!
\file scatra_timint_ost.cpp
\brief One-Step-Theta time integration scheme

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include "scatra_timint_ost.H"
#include <Teuchos_TimeMonitor.hpp>
#include "../drt_inpar/inpar_elch.H"


/*----------------------------------------------------------------------*
 |  Constructor (public)                                      gjb 08/08 |
 *----------------------------------------------------------------------*/
SCATRA::TimIntOneStepTheta::TimIntOneStepTheta(
  RCP<DRT::Discretization>      actdis,
  RCP<LINALG::Solver>           solver,
  RCP<ParameterList>            params,
  RCP<ParameterList>            extraparams,
  RCP<IO::DiscretizationWriter> output)
: ScaTraTimIntImpl(actdis,solver,params,extraparams,output),
  theta_(params_->get<double>("THETA"))
{
  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  // -------------------------------------------------------------------
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // temporal solution derivative at time n
  phidtn_ = LINALG::CreateVector(*dofrowmap,true);

  // ELCH with natural convection
  if (prbtype_ == "elch" &&
      (extraparams_->get<INPAR::ELCH::NatConv>("Natural Convection")!= INPAR::ELCH::natural_convection_no))
  {
    const Epetra_Map* noderowmap = discret_->NodeRowMap();

    // density at time n
    elchdensn_ = LINALG::CreateVector(*noderowmap,true);
    elchdensn_->PutScalar(1.0);

  }

  // fine-scale vector at time n+1
  if (fssgd_ != INPAR::SCATRA::fssugrdiff_no)
    fsphinp_ = LINALG::CreateVector(*dofrowmap,true);

  // initialize time-dependent electrode kinetics variables (galvanostatic mode)
  ElectrodeKineticsTimeUpdate(true);

  return;
}


/*----------------------------------------------------------------------*
| Destructor dtor (public)                                   gjb 08/08 |
*----------------------------------------------------------------------*/
SCATRA::TimIntOneStepTheta::~TimIntOneStepTheta()
{
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

  return;
}


/*----------------------------------------------------------------------*
 | perform an explicit predictor step                         gjb 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::ExplicitPredictor()
{
  phinp_->Update(dta_, *phidtn_,1.0);

  // for the electric potential we just use the 'old' value of last time step
  Teuchos::RCP<Epetra_Vector> onlypot = splitter_.ExtractCondVector(phin_);
  splitter_.InsertCondVector(onlypot, phinp_);

  return;
}


/*----------------------------------------------------------------------*
 | predict thermodynamic pressure and time derivative          vg 12/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::PredictThermPressure()
{
  // same-thermodynamic-pressure predictor (not required to be performed,
  // since we just updated the thermodynamic pressure, and thus,
  // thermpressnp_ = thermpressn_)

  // same-thermodynamic-pressure-derivative predictor (currently not used)
  //thermpressnp_ += dta_*thermpressdtn_;

  return;
}


/*----------------------------------------------------------------------*
 | set time for evaluation of Neumann boundary conditions      vg 12/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::SetTimeForNeumannEvaluation(
  ParameterList& params)
{
  params.set("total time",time_);
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
 | add parameters specific for time-integration scheme         vg 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::AddSpecificTimeIntegrationParameters(
  ParameterList& params)
{
  params.set("using stationary formulation",false);
  params.set("using generalized-alpha time integration",false);
  params.set("total time",time_);
  params.set("time factor",theta_*dta_);
  params.set("alpha_F",1.0);

  if (prbtype_ == "loma")
  {
    params.set("thermodynamic pressure",thermpressnp_);
    params.set("time derivative of thermodynamic pressure",thermpressdtnp_);
  }

  discret_->SetState("hist",hist_);
  discret_->SetState("phinp",phinp_);

  return;
}


/*----------------------------------------------------------------------*
 | compute thermodynamic pressure for low-Mach-number flow     vg 12/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::ComputeThermPressure()
{
  // compute "history" part
  const double hist = thermpressn_ + (1.0-theta_)*dta_*thermpressdtn_;

  // set scalar and density vector values needed by elements
  discret_->ClearState();
  discret_->SetState("phinp",phinp_);

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
  eleparams.set("total time",time_);

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
  const double shr = 1.4;
  const double lhs = theta_*dta_*shr*pardivuint/pardomint;
  const double rhs = theta_*dta_*(shr-1.0)*(pardiffint+parbofint)/pardomint;
  thermpressnp_ = (rhs + hist)/(1.0 + lhs);

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
 | update time derivative                                      vg 09/09 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::UpdateTimeDerivative()
{
  // time derivative of phi:
  // phidt(n+1) = (phi(n+1)-phi(n)) / (theta*dt) + (1-(1/theta))*phidt(n)
  const double fact1 = 1.0/(theta_*dta_);
  const double fact2 = 1.0 - (1.0/theta_);
  phidtnp_->Update(fact2,*phidtn_,0.0);
  phidtnp_->Update(fact1,*phinp_,-fact1,*phin_,1.0);

  // we know the first time derivative on Dirichlet boundaries
  // so we do not need an approximation of these values!
  ApplyDirichletBC(time_,Teuchos::null,phidtnp_);

  return;
}


/*----------------------------------------------------------------------*
 | update time derivative of thermodynamic pressure            vg 09/09 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::UpdateThermPressureTimeDerivative()
{
  // time derivative of thermodynamic pressure:
  // tpdt(n+1) = (tp(n+1)-tp(n))/(theta*dt)+((theta-1)/theta)*tpdt(n)
  double fact1 = 1.0/(theta_*dta_);
  double fact2 = (theta_-1.0)/theta_;
  thermpressdtnp_ = fact1*(thermpressnp_-thermpressn_) + fact2*thermpressdtn_;

  return;
}


/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                            gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::Update()
{
  // solution of this step becomes most recent solution of the last step
  phin_ ->Update(1.0,*phinp_,0.0);

  // time deriv. of this step becomes most recent time derivative of
  // last step
  phidtn_->Update(1.0,*phidtnp_,0.0);

  ElectrodeKineticsTimeUpdate();

  return;
}


/*----------------------------------------------------------------------*
 | update thermodynamic pressure at n for low-Mach-number flow vg 12/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::UpdateThermPressure()
{
  thermpressn_   = thermpressnp_;
  thermpressdtn_ = thermpressdtnp_;

  return;
}


/*----------------------------------------------------------------------*
 | update density at n for ELCH natural convection            gjb 07/09 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::UpdateDensityElch()
{
  elchdensn_->Update(1.0,*elchdensnp_,0.0);

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

  // read state vectors that are needed for One-Step-Theta restart
  reader.ReadVector(phinp_, "phinp");
  reader.ReadVector(phin_,  "phin");
  reader.ReadVector(phidtn_,"phidtn");

  return;
}

/*----------------------------------------------------------------------*
 | Initialization procedure before the first time step        gjb 08/08 |
 -----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::PrepareFirstTimeStep()
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


/*----------------------------------------------------------------------*
 | update of time-dependent variables for electrode kinetics  gjb 11/09 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::ElectrodeKineticsTimeUpdate(const bool init)
{
  if (prbtype_ == "elch")
  {
    if (Teuchos::getIntegralValue<int>(extraparams_->sublist("ELCH CONTROL"),"GALVANOSTATIC"))
    {
      vector<DRT::Condition*> cond;
      discret_->GetCondition("ElectrodeKinetics",cond);
      for (size_t i=0; i < cond.size(); i++) // we update simply every condition!
      {
        double potnp = cond[i]->GetDouble("pot");
        if (init) // create and initialize additional b.c. entries if desired
        {
          cond[i]->Add("potn",potnp);
          cond[i]->Add("potdtn",0.0);
          cond[i]->Add("pothist",0.0);
        }
        double potn = cond[i]->GetDouble("potn");
        double potdtn = cond[i]->GetDouble("potdtn");
        // update the time derivative
        potdtn = (1/(theta_*dta_))*(potnp- potn) + (1 - (1/theta_))*dta_*potdtn;
        cond[i]->Add("potdtn",potdtn);
        // shift status variables
        cond[i]->Add("potn",potnp);
        // prepare old part of rhs for galvanostatic mode
        double pothist = potn + theta_*dta_*potdtn;
        cond[i]->Add("pothist",pothist);
      }
    }
  }
  return;
}

#endif /* CCADISCRET */
