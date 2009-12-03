/*----------------------------------------------------------------------*/
/*!
\file scatra_timint_bdf2.cpp
\brief

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include "scatra_timint_bdf2.H"
#include <Teuchos_TimeMonitor.hpp>
#include "../drt_inpar/inpar_elch.H"


/*----------------------------------------------------------------------*
 |  Constructor (public)                                      gjb 08/08 |
 *----------------------------------------------------------------------*/
SCATRA::TimIntBDF2::TimIntBDF2(
  RCP<DRT::Discretization>      actdis,
  RCP<LINALG::Solver>           solver,
  RCP<ParameterList>            params,
  RCP<ParameterList>            extraparams,
  RCP<IO::DiscretizationWriter> output)
: ScaTraTimIntImpl(actdis,solver,params,extraparams,output),
  theta_(0.0)
{
  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  // -------------------------------------------------------------------
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // state vector for solution at time t_{n-1}
  phinm_ = LINALG::CreateVector(*dofrowmap,true);

  // ELCH with natural convection
  if (prbtype_ == "elch" &&
      (extraparams_->get<INPAR::ELCH::NatConv>("Natural Convection")!= INPAR::ELCH::natural_convection_no))
  {
    const Epetra_Map* noderowmap = discret_->NodeRowMap();

    // density at time n
    elchdensn_  = LINALG::CreateVector(*noderowmap,true);
    elchdensn_->PutScalar(1.0);

    // density at time n-1
    elchdensnm_  = LINALG::CreateVector(*noderowmap,true);
    elchdensnm_->PutScalar(1.0);
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
SCATRA::TimIntBDF2::~TimIntBDF2()
{
  return;
}


/*----------------------------------------------------------------------*
 | set part of the residual vector belonging to the old timestep        |
 |                                                            gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntBDF2::SetOldPartOfRighthandside()
{
  /*
  BDF2: for variable time step:

                 hist_ = (1+omega)^2/(1+ 2*omega) * phin_
                           - omega^2/(1+ 2*omega) * phinm_

  BDF2: for constant time step:

                 hist_ = 4/3 phin_ - 1/3 phinm_
  */
  if (step_>1)
  {
    double omega = dta_/dtp_;
    double fact1 = (1.0 + omega)*(1.0 + omega)/(1.0 + (2.0*omega));
    double fact2 = -(omega*omega)/(1+ (2.0*omega));
    hist_->Update(fact1, *phin_, fact2, *phinm_, 0.0);

    // for BDF2 theta is set by the timestepsizes, 2/3 for const. dt
    theta_ = (dta_+dtp_)/(2.0*dta_ + dtp_);
  }
  else
  {
    // for start-up of BDF2 we do one step with backward Euler
    hist_->Update(1.0, *phin_, 0.0);

    // backward Euler => use theta=1.0
    theta_=1.0;
  }

  return;
}


/*----------------------------------------------------------------------*
 | perform an explicit predictor step                         gjb 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntBDF2::ExplicitPredictor()
{
  if (step_>1) phinp_->Update(-1.0, *phinm_,2.0);
  // for step == 1 phinp_ is already correctly initialized with the
  // initial field phin_

  return;
}


/*----------------------------------------------------------------------*
 | predict thermodynamic pressure and time derivative          vg 12/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntBDF2::PredictThermPressure()
{
  // same-thermodynamic-pressure predictor (not required to be performed,
  // since we just updated the thermodynamic pressure, and thus,
  // thermpressnp_ = thermpressn_)

  // same-thermodynamic-pressure-increment predictor (currently not used)
  //if (step_>1) thermpressnp_ = 2.0*thermpressn_ - thermpressnm_;

  return;
}


/*----------------------------------------------------------------------*
 | set time for evaluation of Neumann boundary conditions      vg 12/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntBDF2::SetTimeForNeumannEvaluation(
  ParameterList& params)
{
  params.set("total time",time_);
  return;
}


/*----------------------------------------------------------------------*
 | add actual Neumann loads                                             |
 | scaled with a factor resulting from time discretization     vg 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntBDF2::AddNeumannToResidual()
{
  residual_->Update(theta_*dta_,*neumann_loads_,1.0);
  return;
}


/*----------------------------------------------------------------------*
 | AVM3-based scale separation                                 vg 03/09 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntBDF2::AVM3Separation()
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
void SCATRA::TimIntBDF2::AddSpecificTimeIntegrationParameters(
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
void SCATRA::TimIntBDF2::ComputeThermPressure()
{
  // compute "history" part (start-up of BDF2: one step backward Euler)
  double hist = 0.0;
  if (step_>1)
  {
    double omega = dta_/dtp_;
    double fact1 = (1.0 + omega)*(1.0 + omega)/(1.0 + (2.0*omega));
    double fact2 = -(omega*omega)/(1+ (2.0*omega));

    hist = fact1*thermpressn_ + fact2*thermpressnm_;
  }
  else hist = thermpressn_;

  // set scalar values needed by elements
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

  // compute time derivative of thermodynamic pressure at n+1
  thermpressdtnp_ = (thermpressnp_-thermpressn_)/dta_;

  return;
}


/*----------------------------------------------------------------------*
 | update time derivative                                      vg 09/09 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntBDF2::UpdateTimeDerivative()
{
  if (step_ == 1)
  {
    // time derivative of phi for first time step:
    // phidt(n+1) = (phi(n+1)-phi(n))/dt
    const double fact = 1.0/dta_;
    phidtnp_->Update(fact,*phinp_,-fact,*phin_,0.0);
  }
  else
  {
    // time derivative of phi:
    // phidt(n+1) = ((3/2)*phi(n+1)-2*phi(n)+(1/2)*phi(n-1))/dt
    const double fact1 = 3.0/(2.0*dta_);
    const double fact2 = -2.0/dta_;
    const double fact3 = 1.0/(2.0*dta_);
    phidtnp_->Update(fact3,*phinm_,0.0);
    phidtnp_->Update(fact1,*phinp_,fact2,*phin_,1.0);
  }

  // we know the first time derivative on Dirichlet boundaries
  // so we do not need an approximation of these values!
  ApplyDirichletBC(time_,Teuchos::null,phidtnp_);

  return;
}


/*----------------------------------------------------------------------*
 | update time derivative of thermodynamic pressure            vg 09/09 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntBDF2::UpdateThermPressureTimeDerivative()
{
  if (step_ == 1)
  {
    // time derivative of thermodynamic pressure for first time step:
    // tpdt(n+1) = (tp(n+1)-tp(n))/dt
    const double fact = 1.0/dta_;
    thermpressdtnp_ = fact*(thermpressnp_-thermpressn_);
  }
  else
  {
    // time derivative of of thermodynamic pressure:
    // tpdt(n+1) = ((3/2)*tp(n+1)-2*tp(n)+(1/2)*tp(n-1))/dt
    const double fact1 = 3.0/(2.0*dta_);
    const double fact2 = -2.0/dta_;
    const double fact3 = 1.0/(2.0*dta_);
    thermpressdtnp_ = fact1*thermpressnp_+fact2*thermpressn_+fact3*thermpressnm_;
  }

  return;
}


/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                            gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntBDF2::Update()
{
  // solution of this step becomes most recent solution of the last step
  phinm_->Update(1.0,*phin_ ,0.0);
  phin_ ->Update(1.0,*phinp_,0.0);

  ElectrodeKineticsTimeUpdate();

  return;
}


/*----------------------------------------------------------------------*
 | update thermodynamic pressure at n for low-Mach-number flow vg 12/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntBDF2::UpdateThermPressure()
{
  thermpressnm_ = thermpressn_;
  thermpressn_  = thermpressnp_;

  return;
}


/*----------------------------------------------------------------------*
 | update density at n-1 and n for ELCH natural convection    gjb 07/09 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntBDF2::UpdateDensityElch()
{
  elchdensnm_->Update(1.0,*elchdensn_ ,0.0);
  elchdensn_->Update(1.0,*elchdensnp_,0.0);

  return;
}


/*----------------------------------------------------------------------*
 | write additional data required for restart                 gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntBDF2::OutputRestart()
{
  // additional state vectors that are needed for BDF2 restart
  output_->WriteVector("phin", phin_);
  output_->WriteVector("phinm", phinm_);

  return;
}


/*----------------------------------------------------------------------*
 |                                                            gjb 08/0 |
 -----------------------------------------------------------------------*/
void SCATRA::TimIntBDF2::ReadRestart(int step)
{
  IO::DiscretizationReader reader(discret_,step);
  time_ = reader.ReadDouble("time");
  step_ = reader.ReadInt("step");

  // read state vectors that are needed for BDF2 restart
  reader.ReadVector(phinp_,"phinp");
  reader.ReadVector(phin_, "phin");
  reader.ReadVector(phinm_,"phinm");

  return;
}


/*----------------------------------------------------------------------*
 | initialization procedure before the first time step         gjb 07/08|
 *----------------------------------------------------------------------*/
void SCATRA::TimIntBDF2::PrepareFirstTimeStep()
{
  ApplyDirichletBC(time_, phin_,Teuchos::null);
  return;
}


/*----------------------------------------------------------------------*
 | update of time-dependent variables for electrode kinetics  gjb 11/09 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntBDF2::ElectrodeKineticsTimeUpdate(const bool init)
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
          cond[i]->Add("potnm",potnp);
          cond[i]->Add("pothist",0.0);
        }
        double potn = cond[i]->GetDouble("potn");
        // shift status variables
        cond[i]->Add("potnm",potn);
        cond[i]->Add("potn",potnp);

        // prepare old part of rhs for galvanostatic mode
        double pothist(0.0);
        if (step_>1)
        {
          double omega = dta_/dtp_;
          double fact1 = (1.0 + omega)*(1.0 + omega)/(1.0 + (2.0*omega));
          double fact2 = -(omega*omega)/(1+ (2.0*omega));
          pothist= fact1*potnp + fact2*potn; // potnp is the potn of the next time step, etc.
        }
        else
        {
          // for start-up of BDF2 we do one step with backward Euler
          pothist=potnp;
        }
        cond[i]->Add("pothist",pothist);
      }
    }
  }
  return;
}


#endif /* CCADISCRET */
