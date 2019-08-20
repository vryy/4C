/*----------------------------------------------------------------------*/
/*! \file
\brief  connecting time-integration schemes (OST, BDF2, GenAlpha, Stationary) with
        elch-specific implementation (class ScaTraTimIntElch)
\level 2

\maintainer Christoph Schmidt

*/
/*----------------------------------------------------------------------*/
#include "scatra_timint_elch_scheme.H"
#include "scatra_timint_meshtying_strategy_base.H"

#include "../drt_io/io.H"

#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*
 |  Constructor (public)                                     ehrl 01/14 |
 *----------------------------------------------------------------------*/
SCATRA::ScaTraTimIntElchOST::ScaTraTimIntElchOST(Teuchos::RCP<DRT::Discretization> actdis,
    Teuchos::RCP<LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
    Teuchos::RCP<Teuchos::ParameterList> extraparams, Teuchos::RCP<IO::DiscretizationWriter> output)
    : ScaTraTimIntImpl(actdis, solver, sctratimintparams, extraparams, output),
      ScaTraTimIntElch(actdis, solver, params, sctratimintparams, extraparams, output),
      TimIntOneStepTheta(actdis, solver, sctratimintparams, extraparams, output)
{
  return;
}


/*----------------------------------------------------------------------*
 |  initialize time integration                              ehrl 01/14 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchOST::Init()
{
  // call Init()-functions of base classes
  // note: this order is important
  TimIntOneStepTheta::Init();
  ScaTraTimIntElch::Init();

  return;
}


/*----------------------------------------------------------------------*
 |  initialize time integration                              ehrl 01/14 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchOST::Setup()
{
  // call Setup()-functions of base classes
  // note: this order is important
  TimIntOneStepTheta::Setup();
  ScaTraTimIntElch::Setup();

  return;
}


/*----------------------------------------------------------------------*
| Destructor dtor (public)                                   ehrl 01/14 |
*-----------------------------------------------------------------------*/
SCATRA::ScaTraTimIntElchOST::~ScaTraTimIntElchOST() { return; }


/*--------------------------------------------------------------------------*
 | calculate initial electric potential field                    fang 09/15 |
 *--------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchOST::CalcInitialPotentialField()
{
  // evaluate Dirichlet boundary conditions at time t=0
  // the values should match your initial field at the boundary!
  ApplyDirichletBC(time_, phin_, Teuchos::null);
  ApplyDirichletBC(time_, phinp_, Teuchos::null);
  ComputeIntermediateValues();

  // evaluate Neumann boundary conditions at time t = 0
  ApplyNeumannBC(neumann_loads_);

  // standard general element parameters without stabilization
  SetElementGeneralParameters(true);

  // we also have to modify the time-parameter list (incremental solve)
  // actually we do not need a time integration scheme for calculating the initial electric
  // potential field, but the rhs of the standard element routine is used as starting point for this
  // special system of equations. Therefore, the rhs vector has to be scaled correctly.
  SetElementTimeParameter(true);

  // deactivate turbulence settings
  SetElementTurbulenceParameters(true);

  // call core algorithm
  ScaTraTimIntElch::CalcInitialPotentialField();

  // and finally undo our temporary settings
  SetElementGeneralParameters(false);
  SetElementTimeParameter(false);
  SetElementTurbulenceParameters(false);

  return;
}


/*----------------------------------------------------------------------*
 | write additional data required for restart                 gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchOST::OutputRestart() const
{
  // output restart information associated with one-step-theta time integration scheme
  TimIntOneStepTheta::OutputRestart();

  // output restart information associated with electrochemistry
  ScaTraTimIntElch::OutputRestart();

  // write additional restart data for galvanostatic applications or simulations including a double
  // layer formulation
  if (DRT::INPUT::IntegralValue<int>(*elchparams_, "GALVANOSTATIC") or dlcapexists_ == true)
  {
    // define a vector with all electrode kinetics BCs
    std::vector<DRT::Condition*> cond;
    discret_->GetCondition("ElchBoundaryKinetics", cond);
    if (!cond.size()) discret_->GetCondition("ElchBoundaryKineticsPoint", cond);

    int condid_cathode = elchparams_->get<int>("GSTATCONDID_CATHODE");

    std::vector<DRT::Condition*>::iterator fool;
    // loop through conditions and find the cathode
    for (fool = cond.begin(); fool != cond.end(); ++fool)
    {
      DRT::Condition* mycond = (*(fool));
      const int condid = mycond->GetInt("ConditionID");
      // galvanostatic mode: only applied potential of cathode is adapted
      if (condid_cathode == condid or dlcapexists_ == true)
      {
        std::stringstream temp;
        temp << condid;

        // electrode potential of the adjusted electrode kinetics BC at time n+1
        double pot = mycond->GetDouble("pot");
        output_->WriteDouble("pot_" + temp.str(), pot);

        // electrode potential of the adjusted electrode kinetics BC at time n
        double pot0n = mycond->GetDouble("pot0n");
        output_->WriteDouble("pot0n_" + temp.str(), pot0n);

        // electrode potential time derivative of the adjusted electrode kinetics BC at time n
        double pot0dtn = mycond->GetDouble("pot0dtn");
        output_->WriteDouble("pot0dtn_" + temp.str(), pot0dtn);

        // history of electrode potential of the adjusted electrode kinetics BC
        double pothist = mycond->GetDouble("pot0hist");
        output_->WriteDouble("pot0hist_" + temp.str(), pothist);
      }
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 |                                                            gjb 08/08 |
 -----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchOST::ReadRestart(const int step, Teuchos::RCP<IO::InputControl> input)
{
  TimIntOneStepTheta::ReadRestart(step, input);

  Teuchos::RCP<IO::DiscretizationReader> reader(Teuchos::null);
  if (input == Teuchos::null)
    reader = Teuchos::rcp(new IO::DiscretizationReader(discret_, step));
  else
    reader = Teuchos::rcp(new IO::DiscretizationReader(discret_, input, step));

  // Initialize Nernst-BC
  InitNernstBC();

  if (DRT::INPUT::IntegralValue<int>(*elchparams_, "GALVANOSTATIC") or dlcapexists_ == true)
  {
    // define a vector with all electrode kinetics BCs
    std::vector<DRT::Condition*> cond;
    discret_->GetCondition("ElchBoundaryKinetics", cond);
    if (!cond.size()) discret_->GetCondition("ElchBoundaryKineticsPoint", cond);

    int condid_cathode = elchparams_->get<int>("GSTATCONDID_CATHODE");
    std::vector<DRT::Condition*>::iterator fool;
    bool read_pot = false;

    // read desired values from the .control file and add/set the value to
    // the electrode kinetics boundary condition representing the cathode
    for (fool = cond.begin(); fool != cond.end(); ++fool)
    {
      DRT::Condition* mycond = (*(fool));
      const int condid = mycond->GetInt("ConditionID");
      // galvanostatic mode: only applied potential of cathode is adapted
      if (condid_cathode == condid or dlcapexists_ == true)
      {
        std::stringstream temp;
        temp << condid;

        double pot = reader->ReadDouble("pot_" + temp.str());
        mycond->Add("pot", pot);
        double pot0n = reader->ReadDouble("pot0n_" + temp.str());
        mycond->Add("pot0n", pot0n);
        double pot0hist = reader->ReadDouble("pot0hist_" + temp.str());
        mycond->Add("pot0hist", pot0hist);
        double pot0dtn = reader->ReadDouble("pot0dtn_" + temp.str());
        mycond->Add("pot0dtn", pot0dtn);
        read_pot = true;
        if (myrank_ == 0)
          std::cout << "Successfully read restart data for galvanostatic mode (condid " << condid
                    << ")" << std::endl;
      }
    }
    if (!read_pot) dserror("Reading of electrode potential for restart not successful.");
  }

  return;
}


/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                            gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchOST::Update(const int num)
{
  TimIntOneStepTheta::Update(num);
  ScaTraTimIntElch::Update(num);

  return;
}


/*----------------------------------------------------------------------*
 | update of time-dependent variables for electrode kinetics  gjb 11/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchOST::ElectrodeKineticsTimeUpdate()
{
  if ((DRT::INPUT::IntegralValue<int>(*elchparams_, "GALVANOSTATIC")) or dlcapexists_ == true)
  {
    ComputeTimeDerivPot0(false);

    std::vector<DRT::Condition*> cond;
    discret_->GetCondition("ElchBoundaryKinetics", cond);
    if (!cond.size()) discret_->GetCondition("ElchBoundaryKineticsPoint", cond);
    for (size_t i = 0; i < cond.size(); i++)  // we update simply every condition!
    {
      {
        double pot0np = cond[i]->GetDouble("pot");
        cond[i]->Add("pot0n", pot0np);

        double pot0dtnp = cond[i]->GetDouble("pot0dtnp");
        cond[i]->Add("pot0dtn", pot0dtnp);
      }
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | explicit predictor for nonlinear solver                   fang 10/15 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchOST::ExplicitPredictor() const
{
  // call base class routine
  TimIntOneStepTheta::ExplicitPredictor();

  // for the electric potential we just use the old values from the previous time step
  splitter_->InsertCondVector(splitter_->ExtractCondVector(phin_), phinp_);

  return;
}


/*-------------------------------------------------------------------------------------*
 | compute time derivative of applied electrode potential                   ehrl 08/13 |
 *-------------------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchOST::ComputeTimeDerivPot0(const bool init)
{
  std::vector<DRT::Condition*> cond;
  discret_->GetCondition("ElchBoundaryKinetics", cond);
  if (!cond.size()) discret_->GetCondition("ElchBoundaryKineticsPoint", cond);
  int numcond = cond.size();

  for (int icond = 0; icond < numcond; icond++)
  {
    double pot0np = cond[icond]->GetDouble("pot");
    const int functnum = cond[icond]->GetInt("funct");
    double dlcap = cond[icond]->GetDouble("dl_spec_cap");

    if (init)
    {
      // create and initialize additional b.c. entries for galvanostatic simulations or
      // simulations including a double layer
      cond[icond]->Add("pot0n", 0.0);
      cond[icond]->Add("pot0dtnp", 0.0);
      cond[icond]->Add("pot0dtn", 0.0);
      cond[icond]->Add("pot0hist", 0.0);

      if (dlcap != 0.0) dlcapexists_ = true;
    }
    else
    {
      // compute time derivative of applied potential
      if (functnum >= 0)
      {
        const double functfac = problem_->Funct(functnum).EvaluateTime(time_);

        // adjust potential at metal side accordingly
        pot0np *= functfac;
      }
      // compute time derivative of applied potential pot0
      // pot0dt(n+1) = (pot0(n+1)-pot0(n)) / (theta*dt) + (1-(1/theta))*pot0dt(n)
      double pot0n = cond[icond]->GetDouble("pot0n");
      double pot0dtn = cond[icond]->GetDouble("pot0dtn");
      double pot0dtnp = (pot0np - pot0n) / (dta_ * theta_) + (1 - (1 / theta_)) * pot0dtn;
      // add time derivative of applied potential pot0dtnp to BC
      cond[icond]->Add("pot0dtnp", pot0dtnp);
    }
  }

  return;
}


/*--------------------------------------------------------------------------*
 | set part of residual vector belonging to previous time step   fang 02/15 |
 *--------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchOST::SetOldPartOfRighthandside()
{
  // call base class routine
  TimIntOneStepTheta::SetOldPartOfRighthandside();

  // contribution from galvanostatic equation
  if ((DRT::INPUT::IntegralValue<int>(*elchparams_, "GALVANOSTATIC")) or dlcapexists_ == true)
  {
    std::vector<DRT::Condition*> cond;
    discret_->GetCondition("ElchBoundaryKinetics", cond);
    if (!cond.size()) discret_->GetCondition("ElchBoundaryKineticsPoint", cond);
    for (size_t i = 0; i < cond.size(); i++)  // we update simply every condition!
    {
      // prepare "old part of rhs" for galvanostatic equation (to be used at this time step)
      {
        // re-read values (just to be really sure no mix-up occurs)
        double pot0n = cond[i]->GetDouble("pot0n");
        double pot0dtn = cond[i]->GetDouble("pot0dtn");
        // prepare old part of rhs for galvanostatic mode
        double pothist = pot0n + (1.0 - theta_) * dta_ * pot0dtn;
        cond[i]->Add("pot0hist", pothist);
      }
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Constructor (public)                                     ehrl 01/14 |
 *----------------------------------------------------------------------*/
SCATRA::ScaTraTimIntElchBDF2::ScaTraTimIntElchBDF2(Teuchos::RCP<DRT::Discretization> actdis,
    Teuchos::RCP<LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
    Teuchos::RCP<Teuchos::ParameterList> extraparams, Teuchos::RCP<IO::DiscretizationWriter> output)
    : ScaTraTimIntImpl(actdis, solver, sctratimintparams, extraparams, output),
      ScaTraTimIntElch(actdis, solver, params, sctratimintparams, extraparams, output),
      TimIntBDF2(actdis, solver, sctratimintparams, extraparams, output)
{
  return;
}


/*----------------------------------------------------------------------*
 |  initialize time integration                              ehrl 01/14 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchBDF2::Init()
{
  // call Init()-functions of base classes
  // note: this order is important
  TimIntBDF2::Init();
  ScaTraTimIntElch::Init();

  return;
}

/*----------------------------------------------------------------------*
 |  initialize time integration                             rauch 09/16 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchBDF2::Setup()
{
  // call Setup()-functions of base classes
  // note: this order is important
  TimIntBDF2::Setup();
  ScaTraTimIntElch::Setup();

  return;
}


/*----------------------------------------------------------------------*
| Destructor dtor (public)                                   ehrl 01/14 |
*-----------------------------------------------------------------------*/
SCATRA::ScaTraTimIntElchBDF2::~ScaTraTimIntElchBDF2() { return; }


/*-----------------------------------------------------------------------------------*
 | calculate initial electric potential field                             fang 09/15 |
 *-----------------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchBDF2::CalcInitialPotentialField()
{
  ApplyDirichletBC(time_, phin_, Teuchos::null);
  ApplyDirichletBC(time_, phinp_, Teuchos::null);
  ApplyNeumannBC(neumann_loads_);

  // call core algorithm
  ScaTraTimIntElch::CalcInitialPotentialField();

  return;
}


/*----------------------------------------------------------------------*
 | write additional data required for restart                 gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchBDF2::OutputRestart() const
{
  // output restart information associated with BDF2 time integration scheme
  TimIntBDF2::OutputRestart();

  // output restart information associated with electrochemistry
  ScaTraTimIntElch::OutputRestart();

  // write additional restart data for galvanostatic applications
  if (DRT::INPUT::IntegralValue<int>(*elchparams_, "GALVANOSTATIC") or dlcapexists_ == true)
  {
    // define a vector with all electrode kinetics BCs
    std::vector<DRT::Condition*> cond;
    discret_->GetCondition("ElchBoundaryKinetics", cond);
    if (!cond.size()) discret_->GetCondition("ElchBoundaryKineticsPoint", cond);

    int condid_cathode = elchparams_->get<int>("GSTATCONDID_CATHODE");

    std::vector<DRT::Condition*>::iterator fool;
    // loop through conditions and find the cathode
    for (fool = cond.begin(); fool != cond.end(); ++fool)
    {
      DRT::Condition* mycond = (*(fool));
      const int condid = mycond->GetInt("ConditionID");
      if (condid_cathode == condid or dlcapexists_ == true)
      {
        // electrode potential of the adjusted electrode kinetics BC at time n+1
        double pot = mycond->GetDouble("pot");
        output_->WriteDouble("pot", pot);

        // electrode potential of the adjusted electrode kinetics BC at time n
        double potn = mycond->GetDouble("pot0n");
        output_->WriteDouble("pot0n", potn);

        // electrode potential of the adjusted electrode kinetics BC at time n -1
        double potnm = mycond->GetDouble("potnm");
        output_->WriteDouble("potnm", potnm);

        // history of electrode potential of the adjusted electrode kinetics BC
        double pothist = mycond->GetDouble("pothist");
        output_->WriteDouble("pothist", pothist);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |                                                            gjb 08/08 |
 -----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchBDF2::ReadRestart(const int step, Teuchos::RCP<IO::InputControl> input)
{
  TimIntBDF2::ReadRestart(step, input);

  Teuchos::RCP<IO::DiscretizationReader> reader(Teuchos::null);
  if (input == Teuchos::null)
    reader = Teuchos::rcp(new IO::DiscretizationReader(discret_, step));
  else
    reader = Teuchos::rcp(new IO::DiscretizationReader(discret_, input, step));

  // Initialize Nernst-BC
  InitNernstBC();

  // restart for galvanostatic applications
  if (DRT::INPUT::IntegralValue<int>(*elchparams_, "GALVANOSTATIC") or dlcapexists_ == true)
  {
    // define a vector with all electrode kinetics BCs
    std::vector<DRT::Condition*> cond;
    discret_->GetCondition("ElchBoundaryKinetics", cond);
    if (!cond.size()) discret_->GetCondition("ElchBoundaryKineticsPoint", cond);

    int condid_cathode = elchparams_->get<int>("GSTATCONDID_CATHODE");
    std::vector<DRT::Condition*>::iterator fool;
    bool read_pot = false;

    // read desired values from the .control file and add/set the value to
    // the electrode kinetics boundary condition representing the cathode
    for (fool = cond.begin(); fool != cond.end(); ++fool)
    {
      DRT::Condition* mycond = (*(fool));
      const int condid = mycond->GetInt("ConditionID");
      if (condid_cathode == condid or dlcapexists_ == true)
      {
        double pot = reader->ReadDouble("pot");
        mycond->Add("pot", pot);
        double potn = reader->ReadDouble("pot0n");
        mycond->Add("pot0n", potn);
        double potnm = reader->ReadDouble("potnm");
        mycond->Add("potnm", potnm);
        double pothist = reader->ReadDouble("pothist");
        mycond->Add("pothist", pothist);
        read_pot = true;
        if (myrank_ == 0)
          std::cout << "Successfully read restart data for galvanostatic mode (condid " << condid
                    << ")" << std::endl;
      }
    }
    if (!read_pot) dserror("Reading of electrode potential for restart not successful.");
  }

  return;
}

/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                            gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchBDF2::Update(const int num)
{
  TimIntBDF2::Update(num);
  ScaTraTimIntElch::Update(num);

  return;
}


/*----------------------------------------------------------------------*
 | update of time-dependent variables for electrode kinetics  gjb 11/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchBDF2::ElectrodeKineticsTimeUpdate()
{
  // The galvanostatic mode and double layer charging has never been tested if it is implemented
  // correctly!! The code have to be checked in detail, if somebody want to use it!!

  if (DRT::INPUT::IntegralValue<int>(*elchparams_, "GALVANOSTATIC") or dlcapexists_ == true)
  {
    // ComputeTimeDerivPot0(false);

    std::vector<DRT::Condition*> cond;
    discret_->GetCondition("ElchBoundaryKinetics", cond);
    if (!cond.size()) discret_->GetCondition("ElchBoundaryKineticsPoint", cond);
    for (size_t i = 0; i < cond.size(); i++)  // we update simply every condition!
    {
      {
        double potnp = cond[i]->GetDouble("pot");
        double potn = cond[i]->GetDouble("potn");
        // shift status variables
        cond[i]->Add("potnm", potn);
        cond[i]->Add("potn", potnp);
      }
    }
  }

  return;
}


/*-------------------------------------------------------------------------------------*
 | compute time derivative of applied electrode potential                   ehrl 08/13 |
 *-------------------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchBDF2::ComputeTimeDerivPot0(const bool init)
{
  std::vector<DRT::Condition*> cond;
  discret_->GetCondition("ElchBoundaryKinetics", cond);
  if (!cond.size()) discret_->GetCondition("ElchBoundaryKineticsPoint", cond);
  int numcond = cond.size();

  for (int icond = 0; icond < numcond; icond++)
  {
    double dlcap = cond[icond]->GetDouble("dl_spec_cap");

    if (init)
    {
      // create and initialize additional b.c. entries for galvanostatic simulations or
      // simulations including a double layer
      cond[icond]->Add("pot0n", 0.0);
      cond[icond]->Add("pot0nm", 0.0);
      cond[icond]->Add("pot0hist", 0.0);

      // The galvanostatic mode and double layer charging has never been tested if it is implemented
      // correctly!! The code have to be checked in detail, if somebody want to use it!!
      if (dlcap != 0.0)
        dserror(
            "Double layer charging and galvanostatic mode are not implemented for BDF2! You have "
            "to use one-step-theta time integration scheme");

      if (DRT::INPUT::IntegralValue<int>(*elchparams_, "GALVANOSTATIC") == true)
        dserror(
            "Double layer charging and galvanostatic mode are not implemented for BDF2! You have "
            "to use one-step-theta time integration scheme");
    }
  }
  return;
}


/*--------------------------------------------------------------------------*
 | set part of residual vector belonging to previous time step   fang 02/15 |
 *--------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchBDF2::SetOldPartOfRighthandside()
{
  // call base class routine
  TimIntBDF2::SetOldPartOfRighthandside();

  // contribution from galvanostatic equation
  if ((DRT::INPUT::IntegralValue<int>(*elchparams_, "GALVANOSTATIC")) or dlcapexists_ == true)
  {
    std::vector<DRT::Condition*> cond;
    discret_->GetCondition("ElchBoundaryKinetics", cond);
    if (!cond.size()) discret_->GetCondition("ElchBoundaryKineticsPoint", cond);
    for (size_t i = 0; i < cond.size(); i++)  // we update simply every condition!
                                              // prepare "old part of rhs" for galvanostatic
                                              // equation (to be used at next time step)
    {
      double pothist(0.0);
      double potn = cond[i]->GetDouble("pot0n");
      double potnm = cond[i]->GetDouble("potnm");
      if (step_ > 1)
      {
        // ?? tpdt(n+1) = ((3/2)*tp(n+1)-2*tp(n)+(1/2)*tp(n-1))/dt
        double fact1 = 4.0 / 3.0;
        double fact2 = -1.0 / 3.0;
        pothist = fact1 * potn + fact2 * potnm;
      }
      else
      {
        // for start-up of BDF2 we do one step with backward Euler
        pothist = potn;
      }
      cond[i]->Add("pothist", pothist);
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Constructor (public)                                     ehrl 01/14 |
 *----------------------------------------------------------------------*/
SCATRA::ScaTraTimIntElchGenAlpha::ScaTraTimIntElchGenAlpha(Teuchos::RCP<DRT::Discretization> actdis,
    Teuchos::RCP<LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
    Teuchos::RCP<Teuchos::ParameterList> extraparams, Teuchos::RCP<IO::DiscretizationWriter> output)
    : ScaTraTimIntImpl(actdis, solver, sctratimintparams, extraparams, output),
      ScaTraTimIntElch(actdis, solver, params, sctratimintparams, extraparams, output),
      TimIntGenAlpha(actdis, solver, sctratimintparams, extraparams, output)
{
  return;
}


/*----------------------------------------------------------------------*
 |  initialize time integration                              ehrl 01/14 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchGenAlpha::Init()
{
  // call Init()-functions of base classes
  // note: this order is important
  TimIntGenAlpha::Init();
  ScaTraTimIntElch::Init();

  return;
}


/*----------------------------------------------------------------------*
 |  initialize time integration                             rauch 09/16 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchGenAlpha::Setup()
{
  // call Setup()-functions of base classes
  // note: this order is important
  TimIntGenAlpha::Setup();
  ScaTraTimIntElch::Setup();

  return;
}


/*----------------------------------------------------------------------*
| Destructor dtor (public)                                   ehrl 01/14 |
*-----------------------------------------------------------------------*/
SCATRA::ScaTraTimIntElchGenAlpha::~ScaTraTimIntElchGenAlpha() { return; }


/*---------------------------------------------------------------------------*
 | calculate initial electric potential field                     fang 09/15 |
 ----------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchGenAlpha::CalcInitialPotentialField()
{
  // evaluate Dirichlet boundary conditions at time t = 0
  // the values should match your initial field at the boundary!
  ApplyDirichletBC(time_, phin_, Teuchos::null);
  ApplyDirichletBC(time_, phinp_, Teuchos::null);
  ComputeIntermediateValues();

  // evaluate Neumann boundary conditions at time t = 0
  ApplyNeumannBC(neumann_loads_);

  // for calculation of initial electric potential field, we have to switch off all stabilization
  // and turbulence modeling terms standard general element parameter without stabilization
  SetElementGeneralParameters(true);

  // we also have to modify the time-parameter list (incremental solve)
  // actually we do not need a time integration scheme for calculating the initial electric
  // potential field, but the rhs of the standard element routine is used as starting point for this
  // special system of equations. Therefore, the rhs vector has to be scaled correctly. Since the
  // genalpha scheme cannot be adapted easily, the backward Euler scheme is used instead.
  SetElementTimeParameterBackwardEuler();

  // deactivate turbulence settings
  SetElementTurbulenceParameters(true);

  // call core algorithm
  ScaTraTimIntElch::CalcInitialPotentialField();

  // and finally undo our temporary settings
  SetElementGeneralParameters();
  SetElementTimeParameter();
  SetElementTurbulenceParameters();

  return;
}


/*----------------------------------------------------------------------*
 | write additional data required for restart                 gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchGenAlpha::OutputRestart() const
{
  // output restart information associated with generalized-alpha time integration scheme
  TimIntGenAlpha::OutputRestart();

  // output restart information associated with electrochemistry
  ScaTraTimIntElch::OutputRestart();

  // write additional restart data for galvanostatic applications
  if (DRT::INPUT::IntegralValue<int>(*elchparams_, "GALVANOSTATIC") or dlcapexists_ == true)
  {
    // define a vector with all electrode kinetics BCs
    std::vector<DRT::Condition*> cond;
    discret_->GetCondition("ElchBoundaryKinetics", cond);
    if (!cond.size()) discret_->GetCondition("ElchBoundaryKineticsPoint", cond);

    int condid_cathode = elchparams_->get<int>("GSTATCONDID_CATHODE");

    std::vector<DRT::Condition*>::iterator fool;
    // loop through conditions and find the cathode
    for (fool = cond.begin(); fool != cond.end(); ++fool)
    {
      DRT::Condition* mycond = (*(fool));
      const int condid = mycond->GetInt("ConditionID");
      if (condid_cathode == condid or dlcapexists_ == true)
      {
        // electrode potential of the adjusted electrode kinetics BC at time n+1
        double pot = mycond->GetDouble("pot");
        output_->WriteDouble("pot", pot);

        // electrode potential of the adjusted electrode kinetics BC at time n
        double potn = mycond->GetDouble("pot0n");
        output_->WriteDouble("pot0n", potn);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |                                                            gjb 08/08 |
 -----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchGenAlpha::ReadRestart(
    const int step, Teuchos::RCP<IO::InputControl> input)
{
  TimIntGenAlpha::ReadRestart(step, input);

  Teuchos::RCP<IO::DiscretizationReader> reader(Teuchos::null);
  if (input == Teuchos::null)
    reader = Teuchos::rcp(new IO::DiscretizationReader(discret_, step));
  else
    reader = Teuchos::rcp(new IO::DiscretizationReader(discret_, input, step));

  // Initialize Nernst-BC
  InitNernstBC();

  if (DRT::INPUT::IntegralValue<int>(*elchparams_, "GALVANOSTATIC") or dlcapexists_ == true)
  {
    // define a vector with all electrode kinetics BCs
    std::vector<DRT::Condition*> cond;
    discret_->GetCondition("ElchBoundaryKinetics", cond);
    if (!cond.size()) discret_->GetCondition("ElchBoundaryKineticsPoint", cond);

    int condid_cathode = elchparams_->get<int>("GSTATCONDID_CATHODE");
    std::vector<DRT::Condition*>::iterator fool;
    bool read_pot = false;

    // read desired values from the .control file and add/set the value to
    // the electrode kinetics boundary condition representing the cathode
    for (fool = cond.begin(); fool != cond.end(); ++fool)
    {
      DRT::Condition* mycond = (*(fool));
      const int condid = mycond->GetInt("ConditionID");
      if (condid_cathode == condid or dlcapexists_ == true)
      {
        double pot = reader->ReadDouble("pot");
        mycond->Add("pot", pot);

        double potn = reader->ReadDouble("pot0n");
        mycond->Add("pot0n", potn);

        read_pot = true;
        if (myrank_ == 0)
          std::cout << "Successfully read restart data for galvanostatic mode (condid " << condid
                    << ")" << std::endl;
      }
    }
    if (!read_pot) dserror("Reading of electrode potential for restart not successful.");
  }

  return;
}

/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                            gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchGenAlpha::Update(const int num)
{
  TimIntGenAlpha::Update(num);
  ScaTraTimIntElch::Update(num);

  return;
}


/*----------------------------------------------------------------------*
 | update of time-dependent variables for electrode kinetics  gjb 11/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchGenAlpha::ElectrodeKineticsTimeUpdate()
{
  if ((DRT::INPUT::IntegralValue<int>(*elchparams_, "GALVANOSTATIC")) or dlcapexists_ == true)
  {
    ComputeTimeDerivPot0(false);

    std::vector<DRT::Condition*> cond;
    discret_->GetCondition("ElchBoundaryKinetics", cond);
    if (!cond.size()) discret_->GetCondition("ElchBoundaryKineticsPoint", cond);
    for (size_t i = 0; i < cond.size(); i++)  // we update simply every condition!
    {
      {
        double pot0np = cond[i]->GetDouble("pot");
        cond[i]->Add("pot0n", pot0np);
      }
    }
  }

  return;
}



/*-------------------------------------------------------------------------------------*
 | compute time derivative of applied electrode potential                   ehrl 08/13 |
 *-------------------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchGenAlpha::ComputeTimeDerivPot0(const bool init)
{
  std::vector<DRT::Condition*> cond;
  discret_->GetCondition("ElchBoundaryKinetics", cond);
  if (!cond.size()) discret_->GetCondition("ElchBoundaryKineticsPoint", cond);
  int numcond = cond.size();

  for (int icond = 0; icond < numcond; icond++)
  {
    double pot0np = cond[icond]->GetDouble("pot");
    const int functnum = cond[icond]->GetInt("funct");
    double dlcap = cond[icond]->GetDouble("dl_spec_cap");

    if (init)
    {
      // create and initialize additional b.c. entries for galvanostatic simulations or
      // simulations including a double layer
      cond[icond]->Add("pot0n", 0.0);
      cond[icond]->Add("pot0dtnp", 0.0);
      cond[icond]->Add("pot0dtn", 0.0);
      cond[icond]->Add("pot0hist", 0.0);

      // Double layer charging can not be integrated into the exiting framework without major
      // restructuring on element level: Derivation is based on a history vector!!
      if (dlcap != 0.0)
        dserror(
            "Double layer charging and galvanostatic mode are not implemented for generalized "
            "alpha time-integration scheme! You have to use one-step-theta time integration "
            "scheme");
      // dlcapexists_=true;

      // if(DRT::INPUT::IntegralValue<int>(*elchparams_,"GALVANOSTATIC")==true)
      //    dserror("Double layer charging and galvanostatic mode are not implemented for
      //    generalized-alpha! You have to use one-step-theta time integration scheme");
    }
    else
    {
      // these values are not used without double layer charging
      // compute time derivative of applied potential
      if (functnum >= 0)
      {
        const double functfac = problem_->Funct(functnum).EvaluateTime(time_);
        // adjust potential at metal side accordingly

        pot0np *= functfac;
      }
      // compute time derivative of applied potential pot0
      // pot0dt(n+1) = (pot0(n+1)-pot0(n)) / (theta*dt) + (1-(1/theta))*pot0dt(n)
      double pot0n = cond[icond]->GetDouble("pot0n");
      double pot0dtn = cond[icond]->GetDouble("pot0dtn");
      double pot0dtnp = (pot0np - pot0n) / (dta_ * gamma_) + (1 - (1 / gamma_)) * pot0dtn;
      // add time derivative of applied potential pot0dtnp to BC
      cond[icond]->Add("pot0dtnp", pot0dtnp);
    }
  }

  return;
}



/*----------------------------------------------------------------------*
 |  Constructor (public)                                     ehrl 01/14 |
 *----------------------------------------------------------------------*/
SCATRA::ScaTraTimIntElchStationary::ScaTraTimIntElchStationary(
    Teuchos::RCP<DRT::Discretization> actdis, Teuchos::RCP<LINALG::Solver> solver,
    Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
    Teuchos::RCP<Teuchos::ParameterList> extraparams, Teuchos::RCP<IO::DiscretizationWriter> output)
    : ScaTraTimIntImpl(actdis, solver, sctratimintparams, extraparams, output),
      ScaTraTimIntElch(actdis, solver, params, sctratimintparams, extraparams, output),
      TimIntStationary(actdis, solver, sctratimintparams, extraparams, output)
{
  return;
}


/*----------------------------------------------------------------------*
 |  initialize time integration                              ehrl 01/14 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchStationary::Init()
{
  // call Init()-functions of base classes
  // note: this order is important
  TimIntStationary::Init();
  ScaTraTimIntElch::Init();

  return;
}



/*----------------------------------------------------------------------*
 |  initialize time integration                             rauch 09/16 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchStationary::Setup()
{
  // call Setup()-functions of base classes
  // note: this order is important
  TimIntStationary::Setup();
  ScaTraTimIntElch::Setup();

  return;
}


/*----------------------------------------------------------------------*
| Destructor dtor (public)                                   ehrl 01/14 |
*-----------------------------------------------------------------------*/
SCATRA::ScaTraTimIntElchStationary::~ScaTraTimIntElchStationary() { return; }


/*---------------------------------------------------------------------------*
 | calculate initial electric potential field                     fang 09/15 |
 *---------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchStationary::CalcInitialPotentialField()
{
  ApplyDirichletBC(time_, phin_, Teuchos::null);
  ApplyDirichletBC(time_, phinp_, Teuchos::null);
  ApplyNeumannBC(neumann_loads_);

  // call core algorithm
  ScaTraTimIntElch::CalcInitialPotentialField();

  return;
}


/*----------------------------------------------------------------------*
 | write additional data required for restart                 gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchStationary::OutputRestart() const
{
  // output restart information associated with stationary time integration scheme
  TimIntStationary::OutputRestart();

  // output restart information associated with electrochemistry
  ScaTraTimIntElch::OutputRestart();

  // write additional restart data for galvanostatic applications
  if (DRT::INPUT::IntegralValue<int>(*elchparams_, "GALVANOSTATIC") or dlcapexists_ == true)
  {
    // define a vector with all electrode kinetics BCs
    std::vector<DRT::Condition*> cond;
    discret_->GetCondition("ElchBoundaryKinetics", cond);
    if (!cond.size()) discret_->GetCondition("ElchBoundaryKineticsPoint", cond);

    int condid_cathode = elchparams_->get<int>("GSTATCONDID_CATHODE");

    std::vector<DRT::Condition*>::iterator fool;
    // loop through conditions and find the cathode
    for (fool = cond.begin(); fool != cond.end(); ++fool)
    {
      DRT::Condition* mycond = (*(fool));
      const int condid = mycond->GetInt("ConditionID");
      if (condid_cathode == condid or dlcapexists_ == true)
      {
        // electrode potential of the adjusted electrode kinetics BC at time n+1
        double pot = mycond->GetDouble("pot");
        output_->WriteDouble("pot", pot);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |                                                            gjb 08/08 |
 -----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchStationary::ReadRestart(
    const int step, Teuchos::RCP<IO::InputControl> input)
{
  TimIntStationary::ReadRestart(step, input);

  Teuchos::RCP<IO::DiscretizationReader> reader(Teuchos::null);
  if (input == Teuchos::null)
    reader = Teuchos::rcp(new IO::DiscretizationReader(discret_, step));
  else
    reader = Teuchos::rcp(new IO::DiscretizationReader(discret_, input, step));

  // Initialize Nernst-BC
  InitNernstBC();

  // restart for galvanostatic applications
  if (DRT::INPUT::IntegralValue<int>(*elchparams_, "GALVANOSTATIC") or dlcapexists_ == true)
  {
    // define a vector with all electrode kinetics BCs
    std::vector<DRT::Condition*> cond;
    discret_->GetCondition("ElchBoundaryKinetics", cond);
    if (!cond.size()) discret_->GetCondition("ElchBoundaryKineticsPoint", cond);

    int condid_cathode = elchparams_->get<int>("GSTATCONDID_CATHODE");
    std::vector<DRT::Condition*>::iterator fool;
    bool read_pot = false;

    // read desired values from the .control file and add/set the value to
    // the electrode kinetics boundary condition representing the cathode
    for (fool = cond.begin(); fool != cond.end(); ++fool)
    {
      DRT::Condition* mycond = (*(fool));
      const int condid = mycond->GetInt("ConditionID");
      if (condid_cathode == condid or dlcapexists_ == true)
      {
        double pot = reader->ReadDouble("pot");
        mycond->Add("pot", pot);
        read_pot = true;
        if (myrank_ == 0)
          std::cout << "Successfully read restart data for galvanostatic mode (condid " << condid
                    << ")" << std::endl;
      }
    }
    if (!read_pot) dserror("Reading of electrode potential for restart not successful.");
  }

  return;
}

/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                            gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchStationary::Update(const int num)
{
  TimIntStationary::Update(num);

  return;
}

/*-------------------------------------------------------------------------------------*
 | compute time derivative of applied electrode potential                   ehrl 08/13 |
 *-------------------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchStationary::ComputeTimeDerivPot0(const bool init)
{
  std::vector<DRT::Condition*> cond;
  discret_->GetCondition("ElchBoundaryKinetics", cond);
  if (!cond.size()) discret_->GetCondition("ElchBoundaryKineticsPoint", cond);
  int numcond = cond.size();

  for (int icond = 0; icond < numcond; icond++)
  {
    double dlcap = cond[icond]->GetDouble("dl_spec_cap");

    if (init)
    {
      if (dlcap != 0.0)
        dserror(
            "Double layer charging and galvanostatic mode are not implemented for stationary time "
            "integration scheme! You have to use one-step-theta time integration scheme!");
      // dlcapexists_=true;

      if (DRT::INPUT::IntegralValue<int>(*elchparams_, "GALVANOSTATIC") == true)
        dserror(
            "Double layer charging and galvanostatic mode are not implemented for stationary time "
            "integration scheme! You have to use one-step-theta time integration scheme!");
    }
  }

  return;
}
