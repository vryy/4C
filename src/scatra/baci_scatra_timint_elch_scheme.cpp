/*----------------------------------------------------------------------*/
/*! \file
\brief  connecting time-integration schemes (OST, BDF2, GenAlpha, Stationary) with
        elch-specific implementation (class ScaTraTimIntElch)
\level 2


*/
/*----------------------------------------------------------------------*/
#include "baci_scatra_timint_elch_scheme.hpp"

#include "baci_global_data.hpp"
#include "baci_io.hpp"
#include "baci_lib_discret.hpp"
#include "baci_scatra_timint_meshtying_strategy_base.hpp"
#include "baci_utils_function_of_time.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  Constructor (public)                                     ehrl 01/14 |
 *----------------------------------------------------------------------*/
SCATRA::ScaTraTimIntElchOST::ScaTraTimIntElchOST(Teuchos::RCP<DRT::Discretization> actdis,
    Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
    Teuchos::RCP<Teuchos::ParameterList> extraparams, Teuchos::RCP<IO::DiscretizationWriter> output)
    : ScaTraTimIntImpl(actdis, solver, sctratimintparams, extraparams, output),
      ScaTraTimIntElch(actdis, solver, params, sctratimintparams, extraparams, output),
      TimIntOneStepTheta(actdis, solver, sctratimintparams, extraparams, output)
{
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
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchOST::PreCalcInitialPotentialField()
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
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchOST::PostCalcInitialPotentialField()
{  // and finally undo our temporary settings
  SetElementGeneralParameters(false);
  SetElementTimeParameter(false);
  SetElementTurbulenceParameters(false);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchOST::WriteRestart() const
{
  // output restart information associated with one-step-theta time integration scheme
  TimIntOneStepTheta::WriteRestart();

  // output restart information associated with electrochemistry
  ScaTraTimIntElch::WriteRestart();

  // write additional restart data for galvanostatic applications or simulations including a double
  // layer formulation
  if (CORE::UTILS::IntegralValue<int>(*elchparams_, "GALVANOSTATIC") or dlcapexists_)
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
      const int condid = *mycond->Get<int>("ConditionID");
      // galvanostatic mode: only applied potential of cathode is adapted
      if (condid_cathode == condid or dlcapexists_)
      {
        std::stringstream temp;
        temp << condid;

        // electrode potential of the adjusted electrode kinetics BC at time n+1
        auto pot = *mycond->Get<double>("pot");
        output_->WriteDouble("pot_" + temp.str(), pot);

        // electrode potential of the adjusted electrode kinetics BC at time n
        auto pot0n = *mycond->Get<double>("pot0n");
        output_->WriteDouble("pot0n_" + temp.str(), pot0n);

        // electrode potential time derivative of the adjusted electrode kinetics BC at time n
        auto pot0dtn = *mycond->Get<double>("pot0dtn");
        output_->WriteDouble("pot0dtn_" + temp.str(), pot0dtn);

        // history of electrode potential of the adjusted electrode kinetics BC
        auto pothist = *mycond->Get<double>("pot0hist");
        output_->WriteDouble("pot0hist_" + temp.str(), pothist);
      }
    }
  }
}


/*----------------------------------------------------------------------*
 |                                                            gjb 08/08 |
 -----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchOST::ReadRestart(const int step, Teuchos::RCP<IO::InputControl> input)
{
  TimIntOneStepTheta::ReadRestart(step, input);

  Teuchos::RCP<IO::DiscretizationReader> reader(Teuchos::null);
  if (input == Teuchos::null)
    reader = Teuchos::rcp(new IO::DiscretizationReader(
        discret_, GLOBAL::Problem::Instance()->InputControlFile(), step));
  else
    reader = Teuchos::rcp(new IO::DiscretizationReader(discret_, input, step));

  // Initialize Nernst-BC
  InitNernstBC();

  if (CORE::UTILS::IntegralValue<int>(*elchparams_, "GALVANOSTATIC") or dlcapexists_)
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
      const int condid = *mycond->Get<int>("ConditionID");
      // galvanostatic mode: only applied potential of cathode is adapted
      if (condid_cathode == condid or dlcapexists_)
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
    if (!read_pot) FOUR_C_THROW("Reading of electrode potential for restart not successful.");
  }
}


/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                            gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchOST::Update()
{
  TimIntOneStepTheta::Update();
  ScaTraTimIntElch::Update();
}


/*----------------------------------------------------------------------*
 | update of time-dependent variables for electrode kinetics  gjb 11/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchOST::ElectrodeKineticsTimeUpdate()
{
  if ((CORE::UTILS::IntegralValue<int>(*elchparams_, "GALVANOSTATIC")) or dlcapexists_)
  {
    ComputeTimeDerivPot0(false);

    std::vector<DRT::Condition*> conditions;
    discret_->GetCondition("ElchBoundaryKinetics", conditions);
    if (!conditions.size()) discret_->GetCondition("ElchBoundaryKineticsPoint", conditions);
    for (auto& condition : conditions)  // we update simply every condition!
    {
      {
        auto pot0np = *condition->Get<double>("pot");
        condition->Add("pot0n", pot0np);

        auto pot0dtnp = *condition->Get<double>("pot0dtnp");
        condition->Add("pot0dtn", pot0dtnp);
      }
    }
  }
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
}


/*-------------------------------------------------------------------------------------*
 | compute time derivative of applied electrode potential                   ehrl 08/13 |
 *-------------------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchOST::ComputeTimeDerivPot0(const bool init)
{
  std::vector<DRT::Condition*> cond;
  discret_->GetCondition("ElchBoundaryKinetics", cond);
  if (!cond.size()) discret_->GetCondition("ElchBoundaryKineticsPoint", cond);
  int numcond = static_cast<int>(cond.size());

  for (int icond = 0; icond < numcond; icond++)
  {
    auto pot0np = *cond[icond]->Get<double>("pot");
    const auto functnum = *cond[icond]->Get<int>("funct");
    auto dlcap = *cond[icond]->Get<double>("dl_spec_cap");

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
        const double functfac =
            problem_->FunctionById<CORE::UTILS::FunctionOfTime>(functnum).Evaluate(time_);

        // adjust potential at metal side accordingly
        pot0np *= functfac;
      }
      // compute time derivative of applied potential pot0
      // pot0dt(n+1) = (pot0(n+1)-pot0(n)) / (theta*dt) + (1-(1/theta))*pot0dt(n)
      auto pot0n = *cond[icond]->Get<double>("pot0n");
      auto pot0dtn = *cond[icond]->Get<double>("pot0dtn");
      double pot0dtnp = (pot0np - pot0n) / (dta_ * theta_) + (1 - (1 / theta_)) * pot0dtn;
      // add time derivative of applied potential pot0dtnp to BC
      cond[icond]->Add("pot0dtnp", pot0dtnp);
    }
  }
}


/*--------------------------------------------------------------------------*
 | set part of residual vector belonging to previous time step   fang 02/15 |
 *--------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchOST::SetOldPartOfRighthandside()
{
  // call base class routine
  TimIntOneStepTheta::SetOldPartOfRighthandside();

  // contribution from galvanostatic equation
  if ((CORE::UTILS::IntegralValue<int>(*elchparams_, "GALVANOSTATIC")) or dlcapexists_)
  {
    std::vector<DRT::Condition*> conditions;
    discret_->GetCondition("ElchBoundaryKinetics", conditions);
    if (!conditions.size()) discret_->GetCondition("ElchBoundaryKineticsPoint", conditions);
    for (auto& condition : conditions)  // we update simply every condition!
    {
      // prepare "old part of rhs" for galvanostatic equation (to be used at this time step)
      {
        // re-read values (just to be really sure no mix-up occurs)
        double pot0n = *condition->Get<double>("pot0n");
        double pot0dtn = *condition->Get<double>("pot0dtn");
        // prepare old part of rhs for galvanostatic mode
        double pothist = pot0n + (1.0 - theta_) * dta_ * pot0dtn;
        condition->Add("pot0hist", pothist);
      }
    }
  }
}


/*----------------------------------------------------------------------*
 |  Constructor (public)                                     ehrl 01/14 |
 *----------------------------------------------------------------------*/
SCATRA::ScaTraTimIntElchBDF2::ScaTraTimIntElchBDF2(Teuchos::RCP<DRT::Discretization> actdis,
    Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
    Teuchos::RCP<Teuchos::ParameterList> extraparams, Teuchos::RCP<IO::DiscretizationWriter> output)
    : ScaTraTimIntImpl(actdis, solver, sctratimintparams, extraparams, output),
      ScaTraTimIntElch(actdis, solver, params, sctratimintparams, extraparams, output),
      TimIntBDF2(actdis, solver, sctratimintparams, extraparams, output)
{
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
}


/*-----------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchBDF2::PreCalcInitialPotentialField()
{
  ApplyDirichletBC(time_, phin_, Teuchos::null);
  ApplyDirichletBC(time_, phinp_, Teuchos::null);
  ApplyNeumannBC(neumann_loads_);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchBDF2::WriteRestart() const
{
  // output restart information associated with BDF2 time integration scheme
  TimIntBDF2::WriteRestart();

  // output restart information associated with electrochemistry
  ScaTraTimIntElch::WriteRestart();

  // write additional restart data for galvanostatic applications
  if (CORE::UTILS::IntegralValue<int>(*elchparams_, "GALVANOSTATIC") or dlcapexists_)
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
      const int condid = *mycond->Get<int>("ConditionID");
      if (condid_cathode == condid or dlcapexists_)
      {
        // electrode potential of the adjusted electrode kinetics BC at time n+1
        auto pot = *mycond->Get<double>("pot");
        output_->WriteDouble("pot", pot);

        // electrode potential of the adjusted electrode kinetics BC at time n
        auto potn = *mycond->Get<double>("pot0n");
        output_->WriteDouble("pot0n", potn);

        // electrode potential of the adjusted electrode kinetics BC at time n -1
        auto potnm = *mycond->Get<double>("potnm");
        output_->WriteDouble("potnm", potnm);

        // history of electrode potential of the adjusted electrode kinetics BC
        auto pothist = *mycond->Get<double>("pothist");
        output_->WriteDouble("pothist", pothist);
      }
    }
  }
}

/*----------------------------------------------------------------------*
 |                                                            gjb 08/08 |
 -----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchBDF2::ReadRestart(const int step, Teuchos::RCP<IO::InputControl> input)
{
  TimIntBDF2::ReadRestart(step, input);

  Teuchos::RCP<IO::DiscretizationReader> reader(Teuchos::null);
  if (input == Teuchos::null)
    reader = Teuchos::rcp(new IO::DiscretizationReader(
        discret_, GLOBAL::Problem::Instance()->InputControlFile(), step));
  else
    reader = Teuchos::rcp(new IO::DiscretizationReader(discret_, input, step));

  // Initialize Nernst-BC
  InitNernstBC();

  // restart for galvanostatic applications
  if (CORE::UTILS::IntegralValue<int>(*elchparams_, "GALVANOSTATIC") or dlcapexists_)
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
      const int condid = *mycond->Get<int>("ConditionID");
      if (condid_cathode == condid or dlcapexists_)
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
    if (!read_pot) FOUR_C_THROW("Reading of electrode potential for restart not successful.");
  }
}

/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                            gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchBDF2::Update()
{
  TimIntBDF2::Update();
  ScaTraTimIntElch::Update();
}


/*----------------------------------------------------------------------*
 | update of time-dependent variables for electrode kinetics  gjb 11/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchBDF2::ElectrodeKineticsTimeUpdate()
{
  // The galvanostatic mode and double layer charging has never been tested if it is implemented
  // correctly!! The code have to be checked in detail, if somebody want to use it!!

  if (CORE::UTILS::IntegralValue<int>(*elchparams_, "GALVANOSTATIC") or dlcapexists_)
  {
    // ComputeTimeDerivPot0(false);

    std::vector<DRT::Condition*> conditions;
    discret_->GetCondition("ElchBoundaryKinetics", conditions);
    if (!conditions.size()) discret_->GetCondition("ElchBoundaryKineticsPoint", conditions);
    for (auto& condition : conditions)  // we update simply every condition!
    {
      {
        double potnp = *condition->Get<double>("pot");
        double potn = *condition->Get<double>("potn");
        // shift status variables
        condition->Add("potnm", potn);
        condition->Add("potn", potnp);
      }
    }
  }
}


/*-------------------------------------------------------------------------------------*
 | compute time derivative of applied electrode potential                   ehrl 08/13 |
 *-------------------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchBDF2::ComputeTimeDerivPot0(const bool init)
{
  std::vector<DRT::Condition*> cond;
  discret_->GetCondition("ElchBoundaryKinetics", cond);
  if (!cond.size()) discret_->GetCondition("ElchBoundaryKineticsPoint", cond);
  int numcond = static_cast<int>(cond.size());

  for (int icond = 0; icond < numcond; icond++)
  {
    auto dlcap = *cond[icond]->Get<double>("dl_spec_cap");

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
      {
        FOUR_C_THROW(
            "Double layer charging and galvanostatic mode are not implemented for BDF2! You have "
            "to use one-step-theta time integration scheme");
      }

      if (CORE::UTILS::IntegralValue<int>(*elchparams_, "GALVANOSTATIC") == true)
      {
        FOUR_C_THROW(
            "Double layer charging and galvanostatic mode are not implemented for BDF2! You have "
            "to use one-step-theta time integration scheme");
      }
    }
  }
}


/*--------------------------------------------------------------------------*
 | set part of residual vector belonging to previous time step   fang 02/15 |
 *--------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchBDF2::SetOldPartOfRighthandside()
{
  // call base class routine
  TimIntBDF2::SetOldPartOfRighthandside();

  // contribution from galvanostatic equation
  if ((CORE::UTILS::IntegralValue<int>(*elchparams_, "GALVANOSTATIC")) or dlcapexists_)
  {
    std::vector<DRT::Condition*> conditions;
    discret_->GetCondition("ElchBoundaryKinetics", conditions);
    if (!conditions.size()) discret_->GetCondition("ElchBoundaryKineticsPoint", conditions);
    for (auto& condition : conditions)  // we update simply every condition!
                                        // prepare "old part of rhs" for galvanostatic
                                        // equation (to be used at next time step)
    {
      double pothist;
      double potn = *condition->Get<double>("pot0n");
      double potnm = *condition->Get<double>("potnm");
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
      condition->Add("pothist", pothist);
    }
  }
}


/*----------------------------------------------------------------------*
 |  Constructor (public)                                     ehrl 01/14 |
 *----------------------------------------------------------------------*/
SCATRA::ScaTraTimIntElchGenAlpha::ScaTraTimIntElchGenAlpha(Teuchos::RCP<DRT::Discretization> actdis,
    Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
    Teuchos::RCP<Teuchos::ParameterList> extraparams, Teuchos::RCP<IO::DiscretizationWriter> output)
    : ScaTraTimIntImpl(actdis, solver, sctratimintparams, extraparams, output),
      ScaTraTimIntElch(actdis, solver, params, sctratimintparams, extraparams, output),
      TimIntGenAlpha(actdis, solver, sctratimintparams, extraparams, output)
{
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
}


/*---------------------------------------------------------------------------*
 ----------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchGenAlpha::PreCalcInitialPotentialField()
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
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchGenAlpha::PostCalcInitialPotentialField()
{
  // and finally undo our temporary settings
  SetElementGeneralParameters();
  SetElementTimeParameter();
  SetElementTurbulenceParameters();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchGenAlpha::WriteRestart() const
{
  // output restart information associated with generalized-alpha time integration scheme
  TimIntGenAlpha::WriteRestart();

  // output restart information associated with electrochemistry
  ScaTraTimIntElch::WriteRestart();

  // write additional restart data for galvanostatic applications
  if (CORE::UTILS::IntegralValue<int>(*elchparams_, "GALVANOSTATIC") or dlcapexists_)
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
      const int condid = *mycond->Get<int>("ConditionID");
      if (condid_cathode == condid or dlcapexists_)
      {
        // electrode potential of the adjusted electrode kinetics BC at time n+1
        double pot = *mycond->Get<double>("pot");
        output_->WriteDouble("pot", pot);

        // electrode potential of the adjusted electrode kinetics BC at time n
        double potn = *mycond->Get<double>("pot0n");
        output_->WriteDouble("pot0n", potn);
      }
    }
  }
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
    reader = Teuchos::rcp(new IO::DiscretizationReader(
        discret_, GLOBAL::Problem::Instance()->InputControlFile(), step));
  else
    reader = Teuchos::rcp(new IO::DiscretizationReader(discret_, input, step));

  // Initialize Nernst-BC
  InitNernstBC();

  if (CORE::UTILS::IntegralValue<int>(*elchparams_, "GALVANOSTATIC") or dlcapexists_)
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
      const int condid = *mycond->Get<int>("ConditionID");
      if (condid_cathode == condid or dlcapexists_)
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
    if (!read_pot) FOUR_C_THROW("Reading of electrode potential for restart not successful.");
  }
}

/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                            gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchGenAlpha::Update()
{
  TimIntGenAlpha::Update();
  ScaTraTimIntElch::Update();
}


/*----------------------------------------------------------------------*
 | update of time-dependent variables for electrode kinetics  gjb 11/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchGenAlpha::ElectrodeKineticsTimeUpdate()
{
  if ((CORE::UTILS::IntegralValue<int>(*elchparams_, "GALVANOSTATIC")) or dlcapexists_)
  {
    ComputeTimeDerivPot0(false);

    std::vector<DRT::Condition*> conditions;
    discret_->GetCondition("ElchBoundaryKinetics", conditions);
    if (!conditions.size()) discret_->GetCondition("ElchBoundaryKineticsPoint", conditions);
    for (auto& condition : conditions)  // we update simply every condition!
    {
      {
        auto pot0np = *condition->Get<double>("pot");
        condition->Add("pot0n", pot0np);
      }
    }
  }
}



/*-------------------------------------------------------------------------------------*
 | compute time derivative of applied electrode potential                   ehrl 08/13 |
 *-------------------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchGenAlpha::ComputeTimeDerivPot0(const bool init)
{
  std::vector<DRT::Condition*> cond;
  discret_->GetCondition("ElchBoundaryKinetics", cond);
  if (!cond.size()) discret_->GetCondition("ElchBoundaryKineticsPoint", cond);
  int numcond = static_cast<int>(cond.size());

  for (int icond = 0; icond < numcond; icond++)
  {
    double pot0np = *cond[icond]->Get<double>("pot");
    const int functnum = *cond[icond]->Get<int>("funct");
    double dlcap = *cond[icond]->Get<double>("dl_spec_cap");

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
      {
        FOUR_C_THROW(
            "Double layer charging and galvanostatic mode are not implemented for generalized "
            "alpha time-integration scheme! You have to use one-step-theta time integration "
            "scheme");
      }
    }
    else
    {
      // these values are not used without double layer charging
      // compute time derivative of applied potential
      if (functnum >= 0)
      {
        const double functfac =
            problem_->FunctionById<CORE::UTILS::FunctionOfTime>(functnum).Evaluate(time_);
        // adjust potential at metal side accordingly

        pot0np *= functfac;
      }
      // compute time derivative of applied potential pot0
      // pot0dt(n+1) = (pot0(n+1)-pot0(n)) / (theta*dt) + (1-(1/theta))*pot0dt(n)
      auto pot0n = *cond[icond]->Get<double>("pot0n");
      auto pot0dtn = *cond[icond]->Get<double>("pot0dtn");
      double pot0dtnp = (pot0np - pot0n) / (dta_ * gamma_) + (1 - (1 / gamma_)) * pot0dtn;
      // add time derivative of applied potential pot0dtnp to BC
      cond[icond]->Add("pot0dtnp", pot0dtnp);
    }
  }
}



/*----------------------------------------------------------------------*
 |  Constructor (public)                                     ehrl 01/14 |
 *----------------------------------------------------------------------*/
SCATRA::ScaTraTimIntElchStationary::ScaTraTimIntElchStationary(
    Teuchos::RCP<DRT::Discretization> actdis, Teuchos::RCP<CORE::LINALG::Solver> solver,
    Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
    Teuchos::RCP<Teuchos::ParameterList> extraparams, Teuchos::RCP<IO::DiscretizationWriter> output)
    : ScaTraTimIntImpl(actdis, solver, sctratimintparams, extraparams, output),
      ScaTraTimIntElch(actdis, solver, params, sctratimintparams, extraparams, output),
      TimIntStationary(actdis, solver, sctratimintparams, extraparams, output)
{
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
}


/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchStationary::PreCalcInitialPotentialField()
{
  ApplyDirichletBC(time_, phin_, Teuchos::null);
  ApplyDirichletBC(time_, phinp_, Teuchos::null);
  ApplyNeumannBC(neumann_loads_);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchStationary::WriteRestart() const
{
  // output restart information associated with stationary time integration scheme
  TimIntStationary::WriteRestart();

  // output restart information associated with electrochemistry
  ScaTraTimIntElch::WriteRestart();

  // write additional restart data for galvanostatic applications
  if (CORE::UTILS::IntegralValue<int>(*elchparams_, "GALVANOSTATIC") or dlcapexists_)
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
      const int condid = *mycond->Get<int>("ConditionID");
      if (condid_cathode == condid or dlcapexists_)
      {
        // electrode potential of the adjusted electrode kinetics BC at time n+1
        double pot = *mycond->Get<double>("pot");
        output_->WriteDouble("pot", pot);
      }
    }
  }
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
    reader = Teuchos::rcp(new IO::DiscretizationReader(
        discret_, GLOBAL::Problem::Instance()->InputControlFile(), step));
  else
    reader = Teuchos::rcp(new IO::DiscretizationReader(discret_, input, step));

  // Initialize Nernst-BC
  InitNernstBC();

  // restart for galvanostatic applications
  if (CORE::UTILS::IntegralValue<int>(*elchparams_, "GALVANOSTATIC") or dlcapexists_)
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
      const int condid = *mycond->Get<int>("ConditionID");
      if (condid_cathode == condid or dlcapexists_)
      {
        double pot = reader->ReadDouble("pot");
        mycond->Add("pot", pot);
        read_pot = true;
        if (myrank_ == 0)
          std::cout << "Successfully read restart data for galvanostatic mode (condid " << condid
                    << ")" << std::endl;
      }
    }
    if (!read_pot) FOUR_C_THROW("Reading of electrode potential for restart not successful.");
  }
}

/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                            gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchStationary::Update() { TimIntStationary::Update(); }

/*-------------------------------------------------------------------------------------*
 | compute time derivative of applied electrode potential                   ehrl 08/13 |
 *-------------------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchStationary::ComputeTimeDerivPot0(const bool init)
{
  std::vector<DRT::Condition*> cond;
  discret_->GetCondition("ElchBoundaryKinetics", cond);
  if (!cond.size()) discret_->GetCondition("ElchBoundaryKineticsPoint", cond);
  int numcond = static_cast<int>(cond.size());

  for (int icond = 0; icond < numcond; icond++)
  {
    auto dlcap = *cond[icond]->Get<double>("dl_spec_cap");

    if (init)
    {
      if (dlcap != 0.0)
      {
        FOUR_C_THROW(
            "Double layer charging and galvanostatic mode are not implemented for stationary time "
            "integration scheme! You have to use one-step-theta time integration scheme!");
      }

      if (CORE::UTILS::IntegralValue<int>(*elchparams_, "GALVANOSTATIC") == true)
      {
        FOUR_C_THROW(
            "Double layer charging and galvanostatic mode are not implemented for stationary time "
            "integration scheme! You have to use one-step-theta time integration scheme!");
      }
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SCATRA::ScaTraTimIntElchSCLOST::ScaTraTimIntElchSCLOST(Teuchos::RCP<DRT::Discretization> actdis,
    Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
    Teuchos::RCP<Teuchos::ParameterList> extraparams, Teuchos::RCP<IO::DiscretizationWriter> output)
    : ScaTraTimIntImpl(actdis, solver, sctratimintparams, extraparams, output),
      ScaTraTimIntElch(actdis, solver, params, sctratimintparams, extraparams, output),
      ScaTraTimIntElchSCL(actdis, solver, params, sctratimintparams, extraparams, output),
      TimIntOneStepTheta(actdis, solver, sctratimintparams, extraparams, output)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchSCLOST::Init()
{
  // call Init()-functions of base classes
  // note: this order is important
  TimIntOneStepTheta::Init();
  ScaTraTimIntElchSCL::Init();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchSCLOST::Setup()
{
  // call Setup()-functions of base classes
  // note: this order is important
  TimIntOneStepTheta::Setup();
  ScaTraTimIntElchSCL::Setup();
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchSCLOST::PreCalcInitialPotentialField()
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
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchSCLOST::PostCalcInitialPotentialField()
{  // and finally undo our temporary settings
  SetElementGeneralParameters(false);
  SetElementTimeParameter(false);
  SetElementTurbulenceParameters(false);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchSCLOST::WriteRestart() const
{
  // output restart information associated with one-step-theta time integration scheme
  TimIntOneStepTheta::WriteRestart();

  // output restart information associated with electrochemistry
  ScaTraTimIntElchSCL::WriteRestart();
}


/*----------------------------------------------------------------------*
 -----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchSCLOST::ReadRestart(
    const int step, Teuchos::RCP<IO::InputControl> input)
{
  FOUR_C_THROW("Restart is not implemented for coupled space-charge layers.");
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchSCLOST::Update()
{
  TimIntOneStepTheta::Update();
  ScaTraTimIntElchSCL::Update();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchSCLOST::ExplicitPredictor() const
{
  // call base class routine
  TimIntOneStepTheta::ExplicitPredictor();

  // for the electric potential we just use the old values from the previous time step
  splitter_->InsertCondVector(splitter_->ExtractCondVector(phin_), phinp_);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchSCLOST::SetOldPartOfRighthandside()
{
  // call base class routine
  TimIntOneStepTheta::SetOldPartOfRighthandside();

  // contribution from galvanostatic equation
  if ((CORE::UTILS::IntegralValue<int>(*elchparams_, "GALVANOSTATIC")) or dlcapexists_)
  {
    std::vector<DRT::Condition*> conditions;
    discret_->GetCondition("ElchBoundaryKinetics", conditions);
    if (!conditions.size()) discret_->GetCondition("ElchBoundaryKineticsPoint", conditions);
    for (auto& condition : conditions)  // we update simply every condition!
    {
      // prepare "old part of rhs" for galvanostatic equation (to be used at this time step)
      {
        // re-read values (just to be really sure no mix-up occurs)
        auto pot0n = *condition->Get<double>("pot0n");
        auto pot0dtn = *condition->Get<double>("pot0dtn");
        // prepare old part of rhs for galvanostatic mode
        double pothist = pot0n + (1.0 - theta_) * dta_ * pot0dtn;
        condition->Add("pot0hist", pothist);
      }
    }
  }
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchSCLOST::AddTimeIntegrationSpecificVectors(bool forcedincrementalsolver)
{
  TimIntOneStepTheta::AddTimeIntegrationSpecificVectors(forcedincrementalsolver);
  ScaTraTimIntElchSCL::AddTimeIntegrationSpecificVectors(forcedincrementalsolver);
}
FOUR_C_NAMESPACE_CLOSE
