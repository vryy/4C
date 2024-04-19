/*----------------------------------------------------------------------*/
/*! \file

\brief service routines for scatra time integration for elch

\level 2


 *------------------------------------------------------------------------------------------------*/

#include "baci_scatra_timint_elch_service.hpp"

#include "baci_io.hpp"

FOUR_C_NAMESPACE_OPEN

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
SCATRA::CCCVCondition::CCCVCondition(const DRT::Condition& cccvcyclingcondition,
    const std::vector<DRT::Condition*>& cccvhalfcycleconditions, const bool adaptivetimestepping,
    const int num_dofs)
    : adaptivetimesteppingonoff_(
          static_cast<bool>(*cccvcyclingcondition.Get<int>("AdaptiveTimeSteppingInitRelax"))),
      beginwithcharge_(static_cast<bool>(*cccvcyclingcondition.Get<int>("BeginWithCharging"))),
      charging_(false),
      ihalfcycle_(-1),
      initrelaxtime_(*cccvcyclingcondition.Get<double>("InitRelaxTime")),
      min_time_steps_during_init_relax_(
          *cccvcyclingcondition.Get<int>("MinTimeStepsDuringInitRelax")),
      nhalfcycles_(*cccvcyclingcondition.Get<int>("NumberOfHalfCycles")),
      num_add_adapt_timesteps_(*cccvcyclingcondition.Get<int>("NumAddAdaptTimeSteps")),
      num_dofs_(num_dofs),
      phasechanged_(false),
      phaseinitialrelaxation_(false),
      steplastphasechange_(-1)
{
  // safety checks
  if (adaptivetimesteppingonoff_ and !adaptivetimestepping)
  {
    FOUR_C_THROW(
        "Must not activate adaptive time stepping for initial relaxation while adaptive time "
        "stepping in scatra is disabled.");
  }
  if (nhalfcycles_ < 1)
  {
    FOUR_C_THROW(
        "Less than one constant-current constant-voltage (CCCV) half-cycle specified in CCCV "
        "cell cycling boundary condition!");
  }

  // loop over all conditions and create half cycles
  for (const auto& condition : cccvhalfcycleconditions)
  {
    if (*condition->Get<int>("ConditionID") < 0)
    {
      FOUR_C_THROW(
          "Constant-current constant-voltage (CCCV) half-cycle boundary condition has invalid "
          "condition ID!");
    }

    if (*condition->Get<int>("ConditionID") ==
        *cccvcyclingcondition.Get<int>("ConditionIDForCharge"))
      halfcycle_charge_ =
          Teuchos::rcp(new SCATRA::CCCVHalfCycleCondition(*condition, adaptivetimestepping));
    if (*condition->Get<int>("ConditionID") ==
        *cccvcyclingcondition.Get<int>("ConditionIDForDischarge"))
      halfcycle_discharge_ =
          Teuchos::rcp(new SCATRA::CCCVHalfCycleCondition(*condition, adaptivetimestepping));
  }
  if (halfcycle_charge_ == Teuchos::null) FOUR_C_THROW("Invalid halfcycle for charge!");
  if (halfcycle_discharge_ == Teuchos::null) FOUR_C_THROW("Invalid halfcycle for discharge!");

  // activate first half cycle depending on initial relaxation
  if (initrelaxtime_ < 0.0)
    SetFirstCCCVHalfCycle(1);
  else if (initrelaxtime_ > 0.0)
    phaseinitialrelaxation_ = true;
  else
    FOUR_C_THROW("Please choose an initial relaxation time larger than 0.0");
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void SCATRA::CCCVCondition::SetFirstCCCVHalfCycle(const int step)
{
  // only the current simulation is not restarted
  if (ihalfcycle_ == -1)
  {
    // if initial relaxation was before
    if (phaseinitialrelaxation_) SetPhaseChangeObserver(step);

    // determine whether simulation starts with charge or discharge half-cycle
    charging_ = beginwithcharge_;

    // set number of first charge or discharge half-cycle
    ihalfcycle_ = 1;

    // set phase of half cycles to constant current
    halfcycle_charge_->ResetPhase();
    halfcycle_discharge_->ResetPhase();

    // initial relaxation is over (or has never been started)
    phaseinitialrelaxation_ = false;
  }
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
INPAR::ELCH::CCCVHalfCyclePhase SCATRA::CCCVCondition::GetCCCVHalfCyclePhase() const
{
  // find current phase. If not initial relaxation look up in active half cycle
  return (phaseinitialrelaxation_ ? INPAR::ELCH::CCCVHalfCyclePhase::initital_relaxation
                                  : (charging_ ? halfcycle_charge_->GetCCCVHalfCyclePhase()
                                               : halfcycle_discharge_->GetCCCVHalfCyclePhase()));
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
bool SCATRA::CCCVCondition::IsEndOfHalfCyclePhase(
    const double cellvoltage, const double cellcrate, const double time) const
{
  bool phasefinished = false;

  // is the condition fulfilled to finish current phase?
  switch (GetCCCVHalfCyclePhase())
  {
    case INPAR::ELCH::CCCVHalfCyclePhase::constant_current:
    {
      phasefinished = ExceedCellVoltage(cellvoltage);
      break;
    }
    case INPAR::ELCH::CCCVHalfCyclePhase::constant_voltage:
    {
      phasefinished = ExceedCellCRate(cellcrate);
      break;
    }
    case INPAR::ELCH::CCCVHalfCyclePhase::relaxation:
    {
      phasefinished = (time >= GetRelaxEndTime() - 1e-14);
      break;
    }
    default:
      FOUR_C_THROW("illegal CC-CV operation mode");
      break;
  }

  return phasefinished;
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void SCATRA::CCCVCondition::NextPhase(const int step, const double time, const bool print)
{
  // store, that phase is changed now.
  SetPhaseChangeObserver(step);

  // print cycling information to improve readability of the output
  if (print)
  {
    std::cout << "\n CCCV cycling: ";
    if (charging_)
      std::cout << "Charging half-cycle: ";
    else
      std::cout << "Discharging half-cycle: ";
  }

  // update phase and check if this half cycle is over?
  if (charging_ ? halfcycle_charge_->IsEndOfHalfCycleNextPhase(time, print)
                : halfcycle_discharge_->IsEndOfHalfCycleNextPhase(time, print))
  {
    // if half cylce is over reset all phases to constant current, flip charge mode and increase
    // counter
    charging_ ? halfcycle_discharge_->ResetPhase() : halfcycle_charge_->ResetPhase();
    charging_ = !charging_;
    ihalfcycle_++;
  }
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
bool SCATRA::CCCVCondition::ExceedCellVoltage(const double expected_cellvoltage) const
{
  return ((charging_ and expected_cellvoltage > halfcycle_charge_->GetCutOffVoltage() - 1e-14) or
          (!charging_ and expected_cellvoltage < halfcycle_discharge_->GetCutOffVoltage() + 1e-14));
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
bool SCATRA::CCCVCondition::ExceedCellCRate(const double expected_cellcrate) const
{
  return ((charging_ and expected_cellcrate < (halfcycle_charge_->GetCutOffCRate() + 1e-14)) or
          (!charging_ and expected_cellcrate < (halfcycle_discharge_->GetCutOffCRate() + 1e-14)));
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
int SCATRA::CCCVCondition::GetHalfCycleConditionID() const
{
  return charging_ ? halfcycle_charge_->GetConditionID() : halfcycle_discharge_->GetConditionID();
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
bool SCATRA::CCCVCondition::ExceedMaxStepsFromLastPhaseChange(const int step)
{
  const bool out = (phasechanged_ && (step >= (steplastphasechange_ + num_add_adapt_timesteps_)));
  if (out) phasechanged_ = false;

  return out;
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
double SCATRA::CCCVCondition::GetRelaxEndTime() const
{
  return charging_ ? halfcycle_charge_->GetRelaxEndTime() : halfcycle_discharge_->GetRelaxEndTime();
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
bool SCATRA::CCCVCondition::IsAdaptiveTimeSteppingPhase() const
{
  return (phaseinitialrelaxation_
              ? static_cast<bool>(adaptivetimesteppingonoff_)
              : (charging_ ? halfcycle_charge_->IsAdaptiveTimeSteppingPhase()
                           : halfcycle_discharge_->IsAdaptiveTimeSteppingPhase()));
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void SCATRA::CCCVCondition::SetPhaseChangeObserver(const int step)
{
  // store, that phase is changed now
  phasechanged_ = true;
  steplastphasechange_ = step;
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void SCATRA::CCCVCondition::ResetPhaseChangeObserver() { phasechanged_ = false; }

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void SCATRA::CCCVCondition::ReadRestart(IO::DiscretizationReader& reader)
{
  // extract number of current charge or discharge half-cycle
  ihalfcycle_ = reader.ReadInt("ihalfcycle");

  // check whether number of current charge or discharge half-cycle is even, i.e., whether
  // current half-cycle is opposite to first one w.r.t. charge/discharge
  charging_ = ihalfcycle_ % 2 == 0 ? !beginwithcharge_ : beginwithcharge_;

  // phase was changed since last time step adaptivity
  phasechanged_ = static_cast<bool>(reader.ReadInt("phasechanged"));

  // we are in initial relaxation
  phaseinitialrelaxation_ = static_cast<bool>(reader.ReadInt("phaseinitialrelaxation"));

  // phase was lastly change at that step
  steplastphasechange_ = reader.ReadInt("steplastphasechange");

  // read restart in half cycles
  halfcycle_charge_->ReadRestart(reader);
  halfcycle_discharge_->ReadRestart(reader);
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
SCATRA::CCCVHalfCycleCondition::CCCVHalfCycleCondition(
    const DRT::Condition& cccvhalfcyclecondition, const bool adaptivetimestepping)
    : adaptivetimesteppingonoff_(
          *cccvhalfcyclecondition.Get<std::vector<int>>("AdaptiveTimeSteppingPhaseOnOff")),
      cutoffcrate_(*cccvhalfcyclecondition.Get<double>("CutoffCRate")),
      cutoffvoltage_(*cccvhalfcyclecondition.Get<double>("CutoffVoltage")),
      halfcycleconditionID_(*cccvhalfcyclecondition.Get<int>("ConditionID")),
      phase_cccv_(INPAR::ELCH::CCCVHalfCyclePhase::undefined),
      relaxendtime_(-1.0),
      relaxtime_(*cccvhalfcyclecondition.Get<double>("RelaxTime"))
{
  // safety check
  for (int i : adaptivetimesteppingonoff_)
  {
    if (i and !adaptivetimestepping)
    {
      FOUR_C_THROW(
          "Must not activate adaptive time stepping for half cycles while adaptive time "
          "stepping in scatra is disabled.");
    }
  }
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
bool SCATRA::CCCVHalfCycleCondition::IsEndOfHalfCycleNextPhase(const double time, const bool print)
{
  bool halfcycle_finished = false;

  // update phase and check if half cycle is over
  switch (phase_cccv_)
  {
    case INPAR::ELCH::CCCVHalfCyclePhase::constant_current:
    {
      phase_cccv_ = INPAR::ELCH::CCCVHalfCyclePhase::constant_voltage;
      if (print) std::cout << "CC-phase finished! \n";
      break;
    }
    case INPAR::ELCH::CCCVHalfCyclePhase::constant_voltage:
    {
      phase_cccv_ = INPAR::ELCH::CCCVHalfCyclePhase::relaxation;
      relaxendtime_ = time + relaxtime_;
      if (print) std::cout << "CV-phase finished! \n";
      break;
    }
    case INPAR::ELCH::CCCVHalfCyclePhase::relaxation:
    {
      halfcycle_finished = true;
      if (print) std::cout << "Relaxation-phase finished! \n";
      break;
    }
    default:
      FOUR_C_THROW("illegal CC-CV mode");
      break;
  }

  return halfcycle_finished;
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void SCATRA::CCCVHalfCycleCondition::ResetPhase()
{
  phase_cccv_ = INPAR::ELCH::CCCVHalfCyclePhase::constant_current;
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
bool SCATRA::CCCVHalfCycleCondition::IsAdaptiveTimeSteppingPhase() const
{
  bool adaptivetimestepping = false;

  // return, if adaptive time stepping is active for current phase
  switch (phase_cccv_)
  {
    case INPAR::ELCH::CCCVHalfCyclePhase::constant_current:
    {
      adaptivetimestepping = static_cast<bool>(adaptivetimesteppingonoff_[0]);
      break;
    }
    case INPAR::ELCH::CCCVHalfCyclePhase::constant_voltage:
    {
      adaptivetimestepping = static_cast<bool>(adaptivetimesteppingonoff_[1]);
      break;
    }
    case INPAR::ELCH::CCCVHalfCyclePhase::relaxation:
    {
      adaptivetimestepping = static_cast<bool>(adaptivetimesteppingonoff_[2]);
      break;
    }
    default:
      FOUR_C_THROW("illegal CC-CV mode");
      break;
  }

  return adaptivetimestepping;
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void SCATRA::CCCVHalfCycleCondition::ReadRestart(IO::DiscretizationReader& reader)
{
  // end time of current relaxation phase
  relaxendtime_ = reader.ReadDouble("relaxendtime");

  // current phase in half cycle
  phase_cccv_ = static_cast<INPAR::ELCH::CCCVHalfCyclePhase>(reader.ReadInt("phase_cccv"));
}

FOUR_C_NAMESPACE_CLOSE
