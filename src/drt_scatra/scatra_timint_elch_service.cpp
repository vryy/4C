/*----------------------------------------------------------------------*/
/*! \file

\brief service routines for scatra time integration for elch

\level 2


 *------------------------------------------------------------------------------------------------*/

#include "scatra_timint_elch_service.H"

#include "../drt_io/io.H"

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
SCATRA::CCCVCondition::CCCVCondition(const DRT::Condition& cccvcyclingcondition,
    std::vector<DRT::Condition*> cccvhalfcycleconditions, bool adaptivetimestepping)
    : adaptivetimesteppingonoff_(
          static_cast<bool>(cccvcyclingcondition.GetInt("AdaptiveTimeSteppingInitRelax"))),
      beginwithcharge_(static_cast<bool>(cccvcyclingcondition.GetInt("BeginWithCharging"))),
      charging_(false),
      halfcycle_charge_(Teuchos::null),
      halfcycle_discharge_(Teuchos::null),
      ihalfcycle_(-1),
      initrelaxtime_(cccvcyclingcondition.GetDouble("InitRelaxTime")),
      nhalfcycles_(cccvcyclingcondition.GetInt("NumberOfHalfCycles")),
      phasechanged_(false),
      phaseinitialrelaxation_(false),
      steplastphasechange_(-1)
{
  // safety checks
  if (adaptivetimesteppingonoff_ and !adaptivetimestepping)
    dserror(
        "Must not activate adaptive time stepping for initial relaxation while adaptive time "
        "stepping in scatra is disabled.");
  if (nhalfcycles_ < 1)
    dserror(
        "Less than one constant-current constant-voltage (CCCV) half-cycle specified in CCCV "
        "cell cycling boundary condition!");

  // loop over all conditions and create half cycles
  for (const auto& condition : cccvhalfcycleconditions)
  {
    if (condition->GetInt("ConditionID") < 0)
      dserror(
          "Constant-current constant-voltage (CCCV) half-cycle boundary condition has invalid "
          "condition ID!");

    if (condition->GetInt("ConditionID") == cccvcyclingcondition.GetInt("ConditionIDForCharge"))
      halfcycle_charge_ =
          Teuchos::rcp(new SCATRA::CCCVHalfCycleCondition(*condition, adaptivetimestepping));
    if (condition->GetInt("ConditionID") == cccvcyclingcondition.GetInt("ConditionIDForDischarge"))
      halfcycle_discharge_ =
          Teuchos::rcp(new SCATRA::CCCVHalfCycleCondition(*condition, adaptivetimestepping));
  }
  if (halfcycle_charge_ == Teuchos::null) dserror("Invalid halfcycle for charge!");
  if (halfcycle_discharge_ == Teuchos::null) dserror("Invalid halfcycle for discharge!");

  // activate first half cycle depending on initial relaxation
  if (initrelaxtime_ < 0.0)
    SetFirstCCCVHalfCycle(1);
  else if (initrelaxtime_ > 0.0 and ihalfcycle_ == -1)
    phaseinitialrelaxation_ = true;
  else
    dserror("Please choose an initial relaxation time larger than 0.0");

  return;
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void SCATRA::CCCVCondition::SetFirstCCCVHalfCycle(int step)
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
    halfcycle_charge_->ResetPhase(step);
    halfcycle_discharge_->ResetPhase(step);

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
    double cellvoltage, double cellcrate, double time) const
{
  bool phasefinished = false;

  // is the condition fulfilled to finish current phase?
  switch (GetCCCVHalfCyclePhase())
  {
    case INPAR::ELCH::CCCVHalfCyclePhase::constant_current:
    {
      phasefinished = IsExceedCellVoltage(cellvoltage);
      break;
    }
    case INPAR::ELCH::CCCVHalfCyclePhase::constant_voltage:
    {
      phasefinished = IsExceedCellCCRate(cellcrate);
      break;
    }
    case INPAR::ELCH::CCCVHalfCyclePhase::relaxation:
    {
      phasefinished = (time >= GetRelaxEndTime() - 1e-14);
      break;
    }
    default:
      dserror("illegal CC-CV operation mode");
      break;
  }

  return phasefinished;
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void SCATRA::CCCVCondition::NextPhase(int step, double time)
{
  // store, that phase is changed now.
  SetPhaseChangeObserver(step);

  // update phase and check if this half cycle is over?
  if (charging_ ? halfcycle_charge_->IsEndOfHalfCycleNextPhase(time)
                : halfcycle_discharge_->IsEndOfHalfCycleNextPhase(time))
  {
    // if half cylce is over reset all phases to constant current, flip charge mode and increase
    // counter
    charging_ ? halfcycle_discharge_->ResetPhase(step) : halfcycle_charge_->ResetPhase(step);
    charging_ = !charging_;
    ihalfcycle_++;
  }
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
bool SCATRA::CCCVCondition::IsExceedCellVoltage(double expected_cellvoltage) const
{
  return ((charging_ and expected_cellvoltage > halfcycle_charge_->GetCutOffVoltage() - 1e-14) or
          (!charging_ and expected_cellvoltage < halfcycle_discharge_->GetCutOffVoltage() + 1e-14));
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
bool SCATRA::CCCVCondition::IsExceedCellCCRate(double expected_cellcrate) const
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
bool SCATRA::CCCVCondition::IsStepsFromLastPhaseChange(int deltastep, int step)
{
  bool out = phasechanged_ ? (step >= (steplastphasechange_ + deltastep)) : false;
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
void SCATRA::CCCVCondition::SetPhaseChangeObserver(int step)
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
  if (ihalfcycle_ % 2 == 0) charging_ = !charging_;

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
    const DRT::Condition cccvhalfcyclecondition, bool adaptivetimestepping)
    : adaptivetimesteppingonoff_(
          *cccvhalfcyclecondition.Get<std::vector<int>>("AdaptiveTimeSteppingPhaseOnOff")),
      current_(0.0),
      cutoffcrate_(cccvhalfcyclecondition.GetDouble("CutoffCRate")),
      cutoffvoltage_(cccvhalfcyclecondition.GetDouble("CutoffVoltage")),
      halfcycleconditionID_(cccvhalfcyclecondition.GetInt("ConditionID")),
      phase_cccv_(INPAR::ELCH::CCCVHalfCyclePhase::undefined),
      relaxendtime_(-1.0),
      relaxtime_(cccvhalfcyclecondition.GetDouble("RelaxTime"))
{
  // safety check
  for (unsigned i = 0; i < adaptivetimesteppingonoff_.size(); ++i)
    if (adaptivetimesteppingonoff_[i] and !adaptivetimestepping)
      dserror(
          "Must not activate adaptive time stepping for half cycles while adaptive time "
          "stepping in scatra is disabled.");

  return;
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
bool SCATRA::CCCVHalfCycleCondition::IsEndOfHalfCycleNextPhase(double time)
{
  bool halfcycle_finished = false;

  // update phase and check if half cycle is over
  switch (phase_cccv_)
  {
    case INPAR::ELCH::CCCVHalfCyclePhase::constant_current:
    {
      phase_cccv_ = INPAR::ELCH::CCCVHalfCyclePhase::constant_voltage;
      break;
    }
    case INPAR::ELCH::CCCVHalfCyclePhase::constant_voltage:
    {
      phase_cccv_ = INPAR::ELCH::CCCVHalfCyclePhase::relaxation;
      relaxendtime_ = time + relaxtime_;
      break;
    }
    case INPAR::ELCH::CCCVHalfCyclePhase::relaxation:
    {
      halfcycle_finished = true;
      break;
    }
    default:
      dserror("illegal CC-CV mode");
      break;
  }

  return halfcycle_finished;
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void SCATRA::CCCVHalfCycleCondition::ResetPhase(int step)
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
      dserror("illegal CC-CV mode");
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
