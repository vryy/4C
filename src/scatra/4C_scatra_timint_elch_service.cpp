/*----------------------------------------------------------------------*/
/*! \file

\brief service routines for scatra time integration for elch

\level 2


 *------------------------------------------------------------------------------------------------*/

#include "4C_scatra_timint_elch_service.hpp"

#include "4C_io.hpp"

FOUR_C_NAMESPACE_OPEN

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
SCATRA::CCCVCondition::CCCVCondition(const CORE::Conditions::Condition& cccvcyclingcondition,
    const std::vector<CORE::Conditions::Condition*>& cccvhalfcycleconditions,
    const bool adaptivetimestepping, const int num_dofs)
    : adaptivetimesteppingonoff_(static_cast<bool>(
          cccvcyclingcondition.parameters().Get<int>("AdaptiveTimeSteppingInitRelax"))),
      beginwithcharge_(
          static_cast<bool>(cccvcyclingcondition.parameters().Get<int>("BeginWithCharging"))),
      charging_(false),
      ihalfcycle_(-1),
      initrelaxtime_(cccvcyclingcondition.parameters().Get<double>("InitRelaxTime")),
      min_time_steps_during_init_relax_(
          cccvcyclingcondition.parameters().Get<int>("min_time_steps_during_init_relax")),
      nhalfcycles_(cccvcyclingcondition.parameters().Get<int>("NumberOfHalfCycles")),
      num_add_adapt_timesteps_(cccvcyclingcondition.parameters().Get<int>("NumAddAdaptTimeSteps")),
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
    if (condition->parameters().Get<int>("ConditionID") < 0)
    {
      FOUR_C_THROW(
          "Constant-current constant-voltage (CCCV) half-cycle boundary condition has invalid "
          "condition ID!");
    }

    if (condition->parameters().Get<int>("ConditionID") ==
        cccvcyclingcondition.parameters().Get<int>("ConditionIDForCharge"))
      halfcycle_charge_ =
          Teuchos::rcp(new SCATRA::CCCVHalfCycleCondition(*condition, adaptivetimestepping));
    if (condition->parameters().Get<int>("ConditionID") ==
        cccvcyclingcondition.parameters().Get<int>("ConditionIDForDischarge"))
      halfcycle_discharge_ =
          Teuchos::rcp(new SCATRA::CCCVHalfCycleCondition(*condition, adaptivetimestepping));
  }
  if (halfcycle_charge_ == Teuchos::null) FOUR_C_THROW("Invalid halfcycle for charge!");
  if (halfcycle_discharge_ == Teuchos::null) FOUR_C_THROW("Invalid halfcycle for discharge!");

  // activate first half cycle depending on initial relaxation
  if (initrelaxtime_ < 0.0)
    set_first_cccv_half_cycle(1);
  else if (initrelaxtime_ > 0.0)
    phaseinitialrelaxation_ = true;
  else
    FOUR_C_THROW("Please choose an initial relaxation time larger than 0.0");
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void SCATRA::CCCVCondition::set_first_cccv_half_cycle(const int step)
{
  // only the current simulation is not restarted
  if (ihalfcycle_ == -1)
  {
    // if initial relaxation was before
    if (phaseinitialrelaxation_) set_phase_change_observer(step);

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
INPAR::ELCH::CCCVHalfCyclePhase SCATRA::CCCVCondition::get_cccv_half_cycle_phase() const
{
  // find current phase. If not initial relaxation look up in active half cycle
  return (phaseinitialrelaxation_
              ? INPAR::ELCH::CCCVHalfCyclePhase::initital_relaxation
              : (charging_ ? halfcycle_charge_->get_cccv_half_cycle_phase()
                           : halfcycle_discharge_->get_cccv_half_cycle_phase()));
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
bool SCATRA::CCCVCondition::is_end_of_half_cycle_phase(
    const double cellvoltage, const double cellcrate, const double time) const
{
  bool phasefinished = false;

  // is the condition fulfilled to finish current phase?
  switch (get_cccv_half_cycle_phase())
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
  set_phase_change_observer(step);

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
  if (charging_ ? halfcycle_charge_->is_end_of_half_cycle_next_phase(time, print)
                : halfcycle_discharge_->is_end_of_half_cycle_next_phase(time, print))
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
int SCATRA::CCCVCondition::get_half_cycle_condition_id() const
{
  return charging_ ? halfcycle_charge_->GetConditionID() : halfcycle_discharge_->GetConditionID();
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
bool SCATRA::CCCVCondition::exceed_max_steps_from_last_phase_change(const int step)
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
bool SCATRA::CCCVCondition::is_adaptive_time_stepping_phase() const
{
  return (phaseinitialrelaxation_
              ? static_cast<bool>(adaptivetimesteppingonoff_)
              : (charging_ ? halfcycle_charge_->is_adaptive_time_stepping_phase()
                           : halfcycle_discharge_->is_adaptive_time_stepping_phase()));
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void SCATRA::CCCVCondition::set_phase_change_observer(const int step)
{
  // store, that phase is changed now
  phasechanged_ = true;
  steplastphasechange_ = step;
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void SCATRA::CCCVCondition::reset_phase_change_observer() { phasechanged_ = false; }

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void SCATRA::CCCVCondition::read_restart(CORE::IO::DiscretizationReader& reader)
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
  halfcycle_charge_->read_restart(reader);
  halfcycle_discharge_->read_restart(reader);
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
SCATRA::CCCVHalfCycleCondition::CCCVHalfCycleCondition(
    const CORE::Conditions::Condition& cccvhalfcyclecondition, const bool adaptivetimestepping)
    : adaptivetimesteppingonoff_(cccvhalfcyclecondition.parameters().Get<std::vector<int>>(
          "AdaptiveTimeSteppingPhaseOnOff")),
      cutoffcrate_(cccvhalfcyclecondition.parameters().Get<double>("CutoffCRate")),
      cutoffvoltage_(cccvhalfcyclecondition.parameters().Get<double>("CutoffVoltage")),
      halfcyclecondition_id_(cccvhalfcyclecondition.parameters().Get<int>("ConditionID")),
      phase_cccv_(INPAR::ELCH::CCCVHalfCyclePhase::undefined),
      relaxendtime_(-1.0),
      relaxtime_(cccvhalfcyclecondition.parameters().Get<double>("RelaxTime"))
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
bool SCATRA::CCCVHalfCycleCondition::is_end_of_half_cycle_next_phase(
    const double time, const bool print)
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
bool SCATRA::CCCVHalfCycleCondition::is_adaptive_time_stepping_phase() const
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
void SCATRA::CCCVHalfCycleCondition::read_restart(CORE::IO::DiscretizationReader& reader)
{
  // end time of current relaxation phase
  relaxendtime_ = reader.ReadDouble("relaxendtime");

  // current phase in half cycle
  phase_cccv_ = static_cast<INPAR::ELCH::CCCVHalfCyclePhase>(reader.ReadInt("phase_cccv"));
}

FOUR_C_NAMESPACE_CLOSE
