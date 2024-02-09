/*----------------------------------------------------------------------*/
/*! \file

\brief service routines for scatra time integration for elch

\level 2


 *------------------------------------------------------------------------------------------------*/

#ifndef BACI_SCATRA_TIMINT_ELCH_SERVICE_HPP
#define BACI_SCATRA_TIMINT_ELCH_SERVICE_HPP

#include "baci_config.hpp"

#include "baci_inpar_elch.hpp"
#include "baci_lib_condition.hpp"
#include "baci_scatra_timint_implicit.hpp"

BACI_NAMESPACE_OPEN

/*==========================================================================*/
// forward declarations
/*==========================================================================*/

namespace SCATRA
{
  /*!
   * \brief Control routine for constant-current, constant-voltage half cycle
   *
   * Holds all values for one constant-current, constant-voltage half cycle, including static
   * values from input file and dynamic values. This class controls the current phase within one
   * half cycle.
   */
  class CCCVHalfCycleCondition
  {
   public:
    //! constructor
    CCCVHalfCycleCondition(const DRT::Condition& cccvhalfcyclecondition, bool adaptivetimestepping);

    //! Get phase of half cycle
    INPAR::ELCH::CCCVHalfCyclePhase GetCCCVHalfCyclePhase() const { return phase_cccv_; };

    //! Get ID of this half cycle condition
    int GetConditionID() const { return halfcycleconditionID_; };

    //! get cut off c-rate during constant voltage
    double GetCutOffCRate() const { return cutoffcrate_; };

    //! get cut off voltage during constant current
    double GetCutOffVoltage() const { return cutoffvoltage_; };

    //! get end time of current relaxation phase
    double GetRelaxEndTime() const { return relaxendtime_; };

    //! does this phase have adaptive time stepping?
    bool IsAdaptiveTimeSteppingPhase() const;

    //! is this half cycle completed after updating the phase?
    bool IsEndOfHalfCycleNextPhase(double time, bool print);

    //! reset phase of this half cycle to constant current
    void ResetPhase();

    //! read restart
    void ReadRestart(IO::DiscretizationReader& reader);

   private:
    //! adaptive time stepping at end of phases?
    std::vector<int> adaptivetimesteppingonoff_;

    //! cut off c-rate during constant voltage
    const double cutoffcrate_;

    //! ut off voltage during constant current
    const double cutoffvoltage_;

    //! ID of this half cycle condition
    const int halfcycleconditionID_;

    //! flag indicating whether cell is currently being operated in constant-current (CC),
    //! constant-voltage (CV), relaxation (RX), or initial relaxation mode
    INPAR::ELCH::CCCVHalfCyclePhase phase_cccv_;

    //! end time of current relaxation phase
    double relaxendtime_;

    //! duration of relaxation phase
    const double relaxtime_;
  };  // class CCCVHalfCycleCondition

  /*========================================================================*/
  /*========================================================================*/

  /*!
   * \brief Control routine for constant-current, constant-voltage condition
   *
   * Holds two half cycles (charge and discharge). Controls, which half cycle is activated. Serves
   * as interface between one half cycle and time integration.
   */
  class CCCVCondition
  {
   public:
    //! constructor
    CCCVCondition(const DRT::Condition& cccvcyclingcondition,
        const std::vector<DRT::Condition*>& cccvhalfcycleconditions, bool adaptivetimestepping,
        int num_dofs);

    //! true, when all half cylces are completed
    bool NotFinished() const { return nhalfcycles_ >= ihalfcycle_; };

    //! phase of active half cycle
    INPAR::ELCH::CCCVHalfCyclePhase GetCCCVHalfCyclePhase() const;

    //! ID of current half cycle
    int GetHalfCycleConditionID() const;

    double GetInitialRelaxTime() const { return initrelaxtime_; };

    //! number of current half cycle
    int GetNumCurrentHalfCycle() const { return ihalfcycle_; };

    //! Step when phase was changed last time
    int GetStepLastPhaseChange() const { return steplastphasechange_; };

    //! get end time of current relaxation phase
    double GetRelaxEndTime() const;

    //! does this phase have adaptive time stepping?
    bool IsAdaptiveTimeSteppingPhase() const;

    //! is this phase finished?
    bool IsEndOfHalfCyclePhase(double cellvoltage, double cellcrate, double time) const;

    //! is cut off c rate exceeded?
    bool ExceedCellCRate(double expected_cellcrate) const;

    //! does cell voltage exceed bounds of current half cycle
    bool ExceedCellVoltage(double expected_cellvoltage) const;

    //! is this condition in initial relaxation?
    bool IsInitialRelaxation(const double time, const double dt) const
    {
      return time <= initrelaxtime_ - dt;
    }

    //! was phase changed since last adaption of time step?
    bool IsPhaseChanged() const { return phasechanged_; };

    //! are we in initial relaxation phase?
    bool IsPhaseInitialRelaxation() const { return phaseinitialrelaxation_; };

    //! true if difference between @p step and the step when the phase was last changed is equal
    //! or more than  num_add_adapt_timesteps_
    bool ExceedMaxStepsFromLastPhaseChange(int step);

    //! return minimum amount of time steps the initial relaxation time of cccv condition has to be
    //! discretized with
    int MinTimeStepsDuringInitRelax() const { return min_time_steps_during_init_relax_; }

    //! go to next phase. If half cycle is finished switch to other half cycle and start with
    //! constant current
    void NextPhase(int step, double time, bool print);

    //! number of dofs of this cccv condition
    int NumDofs() const { return num_dofs_; }

    //! read restart
    void ReadRestart(IO::DiscretizationReader& reader);

    //! reset phasechanged_
    void ResetPhaseChangeObserver();

    //! initialize first half cycle
    void SetFirstCCCVHalfCycle(int step);

    //! start counting steps from now on (phase was changed)
    void SetPhaseChangeObserver(int step);

   private:
    //! adaptive time stepping for initial relaxation?
    const int adaptivetimesteppingonoff_;

    //! how to begin with CCCV condition (charge/discharge)
    const bool beginwithcharge_;

    //! flag indicating whether cell is currently being charged or discharged
    bool charging_;

    //! half cycle of charge
    Teuchos::RCP<SCATRA::CCCVHalfCycleCondition> halfcycle_charge_;

    //! half cycle of discharge
    Teuchos::RCP<SCATRA::CCCVHalfCycleCondition> halfcycle_discharge_;

    //! number of current charge or discharge half-cycle
    int ihalfcycle_;

    //! initial relaxation time of cccv condition
    const double initrelaxtime_;

    //! minimum amount of time steps the initial relaxation time of cccv condition has to be
    //! discretized with
    const int min_time_steps_during_init_relax_;

    //! total number of charge and discharge half-cycles
    const int nhalfcycles_;

    //! number of steps after phase was changed until reset time step to original value
    const int num_add_adapt_timesteps_;

    //! number of dofs of this cccv condition
    const int num_dofs_;

    //! for time step adaptivity: was phase changed since last time step adaptivity?
    bool phasechanged_;

    //! indicating, if currently in initial relaxation
    bool phaseinitialrelaxation_;

    //! step when phase lastly changed
    int steplastphasechange_;

  };  // class CCCVCondition

}  // namespace SCATRA
BACI_NAMESPACE_CLOSE

#endif
