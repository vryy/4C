/*----------------------------------------------------------------------*/
/*! \file

\brief service routines for scatra time integration for elch

\level 2


 *------------------------------------------------------------------------------------------------*/

#ifndef FOUR_C_SCATRA_TIMINT_ELCH_SERVICE_HPP
#define FOUR_C_SCATRA_TIMINT_ELCH_SERVICE_HPP

#include "4C_config.hpp"

#include "4C_discretization_condition.hpp"
#include "4C_inpar_elch.hpp"
#include "4C_scatra_timint_implicit.hpp"

FOUR_C_NAMESPACE_OPEN

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
    CCCVHalfCycleCondition(
        const CORE::Conditions::Condition& cccvhalfcyclecondition, bool adaptivetimestepping);

    //! Get phase of half cycle
    INPAR::ELCH::CCCVHalfCyclePhase get_cccv_half_cycle_phase() const { return phase_cccv_; };

    //! Get ID of this half cycle condition
    int GetConditionID() const { return halfcyclecondition_id_; };

    //! get cut off c-rate during constant voltage
    double GetCutOffCRate() const { return cutoffcrate_; };

    //! get cut off voltage during constant current
    double GetCutOffVoltage() const { return cutoffvoltage_; };

    //! get end time of current relaxation phase
    double GetRelaxEndTime() const { return relaxendtime_; };

    //! does this phase have adaptive time stepping?
    bool is_adaptive_time_stepping_phase() const;

    //! is this half cycle completed after updating the phase?
    bool is_end_of_half_cycle_next_phase(double time, bool print);

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
    const int halfcyclecondition_id_;

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
    CCCVCondition(const CORE::Conditions::Condition& cccvcyclingcondition,
        const std::vector<CORE::Conditions::Condition*>& cccvhalfcycleconditions,
        bool adaptivetimestepping, int num_dofs);

    //! true, when all half cylces are completed
    bool NotFinished() const { return nhalfcycles_ >= ihalfcycle_; };

    //! phase of active half cycle
    INPAR::ELCH::CCCVHalfCyclePhase get_cccv_half_cycle_phase() const;

    //! ID of current half cycle
    int get_half_cycle_condition_id() const;

    double GetInitialRelaxTime() const { return initrelaxtime_; };

    //! number of current half cycle
    int get_num_current_half_cycle() const { return ihalfcycle_; };

    //! Step when phase was changed last time
    int get_step_last_phase_change() const { return steplastphasechange_; };

    //! get end time of current relaxation phase
    double GetRelaxEndTime() const;

    //! does this phase have adaptive time stepping?
    bool is_adaptive_time_stepping_phase() const;

    //! is this phase finished?
    bool is_end_of_half_cycle_phase(double cellvoltage, double cellcrate, double time) const;

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
    bool is_phase_initial_relaxation() const { return phaseinitialrelaxation_; };

    //! true if difference between @p step and the step when the phase was last changed is equal
    //! or more than  num_add_adapt_timesteps_
    bool exceed_max_steps_from_last_phase_change(int step);

    //! return minimum amount of time steps the initial relaxation time of cccv condition has to be
    //! discretized with
    int min_time_steps_during_init_relax() const { return min_time_steps_during_init_relax_; }

    //! go to next phase. If half cycle is finished switch to other half cycle and start with
    //! constant current
    void NextPhase(int step, double time, bool print);

    //! number of dofs of this cccv condition
    int NumDofs() const { return num_dofs_; }

    //! read restart
    void ReadRestart(IO::DiscretizationReader& reader);

    //! reset phasechanged_
    void reset_phase_change_observer();

    //! initialize first half cycle
    void set_first_cccv_half_cycle(int step);

    //! start counting steps from now on (phase was changed)
    void set_phase_change_observer(int step);

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
FOUR_C_NAMESPACE_CLOSE

#endif
