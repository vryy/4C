/*----------------------------------------------------------------------*/
/*! \file

\brief  connecting time-integration schemes (OST, BDF2, GenAlpha, Stationary) with
        elch-specific implementation (class ScaTraTimIntElch)
\level 2



*----------------------------------------------------------------------*/

#ifndef FOUR_C_SCATRA_TIMINT_ELCH_SCHEME_HPP
#define FOUR_C_SCATRA_TIMINT_ELCH_SCHEME_HPP

#include "4C_config.hpp"

#include "4C_scatra_timint_bdf2.hpp"
#include "4C_scatra_timint_elch.hpp"
#include "4C_scatra_timint_elch_scl.hpp"
#include "4C_scatra_timint_genalpha.hpp"
#include "4C_scatra_timint_ost.hpp"
#include "4C_scatra_timint_stat.hpp"

FOUR_C_NAMESPACE_OPEN

namespace ScaTra
{
  class ScaTraTimIntElchOST : public ScaTraTimIntElch, public TimIntOneStepTheta
  {
   public:
    //! Standard Constructor
    ScaTraTimIntElchOST(Teuchos::RCP<Core::FE::Discretization> dis,
        Teuchos::RCP<Core::LinAlg::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<Core::IO::DiscretizationWriter> output);

    void init() override;

    void setup() override;

    void Update() override;

    void read_restart(
        const int step, Teuchos::RCP<Core::IO::InputControl> input = Teuchos::null) override;

    void pre_calc_initial_potential_field() override;

    void post_calc_initial_potential_field() override;

   protected:
    void write_restart() const override;

    void electrode_kinetics_time_update() override;

    void explicit_predictor() const override;

    void compute_time_deriv_pot0(const bool init) override;

    void set_old_part_of_righthandside() override;
  };

  class ScaTraTimIntElchSCLOST : public ScaTraTimIntElchSCL, public TimIntOneStepTheta
  {
   public:
    //! Standard Constructor
    ScaTraTimIntElchSCLOST(Teuchos::RCP<Core::FE::Discretization> dis,
        Teuchos::RCP<Core::LinAlg::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<Core::IO::DiscretizationWriter> output);

    void add_time_integration_specific_vectors(bool forcedincrementalsolver = false) override;

    void init() override;

    void post_calc_initial_potential_field() override;

    void pre_calc_initial_potential_field() override;

    void read_restart(
        const int step, Teuchos::RCP<Core::IO::InputControl> input = Teuchos::null) override;

    void setup() override;

    void Update() override;

   protected:
    void compute_time_deriv_pot0(const bool init) override{};

    void electrode_kinetics_time_update() override{};

    void explicit_predictor() const override;

    void set_old_part_of_righthandside() override;

    void write_restart() const override;
  };


  class ScaTraTimIntElchBDF2 : public ScaTraTimIntElch, public TimIntBDF2
  {
   public:
    //! Standard Constructor
    ScaTraTimIntElchBDF2(Teuchos::RCP<Core::FE::Discretization> dis,
        Teuchos::RCP<Core::LinAlg::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<Core::IO::DiscretizationWriter> output);

    void init() override;

    void setup() override;

    void Update() override;

    void read_restart(
        const int step, Teuchos::RCP<Core::IO::InputControl> input = Teuchos::null) override;

    void pre_calc_initial_potential_field() override;

    void post_calc_initial_potential_field() override{};

   protected:
    void write_restart() const override;

    void electrode_kinetics_time_update() override;

    void compute_time_deriv_pot0(const bool init) override;

    void set_old_part_of_righthandside() override;
  };


  class ScaTraTimIntElchGenAlpha : public ScaTraTimIntElch, public TimIntGenAlpha
  {
   public:
    //! Standard Constructor
    ScaTraTimIntElchGenAlpha(Teuchos::RCP<Core::FE::Discretization> dis,
        Teuchos::RCP<Core::LinAlg::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<Core::IO::DiscretizationWriter> output);

    void init() override;

    void setup() override;

    void Update() override;

    void read_restart(
        const int step, Teuchos::RCP<Core::IO::InputControl> input = Teuchos::null) override;

    void pre_calc_initial_potential_field() override;

    void post_calc_initial_potential_field() override;

   protected:
    void write_restart() const override;

    void electrode_kinetics_time_update() override;

    void compute_time_deriv_pot0(const bool init) override;
  };


  class ScaTraTimIntElchStationary : public ScaTraTimIntElch, public TimIntStationary
  {
   public:
    //! Standard Constructor
    ScaTraTimIntElchStationary(Teuchos::RCP<Core::FE::Discretization> dis,
        Teuchos::RCP<Core::LinAlg::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<Core::IO::DiscretizationWriter> output);

    void init() override;

    void setup() override;

    void Update() override;

    void read_restart(
        const int step, Teuchos::RCP<Core::IO::InputControl> input = Teuchos::null) override;

    void pre_calc_initial_potential_field() override;

    void post_calc_initial_potential_field() override{};

   protected:
    void write_restart() const override;

    void electrode_kinetics_time_update() override
    {
      FOUR_C_THROW(
          "Galvanostatic-BC is not implemented for the stationary time-integration scheme");
    };

    void compute_time_deriv_pot0(const bool init) override;
  };
}  // namespace ScaTra

FOUR_C_NAMESPACE_CLOSE

#endif
