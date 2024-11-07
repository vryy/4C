// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

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
    ScaTraTimIntElchOST(std::shared_ptr<Core::FE::Discretization> dis,
        std::shared_ptr<Core::LinAlg::Solver> solver,
        std::shared_ptr<Teuchos::ParameterList> params,
        std::shared_ptr<Teuchos::ParameterList> sctratimintparams,
        std::shared_ptr<Teuchos::ParameterList> extraparams,
        std::shared_ptr<Core::IO::DiscretizationWriter> output);

    void init() override;

    void setup() override;

    void update() override;

    void read_restart(
        const int step, std::shared_ptr<Core::IO::InputControl> input = nullptr) override;

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
    ScaTraTimIntElchSCLOST(std::shared_ptr<Core::FE::Discretization> dis,
        std::shared_ptr<Core::LinAlg::Solver> solver,
        std::shared_ptr<Teuchos::ParameterList> params,
        std::shared_ptr<Teuchos::ParameterList> sctratimintparams,
        std::shared_ptr<Teuchos::ParameterList> extraparams,
        std::shared_ptr<Core::IO::DiscretizationWriter> output);

    void add_time_integration_specific_vectors(bool forcedincrementalsolver = false) override;

    void init() override;

    void post_calc_initial_potential_field() override;

    void pre_calc_initial_potential_field() override;

    void read_restart(
        const int step, std::shared_ptr<Core::IO::InputControl> input = nullptr) override;

    void setup() override;

    void update() override;

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
    ScaTraTimIntElchBDF2(std::shared_ptr<Core::FE::Discretization> dis,
        std::shared_ptr<Core::LinAlg::Solver> solver,
        std::shared_ptr<Teuchos::ParameterList> params,
        std::shared_ptr<Teuchos::ParameterList> sctratimintparams,
        std::shared_ptr<Teuchos::ParameterList> extraparams,
        std::shared_ptr<Core::IO::DiscretizationWriter> output);

    void init() override;

    void setup() override;

    void update() override;

    void read_restart(
        const int step, std::shared_ptr<Core::IO::InputControl> input = nullptr) override;

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
    ScaTraTimIntElchGenAlpha(std::shared_ptr<Core::FE::Discretization> dis,
        std::shared_ptr<Core::LinAlg::Solver> solver,
        std::shared_ptr<Teuchos::ParameterList> params,
        std::shared_ptr<Teuchos::ParameterList> sctratimintparams,
        std::shared_ptr<Teuchos::ParameterList> extraparams,
        std::shared_ptr<Core::IO::DiscretizationWriter> output);

    void init() override;

    void setup() override;

    void update() override;

    void read_restart(
        const int step, std::shared_ptr<Core::IO::InputControl> input = nullptr) override;

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
    ScaTraTimIntElchStationary(std::shared_ptr<Core::FE::Discretization> dis,
        std::shared_ptr<Core::LinAlg::Solver> solver,
        std::shared_ptr<Teuchos::ParameterList> params,
        std::shared_ptr<Teuchos::ParameterList> sctratimintparams,
        std::shared_ptr<Teuchos::ParameterList> extraparams,
        std::shared_ptr<Core::IO::DiscretizationWriter> output);

    void init() override;

    void setup() override;

    void update() override;

    void read_restart(
        const int step, std::shared_ptr<Core::IO::InputControl> input = nullptr) override;

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
