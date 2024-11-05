// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SCATRA_TIMINT_STAT_HDG_HPP
#define FOUR_C_SCATRA_TIMINT_STAT_HDG_HPP

#include "4C_config.hpp"

#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_scatra_timint_hdg.hpp"

FOUR_C_NAMESPACE_OPEN

namespace ScaTra
{
  class TimIntStationaryHDG : public TimIntHDG
  {
   public:
    /// Standard Constructor
    TimIntStationaryHDG(std::shared_ptr<Core::FE::Discretization> dis,
        std::shared_ptr<Core::LinAlg::Solver> solver,
        std::shared_ptr<Teuchos::ParameterList> params,
        std::shared_ptr<Teuchos::ParameterList> extraparams,
        std::shared_ptr<Core::IO::DiscretizationWriter> output);

    /// initialize time integration scheme
    void init() override;

    /// compute values at intermediate time steps (required for generalized-alpha)
    void compute_intermediate_values() override { return; };

    /// routine to return time integration specific parameters
    std::shared_ptr<Teuchos::ParameterList> scatra_time_parameter_list() override
    {
      FOUR_C_THROW("Not yet implemented!");
      return nullptr;
    }

   protected:
    /// don't want = operator and cctor
    TimIntStationaryHDG operator=(const TimIntStationaryHDG& old);

    /// copy constructor
    TimIntStationaryHDG(const TimIntStationaryHDG& old);

    /// set time parameter for element evaluation
    void set_element_time_parameter(bool forcedincrementalsolver = false) const override;

    /// set time for evaluation of Neumann boundary conditions
    void set_time_for_neumann_evaluation(Teuchos::ParameterList& params) override;

    /// calculate consistent initial conditions in compliance with initial scalar field
    /// this is not necessary for stationary calculations
    void calc_initial_time_derivative() override { return; };

    /// do explicit predictor step (nothing to predict for stationary problems!)
    void explicit_predictor() const override { return; };

    /// dynamic Smagorinsky model
    void dynamic_computation_of_cs() override
    {
      FOUR_C_THROW("no turbulence in stationary flows!");
      return;
    };

    /// dynamic Vreman model
    void dynamic_computation_of_cv() override
    {
      FOUR_C_THROW("no turbulence in stationary flows!");
      return;
    };

  };  // class TimIntStationaryHDG

}  // namespace ScaTra

FOUR_C_NAMESPACE_CLOSE

#endif