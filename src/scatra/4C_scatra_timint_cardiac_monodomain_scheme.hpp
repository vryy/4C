// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SCATRA_TIMINT_CARDIAC_MONODOMAIN_SCHEME_HPP
#define FOUR_C_SCATRA_TIMINT_CARDIAC_MONODOMAIN_SCHEME_HPP

#include "4C_config.hpp"

#include "4C_scatra_timint_bdf2.hpp"
#include "4C_scatra_timint_cardiac_monodomain.hpp"
#include "4C_scatra_timint_genalpha.hpp"
#include "4C_scatra_timint_ost.hpp"

FOUR_C_NAMESPACE_OPEN


namespace ScaTra
{
  class TimIntCardiacMonodomainOST : public virtual TimIntCardiacMonodomain,
                                     public virtual TimIntOneStepTheta
  {
   public:
    //! Standard Constructor
    TimIntCardiacMonodomainOST(std::shared_ptr<Core::FE::Discretization> dis,
        std::shared_ptr<Core::LinAlg::Solver> solver,
        std::shared_ptr<Teuchos::ParameterList> params,
        std::shared_ptr<Teuchos::ParameterList> sctratimintparams,
        std::shared_ptr<Teuchos::ParameterList> extraparams,
        std::shared_ptr<Core::IO::DiscretizationWriter> output);


    //! setup time integration scheme
    void setup() override;

    //! Update the solution after convergence of the nonlinear iteration.
    //! Current solution becomes old solution of next timestep.
    void update() override;

    //! read restart data
    void read_restart(
        const int step, std::shared_ptr<Core::IO::InputControl> input = nullptr) override;

   protected:
    void write_restart() const override;

    /// add parameters specific for time-integration scheme
    void add_time_integration_specific_vectors(bool forcedincrementalsolver = false) override;

    //! do not calculate initial scalar time derivatives for ep
    void calc_initial_time_derivative() override { return; };

  };  // class TimIntCardiacMonodomainOST


  class TimIntCardiacMonodomainBDF2 : public virtual TimIntCardiacMonodomain,
                                      public virtual TimIntBDF2
  {
   public:
    //! Standard Constructor
    TimIntCardiacMonodomainBDF2(std::shared_ptr<Core::FE::Discretization> dis,
        std::shared_ptr<Core::LinAlg::Solver> solver,
        std::shared_ptr<Teuchos::ParameterList> params,
        std::shared_ptr<Teuchos::ParameterList> sctratimintparams,
        std::shared_ptr<Teuchos::ParameterList> extraparams,
        std::shared_ptr<Core::IO::DiscretizationWriter> output);


    //! setup time integration scheme
    void setup() override;

    //! Update the solution after convergence of the nonlinear iteration.
    //! Current solution becomes old solution of next timestep.
    void update() override;

    //! read restart data
    void read_restart(
        const int step, std::shared_ptr<Core::IO::InputControl> input = nullptr) override;

   protected:
    void write_restart() const override;

    //! do not calculate initial scalar time derivatives for ep
    void calc_initial_time_derivative() override { return; };

  };  // class TimIntCardiacMonodomainBDF2


  class TimIntCardiacMonodomainGenAlpha : public virtual TimIntCardiacMonodomain,
                                          public virtual TimIntGenAlpha
  {
   public:
    //! Standard Constructor
    TimIntCardiacMonodomainGenAlpha(std::shared_ptr<Core::FE::Discretization> dis,
        std::shared_ptr<Core::LinAlg::Solver> solver,
        std::shared_ptr<Teuchos::ParameterList> params,
        std::shared_ptr<Teuchos::ParameterList> sctratimintparams,
        std::shared_ptr<Teuchos::ParameterList> extraparams,
        std::shared_ptr<Core::IO::DiscretizationWriter> output);


    //! Setup time integration scheme
    void setup() override;

    //! Update the solution after convergence of the nonlinear iteration.
    //! Current solution becomes old solution of next timestep.
    void update() override;

    //! read restart data
    void read_restart(
        const int step, std::shared_ptr<Core::IO::InputControl> input = nullptr) override;

    /// add parameters specific for time-integration scheme
    void add_time_integration_specific_vectors(bool forcedincrementalsolver = false) override;

   protected:
    void write_restart() const override;

    //! do not calculate initial scalar time derivatives for ep
    void calc_initial_time_derivative() override { return; };

  };  // class TimIntCardiacMonodomainGenAlpha
}  // namespace ScaTra

FOUR_C_NAMESPACE_CLOSE

#endif
