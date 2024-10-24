// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_LUBRICATION_TIMINT_STAT_HPP
#define FOUR_C_LUBRICATION_TIMINT_STAT_HPP


#include "4C_config.hpp"

#include "4C_lubrication_timint_implicit.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Lubrication
{
  class TimIntStationary : public TimIntImpl
  {
   public:
    /// Standard Constructor
    TimIntStationary(Teuchos::RCP<Core::FE::Discretization> dis,
        Teuchos::RCP<Core::LinAlg::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<Core::IO::DiscretizationWriter> output);


    /// initialize time integration scheme
    void init() override;

    /// Update the solution after convergence of the nonlinear iteration.
    /// Current solution becomes old solution of next timestep.
    void update(const int num = 0) override;

    /// read restart data
    void read_restart(int step) override;

    //! Update iteration incrementally
    //!
    //! This update is carried out by computing the new #raten_
    //! from scratch by using the newly updated #prenp_. The method
    //! respects the Dirichlet DOFs which are not touched.
    //! This method is necessary for certain predictors
    //! (like #predict_const_temp_consist_rate)
    void update_iter_incrementally() override;

   protected:
    /// don't want = operator and cctor
    TimIntStationary operator=(const TimIntStationary& old);

    /// copy constructor
    TimIntStationary(const TimIntStationary& old);

    /// set time parameter for element evaluation
    void set_element_time_parameter() const override;

    //! set time for evaluation of Neumann boundary conditions
    void set_time_for_neumann_evaluation(Teuchos::ParameterList& params) override;

    /// add actual Neumann loads with time factor
    void add_neumann_to_residual() override;

    /// add parameters specific for time-integration scheme
    void add_time_integration_specific_vectors(bool forcedincrementalsolver = false) override;

    /// return the right time-scaling-factor for the true residual
    double residual_scaling() const override { return 1.0; }

  };  // class TimIntStationary

}  // namespace Lubrication

FOUR_C_NAMESPACE_CLOSE

#endif
