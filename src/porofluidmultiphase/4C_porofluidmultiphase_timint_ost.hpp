// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_POROFLUIDMULTIPHASE_TIMINT_OST_HPP
#define FOUR_C_POROFLUIDMULTIPHASE_TIMINT_OST_HPP

#include "4C_config.hpp"

#include "4C_linalg_utils_densematrix_communication.hpp"
#include "4C_porofluidmultiphase_timint_implicit.hpp"

FOUR_C_NAMESPACE_OPEN

namespace POROFLUIDMULTIPHASE
{
  class TimIntOneStepTheta : public TimIntImpl
  {
   public:
    /// Standard Constructor
    TimIntOneStepTheta(std::shared_ptr<Core::FE::Discretization> dis,  //!< discretization
        const int linsolvernumber,                                     //!< number of linear solver
        const Teuchos::ParameterList& probparams,  //!< parameter list of global problem
        const Teuchos::ParameterList& poroparams,  //!< parameter list of poro problem
        std::shared_ptr<Core::IO::DiscretizationWriter> output  //!< output writer
    );


    /// Print information about current time step to screen (reimplementation for OST)
    void print_time_step_info() override;

    /// compute values at intermediate time steps (required for generalized-alpha)
    void compute_intermediate_values() override { return; };

    ///  compute scalar time derivative
    void compute_time_derivative() override;

    /// update the solution after convergence of the nonlinear iteration.
    /// current solution becomes old solution of next timestep.
    void update() override;

    /// read restart data
    void read_restart(const int step) override;

   protected:
    /// don't want = operator and cctor
    TimIntOneStepTheta operator=(const TimIntOneStepTheta& old);

    /// copy constructor
    TimIntOneStepTheta(const TimIntOneStepTheta& old);

    /// set time parameter for element evaluation (called before every time step)
    void set_element_time_step_parameter() const override;

    //! set time for evaluation of Neumann boundary conditions
    void set_time_for_neumann_evaluation(Teuchos::ParameterList& params) override;

    //! initialization procedure prior to evaluation of first time step
    void calc_initial_time_derivative() override;

    /// set part of residual vector belonging to previous time step
    void set_old_part_of_righthandside() override;

    /// do explicit predictor step (-> better starting value for nonlinear solver)
    virtual void explicit_predictor();

    /// add actual Neumann loads with time factor
    void add_neumann_to_residual() override;

    /// add parameters specific for time-integration scheme
    void add_time_integration_specific_vectors() override;

    /// write additional data required for restart
    void output_restart() override;

    /// return the right time-scaling-factor for the true residual
    double residual_scaling() const override { return 1.0 / (dt_ * theta_); }

    /// time factor for one-step-theta/BDF2 time integration
    double theta_;

  };  // class TimIntOneStepTheta

}  // namespace POROFLUIDMULTIPHASE



FOUR_C_NAMESPACE_CLOSE

#endif
