// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SCATRA_TIMINT_HDG_HPP
#define FOUR_C_SCATRA_TIMINT_HDG_HPP

#include "4C_config.hpp"

#include "4C_fem_discretization_hdg.hpp"
#include "4C_scatra_timint_genalpha.hpp"

FOUR_C_NAMESPACE_OPEN

namespace ScaTra
{
  /*!
  \brief time integration for HDG scatra
  */
  class TimIntHDG : public TimIntGenAlpha
  {
   public:
    //! standard constructor
    TimIntHDG(const std::shared_ptr<Core::FE::Discretization>& actdis,
        const std::shared_ptr<Core::LinAlg::Solver>& solver,
        const std::shared_ptr<Teuchos::ParameterList>& params,
        const std::shared_ptr<Teuchos::ParameterList>& extraparams,
        std::shared_ptr<Core::IO::DiscretizationWriter> output);

    //! setup
    void setup() override;

    void setup_context_vector() override {};

    //! set theta_ to its value, dependent on integration method for GenAlpha and BDF2
    virtual void set_theta();

    //! set states in the time integration schemes: additional vectors for HDG
    void add_time_integration_specific_vectors(bool forcedincrementalsolver = false) override;

    //! set the part of the right hand side belonging to the last time step
    void set_old_part_of_righthandside() override;

    //! update the solution after convergence of the nonlinear iteration,
    //! current solution becomes old solution of next time step
    void update() override;

    void write_restart() const override;

    void collect_runtime_output_data() override;

    //! read restart
    void read_restart(
        const int step, std::shared_ptr<Core::IO::InputControl> input = nullptr) override;

    //! set the initial scalar field phi
    void set_initial_field(const Inpar::ScaTra::InitialField init,  //!< type of initial field
        const int startfuncno                                       //!< number of spatial function
        ) override;

    //! accessor to interior concentrations
    virtual std::shared_ptr<Core::LinAlg::Vector<double>> return_int_phinp() { return intphinp_; }
    virtual std::shared_ptr<Core::LinAlg::Vector<double>> return_int_phin() { return intphin_; }

    virtual std::shared_ptr<Core::LinAlg::Vector<double>> interpolated_phinp() const
    {
      return interpolatedPhinp_;
    }

    //! prepare time loop
    void prepare_time_loop() override;

    /*!
    \brief Compare the numerical solution to the analytical one.
    */
    virtual std::shared_ptr<Core::LinAlg::SerialDenseVector> compute_error() const;

    std::shared_ptr<Core::Utils::ResultTest> create_scatra_field_test() override;

   protected:
    //! copy constructor
    TimIntHDG(const TimIntHDG& old);

    //! update time derivative for generalized-alpha time integration
    virtual void gen_alpha_compute_time_derivative();

    //! compute values at intermediate time steps for gen.-alpha
    virtual void gen_alpha_intermediate_values();

    //! number of dofset for interior variables
    int nds_intvar_;

    //! @name concentration and concentration gradient at different times for element interior for
    //! HDG
    //@{
    std::shared_ptr<Core::LinAlg::Vector<double>>
        intphinp_;  //!< concentration at time \f$t^{n+1}\f$
    std::shared_ptr<Core::LinAlg::Vector<double>> intphin_;  //!< concentration at time \f$t^{n}\f$
    //@}

    //! @name other HDG-specific auxiliary vectors
    //@{
    std::shared_ptr<Core::LinAlg::Vector<double>>
        interpolatedPhinp_;  //!< concentrations for output at time \f$t^{n+1}\f$
    //@}

    //! calculate intermediate solution
    void compute_intermediate_values() override;

    //! compute values at the interior of the elements
    void compute_interior_values() override;

    //! update interior variables
    virtual void update_interior_variables(
        std::shared_ptr<Core::LinAlg::Vector<double>> updatevector);

    //! write problem specific output
    virtual void write_problem_specific_output(
        std::shared_ptr<Core::LinAlg::Vector<double>> interpolatedPhi)
    {
    }

    virtual void collect_problem_specific_runtime_output_data(
        std::shared_ptr<Core::LinAlg::Vector<double>> interpolatedPhi)
    {
    }

    //! calculate consistent initial scalar time derivatives in compliance with initial scalar field
    void calc_initial_time_derivative() override {}

    //! fd check
    void fd_check() override;

    //! calculation of error with reference to analytical solution during the simulation
    void evaluate_error_compared_to_analytical_sol() override;

    //! adapt degree of test functions and change dofsets accordingly
    virtual void adapt_degree();

    //! adapt variable vectors required due to the change of the degrees of the test functions
    virtual void adapt_variable_vector(std::shared_ptr<Core::LinAlg::Vector<double>> phi_new,
        std::shared_ptr<Core::LinAlg::Vector<double>> phi_old,
        std::shared_ptr<Core::LinAlg::Vector<double>> intphi_new,
        std::shared_ptr<Core::LinAlg::Vector<double>> intphi_old, int nds_var_old,
        int nds_intvar_old, std::vector<Core::Elements::LocationArray> la_old);

    //! calculate matrices on element
    virtual void calc_mat_initial();

    //! chooses the assemble process (assemble matrix and rhs or only rhs)
    void assemble_mat_and_rhs() override;

    //! contains the assembly process only for rhs
    void assemble_rhs();

    //! pack material
    virtual void pack_material() {}

    //! adapt material
    virtual void unpack_material() {}

    //! project material field
    virtual void project_material() {}

   private:
    //! time algorithm flag actually set (we internally reset it)
    Inpar::ScaTra::TimeIntegrationScheme timealgoset_;

    //! @name time stepping variable
    bool startalgo_;  //!< flag for starting algorithm

    //! time-integration-scheme factor theta
    double theta_;

    //! activation_time at times n+1
    std::shared_ptr<Core::LinAlg::Vector<double>> activation_time_interpol_np_;

    //! HDG discretization
    Core::FE::DiscretizationHDG* hdgdis_;

    //! p-adativitity
    bool padaptivity_;

    //! error tolerance for p-adativitity
    double padapterrortol_;

    //! error tolerance base for p-adativitity (delta p=log_padapterrorbase_(E/padapterrortol_))
    double padapterrorbase_;


    //! max. degree of shape functions for p-adativitity
    double padaptdegreemax_;

    //! element degree
    std::shared_ptr<Core::LinAlg::Vector<double>> elementdegree_;

  };  // class TimIntHDG
}  // namespace ScaTra

FOUR_C_NAMESPACE_CLOSE

#endif
