/*----------------------------------------------------------------------*/
/*! \file
\brief Generalized-alpha time-integration scheme

\level 1


*----------------------------------------------------------------------*/

#ifndef FOUR_C_SCATRA_TIMINT_GENALPHA_HPP
#define FOUR_C_SCATRA_TIMINT_GENALPHA_HPP

#include "4C_config.hpp"

#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_scatra_timint_implicit.hpp"

FOUR_C_NAMESPACE_OPEN

namespace ScaTra
{
  class TimIntGenAlpha : public virtual ScaTraTimIntImpl
  {
   public:
    /// Standard Constructor
    TimIntGenAlpha(Teuchos::RCP<Core::FE::Discretization> dis,
        Teuchos::RCP<Core::LinAlg::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<Core::IO::DiscretizationWriter> output);

    void setup() override;

    void print_time_step_info() override
    {
      if (myrank_ == 0)
      {
        printf(
            "\nTIME: %11.4E/%11.4E  DT = %11.4E  %s(a_F=%3.2f | a_M=%3.2f | gamma=%3.2f) STEP = "
            "%4d/%4d\n",
            time_, maxtime_, dta_, method_title().c_str(), alphaF_, alphaM_, gamma_, step_,
            stepmax_);
      }
    }

    void compute_intermediate_values() override;

    void compute_interior_values() override{};

    void compute_time_derivative() override;

    void compute_time_deriv_pot0(const bool init) override{};

    void update() override;

    void read_restart(
        const int step, Teuchos::RCP<Core::IO::InputControl> input = Teuchos::null) override;

    Teuchos::RCP<Epetra_Vector> phiaf() override { return phiaf_; }

    Teuchos::RCP<Epetra_Vector> phiafnp() override { return phiaf_; }

    Teuchos::RCP<Epetra_Vector> phiam() override { return phiam_; }

    Teuchos::RCP<Epetra_Vector> phidtam() override { return phidtam_; }

    Teuchos::RCP<Epetra_Vector> fs_phi() override
    {
      if (Sep_ != Teuchos::null) Sep_->multiply(false, *phiaf_, *fsphiaf_);
      return fsphiaf_;
    };

    Teuchos::RCP<Teuchos::ParameterList> scatra_time_parameter_list() override
    {
      Teuchos::RCP<Teuchos::ParameterList> timeparams;
      timeparams = Teuchos::rcp(new Teuchos::ParameterList());
      timeparams->set("using stationary formulation", false);
      timeparams->set("using generalized-alpha time integration", true);
      timeparams->set("total time", time_ - (1 - alphaF_) * dta_);
      timeparams->set("time factor", genalphafac_ * dta_);
      timeparams->set("alpha_F", alphaF_);
      return timeparams;
    }

    void pre_calc_initial_time_derivative() override;

    void post_calc_initial_time_derivative() override;


   protected:
    /// don't want = operator and cctor
    TimIntGenAlpha operator=(const TimIntGenAlpha& old);

    /// copy constructor
    TimIntGenAlpha(const TimIntGenAlpha& old);

    void set_element_time_parameter(bool forcedincrementalsolver = false) const override;

    void set_time_for_neumann_evaluation(Teuchos::ParameterList& params) override;

    void set_element_time_parameter_backward_euler() const override;

    void calc_initial_time_derivative() override;

    void set_old_part_of_righthandside() override;

    void explicit_predictor() const override;

    void add_neumann_to_residual() override;

    void av_m3_separation() override;

    void dynamic_computation_of_cs() override;

    void dynamic_computation_of_cv() override;

    void add_time_integration_specific_vectors(bool forcedincrementalsolver = false) override;

    void write_restart() const override;

    double residual_scaling() const override { return 1.0 / (dta_ * genalphafac_); }

    /// scalar at time n+alpha_F and n+alpha_M
    Teuchos::RCP<Epetra_Vector> phiaf_;
    Teuchos::RCP<Epetra_Vector> phiam_;

    /// scalar time derivative at time n+alpha_M
    Teuchos::RCP<Epetra_Vector> phidtam_;

    /// fine-scale part at time n+alpha_F
    Teuchos::RCP<Epetra_Vector> fsphiaf_;

    /// time factors for generalized-alpha time integration
    double alphaM_;
    double alphaF_;
    double gamma_;
    double genalphafac_;
  };  // class TimIntGenAlpha

}  // namespace ScaTra

FOUR_C_NAMESPACE_CLOSE

#endif
