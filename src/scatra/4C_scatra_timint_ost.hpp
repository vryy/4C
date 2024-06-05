/*----------------------------------------------------------------------*/
/*! \file
\brief One-Step-Theta time-integration scheme

\level 1


*----------------------------------------------------------------------*/

#ifndef FOUR_C_SCATRA_TIMINT_OST_HPP
#define FOUR_C_SCATRA_TIMINT_OST_HPP

#include "4C_config.hpp"

#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_scatra_timint_implicit.hpp"

FOUR_C_NAMESPACE_OPEN

namespace SCATRA
{
  class TimIntOneStepTheta : public virtual ScaTraTimIntImpl
  {
   public:
    /// Standard Constructor
    TimIntOneStepTheta(Teuchos::RCP<DRT::Discretization> actdis,  //!< discretization
        Teuchos::RCP<CORE::LINALG::Solver> solver,                //!< linear solver
        Teuchos::RCP<Teuchos::ParameterList> params,              //!< parameter list
        Teuchos::RCP<Teuchos::ParameterList> extraparams,         //!< supplementary parameter list
        Teuchos::RCP<CORE::IO::DiscretizationWriter> output,      //!< output writer
        const int probnum = 0                                     //!< global problem number
    );

    /// don't want = operator and cctor
    TimIntOneStepTheta operator=(const TimIntOneStepTheta& old) = delete;

    /// copy constructor
    TimIntOneStepTheta(const TimIntOneStepTheta& old) = delete;

    void Init() override;

    void Setup() override;

    void PreSolve() override{};

    void PostSolve() override{};

    void print_time_step_info() override;

    void compute_intermediate_values() override{};

    void compute_interior_values() override{};

    void compute_time_derivative() override;

    void compute_time_deriv_pot0(const bool init) override{};

    void Update() override;

    void read_restart(
        const int step, Teuchos::RCP<CORE::IO::InputControl> input = Teuchos::null) override;

    Teuchos::RCP<Epetra_Vector> Phiaf() override { return Teuchos::null; }

    Teuchos::RCP<Epetra_Vector> Phiam() override { return Teuchos::null; }

    Teuchos::RCP<Epetra_Vector> Phidtam() override { return Teuchos::null; }

    Teuchos::RCP<Epetra_Vector> FsPhi() override
    {
      if (Sep_ != Teuchos::null) Sep_->Multiply(false, *phinp_, *fsphinp_);
      return fsphinp_;
    }

    Teuchos::RCP<Teuchos::ParameterList> scatra_time_parameter_list() override
    {
      Teuchos::RCP<Teuchos::ParameterList> timeparams;
      timeparams = Teuchos::rcp(new Teuchos::ParameterList());
      timeparams->set("using stationary formulation", false);
      timeparams->set("using generalized-alpha time integration", false);
      timeparams->set("total time", time_);
      timeparams->set("time factor", theta_ * dta_);
      timeparams->set("alpha_F", 1.0);
      return timeparams;
    }

    //! set state on micro scale in multi-scale simulations
    void set_state(Teuchos::RCP<Epetra_Vector> phin,  //!< micro-scale state vector at old time step
        Teuchos::RCP<Epetra_Vector> phinp,            //!< micro-scale state vector at new time step
        Teuchos::RCP<Epetra_Vector>
            phidtn,  //!< time derivative of micro-scale state vector at old time step
        Teuchos::RCP<Epetra_Vector>
            phidtnp,  //!< time derivative of micro-scale state vector at new time step
        Teuchos::RCP<Epetra_Vector> hist,                     //!< micro-scale history vector
        Teuchos::RCP<CORE::IO::DiscretizationWriter> output,  //!< micro-scale discretization writer
        const std::vector<double>&
            phinp_macro,   //!< values of state variables at macro-scale Gauss point
        const int step,    //!< time step
        const double time  //!< time
    );

    //! clear state on micro scale in multi-scale simulations
    void ClearState();

    void pre_calc_initial_time_derivative() override;

    void post_calc_initial_time_derivative() override;

    void write_restart() const override;

   protected:
    void set_element_time_parameter(bool forcedincrementalsolver = false) const override;

    void set_time_for_neumann_evaluation(Teuchos::ParameterList& params) override;

    void calc_initial_time_derivative() override;

    void set_old_part_of_righthandside() override;

    void explicit_predictor() const override;

    void add_neumann_to_residual() override;

    void av_m3_separation() override;

    void dynamic_computation_of_cs() override;

    void dynamic_computation_of_cv() override;

    void add_time_integration_specific_vectors(bool forcedincrementalsolver = false) override;

    double residual_scaling() const override { return 1.0 / (dta_ * theta_); }

    /// time factor for one-step-theta/BDF2 time integration
    double theta_;

    /// fine-scale solution vector at time n+1
    Teuchos::RCP<Epetra_Vector> fsphinp_;
  };

}  // namespace SCATRA

FOUR_C_NAMESPACE_CLOSE

#endif
