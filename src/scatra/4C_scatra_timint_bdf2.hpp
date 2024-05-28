/*----------------------------------------------------------------------*/
/*! \file
\brief BDF2 time-integration scheme

\level 1



*----------------------------------------------------------------------*/

#ifndef FOUR_C_SCATRA_TIMINT_BDF2_HPP
#define FOUR_C_SCATRA_TIMINT_BDF2_HPP

#include "4C_config.hpp"

#include "4C_linalg_sparsematrix.hpp"
#include "4C_scatra_timint_implicit.hpp"

FOUR_C_NAMESPACE_OPEN

namespace SCATRA
{
  class TimIntBDF2 : public virtual ScaTraTimIntImpl
  {
   public:
    /// Standard Constructor
    TimIntBDF2(Teuchos::RCP<DRT::Discretization> dis, Teuchos::RCP<CORE::LINALG::Solver> solver,
        Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);

    /// Setup time integration scheme
    void Setup() override;

    /// compute values at intermediate time steps (required for generalized-alpha)
    void compute_intermediate_values() override{};

    /// compute values at the interior of the elements (required for hdg)
    void compute_interior_values() override{};

    ///  compute scalar time derivative
    void compute_time_derivative() override;

    ///  compute scalar time derivate parameters of the input voltage to compute double layer
    ///  current densities
    void compute_time_deriv_pot0(const bool init) override{};

    /// Update the solution after convergence of the nonlinear iteration.
    /// Current solution becomes old solution of next timestep.
    void Update() override;

    /// read restart data
    void read_restart(
        const int step, Teuchos::RCP<IO::InputControl> input = Teuchos::null) override;

    /// routine to return scalar field phi at time step n+alpha_F
    Teuchos::RCP<Epetra_Vector> Phiaf() override { return Teuchos::null; }

    /// routine to return scalar field phi at time step n+alpha_M
    Teuchos::RCP<Epetra_Vector> Phiam() override { return Teuchos::null; }

    /// routine to return time derivative of scalar field phi at time step n+alpha_M
    Teuchos::RCP<Epetra_Vector> Phidtam() override { return Teuchos::null; }

    /// routine to return fine-scale scalar field fsphi at time step n+1
    Teuchos::RCP<Epetra_Vector> FsPhi() override
    {
      if (Sep_ != Teuchos::null) Sep_->Multiply(false, *phinp_, *fsphinp_);
      return fsphinp_;
    };

    /// routine to return time integration specific parameters
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

    /// don't want = operator and cctor
    TimIntBDF2 operator=(const TimIntBDF2& old) = delete;

    /// copy constructor
    TimIntBDF2(const TimIntBDF2& old) = delete;

   protected:
    /// set time parameter for element evaluation
    void set_element_time_parameter(bool forcedincrementalsolver = false) const override;

    //! set time for evaluation of Neumann boundary conditions
    void set_time_for_neumann_evaluation(Teuchos::ParameterList& params) override;

    //! calculate consistent initial conditions in compliance with initial scalar field
    //! this is not necessary for BDF2 time integration scheme
    void calc_initial_time_derivative() override{};

    /// Set the part of the righthandside belonging to the last timestep.
    void set_old_part_of_righthandside() override;

    /// do explicit predictor step (-> better starting value for nonlinear solver)
    void explicit_predictor() const override;

    /// add actual Neumann loads with time factor
    void add_neumann_to_residual() override;

    /// AVM3-based scale separation
    void av_m3_separation() override;

    /// dynamic Smagorinsky model
    void dynamic_computation_of_cs() override;

    /// dynamic Vreman model
    void dynamic_computation_of_cv() override;

    /// add parameters specific for time-integration scheme
    void add_time_integration_specific_vectors(bool forcedincrementalsolver = false) override;

    void write_restart() const override;

    /// return the right time-scaling-factor for the true residual
    double residual_scaling() const override { return 1.0 / (dta_ * theta_); }

    /// time factor for one-step-theta/BDF2 time integration
    double theta_;

    /// solution vector phi at time n-1
    Teuchos::RCP<Epetra_Vector> phinm_;

    /// fine-scale solution vector at time n+1
    Teuchos::RCP<Epetra_Vector> fsphinp_;
  };  // class TimIntBDF2
}  // namespace SCATRA
FOUR_C_NAMESPACE_CLOSE

#endif
