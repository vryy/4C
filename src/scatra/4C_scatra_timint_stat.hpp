/*----------------------------------------------------------------------*/
/*! \file
\brief solution algorithm for stationary problems

\level 1


*----------------------------------------------------------------------*/

#ifndef FOUR_C_SCATRA_TIMINT_STAT_HPP
#define FOUR_C_SCATRA_TIMINT_STAT_HPP

#include "4C_config.hpp"

#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_scatra_timint_implicit.hpp"

FOUR_C_NAMESPACE_OPEN

namespace ScaTra
{
  class TimIntStationary : public virtual ScaTraTimIntImpl
  {
   public:
    /// Standard Constructor
    TimIntStationary(Teuchos::RCP<Core::FE::Discretization> dis,
        Teuchos::RCP<Core::LinAlg::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<Core::IO::DiscretizationWriter> output);


    /// initialize time integration scheme
    void init() override;

    /// compute values at intermediate time steps (required for generalized-alpha)
    void compute_intermediate_values() override { return; };

    /// compute values at the interior of the elements (required for hdg)
    void compute_interior_values() override { return; };

    ///  compute scalar time derivate parameters of the input voltage
    void compute_time_deriv_pot0(const bool init) override { return; };

    void setup() override;

    /// Update the solution after convergence of the nonlinear iteration.
    /// Current solution becomes old solution of next timestep.
    void update() override;

    /// read restart data
    void read_restart(
        const int step, Teuchos::RCP<Core::IO::InputControl> input = Teuchos::null) override;

    // routine to return scalar field phi at time step n-1
    Teuchos::RCP<Epetra_Vector> phinm() { return Teuchos::null; }

    /// routine to return scalar field phi at time step n+alpha_F
    Teuchos::RCP<Epetra_Vector> phiaf() override { return Teuchos::null; }

    /// routine to return scalar field phi at time step n+alpha_M
    Teuchos::RCP<Epetra_Vector> phiam() override { return Teuchos::null; }

    /// routine to return time derivative of scalar field phi at time step n+alpha_M
    Teuchos::RCP<Epetra_Vector> phidtam() override { return Teuchos::null; }

    /// routine to return fine-scale scalar field fsphi
    Teuchos::RCP<Epetra_Vector> fs_phi() override
    {
      if (Sep_ != Teuchos::null) Sep_->multiply(false, *phinp_, *fsphinp_);
      return fsphinp_;
    };

    /// routine to return time integration specific parameters
    Teuchos::RCP<Teuchos::ParameterList> scatra_time_parameter_list() override
    {
      FOUR_C_THROW("Not yet implemented!");
      return Teuchos::null;
    }


   protected:
    /// don't want = operator and cctor
    TimIntStationary operator=(const TimIntStationary& old);

    /// copy constructor
    TimIntStationary(const TimIntStationary& old);

    /// set time parameter for element evaluation
    void set_element_time_parameter(bool forcedincrementalsolver = false) const override;

    //! set time for evaluation of Neumann boundary conditions
    void set_time_for_neumann_evaluation(Teuchos::ParameterList& params) override;

    //! calculate consistent initial conditions in compliance with initial scalar field
    //! this is not necessary for stationary calculations
    void calc_initial_time_derivative() override { return; };

    /// Set the part of the righthandside belonging to the last timestep.
    void set_old_part_of_righthandside() override;

    /// do explicit predictor step (nothing to predict for stationary problems!)
    void explicit_predictor() const override { return; };

    /// add actual Neumann loads with time factor
    void add_neumann_to_residual() override;

    /// AVM3-based scale separation
    void av_m3_separation() override;

    /// add parameters specific for time-integration scheme
    void add_time_integration_specific_vectors(bool forcedincrementalsolver = false) override;

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

    void write_restart() const override;

    /// return the right time-scaling-factor for the true residual
    double residual_scaling() const override { return 1.0; }

   private:
    /// fine-scale solution vector at time n+1
    Teuchos::RCP<Epetra_Vector> fsphinp_;


  };  // class TimIntStationary

}  // namespace ScaTra

FOUR_C_NAMESPACE_CLOSE

#endif
