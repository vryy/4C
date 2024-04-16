/*----------------------------------------------------------------------*/
/*! \file
\brief BDF2 time-integration scheme

\level 1



*----------------------------------------------------------------------*/

#ifndef FOUR_C_SCATRA_TIMINT_BDF2_HPP
#define FOUR_C_SCATRA_TIMINT_BDF2_HPP

#include "baci_config.hpp"

#include "baci_linalg_sparsematrix.hpp"
#include "baci_scatra_timint_implicit.hpp"

BACI_NAMESPACE_OPEN

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
    void ComputeIntermediateValues() override{};

    /// compute values at the interior of the elements (required for hdg)
    void ComputeInteriorValues() override{};

    ///  compute scalar time derivative
    void ComputeTimeDerivative() override;

    ///  compute scalar time derivate parameters of the input voltage to compute double layer
    ///  current densities
    void ComputeTimeDerivPot0(const bool init) override{};

    /// Update the solution after convergence of the nonlinear iteration.
    /// Current solution becomes old solution of next timestep.
    void Update() override;

    /// read restart data
    void ReadRestart(const int step, Teuchos::RCP<IO::InputControl> input = Teuchos::null) override;

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
    Teuchos::RCP<Teuchos::ParameterList> ScatraTimeParameterList() override
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
    void SetElementTimeParameter(bool forcedincrementalsolver = false) const override;

    //! set time for evaluation of Neumann boundary conditions
    void SetTimeForNeumannEvaluation(Teuchos::ParameterList& params) override;

    //! calculate consistent initial conditions in compliance with initial scalar field
    //! this is not necessary for BDF2 time integration scheme
    void CalcInitialTimeDerivative() override{};

    /// Set the part of the righthandside belonging to the last timestep.
    void SetOldPartOfRighthandside() override;

    /// do explicit predictor step (-> better starting value for nonlinear solver)
    void ExplicitPredictor() const override;

    /// add actual Neumann loads with time factor
    void AddNeumannToResidual() override;

    /// AVM3-based scale separation
    void AVM3Separation() override;

    /// dynamic Smagorinsky model
    void DynamicComputationOfCs() override;

    /// dynamic Vreman model
    void DynamicComputationOfCv() override;

    /// add parameters specific for time-integration scheme
    void AddTimeIntegrationSpecificVectors(bool forcedincrementalsolver = false) override;

    void WriteRestart() const override;

    /// return the right time-scaling-factor for the true residual
    double ResidualScaling() const override { return 1.0 / (dta_ * theta_); }

    /// time factor for one-step-theta/BDF2 time integration
    double theta_;

    /// solution vector phi at time n-1
    Teuchos::RCP<Epetra_Vector> phinm_;

    /// fine-scale solution vector at time n+1
    Teuchos::RCP<Epetra_Vector> fsphinp_;
  };  // class TimIntBDF2
}  // namespace SCATRA
BACI_NAMESPACE_CLOSE

#endif
