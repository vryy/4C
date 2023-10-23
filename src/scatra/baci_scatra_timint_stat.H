/*----------------------------------------------------------------------*/
/*! \file
\brief solution algorithm for stationary problems

\level 1


*----------------------------------------------------------------------*/

#ifndef BACI_SCATRA_TIMINT_STAT_H
#define BACI_SCATRA_TIMINT_STAT_H

#include "baci_linalg_utils_sparse_algebra_create.H"
#include "baci_scatra_timint_implicit.H"

namespace SCATRA
{
  class TimIntStationary : public virtual ScaTraTimIntImpl
  {
   public:
    /// Standard Constructor
    TimIntStationary(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);


    /// initialize time integration scheme
    void Init() override;

    /// compute values at intermediate time steps (required for generalized-alpha)
    void ComputeIntermediateValues() override { return; };

    /// compute values at the interior of the elements (required for hdg)
    void ComputeInteriorValues() override { return; };

    ///  compute scalar time derivate parameters of the input voltage
    void ComputeTimeDerivPot0(const bool init) override { return; };

    void Setup() override;

    /// Update the solution after convergence of the nonlinear iteration.
    /// Current solution becomes old solution of next timestep.
    void Update(const int num = 0) override;

    /// read restart data
    void ReadRestart(const int step, Teuchos::RCP<IO::InputControl> input = Teuchos::null) override;

    // routine to return scalar field phi at time step n-1
    Teuchos::RCP<Epetra_Vector> Phinm() { return Teuchos::null; }

    /// routine to return scalar field phi at time step n+alpha_F
    Teuchos::RCP<Epetra_Vector> Phiaf() override { return Teuchos::null; }

    /// routine to return scalar field phi at time step n+alpha_M
    Teuchos::RCP<Epetra_Vector> Phiam() override { return Teuchos::null; }

    /// routine to return time derivative of scalar field phi at time step n+alpha_M
    Teuchos::RCP<Epetra_Vector> Phidtam() override { return Teuchos::null; }

    /// routine to return fine-scale scalar field fsphi
    Teuchos::RCP<Epetra_Vector> FsPhi() override
    {
      if (Sep_ != Teuchos::null) Sep_->Multiply(false, *phinp_, *fsphinp_);
      return fsphinp_;
    };

    /// routine to return time integration specific parameters
    Teuchos::RCP<Teuchos::ParameterList> ScatraTimeParameterList() override
    {
      dserror("Not yet implemented!");
      return Teuchos::null;
    }


   protected:
    /// don't want = operator and cctor
    TimIntStationary operator=(const TimIntStationary& old);

    /// copy constructor
    TimIntStationary(const TimIntStationary& old);

    /// set time parameter for element evaluation
    void SetElementTimeParameter(bool forcedincrementalsolver = false) const override;

    //! set time for evaluation of Neumann boundary conditions
    void SetTimeForNeumannEvaluation(Teuchos::ParameterList& params) override;

    //! calculate consistent initial conditions in compliance with initial scalar field
    //! this is not necessary for stationary calculations
    void CalcInitialTimeDerivative() override { return; };

    /// Set the part of the righthandside belonging to the last timestep.
    void SetOldPartOfRighthandside() override;

    /// do explicit predictor step (nothing to predict for stationary problems!)
    void ExplicitPredictor() const override { return; };

    /// add actual Neumann loads with time factor
    void AddNeumannToResidual() override;

    /// AVM3-based scale separation
    void AVM3Separation() override;

    /// add parameters specific for time-integration scheme
    void AddTimeIntegrationSpecificVectors(bool forcedincrementalsolver = false) override;

    /// dynamic Smagorinsky model
    void DynamicComputationOfCs() override
    {
      dserror("no turbulence in stationary flows!");
      return;
    };

    /// dynamic Vreman model
    void DynamicComputationOfCv() override
    {
      dserror("no turbulence in stationary flows!");
      return;
    };

    /// write additional data required for restart
    void OutputRestart() const override;

    /// return the right time-scaling-factor for the true residual
    double ResidualScaling() const override { return 1.0; }

   private:
    /// fine-scale solution vector at time n+1
    Teuchos::RCP<Epetra_Vector> fsphinp_;


  };  // class TimIntStationary

}  // namespace SCATRA

#endif  // SCATRA_TIMINT_STAT_H
