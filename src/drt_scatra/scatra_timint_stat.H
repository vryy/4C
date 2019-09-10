/*----------------------------------------------------------------------*/
/*! \file
\brief solution algorithm for stationary problems

\level 1

\maintainer Anh-Tu Vuong

*----------------------------------------------------------------------*/

#ifndef SCATRA_TIMINT_STAT_H
#define SCATRA_TIMINT_STAT_H

#include "scatra_timint_implicit.H"

#include "../linalg/linalg_utils.H"

namespace SCATRA
{
  class TimIntStationary : public virtual ScaTraTimIntImpl
  {
   public:
    /// Standard Constructor
    TimIntStationary(Teuchos::RCP<DRT::Discretization> dis, Teuchos::RCP<LINALG::Solver> solver,
        Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);

    /// Destructor
    virtual ~TimIntStationary();

    /// initialize time integration scheme
    virtual void Init();

    /// compute values at intermediate time steps (required for generalized-alpha)
    void ComputeIntermediateValues() { return; };

    /// compute values at the interior of the elements (required for hdg)
    virtual void ComputeInteriorValues() { return; };

    ///  compute scalar time derivate parameters of the input voltage
    void ComputeTimeDerivPot0(const bool init) { return; };

    /// Update the solution after convergence of the nonlinear iteration.
    /// Current solution becomes old solution of next timestep.
    virtual void Update(const int num = 0);

    /// read restart data
    virtual void ReadRestart(const int step, Teuchos::RCP<IO::InputControl> input = Teuchos::null);

    // routine to return scalar field phi at time step n-1
    Teuchos::RCP<Epetra_Vector> Phinm() { return Teuchos::null; }

    /// routine to return scalar field phi at time step n+alpha_F
    Teuchos::RCP<Epetra_Vector> Phiaf() { return Teuchos::null; }

    /// routine to return scalar field phi at time step n+alpha_M
    Teuchos::RCP<Epetra_Vector> Phiam() { return Teuchos::null; }

    /// routine to return time derivative of scalar field phi at time step n+alpha_M
    Teuchos::RCP<Epetra_Vector> Phidtam() { return Teuchos::null; }

    /// routine to return fine-scale scalar field fsphi
    Teuchos::RCP<Epetra_Vector> FsPhi()
    {
      if (Sep_ != Teuchos::null) Sep_->Multiply(false, *phinp_, *fsphinp_);
      return fsphinp_;
    };

    /// routine to return time integration specific parameters
    Teuchos::RCP<Teuchos::ParameterList> ScatraTimeParameterList()
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
    void SetElementTimeParameter(bool forcedincrementalsolver = false) const;

    //! set time for evaluation of Neumann boundary conditions
    void SetTimeForNeumannEvaluation(Teuchos::ParameterList& params);

    //! calculate consistent initial conditions in compliance with initial scalar field
    //! this is not necessary for stationary calculations
    void CalcInitialTimeDerivative() { return; };

    /// Set the part of the righthandside belonging to the last timestep.
    void SetOldPartOfRighthandside();

    /// do explicit predictor step (nothing to predict for stationary problems!)
    void ExplicitPredictor() const { return; };

    /// add actual Neumann loads with time factor
    void AddNeumannToResidual();

    /// AVM3-based scale separation
    void AVM3Separation();

    /// add parameters specific for time-integration scheme
    void AddTimeIntegrationSpecificVectors(bool forcedincrementalsolver = false);

    /// dynamic Smagorinsky model
    void DynamicComputationOfCs()
    {
      dserror("no turbulence in stationary flows!");
      return;
    };

    /// dynamic Vreman model
    void DynamicComputationOfCv()
    {
      dserror("no turbulence in stationary flows!");
      return;
    };

    /// write additional data required for restart
    virtual void OutputRestart() const;

    /// return the right time-scaling-factor for the true residual
    double ResidualScaling() const { return 1.0; }

   private:
    /// fine-scale solution vector at time n+1
    Teuchos::RCP<Epetra_Vector> fsphinp_;


  };  // class TimIntStationary

}  // namespace SCATRA

#endif  // #ifndef SCATRA_TIMINT_STAT_H
