/*!----------------------------------------------------------------------
\file scatra_timint_genalpha.H
\brief

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>

*----------------------------------------------------------------------*/

#ifdef CCADISCRET
#ifndef SCATRA_TIMINT_GENALPHA_H
#define SCATRA_TIMINT_GENALPHA_H

#include "scatra_timint_implicit.H"


namespace SCATRA
{

class TimIntGenAlpha : public ScaTraTimIntImpl
{

public:

  /// Standard Constructor
  TimIntGenAlpha(RCP<DRT::Discretization>       dis,
                 RCP<LINALG::Solver>            solver,
                 RCP<ParameterList>             params,
                 RCP<ParameterList>             extraparams,
                 RCP<IO::DiscretizationWriter>  output);

  /// Destructor
  virtual ~TimIntGenAlpha();

  /// predict thermodynamic pressure and time derivative for low-Mach-number flow
  void PredictThermPressure();

  /// compute values at intermediate time steps
  void ComputeIntermediateValues();

  /// compute values of thermodynamic pressure at intermediate time steps
  void ComputeThermPressureIntermediateValues();

  /// compute thermodynamic pressure and time derivative for low-Mach-number flow
  void ComputeThermPressure();

  ///  update time derivative after solution
  void UpdateTimeDerivative();

  ///  update time derivative of thermodynamic pressure after solution
  void UpdateThermPressureTimeDerivative();

  /// Update the solution after convergence of the nonlinear iteration.
  /// Current solution becomes old solution of next timestep.
  void Update();

  /// update thermodynamic pressure and time derivative for low-Mach-number flow
  void UpdateThermPressure();

  /// update density field for ELCH natural convection
  void UpdateDensityElch(){dserror("not implemented"); return;};

  /// read restart data
  void ReadRestart(int step);

  /// routine to return scalar field phi at time step n+alpha_F
  Teuchos::RCP<Epetra_Vector> Phiaf() { return phiaf_; }

  /// routine to return scalar field phi at time step n+alpha_M
  Teuchos::RCP<Epetra_Vector> Phiam() { return phiam_; }

  /// routine to return time derivative of scalar field phi at time step n+alpha_M
  Teuchos::RCP<Epetra_Vector> Phidtam() { return phidtam_; }

  /// routine to return thermo. press. at time step n+alpha_F for low-Mach-number flow
  double ThermPressAf() { return thermpressaf_; }

  /// routine to return thermo. press. at time step n+alpha_M for low-Mach-number flow
  double ThermPressAm() { return thermpressam_; }

  /// routine to return thermo. press. at time step n+alpha_M for low-Mach-number flow
  double ThermPressDtAm() { return thermpressdtam_; }


protected:

  /// don't want = operator and cctor
  TimIntGenAlpha operator = (const TimIntGenAlpha& old);

  /// copy constructor
  TimIntGenAlpha (const TimIntGenAlpha& old);

  /// Initialization procedure before the first time step is done
  void PrepareFirstTimeStep();

  /// Set the part of the righthandside belonging to the last timestep.
  void SetOldPartOfRighthandside();

  /// do explicit predictor step (-> better starting value for nonlinear solver)
  void ExplicitPredictor();

  /// set time for evaluation of Neumann boundary conditions
  void SetTimeForNeumannEvaluation(ParameterList& params);

  /// add actual Neumann loads with time factor
  void AddNeumannToResidual();

  /// AVM3-based scale separation
  void AVM3Separation();

  /// add parameters specific for time-integration scheme
  void AddSpecificTimeIntegrationParameters(ParameterList& params);

  /// write additional data required for restart
  void OutputRestart();

  /// return the right time-scaling-factor for the true residual
  double ResidualScaling() const { return 1.0/(dta_*genalphafac_); }


private:

  /// scalar at time n+alpha_F and n+alpha_M
  RCP<Epetra_Vector>  phiaf_;
  RCP<Epetra_Vector>  phiam_;

  /// scalar time derivative at time n+alpha_M
  RCP<Epetra_Vector>  phidtam_;

  /// fine-scale part at time n+alpha_F
  RCP<Epetra_Vector>  fsphiaf_;

  /// time factors for generalized-alpha time integration
  double alphaM_;
  double alphaF_;
  double gamma_;
  double genalphafac_;

  /// LOMA-specific parameter: thermodynamic pressure at n+alpha_F and n+alpha_M
  /// and its time derivative at n+alpha_F and n+alpha_M
  double thermpressaf_;
  double thermpressam_;
  double thermpressdtaf_;
  double thermpressdtam_;

}; // class TimIntGenAlpha

} // namespace SCATRA

#endif  // #ifndef SCATRA_TIMINT_GENALPHA_H
#endif  // #ifdef CCADISCRET
