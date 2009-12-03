/*!----------------------------------------------------------------------
\file scatra_timint_stat.H
\brief

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>

*----------------------------------------------------------------------*/

#ifdef CCADISCRET
#ifndef SCATRA_TIMINT_STAT_H
#define SCATRA_TIMINT_STAT_H

#include "scatra_timint_implicit.H"


namespace SCATRA
{

class TimIntStationary : public ScaTraTimIntImpl
{

public:

  /// Standard Constructor
  TimIntStationary(RCP<DRT::Discretization>      dis,
                  RCP<LINALG::Solver>            solver,
                  RCP<ParameterList>             params,
                  RCP<ParameterList>             extraparams,
                  RCP<IO::DiscretizationWriter>  output);

  /// Destructor
  virtual ~TimIntStationary();

  /// predict thermodynamic pressure and time derivative for low-Mach-number flow
  void PredictThermPressure(){return;};

  /// compute values at intermediate time steps (required for generalized-alpha)
  void ComputeIntermediateValues(){return;};

  /// compute values of thermodynamic pressure at intermediate time steps (required for generalized-alpha)
  void ComputeThermPressureIntermediateValues(){return;};

  /// compute thermodynamic pressure and time derivative for low-Mach-number flow
  void ComputeThermPressure(){return;};

  ///  update time derivative after solution
  void UpdateTimeDerivative(){return;};

  ///  update time derivative of thermodynamic pressure after solution
  void UpdateThermPressureTimeDerivative(){return;};

  /// Update the solution after convergence of the nonlinear iteration.
  /// Current solution becomes old solution of next timestep.
  void Update(){return;};

  /// update thermodynamic pressure and time derivative for low-Mach-number flow
  void UpdateThermPressure(){return;};

  /// update density field for ELCH natural convection
  void UpdateDensityElch(){return;};

  /// read restart data
  void ReadRestart(int step);

  /// routine to return scalar field phi at time step n+alpha_F
  Teuchos::RCP<Epetra_Vector> Phiaf(){return null;}

  /// routine to return scalar field phi at time step n+alpha_M
  Teuchos::RCP<Epetra_Vector> Phiam(){return null;}

  /// routine to return time derivative of scalar field phi at time step n+alpha_M
  Teuchos::RCP<Epetra_Vector> Phidtam(){return null;}

  /// routine to return thermo. press. at time step n+alpha_F for low-Mach-number flow
  double ThermPressAf(){return 1.0;}

  /// routine to return thermo. press. at time step n+alpha_M for low-Mach-number flow
  double ThermPressAm(){return 1.0;}

  /// routine to return thermo. press. at time step n+alpha_M for low-Mach-number flow
  double ThermPressDtAm(){return 0.0;}


protected:

  /// don't want = operator and cctor
  TimIntStationary operator = (const TimIntStationary& old);

  /// copy constructor
  TimIntStationary (const TimIntStationary& old);

  /// Initialization procedure before the first time step is done
  void PrepareFirstTimeStep(){return;};

  /// Set the part of the righthandside belonging to the last timestep.
  void SetOldPartOfRighthandside();

  /// do explicit predictor step (nothing to predict for stationary problems!)
  void ExplicitPredictor(){return;};

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
  double ResidualScaling() const { return 1.0; }

  /// update time-dependent electrode state variables at the end of an time step
  void ElectrodeKineticsTimeUpdate(const bool init = false){return;};

private:

  /// fine-scale solution vector at time n+1
  RCP<Epetra_Vector>  fsphinp_;


}; // class TimIntStationary

} // namespace SCATRA

#endif  // #ifndef SCATRA_TIMINT_STAT_H
#endif  // #ifdef CCADISCRET
