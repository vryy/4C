/*!----------------------------------------------------------------------
\file scatra_timint_bdf2.H
\brief

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>

*----------------------------------------------------------------------*/

#ifdef CCADISCRET
#ifndef SCATRA_TIMINT_BDF2_H
#define SCATRA_TIMINT_BDF2_H


#include "scatra_timint_implicit.H"


namespace SCATRA
{

class TimIntBDF2 : public ScaTraTimIntImpl
{
public:
  /// Standard Constructor
  TimIntBDF2(RCP<DRT::Discretization>       dis,
             RCP<LINALG::Solver>            solver,
             RCP<ParameterList>             params,
             RCP<IO::DiscretizationWriter>  output);

  /// Destructor
  virtual ~TimIntBDF2();

  /// predict thermodynamic pressure for low-Mach-number flow
  void PredictThermPressure();

  /// predict density field for low-Mach-number flow
  void PredictDensity();

  /// compute values at intermediate time steps (required for generalized-alpha)
  void ComputeIntermediateValues(){return;};

  /// compute thermodynamic pressure and time derivative for low-Mach-number flow
  void ComputeThermPressure();

  /// Update the solution after convergence of the nonlinear iteration.
  /// Current solution becomes old solution of next timestep.
  void Update();

  /// update thermodynamic pressure for low-Mach-number flow
  void UpdateThermPressure();

  /// update density field for ELCH natural convection
  void UpdateDensityElch();

  /// read restart data
  void ReadRestart(int step);


protected:

  /// don't want = operator and cctor
  TimIntBDF2 operator = (const TimIntBDF2& old);

  /// copy constructor
  TimIntBDF2 (const TimIntBDF2& old);

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
  double ResidualScaling() const { return 1.0/(dta_*theta_); }

private:

  /// solution vector phi at time n-1
  RCP<Epetra_Vector>    phinm_;

  /// fine-scale solution vector at time n+1
  RCP<Epetra_Vector>  fsphinp_;

  /// time factor for one-step-theta/BDF2 time integration
  double theta_;

  /// LOMA-specific parameters: thermodynamic pressure at n-1
  double thermpressnm_;

}; // class TimIntBDF2

} // namespace SCATRA


#endif  // #ifndef SCATRA_TIMINT_BDF2_H
#endif  // #ifdef CCADISCRET
