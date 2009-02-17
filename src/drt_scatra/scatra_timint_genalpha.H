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
                 RCP<IO::DiscretizationWriter>  output);

  /// Destructor
  virtual ~TimIntGenAlpha();

  /// set velocity field for low-Mach-number flow
  void SetLomaVelocity(RCP<const Epetra_Vector> extvel,
                       RCP<DRT::Discretization> fluiddis);

  /// compute initial time derivative of density field for low-Mach-number flow
  void ComputeInitialDensityDerivative();

  /// predict thermodynamic pressure and time derivative for low-Mach-number flow
  void PredictThermPressure();

  /// predict density field and time derivative for low-Mach-number flow
  void PredictDensity();

  /// compute values at intermediate time steps
  void ComputeIntermediateValues();

  /// compute thermodynamic pressure and time derivative for low-Mach-number flow
  double ComputeThermPressure();

  /// compute time derivative of density field for low-Mach-number flow
  void ComputeDensityDerivative();

  /// Update the solution after convergence of the nonlinear iteration.
  /// Current solution becomes old solution of next timestep.
  void Update();

  /// update thermodynamic pressure and time derivative for low-Mach-number flow
  void UpdateThermPressure();

  /// update density field and time derivative for low-Mach-number flow
  void UpdateDensity();

  /// read restart data
  void ReadRestart(int step);


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

  /// reset the residual vector and add actual Neumann loads with time factor
  void AddNeumannToResidual();

  /// add parameters specific for time-integration scheme
  void AddSpecificTimeIntegrationParameters(ParameterList& params);

  /// write additional data required for restart
  void OutputRestart();

  /// time derivative of solution vector phi at times t_{n+1} and t_{n}
  RCP<Epetra_Vector>  phidtnp_;
  RCP<Epetra_Vector>  phidtn_;


private:

  /// calculate initial value for phidtn_
  void CalcInitialPhidt();

  /// solution vector at time n+alpha_F
  RCP<Epetra_Vector>  phiaf_;

  /// density at times n+alpha_M and n+alpha_F
  RCP<Epetra_Vector>  densam_;
  RCP<Epetra_Vector>  densaf_;

  /// time factors for generalized-alpha time integration
  double alphaM_;
  double alphaF_;
  double gamma_;
  double genalphafac_;

  /// LOMA-specific parameter: time derivative of thermodynamic pressure
  /// at n+alpha_F
  double thermpressdtaf_;

}; // class TimIntGenAlpha

} // namespace SCATRA

#endif  // #ifndef SCATRA_TIMINT_GENALPHA_H
#endif  // #ifdef CCADISCRET
