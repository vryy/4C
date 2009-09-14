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

  /// predict thermodynamic pressure and time derivative for low-Mach-number flow
  void PredictThermPressure();

  /// compute values at intermediate time steps
  void ComputeIntermediateValues();

  /// compute thermodynamic pressure and time derivative for low-Mach-number flow
  void ComputeThermPressure();

  /// Update the solution after convergence of the nonlinear iteration.
  /// Current solution becomes old solution of next timestep.
  void Update();

  /// update thermodynamic pressure and time derivative for low-Mach-number flow
  void UpdateThermPressure();

  /// update density field for ELCH natural convection
  void UpdateDensityElch(){dserror("not implemented"); return;};

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

  /// solution vectors at time n+alpha_F and n+alpha_M
  RCP<Epetra_Vector>  phiaf_;
  RCP<Epetra_Vector>  phiam_;

  /// fine-scale solution vector at time n+alpha_F
  RCP<Epetra_Vector>  fsphiaf_;

  /// time factors for generalized-alpha time integration
  double alphaM_;
  double alphaF_;
  double gamma_;
  double genalphafac_;

  /// LOMA-specific parameter: thermodynamic pressure at n+alpha_F and n+alpha_M
  /// and its time derivative at n+alpha_F
  double thermpressaf_;
  double thermpressam_;
  double thermpressdtaf_;

}; // class TimIntGenAlpha

} // namespace SCATRA

#endif  // #ifndef SCATRA_TIMINT_GENALPHA_H
#endif  // #ifdef CCADISCRET
