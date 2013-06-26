/*!----------------------------------------------------------------------
\file scatra_timint_bdf2.H
\brief

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>

*----------------------------------------------------------------------*/

#ifndef SCATRA_TIMINT_BDF2_H
#define SCATRA_TIMINT_BDF2_H


#include "scatra_timint_implicit.H"

#include "../linalg/linalg_utils.H"


namespace SCATRA
{

class TimIntBDF2 : public ScaTraTimIntImpl
{
public:
  /// Standard Constructor
  TimIntBDF2(Teuchos::RCP<DRT::Discretization> dis,
      Teuchos::RCP<LINALG::Solver>            solver,
      Teuchos::RCP<Teuchos::ParameterList>    params,
      Teuchos::RCP<Teuchos::ParameterList>    extraparams,
      Teuchos::RCP<IO::DiscretizationWriter>  output);

  /// Destructor
  virtual ~TimIntBDF2();

  /// predict thermodynamic pressure for low-Mach-number flow
  void PredictThermPressure();

  /// predict density field for low-Mach-number flow
  void PredictDensity();

  /// compute values at intermediate time steps (required for generalized-alpha)
  void ComputeIntermediateValues(){return;};

  /// compute values of thermodynamic pressure at intermediate time steps (required for generalized-alpha)
  void ComputeThermPressureIntermediateValues(){return;};

  /// compute thermodynamic pressure and time derivative for low-Mach-number flow
  void ComputeThermPressure();

  ///  compute scalar time derivative
  void ComputeTimeDerivative();

  ///  compute time derivative of thermodynamic pressure
  void ComputeThermPressureTimeDerivative();

  /// Update the solution after convergence of the nonlinear iteration.
  /// Current solution becomes old solution of next timestep.
  void Update(const int num=0);

  /// update thermodynamic pressure for low-Mach-number flow
  void UpdateThermPressure();

  /// update density field for ELCH natural convection
  void UpdateDensityElch();

  /// update solution after reinitialization for OST scheme
  void UpdateReinit() { dserror("only implemented for OST"); return; };

  /// read restart data
  void ReadRestart(int step);

  // routine to return scalar field phi at time step n-1
  Teuchos::RCP<Epetra_Vector> Phinm() { return Teuchos::null; }

  /// routine to return scalar field phi at time step n+alpha_F
  Teuchos::RCP<Epetra_Vector> Phiaf(){return Teuchos::null;}

  /// routine to return scalar field phi at time step n+alpha_M
  Teuchos::RCP<Epetra_Vector> Phiam(){return Teuchos::null;}

  /// routine to return time derivative of scalar field phi at time step n+alpha_M
  Teuchos::RCP<Epetra_Vector> Phidtam(){return Teuchos::null;}

  /// routine to return fine-scale scalar field fsphi at time step n+1
  virtual Teuchos::RCP<Epetra_Vector> FsPhi() {
    if (Sep_ != Teuchos::null)
    Sep_->Multiply(false,*phinp_,*fsphinp_);
    return fsphinp_; };

  /// routine to return thermo. press. at time step n+alpha_F for low-Mach-number flow
  double ThermPressAf(){return 1.0;}

  /// routine to return thermo. press. at time step n+alpha_M for low-Mach-number flow
  double ThermPressAm(){return 1.0;}

  /// routine to return time derivative of thermo. press. at time step n+alpha_F for low-Mach-number flow
  double ThermPressDtAf(){return 0.0;}

  /// routine to return time derivative of thermo. press. at time step n+alpha_M for low-Mach-number flow
  double ThermPressDtAm(){return 0.0;}

  /// routine to return time integration specific parameters
  Teuchos::RCP<Teuchos::ParameterList> ScatraTimeParameterList()
  {
    Teuchos::RCP<Teuchos::ParameterList> timeparams;
    timeparams = Teuchos::rcp(new Teuchos::ParameterList());
    timeparams->set("using stationary formulation",false);
    timeparams->set("using generalized-alpha time integration",false);
    timeparams->set("total time",time_);
    timeparams->set("time factor",theta_*dta_);
    timeparams->set("alpha_F",1.0);
    return timeparams;
  }


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
  void SetTimeForNeumannEvaluation(Teuchos::ParameterList& params);

  /// add actual Neumann loads with time factor
  void AddNeumannToResidual();

  /// AVM3-based scale separation
  void AVM3Separation();

  /// dynamic Smagorinsky model
  void DynamicComputationOfCs();

  /// add parameters specific for time-integration scheme
  void AddSpecificTimeIntegrationParameters(Teuchos::ParameterList& params);

  /// write additional data required for restart
  void OutputRestart();

  /// return the right time-scaling-factor for the true residual
  double ResidualScaling() const { return 1.0/(dta_*theta_); }

  /// update time-dependent electrode state variables at the end of an time step
  void ElectrodeKineticsTimeUpdate(const bool init = false);

  /// set old part of RHS for galvanostatic equation
  void ElectrodeKineticsSetOldPartOfRHS();

private:

  /// solution vector phi at time n-1
  Teuchos::RCP<Epetra_Vector>  phinm_;

  /// fine-scale solution vector at time n+1
  Teuchos::RCP<Epetra_Vector>  fsphinp_;

  /// time factor for one-step-theta/BDF2 time integration
  double theta_;

  /// LOMA-specific parameters: thermodynamic pressure at n-1
  double thermpressnm_;

}; // class TimIntBDF2

} // namespace SCATRA


#endif  // #ifndef SCATRA_TIMINT_BDF2_H
