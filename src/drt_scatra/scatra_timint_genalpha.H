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

#ifndef SCATRA_TIMINT_GENALPHA_H
#define SCATRA_TIMINT_GENALPHA_H

#include "scatra_timint_implicit.H"


namespace SCATRA
{

class TimIntGenAlpha : public ScaTraTimIntImpl
{

public:

  /// Standard Constructor
  TimIntGenAlpha(Teuchos::RCP<DRT::Discretization>       dis,
      Teuchos::RCP<LINALG::Solver>            solver,
      Teuchos::RCP<Teuchos::ParameterList>    params,
      Teuchos::RCP<Teuchos::ParameterList>    extraparams,
      Teuchos::RCP<IO::DiscretizationWriter>  output);

  /// Destructor
  virtual ~TimIntGenAlpha();

  /// Print information about current time step to screen (reimplementation for Gen. Alpha)
  virtual void PrintTimeStepInfo()
  {
    if (myrank_==0)
        printf("TIME: %11.4E/%11.4E  DT = %11.4E  %s(a_F=%3.2f | a_M=%3.2f | gamma=%3.2f) STEP = %4d/%4d\n",
               time_,maxtime_,dta_,MethodTitle().c_str(),alphaF_,alphaM_,gamma_,step_,stepmax_);
    return;
  }

  /// predict thermodynamic pressure and time derivative for low-Mach-number flow
  void PredictThermPressure();

  /// compute values at intermediate time steps
  void ComputeIntermediateValues();

  /// compute values of thermodynamic pressure at intermediate time steps
  void ComputeThermPressureIntermediateValues();

  /// compute thermodynamic pressure and time derivative for low-Mach-number flow
  void ComputeThermPressure();

  ///  compute scalar time derivative
  void ComputeTimeDerivative();

  ///  compute time derivative of thermodynamic pressure
  void ComputeThermPressureTimeDerivative();

  /// Update the solution after convergence of the nonlinear iteration.
  /// Current solution becomes old solution of next timestep.
  void Update();

  /// update thermodynamic pressure and time derivative for low-Mach-number flow
  void UpdateThermPressure();

  /// update density field for ELCH natural convection
  void UpdateDensityElch(){dserror("not implemented"); return;};

  /// update solution after reinitialization for OST scheme
  void UpdateReinit();

  /// read restart data
  void ReadRestart(int step);

  // routine to return scalar field phi at time step n-1
  Teuchos::RCP<Epetra_Vector> Phinm() { return Teuchos::null; }

  /// routine to return scalar field phi at time step n+alpha_F
  Teuchos::RCP<Epetra_Vector> Phiaf() { return phiaf_; }

  /// routine to return scalar field phi at time step n+alpha_M
  Teuchos::RCP<Epetra_Vector> Phiam() { return phiam_; }

  /// routine to return time derivative of scalar field phi at time step n+alpha_M
  Teuchos::RCP<Epetra_Vector> Phidtam() { return phidtam_; }

  /// routine to return fine-scale scalar field fsphi at time step n+alpha_M
  virtual Teuchos::RCP<Epetra_Vector> FsPhi() { return fsphiaf_; };

  /// routine to return thermo. press. at time step n+alpha_F for low-Mach-number flow
  double ThermPressAf() { return thermpressaf_; }

  /// routine to return thermo. press. at time step n+alpha_M for low-Mach-number flow
  double ThermPressAm() { return thermpressam_; }

  /// routine to return time derivative of thermo. press. at time step n+alpha_F for low-Mach-number flow
  double ThermPressDtAf() { return thermpressdtaf_; }

  /// routine to return time derivative of thermo. press. at time step n+alpha_M for low-Mach-number flow
  double ThermPressDtAm() { return thermpressdtam_; }

  /// redistribute the scatra discretization and vectors according to nodegraph
  void Redistribute(const Teuchos::RCP<Epetra_CrsGraph> nodegraph);

  /// routine to return time integration specific parameters
  Teuchos::RCP<Teuchos::ParameterList> ScatraTimeParameterList()
  {
    Teuchos::RCP<Teuchos::ParameterList> timeparams;
    timeparams = Teuchos::rcp(new Teuchos::ParameterList());
    timeparams->set("using stationary formulation",false);
    timeparams->set("using generalized-alpha time integration",true);
    timeparams->set("total time",time_-(1-alphaF_)*dta_);
    timeparams->set("time factor",genalphafac_*dta_);
    timeparams->set("alpha_F",alphaF_);
    return timeparams;
  }


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
  double ResidualScaling() const { return 1.0/(dta_*genalphafac_); }

  /// update time-dependent electrode state variables at the end of an time step
  void ElectrodeKineticsTimeUpdate(const bool init = false);

private:

  /// scalar at time n+alpha_F and n+alpha_M
  Teuchos::RCP<Epetra_Vector>  phiaf_;
  Teuchos::RCP<Epetra_Vector>  phiam_;

  /// scalar time derivative at time n+alpha_M
  Teuchos::RCP<Epetra_Vector>  phidtam_;

  /// fine-scale part at time n+alpha_F
  Teuchos::RCP<Epetra_Vector>  fsphiaf_;

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
