/*!----------------------------------------------------------------------
\file scatra_timint_ost.H
\brief

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>

*----------------------------------------------------------------------*/

#ifdef CCADISCRET
#ifndef SCATRA_TIMINT_OST_H
#define SCATRA_TIMINT_OST_H

#include "scatra_timint_implicit.H"


namespace SCATRA
{

class TimIntOneStepTheta : public ScaTraTimIntImpl
{

public:

  /// Standard Constructor
  TimIntOneStepTheta(RCP<DRT::Discretization>       dis,
                     RCP<LINALG::Solver>            solver,
                     RCP<ParameterList>             params,
                     RCP<ParameterList>             extraparams,
                     RCP<IO::DiscretizationWriter>  output);

  /// Destructor
  virtual ~TimIntOneStepTheta();

  /// Print information about current time step to screen (reimplementation for OST)
  virtual void PrintTimeStepInfo()
  {
    if (myrank_==0)
        printf("TIME: %11.4E/%11.4E  DT = %11.4E  %s (theta = %3.2f)  STEP = %4d/%4d \n",
               time_,maxtime_,dta_,MethodTitle().c_str(),theta_,step_,stepmax_);
    return;
  }

  /// predict thermodynamic pressure and time derivative for low-Mach-number flow
  void PredictThermPressure();

  /// compute values at intermediate time steps (required for generalized-alpha)
  void ComputeIntermediateValues(){return;};

  /// compute values of thermodynamic pressure at intermediate time steps (required for generalized-alpha)
  void ComputeThermPressureIntermediateValues(){return;};

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
  void UpdateDensityElch();

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
  TimIntOneStepTheta operator = (const TimIntOneStepTheta& old);

  /// copy constructor
  TimIntOneStepTheta (const TimIntOneStepTheta& old);

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

  /// update time-dependent electrode state variables at the end of an time step
  void ElectrodeKineticsTimeUpdate(const bool init = false);

private:

  /// time factor for one-step-theta/BDF2 time integration
  double theta_;

  /// fine-scale solution vector at time n+1
  RCP<Epetra_Vector>  fsphinp_;

}; // class TimIntOneStepTheta

} // namespace SCATRA

#endif  // #ifndef SCATRA_TIMINT_OST_H
#endif  // #ifdef CCADISCRET
