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

  /// set velocity field for low-Mach-number flow
  void SetLomaVelocity(RCP<const Epetra_Vector> extvel,
                       RCP<DRT::Discretization> fluiddis);

  /// compute initial time derivative of density field for low-Mach-number flow
  void ComputeInitialDensityDerivative();

  /// predict thermodynamic pressure and time derivative for low-Mach-number flow
  void PredictThermPressure();

  /// predict density field and time derivative for low-Mach-number flow
  void PredictDensity();

  /// compute values at intermediate time steps (required for generalized-alpha)
  void ComputeIntermediateValues(){return;};

  /// compute thermodynamic pressure and time derivative for low-Mach-number flow
  double ComputeThermPressure();

  /// compute time derivative of density field for low-Mach-number flow
  void ComputeDensityDerivative(){return;};

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

  /// reset the residual vector and add actual Neumann loads with time factor
  void AddNeumannToResidual();

  /// add parameters specific for time-integration scheme
  void AddSpecificTimeIntegrationParameters(ParameterList& params);

  /// write additional data required for restart
  void OutputRestart();

  /// time derivative of solution vector phi at time t_{n}
  RCP<Epetra_Vector>    phidtn_;


private:

  /// calculate initial value for phidtn_
  void CalcInitialPhidt();

  /// time factor for one-step-theta/BDF2 time integration
  double theta_;

}; // class TimIntOneStepTheta

} // namespace SCATRA

#endif  // #ifndef SCATRA_TIMINT_OST_H
#endif  // #ifdef CCADISCRET
