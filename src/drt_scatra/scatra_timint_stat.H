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
                  RCP<IO::DiscretizationWriter>  output);

  /// Destructor
  virtual ~TimIntStationary();

  //! set velocity field for low-Mach-number flow
  void SetLomaVelocity(RCP<const Epetra_Vector> extvel,
                       RCP<DRT::Discretization> fluiddis);

  //! predict density field for low-Mach-number flow
  void PredictDensity(){return;};

  //! compute time derivative of density field for low-Mach-number flow
  void ComputeDensityDerivative(){return;};

  /// Update the solution after convergence of the nonlinear iteration.
  /// Current solution becomes old solution of next timestep.
  void Update(){return;};

  //! update density field for low-Mach-number flow
  void UpdateDensity(){return;};

  /// read restart data
  void ReadRestart(int step);


protected:

  // don't want = operator and cctor
  TimIntStationary operator = (const TimIntStationary& old);

  //! copy constructor
  TimIntStationary (const TimIntStationary& old);

  /// Initialization procedure before the first time step is done
  void PrepareFirstTimeStep(){return;};

  /// Set the part of the righthandside belonging to the last timestep.
  void SetOldPartOfRighthandside();

  /// do explicit predictor step (nothing to predict for stationary problems!)
  void ExplicitPredictor(){return;};

  /// set time for evaluation of Neumann boundary conditions
  void SetTimeForNeumannEvaluation(ParameterList& params);

  /// reset the residual vector and add actual Neumann loads with time factor
  void AddNeumannToResidual();

  /// add parameters specific for time-integration scheme
  void AddSpecificTimeIntegrationParameters(ParameterList& params);

  /// write additional data required for restart
  void OutputRestart(){return;};


private:


}; // class TimIntStationary

} // namespace SCATRA

#endif  // #ifndef SCATRA_TIMINT_STAT_H
#endif  // #ifdef CCADISCRET
