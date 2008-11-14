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

  /// Update the solution after convergence of the nonlinear iteration.
  /// Current solution becomes old solution of next timestep.
  void Update();

  /// read restart data
  void ReadRestart(int step);


protected:

  // don't want = operator and cctor
  TimIntGenAlpha operator = (const TimIntGenAlpha& old);

  //! copy constructor
  TimIntGenAlpha (const TimIntGenAlpha& old);

  /// Initialization procedure before the first time step is done
  void PrepareFirstTimeStep();

  /// Set the part of the righthandside belonging to the last timestep.
  void SetOldPartOfRighthandside();

  /// do explicit predictor step (-> better starting value for nonlinear solver)
  void ExplicitPredictor();

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

  /// time factors for generalized-alpha time integration
  double alphaM_;
  double alphaF_;
  double gamma_;
  double genalphafac_;

}; // class TimIntGenAlpha

} // namespace SCATRA

#endif  // #ifndef SCATRA_TIMINT_GENALPHA_H
#endif  // #ifdef CCADISCRET
