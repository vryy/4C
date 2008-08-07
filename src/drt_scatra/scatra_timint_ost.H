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


  /// Set the part of the righthandside belonging to the last timestep.
  void SetOldPartOfRighthandside();

  /// Update the solution after convergence of the nonlinear iteration.
  /// Current solution becomes old solution of next timestep.
  void Update();

  /// write additional data required for restart
  void OutputRestart();

  /// read restart data
  void ReadRestart(int step);

protected:

  // don't want = operator and cctor
  TimIntOneStepTheta operator = (const TimIntOneStepTheta& old);
  //! copy constructor
  TimIntOneStepTheta (const TimIntOneStepTheta& old);

private:

  /// time derivative of slution vector phi at time t_{n}
  RCP<Epetra_Vector>    phidtn_;

}; // class TimIntOneStepTheta

} // namespace SCATRA

#endif  // #ifndef SCATRA_TIMINT_OST_H
#endif  // #ifdef CCADISCRET
