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

  /// Update the solution after convergence of the nonlinear iteration.
  /// Current solution becomes old solution of next timestep.
  void Update();

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

  /// write additional data required for restart
  void OutputRestart();

private:

  /// solution vector phi at time t_{n-1}
  RCP<Epetra_Vector>    phinm_;

}; // class TimIntBDF2

} // namespace SCATRA


#endif  // #ifndef SCATRA_TIMINT_BDF2_H
#endif  // #ifdef CCADISCRET
