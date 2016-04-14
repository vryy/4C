/*!----------------------------------------------------------------------
\file tutorial_fixedpoint_scheme.H

\brief student's c++/baci simple fixed point problem

\maintainer  Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240

*----------------------------------------------------------------------*/
#ifndef TUTORIAL_FIXEDPOINT_SCHEME_H
#define TUTORIAL_FIXEDPOINT_SCHEME_H


namespace TUTORIAL {

class FixedPointScheme
{

public:

  /// constructor (ctor)
  explicit FixedPointScheme();

  /// execute the time loop
  void TimeLoop();

  /// return the relaxation parameter
  double Omega(){return omega_;};

  /// return the current simulation time
  double Time(){return time_;};

  /// return the current simulation time
  void IncrementTime(double dt){time_=time_+dt;};

  /// return the time step size
  double Dt(){return dt_;};

  /// return max time
  double MaxTime(){return Tend_;};

private:

  /// print the chosen tutorial type to screen
  void PrintTutorialType();

  /// define parameters for the toy problem
  void ProblemDefinition();

  /// call the partitioned fixed point iteration to equilibrate the time step
  void IterateFixedPoint();

  /// Solve Field 1
  double Operator1(double Fadh);

  /// Solve Field 2
  double Operator2(double disp);

  /// check if fixed-point scheme is converged
  bool ConvergenceCheck(double inc, double tol);

  /// predict the new state
  double InitialGuess();

  /// relaxation parameter
  double omega_;

  /// current simulation time
  double time_;

  /// timestep
  double dt_;

  /// max time
  double Tend_;

  /// convergence tolerance
  double convtol_;

  /// active force of body 2
  double Fact_;

  /// maximum active force
  double Fmax_;

  /// time when maximum active force is reached
  double Tm_;

  /// stiffnesses
  double K_cell_;
  double E_1_;
  double E_2_;

  // displacement state
  double disp_;

  /// current state of partitioned scheme
  double x_;

};// class FixedPointScheme
} // namespace TUTORIAL

#endif
