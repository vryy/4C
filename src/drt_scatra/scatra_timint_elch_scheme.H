/*!----------------------------------------------------------------------
\file scatra_timint_elch_scheme.H

\brief  connecting time-integration schemes (OST, BDF2, GenAlpha, Stationary) with
        elch-specific implementation (class ScaTraTimIntElch)

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>

*----------------------------------------------------------------------*/

#ifndef SCATRA_TIMINT_ELCH_SCHEME_H
#define SCATRA_TIMINT_ELCH_SCHEME_H

#include "scatra_timint_ost.H"
#include "scatra_timint_bdf2.H"
#include "scatra_timint_genalpha.H"
#include "scatra_timint_stat.H"
#include "scatra_timint_elch.H"

namespace SCATRA
{

class ScaTraTimIntElchOST : public ScaTraTimIntElch, public TimIntOneStepTheta
{

public:

  //! Standard Constructor
  ScaTraTimIntElchOST(
    Teuchos::RCP<DRT::Discretization>       dis,
    Teuchos::RCP<LINALG::Solver>            solver,
    Teuchos::RCP<Teuchos::ParameterList>    params,
    Teuchos::RCP<Teuchos::ParameterList>    sctratimintparams,
    Teuchos::RCP<Teuchos::ParameterList>    extraparams,
    Teuchos::RCP<IO::DiscretizationWriter>  output);

  //! Destructor
  virtual ~ScaTraTimIntElchOST();

  //! initialize time integration scheme
  void Init();

  //! Update the solution after convergence of the nonlinear iteration.
  //! Current solution becomes old solution of next timestep.
  void Update(const int num=0);

  //! read restart data
  void ReadRestart(int step);

protected:

  //! write additional data required for restart
  void OutputRestart();

  //! update time-dependent electrode state variables at the end of an time step
  void ElectrodeKineticsTimeUpdate();

  //!  compute scalar time derivate parameters of the input voltage to compute double layer current densities
  void ComputeTimeDerivPot0(const bool init);

  //! set old part of RHS for galvanostatic equation
  void ElectrodeKineticsSetOldPartOfRHS();
}; // class TimIntElchOST


class ScaTraTimIntElchBDF2 : public ScaTraTimIntElch, public TimIntBDF2
{

public:

  //! Standard Constructor
  ScaTraTimIntElchBDF2(
    Teuchos::RCP<DRT::Discretization>       dis,
    Teuchos::RCP<LINALG::Solver>            solver,
    Teuchos::RCP<Teuchos::ParameterList>    params,
    Teuchos::RCP<Teuchos::ParameterList>    sctratimintparams,
    Teuchos::RCP<Teuchos::ParameterList>    extraparams,
    Teuchos::RCP<IO::DiscretizationWriter>  output);

  //! Destructor
  virtual ~ScaTraTimIntElchBDF2();

  //! initialize time integration scheme
  void Init();

  //! Update the solution after convergence of the nonlinear iteration.
  //! Current solution becomes old solution of next timestep.
  void Update(const int num=0);

  //! read restart data
  void ReadRestart(int step);

protected:

  //! write additional data required for restart
  void OutputRestart();

  //! update time-dependent electrode state variables at the end of an time step
  void ElectrodeKineticsTimeUpdate();

  //!  compute scalar time derivate parameters of the input voltage to compute double layer current densities
  void ComputeTimeDerivPot0(const bool init);

  //! set old part of RHS for galvanostatic equation
  void ElectrodeKineticsSetOldPartOfRHS();
}; // class TimIntElchBDF2


class ScaTraTimIntElchGenAlpha : public ScaTraTimIntElch, public TimIntGenAlpha
{

public:

  //! Standard Constructor
  ScaTraTimIntElchGenAlpha(
    Teuchos::RCP<DRT::Discretization>       dis,
    Teuchos::RCP<LINALG::Solver>            solver,
    Teuchos::RCP<Teuchos::ParameterList>    params,
    Teuchos::RCP<Teuchos::ParameterList>    sctratimintparams,
    Teuchos::RCP<Teuchos::ParameterList>    extraparams,
    Teuchos::RCP<IO::DiscretizationWriter>  output);

  //! Destructor
  virtual ~ScaTraTimIntElchGenAlpha();

  //! initialize time integration scheme
  void Init();

  //! Update the solution after convergence of the nonlinear iteration.
  //! Current solution becomes old solution of next timestep.
  void Update(const int num=0);

  //! read restart data
  void ReadRestart(int step);

protected:

  //! write additional data required for restart
  void OutputRestart();

  //! update time-dependent electrode state variables at the end of an time step
  void ElectrodeKineticsTimeUpdate();

  //!  compute scalar time derivate parameters of the input voltage to compute double layer current densities
  void ComputeTimeDerivPot0(const bool init);

  //! set old part of RHS for galvanostatic equation
  void ElectrodeKineticsSetOldPartOfRHS(){dserror(""); return;};
}; // class TimIntElchGenAlpha


class ScaTraTimIntElchStationary : public ScaTraTimIntElch, public TimIntStationary
{

public:

  //! Standard Constructor
  ScaTraTimIntElchStationary(
    Teuchos::RCP<DRT::Discretization>       dis,
    Teuchos::RCP<LINALG::Solver>            solver,
    Teuchos::RCP<Teuchos::ParameterList>    params,
    Teuchos::RCP<Teuchos::ParameterList>    sctratimintparams,
    Teuchos::RCP<Teuchos::ParameterList>    extraparams,
    Teuchos::RCP<IO::DiscretizationWriter>  output);

  //! Destructor
  virtual ~ScaTraTimIntElchStationary();

  //! initialize time integration scheme
  void Init();

  //! Update the solution after convergence of the nonlinear iteration.
  //! Current solution becomes old solution of next timestep.
  void Update(const int num=0);

  //! read restart data
  void ReadRestart(int step);

protected:

  //! write additional data required for restart
  void OutputRestart();

  //! update time-dependent electrode state variables at the end of an time step
  void ElectrodeKineticsTimeUpdate(){dserror("Galvanostatic-BC is not implemented for the stationary time-integration scheme"); return;};

  //!  compute scalar time derivate parameters of the input voltage to compute double layer current densities
  void ComputeTimeDerivPot0(const bool init);

  //! set old part of RHS for galvanostatic equation
  void ElectrodeKineticsSetOldPartOfRHS(){return;};
}; // class TimIntElchStatationary
} // namespace SCATRA

#endif  // #ifndef SCATRA_TIMINT_ELCH_SCHEME_H
