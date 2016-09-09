/*!----------------------------------------------------------------------
\file scatra_timint_loma_ost.H

\brief One-step-theta time-integration scheme with extensions for
       loma problems

\level 2

\maintainer Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245

*----------------------------------------------------------------------*/

#ifndef SCATRA_TIMINT_LOMA_OST_H
#define SCATRA_TIMINT_LOMA_OST_H

#include "scatra_timint_loma.H"
#include "scatra_timint_ost.H"

#include "../linalg/linalg_utils.H"


namespace SCATRA
{

class TimIntLomaOST : public virtual ScaTraTimIntLoma, public virtual TimIntOneStepTheta
{

public:

  /// Standard Constructor
  TimIntLomaOST(
    Teuchos::RCP<DRT::Discretization>       dis,
    Teuchos::RCP<LINALG::Solver>            solver,
    Teuchos::RCP<Teuchos::ParameterList>    params,
    Teuchos::RCP<Teuchos::ParameterList>    sctratimintparams,
    Teuchos::RCP<Teuchos::ParameterList>    extraparams,
    Teuchos::RCP<IO::DiscretizationWriter>  output);

  /// Destructor
  virtual ~TimIntLomaOST();

  /// initialize time integration scheme
  virtual void Init();

  /// setup time integration scheme
  virtual void Setup();

  /// predict thermodynamic pressure and time derivative for low-Mach-number flow
  void PredictThermPressure();

  /// compute values of thermodynamic pressure at intermediate time steps
  void ComputeThermPressureIntermediateValues(){return;};

  /// compute thermodynamic pressure and time derivative for low-Mach-number flow
  void ComputeThermPressure();

  ///  compute time derivative of thermodynamic pressure
  void ComputeThermPressureTimeDerivative();

  /// update thermodynamic pressure and time derivative for low-Mach-number flow
  void UpdateThermPressure();

  /// read restart data
  void ReadRestart(const int step,Teuchos::RCP<IO::InputControl> input = Teuchos::null);

  /// routine to return thermo. press. at time step n+alpha_F for low-Mach-number flow
  double ThermPressAf() { return 1.0; }

  /// routine to return thermo. press. at time step n+alpha_M for low-Mach-number flow
  double ThermPressAm() { return 1.0; }

  /// routine to return time derivative of thermo. press. at time step n+alpha_F for low-Mach-number flow
  double ThermPressDtAf() { return 0.0; }

  /// routine to return time derivative of thermo. press. at time step n+alpha_M for low-Mach-number flow
  double ThermPressDtAm() { return 0.0; }

protected:

  /// write additional data required for restart
  void OutputRestart();

  /// dynamic Smagorinsky model
  void DynamicComputationOfCs();

  /// dynamic Vreman model
  void DynamicComputationOfCv();

  /// add thermodynamic pressure to parameter list for element evaluation
  void AddThermPressToParameterList(
    Teuchos::ParameterList& params //!< parameter list
  );

}; // class TimIntLomaOST

} // namespace SCATRA

#endif  // #ifndef SCATRA_TIMINT_LOMA_OST_H
