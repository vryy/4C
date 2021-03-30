/*----------------------------------------------------------------------*/
/*! \file

\brief One-step-theta time-integration scheme with extensions for
       loma problems

\level 2


*----------------------------------------------------------------------*/

#ifndef SCATRA_TIMINT_LOMA_OST_H
#define SCATRA_TIMINT_LOMA_OST_H

#include "scatra_timint_loma.H"
#include "scatra_timint_ost.H"

#include "../linalg/linalg_utils_sparse_algebra_math.H"


namespace SCATRA
{
  class TimIntLomaOST : public virtual ScaTraTimIntLoma, public virtual TimIntOneStepTheta
  {
   public:
    /// Standard Constructor
    TimIntLomaOST(Teuchos::RCP<DRT::Discretization> dis, Teuchos::RCP<LINALG::Solver> solver,
        Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);

    /// initialize time integration scheme
    void Init() override;

    /// setup time integration scheme
    void Setup() override;

    /// predict thermodynamic pressure and time derivative for low-Mach-number flow
    void PredictThermPressure() override{};

    /// compute values of thermodynamic pressure at intermediate time steps
    void ComputeThermPressureIntermediateValues() override{};

    /// compute thermodynamic pressure and time derivative for low-Mach-number flow
    void ComputeThermPressure() override;

    ///  compute time derivative of thermodynamic pressure
    void ComputeThermPressureTimeDerivative() override;

    /// update thermodynamic pressure and time derivative for low-Mach-number flow
    void UpdateThermPressure() override;

    /// read restart data
    void ReadRestart(const int step, Teuchos::RCP<IO::InputControl> input = Teuchos::null) override;

    /// routine to return thermo. press. at time step n+alpha_F for low-Mach-number flow
    double ThermPressAf() override { return 1.0; }

    /// routine to return thermo. press. at time step n+alpha_M for low-Mach-number flow
    double ThermPressAm() override { return 1.0; }

    /// routine to return time derivative of thermo. press. at time step n+alpha_F for
    /// low-Mach-number flow
    double ThermPressDtAf() override { return 0.0; }

    /// routine to return time derivative of thermo. press. at time step n+alpha_M for
    /// low-Mach-number flow
    double ThermPressDtAm() override { return 0.0; }

   protected:
    /// write additional data required for restart
    void OutputRestart() const override;

    /// dynamic Smagorinsky model
    void DynamicComputationOfCs() override;

    /// dynamic Vreman model
    void DynamicComputationOfCv() override;

    /*!
     * @brief add thermodynamic pressure to parameter list for element evaluation
     *
     * @param params parameter list
     */
    void AddThermPressToParameterList(Teuchos::ParameterList& params) override;
  };  // class TimIntLomaOST
}  // namespace SCATRA

#endif  // #ifndef SCATRA_TIMINT_LOMA_OST_H
