/*----------------------------------------------------------------------*/
/*! \file

\brief Generalized-alpha time-integration scheme with extensions for
       loma problems

\level 2


*----------------------------------------------------------------------*/

#ifndef FOUR_C_SCATRA_TIMINT_LOMA_GENALPHA_HPP
#define FOUR_C_SCATRA_TIMINT_LOMA_GENALPHA_HPP

#include "baci_config.hpp"

#include "baci_linalg_utils_sparse_algebra_math.hpp"
#include "baci_scatra_timint_genalpha.hpp"
#include "baci_scatra_timint_loma.hpp"

BACI_NAMESPACE_OPEN


namespace SCATRA
{
  class TimIntLomaGenAlpha : public virtual ScaTraTimIntLoma, public virtual TimIntGenAlpha
  {
   public:
    /// Standard Constructor
    TimIntLomaGenAlpha(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);


    /// initialize time integration scheme
    void Init() override;

    /// setup time integration scheme
    void Setup() override;

    /// predict thermodynamic pressure and time derivative for low-Mach-number flow
    void PredictThermPressure() override;

    /// compute values of thermodynamic pressure at intermediate time steps
    void ComputeThermPressureIntermediateValues() override;

    /// compute thermodynamic pressure and time derivative for low-Mach-number flow
    void ComputeThermPressure() override;

    ///  compute time derivative of thermodynamic pressure
    void ComputeThermPressureTimeDerivative() override;

    /// update thermodynamic pressure and time derivative for low-Mach-number flow
    void UpdateThermPressure() override;

    /// read restart data
    void ReadRestart(const int step, Teuchos::RCP<IO::InputControl> input = Teuchos::null) override;

    /// routine to return thermo. press. at time step n+alpha_F for low-Mach-number flow
    double ThermPressAf() override { return thermpressaf_; }

    /// routine to return thermo. press. at time step n+alpha_M for low-Mach-number flow
    double ThermPressAm() override { return thermpressam_; }

    /// routine to return time derivative of thermo. press. at time step n+alpha_F for
    /// low-Mach-number flow
    double ThermPressDtAf() override { return thermpressdtaf_; }

    /// routine to return time derivative of thermo. press. at time step n+alpha_M for
    /// low-Mach-number flow
    double ThermPressDtAm() override { return thermpressdtam_; }

   protected:
    void WriteRestart() const override;

    /// dynamic Smagorinsky model
    void DynamicComputationOfCs() override;

    /// dynamic Vreman model
    void DynamicComputationOfCv() override;

    /// add thermodynamic pressure to parameter list for element evaluation
    void AddThermPressToParameterList(Teuchos::ParameterList& params  //!< parameter list
        ) override;

   private:
    /// LOMA-specific parameter: thermodynamic pressure at n+alpha_F and n+alpha_M
    /// and its time derivative at n+alpha_F and n+alpha_M
    double thermpressaf_;
    double thermpressam_;
    double thermpressdtaf_;
    double thermpressdtam_;


  };  // class TimIntLomaGenAlpha

}  // namespace SCATRA

BACI_NAMESPACE_CLOSE

#endif  // BACI_SCATRA_TIMINT_LOMA_GENALPHA_H
