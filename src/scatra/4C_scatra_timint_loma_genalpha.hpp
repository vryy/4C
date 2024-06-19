/*----------------------------------------------------------------------*/
/*! \file

\brief Generalized-alpha time-integration scheme with extensions for
       loma problems

\level 2


*----------------------------------------------------------------------*/

#ifndef FOUR_C_SCATRA_TIMINT_LOMA_GENALPHA_HPP
#define FOUR_C_SCATRA_TIMINT_LOMA_GENALPHA_HPP

#include "4C_config.hpp"

#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_scatra_timint_genalpha.hpp"
#include "4C_scatra_timint_loma.hpp"

FOUR_C_NAMESPACE_OPEN


namespace ScaTra
{
  class TimIntLomaGenAlpha : public virtual ScaTraTimIntLoma, public virtual TimIntGenAlpha
  {
   public:
    /// Standard Constructor
    TimIntLomaGenAlpha(Teuchos::RCP<Core::FE::Discretization> dis,
        Teuchos::RCP<Core::LinAlg::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<Core::IO::DiscretizationWriter> output);


    /// initialize time integration scheme
    void init() override;

    /// setup time integration scheme
    void setup() override;

    /// predict thermodynamic pressure and time derivative for low-Mach-number flow
    void predict_therm_pressure() override;

    /// compute values of thermodynamic pressure at intermediate time steps
    void compute_therm_pressure_intermediate_values() override;

    /// compute thermodynamic pressure and time derivative for low-Mach-number flow
    void compute_therm_pressure() override;

    ///  compute time derivative of thermodynamic pressure
    void compute_therm_pressure_time_derivative() override;

    /// update thermodynamic pressure and time derivative for low-Mach-number flow
    void UpdateThermPressure() override;

    /// read restart data
    void read_restart(
        const int step, Teuchos::RCP<Core::IO::InputControl> input = Teuchos::null) override;

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
    void write_restart() const override;

    /// dynamic Smagorinsky model
    void dynamic_computation_of_cs() override;

    /// dynamic Vreman model
    void dynamic_computation_of_cv() override;

    /// add thermodynamic pressure to parameter list for element evaluation
    void add_therm_press_to_parameter_list(Teuchos::ParameterList& params  //!< parameter list
        ) override;

   private:
    /// LOMA-specific parameter: thermodynamic pressure at n+alpha_F and n+alpha_M
    /// and its time derivative at n+alpha_F and n+alpha_M
    double thermpressaf_;
    double thermpressam_;
    double thermpressdtaf_;
    double thermpressdtam_;


  };  // class TimIntLomaGenAlpha

}  // namespace ScaTra

FOUR_C_NAMESPACE_CLOSE

#endif
