/*----------------------------------------------------------------------*/
/*! \file
 \brief Managing parameters of (reaction) functions for poromultiphase_scatra
        problems

\level 3

*----------------------------------------------------------------------*/

#ifndef BACI_POROMULTIPHASE_SCATRA_FUNCTION_PARAMETERS_HPP
#define BACI_POROMULTIPHASE_SCATRA_FUNCTION_PARAMETERS_HPP


#include "baci_config.hpp"

#include <array>
#include <string>
#include <vector>

BACI_NAMESPACE_OPEN


namespace POROMULTIPHASESCATRA
{

  /*!
   *  @brief struct for managing the parameters of the POROMULTIPHASESCATRA_FUNCTION
   *  TUMOR_GROWTH_LAW_HEAVISIDE
   */
  struct TumorGrowthLawHeavisideParameters
  {
    TumorGrowthLawHeavisideParameters() = default;

    TumorGrowthLawHeavisideParameters(
        const std::vector<std::pair<std::string, double>>& funct_params);

    double gamma_T_growth{};
    double w_nl_crit{};
    double w_nl_env{};
    double lambda{};
    double p_t_crit{};
  };

  /*!
   *  @brief struct for managing the parameters of the POROMULTIPHASESCATRA_FUNCTION
   *  NECROSIS_LAW_HEAVISIDE
   */
  struct NecrosisLawHeavisideParameters
  {
    NecrosisLawHeavisideParameters() = default;

    NecrosisLawHeavisideParameters(const std::vector<std::pair<std::string, double>>& funct_params);

    double gamma_t_necr{};
    double w_nl_crit{};
    double w_nl_env{};
    double delta_a_t{};
    double p_t_crit{};
  };

  /*!
   *  @brief struct for managing the parameters of the POROMULTIPHASESCATRA_FUNCTION
   *  OXYGEN_CONSUMPTION_LAW_HEAVISIDE
   */
  struct OxygenConsumptionLawHeavisideParameters
  {
    OxygenConsumptionLawHeavisideParameters() = default;

    OxygenConsumptionLawHeavisideParameters(
        const std::vector<std::pair<std::string, double>>& funct_params);

    double gamma_nl_growth{};
    double gamma_0_nl{};
    double w_nl_crit{};
    double w_nl_env{};
    double p_t_crit{};
  };


  /*!
   *  @brief struct for managing the parameters of the POROMULTIPHASESCATRA_FUNCTION
   *  TUMOR_GROWTH_LAW_HEAVISIDE_OXY and TUMOR_GROWTH_LAW_HEAVISIDE_NECRO
   */
  struct TumorGrowthLawHeavisideNecroOxyParameters
  {
    TumorGrowthLawHeavisideNecroOxyParameters() = default;

    TumorGrowthLawHeavisideNecroOxyParameters(
        const std::vector<std::pair<std::string, double>>& funct_params);

    double gamma_T_growth{};
    double w_nl_crit{};
    double w_nl_env{};
    double lambda{};
    double p_t_crit{};
  };

  /*!
   *  @brief struct for managing the parameters of the POROMULTIPHASESCATRA_FUNCTION
   *  OXYGEN_TRANSVASCULAR_EXCHANGE_LAW_CONT
   */
  struct OxygenTransvascularExchangeLawContParameters
  {
    OxygenTransvascularExchangeLawContParameters() = default;

    OxygenTransvascularExchangeLawContParameters(
        const std::vector<std::pair<std::string, double>>& funct_params);

    double n{};
    double Pb50{};
    double CaO2_max{};
    double alpha_bl_eff{};
    double gammarhoSV{};
    double rho_oxy{};
    double rho_if{};
    double rho_bl{};
    double alpha_IF{};
  };

  /*!
   *  @brief struct for managing the parameters of the POROMULTIPHASESCATRA_FUNCTION
   *  OXYGEN_TRANSVASCULAR_EXCHANGE_LAW_DISC
   */
  struct OxygenTransvascularExchangeLawDiscParameters
  {
    OxygenTransvascularExchangeLawDiscParameters() = default;

    OxygenTransvascularExchangeLawDiscParameters(
        const std::vector<std::pair<std::string, double>>& funct_params);

    double n{};
    double Pb50{};
    double CaO2_max{};
    double alpha_bl_eff{};
    double gammarho{};
    double rho_oxy{};
    double rho_if{};
    double rho_bl{};
    double S2_max{};
    double alpha_IF{};
  };

  /*!
   *  @brief struct for managing the parameters of the POROMULTIPHASESCATRA_FUNCTION
   *  LUNG_OXYGEN_EXCHANGE_LAW
   */
  struct LungOxygenExchangeLawParameters
  {
    LungOxygenExchangeLawParameters() = default;

    LungOxygenExchangeLawParameters(
        const std::vector<std::pair<std::string, double>>& funct_params);

    double rho_oxy{};
    double DiffAdVTLC{};
    double alpha_oxy{};
    double rho_air{};
    double rho_bl{};
    double n{};
    double P_oB50{};
    double NC_Hb{};
    double P_atmospheric{};
    double volfrac_blood_ref{};
  };

  /*!
   *  @brief struct for managing the parameters of the POROMULTIPHASESCATRA_FUNCTION
   *  LUNG_CARBONDIOXIDE_EXCHANGE_LAW
   */
  struct LungCarbonDioxideExchangeLawParameters
  {
    LungCarbonDioxideExchangeLawParameters() = default;

    LungCarbonDioxideExchangeLawParameters(
        const std::vector<std::pair<std::string, double>>& funct_params);

    double rho_CO2{};
    double DiffsolAdVTLC{};
    double pH{};
    double rho_air{};
    double rho_bl{};
    double rho_oxy{};
    double n{};
    double P_oB50{};
    double C_Hb{};
    double NC_Hb{};
    double alpha_oxy{};
    double P_atmospheric{};
    double ScalingFormmHg{};
    double volfrac_blood_ref{};
  };
}  // namespace POROMULTIPHASESCATRA


BACI_NAMESPACE_CLOSE

#endif