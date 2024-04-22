/*----------------------------------------------------------------------*/
/*! \file
 \brief Managing and evaluating of (reaction) functions for poromultiphase_scatra
        problems

\level 3

    *----------------------------------------------------------------------*/

#include "baci_poromultiphase_scatra_function.hpp"

#include "baci_global_data.hpp"
#include "baci_io_linedefinition.hpp"
#include "baci_poromultiphase_scatra_utils.hpp"
#include "baci_utils_fad.hpp"
#include "baci_utils_function_manager.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace
{

  template <int dim>
  Teuchos::RCP<CORE::UTILS::FunctionOfAnything> CreatePoroFunction(
      const std::string& type, const std::vector<std::pair<std::string, double>>& params)
  {
    if (type == "TUMOR_GROWTH_LAW_HEAVISIDE")
      return Teuchos::rcp(new POROMULTIPHASESCATRA::TumorGrowthLawHeaviside<dim>(params));
    else if (type == "NECROSIS_LAW_HEAVISIDE")
      return Teuchos::rcp(new POROMULTIPHASESCATRA::NecrosisLawHeaviside<dim>(params));
    else if (type == "OXYGEN_CONSUMPTION_LAW_HEAVISIDE")
      return Teuchos::rcp(new POROMULTIPHASESCATRA::OxygenConsumptionLawHeaviside<dim>(params));
    else if (type == "TUMOR_GROWTH_LAW_HEAVISIDE_OXY")
      return Teuchos::rcp(new POROMULTIPHASESCATRA::TumorGrowthLawHeavisideOxy<dim>(params));
    else if (type == "TUMOR_GROWTH_LAW_HEAVISIDE_NECRO")
      return Teuchos::rcp(new POROMULTIPHASESCATRA::TumorGrowthLawHeavisideNecro<dim>(params));
    else if (type == "OXYGEN_TRANSVASCULAR_EXCHANGE_LAW_CONT")
    {
      return Teuchos::rcp(
          new POROMULTIPHASESCATRA::OxygenTransvascularExchangeLawCont<dim>(params));
    }
    else if (type == "OXYGEN_TRANSVASCULAR_EXCHANGE_LAW_DISC")
    {
      return Teuchos::rcp(
          new POROMULTIPHASESCATRA::OxygenTransvascularExchangeLawDisc<dim>(params));
    }
    else if (type == "LUNG_OXYGEN_EXCHANGE_LAW")
    {
      return Teuchos::rcp(new POROMULTIPHASESCATRA::LungOxygenExchangeLaw<dim>(params));
    }
    else if (type == "LUNG_CARBONDIOXIDE_EXCHANGE_LAW")
    {
      return Teuchos::rcp(new POROMULTIPHASESCATRA::LungCarbonDioxideExchangeLaw<dim>(params));
    }
    else
    {
      FOUR_C_THROW("Wrong type of POROMULTIPHASESCATRA_FUNCTION");
      return Teuchos::RCP<CORE::UTILS::FunctionOfAnything>(nullptr);
    }
  }



  template <int dim>
  Teuchos::RCP<CORE::UTILS::FunctionOfAnything> TryCreatePoroFunction(
      const std::vector<INPUT::LineDefinition>& function_line_defs)
  {
    if (function_line_defs.size() != 1) return Teuchos::null;

    const auto& function_lin_def = function_line_defs.front();

    if (function_lin_def.HaveNamed("POROMULTIPHASESCATRA_FUNCTION"))
    {
      std::string type;
      function_lin_def.ExtractString("POROMULTIPHASESCATRA_FUNCTION", type);

      std::vector<std::pair<std::string, double>> params;
      if (function_lin_def.HaveNamed("PARAMS"))
        function_lin_def.ExtractPairOfStringAndDoubleVector("PARAMS", params);

      return CreatePoroFunction<dim>(type, params);
    }
    else
    {
      return Teuchos::RCP<CORE::UTILS::FunctionOfAnything>(nullptr);
    }
  }

  auto TryCreatePoroFunctionDispatch(const std::vector<INPUT::LineDefinition>& function_line_defs)
  {
    switch (GLOBAL::Problem::Instance()->NDim())
    {
      case 1:
        return TryCreatePoroFunction<1>(function_line_defs);
      case 2:
        return TryCreatePoroFunction<2>(function_line_defs);
      case 3:
        return TryCreatePoroFunction<3>(function_line_defs);
      default:
        FOUR_C_THROW("Unsupported dimension %d.", GLOBAL::Problem::Instance()->NDim());
    }
  }
}  // namespace

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
POROMULTIPHASESCATRA::PoroMultiPhaseScaTraFunction<dim>::PoroMultiPhaseScaTraFunction()
    : order_checked_(false)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::AddValidPoroFunctions(CORE::UTILS::FunctionManager& function_manager)
{
  function_manager.AddFunctionDefinition(
      {INPUT::LineDefinition::Builder()
              .AddNamedString("POROMULTIPHASESCATRA_FUNCTION")
              .AddOptionalNamedInt("NUMPARAMS")
              .AddOptionalNamedPairOfStringAndDoubleVector(
                  "PARAMS", INPUT::LengthFromIntNamed("NUMPARAMS"))
              .Build()},
      TryCreatePoroFunctionDispatch);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
POROMULTIPHASESCATRA::TumorGrowthLawHeaviside<dim>::TumorGrowthLawHeaviside(
    const std::vector<std::pair<std::string, double>>& funct_params)
    : PoroMultiPhaseScaTraFunction<dim>(), parameter(funct_params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
void POROMULTIPHASESCATRA::TumorGrowthLawHeaviside<dim>::CheckOrder(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants) const
{
  // safety check for correct ordering of variables and constants
  // they should have been added in exactly the same way in fluidporo_multiphase_singlereaction, but
  // order might be different if we do not use exactly three fluid phases
  if (variables[1].first != "p2")
    FOUR_C_THROW("wrong order in variable vector, p2 not at position 2");
  if (variables[4].first != "S2")
    FOUR_C_THROW("wrong order in variable vector, S2 not at position 5");
  if (variables[6].first != "porosity")
    FOUR_C_THROW("wrong order in variable vector, porosity not at position 7");
  if (variables[7].first != "phi1")
    FOUR_C_THROW("wrong order in variable vector, phi1 (oxygen mass fraction) not at position 8");
  if (variables[8].first != "phi2")
    FOUR_C_THROW("wrong order in variable vector, phi2 (necrotic mass fraction) not at position 9");

  // order is correct
  this->order_checked_ = true;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
double POROMULTIPHASESCATRA::TumorGrowthLawHeaviside<dim>::Evaluate(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants, const size_t component) const
{
  // Check order (only once since it does not change)
  if (not this->order_checked_) CheckOrder(variables, constants);

  // read variables and constants (order is crucial)
  const double p2 = variables[1].second;
  const double S2 = variables[4].second;
  const double porosity = variables[6].second;
  const double oxy_mass_frac = variables[7].second;
  const double necr_frac = variables[8].second;

  // evaluate heaviside
  const double heaviside_oxy((oxy_mass_frac - parameter.w_nl_crit) > 0. ? 1. : 0.);
  const double heaviside_pres((parameter.p_t_crit - p2) > 0. ? 1. : 0.);
  const double macaulay = parameter.gamma_T_growth * (oxy_mass_frac - parameter.w_nl_crit) /
                          (parameter.w_nl_env - parameter.w_nl_crit) * heaviside_oxy;

  // evaluate function
  const double functval = (macaulay * heaviside_pres) * (1 - necr_frac) * porosity * S2 -
                          parameter.lambda * porosity * oxy_mass_frac * S2;
  return functval;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
std::vector<double> POROMULTIPHASESCATRA::TumorGrowthLawHeaviside<dim>::EvaluateDerivative(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants, const size_t component) const
{
  // create derivative vector (should have size of variables)
  std::vector<double> deriv(variables.size(), 0.0);

  // read variables and constants
  const double p2 = variables[1].second;
  const double S2 = variables[4].second;
  const double porosity = variables[6].second;
  const double oxy_mass_frac = variables[7].second;
  const double necr_frac = variables[8].second;

  // evaluate heaviside
  const double heaviside_oxy((oxy_mass_frac - parameter.w_nl_crit) > 0. ? 1. : 0.);
  const double heaviside_pres((parameter.p_t_crit - p2) > 0. ? 1. : 0.);
  const double macaulay = parameter.gamma_T_growth * (oxy_mass_frac - parameter.w_nl_crit) /
                          (parameter.w_nl_env - parameter.w_nl_crit) * heaviside_oxy;

  // set saturation derivs w.r.t. S2
  const double saturationderiv = (macaulay * heaviside_pres) * (1 - necr_frac) * porosity -
                                 parameter.lambda * necr_frac * porosity;
  deriv[4] = saturationderiv;

  // set porosity derivs
  const double porosityderiv =
      (macaulay * heaviside_pres) * (1 - necr_frac) * S2 - parameter.lambda * necr_frac * S2;
  deriv[6] = porosityderiv;

  // set scalar derivs w.r.t. phi1
  const double oxygenderiv = parameter.gamma_T_growth * heaviside_oxy * heaviside_pres * 1.0 /
                             (parameter.w_nl_env - parameter.w_nl_crit) * (1 - necr_frac) *
                             porosity * S2;
  deriv[7] = oxygenderiv;

  // set scalar derivs w.r.t. phi2
  const double necroticderiv =
      -(macaulay * heaviside_pres) * porosity * S2 - parameter.lambda * porosity * S2;
  deriv[8] = necroticderiv;

  return deriv;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
POROMULTIPHASESCATRA::NecrosisLawHeaviside<dim>::NecrosisLawHeaviside(
    const std::vector<std::pair<std::string, double>>& funct_params)
    : PoroMultiPhaseScaTraFunction<dim>(), parameter(funct_params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
void POROMULTIPHASESCATRA::NecrosisLawHeaviside<dim>::CheckOrder(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants) const
{
  // safety check for correct ordering of variables and constants
  // they should have been added in exactly the same way in scatra_ele_calc_multiporo_reaction, but
  // order might be different if we do not use exactly three fluid phases
  if (constants[1].first != "p2")
    FOUR_C_THROW("wrong order in variable vector, p2 not at position 2");
  if (constants[4].first != "S2")
    FOUR_C_THROW("wrong order in variable vector, S2 not at position 5");
  if (constants[6].first != "porosity")
    FOUR_C_THROW("wrong order in variable vector, porosity not at position 7");
  if (variables[0].first != "phi1")
    FOUR_C_THROW("wrong order in variable vector, phi1 (oxygen mass fraction) not at position 1");
  if (variables[1].first != "phi2")
    FOUR_C_THROW("wrong order in variable vector, phi1 (necrotic mass fraction) not at position 2");

  // order is correct
  this->order_checked_ = true;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
double POROMULTIPHASESCATRA::NecrosisLawHeaviside<dim>::Evaluate(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants, const size_t component) const
{
  // Check order (only once since it does not change)
  if (not this->order_checked_) CheckOrder(variables, constants);

  // read variables and constants (order is crucial)
  const double p2 = constants[1].second;
  const double S2 = constants[4].second;
  const double porosity = constants[6].second;
  const double oxy_mass_frac = variables[0].second;
  const double necr_frac = variables[1].second;

  // evaluate heaviside
  const double heaviside_oxy((-(oxy_mass_frac - parameter.w_nl_crit)) > 0. ? 1. : 0.);
  const double macaulay = -parameter.gamma_t_necr * (oxy_mass_frac - parameter.w_nl_crit) /
                          (parameter.w_nl_env - parameter.w_nl_crit) * heaviside_oxy;
  const double heaviside_pres((p2 - parameter.p_t_crit) > 0. ? 1. : 0.);

  // evaluate the function
  const double functval =
      (1.0 - necr_frac) * S2 * porosity * (macaulay + parameter.delta_a_t * heaviside_pres);

  return functval;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
std::vector<double> POROMULTIPHASESCATRA::NecrosisLawHeaviside<dim>::EvaluateDerivative(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants, const size_t component) const
{
  // create derivative vector (should have size of variables)
  std::vector<double> deriv(variables.size(), 0.0);

  if (variables[0].first == "phi1")  // maindiag-derivative
  {
    // read variables and constants (order is crucial)
    const double p2 = constants[1].second;
    const double S2 = constants[4].second;
    const double porosity = constants[6].second;
    const double oxy_mass_frac = variables[0].second;
    const double necr_frac = variables[1].second;

    // evaluate heaviside
    const double heaviside_oxy((-(oxy_mass_frac - parameter.w_nl_crit)) > 0. ? 1. : 0.);
    const double macaulay = -parameter.gamma_t_necr * (oxy_mass_frac - parameter.w_nl_crit) /
                            (parameter.w_nl_env - parameter.w_nl_crit) * heaviside_oxy;
    const double heaviside_pres((p2 - parameter.p_t_crit) > 0. ? 1. : 0.);

    // derivative w.r.t. oxygen mass fraction
    const double oxy_deriv = (1.0 - necr_frac) * porosity * S2 *
                             (-parameter.gamma_t_necr * heaviside_oxy * 1.0 /
                                 (parameter.w_nl_env - parameter.w_nl_crit));
    deriv[0] = oxy_deriv;

    // derivative w.r.t. necrotic cell mass fraction
    const double necro_deriv =
        porosity * S2 * (macaulay + parameter.delta_a_t * heaviside_pres) * (-1.0);
    deriv[1] = necro_deriv;
  }
  else if (variables[0].first == "p1")  // OD-derivative
  {
    // read variables and constants (order is crucial)
    const double S2 = variables[4].second;
    const double porosity = variables[6].second;
    const double p2 = variables[1].second;
    const double oxy_mass_frac = constants[0].second;
    const double necr_frac = constants[1].second;

    // evaluate heaviside
    const double heaviside_oxy((-(oxy_mass_frac - parameter.w_nl_crit)) > 0. ? 1. : 0.);
    const double macaulay = -parameter.gamma_t_necr * (oxy_mass_frac - parameter.w_nl_crit) /
                            (parameter.w_nl_env - parameter.w_nl_crit) * heaviside_oxy;
    const double heaviside_pres((p2 - parameter.p_t_crit) > 0. ? 1. : 0.);

    // derivative w.r.t. tumor cell saturation S2
    const double tc_deriv =
        (1.0 - necr_frac) * porosity * (macaulay + parameter.delta_a_t * heaviside_pres);
    deriv[4] = tc_deriv;

    // derivative w.r.t. porosity
    const double poro_deriv =
        (1.0 - necr_frac) * S2 * (macaulay + parameter.delta_a_t * heaviside_pres);
    deriv[6] = poro_deriv;

    // Note: no pressure derivative, only coupling is with heaviside --> derivative zero
  }
  else
    FOUR_C_THROW("Something went wrong in derivative evaluation of NECROSIS_LAW_HEAVISIDE");

  return deriv;
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
POROMULTIPHASESCATRA::OxygenConsumptionLawHeaviside<dim>::OxygenConsumptionLawHeaviside(
    const std::vector<std::pair<std::string, double>>& funct_params)
    : PoroMultiPhaseScaTraFunction<dim>(), parameter(funct_params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
void POROMULTIPHASESCATRA::OxygenConsumptionLawHeaviside<dim>::CheckOrder(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants) const
{
  // safety check for correct ordering of variables and constants
  // they should have been added in exactly the same way in scatra_ele_calc_multiporo_reaction, but
  // order might be different if we do not use exactly three fluid phases
  if (constants[1].first != "p2")
    FOUR_C_THROW("wrong order in variable vector, p2 not at position 2");
  if (constants[4].first != "S2")
    FOUR_C_THROW("wrong order in variable vector, S2 not at position 5");
  if (constants[6].first != "porosity")
    FOUR_C_THROW("wrong order in variable vector, porosity not at position 7");
  if (variables[0].first != "phi1")
    FOUR_C_THROW("wrong order in variable vector, phi1 (oxygen mass fraction) not at position 1");
  if (variables[1].first != "phi2")
    FOUR_C_THROW("wrong order in variable vector, phi1 (necrotic mass fraction) not at position 2");

  // order is correct
  this->order_checked_ = true;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
double POROMULTIPHASESCATRA::OxygenConsumptionLawHeaviside<dim>::Evaluate(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants, const size_t component) const
{
  // Check order (only once since it does not change)
  if (not this->order_checked_) CheckOrder(variables, constants);

  // read variables and constants (order is crucial)
  const double S2 = constants[4].second;
  const double porosity = constants[6].second;
  const double p2 = constants[1].second;
  const double oxy_mass_frac = variables[0].second;
  const double necr_frac = variables[1].second;

  // evaluate heaviside
  const double heaviside_oxy((oxy_mass_frac - parameter.w_nl_crit) > 0. ? 1. : 0.);
  const double macaulay = parameter.gamma_nl_growth * (oxy_mass_frac - parameter.w_nl_crit) /
                          (parameter.w_nl_env - parameter.w_nl_crit) * heaviside_oxy;
  const double heaviside_pres((parameter.p_t_crit - p2) > 0. ? 1. : 0.);

  // evaluate the function
  const double functval =
      (1.0 - necr_frac) * S2 * porosity *
      (macaulay * heaviside_pres +
          parameter.gamma_0_nl * sin(M_PI / 2.0 * oxy_mass_frac / parameter.w_nl_env));

  return functval;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
std::vector<double> POROMULTIPHASESCATRA::OxygenConsumptionLawHeaviside<dim>::EvaluateDerivative(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants, const size_t component) const
{
  // create derivative vector (should have size of variables)
  std::vector<double> deriv(variables.size(), 0.0);

  if (variables[0].first == "phi1")  // maindiag-derivative
  {
    // read variables and constants (order is crucial)
    const double S2 = constants[4].second;
    const double porosity = constants[6].second;
    const double p2 = constants[1].second;
    const double oxy_mass_frac = variables[0].second;
    const double necr_frac = variables[1].second;

    // evaluate heaviside
    const double heaviside_oxy((oxy_mass_frac - parameter.w_nl_crit) > 0. ? 1. : 0.);
    const double macaulay = parameter.gamma_nl_growth * (oxy_mass_frac - parameter.w_nl_crit) /
                            (parameter.w_nl_env - parameter.w_nl_crit) * heaviside_oxy;
    const double heaviside_pres((parameter.p_t_crit - p2) > 0. ? 1. : 0.);

    // derivative w.r.t. oxygen mass fraction
    const double oxy_deriv = (1.0 - necr_frac) * S2 * porosity *
                             (parameter.gamma_nl_growth * heaviside_oxy * heaviside_pres * 1.0 /
                                     (parameter.w_nl_env - parameter.w_nl_crit) +
                                 parameter.gamma_0_nl * M_PI / 2.0 / parameter.w_nl_env *
                                     cos(M_PI / 2.0 * oxy_mass_frac / parameter.w_nl_env));
    deriv[0] = oxy_deriv;

    // derivative w.r.t. necrotic cell mass fraction
    const double necro_deriv =
        S2 * porosity * (-1.0) *
        (macaulay * heaviside_pres +
            parameter.gamma_0_nl * sin(M_PI / 2.0 * oxy_mass_frac / parameter.w_nl_env));
    deriv[1] = necro_deriv;
  }
  else if (variables[0].first == "p1")  // OD-derivative
  {
    // read variables and constants (order is crucial)
    const double S2 = variables[4].second;
    const double porosity = variables[6].second;
    const double p2 = variables[1].second;
    const double oxy_mass_frac = constants[0].second;
    const double necr_frac = constants[1].second;

    // evaluate heaviside
    const double heaviside_oxy((oxy_mass_frac - parameter.w_nl_crit) > 0. ? 1. : 0.);
    const double macaulay = parameter.gamma_nl_growth * (oxy_mass_frac - parameter.w_nl_crit) /
                            (parameter.w_nl_env - parameter.w_nl_crit) * heaviside_oxy;
    const double heaviside_pres((parameter.p_t_crit - p2) > 0. ? 1. : 0.);

    // derivative w.r.t. tumor cell saturation S2
    const double tc_deriv =
        (1.0 - necr_frac) * porosity *
        (macaulay * heaviside_pres +
            parameter.gamma_0_nl * sin(M_PI / 2.0 * oxy_mass_frac / parameter.w_nl_env));
    deriv[4] = tc_deriv;

    // derivative w.r.t. porosity
    const double poro_deriv =
        (1.0 - necr_frac) * S2 *
        (macaulay * heaviside_pres +
            parameter.gamma_0_nl * sin(M_PI / 2.0 * oxy_mass_frac / parameter.w_nl_env));
    deriv[6] = poro_deriv;
  }
  else
    FOUR_C_THROW(
        "Something went wrong in derivative evaluation of OXYGEN_CONSUMPTION_LAW_HEAVISIDE");

  return deriv;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
POROMULTIPHASESCATRA::TumorGrowthLawHeavisideOxy<dim>::TumorGrowthLawHeavisideOxy(
    const std::vector<std::pair<std::string, double>>& funct_params)
    : PoroMultiPhaseScaTraFunction<dim>(), parameter(funct_params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
void POROMULTIPHASESCATRA::TumorGrowthLawHeavisideOxy<dim>::CheckOrder(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants) const
{
  // safety check for correct ordering of variables and constants
  // they should have been added in exactly the same way in scatra_ele_calc_multiporo_reaction, but
  // order might be different if we do not use exactly three fluid phases
  if (constants[1].first != "p2")
    FOUR_C_THROW("wrong order in variable vector, p2 not at position 2");
  if (constants[4].first != "S2")
    FOUR_C_THROW("wrong order in variable vector, S2 not at position 5");
  if (constants[6].first != "porosity")
    FOUR_C_THROW("wrong order in variable vector, porosity not at position 7");
  if (variables[0].first != "phi1")
    FOUR_C_THROW("wrong order in variable vector, phi1 (oxygen mass fraction) not at position 1");
  if (variables[1].first != "phi2")
    FOUR_C_THROW("wrong order in variable vector, phi1 (necrotic mass fraction) not at position 2");

  // order is correct
  this->order_checked_ = true;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
double POROMULTIPHASESCATRA::TumorGrowthLawHeavisideOxy<dim>::Evaluate(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants, const size_t component) const
{
  // Check order (only once since it does not change)
  if (not this->order_checked_) CheckOrder(variables, constants);

  // read variables and constants (order is crucial)
  const double p2 = constants[1].second;
  const double S2 = constants[4].second;
  const double porosity = constants[6].second;
  const double oxy_mass_frac = variables[0].second;
  const double necr_frac = variables[1].second;

  // evaluate heaviside
  const double heaviside_oxy((oxy_mass_frac - parameter.w_nl_crit) > 0. ? 1. : 0.);
  const double heaviside_pres((parameter.p_t_crit - p2) > 0. ? 1. : 0.);
  const double macaulay = parameter.gamma_T_growth * (oxy_mass_frac - parameter.w_nl_crit) /
                          (parameter.w_nl_env - parameter.w_nl_crit) * heaviside_oxy;

  // evaluate function
  const double functval =
      oxy_mass_frac * S2 * porosity *
      ((macaulay * heaviside_pres) * (1 - necr_frac) - parameter.lambda * necr_frac);

  return functval;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
std::vector<double> POROMULTIPHASESCATRA::TumorGrowthLawHeavisideOxy<dim>::EvaluateDerivative(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants, const size_t component) const
{
  // create derivative vector (should have size of variables)
  std::vector<double> deriv(variables.size(), 0.0);

  if (variables[0].first == "phi1")  // maindiag-derivative
  {
    // read variables and constants (order is crucial)
    const double S2 = constants[4].second;
    const double porosity = constants[6].second;
    const double p2 = constants[1].second;
    const double oxy_mass_frac = variables[0].second;
    const double necr_frac = variables[1].second;

    // evaluate heaviside
    const double heaviside_oxy((oxy_mass_frac - parameter.w_nl_crit) > 0. ? 1. : 0.);
    const double heaviside_pres((parameter.p_t_crit - p2) > 0. ? 1. : 0.);
    const double macaulay = parameter.gamma_T_growth * (oxy_mass_frac - parameter.w_nl_crit) /
                            (parameter.w_nl_env - parameter.w_nl_crit) * heaviside_oxy;

    // derivative w.r.t. oxygen mass fraction
    const double oxy_deriv =
        (parameter.gamma_T_growth * heaviside_oxy * heaviside_pres * 1.0 /
            (parameter.w_nl_env - parameter.w_nl_crit) * (1 - necr_frac)) *
            oxy_mass_frac * S2 * porosity +
        ((macaulay * heaviside_pres) * (1 - necr_frac) - parameter.lambda * necr_frac) * S2 *
            porosity;
    deriv[0] = oxy_deriv;

    // derivative w.r.t. necrotic cell mass fraction
    const double necro_deriv =
        ((macaulay * heaviside_pres) * (-1.0) - parameter.lambda) * oxy_mass_frac * S2 * porosity;
    deriv[1] = necro_deriv;
  }
  else if (variables[0].first == "p1")  // OD-derivative
  {
    // read variables and constants (order is crucial)
    const double S2 = variables[4].second;
    const double porosity = variables[6].second;
    const double p2 = variables[1].second;
    const double oxy_mass_frac = constants[0].second;
    const double necr_frac = constants[1].second;

    // evaluate heaviside
    const double heaviside_oxy((oxy_mass_frac - parameter.w_nl_crit) > 0. ? 1. : 0.);
    const double heaviside_pres((parameter.p_t_crit - p2) > 0. ? 1. : 0.);
    const double macaulay = parameter.gamma_T_growth * (oxy_mass_frac - parameter.w_nl_crit) /
                            (parameter.w_nl_env - parameter.w_nl_crit) * heaviside_oxy;

    // derivative w.r.t. tumor cell saturation S2
    const double tc_deriv =
        ((macaulay * heaviside_pres) * (1.0 - necr_frac) - parameter.lambda * necr_frac) *
        oxy_mass_frac * porosity;
    deriv[4] = tc_deriv;

    // derivative w.r.t. porosity
    const double poro_deriv =
        (macaulay * heaviside_pres * (1.0 - necr_frac) - parameter.lambda * necr_frac) *
        oxy_mass_frac * S2;
    deriv[6] = poro_deriv;
  }
  else
    FOUR_C_THROW("Something went wrong in derivative evaluation of TUMOR_GROWTH_LAW_HEAVISIDE_OXY");

  return deriv;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
POROMULTIPHASESCATRA::TumorGrowthLawHeavisideNecro<dim>::TumorGrowthLawHeavisideNecro(
    const std::vector<std::pair<std::string, double>>& funct_params)
    : PoroMultiPhaseScaTraFunction<dim>(), parameter(funct_params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
void POROMULTIPHASESCATRA::TumorGrowthLawHeavisideNecro<dim>::CheckOrder(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants) const
{
  // safety check for correct ordering of variables and constants
  // they should have been added in exactly the same way in scatra_ele_calc_multiporo_reaction, but
  // order might be different if we do not use exactly three fluid phases
  if (constants[1].first != "p2")
    FOUR_C_THROW("wrong order in variable vector, p2 not at position 2");
  if (constants[4].first != "S2")
    FOUR_C_THROW("wrong order in variable vector, S2 not at position 5");
  if (constants[6].first != "porosity")
    FOUR_C_THROW("wrong order in variable vector, porosity not at position 7");
  if (variables[0].first != "phi1")
    FOUR_C_THROW("wrong order in variable vector, phi1 (oxygen mass fraction) not at position 1");
  if (variables[1].first != "phi2")
    FOUR_C_THROW("wrong order in variable vector, phi1 (necrotic mass fraction) not at position 2");

  // order is correct
  this->order_checked_ = true;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
double POROMULTIPHASESCATRA::TumorGrowthLawHeavisideNecro<dim>::Evaluate(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants, const size_t component) const
{
  // Check order (only once since it does not change)
  if (not this->order_checked_) CheckOrder(variables, constants);

  // read variables and constants (order is crucial)
  const double p2 = constants[1].second;
  const double S2 = constants[4].second;
  const double porosity = constants[6].second;
  const double oxy_mass_frac = variables[0].second;
  const double necr_frac = variables[1].second;

  // evaluate heaviside
  const double heaviside_oxy((oxy_mass_frac - parameter.w_nl_crit) > 0. ? 1. : 0.);
  const double heaviside_pres((parameter.p_t_crit - p2) > 0. ? 1. : 0.);
  const double macaulay = parameter.gamma_T_growth * (oxy_mass_frac - parameter.w_nl_crit) /
                          (parameter.w_nl_env - parameter.w_nl_crit) * heaviside_oxy;

  // evaluate function
  const double functval =
      porosity * S2 *
      (((macaulay * heaviside_pres) * (1 - necr_frac) - parameter.lambda * necr_frac) * necr_frac +
          parameter.lambda * necr_frac);

  return functval;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
std::vector<double> POROMULTIPHASESCATRA::TumorGrowthLawHeavisideNecro<dim>::EvaluateDerivative(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants, const size_t component) const
{
  // create derivative vector (should have size of variables)
  std::vector<double> deriv(variables.size(), 0.0);

  if (variables[0].first == "phi1")  // maindiag-derivative
  {
    // read variables and constants (order is crucial)
    const double p2 = constants[1].second;
    const double S2 = constants[4].second;
    const double porosity = constants[6].second;
    const double oxy_mass_frac = variables[0].second;
    const double necr_frac = variables[1].second;

    // evaluate heaviside
    const double heaviside_oxy((oxy_mass_frac - parameter.w_nl_crit) > 0. ? 1. : 0.);
    const double heaviside_pres((parameter.p_t_crit - p2) > 0. ? 1. : 0.);
    const double macaulay = parameter.gamma_T_growth * (oxy_mass_frac - parameter.w_nl_crit) /
                            (parameter.w_nl_env - parameter.w_nl_crit) * heaviside_oxy;

    // derivative w.r.t. oxygen mass fraction
    const double oxy_deriv = (parameter.gamma_T_growth * heaviside_oxy * heaviside_pres * 1.0 /
                                 (parameter.w_nl_env - parameter.w_nl_crit)) *
                             (1 - necr_frac) * necr_frac * porosity * S2;
    deriv[0] = oxy_deriv;

    // derivative w.r.t. necrotic cell mass fraction
    const double necro_deriv = ((macaulay * heaviside_pres) * (1 - 2.0 * necr_frac) -
                                   2.0 * parameter.lambda * necr_frac + parameter.lambda) *
                               porosity * S2;
    deriv[1] = necro_deriv;
  }
  else if (variables[0].first == "p1")  // OD-derivative
  {
    // read variables and constants (order is crucial)
    const double S2 = variables[4].second;
    const double porosity = variables[6].second;
    const double p2 = variables[1].second;
    const double oxy_mass_frac = constants[0].second;
    const double necr_frac = constants[1].second;

    // evaluate heaviside
    const double heaviside_oxy((oxy_mass_frac - parameter.w_nl_crit) > 0. ? 1. : 0.);
    const double heaviside_pres((parameter.p_t_crit - p2) > 0. ? 1. : 0.);
    const double macaulay = parameter.gamma_T_growth * (oxy_mass_frac - parameter.w_nl_crit) /
                            (parameter.w_nl_env - parameter.w_nl_crit) * heaviside_oxy;

    // derivative w.r.t. tumor cell saturation S2
    const double tc_deriv =
        porosity * (((macaulay * heaviside_pres) * (1 - necr_frac) - parameter.lambda * necr_frac) *
                           necr_frac +
                       parameter.lambda * necr_frac);
    deriv[4] = tc_deriv;

    // derivative w.r.t. porosity
    const double poro_deriv =
        S2 * (((macaulay * heaviside_pres) * (1 - necr_frac) - parameter.lambda * necr_frac) *
                     necr_frac +
                 parameter.lambda * necr_frac);
    deriv[6] = poro_deriv;

    // Note: no pressure derivative, only coupling is with heaviside --> derivative zero
  }
  else
    FOUR_C_THROW(
        "Something went wrong in derivative evaluation of TUMOR_GROWTH_LAW_HEAVISIDE_NECRO");

  return deriv;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
POROMULTIPHASESCATRA::OxygenTransvascularExchangeLawCont<dim>::OxygenTransvascularExchangeLawCont(
    const std::vector<std::pair<std::string, double>>& funct_params)
    : PoroMultiPhaseScaTraFunction<dim>(), parameter(funct_params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
void POROMULTIPHASESCATRA::OxygenTransvascularExchangeLawCont<dim>::CheckOrder(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants) const
{
  // safety check for correct ordering of variables and constants
  // they should have been added in exactly the same way in scatra_ele_calc_multiporo_reaction, but
  // order might be different if we do not use exactly three fluid phases
  if (constants[7].first != "VF1")
    FOUR_C_THROW("wrong order in variable vector, porosity not at position 8");
  if (variables[0].first != "phi1")
    FOUR_C_THROW("wrong order in variable vector, phi1 (oxygen mass fraction) not at position 1");
  if (variables[1].first != "phi2")
    FOUR_C_THROW("wrong order in variable vector, phi2 (necrotic mass fraction) not at position 2");
  if (variables[2].first != "phi3")
    FOUR_C_THROW("wrong order in variable vector, phi3 (necrotic mass fraction) not at position 3");
  if (variables[3].first != "phi4")
    FOUR_C_THROW("wrong order in variable vector, phi4 (necrotic mass fraction) not at position 4");

  // order is correct
  this->order_checked_ = true;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
double POROMULTIPHASESCATRA::OxygenTransvascularExchangeLawCont<dim>::Evaluate(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants, const size_t component) const
{
  // Check order (only once since it does not change)
  if (not this->order_checked_) CheckOrder(variables, constants);

  const double fac_if = parameter.rho_oxy / parameter.rho_if * parameter.alpha_IF;

  // read variables and constants (order is crucial)
  double VF1 = constants[7].second;
  double oxy_mass_frac_if = variables[0].second;
  double oxy_mass_frac_nv = variables[3].second;

  double Pb = 0.0;
  double CaO2 = oxy_mass_frac_nv * parameter.rho_bl / parameter.rho_oxy;
  // safety check --> should not be larger than CaO2_max, which already correponds to partial
  // pressures of ~250Pa
  CaO2 = std::max(0.0, std::min(CaO2, 1.0 * parameter.CaO2_max));
  POROMULTIPHASESCATRA::UTILS::GetOxyPartialPressureFromConcentration<double>(
      Pb, CaO2, parameter.CaO2_max, parameter.Pb50, parameter.n, parameter.alpha_bl_eff);

  // evaluate function
  const double heaviside_oxy((Pb - oxy_mass_frac_if / fac_if) > 0. ? 1. : 0.);
  const double functval =
      parameter.gammarhoSV * heaviside_oxy * (Pb - oxy_mass_frac_if / fac_if) * VF1;

  return functval;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
std::vector<double>
POROMULTIPHASESCATRA::OxygenTransvascularExchangeLawCont<dim>::EvaluateDerivative(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants, const size_t component) const
{
  // create derivative vector (should have size of variables)
  std::vector<double> deriv(variables.size(), 0.0);

  const double fac_if = parameter.rho_oxy / parameter.rho_if * parameter.alpha_IF;

  double VF1 = 0.0, oxy_mass_frac_if = 0.0;

  // define Fad object for evaluation
  using FAD = Sacado::Fad::DFad<double>;
  FAD oxy_mass_frac_nv = 0.0;
  oxy_mass_frac_nv.diff(0, 1);       // independent variable 0 out of a total of 1
  if (variables[0].first == "phi1")  // maindiag-derivative
  {
    // read variables and constants (order is crucial)
    VF1 = constants[7].second;
    oxy_mass_frac_if = variables[0].second;
    oxy_mass_frac_nv.val() = variables[3].second;
  }
  else if (variables[0].first == "p1")  // OD-derivative
  {
    // read variables and constants (order is crucial)
    VF1 = variables[7].second;
    oxy_mass_frac_if = constants[0].second;
    oxy_mass_frac_nv.val() = constants[3].second;
  }
  else
    FOUR_C_THROW(
        "Something went wrong in derivative evaluation of OXYGEN_TRANSVASCULAR_EXCHANGE_LAW_CONT");

  FAD Pb = 0.0;
  FAD CaO2 = oxy_mass_frac_nv * parameter.rho_bl / parameter.rho_oxy;
  // safety check --> should not be larger than CaO2_max, which already correponds to partial
  // pressures of ~250Pa
  CaO2 = std::max(0.0, std::min(CaO2, 1.0 * parameter.CaO2_max));
  POROMULTIPHASESCATRA::UTILS::GetOxyPartialPressureFromConcentration<FAD>(
      Pb, CaO2, parameter.CaO2_max, parameter.Pb50, parameter.n, parameter.alpha_bl_eff);
  const double heaviside_oxy((Pb - oxy_mass_frac_if / fac_if) > 0. ? 1. : 0.);

  if (variables[0].first == "phi1")  // maindiag-derivative
  {
    deriv[0] = parameter.gammarhoSV * VF1 * (-1.0 / fac_if) * heaviside_oxy;
    deriv[3] = parameter.gammarhoSV * VF1 * (Pb.fastAccessDx(0)) * heaviside_oxy;
  }
  else if (variables[0].first == "p1")  // OD-derivative
  {
    deriv[7] = parameter.gammarhoSV * (Pb.val() - oxy_mass_frac_if / fac_if) * heaviside_oxy;
  }
  else
    FOUR_C_THROW(
        "Something went wrong in derivative evaluation of OXYGEN_TRANSVASCULAR_EXCHANGE_LAW_CONT");

  return deriv;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
POROMULTIPHASESCATRA::OxygenTransvascularExchangeLawDisc<dim>::OxygenTransvascularExchangeLawDisc(
    const std::vector<std::pair<std::string, double>>& funct_params)
    : PoroMultiPhaseScaTraFunction<dim>(), parameter(funct_params), pos_oxy_art_(-1), pos_diam_(-1)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
void POROMULTIPHASESCATRA::OxygenTransvascularExchangeLawDisc<dim>::CheckOrder(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants) const
{
  // safety check for correct ordering of variables and constants
  // they should have been added in exactly the same way in scatra_ele_calc_multiporo_reaction, but
  // order might be different if we do not use exactly three fluid phases
  if (variables[0].first != "phi1")
    FOUR_C_THROW("wrong order in variable vector, phi1 (oxygen mass fraction) not at position 1");

  // oxygen in artery
  // we have no neo-vasculature --> at position 2, species 1: oxy in IF, species 2: NTC
  if (variables[2].first == "phi_art1") pos_oxy_art_ = 2;
  // we have no neo-vasculature --> at position 4, species 1: oxy in IF, species 2: NTC,
  // species 3: TAF, species 4: oxy in NV
  if (variables.size() >= 5 && variables[4].first == "phi_art1") pos_oxy_art_ = 4;
  if (pos_oxy_art_ == -1) FOUR_C_THROW("cannot find position of oxygen in arteries");

  // fluid variables
  if (constants[1].first != "p2")
    FOUR_C_THROW("wrong order in constants vector, p2 (pressure of tumor cells) not at position 2");
  if (constants[4].first != "S2")
    FOUR_C_THROW(
        "wrong order in constants vector, S2 (saturation of tumor cells) not at position 5");

  // diameter
  if (constants[8].first == "D") pos_diam_ = 8;
  if (constants.size() >= 11 && constants[10].first == "D") pos_diam_ = 10;
  if (pos_diam_ == -1) FOUR_C_THROW("cannot find position of artery diameter");

  // order is correct
  this->order_checked_ = true;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
double POROMULTIPHASESCATRA::OxygenTransvascularExchangeLawDisc<dim>::Evaluate(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants, const size_t component) const
{
  // Check order (only once since it does not change)
  if (not this->order_checked_) CheckOrder(variables, constants);

  const double fac_if = parameter.rho_oxy / parameter.rho_if * parameter.alpha_IF;

  // read variables and constants (order is crucial)
  double oxy_mass_frac_if = variables[0].second;
  double oxy_mass_frac_nv = variables[pos_oxy_art_].second;
  const double D = constants[pos_diam_].second;
  const double S2 = constants[4].second;

  double Pb = 0.0;
  double CaO2 = oxy_mass_frac_nv * parameter.rho_bl / parameter.rho_oxy;
  // safety check --> should not be larger than CaO2_max, which already correponds to partial
  // pressures of ~250Pa
  CaO2 = std::max(0.0, std::min(CaO2, 1.0 * parameter.CaO2_max));
  POROMULTIPHASESCATRA::UTILS::GetOxyPartialPressureFromConcentration<double>(
      Pb, CaO2, parameter.CaO2_max, parameter.Pb50, parameter.n, parameter.alpha_bl_eff);

  // evaluate function
  const double heaviside_S2 = ((S2 - parameter.S2_max) > 0. ? 0. : 1.);
  const double functval =
      parameter.gammarho * M_PI * D * (Pb - oxy_mass_frac_if / fac_if) * heaviside_S2;

  return functval;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
std::vector<double>
POROMULTIPHASESCATRA::OxygenTransvascularExchangeLawDisc<dim>::EvaluateDerivative(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants, const size_t component) const
{
  // create derivative vector (should have size of variables)
  std::vector<double> deriv(variables.size(), 0.0);

  const double fac_if = parameter.rho_oxy / parameter.rho_if * parameter.alpha_IF;

  // define Fad object for evaluation
  using FAD = Sacado::Fad::DFad<double>;
  FAD oxy_mass_frac_nv = 0.0;
  oxy_mass_frac_nv.diff(0, 1);  // independent variable 0 out of a total of 1

  // read variables and constants (order is crucial)
  oxy_mass_frac_nv.val() = variables[pos_oxy_art_].second;
  const double D = constants[pos_diam_].second;
  const double S2 = constants[4].second;

  FAD Pb = 0.0;
  FAD CaO2 = oxy_mass_frac_nv * parameter.rho_bl / parameter.rho_oxy;
  // safety check --> should not be larger than CaO2_max, which already correponds to partial
  // pressures of ~250Pa
  CaO2 = std::max(0.0, std::min(CaO2, 1.0 * parameter.CaO2_max));
  POROMULTIPHASESCATRA::UTILS::GetOxyPartialPressureFromConcentration<FAD>(
      Pb, CaO2, parameter.CaO2_max, parameter.Pb50, parameter.n, parameter.alpha_bl_eff);

  // evaluate function
  const double heaviside_S2 = ((S2 - parameter.S2_max) > 0. ? 0. : 1.);

  deriv[0] = parameter.gammarho * M_PI * D * (-1.0 / fac_if) * heaviside_S2;
  deriv[pos_oxy_art_] = parameter.gammarho * M_PI * D * (Pb.fastAccessDx(0)) * heaviside_S2;

  return deriv;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
POROMULTIPHASESCATRA::LungOxygenExchangeLaw<dim>::LungOxygenExchangeLaw(
    const std::vector<std::pair<std::string, double>>& funct_params)
    : PoroMultiPhaseScaTraFunction<dim>(), parameter(funct_params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
void POROMULTIPHASESCATRA::LungOxygenExchangeLaw<dim>::CheckOrder(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants) const
{
  // safety check for correct ordering of variables and constants
  if (variables[0].first == "phi1")
  {
    if (constants[0].first != "p1")
      FOUR_C_THROW("wrong order in constants vector, P1 (Pressure of air) not at position 1");
    if (constants[3].first != "VF1")
      FOUR_C_THROW(
          "wrong order in constants vector, VF1 (volume fraction of additional porous network "
          "(blood phase)) not at position 4");
    if (variables[1].first != "phi2")
      FOUR_C_THROW(
          "wrong order in variable vector, phi2 (oxygen mass fraction in blood) not at position 2");
  }
  else if (variables[0].first == "p1")
  {
    if (variables[3].first != "VF1")
    {
      FOUR_C_THROW(
          "wrong order in variable vector, VF1 (volume fraction of additional porous network "
          "(blood)) not at position 4");
    }
    if (constants[0].first != "phi1")
      FOUR_C_THROW(
          "wrong order in variable vector, phi1 (oxygen mass fraction in air) not at position 1");
    if (constants[1].first != "phi2")
      FOUR_C_THROW(
          "wrong order in variable vector, phi2 (oxygen mass fraction in blood) not at position 2");
  }
  else
  {
    FOUR_C_THROW("Variable <%s> not supported on position 0. Wrong order in variable vector! ",
        variables[0].first.c_str());
  }

  // order is correct
  this->order_checked_ = true;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
double POROMULTIPHASESCATRA::LungOxygenExchangeLaw<dim>::Evaluate(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants, const size_t component) const
{
  // Check order of variables and constants vector only once (since it does not change)
  if (not this->order_checked_) CheckOrder(variables, constants);

    // In debug mode, check order of variables and constants vector on every call
#ifdef FOUR_C_ENABLE_ASSERTIONS
  CheckOrder(variables, constants);
#endif

  // read variables (order is crucial)
  const double oxy_mass_frac_air = variables[0].second;
  const double oxy_mass_frac_bl = variables[1].second;

  // read constants (order is crucial)
  const double P_air = constants[0].second;
  const double volfrac_blood = constants[3].second;

  // partial pressure of oxygen in air
  const double P_oA =
      oxy_mass_frac_air * (P_air + parameter.P_atmospheric) * parameter.rho_air / parameter.rho_oxy;

  // CoB_total is total concentration of oxygen in blood (physically dissolved and bound to
  // hemoglobin)
  const double CoB_total = oxy_mass_frac_bl * parameter.rho_bl / parameter.rho_oxy;

  // partial pressure of oxygen in blood
  double P_oB = 0.0;

  // Calculate partial pressure of oxygen in blood
  POROMULTIPHASESCATRA::UTILS::GetOxyPartialPressureFromConcentration<double>(
      P_oB, CoB_total, parameter.NC_Hb, parameter.P_oB50, parameter.n, parameter.alpha_oxy);

  // evaluate function
  const double functval = parameter.rho_oxy * parameter.DiffAdVTLC * parameter.alpha_oxy *
                          (volfrac_blood / parameter.volfrac_blood_ref) * (P_oA - P_oB);

  return functval;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
std::vector<double> POROMULTIPHASESCATRA::LungOxygenExchangeLaw<dim>::EvaluateDerivative(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants, const size_t component) const
{
// In debug mode, check order of variables and constants vector on every call
#ifdef FOUR_C_ENABLE_ASSERTIONS
  CheckOrder(variables, constants);
#endif
  // Check order of variables and constants vector only once (since it does not change)
  if (not this->order_checked_) CheckOrder(variables, constants);

  // create derivative vector (should have size of variables)
  std::vector<double> deriv(variables.size(), 0.0);

  // define Fad object for evaluation
  using FAD = Sacado::Fad::DFad<double>;
  FAD oxy_mass_frac_bl = 0.0;
  oxy_mass_frac_bl.diff(0, 1);  // independent variable 0 out of a total of 1

  double oxy_mass_frac_air = 0.0, P_air = 0.0, volfrac_blood = 0.0;

  if (variables[0].first == "phi1")  // maindiag-derivative
  {
    // read variables and constants (order is crucial)
    oxy_mass_frac_air = variables[0].second;
    oxy_mass_frac_bl.val() = variables[1].second;
    P_air = constants[0].second;
    volfrac_blood = constants[3].second;
  }
  else if (variables[0].first == "p1")  // OD-derivative
  {
    // read variables and constants (order is crucial)
    oxy_mass_frac_air = constants[0].second;
    oxy_mass_frac_bl.val() = constants[1].second;
    P_air = variables[0].second;
    volfrac_blood = variables[3].second;
  }
  else
    FOUR_C_THROW("Derivative w.r.t. <%s> not supported in LUNG_OXYGEN_EXCHANGE_LAW.",
        variables[0].first.c_str());

  // volfrac relation
  const double volfrac_relation = (volfrac_blood / parameter.volfrac_blood_ref);

  FAD P_oB = 0.0;
  FAD C_oB_total = oxy_mass_frac_bl * parameter.rho_bl / parameter.rho_oxy;

  POROMULTIPHASESCATRA::UTILS::GetOxyPartialPressureFromConcentration<FAD>(
      P_oB, C_oB_total, parameter.NC_Hb, parameter.P_oB50, parameter.n, parameter.alpha_oxy);


  if (variables[0].first == "phi1")  // maindiag-derivative
  {
    deriv[0] = parameter.rho_oxy * parameter.DiffAdVTLC * volfrac_relation * parameter.alpha_oxy *
               ((P_air + parameter.P_atmospheric) * parameter.rho_air / parameter.rho_oxy);
    deriv[1] = parameter.rho_oxy * parameter.DiffAdVTLC * volfrac_relation * parameter.alpha_oxy *
               (-1.0) * P_oB.fastAccessDx(0);
  }
  else if (variables[0].first == "p1")  // OD-derivative
  {
    deriv[0] = (parameter.rho_oxy * parameter.DiffAdVTLC * volfrac_relation * parameter.alpha_oxy) *
               ((oxy_mass_frac_air * parameter.rho_air) /
                   parameter.rho_oxy);  // derivative wrt P_air (dFunc/dP_oA * dP_oA/P_air)
                                        // partial pressure of oxygen in air
    const double P_oA = oxy_mass_frac_air * (P_air + parameter.P_atmospheric) * parameter.rho_air /
                        parameter.rho_oxy;
    deriv[3] = parameter.rho_oxy * parameter.DiffAdVTLC * parameter.alpha_oxy *
               (1 / parameter.volfrac_blood_ref) * (P_oA - P_oB.val());
  }
  else
    FOUR_C_THROW("Derivative w.r.t. <%s> not supported in LUNG_OXYGEN_EXCHANGE_LAW.",
        variables[0].first.c_str());

  return deriv;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
POROMULTIPHASESCATRA::LungCarbonDioxideExchangeLaw<dim>::LungCarbonDioxideExchangeLaw(
    const std::vector<std::pair<std::string, double>>& funct_params)
    : PoroMultiPhaseScaTraFunction<dim>(), parameter(funct_params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
void POROMULTIPHASESCATRA::LungCarbonDioxideExchangeLaw<dim>::CheckOrder(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants) const
{
  // safety check for correct ordering of variables and constants
  if (variables[0].first == "phi1")
  {
    if (constants[0].first != "p1")
      FOUR_C_THROW("wrong order in constants vector, P1 (Pressure of air) not at position 1");
    if (constants[1].first != "S1")
      FOUR_C_THROW("wrong order in constants vector, S1 (Saturation of air) not at position 2");
    if (constants[3].first != "VF1")
      FOUR_C_THROW("wrong order in constants vector, VF1 (volfrac 1) not at position 4");
    if (variables[1].first != "phi2")
      FOUR_C_THROW(
          "wrong order in variable vector, phi2 (oxygen mass fraction in blood) not at position 2");
  }
  else if (variables[0].first == "p1")
  {
    if (constants[0].first != "phi1")
      FOUR_C_THROW(
          "wrong order in variable vector, phi1 (oxygen mass fraction in air) not at position 1");
    if (constants[1].first != "phi2")
      FOUR_C_THROW(
          "wrong order in variable vector, phi2 (oxygen mass fraction in blood) not at position 2");
    if (constants[2].first != "phi3")
      FOUR_C_THROW(
          "wrong order in variable vector, phi3 (oxygen mass fraction in blood) not at position 3");
    if (constants[3].first != "phi4")
      FOUR_C_THROW(
          "wrong order in variable vector, phi4 (oxygen mass fraction in blood) not at position 4");
  }
  else
  {
    FOUR_C_THROW("Variable <%s> not supported on position 0. Wrong order in variable vector! ",
        variables[0].first.c_str());
  }

  // order is correct
  this->order_checked_ = true;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
double POROMULTIPHASESCATRA::LungCarbonDioxideExchangeLaw<dim>::Evaluate(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants, const size_t component) const
{
  // Check order of variables and constants vector only once (since it does not change)
  if (not this->order_checked_) CheckOrder(variables, constants);

    // In debug mode, check order of variables and constants vector on every call
#ifdef FOUR_C_ENABLE_ASSERTIONS
  CheckOrder(variables, constants);
#endif

  // read variables (order is crucial)
  const double O2_mass_frac_bl = variables[1].second;
  const double CO2_mass_frac_air = variables[2].second;
  const double CO2_mass_frac_bl = variables[3].second;

  // read constants (order is crucial)
  const double P_air = constants[0].second;
  const double volfrac_blood = constants[3].second;

  // partial pressure of carbon dioxide in air
  const double P_CO2A =
      CO2_mass_frac_air * (P_air + parameter.P_atmospheric) * parameter.rho_air / parameter.rho_CO2;

  // CoB_total is total concentration of oxygen in blood (physically dissolved and bound to
  // hemoglobin)
  const double CoB_total = O2_mass_frac_bl * parameter.rho_bl / parameter.rho_oxy;

  // partial pressure of oxygen in blood
  double P_O2B = 0.0;

  // Calculate partial pressure of oxygen in blood
  POROMULTIPHASESCATRA::UTILS::GetOxyPartialPressureFromConcentration<double>(
      P_O2B, CoB_total, parameter.NC_Hb, parameter.P_oB50, parameter.n, parameter.alpha_oxy);

  // saturation of hemoglobin with oxygen from hill equation
  const double SO2 =
      pow(P_O2B, parameter.n) / (pow(P_O2B, parameter.n) + pow(parameter.P_oB50, parameter.n));

  // temporary help variable for calculating partial pressure of carbon dioxide in blood
  const double temp =
      (1.0 - (0.02924 * parameter.C_Hb) / ((2.244 - 0.422 * SO2) * (8.740 - parameter.pH))) *
      0.0301 * 2.226 * (1 + pow(10, parameter.pH - 6.1));

  // partial pressure of carbon dioxide in blood
  double P_CO2B = (CO2_mass_frac_bl * parameter.rho_bl) / (parameter.rho_CO2 * temp);

  // scaling of P_CO2B to get from mmHg to the used pressure unit in the input file
  P_CO2B *= parameter.ScalingFormmHg;

  // evaluate function
  const double functval = parameter.rho_CO2 * parameter.DiffsolAdVTLC *
                          (volfrac_blood / parameter.volfrac_blood_ref) * (P_CO2B - P_CO2A);

  return functval;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
std::vector<double> POROMULTIPHASESCATRA::LungCarbonDioxideExchangeLaw<dim>::EvaluateDerivative(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants, const size_t component) const
{
// In debug mode, check order of variables and constants vector on every call
#ifdef FOUR_C_ENABLE_ASSERTIONS
  CheckOrder(variables, constants);
#endif

  // Check order of variables and constants vector only once (since it does not change)
  if (not this->order_checked_) CheckOrder(variables, constants);

  // create derivative vector (should have size of variables)
  std::vector<double> deriv(variables.size(), 0.0);

  // define Fad object for evaluation
  using FAD = Sacado::Fad::DFad<double>;
  FAD O2_mass_frac_bl = 0.0;
  O2_mass_frac_bl.diff(0, 1);

  double CO2_mass_frac_air = 0.0, P_air = 0.0, CO2_mass_frac_bl = 0.0, volfrac_blood = 0.0;

  if (variables[0].first == "phi1")  // maindiag-derivative
  {
    // read variables and constants (order is crucial)
    O2_mass_frac_bl.val() = variables[1].second;
    CO2_mass_frac_air = variables[2].second;
    CO2_mass_frac_bl = variables[3].second;
    P_air = constants[0].second;
    volfrac_blood = constants[3].second;
  }
  else if (variables[0].first == "p1")  // OD-derivative
  {
    // read variables and constants (order is crucial)
    O2_mass_frac_bl.val() = constants[1].second;
    CO2_mass_frac_air = constants[2].second;
    CO2_mass_frac_bl = constants[3].second;
    P_air = variables[0].second;
    volfrac_blood = variables[3].second;
  }
  else
    FOUR_C_THROW("Derivative w.r.t. <%s> not supported in LUNG_CARBONDIOXIDE_EXCHANGE_LAW.",
        variables[0].first.c_str());

  // volfrac relation
  const double volfrac_relation = (volfrac_blood / parameter.volfrac_blood_ref);

  FAD P_O2B = 0.0;
  FAD C_oB_total = O2_mass_frac_bl * parameter.rho_bl / parameter.rho_oxy;

  POROMULTIPHASESCATRA::UTILS::GetOxyPartialPressureFromConcentration<FAD>(
      P_O2B, C_oB_total, parameter.NC_Hb, parameter.P_oB50, parameter.n, parameter.alpha_oxy);

  // saturation of hemoglobin with oxygen from hill equation
  const double SO2 = pow(P_O2B.val(), parameter.n) /
                     (pow(P_O2B.val(), parameter.n) + pow(parameter.P_oB50, parameter.n));

  // temporary help variable for calculating partial pressure of carbon dioxide in blood
  const double temp =
      (1.0 - (0.02924 * parameter.C_Hb) / ((2.244 - 0.422 * SO2) * (8.740 - parameter.pH))) *
      0.0301 * 2.226 * (1.0 + pow(10.0, parameter.pH - 6.1));


  if (variables[0].first == "phi1")  // maindiag-derivative
  {
    // linearization w.r.t. phi2 (oxygen in blood) = dMassexchangeCO2/dwO2B = dMassexchangeCO2/dPCO2
    // * dPCO2/dSO2 * dSO2/dPO2B * dPO2B/dwO2B
    double dMassexchangeCO2dPCO2 =
        parameter.rho_CO2 * parameter.DiffsolAdVTLC * volfrac_relation * parameter.ScalingFormmHg;
    double dPCO2dSO2 = CO2_mass_frac_bl * (parameter.rho_bl / parameter.rho_CO2) * pow(temp, -2.0) *
                       0.0301 * (1.0 + pow(10.0, parameter.pH - 6.10)) * 2.226 *
                       (0.02924 * parameter.C_Hb) / ((8.740 - parameter.pH)) *
                       pow(2.244 - 0.422 * SO2, -2.0) * 0.422;
    double dSO2dPO2B =
        parameter.n * pow(P_O2B.val(), parameter.n - 1.0) *
            pow(pow(P_O2B.val(), parameter.n) + pow(parameter.P_oB50, parameter.n), -1.0) -
        pow(pow(P_O2B.val(), parameter.n) + pow(parameter.P_oB50, parameter.n), -2.0) *
            pow(P_O2B.val(), parameter.n) * parameter.n * pow(P_O2B.val(), parameter.n - 1.0);
    double dPO2BdwO2B = P_O2B.fastAccessDx(0);
    deriv[1] = dMassexchangeCO2dPCO2 * dPCO2dSO2 * dSO2dPO2B * dPO2BdwO2B;

    deriv[2] = (-1.0) * parameter.rho_CO2 * parameter.DiffsolAdVTLC * volfrac_relation *
               ((P_air + parameter.P_atmospheric) * parameter.rho_air / parameter.rho_CO2);

    deriv[3] = parameter.rho_CO2 * parameter.DiffsolAdVTLC * volfrac_relation *
               parameter.ScalingFormmHg * parameter.rho_bl / (parameter.rho_CO2 * temp);
  }
  else if (variables[0].first == "p1")  // OD-derivative
  {
    // derivative w.r.t. P_air (dFunc/dP_CO2A * dP_CO2A/P_air)
    deriv[0] = (-1.0) * (parameter.rho_CO2 * parameter.DiffsolAdVTLC * volfrac_relation) *
               ((CO2_mass_frac_air * parameter.rho_air) / parameter.rho_CO2);

    // partial pressure of carbon dioxide in air
    const double P_CO2A = CO2_mass_frac_air * (P_air + parameter.P_atmospheric) *
                          parameter.rho_air / parameter.rho_CO2;
    // partial pressure of carbon dioxide in blood
    double P_CO2B = ((CO2_mass_frac_bl * parameter.rho_bl / parameter.rho_CO2) / temp) *
                    parameter.ScalingFormmHg;
    deriv[3] = parameter.rho_CO2 * parameter.DiffsolAdVTLC * (1 / parameter.volfrac_blood_ref) *
               (P_CO2B - P_CO2A);
  }
  else
    FOUR_C_THROW("Derivative w.r.t. <%s> not supported in LUNG_CARBONDIOXIDE_EXCHANGE_LAW.",
        variables[0].first.c_str());

  return deriv;
}

// explicit instantiations

template class POROMULTIPHASESCATRA::TumorGrowthLawHeaviside<1>;
template class POROMULTIPHASESCATRA::TumorGrowthLawHeaviside<2>;
template class POROMULTIPHASESCATRA::TumorGrowthLawHeaviside<3>;

template class POROMULTIPHASESCATRA::NecrosisLawHeaviside<1>;
template class POROMULTIPHASESCATRA::NecrosisLawHeaviside<2>;
template class POROMULTIPHASESCATRA::NecrosisLawHeaviside<3>;

template class POROMULTIPHASESCATRA::OxygenConsumptionLawHeaviside<1>;
template class POROMULTIPHASESCATRA::OxygenConsumptionLawHeaviside<2>;
template class POROMULTIPHASESCATRA::OxygenConsumptionLawHeaviside<3>;

template class POROMULTIPHASESCATRA::TumorGrowthLawHeavisideOxy<1>;
template class POROMULTIPHASESCATRA::TumorGrowthLawHeavisideOxy<2>;
template class POROMULTIPHASESCATRA::TumorGrowthLawHeavisideOxy<3>;

template class POROMULTIPHASESCATRA::TumorGrowthLawHeavisideNecro<1>;
template class POROMULTIPHASESCATRA::TumorGrowthLawHeavisideNecro<2>;
template class POROMULTIPHASESCATRA::TumorGrowthLawHeavisideNecro<3>;

template class POROMULTIPHASESCATRA::OxygenTransvascularExchangeLawCont<1>;
template class POROMULTIPHASESCATRA::OxygenTransvascularExchangeLawCont<2>;
template class POROMULTIPHASESCATRA::OxygenTransvascularExchangeLawCont<3>;

template class POROMULTIPHASESCATRA::OxygenTransvascularExchangeLawDisc<1>;
template class POROMULTIPHASESCATRA::OxygenTransvascularExchangeLawDisc<2>;
template class POROMULTIPHASESCATRA::OxygenTransvascularExchangeLawDisc<3>;

template class POROMULTIPHASESCATRA::LungOxygenExchangeLaw<1>;
template class POROMULTIPHASESCATRA::LungOxygenExchangeLaw<2>;
template class POROMULTIPHASESCATRA::LungOxygenExchangeLaw<3>;

template class POROMULTIPHASESCATRA::LungCarbonDioxideExchangeLaw<1>;
template class POROMULTIPHASESCATRA::LungCarbonDioxideExchangeLaw<2>;
template class POROMULTIPHASESCATRA::LungCarbonDioxideExchangeLaw<3>;

FOUR_C_NAMESPACE_CLOSE
