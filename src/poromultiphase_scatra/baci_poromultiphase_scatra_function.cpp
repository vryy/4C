/*----------------------------------------------------------------------*/
/*! \file
 \brief Managing and evaluating of (reaction) functions for poromultiphase_scatra
        problems

\level 3

    *----------------------------------------------------------------------*/

#include "baci_poromultiphase_scatra_function.H"
#include "baci_poromultiphase_scatra_utils.H"
#include "baci_utils_fad.H"
#include <Teuchos_RCP.hpp>
#include "baci_lib_linedefinition.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
POROMULTIPHASESCATRA::PoroMultiPhaseScaTraFunction<dim>::PoroMultiPhaseScaTraFunction()
    : order_checked_(false)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::AddValidPoroFunctionLines(Teuchos::RCP<DRT::INPUT::Lines> lines)
{
  DRT::INPUT::LineDefinition poromultiphasescatra_funct;
  poromultiphasescatra_funct.AddNamedString("POROMULTIPHASESCATRA_FUNCTION")
      .AddOptionalNamedInt("NUMPARAMS")
      .AddOptionalNamedPairOfStringAndDoubleVector("PARAMS", "NUMPARAMS");

  lines->Add(poromultiphasescatra_funct);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
Teuchos::RCP<DRT::UTILS::FunctionOfAnything> POROMULTIPHASESCATRA::TryCreatePoroFunction(
    const std::vector<Teuchos::RCP<DRT::INPUT::LineDefinition>>& function_line_defs)
{
  if (function_line_defs.size() != 1) return Teuchos::null;

  const auto& function_lin_def = function_line_defs.front();

  if (function_lin_def->HaveNamed("POROMULTIPHASESCATRA_FUNCTION"))
  {
    std::string type;
    function_lin_def->ExtractString("POROMULTIPHASESCATRA_FUNCTION", type);

    std::vector<std::pair<std::string, double>> params;
    if (function_lin_def->HaveNamed("PARAMS"))
      function_lin_def->ExtractPairOfStringAndDoubleVector("PARAMS", params);

    return CreatePoroFunction<dim>(type, params);
  }
  else
  {
    return Teuchos::RCP<DRT::UTILS::FunctionOfAnything>(nullptr);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
Teuchos::RCP<DRT::UTILS::FunctionOfAnything> POROMULTIPHASESCATRA::CreatePoroFunction(
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
    return Teuchos::rcp(new POROMULTIPHASESCATRA::OxygenTransvascularExchangeLawCont<dim>(params));
  }
  else if (type == "OXYGEN_TRANSVASCULAR_EXCHANGE_LAW_DISC")
  {
    return Teuchos::rcp(new POROMULTIPHASESCATRA::OxygenTransvascularExchangeLawDisc<dim>(params));
  }
  else if (type == "LUNG_OXYGEN_EXCHANGE_LAW")
  {
    return Teuchos::rcp(new POROMULTIPHASESCATRA::LungOxygenExchangeLaw<dim>(params));
  }
  else
  {
    dserror("Wrong type of POROMULTIPHASESCATRA_FUNCTION");
    return Teuchos::RCP<DRT::UTILS::FunctionOfAnything>(nullptr);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
POROMULTIPHASESCATRA::TumorGrowthLawHeaviside<dim>::TumorGrowthLawHeaviside(
    std::vector<std::pair<std::string, double>> funct_params)
    : PoroMultiPhaseScaTraFunction<dim>()
{
  // Check size
  if (funct_params.size() != 5)
  {
    dserror(
        "Wrong size of funct_params for TUMOR_GROWTH_LAW_HEAVISIDE, it should have exactly\n"
        "5 funct_params (in this order) gamma_T_growth, w_nl_crit, w_nl_env, lambda and p_t_crit");
  }

  // Check correct naming and order of funct_params
  if (funct_params[0].first != "gamma_T_growth")
    dserror("First parameter for TUMOR_GROWTH_LAW_HEAVISIDE has to be gamma_T_growth");

  if (funct_params[1].first != "w_nl_crit")
    dserror("Second parameter for TUMOR_GROWTH_LAW_HEAVISIDE has to be w_nl_crit");

  if (funct_params[2].first != "w_nl_env")
    dserror("Third parameter for TUMOR_GROWTH_LAW_HEAVISIDE has to be w_nl_env");

  if (funct_params[3].first != "lambda")
    dserror("Fourth parameter for TUMOR_GROWTH_LAW_HEAVISIDE has to be lambda");

  if (funct_params[4].first != "p_t_crit")
    dserror("Fifth parameter for TUMOR_GROWTH_LAW_HEAVISIDE has to be p_t_crit");

  // save funct_params in class variable
  this->myfunct_params_.resize(5);
  this->myfunct_params_[0] = funct_params[0].second;
  this->myfunct_params_[1] = funct_params[1].second;
  this->myfunct_params_[2] = funct_params[2].second;
  this->myfunct_params_[3] = funct_params[3].second;
  this->myfunct_params_[4] = funct_params[4].second;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
void POROMULTIPHASESCATRA::TumorGrowthLawHeaviside<dim>::CheckOrder(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants)
{
  // safety check for correct ordering of variables and constants
  // they should have been added in exactly the same way in fluidporo_multiphase_singlereaction, but
  // order might be different if we do not use exactly three fluid phases
  if (variables[1].first != "p2") dserror("wrong order in variable vector, p2 not at position 2");
  if (variables[4].first != "S2") dserror("wrong order in variable vector, S2 not at position 5");
  if (variables[6].first != "porosity")
    dserror("wrong order in variable vector, porosity not at position 7");
  if (variables[7].first != "phi1")
    dserror("wrong order in variable vector, phi1 (oxygen mass fraction) not at position 8");
  if (variables[8].first != "phi2")
    dserror("wrong order in variable vector, phi2 (necrotic mass fraction) not at position 9");

  // order is correct
  this->order_checked_ = true;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
double POROMULTIPHASESCATRA::TumorGrowthLawHeaviside<dim>::Evaluate(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants, const size_t component)
{
  // Check order (only once since it does not change)
  if (not this->order_checked_) CheckOrder(variables, constants);

  // read function params
  const double gamma_T_growth = this->myfunct_params_[0];
  const double w_nl_crit = this->myfunct_params_[1];
  const double w_nl_env = this->myfunct_params_[2];
  const double lambda = this->myfunct_params_[3];
  const double p_t_crit = this->myfunct_params_[4];

  // read variables and constants (order is crucial)
  const double p2 = variables[1].second;
  const double S2 = variables[4].second;
  const double porosity = variables[6].second;
  const double oxy_mass_frac = variables[7].second;
  const double necr_frac = variables[8].second;

  // evaluate heaviside
  const double heaviside_oxy((oxy_mass_frac - w_nl_crit) > 0. ? 1. : 0.);
  const double heaviside_pres((p_t_crit - p2) > 0. ? 1. : 0.);
  const double macaulay =
      gamma_T_growth * (oxy_mass_frac - w_nl_crit) / (w_nl_env - w_nl_crit) * heaviside_oxy;

  // evaluate function
  const double functval = (macaulay * heaviside_pres) * (1 - necr_frac) * porosity * S2 -
                          lambda * porosity * oxy_mass_frac * S2;
  return functval;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
std::vector<double> POROMULTIPHASESCATRA::TumorGrowthLawHeaviside<dim>::EvaluateDerivative(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants, const size_t component)
{
  // create derivative vector (should have size of variables)
  std::vector<double> deriv(variables.size(), 0.0);

  // read function params
  const double gamma_T_growth = this->myfunct_params_[0];
  const double w_nl_crit = this->myfunct_params_[1];
  const double w_nl_env = this->myfunct_params_[2];
  const double lambda = this->myfunct_params_[3];
  const double p_t_crit = this->myfunct_params_[4];

  // read variables and constants
  const double p2 = variables[1].second;
  const double S2 = variables[4].second;
  const double porosity = variables[6].second;
  const double oxy_mass_frac = variables[7].second;
  const double necr_frac = variables[8].second;

  // evaluate heaviside
  const double heaviside_oxy((oxy_mass_frac - w_nl_crit) > 0. ? 1. : 0.);
  const double heaviside_pres((p_t_crit - p2) > 0. ? 1. : 0.);
  const double macaulay =
      gamma_T_growth * (oxy_mass_frac - w_nl_crit) / (w_nl_env - w_nl_crit) * heaviside_oxy;

  // set saturation derivs w.r.t. S2
  const double saturationderiv =
      (macaulay * heaviside_pres) * (1 - necr_frac) * porosity - lambda * necr_frac * porosity;
  deriv[4] = saturationderiv;

  // set porosity derivs
  const double porosityderiv =
      (macaulay * heaviside_pres) * (1 - necr_frac) * S2 - lambda * necr_frac * S2;
  deriv[6] = porosityderiv;

  // set scalar derivs w.r.t. phi1
  const double oxygenderiv = gamma_T_growth * heaviside_oxy * heaviside_pres * 1.0 /
                             (w_nl_env - w_nl_crit) * (1 - necr_frac) * porosity * S2;
  deriv[7] = oxygenderiv;

  // set scalar derivs w.r.t. phi2
  const double necroticderiv =
      -(macaulay * heaviside_pres) * porosity * S2 - lambda * porosity * S2;
  deriv[8] = necroticderiv;

  return deriv;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
POROMULTIPHASESCATRA::NecrosisLawHeaviside<dim>::NecrosisLawHeaviside(
    std::vector<std::pair<std::string, double>> funct_params)
    : PoroMultiPhaseScaTraFunction<dim>()
{
  // Check size
  if (funct_params.size() != 5)
  {
    dserror(
        "Wrong size of funct_params for NECROSIS_LAW_HEAVISIDE, it should have exactly\n"
        "5 funct_params (in this order) gamma_t_necr, w_nl_crit, w_nl_env, delta_a_t and p_t_crit");
  }

  // Check correct naming and order of funct_params
  if (funct_params[0].first != "gamma_t_necr")
    dserror("First parameter for NECROSIS_LAW_HEAVISIDE has to be gamma_t_necr");

  if (funct_params[1].first != "w_nl_crit")
    dserror("Second parameter for NECROSIS_LAW_HEAVISIDE has to be w_nl_crit");

  if (funct_params[2].first != "w_nl_env")
    dserror("Third parameter for NECROSIS_LAW_HEAVISIDE has to be w_nl_env");

  if (funct_params[3].first != "delta_a_t")
    dserror("Fourth parameter for NECROSIS_LAW_HEAVISIDE has to be delta_a_t");

  if (funct_params[4].first != "p_t_crit")
    dserror("Fifth parameter for NECROSIS_LAW_HEAVISIDE has to be p_t_crit");

  // save funct_params in class variable
  this->myfunct_params_.resize(5);
  this->myfunct_params_[0] = funct_params[0].second;
  this->myfunct_params_[1] = funct_params[1].second;
  this->myfunct_params_[2] = funct_params[2].second;
  this->myfunct_params_[3] = funct_params[3].second;
  this->myfunct_params_[4] = funct_params[4].second;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
void POROMULTIPHASESCATRA::NecrosisLawHeaviside<dim>::CheckOrder(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants)
{
  // safety check for correct ordering of variables and constants
  // they should have been added in exactly the same way in scatra_ele_calc_multiporo_reaction, but
  // order might be different if we do not use exactly three fluid phases
  if (constants[1].first != "p2") dserror("wrong order in variable vector, p2 not at position 2");
  if (constants[4].first != "S2") dserror("wrong order in variable vector, S2 not at position 5");
  if (constants[6].first != "porosity")
    dserror("wrong order in variable vector, porosity not at position 7");
  if (variables[0].first != "phi1")
    dserror("wrong order in variable vector, phi1 (oxygen mass fraction) not at position 1");
  if (variables[1].first != "phi2")
    dserror("wrong order in variable vector, phi1 (necrotic mass fraction) not at position 2");

  // order is correct
  this->order_checked_ = true;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
double POROMULTIPHASESCATRA::NecrosisLawHeaviside<dim>::Evaluate(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants, const size_t component)
{
  // Check order (only once since it does not change)
  if (not this->order_checked_) CheckOrder(variables, constants);

  // read function params
  const double gamma_t_necr = this->myfunct_params_[0];
  const double w_nl_crit = this->myfunct_params_[1];
  const double w_nl_env = this->myfunct_params_[2];
  const double delta_a_t = this->myfunct_params_[3];
  const double p_t_crit = this->myfunct_params_[4];

  // read variables and constants (order is crucial)
  const double p2 = constants[1].second;
  const double S2 = constants[4].second;
  const double porosity = constants[6].second;
  const double oxy_mass_frac = variables[0].second;
  const double necr_frac = variables[1].second;

  // evaluate heaviside
  const double heaviside_oxy((-(oxy_mass_frac - w_nl_crit)) > 0. ? 1. : 0.);
  const double macaulay =
      -gamma_t_necr * (oxy_mass_frac - w_nl_crit) / (w_nl_env - w_nl_crit) * heaviside_oxy;
  const double heaviside_pres((p2 - p_t_crit) > 0. ? 1. : 0.);

  // evaluate the function
  const double functval =
      (1.0 - necr_frac) * S2 * porosity * (macaulay + delta_a_t * heaviside_pres);

  return functval;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
std::vector<double> POROMULTIPHASESCATRA::NecrosisLawHeaviside<dim>::EvaluateDerivative(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants, const size_t component)
{
  // create derivative vector (should have size of variables)
  std::vector<double> deriv(variables.size(), 0.0);

  // read function params
  const double gamma_t_necr = this->myfunct_params_[0];
  const double w_nl_crit = this->myfunct_params_[1];
  const double w_nl_env = this->myfunct_params_[2];
  const double delta_a_t = this->myfunct_params_[3];
  const double p_t_crit = this->myfunct_params_[4];

  if (variables[0].first == "phi1")  // maindiag-derivative
  {
    // read variables and constants (order is crucial)
    const double p2 = constants[1].second;
    const double S2 = constants[4].second;
    const double porosity = constants[6].second;
    const double oxy_mass_frac = variables[0].second;
    const double necr_frac = variables[1].second;

    // evaluate heaviside
    const double heaviside_oxy((-(oxy_mass_frac - w_nl_crit)) > 0. ? 1. : 0.);
    const double macaulay =
        -gamma_t_necr * (oxy_mass_frac - w_nl_crit) / (w_nl_env - w_nl_crit) * heaviside_oxy;
    const double heaviside_pres((p2 - p_t_crit) > 0. ? 1. : 0.);

    // derivative w.r.t. oxygen mass fraction
    const double oxy_deriv = (1.0 - necr_frac) * porosity * S2 *
                             (-gamma_t_necr * heaviside_oxy * 1.0 / (w_nl_env - w_nl_crit));
    deriv[0] = oxy_deriv;

    // derivative w.r.t. necrotic cell mass fraction
    const double necro_deriv = porosity * S2 * (macaulay + delta_a_t * heaviside_pres) * (-1.0);
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
    const double heaviside_oxy((-(oxy_mass_frac - w_nl_crit)) > 0. ? 1. : 0.);
    const double macaulay =
        -gamma_t_necr * (oxy_mass_frac - w_nl_crit) / (w_nl_env - w_nl_crit) * heaviside_oxy;
    const double heaviside_pres((p2 - p_t_crit) > 0. ? 1. : 0.);

    // derivative w.r.t. tumor cell saturation S2
    const double tc_deriv = (1.0 - necr_frac) * porosity * (macaulay + delta_a_t * heaviside_pres);
    deriv[4] = tc_deriv;

    // derivative w.r.t. porosity
    const double poro_deriv = (1.0 - necr_frac) * S2 * (macaulay + delta_a_t * heaviside_pres);
    deriv[6] = poro_deriv;

    // Note: no pressure derivative, only coupling is with heaviside --> derivative zero
  }
  else
    dserror("Something went wrong in derivative evaluation of NECROSIS_LAW_HEAVISIDE");

  return deriv;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
POROMULTIPHASESCATRA::OxygenConsumptionLawHeaviside<dim>::OxygenConsumptionLawHeaviside(
    std::vector<std::pair<std::string, double>> funct_params)
    : PoroMultiPhaseScaTraFunction<dim>()
{
  // Check size
  if (funct_params.size() != 5)
  {
    dserror(
        "Wrong size of funct_params for NECROSIS_LAW_HEAVISIDE, it should have exactly\n"
        "5 funct_params (in this order) gamma_nl_growth, gamma_0_nl, w_nl_crit, w_nl_env and "
        "p_t_crit");
  }

  // Check correct naming and order of funct_params
  if (funct_params[0].first != "gamma_nl_growth")
    dserror("First parameter for OXYGEN_CONSUMPTION_LAW_HEAVISIDE has to be gamma_nl_growth");

  if (funct_params[1].first != "gamma_0_nl")
    dserror("Second parameter for OXYGEN_CONSUMPTION_LAW_HEAVISIDE has to be gamma_0_nl");

  if (funct_params[2].first != "w_nl_crit")
    dserror("Third parameter for OXYGEN_CONSUMPTION_LAW_HEAVISIDE has to be w_nl_crit");

  if (funct_params[3].first != "w_nl_env")
    dserror("Fourth parameter for OXYGEN_CONSUMPTION_LAW_HEAVISIDE has to be w_nl_env");

  if (funct_params[4].first != "p_t_crit")
    dserror("Fifth parameter for OXYGEN_CONSUMPTION_LAW_HEAVISIDE has to be p_t_crit");

  // save funct_params in class variable
  this->myfunct_params_.resize(5);
  this->myfunct_params_[0] = funct_params[0].second;
  this->myfunct_params_[1] = funct_params[1].second;
  this->myfunct_params_[2] = funct_params[2].second;
  this->myfunct_params_[3] = funct_params[3].second;
  this->myfunct_params_[4] = funct_params[4].second;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
void POROMULTIPHASESCATRA::OxygenConsumptionLawHeaviside<dim>::CheckOrder(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants)
{
  // safety check for correct ordering of variables and constants
  // they should have been added in exactly the same way in scatra_ele_calc_multiporo_reaction, but
  // order might be different if we do not use exactly three fluid phases
  if (constants[1].first != "p2") dserror("wrong order in variable vector, p2 not at position 2");
  if (constants[4].first != "S2") dserror("wrong order in variable vector, S2 not at position 5");
  if (constants[6].first != "porosity")
    dserror("wrong order in variable vector, porosity not at position 7");
  if (variables[0].first != "phi1")
    dserror("wrong order in variable vector, phi1 (oxygen mass fraction) not at position 1");
  if (variables[1].first != "phi2")
    dserror("wrong order in variable vector, phi1 (necrotic mass fraction) not at position 2");

  // order is correct
  this->order_checked_ = true;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
double POROMULTIPHASESCATRA::OxygenConsumptionLawHeaviside<dim>::Evaluate(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants, const size_t component)
{
  // Check order (only once since it does not change)
  if (not this->order_checked_) CheckOrder(variables, constants);

  // read function params
  const double gamma_nl_growth = this->myfunct_params_[0];
  const double gamma_0_nl = this->myfunct_params_[1];
  const double w_nl_crit = this->myfunct_params_[2];
  const double w_nl_env = this->myfunct_params_[3];
  const double p_t_crit = this->myfunct_params_[4];

  // read variables and constants (order is crucial)
  const double S2 = constants[4].second;
  const double porosity = constants[6].second;
  const double p2 = constants[1].second;
  const double oxy_mass_frac = variables[0].second;
  const double necr_frac = variables[1].second;

  // evaluate heaviside
  const double heaviside_oxy((oxy_mass_frac - w_nl_crit) > 0. ? 1. : 0.);
  const double macaulay =
      gamma_nl_growth * (oxy_mass_frac - w_nl_crit) / (w_nl_env - w_nl_crit) * heaviside_oxy;
  const double heaviside_pres((p_t_crit - p2) > 0. ? 1. : 0.);

  // evaluate the function
  const double functval =
      (1.0 - necr_frac) * S2 * porosity *
      (macaulay * heaviside_pres + gamma_0_nl * sin(M_PI / 2.0 * oxy_mass_frac / w_nl_env));

  return functval;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
std::vector<double> POROMULTIPHASESCATRA::OxygenConsumptionLawHeaviside<dim>::EvaluateDerivative(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants, const size_t component)
{
  // create derivative vector (should have size of variables)
  std::vector<double> deriv(variables.size(), 0.0);

  // read function params
  const double gamma_nl_growth = this->myfunct_params_[0];
  const double gamma_0_nl = this->myfunct_params_[1];
  const double w_nl_crit = this->myfunct_params_[2];
  const double w_nl_env = this->myfunct_params_[3];
  const double p_t_crit = this->myfunct_params_[4];

  if (variables[0].first == "phi1")  // maindiag-derivative
  {
    // read variables and constants (order is crucial)
    const double S2 = constants[4].second;
    const double porosity = constants[6].second;
    const double p2 = constants[1].second;
    const double oxy_mass_frac = variables[0].second;
    const double necr_frac = variables[1].second;

    // evaluate heaviside
    const double heaviside_oxy((oxy_mass_frac - w_nl_crit) > 0. ? 1. : 0.);
    const double macaulay =
        gamma_nl_growth * (oxy_mass_frac - w_nl_crit) / (w_nl_env - w_nl_crit) * heaviside_oxy;
    const double heaviside_pres((p_t_crit - p2) > 0. ? 1. : 0.);

    // derivative w.r.t. oxygen mass fraction
    const double oxy_deriv =
        (1.0 - necr_frac) * S2 * porosity *
        (gamma_nl_growth * heaviside_oxy * heaviside_pres * 1.0 / (w_nl_env - w_nl_crit) +
            gamma_0_nl * M_PI / 2.0 / w_nl_env * cos(M_PI / 2.0 * oxy_mass_frac / w_nl_env));
    deriv[0] = oxy_deriv;

    // derivative w.r.t. necrotic cell mass fraction
    const double necro_deriv =
        S2 * porosity * (-1.0) *
        (macaulay * heaviside_pres + gamma_0_nl * sin(M_PI / 2.0 * oxy_mass_frac / w_nl_env));
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
    const double heaviside_oxy((oxy_mass_frac - w_nl_crit) > 0. ? 1. : 0.);
    const double macaulay =
        gamma_nl_growth * (oxy_mass_frac - w_nl_crit) / (w_nl_env - w_nl_crit) * heaviside_oxy;
    const double heaviside_pres((p_t_crit - p2) > 0. ? 1. : 0.);

    // derivative w.r.t. tumor cell saturation S2
    const double tc_deriv =
        (1.0 - necr_frac) * porosity *
        (macaulay * heaviside_pres + gamma_0_nl * sin(M_PI / 2.0 * oxy_mass_frac / w_nl_env));
    deriv[4] = tc_deriv;

    // derivative w.r.t. porosity
    const double poro_deriv =
        (1.0 - necr_frac) * S2 *
        (macaulay * heaviside_pres + gamma_0_nl * sin(M_PI / 2.0 * oxy_mass_frac / w_nl_env));
    deriv[6] = poro_deriv;
  }
  else
    dserror("Something went wrong in derivative evaluation of OXYGEN_CONSUMPTION_LAW_HEAVISIDE");

  return deriv;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
POROMULTIPHASESCATRA::TumorGrowthLawHeavisideOxy<dim>::TumorGrowthLawHeavisideOxy(
    std::vector<std::pair<std::string, double>> funct_params)
    : PoroMultiPhaseScaTraFunction<dim>()
{
  // Check size
  if (funct_params.size() != 5)
  {
    dserror(
        "Wrong size of funct_params for TUMOR_GROWTH_LAW_HEAVISIDE_OXY, it should have exactly\n"
        "5 funct_params (in this order) gamma_T_growth, w_nl_crit, w_nl_env, lambda and p_t_crit");
  }

  // Check correct naming and order of funct_params
  if (funct_params[0].first != "gamma_T_growth")
    dserror("First parameter for TUMOR_GROWTH_LAW_HEAVISIDE_OXY has to be gamma_T_growth");

  if (funct_params[1].first != "w_nl_crit")
    dserror("Second parameter for TUMOR_GROWTH_LAW_HEAVISIDE_OXY has to be w_nl_crit");

  if (funct_params[2].first != "w_nl_env")
    dserror("Third parameter for TUMOR_GROWTH_LAW_HEAVISIDE_OXY has to be w_nl_env");

  if (funct_params[3].first != "lambda")
    dserror("Fourth parameter for TUMOR_GROWTH_LAW_HEAVISIDE_OXY has to be lambda");

  if (funct_params[4].first != "p_t_crit")
    dserror("Fifth parameter for TUMOR_GROWTH_LAW_HEAVISIDE_OXY has to be p_t_crit");

  // save funct_params in class variable
  this->myfunct_params_.resize(5);
  this->myfunct_params_[0] = funct_params[0].second;
  this->myfunct_params_[1] = funct_params[1].second;
  this->myfunct_params_[2] = funct_params[2].second;
  this->myfunct_params_[3] = funct_params[3].second;
  this->myfunct_params_[4] = funct_params[4].second;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
void POROMULTIPHASESCATRA::TumorGrowthLawHeavisideOxy<dim>::CheckOrder(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants)
{
  // safety check for correct ordering of variables and constants
  // they should have been added in exactly the same way in scatra_ele_calc_multiporo_reaction, but
  // order might be different if we do not use exactly three fluid phases
  if (constants[1].first != "p2") dserror("wrong order in variable vector, p2 not at position 2");
  if (constants[4].first != "S2") dserror("wrong order in variable vector, S2 not at position 5");
  if (constants[6].first != "porosity")
    dserror("wrong order in variable vector, porosity not at position 7");
  if (variables[0].first != "phi1")
    dserror("wrong order in variable vector, phi1 (oxygen mass fraction) not at position 1");
  if (variables[1].first != "phi2")
    dserror("wrong order in variable vector, phi1 (necrotic mass fraction) not at position 2");

  // order is correct
  this->order_checked_ = true;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
double POROMULTIPHASESCATRA::TumorGrowthLawHeavisideOxy<dim>::Evaluate(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants, const size_t component)
{
  // Check order (only once since it does not change)
  if (not this->order_checked_) CheckOrder(variables, constants);

  // read function params
  const double gamma_T_growth = this->myfunct_params_[0];
  const double w_nl_crit = this->myfunct_params_[1];
  const double w_nl_env = this->myfunct_params_[2];
  const double lambda = this->myfunct_params_[3];
  const double p_t_crit = this->myfunct_params_[4];

  // read variables and constants (order is crucial)
  const double p2 = constants[1].second;
  const double S2 = constants[4].second;
  const double porosity = constants[6].second;
  const double oxy_mass_frac = variables[0].second;
  const double necr_frac = variables[1].second;

  // evaluate heaviside
  const double heaviside_oxy((oxy_mass_frac - w_nl_crit) > 0. ? 1. : 0.);
  const double heaviside_pres((p_t_crit - p2) > 0. ? 1. : 0.);
  const double macaulay =
      gamma_T_growth * (oxy_mass_frac - w_nl_crit) / (w_nl_env - w_nl_crit) * heaviside_oxy;

  // evaluate function
  const double functval = oxy_mass_frac * S2 * porosity *
                          ((macaulay * heaviside_pres) * (1 - necr_frac) - lambda * necr_frac);

  return functval;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
std::vector<double> POROMULTIPHASESCATRA::TumorGrowthLawHeavisideOxy<dim>::EvaluateDerivative(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants, const size_t component)
{
  // create derivative vector (should have size of variables)
  std::vector<double> deriv(variables.size(), 0.0);

  // read function params
  const double gamma_T_growth = this->myfunct_params_[0];
  const double w_nl_crit = this->myfunct_params_[1];
  const double w_nl_env = this->myfunct_params_[2];
  const double lambda = this->myfunct_params_[3];
  const double p_t_crit = this->myfunct_params_[4];

  if (variables[0].first == "phi1")  // maindiag-derivative
  {
    // read variables and constants (order is crucial)
    const double S2 = constants[4].second;
    const double porosity = constants[6].second;
    const double p2 = constants[1].second;
    const double oxy_mass_frac = variables[0].second;
    const double necr_frac = variables[1].second;

    // evaluate heaviside
    const double heaviside_oxy((oxy_mass_frac - w_nl_crit) > 0. ? 1. : 0.);
    const double heaviside_pres((p_t_crit - p2) > 0. ? 1. : 0.);
    const double macaulay =
        gamma_T_growth * (oxy_mass_frac - w_nl_crit) / (w_nl_env - w_nl_crit) * heaviside_oxy;

    // derivative w.r.t. oxygen mass fraction
    const double oxy_deriv =
        (gamma_T_growth * heaviside_oxy * heaviside_pres * 1.0 / (w_nl_env - w_nl_crit) *
            (1 - necr_frac)) *
            oxy_mass_frac * S2 * porosity +
        ((macaulay * heaviside_pres) * (1 - necr_frac) - lambda * necr_frac) * S2 * porosity;
    deriv[0] = oxy_deriv;

    // derivative w.r.t. necrotic cell mass fraction
    const double necro_deriv =
        ((macaulay * heaviside_pres) * (-1.0) - lambda) * oxy_mass_frac * S2 * porosity;
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
    const double heaviside_oxy((oxy_mass_frac - w_nl_crit) > 0. ? 1. : 0.);
    const double heaviside_pres((p_t_crit - p2) > 0. ? 1. : 0.);
    const double macaulay =
        gamma_T_growth * (oxy_mass_frac - w_nl_crit) / (w_nl_env - w_nl_crit) * heaviside_oxy;

    // derivative w.r.t. tumor cell saturation S2
    const double tc_deriv = ((macaulay * heaviside_pres) * (1.0 - necr_frac) - lambda * necr_frac) *
                            oxy_mass_frac * porosity;
    deriv[4] = tc_deriv;

    // derivative w.r.t. porosity
    const double poro_deriv =
        (macaulay * heaviside_pres * (1.0 - necr_frac) - lambda * necr_frac) * oxy_mass_frac * S2;
    deriv[6] = poro_deriv;
  }
  else
    dserror("Something went wrong in derivative evaluation of TUMOR_GROWTH_LAW_HEAVISIDE_OXY");

  return deriv;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
POROMULTIPHASESCATRA::TumorGrowthLawHeavisideNecro<dim>::TumorGrowthLawHeavisideNecro(
    std::vector<std::pair<std::string, double>> funct_params)
    : PoroMultiPhaseScaTraFunction<dim>()
{
  // Check size
  if (funct_params.size() != 5)
  {
    dserror(
        "Wrong size of funct_params for TUMOR_GROWTH_LAW_HEAVISIDE_OXY, it should have exactly\n"
        "5 funct_params (in this order) gamma_T_growth, w_nl_crit, w_nl_env, lambda and p_t_crit");
  }

  // Check correct naming and order of funct_params
  if (funct_params[0].first != "gamma_T_growth")
    dserror("First parameter for TUMOR_GROWTH_LAW_HEAVISIDE_OXY has to be gamma_T_growth");

  if (funct_params[1].first != "w_nl_crit")
    dserror("Second parameter for TUMOR_GROWTH_LAW_HEAVISIDE_OXY has to be w_nl_crit");

  if (funct_params[2].first != "w_nl_env")
    dserror("Third parameter for TUMOR_GROWTH_LAW_HEAVISIDE_OXY has to be w_nl_env");

  if (funct_params[3].first != "lambda")
    dserror("Fourth parameter for TUMOR_GROWTH_LAW_HEAVISIDE_OXY has to be lambda");

  if (funct_params[4].first != "p_t_crit")
    dserror("Fifth parameter for TUMOR_GROWTH_LAW_HEAVISIDE_OXY has to be p_t_crit");

  // save funct_params in class variable
  this->myfunct_params_.resize(5);
  this->myfunct_params_[0] = funct_params[0].second;
  this->myfunct_params_[1] = funct_params[1].second;
  this->myfunct_params_[2] = funct_params[2].second;
  this->myfunct_params_[3] = funct_params[3].second;
  this->myfunct_params_[4] = funct_params[4].second;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
void POROMULTIPHASESCATRA::TumorGrowthLawHeavisideNecro<dim>::CheckOrder(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants)
{
  // safety check for correct ordering of variables and constants
  // they should have been added in exactly the same way in scatra_ele_calc_multiporo_reaction, but
  // order might be different if we do not use exactly three fluid phases
  if (constants[1].first != "p2") dserror("wrong order in variable vector, p2 not at position 2");
  if (constants[4].first != "S2") dserror("wrong order in variable vector, S2 not at position 5");
  if (constants[6].first != "porosity")
    dserror("wrong order in variable vector, porosity not at position 7");
  if (variables[0].first != "phi1")
    dserror("wrong order in variable vector, phi1 (oxygen mass fraction) not at position 1");
  if (variables[1].first != "phi2")
    dserror("wrong order in variable vector, phi1 (necrotic mass fraction) not at position 2");

  // order is correct
  this->order_checked_ = true;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
double POROMULTIPHASESCATRA::TumorGrowthLawHeavisideNecro<dim>::Evaluate(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants, const size_t component)
{
  // Check order (only once since it does not change)
  if (not this->order_checked_) CheckOrder(variables, constants);

  // read function params
  const double gamma_T_growth = this->myfunct_params_[0];
  const double w_nl_crit = this->myfunct_params_[1];
  const double w_nl_env = this->myfunct_params_[2];
  const double lambda = this->myfunct_params_[3];
  const double p_t_crit = this->myfunct_params_[4];

  // read variables and constants (order is crucial)
  const double p2 = constants[1].second;
  const double S2 = constants[4].second;
  const double porosity = constants[6].second;
  const double oxy_mass_frac = variables[0].second;
  const double necr_frac = variables[1].second;

  // evaluate heaviside
  const double heaviside_oxy((oxy_mass_frac - w_nl_crit) > 0. ? 1. : 0.);
  const double heaviside_pres((p_t_crit - p2) > 0. ? 1. : 0.);
  const double macaulay =
      gamma_T_growth * (oxy_mass_frac - w_nl_crit) / (w_nl_env - w_nl_crit) * heaviside_oxy;

  // evaluate function
  const double functval =
      porosity * S2 *
      (((macaulay * heaviside_pres) * (1 - necr_frac) - lambda * necr_frac) * necr_frac +
          lambda * necr_frac);

  return functval;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
std::vector<double> POROMULTIPHASESCATRA::TumorGrowthLawHeavisideNecro<dim>::EvaluateDerivative(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants, const size_t component)
{
  // create derivative vector (should have size of variables)
  std::vector<double> deriv(variables.size(), 0.0);

  // read function params
  const double gamma_T_growth = this->myfunct_params_[0];
  const double w_nl_crit = this->myfunct_params_[1];
  const double w_nl_env = this->myfunct_params_[2];
  const double lambda = this->myfunct_params_[3];
  const double p_t_crit = this->myfunct_params_[4];

  if (variables[0].first == "phi1")  // maindiag-derivative
  {
    // read variables and constants (order is crucial)
    const double p2 = constants[1].second;
    const double S2 = constants[4].second;
    const double porosity = constants[6].second;
    const double oxy_mass_frac = variables[0].second;
    const double necr_frac = variables[1].second;

    // evaluate heaviside
    const double heaviside_oxy((oxy_mass_frac - w_nl_crit) > 0. ? 1. : 0.);
    const double heaviside_pres((p_t_crit - p2) > 0. ? 1. : 0.);
    const double macaulay =
        gamma_T_growth * (oxy_mass_frac - w_nl_crit) / (w_nl_env - w_nl_crit) * heaviside_oxy;

    // derivative w.r.t. oxygen mass fraction
    const double oxy_deriv =
        (gamma_T_growth * heaviside_oxy * heaviside_pres * 1.0 / (w_nl_env - w_nl_crit)) *
        (1 - necr_frac) * necr_frac * porosity * S2;
    deriv[0] = oxy_deriv;

    // derivative w.r.t. necrotic cell mass fraction
    const double necro_deriv =
        ((macaulay * heaviside_pres) * (1 - 2.0 * necr_frac) - 2.0 * lambda * necr_frac + lambda) *
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
    const double heaviside_oxy((oxy_mass_frac - w_nl_crit) > 0. ? 1. : 0.);
    const double heaviside_pres((p_t_crit - p2) > 0. ? 1. : 0.);
    const double macaulay =
        gamma_T_growth * (oxy_mass_frac - w_nl_crit) / (w_nl_env - w_nl_crit) * heaviside_oxy;

    // derivative w.r.t. tumor cell saturation S2
    const double tc_deriv =
        porosity *
        (((macaulay * heaviside_pres) * (1 - necr_frac) - lambda * necr_frac) * necr_frac +
            lambda * necr_frac);
    deriv[4] = tc_deriv;

    // derivative w.r.t. porosity
    const double poro_deriv =
        S2 * (((macaulay * heaviside_pres) * (1 - necr_frac) - lambda * necr_frac) * necr_frac +
                 lambda * necr_frac);
    deriv[6] = poro_deriv;

    // Note: no pressure derivative, only coupling is with heaviside --> derivative zero
  }
  else
    dserror("Something went wrong in derivative evaluation of TUMOR_GROWTH_LAW_HEAVISIDE_NECRO");

  return deriv;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
POROMULTIPHASESCATRA::OxygenTransvascularExchangeLawCont<dim>::OxygenTransvascularExchangeLawCont(
    std::vector<std::pair<std::string, double>> funct_params)
    : PoroMultiPhaseScaTraFunction<dim>()
{
  // Check size
  if (funct_params.size() != 9)
  {
    dserror(
        "Wrong size of funct_params for OXYGEN_TRANSVASCULAR_EXCHANGE_LAW_CONT, it should have "
        "exactly\n"
        "9 funct_params (in this order) n, Pb50, CaO2_max, alpha_bl_eff, gamma*rho*S/V, rho_oxy, "
        "rho_IF, rho_bl, alpha_IF");
  }

  // Check correct naming and order of funct_params
  if (funct_params[0].first != "n")
    dserror("First parameter for OXYGEN_TRANSVASCULAR_EXCHANGE_LAW_CONT has to be n");

  if (funct_params[1].first != "Pb50")
    dserror("Second parameter for OXYGEN_TRANSVASCULAR_EXCHANGE_LAW_CONT has to be Pb50");

  if (funct_params[2].first != "CaO2_max")
    dserror("Third parameter for OXYGEN_TRANSVASCULAR_EXCHANGE_LAW_CONT has to be CaO2_max");

  if (funct_params[3].first != "alpha_bl_eff")
    dserror("Fourth parameter for OXYGEN_TRANSVASCULAR_EXCHANGE_LAW_CONT has to be alpha_bl_eff");

  if (funct_params[4].first != "gamma*rho*S/V")
    dserror("Fifth parameter for OXYGEN_TRANSVASCULAR_EXCHANGE_LAW_CONT has to be gamma*rho*S/V");

  if (funct_params[5].first != "rho_oxy")
    dserror("Sixth parameter for OXYGEN_TRANSVASCULAR_EXCHANGE_LAW_CONT has to be rho_oxy");

  if (funct_params[6].first != "rho_IF")
    dserror("Sixth parameter for OXYGEN_TRANSVASCULAR_EXCHANGE_LAW_CONT has to be rho_IF");

  if (funct_params[7].first != "rho_bl")
    dserror("Sixth parameter for OXYGEN_TRANSVASCULAR_EXCHANGE_LAW_CONT has to be rho_bl");

  if (funct_params[8].first != "alpha_IF")
    dserror("Sixth parameter for OXYGEN_TRANSVASCULAR_EXCHANGE_LAW_CONT has to be alpha_IF");

  // save funct_params in class variable
  this->myfunct_params_.resize(9);
  this->myfunct_params_[0] = funct_params[0].second;
  this->myfunct_params_[1] = funct_params[1].second;
  this->myfunct_params_[2] = funct_params[2].second;
  this->myfunct_params_[3] = funct_params[3].second;
  this->myfunct_params_[4] = funct_params[4].second;
  this->myfunct_params_[5] = funct_params[5].second;
  this->myfunct_params_[6] = funct_params[6].second;
  this->myfunct_params_[7] = funct_params[7].second;
  this->myfunct_params_[8] = funct_params[8].second;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
void POROMULTIPHASESCATRA::OxygenTransvascularExchangeLawCont<dim>::CheckOrder(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants)
{
  // safety check for correct ordering of variables and constants
  // they should have been added in exactly the same way in scatra_ele_calc_multiporo_reaction, but
  // order might be different if we do not use exactly three fluid phases
  if (constants[7].first != "VF1")
    dserror("wrong order in variable vector, porosity not at position 7");
  if (variables[0].first != "phi1")
    dserror("wrong order in variable vector, phi1 (oxygen mass fraction) not at position 1");
  if (variables[1].first != "phi2")
    dserror("wrong order in variable vector, phi2 (necrotic mass fraction) not at position 2");
  if (variables[2].first != "phi3")
    dserror("wrong order in variable vector, phi3 (necrotic mass fraction) not at position 3");
  if (variables[3].first != "phi4")
    dserror("wrong order in variable vector, phi4 (necrotic mass fraction) not at position 4");

  // order is correct
  this->order_checked_ = true;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
double POROMULTIPHASESCATRA::OxygenTransvascularExchangeLawCont<dim>::Evaluate(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants, const size_t component)
{
  // Check order (only once since it does not change)
  if (not this->order_checked_) CheckOrder(variables, constants);

  // read function params
  const double n = this->myfunct_params_[0];
  const double Pb50 = this->myfunct_params_[1];
  const double CaO2_max = this->myfunct_params_[2];
  const double alpha_bl_eff = this->myfunct_params_[3];
  const double gammarhoSV = this->myfunct_params_[4];
  const double rho_oxy = this->myfunct_params_[5];
  const double rho_if = this->myfunct_params_[6];
  const double rho_bl = this->myfunct_params_[7];
  const double alpha_IF = this->myfunct_params_[8];

  const double fac_if = rho_oxy / rho_if * alpha_IF;

  // read variables and constants (order is crucial)
  double VF1 = constants[7].second;
  double oxy_mass_frac_if = variables[0].second;
  double oxy_mass_frac_nv = variables[3].second;

  double Pb = 0.0;
  double CaO2 = oxy_mass_frac_nv * rho_bl / rho_oxy;
  // safety check --> should not be larger than CaO2_max, which already correponds to partial
  // pressures of ~250Pa
  CaO2 = std::max(0.0, std::min(CaO2, 1.0 * CaO2_max));
  POROMULTIPHASESCATRA::UTILS::GetOxyPartialPressureFromConcentration<double>(
      Pb, CaO2, CaO2_max, Pb50, n, alpha_bl_eff);

  // evaluate function
  const double heaviside_oxy((Pb - oxy_mass_frac_if / fac_if) > 0. ? 1. : 0.);
  const double functval = gammarhoSV * heaviside_oxy * (Pb - oxy_mass_frac_if / fac_if) * VF1;

  return functval;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
std::vector<double>
POROMULTIPHASESCATRA::OxygenTransvascularExchangeLawCont<dim>::EvaluateDerivative(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants, const size_t component)
{
  // create derivative vector (should have size of variables)
  std::vector<double> deriv(variables.size(), 0.0);

  // read function params
  const double n = this->myfunct_params_[0];
  const double Pb50 = this->myfunct_params_[1];
  const double CaO2_max = this->myfunct_params_[2];
  const double alpha_bl_eff = this->myfunct_params_[3];
  const double gammarhoSV = this->myfunct_params_[4];
  const double rho_oxy = this->myfunct_params_[5];
  const double rho_if = this->myfunct_params_[6];
  const double rho_bl = this->myfunct_params_[7];
  const double alpha_IF = this->myfunct_params_[8];

  const double fac_if = rho_oxy / rho_if * alpha_IF;

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
    dserror(
        "Something went wrong in derivative evaluation of OXYGEN_TRANSVASCULAR_EXCHANGE_LAW_CONT");

  FAD Pb = 0.0;
  FAD CaO2 = oxy_mass_frac_nv * rho_bl / rho_oxy;
  // safety check --> should not be larger than CaO2_max, which already correponds to partial
  // pressures of ~250Pa
  CaO2 = std::max(0.0, std::min(CaO2, 1.0 * CaO2_max));
  POROMULTIPHASESCATRA::UTILS::GetOxyPartialPressureFromConcentration<FAD>(
      Pb, CaO2, CaO2_max, Pb50, n, alpha_bl_eff);
  const double heaviside_oxy((Pb - oxy_mass_frac_if / fac_if) > 0. ? 1. : 0.);

  if (variables[0].first == "phi1")  // maindiag-derivative
  {
    deriv[0] = gammarhoSV * VF1 * (-1.0 / fac_if) * heaviside_oxy;
    deriv[3] = gammarhoSV * VF1 * (Pb.fastAccessDx(0)) * heaviside_oxy;
  }
  else if (variables[0].first == "p1")  // OD-derivative
  {
    deriv[7] = gammarhoSV * (Pb.val() - oxy_mass_frac_if / fac_if) * heaviside_oxy;
  }
  else
    dserror(
        "Something went wrong in derivative evaluation of OXYGEN_TRANSVASCULAR_EXCHANGE_LAW_CONT");

  return deriv;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
POROMULTIPHASESCATRA::OxygenTransvascularExchangeLawDisc<dim>::OxygenTransvascularExchangeLawDisc(
    std::vector<std::pair<std::string, double>> funct_params)
    : PoroMultiPhaseScaTraFunction<dim>()
{
  // Check size
  if (funct_params.size() != 10)
  {
    dserror(
        "Wrong size of funct_params for OXYGEN_TRANSVASCULAR_EXCHANGE_LAW_DISC, it should have "
        "exactly\n"
        "10 funct_params (in this order) n, Pb50, CaO2_max, alpha_bl_eff, gamma*rho,rho_oxy, "
        "rho_IF, rho_bl, S2_max, alpha_IF");
  }

  // Check correct naming and order of funct_params
  if (funct_params[0].first != "n")
    dserror("First parameter for OXYGEN_TRANSVASCULAR_EXCHANGE_LAW_DISC has to be n");

  if (funct_params[1].first != "Pb50")
    dserror("Second parameter for OXYGEN_TRANSVASCULAR_EXCHANGE_LAW_DISC has to be Pb50");

  if (funct_params[2].first != "CaO2_max")
    dserror("Third parameter for OXYGEN_TRANSVASCULAR_EXCHANGE_LAW_DISC has to be CaO2_max");

  if (funct_params[3].first != "alpha_bl_eff")
    dserror("Fourth parameter for OXYGEN_TRANSVASCULAR_EXCHANGE_LAW_DISC has to be alpha_bl_eff");

  if (funct_params[4].first != "gamma*rho")
    dserror("Fifth parameter for OXYGEN_TRANSVASCULAR_EXCHANGE_LAW_DISC has to be gamma*rho");

  if (funct_params[5].first != "rho_oxy")
    dserror("Sixth parameter for OXYGEN_TRANSVASCULAR_EXCHANGE_LAW_DISC has to be rho_oxy");

  if (funct_params[6].first != "rho_IF")
    dserror("Sixth parameter for OXYGEN_TRANSVASCULAR_EXCHANGE_LAW_DISC has to be rho_IF");

  if (funct_params[7].first != "rho_bl")
    dserror("Sixth parameter for OXYGEN_TRANSVASCULAR_EXCHANGE_LAW_DISC has to be rho_bl");

  if (funct_params[8].first != "S2_max")
    dserror("Sixth parameter for OXYGEN_TRANSVASCULAR_EXCHANGE_LAW_DISC has to be S2_max");

  if (funct_params[9].first != "alpha_IF")
    dserror("Sixth parameter for OXYGEN_TRANSVASCULAR_EXCHANGE_LAW_DISC has to be alpha_IF");

  // save funct_params in class variable
  this->myfunct_params_.resize(10);
  this->myfunct_params_[0] = funct_params[0].second;
  this->myfunct_params_[1] = funct_params[1].second;
  this->myfunct_params_[2] = funct_params[2].second;
  this->myfunct_params_[3] = funct_params[3].second;
  this->myfunct_params_[4] = funct_params[4].second;
  this->myfunct_params_[5] = funct_params[5].second;
  this->myfunct_params_[6] = funct_params[6].second;
  this->myfunct_params_[7] = funct_params[7].second;
  this->myfunct_params_[8] = funct_params[8].second;
  this->myfunct_params_[9] = funct_params[9].second;

  pos_oxy_art_ = -1;
  pos_diam_ = -1;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
void POROMULTIPHASESCATRA::OxygenTransvascularExchangeLawDisc<dim>::CheckOrder(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants)
{
  // safety check for correct ordering of variables and constants
  // they should have been added in exactly the same way in scatra_ele_calc_multiporo_reaction, but
  // order might be different if we do not use exactly three fluid phases
  if (variables[0].first != "phi1")
    dserror("wrong order in variable vector, phi1 (oxygen mass fraction) not at position 1");

  // oxygen in artery
  // we have no neo-vasculature --> at position 2, species 1: oxy in IF, species 2: NTC
  if (variables[2].first == "phi_art1") pos_oxy_art_ = 2;
  // we have no neo-vasculature --> at position 4, species 1: oxy in IF, species 2: NTC,
  // species 3: TAF, species 4: oxy in NV
  if (variables.size() >= 5 && variables[4].first == "phi_art1") pos_oxy_art_ = 4;
  if (pos_oxy_art_ == -1) dserror("cannot find position of oxygen in arteries");

  // fluid variables
  if (constants[1].first != "p2")
    dserror("wrong order in constants vector, p2 (pressure of tumor cells) not at position 2");
  if (constants[4].first != "S2")
    dserror("wrong order in constants vector, S2 (saturation of tumor cells) not at position 5");

  // diameter
  if (constants[8].first == "D") pos_diam_ = 8;
  if (constants.size() >= 11 && constants[10].first == "D") pos_diam_ = 10;
  if (pos_diam_ == -1) dserror("cannot find position of artery diameter");

  // order is correct
  this->order_checked_ = true;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
double POROMULTIPHASESCATRA::OxygenTransvascularExchangeLawDisc<dim>::Evaluate(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants, const size_t component)
{
  // Check order (only once since it does not change)
  if (not this->order_checked_) CheckOrder(variables, constants);

  // read function params
  const double n = this->myfunct_params_[0];
  const double Pb50 = this->myfunct_params_[1];
  const double CaO2_max = this->myfunct_params_[2];
  const double alpha_bl_eff = this->myfunct_params_[3];
  const double gammarho = this->myfunct_params_[4];
  const double rho_oxy = this->myfunct_params_[5];
  const double rho_if = this->myfunct_params_[6];
  const double rho_bl = this->myfunct_params_[7];
  const double S2_max = this->myfunct_params_[8];
  const double alpha_IF = this->myfunct_params_[9];

  const double fac_if = rho_oxy / rho_if * alpha_IF;

  // read variables and constants (order is crucial)
  double oxy_mass_frac_if = variables[0].second;
  double oxy_mass_frac_nv = variables[pos_oxy_art_].second;
  const double D = constants[pos_diam_].second;
  const double S2 = constants[4].second;

  double Pb = 0.0;
  double CaO2 = oxy_mass_frac_nv * rho_bl / rho_oxy;
  // safety check --> should not be larger than CaO2_max, which already correponds to partial
  // pressures of ~250Pa
  CaO2 = std::max(0.0, std::min(CaO2, 1.0 * CaO2_max));
  POROMULTIPHASESCATRA::UTILS::GetOxyPartialPressureFromConcentration<double>(
      Pb, CaO2, CaO2_max, Pb50, n, alpha_bl_eff);

  // evaluate function
  const double heaviside_S2 = ((S2 - S2_max) > 0. ? 0. : 1.);
  const double functval = gammarho * M_PI * D * (Pb - oxy_mass_frac_if / fac_if) * heaviside_S2;

  return functval;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
std::vector<double>
POROMULTIPHASESCATRA::OxygenTransvascularExchangeLawDisc<dim>::EvaluateDerivative(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants, const size_t component)
{
  // create derivative vector (should have size of variables)
  std::vector<double> deriv(variables.size(), 0.0);

  // read function params
  const double n = this->myfunct_params_[0];
  const double Pb50 = this->myfunct_params_[1];
  const double CaO2_max = this->myfunct_params_[2];
  const double alpha_bl_eff = this->myfunct_params_[3];
  const double gammarho = this->myfunct_params_[4];
  const double rho_oxy = this->myfunct_params_[5];
  const double rho_if = this->myfunct_params_[6];
  const double rho_bl = this->myfunct_params_[7];
  const double S2_max = this->myfunct_params_[8];
  const double alpha_IF = this->myfunct_params_[9];

  const double fac_if = rho_oxy / rho_if * alpha_IF;

  // define Fad object for evaluation
  using FAD = Sacado::Fad::DFad<double>;
  FAD oxy_mass_frac_nv = 0.0;
  oxy_mass_frac_nv.diff(0, 1);  // independent variable 0 out of a total of 1

  // read variables and constants (order is crucial)
  oxy_mass_frac_nv.val() = variables[pos_oxy_art_].second;
  const double D = constants[pos_diam_].second;
  const double S2 = constants[4].second;

  FAD Pb = 0.0;
  FAD CaO2 = oxy_mass_frac_nv * rho_bl / rho_oxy;
  // safety check --> should not be larger than CaO2_max, which already correponds to partial
  // pressures of ~250Pa
  CaO2 = std::max(0.0, std::min(CaO2, 1.0 * CaO2_max));
  POROMULTIPHASESCATRA::UTILS::GetOxyPartialPressureFromConcentration<FAD>(
      Pb, CaO2, CaO2_max, Pb50, n, alpha_bl_eff);

  // evaluate function
  const double heaviside_S2 = ((S2 - S2_max) > 0. ? 0. : 1.);

  deriv[0] = gammarho * M_PI * D * (-1.0 / fac_if) * heaviside_S2;
  deriv[pos_oxy_art_] = gammarho * M_PI * D * (Pb.fastAccessDx(0)) * heaviside_S2;

  return deriv;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
POROMULTIPHASESCATRA::LungOxygenExchangeLaw<dim>::LungOxygenExchangeLaw(
    const std::vector<std::pair<std::string, double>>& funct_params)
    : PoroMultiPhaseScaTraFunction<dim>()
{
  // Check size
  if (funct_params.size() != 9)
  {
    dserror(
        "Wrong size of funct_params for LUNG_OXYGEN_EXCHANGE_LAW, it should have "
        "exactly\n"
        "9 funct_params (in this order) rho_oxy, DiffAdVTLC, alpha_oxy, rho_air, rho_bl, "
        "n, P_oB50, NC_Hb, P_atmospheric");
  }

  // Check correct naming and order of funct_params
  if (funct_params[0].first != "rho_oxy")
    dserror("First parameter for LUNG_OXYGEN_EXCHANGE_LAW has to be rho_oxy");

  if (funct_params[1].first != "DiffAdVTLC")
    dserror("Second parameter for LUNG_OXYGEN_EXCHANGE_LAW has to be DiffAdVTLC");

  if (funct_params[2].first != "alpha_oxy")
    dserror("Third parameter for LUNG_OXYGEN_EXCHANGE_LAW has to be alpha_oxy");

  if (funct_params[3].first != "rho_air")
    dserror("Third parameter for LUNG_OXYGEN_EXCHANGE_LAW has to be rho_air");

  if (funct_params[4].first != "rho_bl")
    dserror("Fourth parameter for LUNG_OXYGEN_EXCHANGE_LAW has to be rho_bl");

  if (funct_params[5].first != "n")
    dserror("Fifth parameter for LUNG_OXYGEN_EXCHANGE_LAW has to be n");

  if (funct_params[6].first != "P_oB50")
    dserror("Sixth parameter for LUNG_OXYGEN_EXCHANGE_LAW has to be P_oB50");

  if (funct_params[7].first != "NC_Hb")
    dserror("Seventh parameter for LUNG_OXYGEN_EXCHANGE_LAW has to be NC_Hb");

  if (funct_params[8].first != "P_atmospheric")
  {
    dserror(
        "Eighth parameter for LUNG_OXYGEN_EXCHANGE_LAW has to be P_atmospheric, which should be "
        "1.013 bar");
  }

  // save funct_params in class variable
  this->myfunct_params_.resize(11);
  this->myfunct_params_[0] = funct_params[0].second;
  this->myfunct_params_[1] = funct_params[1].second;
  this->myfunct_params_[2] = funct_params[2].second;
  this->myfunct_params_[3] = funct_params[3].second;
  this->myfunct_params_[4] = funct_params[4].second;
  this->myfunct_params_[5] = funct_params[5].second;
  this->myfunct_params_[6] = funct_params[6].second;
  this->myfunct_params_[7] = funct_params[7].second;
  this->myfunct_params_[8] = funct_params[8].second;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
void POROMULTIPHASESCATRA::LungOxygenExchangeLaw<dim>::CheckOrder(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants)
{
  // safety check for correct ordering of variables and constants
  if (variables[0].first == "phi1")
  {
    if (constants[0].first != "p1")
      dserror("wrong order in constants vector, P1 (Pressure of air) not at position 0");
    if (variables[1].first != "phi2")
      dserror(
          "wrong order in variable vector, phi2 (oxygen mass fraction in blood) not at position 2");
  }
  else if (variables[0].first == "p1")
  {
    if (constants[0].first != "phi1")
      dserror(
          "wrong order in variable vector, phi1 (oxygen mass fraction in air) not at position 1");
    if (constants[1].first != "phi2")
      dserror(
          "wrong order in variable vector, phi2 (oxygen mass fraction in blood) not at position 2");
  }
  else
  {
    dserror("Variable <%s> not supported on position 0. Wrong order in variable vector! ",
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
    const std::vector<std::pair<std::string, double>>& constants, const size_t component)
{
  // Check order of variables and constants vector only once (since it does not change)
  if (not this->order_checked_) CheckOrder(variables, constants);

    // In debug mode, check order of variables and constants vector on every call
#ifdef DEBUG
  CheckOrder(variables, constants);
#endif

  // read function params
  const double rho_oxy = this->myfunct_params_[0];
  const double DiffAdVTLC = this->myfunct_params_[1];
  const double alpha_oxy = this->myfunct_params_[2];
  const double rho_air = this->myfunct_params_[3];
  const double rho_bl = this->myfunct_params_[4];
  const double n = this->myfunct_params_[5];
  const double P_oB50 = this->myfunct_params_[6];
  const double NC_Hb = this->myfunct_params_[7];
  const double P_atmospheric = this->myfunct_params_[8];

  // read variables (order is crucial)
  const double oxy_mass_frac_air = variables[0].second;
  const double oxy_mass_frac_bl = variables[1].second;

  // read constants (order is crucial)
  const double P_air = constants[0].second;

  // partial pressure of oxygen in air
  const double P_oA = oxy_mass_frac_air * (P_air + P_atmospheric) * rho_air / rho_oxy;

  // CoB_total is total concentration of oxygen in blood (physically dissolved and bound to
  // hemoglobin)
  const double CoB_total = oxy_mass_frac_bl * rho_bl / rho_oxy;

  // partial pressure of oxygen in blood
  double P_oB = 0.0;

  // Calculate partial pressure of oxygen in blood
  POROMULTIPHASESCATRA::UTILS::GetOxyPartialPressureFromConcentration<double>(
      P_oB, CoB_total, NC_Hb, P_oB50, n, alpha_oxy);

  // evaluate function
  const double functval = rho_oxy * DiffAdVTLC * alpha_oxy * (P_oA - P_oB);

  return functval;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int dim>
std::vector<double> POROMULTIPHASESCATRA::LungOxygenExchangeLaw<dim>::EvaluateDerivative(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants, const size_t component)
{
// In debug mode, check order of variables and constants vector on every call
#ifdef DEBUG
  CheckOrder(variables, constants);
#endif

  // create derivative vector (should have size of variables)
  std::vector<double> deriv(variables.size(), 0.0);

  // read function params
  const double rho_oxy = this->myfunct_params_[0];
  const double DiffAdVTLC = this->myfunct_params_[1];
  const double alpha_oxy = this->myfunct_params_[2];
  const double rho_air = this->myfunct_params_[3];
  const double rho_bl = this->myfunct_params_[4];
  const double n = this->myfunct_params_[5];
  const double P_oB50 = this->myfunct_params_[6];
  const double NC_Hb = this->myfunct_params_[7];
  const double P_atmospheric = this->myfunct_params_[8];

  // define Fad object for evaluation
  using FAD = Sacado::Fad::DFad<double>;
  FAD oxy_mass_frac_bl = 0.0;
  oxy_mass_frac_bl.diff(0, 1);  // independent variable 0 out of a total of 1

  double oxy_mass_frac_air = 0.0, P_air = 0.0;

  if (variables[0].first == "phi1")  // maindiag-derivative
  {
    // read variables and constants (order is crucial)
    oxy_mass_frac_air = variables[0].second;
    oxy_mass_frac_bl.val() = variables[1].second;
    P_air = constants[0].second;
  }
  else if (variables[0].first == "p1")  // OD-derivative
  {
    // read variables and constants (order is crucial)
    oxy_mass_frac_air = constants[0].second;
    oxy_mass_frac_bl = constants[1].second;
    P_air = variables[0].second;
  }
  else
    dserror("Derivative w.r.t. <%s> not supported in LUNG_OXYGEN_EXCHANGE_LAW.",
        variables[0].first.c_str());

  FAD P_oB = 0.0;
  FAD C_oB_total = oxy_mass_frac_bl * rho_bl / rho_oxy;

  POROMULTIPHASESCATRA::UTILS::GetOxyPartialPressureFromConcentration<FAD>(
      P_oB, C_oB_total, NC_Hb, P_oB50, n, alpha_oxy);

  if (variables[0].first == "phi1")  // maindiag-derivative
  {
    deriv[0] = rho_oxy * DiffAdVTLC * alpha_oxy * ((P_air + P_atmospheric) * rho_air / rho_oxy);
    deriv[1] = rho_oxy * DiffAdVTLC * alpha_oxy * (-1.0) * P_oB.fastAccessDx(0);
  }
  else if (variables[0].first == "p1")  // OD-derivative
  {
    deriv[0] = (rho_oxy * DiffAdVTLC * alpha_oxy) *
               ((oxy_mass_frac_air * rho_air) /
                   rho_oxy);  // derivative wrt P_air (dFunc/dP_oA * dP_oA/P_air)
  }
  else
    dserror("Derivative w.r.t. <%s> not supported in LUNG_OXYGEN_EXCHANGE_LAW.",
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

template Teuchos::RCP<DRT::UTILS::FunctionOfAnything>
POROMULTIPHASESCATRA::TryCreatePoroFunction<1>(
    const std::vector<Teuchos::RCP<DRT::INPUT::LineDefinition>>& function_line_defs);
template Teuchos::RCP<DRT::UTILS::FunctionOfAnything>
POROMULTIPHASESCATRA::TryCreatePoroFunction<2>(
    const std::vector<Teuchos::RCP<DRT::INPUT::LineDefinition>>& function_line_defs);
template Teuchos::RCP<DRT::UTILS::FunctionOfAnything>
POROMULTIPHASESCATRA::TryCreatePoroFunction<3>(
    const std::vector<Teuchos::RCP<DRT::INPUT::LineDefinition>>& function_line_defs);
