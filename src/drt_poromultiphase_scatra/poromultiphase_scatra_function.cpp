/*----------------------------------------------------------------------*/
/*!
 \file poromultiphase_scatra_function.cpp

 \brief Managing and evaluating of (reaction) functions for poromultiphase_scatra
        problems

   \level 3

   \maintainer  Johannes Kremheller
                kremheller@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/

#include "poromultiphase_scatra_function.H"
#include "poromultiphase_scatra_utils.H"
#include "../headers/FAD_utils.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
POROMULTIPHASESCATRA::PoroMultiPhaseScaTraFunction::PoroMultiPhaseScaTraFunction()
    : order_checked_(false)
{
}


// standard growth law for tumor cells <--> IF (with lysis) and pressure dependency:
// +-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
// | TUMOR_GROWTH_LAW_HEAVISIDE | |
// (gamma_T_growth*(phi1-w_nl_crit)/(w_nl_env-w_nl_crit)*heaviside((phi1-w_nl_crit)/(w_nl_env-w_nl_crit))*heaviside(p_t_crit-p2))*(1-phi2)*porosity*S2
// - lambda*phi2*porosity*S2 | | with phi1: mass fraction of oxygen, phi2: mass fraction of necrotic
// tumor cells, S2: volume fraction of tumor cells | | Furthermore, we assume that phase 1: healthy
// cells, phase2: tumor cells, phase3: IF | | | | INPUT DEFINITION: | |
// POROMULTIPHASESCATRA_FUNCTION TUMOR_GROWTH_LAW_HEAVISIDE NUMPARAMS 5 PARAMS gamma_T_growth 9.6e-6
// w_nl_crit 2.0e-6 w_nl_env 4.2e-6 lambda 0.0 p_t_crit 1.0e9                  |
// +-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

POROMULTIPHASESCATRA::TumorGrowthLawHeaviside::TumorGrowthLawHeaviside(
    std::vector<std::pair<std::string, double>> funct_params)
    : PoroMultiPhaseScaTraFunction()
{
  // Check size
  if (funct_params.size() != 5)
    dserror(
        "Wrong size of funct_params for TUMOR_GROWTH_LAW_HEAVISIDE, it should have exactly\n"
        "5 funct_params (in this order) gamma_T_growth, w_nl_crit, w_nl_env, lambda and p_t_crit");

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
  myfunct_params_.resize(5);
  myfunct_params_[0] = funct_params[0].second;
  myfunct_params_[1] = funct_params[1].second;
  myfunct_params_[2] = funct_params[2].second;
  myfunct_params_[3] = funct_params[3].second;
  myfunct_params_[4] = funct_params[4].second;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::TumorGrowthLawHeaviside::CheckOrder(
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
  order_checked_ = true;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double POROMULTIPHASESCATRA::TumorGrowthLawHeaviside::Evaluate(const int index,
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants)
{
  // Check order (only once since it does not change)
  if (not order_checked_) CheckOrder(variables, constants);

  // read function params
  const double gamma_T_growth = myfunct_params_[0];
  const double w_nl_crit = myfunct_params_[1];
  const double w_nl_env = myfunct_params_[2];
  const double lambda = myfunct_params_[3];
  const double p_t_crit = myfunct_params_[4];

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
std::vector<double> POROMULTIPHASESCATRA::TumorGrowthLawHeaviside::EvaluateDerivative(int index,
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants)
{
  // create derivative vector (should have size of variables)
  std::vector<double> deriv(variables.size(), 0.0);

  // read function params
  const double gamma_T_growth = myfunct_params_[0];
  const double w_nl_crit = myfunct_params_[1];
  const double w_nl_env = myfunct_params_[2];
  const double lambda = myfunct_params_[3];
  const double p_t_crit = myfunct_params_[4];

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

// standard necrosis law for tumor growth model:
// +-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
// | NECROSIS_LAW_HEAVISIDE | |
// (1-phi2)*porosity*S2*(-gamma_t_necr*(phi1-w_nl_crit)/(w_nl_env-w_nl_crit)*heaviside(-(phi1-w_nl_crit)/(w_nl_env-w_nl_crit))+
// delta_a_t*heaviside(p2-p_t_crit) )               | | with phi1: mass fraction of oxygen, phi2:
// mass fraction of necrotic tumor cells, S2: volume fraction of tumor cells | | Furthermore, we
// assume that phase 1: healthy cells, phase2: tumor cells, phase3: IF | | | | (possible) INPUT
// DEFINITION with exactly 5 function parameters: | | POROMULTIPHASESCATRA_FUNCTION
// NECROSIS_LAW_HEAVISIDE NUMPARAMS 5 PARAMS gamma_t_necr 9.6e-6 w_nl_crit 2.0e-6 w_nl_env 4.2e-6
// delta_a_t 0.0 p_t_crit 1.0e9                     |
// +-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

POROMULTIPHASESCATRA::NecrosisLawHeaviside::NecrosisLawHeaviside(
    std::vector<std::pair<std::string, double>> funct_params)
    : PoroMultiPhaseScaTraFunction()
{
  // Check size
  if (funct_params.size() != 5)
    dserror(
        "Wrong size of funct_params for NECROSIS_LAW_HEAVISIDE, it should have exactly\n"
        "5 funct_params (in this order) gamma_t_necr, w_nl_crit, w_nl_env, delta_a_t and p_t_crit");

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
  myfunct_params_.resize(5);
  myfunct_params_[0] = funct_params[0].second;
  myfunct_params_[1] = funct_params[1].second;
  myfunct_params_[2] = funct_params[2].second;
  myfunct_params_[3] = funct_params[3].second;
  myfunct_params_[4] = funct_params[4].second;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::NecrosisLawHeaviside::CheckOrder(
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
  order_checked_ = true;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double POROMULTIPHASESCATRA::NecrosisLawHeaviside::Evaluate(const int index,
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants)
{
  // Check order (only once since it does not change)
  if (not order_checked_) CheckOrder(variables, constants);

  // read function params
  const double gamma_t_necr = myfunct_params_[0];
  const double w_nl_crit = myfunct_params_[1];
  const double w_nl_env = myfunct_params_[2];
  const double delta_a_t = myfunct_params_[3];
  const double p_t_crit = myfunct_params_[4];

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
std::vector<double> POROMULTIPHASESCATRA::NecrosisLawHeaviside::EvaluateDerivative(int index,
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants)
{
  // create derivative vector (should have size of variables)
  std::vector<double> deriv(variables.size(), 0.0);

  // read function params
  const double gamma_t_necr = myfunct_params_[0];
  const double w_nl_crit = myfunct_params_[1];
  const double w_nl_env = myfunct_params_[2];
  const double delta_a_t = myfunct_params_[3];
  const double p_t_crit = myfunct_params_[4];

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

// standard oxygen consumption law for tumor growth model:
// +--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
// | OXYGEN_CONSUMPTION_LAW_HEAVISIDE | |
// porosity*(1-phi2)*S2*(gamma_nl_growth*(phi1-w_nl_crit)/(w_nl_env-w_nl_crit)*heaviside((phi1-w_nl_crit)/(w_nl_env-w_nl_crit))*heaviside(p_t_crit-p2)+
// gamma_0_nl*sin(pi/2.0*phi1/w_nl_env)) | | with phi1: mass fraction of oxygen, phi2: mass fraction
// of necrotic tumor cells, S2: volume fraction of tumor cells | | Furthermore, we assume that phase
// 1: healthy cells, phase2: tumor cells, phase3: IF | | | | (possible) INPUT DEFINITION with
// exactly 5 function parameters: | | POROMULTIPHASESCATRA_FUNCTION OXYGEN_CONSUMPTION_LAW_HEAVISIDE
// NUMPARAMS 5 PARAMS gamma_nl_growth 2.4e-7 gamma_0_nl 6e-7 w_nl_crit 2.0e-6 w_nl_env 4.2e-6
// p_t_crit 1.0e9                   |
// +--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

POROMULTIPHASESCATRA::OxygenConsumptionLawHeaviside::OxygenConsumptionLawHeaviside(
    std::vector<std::pair<std::string, double>> funct_params)
    : PoroMultiPhaseScaTraFunction()
{
  // Check size
  if (funct_params.size() != 5)
    dserror(
        "Wrong size of funct_params for NECROSIS_LAW_HEAVISIDE, it should have exactly\n"
        "5 funct_params (in this order) gamma_nl_growth, gamma_0_nl, w_nl_crit, w_nl_env and "
        "p_t_crit");

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
  myfunct_params_.resize(5);
  myfunct_params_[0] = funct_params[0].second;
  myfunct_params_[1] = funct_params[1].second;
  myfunct_params_[2] = funct_params[2].second;
  myfunct_params_[3] = funct_params[3].second;
  myfunct_params_[4] = funct_params[4].second;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::OxygenConsumptionLawHeaviside::CheckOrder(
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
  order_checked_ = true;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double POROMULTIPHASESCATRA::OxygenConsumptionLawHeaviside::Evaluate(const int index,
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants)
{
  // Check order (only once since it does not change)
  if (not order_checked_) CheckOrder(variables, constants);

  // read function params
  const double gamma_nl_growth = myfunct_params_[0];
  const double gamma_0_nl = myfunct_params_[1];
  const double w_nl_crit = myfunct_params_[2];
  const double w_nl_env = myfunct_params_[3];
  const double p_t_crit = myfunct_params_[4];

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
std::vector<double> POROMULTIPHASESCATRA::OxygenConsumptionLawHeaviside::EvaluateDerivative(
    int index, const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants)
{
  // create derivative vector (should have size of variables)
  std::vector<double> deriv(variables.size(), 0.0);

  // read function params
  const double gamma_nl_growth = myfunct_params_[0];
  const double gamma_0_nl = myfunct_params_[1];
  const double w_nl_crit = myfunct_params_[2];
  const double w_nl_env = myfunct_params_[3];
  const double p_t_crit = myfunct_params_[4];

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

// standard growth law for tumor cells <--> IF (with lysis) and pressure dependency as introduced
// into balance of mass of oxygen:
// +-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
// | TUMOR_GROWTH_LAW_HEAVISIDE_OXY | | phi1*S2*porosity*((
// gamma_T_growth*(phi1-w_nl_crit)/(w_nl_env-w_nl_crit)*heaviside((phi1-w_nl_crit)/(w_nl_env-w_nl_crit))*heaviside(p_t_crit-p2))*(1-phi2)-
// lambda*phi2)      | | with phi1: mass fraction of oxygen, phi2: mass fraction of necrotic tumor
// cells, S2: volume fraction of tumor cells | | Furthermore, we assume that phase 1: healthy cells,
// phase2: tumor cells, phase3: IF | | | | (possible) INPUT DEFINITION: | |
// POROMULTIPHASESCATRA_FUNCTION TUMOR_GROWTH_LAW_HEAVISIDE_OXY NUMPARAMS 5 PARAMS
// gamma_T_growth 9.6e-6 w_nl_crit 2.0e-6 w_nl_env 4.2e-6 lambda 0.0 p_t_crit 1.0e9              |
// +-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

POROMULTIPHASESCATRA::TumorGrowthLawHeavisideOxy::TumorGrowthLawHeavisideOxy(
    std::vector<std::pair<std::string, double>> funct_params)
    : PoroMultiPhaseScaTraFunction()
{
  // Check size
  if (funct_params.size() != 5)
    dserror(
        "Wrong size of funct_params for TUMOR_GROWTH_LAW_HEAVISIDE_OXY, it should have exactly\n"
        "5 funct_params (in this order) gamma_T_growth, w_nl_crit, w_nl_env, lambda and p_t_crit");

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
  myfunct_params_.resize(5);
  myfunct_params_[0] = funct_params[0].second;
  myfunct_params_[1] = funct_params[1].second;
  myfunct_params_[2] = funct_params[2].second;
  myfunct_params_[3] = funct_params[3].second;
  myfunct_params_[4] = funct_params[4].second;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::TumorGrowthLawHeavisideOxy::CheckOrder(
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
  order_checked_ = true;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double POROMULTIPHASESCATRA::TumorGrowthLawHeavisideOxy::Evaluate(const int index,
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants)
{
  // Check order (only once since it does not change)
  if (not order_checked_) CheckOrder(variables, constants);

  // read function params
  const double gamma_T_growth = myfunct_params_[0];
  const double w_nl_crit = myfunct_params_[1];
  const double w_nl_env = myfunct_params_[2];
  const double lambda = myfunct_params_[3];
  const double p_t_crit = myfunct_params_[4];

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
std::vector<double> POROMULTIPHASESCATRA::TumorGrowthLawHeavisideOxy::EvaluateDerivative(int index,
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants)
{
  // create derivative vector (should have size of variables)
  std::vector<double> deriv(variables.size(), 0.0);

  // read function params
  const double gamma_T_growth = myfunct_params_[0];
  const double w_nl_crit = myfunct_params_[1];
  const double w_nl_env = myfunct_params_[2];
  const double lambda = myfunct_params_[3];
  const double p_t_crit = myfunct_params_[4];

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

// standard growth law for tumor cells <--> IF (with lysis) and pressure dependency as introduced
// into balance of mass of necrotic cells:
// +---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
// | TUMOR_GROWTH_LAW_HEAVISIDE_NECRO | | porosity*S2*(((
// gamma_T_growth*(phi1-w_nl_crit)/(w_nl_env-w_nl_crit)*heaviside((phi1-w_nl_crit)/(w_nl_env-w_nl_crit))*heaviside(p_t_crit-p2))*(1-phi2)
// - lambda*phi2)*phi2 + lambda*phi2)   | | with phi1: mass fraction of oxygen, phi2: mass fraction
// of necrotic tumor cells | | Furthermore, we assume that phase 1: healthy cells, phase2: tumor
// cells, phase3: IF | | | | (possible) INPUT DEFINITION: | | POROMULTIPHASESCATRA_FUNCTION
// TUMOR_GROWTH_LAW_HEAVISIDE_NECRO NUMPARAMS 5 PARAMS gamma_T_growth 9.6e-6 w_nl_crit 2.0e-6
// w_nl_env 4.2e-6 lambda 0.0 p_t_crit 1.0e9                          |
// +---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

POROMULTIPHASESCATRA::TumorGrowthLawHeavisideNecro::TumorGrowthLawHeavisideNecro(
    std::vector<std::pair<std::string, double>> funct_params)
    : PoroMultiPhaseScaTraFunction()
{
  // Check size
  if (funct_params.size() != 5)
    dserror(
        "Wrong size of funct_params for TUMOR_GROWTH_LAW_HEAVISIDE_OXY, it should have exactly\n"
        "5 funct_params (in this order) gamma_T_growth, w_nl_crit, w_nl_env, lambda and p_t_crit");

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
  myfunct_params_.resize(5);
  myfunct_params_[0] = funct_params[0].second;
  myfunct_params_[1] = funct_params[1].second;
  myfunct_params_[2] = funct_params[2].second;
  myfunct_params_[3] = funct_params[3].second;
  myfunct_params_[4] = funct_params[4].second;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::TumorGrowthLawHeavisideNecro::CheckOrder(
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
  order_checked_ = true;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double POROMULTIPHASESCATRA::TumorGrowthLawHeavisideNecro::Evaluate(const int index,
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants)
{
  // Check order (only once since it does not change)
  if (not order_checked_) CheckOrder(variables, constants);

  // read function params
  const double gamma_T_growth = myfunct_params_[0];
  const double w_nl_crit = myfunct_params_[1];
  const double w_nl_env = myfunct_params_[2];
  const double lambda = myfunct_params_[3];
  const double p_t_crit = myfunct_params_[4];

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
std::vector<double> POROMULTIPHASESCATRA::TumorGrowthLawHeavisideNecro::EvaluateDerivative(
    int index, const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants)
{
  // create derivative vector (should have size of variables)
  std::vector<double> deriv(variables.size(), 0.0);

  // read function params
  const double gamma_T_growth = myfunct_params_[0];
  const double w_nl_crit = myfunct_params_[1];
  const double w_nl_env = myfunct_params_[2];
  const double lambda = myfunct_params_[3];
  const double p_t_crit = myfunct_params_[4];

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

// transvascular exchange of oxygen from neovasculature into interstitial fluid
// +-------------------------------------------------------------------------------------------+
// | OXYGEN_TRANSVASCULAR_EXCHANGE_LAW_CONT                                                    |
// | S/V*gamma*rho*(P_nv-P_if)*heaviside(P_nv-P_if)*VF1                                        |
// | with partial pressures of oxygen in neovasculature P_nv and in IF P_if                    |
// +-------------------------------------------------------------------------------------------+
POROMULTIPHASESCATRA::OxygenTransvascularExchangeLawCont::OxygenTransvascularExchangeLawCont(
    std::vector<std::pair<std::string, double>> funct_params)
    : PoroMultiPhaseScaTraFunction()
{
  // Check size
  if (funct_params.size() != 9)
    dserror(
        "Wrong size of funct_params for OXYGEN_TRANSVASCULAR_EXCHANGE_LAW_CONT, it should have "
        "exactly\n"
        "9 funct_params (in this order) n, Pb50, CaO2_max, alpha_bl_eff, gamma*rho*S/V, rho_oxy, "
        "rho_IF, rho_bl, alpha_IF");

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
  myfunct_params_.resize(9);
  myfunct_params_[0] = funct_params[0].second;
  myfunct_params_[1] = funct_params[1].second;
  myfunct_params_[2] = funct_params[2].second;
  myfunct_params_[3] = funct_params[3].second;
  myfunct_params_[4] = funct_params[4].second;
  myfunct_params_[5] = funct_params[5].second;
  myfunct_params_[6] = funct_params[6].second;
  myfunct_params_[7] = funct_params[7].second;
  myfunct_params_[8] = funct_params[8].second;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::OxygenTransvascularExchangeLawCont::CheckOrder(
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
  order_checked_ = true;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double POROMULTIPHASESCATRA::OxygenTransvascularExchangeLawCont::Evaluate(const int index,
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants)
{
  // Check order (only once since it does not change)
  if (not order_checked_) CheckOrder(variables, constants);

  // read function params
  const double n = myfunct_params_[0];
  const double Pb50 = myfunct_params_[1];
  const double CaO2_max = myfunct_params_[2];
  const double alpha_bl_eff = myfunct_params_[3];
  const double gammarhoSV = myfunct_params_[4];
  const double rho_oxy = myfunct_params_[5];
  const double rho_if = myfunct_params_[6];
  const double rho_bl = myfunct_params_[7];
  const double alpha_IF = myfunct_params_[8];

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
  POROMULTIPHASESCATRA::UTILS::GetOxyPartialPressure<double>(
      Pb, CaO2, CaO2_max, Pb50, n, alpha_bl_eff);

  // evaluate function
  const double heaviside_oxy((Pb - oxy_mass_frac_if / fac_if) > 0. ? 1. : 0.);
  const double functval = gammarhoSV * heaviside_oxy * (Pb - oxy_mass_frac_if / fac_if) * VF1;

  return functval;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> POROMULTIPHASESCATRA::OxygenTransvascularExchangeLawCont::EvaluateDerivative(
    int index, const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants)
{
  // create derivative vector (should have size of variables)
  std::vector<double> deriv(variables.size(), 0.0);

  // read function params
  const double n = myfunct_params_[0];
  const double Pb50 = myfunct_params_[1];
  const double CaO2_max = myfunct_params_[2];
  const double alpha_bl_eff = myfunct_params_[3];
  const double gammarhoSV = myfunct_params_[4];
  const double rho_oxy = myfunct_params_[5];
  const double rho_if = myfunct_params_[6];
  const double rho_bl = myfunct_params_[7];
  const double alpha_IF = myfunct_params_[8];

  const double fac_if = rho_oxy / rho_if * alpha_IF;

  double VF1, oxy_mass_frac_if;
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
  POROMULTIPHASESCATRA::UTILS::GetOxyPartialPressure<FAD>(
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

// transvascular exchange of oxygen from pre-existing vasculature into interstitial fluid
// +-------------------------------------------------------------------------------------------+
// | OXYGEN_TRANSVASCULAR_EXCHANGE_LAW_DISC                                                    |
// | pi*D*gamma*rho*(P_v-P_if)*VF1*heaviside(S2_max-S2)                                        |
// | with partial pressures of oxygen in vasculature P_v and in IF P_if                        |
// +-------------------------------------------------------------------------------------------+
POROMULTIPHASESCATRA::OxygenTransvascularExchangeLawDisc::OxygenTransvascularExchangeLawDisc(
    std::vector<std::pair<std::string, double>> funct_params)
    : PoroMultiPhaseScaTraFunction()
{
  // Check size
  if (funct_params.size() != 10)
    dserror(
        "Wrong size of funct_params for OXYGEN_TRANSVASCULAR_EXCHANGE_LAW_DISC, it should have "
        "exactly\n"
        "10 funct_params (in this order) n, Pb50, CaO2_max, alpha_bl_eff, gamma*rho,rho_oxy, "
        "rho_IF, rho_bl, S2_max, alpha_IF");

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
  myfunct_params_.resize(10);
  myfunct_params_[0] = funct_params[0].second;
  myfunct_params_[1] = funct_params[1].second;
  myfunct_params_[2] = funct_params[2].second;
  myfunct_params_[3] = funct_params[3].second;
  myfunct_params_[4] = funct_params[4].second;
  myfunct_params_[5] = funct_params[5].second;
  myfunct_params_[6] = funct_params[6].second;
  myfunct_params_[7] = funct_params[7].second;
  myfunct_params_[8] = funct_params[8].second;
  myfunct_params_[9] = funct_params[9].second;

  pos_oxy_art_ = -1;
  pos_diam_ = -1;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::OxygenTransvascularExchangeLawDisc::CheckOrder(
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
  order_checked_ = true;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double POROMULTIPHASESCATRA::OxygenTransvascularExchangeLawDisc::Evaluate(const int index,
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants)
{
  // Check order (only once since it does not change)
  if (not order_checked_) CheckOrder(variables, constants);

  // read function params
  const double n = myfunct_params_[0];
  const double Pb50 = myfunct_params_[1];
  const double CaO2_max = myfunct_params_[2];
  const double alpha_bl_eff = myfunct_params_[3];
  const double gammarho = myfunct_params_[4];
  const double rho_oxy = myfunct_params_[5];
  const double rho_if = myfunct_params_[6];
  const double rho_bl = myfunct_params_[7];
  const double S2_max = myfunct_params_[8];
  const double alpha_IF = myfunct_params_[9];

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
  POROMULTIPHASESCATRA::UTILS::GetOxyPartialPressure<double>(
      Pb, CaO2, CaO2_max, Pb50, n, alpha_bl_eff);

  // evaluate function
  const double heaviside_S2 = ((S2 - S2_max) > 0. ? 0. : 1.);
  const double functval = gammarho * M_PI * D * (Pb - oxy_mass_frac_if / fac_if) * heaviside_S2;

  return functval;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> POROMULTIPHASESCATRA::OxygenTransvascularExchangeLawDisc::EvaluateDerivative(
    int index, const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants)
{
  // create derivative vector (should have size of variables)
  std::vector<double> deriv(variables.size(), 0.0);

  // read function params
  const double n = myfunct_params_[0];
  const double Pb50 = myfunct_params_[1];
  const double CaO2_max = myfunct_params_[2];
  const double alpha_bl_eff = myfunct_params_[3];
  const double gammarho = myfunct_params_[4];
  const double rho_oxy = myfunct_params_[5];
  const double rho_if = myfunct_params_[6];
  const double rho_bl = myfunct_params_[7];
  const double S2_max = myfunct_params_[8];
  const double alpha_IF = myfunct_params_[9];

  const double fac_if = rho_oxy / rho_if * alpha_IF;

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
  POROMULTIPHASESCATRA::UTILS::GetOxyPartialPressure<FAD>(
      Pb, CaO2, CaO2_max, Pb50, n, alpha_bl_eff);

  // evaluate function
  const double heaviside_S2 = ((S2 - S2_max) > 0. ? 0. : 1.);

  deriv[0] = gammarho * M_PI * D * (-1.0 / fac_if) * heaviside_S2;
  deriv[pos_oxy_art_] = gammarho * M_PI * D * (Pb.fastAccessDx(0)) * heaviside_S2;

  return deriv;
}
