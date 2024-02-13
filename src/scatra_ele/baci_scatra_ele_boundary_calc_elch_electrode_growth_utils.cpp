/*----------------------------------------------------------------------*/
/*! \file

\brief Utility class supporting evaluation of electrode materials

\level 2

*/
/*----------------------------------------------------------------------*/

#include "baci_scatra_ele_boundary_calc_elch_electrode_growth_utils.hpp"

#include "baci_scatra_ele_parameter_boundary.hpp"
#include "baci_scatra_ele_parameter_std.hpp"
#include "baci_utils_exceptions.hpp"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double DRT::ELEMENTS::CalculateGrowthExchangeMassFluxDensity(const double kr, const double alpha_a,
    const double c_el, const int kinetic_model,
    const DRT::Condition::ConditionType& s2i_condition_type)
{
  dsassert(s2i_condition_type == DRT::Condition::S2IKineticsGrowth,
      "This method is called with the wrong condition type. Check the implementation!");

  switch (kinetic_model)
  {
    case INPAR::S2I::growth_kinetics_butlervolmer:
    {
      return kr * std::pow(c_el, alpha_a);
    }
    default:
      dserror("Did not recognize kinetic model of S2IKineticsGrowth condition!");
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::CalculateS2IGrowthElchLinearizations(const double j0, const double frt,
    const double epdderiv, const double eta, const double resistance, const double regfac,
    const double emasterphiint, const double eslavephiint, const double cmax,
    const DRT::ELEMENTS::ScaTraEleParameterBoundary* const scatraeleparamsboundary,
    double& dj_dc_slave, double& dj_dc_master, double& dj_dpot_slave, double& dj_dpot_master)
{
  const double kr = scatraeleparamsboundary->ChargeTransferConstant();
  const double alphaa = scatraeleparamsboundary->AlphaA();
  const double alphac = scatraeleparamsboundary->AlphaC();

  const double expterm1 = std::exp(alphaa * frt * eta);
  const double expterm2 = std::exp(-alphac * frt * eta);

  // define linearization terms
  double dF_dc_slave(0.0), dF_dc_master(0.0), dF_dpot_slave(0.0), dF_dpot_master(0.0);
  double dF_dj_inverse(0.0);

  switch (scatraeleparamsboundary->RegularizationType())
  {
    case INPAR::S2I::RegularizationType::regularization_none:
    {
      const double expterm = expterm1 - expterm2;

      // compute linearizations of Butler-Volmer mass flux density via implicit differentiation,
      // where F = j0*expterm - j = 0
      dF_dc_slave =
          scatraeleparamsboundary->ConditionType() == DRT::Condition::S2IKinetics
              ? kr * std::pow(emasterphiint, alphaa) * std::pow(cmax - eslavephiint, alphaa - 1.0) *
                        std::pow(eslavephiint, alphac - 1.0) *
                        (-alphaa * eslavephiint + alphac * (cmax - eslavephiint)) * expterm -
                    j0 * (alphaa * frt * epdderiv * expterm1 + alphac * frt * epdderiv * expterm2)
              : 0.0;
      dF_dc_master = j0 * alphaa / emasterphiint * expterm;
      dF_dpot_slave = j0 * frt * (alphaa * expterm1 + alphac * expterm2);
      dF_dpot_master = -dF_dpot_slave;
      dF_dj_inverse =
          -1.0 / (j0 * frt * resistance * (alphaa * expterm1 + alphac * expterm2) + 1.0);

      break;
    }
    case INPAR::S2I::RegularizationType::regularization_trigonometrical:
    case INPAR::S2I::RegularizationType::regularization_polynomial:
    {
      const double expterm = regfac * (expterm1 - expterm2);

      // compute linearizations of Butler-Volmer mass flux density via implicit differentiation,
      // where F = j0*expterm - j = 0
      dF_dc_slave = 0.0;
      dF_dc_master = j0 * alphaa / emasterphiint * expterm;
      dF_dpot_slave = j0 * frt * regfac * (alphaa * expterm1 + alphac * expterm2);
      dF_dpot_master = -dF_dpot_slave;
      dF_dj_inverse =
          -1.0 / (j0 * frt * resistance * regfac * (alphaa * expterm1 + alphac * expterm2) + 1.0);

      break;
    }
    case INPAR::S2I::RegularizationType::regularization_hein:
    {
      const double expterm = regfac * expterm1 - expterm2;

      // compute linearizations of Butler-Volmer mass flux density via implicit differentiation,
      // where F = j0*expterm - j = 0
      dF_dc_slave = 0.0;
      dF_dc_master = j0 * alphaa / emasterphiint * expterm;
      dF_dpot_slave = j0 * frt * (alphaa * regfac * expterm1 + alphac * expterm2);
      dF_dpot_master = -dF_dpot_slave;
      dF_dj_inverse =
          -1.0 / (j0 * frt * resistance * (alphaa * regfac * expterm1 + alphac * expterm2) + 1.0);

      break;
    }
    default:
    {
      dserror("Regularization type %i not recognized!",
          static_cast<int>(scatraeleparamsboundary->RegularizationType()));
    }
  }

  dj_dc_slave = -dF_dc_slave * dF_dj_inverse;
  dj_dc_master = -dF_dc_master * dF_dj_inverse;
  dj_dpot_slave = -dF_dpot_slave * dF_dj_inverse;
  dj_dpot_master = -dF_dpot_master * dF_dj_inverse;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double DRT::ELEMENTS::CalculateS2IElchGrowthLinearizations(const double j0, const double j,
    const double frt, const double resistivity, const double resistance, const double regfac,
    const double regfacderiv, const double expterm1, const double expterm2,
    const DRT::ELEMENTS::ScaTraEleParameterBoundary* const scatraeleparamsboundary)
{
  const double alphaa = scatraeleparamsboundary->AlphaA();
  const double alphac = scatraeleparamsboundary->AlphaC();

  double dF_dgrowth(0.0), dF_dj_inverse(0.0);

  switch (scatraeleparamsboundary->RegularizationType())
  {
    case INPAR::S2I::RegularizationType::regularization_none:
    {
      // compute linearization of Butler-Volmer mass flux density w.r.t. scatra-scatra interface
      // layer thickness via implicit differentiation, where F = j0*expterm - j = 0
      dF_dgrowth = -j0 * j * frt * resistivity * (alphaa * expterm1 + alphac * expterm2);
      dF_dj_inverse =
          -1.0 / (j0 * frt * resistance * (alphaa * expterm1 + alphac * expterm2) + 1.0);

      break;
    }
    case INPAR::S2I::RegularizationType::regularization_trigonometrical:
    case INPAR::S2I::RegularizationType::regularization_polynomial:
    {
      // compute linearization of Butler-Volmer mass flux density w.r.t. scatra-scatra interface
      // layer thickness via implicit differentiation, where F = j0*expterm - j = 0
      dF_dgrowth = -j0 * j * frt * regfac * resistivity * (alphaa * expterm1 + alphac * expterm2) +
                   regfacderiv * j0 * (expterm1 - expterm2);
      dF_dj_inverse =
          -1.0 / (j0 * frt * resistance * regfac * (alphaa * expterm1 + alphac * expterm2) + 1.0);

      break;
    }
    case INPAR::S2I::RegularizationType::regularization_hein:
    {
      // compute linearization of Butler-Volmer mass flux density w.r.t. scatra-scatra interface
      // layer thickness via implicit differentiation, where F = j0*expterm - j = 0
      dF_dgrowth = -j0 * j * frt * resistivity * (alphaa * regfac * expterm1 + alphac * expterm2) +
                   j0 * regfacderiv * expterm1;
      dF_dj_inverse =
          -1.0 / (j0 * frt * resistance * (alphaa * regfac * expterm1 + alphac * expterm2) + 1.0);

      break;
    }
    default:
    {
      dserror("Regularization type %i not recognized!",
          static_cast<int>(scatraeleparamsboundary->RegularizationType()));
    }
  }

  return -dF_dgrowth * dF_dj_inverse;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double DRT::ELEMENTS::CalculateGrowthMassFluxDensity(const double j0, const double frt,
    const double pot_ed, const double pot_el, const double epd, const double resistance,
    const double thickness, const double faraday,
    const DRT::ELEMENTS::ScaTraEleParameterStd* const scatraparameterstd,
    const DRT::ELEMENTS::ScaTraEleParameterBoundary* const scatraeleparamsboundary)
{
  // Iterations are conducted over current density i which is scaled down to mass flux
  // density j by j = i / faraday at the end of the function in order to reduce the effect
  // of numerical error introduced into global problem since i is roughly 10^5 times bigger than j

  // initialize Butler-Volmer current density
  double i(0.0);
  const double i0 = j0 * faraday;

  const double alphaa = scatraeleparamsboundary->AlphaA();
  const double alphac = scatraeleparamsboundary->AlphaC();

  // compute Butler-Volmer current density in case of physically reasonable half-cell open-circuit
  // potential
  if (not std::isinf(epd))
  {
    // compute initial guess of Butler-Volmer current density, neglecting overpotential due to
    // scatra-scatra interface layer resistance
    double eta = pot_ed - pot_el - epd;
    double regfac = GetRegularizationFactor(thickness, eta, scatraeleparamsboundary);
    i = (scatraeleparamsboundary->RegularizationType() ==
            INPAR::S2I::RegularizationType::regularization_hein)
            ? i0 * (regfac * std::exp(alphaa * frt * eta) - std::exp(-alphac * frt * eta))
            : i0 * regfac * (std::exp(alphaa * frt * eta) - std::exp(-alphac * frt * eta));

    // initialize Newton-Raphson iteration counter
    unsigned iternum(0);

    // apply Newton-Raphson method to compute Butler-Volmer current density, involving overpotential
    // due to scatra-scatra interface layer resistance
    while (true)
    {
      // increment counter
      ++iternum;

      // compute current Newton-Raphson residual
      eta = pot_ed - pot_el - epd - resistance * i;
      regfac = GetRegularizationFactor(thickness, eta, scatraeleparamsboundary);
      const double expterm1 = std::exp(alphaa * frt * eta);
      const double expterm2 = std::exp(-alphac * frt * eta);
      const double residual = (scatraeleparamsboundary->RegularizationType() ==
                                  INPAR::S2I::RegularizationType::regularization_hein)
                                  ? i0 * (regfac * expterm1 - expterm2) - i
                                  : i0 * regfac * (expterm1 - expterm2) - i;

      // convergence check
      if (std::abs(residual) < scatraparameterstd->IntLayerGrowthConvTol())
        break;
      else if (iternum == scatraparameterstd->IntLayerGrowthIteMax())
        dserror(
            "Local Newton-Raphson iteration for Butler-Volmer current density did not converge!");

      // compute linearization of current Newton-Raphson residual w.r.t. Butler-Volmer current
      // density
      const double linearization =
          (scatraeleparamsboundary->RegularizationType() ==
              INPAR::S2I::RegularizationType::regularization_hein)
              ? -i0 * resistance * frt * (regfac * alphaa * expterm1 + alphac * expterm2) - 1.0
              : -i0 * resistance * frt * regfac * (alphaa * expterm1 + alphac * expterm2) - 1.0;

      // update Butler-Volmer current density
      i -= residual / linearization;
    }

    // enforce plating condition, i.e., consider initial lithium plating only in case of negative
    // overpotential
    if (std::abs(regfac) < 1.0e-16 and eta >= 0.0) i = 0.0;
  }

  return i / faraday;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double DRT::ELEMENTS::GetRegularizationFactor(const double thickness, const double eta,
    const DRT::ELEMENTS::ScaTraEleParameterBoundary* const scatraeleparamsboundary)
{
  // initialize regularization factor
  double regfac(1.0);

  // get the S2I condition type
  const DRT::Condition::ConditionType conditiontype = scatraeleparamsboundary->ConditionType();

  // get the S2I coupling growth regularization type
  const INPAR::S2I::RegularizationType regularizationtype =
      scatraeleparamsboundary->RegularizationType();

  // actually compute regularization factor if lithium stripping is relevant
  if (conditiontype == DRT::Condition::S2IKineticsGrowth and
      (eta > 0.0 or regularizationtype == INPAR::S2I::RegularizationType::regularization_hein))
  {
    // get the S2I coupling growth regularization parameter
    const double regularizationparameter = scatraeleparamsboundary->RegularizationParameter();
    if (regularizationparameter < 0.0)
      dserror("Regularization parameter for lithium stripping must not be negative!");

    // evaluate dependent on the regularization type
    switch (regularizationtype)
    {
      case INPAR::S2I::RegularizationType::regularization_polynomial:
      // polynomial regularization, cf. Hein, Latz, Electrochimica Acta 201 (2016) 354-365
      case INPAR::S2I::RegularizationType::regularization_hein:
      {
        // use regularization parameter if specified, otherwise take default value according to
        // reference
        const double thickness0 = regularizationparameter > 0.0 ? regularizationparameter : 4.8e-7;

        // compute regularization factor
        regfac = thickness <= 0.0
                     ? 0.0
                     : std::pow(thickness, 4) / (std::pow(thickness, 4) + std::pow(thickness0, 4));

        break;
      }
      // trigonometrical regularization involving (co)sine half-wave
      case INPAR::S2I::RegularizationType::regularization_trigonometrical:
      {
        // use regularization parameter if specified, otherwise take lithium atom diameter as
        // default value
        const double thickness_regend =
            regularizationparameter > 0.0 ? regularizationparameter : 2.9e-7;

        // compute regularization factor
        if (thickness <= 0.0)
          regfac = 0.0;
        else if (thickness < thickness_regend)
          regfac = 0.5 * std::cos(thickness / thickness_regend * M_PI - M_PI) + 0.5;

        break;
      }
      // non-regularized Heaviside function
      case INPAR::S2I::RegularizationType::regularization_none:
      {
        if (thickness <= 0.0) regfac = 0.0;

        break;
      }
      // safety check
      default:
        dserror("Invalid type of regularization: %i for lithium stripping!",
            static_cast<int>(regularizationtype));
    }
  }

  return regfac;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double DRT::ELEMENTS::GetRegularizationFactorDerivative(const double thickness, const double eta,
    const DRT::ELEMENTS::ScaTraEleParameterBoundary* const scatraeleparamsboundary)
{
  // initialize derivative of regularization factor
  double regfacderiv(0.0);

  // get the S2I condition type
  const DRT::Condition::ConditionType conditiontype = scatraeleparamsboundary->ConditionType();

  // get the S2I coupling growth regularization type
  const INPAR::S2I::RegularizationType regularizationtype =
      scatraeleparamsboundary->RegularizationType();

  // actually compute derivative of regularization factor if lithium stripping is relevant
  if (conditiontype == DRT::Condition::S2IKineticsGrowth and thickness > 0.0 and
      (eta > 0.0 or regularizationtype == INPAR::S2I::RegularizationType::regularization_hein))
  {
    // get the S2I coupling growth regularization parameter
    const double regularizationparameter = scatraeleparamsboundary->RegularizationParameter();

    // evaluate dependent on the regularization type
    switch (regularizationtype)
    {
      case INPAR::S2I::RegularizationType::regularization_polynomial:
      // polynomial regularization, cf. Hein, Latz, Electrochimica Acta 201 (2016) 354-365
      case INPAR::S2I::RegularizationType::regularization_hein:
      {
        // use regularization parameter if specified, otherwise take default value according to
        // reference
        const double thickness0 = regularizationparameter > 0.0 ? regularizationparameter : 4.8e-7;

        // compute derivative of regularization factor
        const double thickness0_pow4 = std::pow(thickness0, 4);
        regfacderiv = 4. * std::pow(thickness, 3) * thickness0_pow4 /
                      std::pow(thickness0_pow4 + std::pow(thickness, 4), 2);

        break;
      }
      // trigonometrical regularization involving (co)sine half-wave
      case INPAR::S2I::RegularizationType::regularization_trigonometrical:
      {
        // use regularization parameter if specified, otherwise take lithium atom diameter as
        // default value
        const double thickness_regend =
            regularizationparameter > 0.0 ? regularizationparameter : 2.9e-7;

        // compute derivative of regularization factor
        const double thickness_regend_inverse = 1.0 / thickness_regend;
        if (thickness < thickness_regend)
          regfacderiv = 0.5 * std::sin(thickness * thickness_regend_inverse * M_PI) * M_PI *
                        thickness_regend_inverse;

        break;
      }
      case INPAR::S2I::RegularizationType::regularization_none:
      {
        // do nothing and retain derivative as initialized, since non-regularized Heaviside function
        // cannot be properly differentiated
        break;
      }
      // safety check
      default:
      {
        dserror("Invalid type of regularization: %i for lithium stripping!",
            static_cast<int>(regularizationtype));
      }
    }
  }

  return regfacderiv;
}

BACI_NAMESPACE_CLOSE
