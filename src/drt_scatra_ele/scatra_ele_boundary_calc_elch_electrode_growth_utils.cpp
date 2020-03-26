/*----------------------------------------------------------------------*/
/*! \file

\brief Utility class supporting evaluation of electrode materials

\level 2

\maintainer Christoph Schmidt
*/
/*----------------------------------------------------------------------*/

#include "scatra_ele_boundary_calc_elch_electrode_growth_utils.H"
#include "scatra_ele_parameter_boundary.H"
#include "scatra_ele_parameter_std.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowthUtils::
    CalculateS2IElchElchLinearizations(const double i0, const double frt, const double epdderiv,
        const double resistance, const double regfac, const double expterm1, const double expterm2,
        const double faraday, const double emasterphiint, const double eslavephiint,
        const double cmax,
        const DRT::ELEMENTS::ScaTraEleParameterBoundary* const scatraeleparamsboundary,
        double& di_dc_slave, double& di_dc_master, double& di_dpot_slave, double& di_dpot_master)
{
  const double kr = scatraeleparamsboundary->Kr();
  const double alphaa = scatraeleparamsboundary->AlphaA();
  const double alphac = scatraeleparamsboundary->AlphaC();

  // define linearization terms
  double dF_dc_slave(0.0), dF_dc_master(0.0), dF_dpot_slave(0.0), dF_dpot_master(0.0);
  double dF_di_inverse(0.0);

  switch (scatraeleparamsboundary->RegularizationType())
  {
    case INPAR::S2I::RegularizationType::regularization_none:
    case INPAR::S2I::RegularizationType::regularization_trigonometrical:
    case INPAR::S2I::RegularizationType::regularization_polynomial:
    {
      const double expterm = regfac * (expterm1 - expterm2);

      // safety check
      if (abs(expterm) > 1.0e5)
        dserror("Overflow of exponential term in Butler-Volmer formulation detected! Value: %lf",
            expterm);

      // compute linearizations of Butler-Volmer current density via implicit differentiation, where
      // F = i0*expterm - i = 0
      // also true for this DRT::Condition::S2ICoupling, as regfac is equal to 1 in this case
      dF_dc_slave =
          scatraeleparamsboundary->ConditionType() == DRT::Condition::S2ICoupling
              ? kr * faraday * pow(emasterphiint, alphaa) * pow(cmax - eslavephiint, alphaa - 1.0) *
                        pow(eslavephiint, alphac - 1.0) *
                        (-alphaa * eslavephiint + alphac * (cmax - eslavephiint)) * expterm -
                    i0 * (alphaa * frt * epdderiv * expterm1 + alphac * frt * epdderiv * expterm2)
              : 0.0;
      dF_dc_master = i0 * alphaa / emasterphiint * expterm;
      dF_dpot_slave = i0 * frt * regfac * (alphaa * expterm1 + alphac * expterm2);
      dF_dpot_master = -dF_dpot_slave;
      dF_di_inverse =
          -1.0 / (i0 * frt * resistance * regfac * (alphaa * expterm1 + alphac * expterm2) + 1.0);

      break;
    }
    case INPAR::S2I::RegularizationType::regularization_hein:
    {
      const double expterm = regfac * expterm1 - expterm2;

      // safety check
      if (abs(expterm) > 1.0e5)
        dserror("Overflow of exponential term in Butler-Volmer formulation detected! Value: %lf",
            expterm);

      // compute linearizations of Butler-Volmer current density via implicit differentiation, where
      // F = i0*expterm - i = 0
      // also true for this DRT::Condition::S2ICoupling, as regfac is equal to 1 in this case
      dF_dc_slave =
          scatraeleparamsboundary->ConditionType() == DRT::Condition::S2ICoupling
              ? kr * faraday * pow(emasterphiint, alphaa) * pow(cmax - eslavephiint, alphaa - 1.0) *
                        pow(eslavephiint, alphac - 1.0) *
                        (-alphaa * eslavephiint + alphac * (cmax - eslavephiint)) * expterm -
                    i0 * (alphaa * frt * epdderiv * expterm1 + alphac * frt * epdderiv * expterm2)
              : 0.0;
      dF_dc_master = i0 * alphaa / emasterphiint * expterm;
      dF_dpot_slave = i0 * frt * (alphaa * regfac * expterm1 + alphac * expterm2);
      dF_dpot_master = -dF_dpot_slave;
      dF_di_inverse =
          -1.0 / (i0 * frt * resistance * (alphaa * regfac * expterm1 + alphac * expterm2) + 1.0);

      break;
    }
    default:
    {
      dserror("Regularization type %i not recognized!",
          static_cast<int>(scatraeleparamsboundary->RegularizationType()));
      break;
    }
  }

  di_dc_slave = -dF_dc_slave * dF_di_inverse;
  di_dc_master = -dF_dc_master * dF_di_inverse;
  di_dpot_slave = -dF_dpot_slave * dF_di_inverse;
  di_dpot_master = -dF_dpot_master * dF_di_inverse;

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double
DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowthUtils::CalculateS2IElchGrowthLinearizations(
    const double i0, const double i, const double frt, const double resistivity,
    const double resistance, const double regfac, const double regfacderiv, const double expterm1,
    const double expterm2,
    const DRT::ELEMENTS::ScaTraEleParameterBoundary* const scatraeleparamsboundary)
{
  const double alphaa = scatraeleparamsboundary->AlphaA();
  const double alphac = scatraeleparamsboundary->AlphaC();

  double dF_dgrowth(0.0), dF_di_inverse(0.0), di_dgrowth(0.0);

  switch (scatraeleparamsboundary->RegularizationType())
  {
    case INPAR::S2I::RegularizationType::regularization_none:
    case INPAR::S2I::RegularizationType::regularization_trigonometrical:
    case INPAR::S2I::RegularizationType::regularization_polynomial:
    {
      // compute linearization of Butler-Volmer current density w.r.t. scatra-scatra interface layer
      // thickness via implicit differentiation, where F = i0*expterm - i = 0
      dF_dgrowth = -i0 * i * frt * regfac * resistivity * (alphaa * expterm1 + alphac * expterm2) +
                   regfacderiv * i0 * (expterm1 - expterm2);
      dF_di_inverse =
          -1.0 / (i0 * frt * resistance * regfac * (alphaa * expterm1 + alphac * expterm2) + 1.0);

      break;
    }
    case INPAR::S2I::RegularizationType::regularization_hein:
    {
      // compute linearization of Butler-Volmer current density w.r.t. scatra-scatra interface layer
      // thickness via implicit differentiation, where F = i0*expterm - i = 0
      dF_dgrowth = -i0 * i * frt * resistivity * (alphaa * regfac * expterm1 + alphac * expterm2) +
                   i0 * regfacderiv * expterm1;
      dF_di_inverse =
          -1.0 / (i0 * frt * resistance * (alphaa * regfac * expterm1 + alphac * expterm2) + 1.0);

      break;
    }
    default:
    {
      dserror("Regularization type %i not recognized!",
          static_cast<int>(scatraeleparamsboundary->RegularizationType()));
      break;
    }
  }

  di_dgrowth = -dF_dgrowth * dF_di_inverse;

  return di_dgrowth;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowthUtils::GetButlerVolmerCurrentDensity(
    const double i0, const double frt, const double pot_ed, const double pot_el, const double epd,
    const double resistance, const double thickness,
    const DRT::ELEMENTS::ScaTraEleParameterStd* const scatraparameterstd,
    const DRT::ELEMENTS::ScaTraEleParameterBoundary* const scatraeleparamsboundary)
{
  // initialize Butler-Volmer current density
  double i(0.0);

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
            ? i0 * (regfac * exp(alphaa * frt * eta) - exp(-alphac * frt * eta))
            : i0 * regfac * (exp(alphaa * frt * eta) - exp(-alphac * frt * eta));

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
      const double expterm1 = exp(alphaa * frt * eta);
      const double expterm2 = exp(-alphac * frt * eta);
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

  return i;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowthUtils::GetRegularizationFactor(
    const double thickness, const double eta,
    const DRT::ELEMENTS::ScaTraEleParameterBoundary* const scatraeleboundary)
{
  // initialize regularization factor
  double regfac(1.0);

  // get the S2I condition type
  const DRT::Condition::ConditionType conditiontype = scatraeleboundary->ConditionType();

  // get the S2I coupling growth regularization type
  const INPAR::S2I::RegularizationType regularizationtype = scatraeleboundary->RegularizationType();

  // actually compute regularization factor if lithium stripping is relevant
  if (conditiontype == DRT::Condition::S2ICouplingGrowth and
      (eta > 0.0 or regularizationtype == INPAR::S2I::RegularizationType::regularization_hein))
  {
    // get the S2I coupling growth regularization parameter
    const double regularizationparameter = scatraeleboundary->RegularizationParameter();
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
        regfac =
            thickness <= 0.0 ? 0.0 : pow(thickness, 4) / (pow(thickness, 4) + pow(thickness0, 4));

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
          regfac = 0.5 * cos(thickness / thickness_regend * M_PI - M_PI) + 0.5;

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
        break;
    }
  }

  return regfac;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double
DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowthUtils::GetRegularizationFactorDerivative(
    const double thickness, const double eta,
    const DRT::ELEMENTS::ScaTraEleParameterBoundary* const scatraeleboundary)
{
  // initialize derivative of regularization factor
  double regfacderiv(0.0);

  // get the S2I condition type
  const DRT::Condition::ConditionType conditiontype = scatraeleboundary->ConditionType();

  // get the S2I coupling growth regularization type
  const INPAR::S2I::RegularizationType regularizationtype = scatraeleboundary->RegularizationType();

  // actually compute derivative of regularization factor if lithium stripping is relevant
  if (conditiontype == DRT::Condition::S2ICouplingGrowth and thickness > 0.0 and
      (eta > 0.0 or regularizationtype == INPAR::S2I::RegularizationType::regularization_hein))
  {
    // get the S2I coupling growth regularization parameter
    const double regularizationparameter = scatraeleboundary->RegularizationParameter();
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

        // compute derivative of regularization factor
        const double thickness0_pow4 = pow(thickness0, 4);
        regfacderiv =
            4. * pow(thickness, 3) * thickness0_pow4 / pow(thickness0_pow4 + pow(thickness, 4), 2);

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
          regfacderiv = 0.5 * sin(thickness * thickness_regend_inverse * M_PI) * M_PI *
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
        break;
      }
    }
  }

  return regfacderiv;
}
