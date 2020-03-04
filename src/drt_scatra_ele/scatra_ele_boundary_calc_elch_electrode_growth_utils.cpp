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
double DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowthUtils::GetButlerVolmerCurrentDensity(
    const double i0, const double alphaa, const double alphac, const double frt,
    const double pot_ed, const double pot_el, const double epd, const double resistance,
    const double thickness, const DRT::ELEMENTS::ScaTraEleParameterStd* const scatraparameterstd,
    const DRT::ELEMENTS::ScaTraEleParameterBoundary* const scatraeleparamsboundary)
{
  // initialize Butler-Volmer current density
  double i(0.0);

  // compute Butler-Volmer current density in case of physically reasonable half-cell open-circuit
  // potential
  if (not std::isinf(epd))
  {
    // compute initial guess of Butler-Volmer current density, neglecting overpotential due to
    // scatra-scatra interface layer resistance
    double eta = pot_ed - pot_el - epd;
    double regfac = GetRegularizationFactor(thickness, eta, scatraeleparamsboundary);
    i = i0 * regfac * (exp(alphaa * frt * eta) - exp(-alphac * frt * eta));

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
      const double residual = i0 * regfac * (expterm1 - expterm2) - i;

      // convergence check
      if (std::abs(residual) < scatraparameterstd->IntLayerGrowthConvTol())
        break;
      else if (iternum == scatraparameterstd->IntLayerGrowthIteMax())
        dserror(
            "Local Newton-Raphson iteration for Butler-Volmer current density did not converge!");

      // compute linearization of current Newton-Raphson residual w.r.t. Butler-Volmer current
      // density
      const double linearization =
          -i0 * resistance * frt * regfac * (alphaa * expterm1 + alphac * expterm2) - 1.0;

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

  // actually compute regularization factor if lithium stripping is relevant
  if (conditiontype == DRT::Condition::S2ICouplingGrowth and eta > 0.0)
  {
    // get the S2I coupling growth regularization parameter
    const double regularizationparameter = scatraeleboundary->RegularizationParameter();
    if (regularizationparameter < 0.0)
      dserror("Regularization parameter for lithium stripping must not be negative!");

    // get the S2I coupling growth regularization type
    const INPAR::S2I::RegularizationType regularizationtype =
        scatraeleboundary->RegularizationType();

    // evaluate dependent on the regularization type
    switch (regularizationtype)
    {
      // polynomial regularization, cf. Hein, Latz, Electrochimica Acta 201 (2016) 354-365
      case INPAR::S2I::RegularizationType::regularization_polynomial:
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

  // actually compute derivative of regularization factor if lithium stripping is relevant
  if (conditiontype == DRT::Condition::S2ICouplingGrowth and thickness > 0.0 and eta > 0.0)
  {
    // get the S2I coupling growth regularization parameter
    const double regularizationparameter = scatraeleboundary->RegularizationParameter();
    if (regularizationparameter < 0.0)
      dserror("Regularization parameter for lithium stripping must not be negative!");

    // get the S2I coupling growth regularization type
    const INPAR::S2I::RegularizationType regularizationtype =
        scatraeleboundary->RegularizationType();

    // evaluate dependent on the regularization type
    switch (regularizationtype)
    {
      // polynomial regularization, cf. Hein, Latz, Electrochimica Acta 201 (2016) 354-365
      case INPAR::S2I::RegularizationType::regularization_polynomial:
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
