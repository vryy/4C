/*----------------------------------------------------------------------*/
/*! \file

\brief Utility class supporting evaluation of electrode materials

\level 2

\maintainer Christoph Schmidt
*/
/*----------------------------------------------------------------------*/

#include "scatra_ele_boundary_calc_elch_electrode_utils.H"

#include "../drt_inpar/inpar_s2i.H"
#include "../drt_lib/drt_dserror.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeUtils::CalculateCoreLinearizations(
    const int kineticmodel, const double timefacfac, const double timefacrhsfac, const double j0,
    const double frt, const double epdderiv, const double alphaa, const double alphac,
    const double resistance, const double expterm1, const double expterm2, const double kr,
    const double faraday, const double emasterphiint, const double eslavephiint, const double cmax,
    double& dj_dc_slave, double& dj_dc_master, double& dj_dpot_slave, double& dj_dpot_master)
{
  const double expterm = expterm1 - expterm2;
  // core linearizations associated with Butler-Volmer mass flux density
  switch (kineticmodel)
  {
    case INPAR::S2I::kinetics_butlervolmerreduced:
    {
      dj_dc_slave = timefacfac * j0 * frt * epdderiv * (-alphaa * expterm1 - alphac * expterm2);
      dj_dc_master = 0.0;
      dj_dpot_slave = timefacfac * j0 * (alphaa * frt * expterm1 + alphac * frt * expterm2);
      dj_dpot_master = -dj_dpot_slave;
      break;
    }
    case INPAR::S2I::kinetics_butlervolmer:
    case INPAR::S2I::kinetics_butlervolmerpeltier:
    {
      dj_dc_slave =
          timefacfac * (kr * pow(emasterphiint, alphaa) * pow(cmax - eslavephiint, alphaa - 1.) *
                               pow(eslavephiint, alphac - 1.) *
                               (-alphaa * eslavephiint + alphac * (cmax - eslavephiint)) * expterm +
                           j0 * frt * epdderiv * (-alphaa * expterm1 - alphac * expterm2));
      dj_dc_master = timefacfac * j0 * alphaa / emasterphiint * expterm;
      dj_dpot_slave = timefacfac * j0 * (alphaa * frt * expterm1 + alphac * frt * expterm2);
      dj_dpot_master = -dj_dpot_slave;
      break;
    }
    case INPAR::S2I::kinetics_butlervolmerresistance:
    {
      // core linearizations associated with Butler-Volmer current density according to MA Schmidt
      // 2016 via implicit differentiation where F(x,i) = i - i0 * expterm
      const double dF_di_inverse =
          1. / (1.0 + j0 * faraday * resistance * frt * (alphaa * expterm1 + alphac * expterm2));
      const double dF_dc_slave =
          timefacfac * (-kr * pow(emasterphiint, alphaa) * pow(cmax - eslavephiint, alphaa - 1.) *
                               pow(eslavephiint, alphac - 1.) *
                               (-alphaa * eslavephiint + alphac * (cmax - eslavephiint)) * expterm +
                           j0 * frt * epdderiv * (alphaa * expterm1 + alphac * expterm2));
      const double dF_dc_master = -timefacfac * j0 * alphaa / emasterphiint * expterm;
      const double dF_dpot_slave = -timefacfac * j0 * frt * (alphaa * expterm1 + alphac * expterm2);
      const double dF_dpot_master = -dF_dpot_slave;

      // rule of implicit differentiation
      dj_dc_slave = -dF_dc_slave * dF_di_inverse;
      dj_dc_master = -dF_dc_master * dF_di_inverse;
      dj_dpot_slave = -dF_dpot_slave * dF_di_inverse;
      dj_dpot_master = -dF_dpot_master * dF_di_inverse;
      break;
    }  // case INPAR::S2I::kinetics_butlervolmerresistance
    case INPAR::S2I::kinetics_butlervolmerreducedwithresistance:
    {
      // core linearizations associated with Butler-Volmer current density according to MA Schmidt
      // 2016 via implicit differentiation where F(x,i) = i - i0 * expterm
      const double dF_di_inverse =
          1. / (1.0 + j0 * faraday * resistance * frt * (alphaa * expterm1 + alphac * expterm2));
      const double dF_dc_slave =
          timefacfac * j0 * frt * epdderiv * (alphaa * expterm1 + alphac * expterm2);
      const double dF_dc_master = 0.0;
      const double dF_dpot_slave = -timefacfac * j0 * frt * (alphaa * expterm1 + alphac * expterm2);
      const double dF_dpot_master = -dF_dpot_slave;

      // rule of implicit differentiation
      dj_dc_slave = -dF_dc_slave * dF_di_inverse;
      dj_dc_master = -dF_dc_master * dF_di_inverse;
      dj_dpot_slave = -dF_dpot_slave * dF_di_inverse;
      dj_dpot_master = -dF_dpot_master * dF_di_inverse;
      break;
    }  // case INPAR::S2I::kinetics_butlervolmerreducedwithresistance
    default:
    {
      dserror("Unknown scatra-scatra interface kinetic model: %i", kineticmodel);
      break;
    }
  }  // switch(kineticmodel)
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeUtils::
    CalculateModifiedButlerVolmerMassFluxDensity(const double j0, const double alphaa,
        const double alphac, const double frt, const double pot_ed, const double pot_el,
        const double epd, const double resistance, const double itemax, const double convtol,
        const double faraday)
{
  // Iterations are conducted over current density i which is scaled down to mass flux
  // density j by j = i / faraday at the end of the function in order to reduce the effect
  // of numerical error introduced into global problem since i is roughly 10^5 times bigger than j

  // initialize Butler-Volmer current density
  double i(0.0);

  const double i0 = j0 * faraday;
  // compute Butler-Volmer current density in case of physically reasonable half-cell open-circuit
  // potential
  if (not std::isinf(epd))
  {
    // initialize Newton-Raphson iteration counter
    unsigned iternum(0);

    // apply Newton-Raphson method to compute Butler-Volmer current density, involving overpotential
    // due to scatra-scatra interface layer resistance
    while (true)
    {
      // increment counter
      ++iternum;

      // compute current Newton-Raphson residual
      const double eta = pot_ed - pot_el - epd - resistance * i;
      const double expterm1 = exp(alphaa * frt * eta);
      const double expterm2 = exp(-alphac * frt * eta);
      const double residual = i0 * (expterm1 - expterm2) - i;

      // convergence check
      if (std::abs(residual) < convtol)
        break;

      else if (iternum == itemax)
        dserror(
            "Local Newton-Raphson iteration for Butler-Volmer current density did not converge!");

      // compute linearization of current Newton-Raphson residual w.r.t. Butler-Volmer current
      // density
      const double linearization =
          -i0 * resistance * frt * (alphaa * expterm1 + alphac * expterm2) - 1.0;

      // update Butler-Volmer current density
      i -= residual / linearization;
    }
  }
  // final scaling
  return i / faraday;
}
