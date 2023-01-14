/*----------------------------------------------------------------------*/
/*! \file

\brief Utility class supporting evaluation of electrode growth materials

\level 2

*/
/*----------------------------------------------------------------------*/

#ifndef SCATRA_ELE_BOUNDARY_CALC_ELCH_ELECTRODE_GROWTH_UTILS_H
#define SCATRA_ELE_BOUNDARY_CALC_ELCH_ELECTRODE_GROWTH_UTILS_H

// forward declaration
namespace DRT
{
  namespace ELEMENTS
  {
    class ScaTraEleParameterStd;
    class ScaTraEleParameterBoundary;
  }  // namespace ELEMENTS
}  // namespace DRT

namespace DRT
{
  namespace ELEMENTS
  {
    // class implementation
    class ScaTraEleBoundaryCalcElchElectrodeGrowthUtils
    {
     public:
      /*!
       * \brief calculate linearizations associated with Butler-Volmer current density w.r.t. the
       * electrochemical primary variables
       *
       * @param[in] i0               exchange current density
       * @param[in] frt              factor F/(RT)
       * @param[in] epdderiv         derivative of equilibrium electric potential difference w.r.t.
       *                             concentration at electrode surface
       * @param[in] resistance       ohmic resistance on the interface
       * @param[in] regfac           regularization factor for scatra-scatra interface growth
       * @param[in] expterm1         first exponential term of Butler-Volmer equation
       * @param[in] expterm2         second exponential term of Butler-Volmer equation
       * @param[in] faraday          faraday constant
       * @param[in] emasterphiint    state variables on master-side integration points
       * @param[in] eslavephiint     state variables on slave-side integration points
       * @param[in] cmax             saturation value of intercalated lithium concentration from
       *                             electrode material
       * @param[in] scatraeleparamsboundary  scatra ele boundary parameter class
       * @param[out] dj_dc_slave     linearization of Butler-Volmer mass flux density w.r.t.
       *                             concentration on slave-side
       * @param[out] dj_dc_master    linearization of Butler-Volmer mass flux density w.r.t.
       *                             concentration on master-side
       * @param[out] dj_dpot_slave   linearization of Butler-Volmer mass flux density w.r.t.
       *                             electric potential on slave-side
       * @param[out] dj_dpot_master  linearization of Butler-Volmer mass flux density w.r.t.
       *                             electric potential on master-side
       *
       * \author @cschmidt
       *
       * \date 02/20
       */
      static void CalculateS2IElchElchLinearizations(const double i0, const double frt,
          const double epdderiv, const double resistance, const double regfac,
          const double expterm1, const double expterm2, const double faraday,
          const double emasterphiint, const double eslavephiint, const double cmax,
          const DRT::ELEMENTS::ScaTraEleParameterBoundary* const scatraeleparamsboundary,
          double& di_dc_slave, double& di_dc_master, double& di_dpot_slave, double& di_dpot_master);

      /*!
       * \brief calculate linearizations associated with Butler-Volmer current density w.r.t. the
       * interface film thickness
       *
       * @param[in] i0            exchange current density
       * @param[in] i             interface current denisty
       * @param[in] frt           factor F/(RT)
       * @param[in] resistivity   resistivity associated to the material of the interface film
       * @param[in] resistance    ohmic resistance of the interface
       * @param[in] regfac        regularization factor for scatra-scatra interface growth
       * @param[in] regfacderiv   derivative of the regularization factor for scatra-scatra
       *                          interface growth w.r.t. the interface film thickness
       * @param[in] expterm1      first exponential term of Butler-Volmer equation
       * @param[in] expterm2      second exponential term of Butler-Volmer equation
       * @param[in] scatraeleparamsboundary   scatra ele boundary parameter class
       * @return  linearization of Butler-Volmer current density w.r.t. interface film thickness
       *
       * \author @cschmidt
       *
       * \date 02/20
       */
      static double CalculateS2IElchGrowthLinearizations(const double i0, const double i,
          const double frt, const double resistivity, const double resistance, const double regfac,
          const double regfacderiv, const double expterm1, const double expterm2,
          const DRT::ELEMENTS::ScaTraEleParameterBoundary* const scatraeleparamsboundary);

      /*!
       * \brief compute Butler-Volmer current density via local Newton-Raphson iteration
       *
       * @param[in] i0           exchange current density
       * @param[in] frt          factor F/(RT)
       * @param[in] pot_ed       electrode-side electric potential
       * @param[in] pot_el       electrolyte-side electric potential
       * @param[in] epd          half-cell open-circuit potential
       * @param[in] resistance   scatra-scatra interface layer resistance
       * @param[in] thickness    scatra-scatra interface layer thickness
       * @param[in] scatraparameterstd       scatra ele std parameter class
       * @param[in] scatraeleparamsboundary  scatra ele boundary parameter class
       * @return Butler-Volmer current density
       *
       * \author Rui Fang
       *
       * \date 01/17
       */
      static double GetButlerVolmerCurrentDensity(const double i0, const double frt,
          const double pot_ed, const double pot_el, const double epd, const double resistance,
          const double thickness,
          const DRT::ELEMENTS::ScaTraEleParameterStd* const scatraparameterstd,
          const DRT::ELEMENTS::ScaTraEleParameterBoundary* const scatraeleparamsboundary);


      /*!
       * \brief compute regularization factor for lithium stripping
       *
       * @param[in] thickness    scatra-scatra interface layer thickness
       * @param[in] eta          electrode-electrolyte overpotential at integration point
       * @param[in] scatraeleparamsboundary   scatra ele boundary parameter class
       * @return  return the regularization factor if regularization is applied
       *
       * \author Rui Fang
       *
       * \date 01/17
       *
       */
      static double GetRegularizationFactor(const double thickness, const double eta,
          const DRT::ELEMENTS::ScaTraEleParameterBoundary* const scatraeleparamsboundary);


      /*!
       * \brief compute derivative of regularization factor for lithium stripping w.r.t. thickness
       * of plated lithium
       *
       * @param[in] thickness     scatra-scatra interface layer thickness
       * @param[in] eta           electrode-electrolyte overpotential at integration point
       * @param[in] scatraeleparamsboundary   scatra ele boundary parameter class
       * @return  return the derivative of the regularization factor w.r.t. thickness of deposited
       *          material if regularization is applied
       *
       * \author Rui Fang
       *
       * \date 01/17
       */
      static double GetRegularizationFactorDerivative(const double thickness, const double eta,
          const DRT::ELEMENTS::ScaTraEleParameterBoundary* const scatraeleparamsboundary);
    };
  }  // namespace ELEMENTS
}  // namespace DRT

#endif
