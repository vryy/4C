/*----------------------------------------------------------------------*/
/*! \file

\brief Utility class supporting evaluation of electrode growth materials

\level 2

\maintainer Christoph Schmidt
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
       * \brief compute Butler-Volmer current density via Newton-Raphson iteration
       *
       * @param[in] i0         exchange current density
       * @param[in] alphaa     symmetry coefficient of anodic intercalation reaction
       * @param[in] alphac     symmetry coefficient of cathodic intercalation reaction
       * @param[in] frt        factor F/RT
       * @param[in] pot_ed     electrode-side electric potential
       * @param[in] pot_el     electrolyte-side electric potential
       * @param[in] epd        half-cell open-circuit potential
       * @param[in] resistance scatra-scatra interface layer resistance
       * @param[in] thickness  scatra-scatra interface layer thickness
       * @param[in] scatraparameterstd       scatra ele std parameter class
       * @param[in] scatraeleparamsboundary  scatra ele boundary parameter class
       * @return Butler-Volmer current density
       *
       * \author Rui Fang
       *
       * \date 01/17
       */
      static double GetButlerVolmerCurrentDensity(const double i0, const double alphaa,
          const double alphac, const double frt, const double pot_ed, const double pot_el,
          const double epd, const double resistance, const double thickness,
          DRT::ELEMENTS::ScaTraEleParameterStd* const scatraparameterstd,
          DRT::ELEMENTS::ScaTraEleParameterBoundary* const scatraeleparamsboundary);


      /*!
       * \brief compute regularization factor for lithium stripping
       *
       * @param[in] thickness          scatra-scatra interface layer thickness
       * @param[in] eta                electrode-electrolyte overpotential at integration point
       * @param[in] scatraeleparamsboundary  scatra ele boundary parameter class
       * @return  return the regularization factor if regularization is applied
       *
       * \author Rui Fang
       *
       * \date 01/17
       *
       */
      static double GetRegularizationFactor(const double thickness, const double eta,
          DRT::ELEMENTS::ScaTraEleParameterBoundary* const scatraeleparamsboundary);


      /*!
       * \brief compute derivative of regularization factor for lithium stripping w.r.t. thickness
       * of plated lithium
       *
       * @param[in] thickness          scatra-scatra interface layer thickness
       * @param[in] eta                electrode-electrolyte overpotential at integration point
       * @param[in] scatraeleparamsboundary  scatra ele boundary parameter class
       * @return  return the derivative of the regularization factor w.r.t. thickness of deposited
       *          material if regularization is applied
       *
       * \author Rui Fang
       *
       * \date 01/17
       */
      static double GetRegularizationFactorDerivative(const double thickness, const double eta,
          DRT::ELEMENTS::ScaTraEleParameterBoundary* const scatraeleparamsboundary);
    };
  }  // namespace ELEMENTS
}  // namespace DRT

#endif
