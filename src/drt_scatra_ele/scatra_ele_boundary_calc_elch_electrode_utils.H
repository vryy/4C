/*----------------------------------------------------------------------*/
/*! \file

\brief Utility class supporting evaluation of electrode materials

\level 2

\maintainer Christoph Schmidt
*/
/*----------------------------------------------------------------------*/

#ifndef SCATRA_ELE_BOUNDARY_CALC_ELCH_ELECTRODE_UTILS_H
#define SCATRA_ELE_BOUNDARY_CALC_ELCH_ELECTRODE_UTILS_H

namespace DRT
{
  namespace ELEMENTS
  {
    // forward declaration

    // class implementation
    class ScaTraEleBoundaryCalcElchElectrodeUtils
    {
     public:
      /*!
       * \brief calculate core linearizations associated with Butler-Volmer mass flux density
       *
       * \remark This is a static method as it is also called from static method
       *         EvaluateS2ICouplingAtIntegrationPoint.
       *
       * @param[in] kineticmodel     kinetic model of scatra-scatra interface condition
       * @param[in] j0               exchange mass flux density
       * @param[in] frt              factor F/(RT)
       * @param[in] epdderiv         derivative of equilibrium electric potential difference w.r.t.
       *                             concentration at electrode surface
       * @param[in] alphaa           symmetry coefficient of anodic intercalation reaction
       * @param[in] alphac           symmetry coefficient of cathodic intercalation reaction
       * @param[in] resistance       ohmic resistance on the interface
       * @param[in] expterm1         first exponential term of Butler-Volmer equation
       * @param[in] expterm2         second exponential term of Butler-Volmer equation
       * @param[in] kr               charge transfer constant
       * @param[in] faraday          faraday constant
       * @param[in] emasterphiint    state variables on master-side integration points
       * @param[in] eslavephiint     state variables on slave-side integration points
       * @param[in] cmax             saturation value of intercalated lithium concentration from
       *                             electrode material
       * @param[out] dj_dc_slave     linearization of Butler-Volmer mass flux density w.r.t.
       *                             concentration on slave-side
       * @param[out] dj_dc_master    linearization of Butler-Volmer mass flux density w.r.t.
       *                             concentration on master-side
       * @param[out] dj_dpot_slave   linearization of Butler-Volmer mass flux density w.r.t.
       *                             electric potential on slave-side
       * @param[out] dj_dpot_master  linearization of Butler-Volmer mass flux density w.r.t.
       *                             electric potential on master-side
       *
       * \author @civaner
       *
       * \date 09/19
       */
      static void CalculateCoreLinearizations(const int kineticmodel, const double j0,
          const double frt, const double epdderiv, const double alphaa, const double alphac,
          const double resistance, const double expterm1, const double expterm2, const double kr,
          const double faraday, const double emasterphiint, const double eslavephiint,
          const double cmax, double& dj_dc_slave, double& dj_dc_master, double& dj_dpot_slave,
          double& dj_dpot_master);

      /*!
       * \brief calculate modified Butler-Volmer mass flux density via Newton method
       *
       * \remark This is a static method as it is called from static method
       *         EvaluateS2ICouplingAtIntegrationPoint.
       *
       * @param[in] j0               exchange mass flux density
       * @param[in] alphaa           symmetry coefficient of anodic intercalation reaction
       * @param[in] alphac           symmetry coefficient of cathodic intercalation reaction
       * @param[in] frt              factor F/(RT)
       * @param[in] pot_ed           electric potential on the electrode-side (slave-side)
       * @param[in] pot_el           electric potential on the electrolyte-side (master-side)
       * @param[in] epd              equilibrium electric potential difference
       * @param[in] resistance       ohmic resistance on the interface
       * @param[in] itemax           max. number of iterations for implicit Butler-Volmer equation
       * @param[in] convtol          convergence tolerance for implicit Butler-Volmer equation
       * @param[in] faraday          faraday constant
       * @param[out] j               Butler-Volmer mass flux density
       *
       * \author @civaner
       *
       * \date 08/19
       */
      static double CalculateModifiedButlerVolmerMassFluxDensity(const double j0,
          const double alphaa, const double alphac, const double frt, const double pot_ed,
          const double pot_el, const double epd, const double resistance, const double itemax,
          const double convtol, const double faraday);
    };
  }  // namespace ELEMENTS
}  // namespace DRT

#endif
