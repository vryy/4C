/*----------------------------------------------------------------------*/
/*! \file

\brief Utility class supporting evaluation of electrode materials

\level 2

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_SCATRA_ELE_BOUNDARY_CALC_ELCH_ELECTRODE_UTILS_H
#define BACI_SCATRA_ELE_BOUNDARY_CALC_ELCH_ELECTRODE_UTILS_H

#include "baci_lib_condition.H"

namespace DRT::ELEMENTS
{
  /*!
   * \brief calculate core linearizations of Butler-Volmer mass flux density w.r.t. elch dofs
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
   * @param[in] eta              overpotential in Butler-Volmer equation
   * @param[out] dj_dc_slave     linearization of Butler-Volmer mass flux density w.r.t.
   *                             concentration on slave-side
   * @param[out] dj_dc_master    linearization of Butler-Volmer mass flux density w.r.t.
   *                             concentration on master-side
   * @param[out] dj_dpot_slave   linearization of Butler-Volmer mass flux density w.r.t.
   *                             electric potential on slave-side
   * @param[out] dj_dpot_master  linearization of Butler-Volmer mass flux density w.r.t.
   *                             electric potential on master-side
   */
  void CalculateButlerVolmerElchLinearizations(int kineticmodel, double j0, double frt,
      double epdderiv, double alphaa, double alphac, double resistance, double expterm1,
      double expterm2, double kr, double faraday, double emasterphiint, double eslavephiint,
      double cmax, double eta, double& dj_dc_slave, double& dj_dc_master, double& dj_dpot_slave,
      double& dj_dpot_master);

  /*!
   * \brief calculate core linearizations of Butler-Volmer mass flux density w.r.t. temperature
   * dofs
   *
   * @param[in] alphaa           symmetry coefficient of anodic intercalation reaction
   * @param[in] alphac           symmetry coefficient of cathodic intercalation reaction
   * @param[in] depddT           equilibrium electric potential difference at electrode surface
   *                             w.r.t. temperature
   * @param[in] eta              electrode-electrolyte overpotential
   * @param[in] etempint         average temperature (master and slave side) at interface
   * @param[in] faraday          Faraday constant
   * @param[in] frt              Faraday/(gasconstant * temperature)
   * @param[in] gasconstant      gasconstant
   * @param[in] j0               exchange mass flux density
   * @param[out] dj_dT_slave     linearization of Butler-Volmer mass flux density w.r.t.
   *                             temperature on slave-side
   */
  void CalculateButlerVolmerTempLinearizations(double alphaa, double alphac, double depddT,
      double eta, double etempint, double faraday, double frt, double gasconstant, double j0,
      double& dj_dT_slave);

  /*!
   *
   * @param[in] kineticmodel  kinetic model of scatra-scatra interface condition
   * @param [in] alphaa       symmetry coefficient of anodic intercalation reaction
   * @param [in] alphac       symmetry coefficient of cathodic intercalation reaction
   * @param [in] frt          Faraday/(gasconstant * temperature)
   * @param [in] j0           exchange mass flux density
   * @param [in] eta          electrode-electrolyte overpotential
   * @param [in] depd_ddetF   derivative of equilibrium potential w.r.t. determinant of the
   *                          deformation gradient at the current Gauss point
   * @param [out] dj_dsqrtdetg linearization of Butler-Volmer mass flux density w.r.t. square
   *                           root of the determinant of the metric tensor
   * @param [out] dj_ddetF     linearization of Butler-Volmer mass flux density w.r.t.
   *                           determinant of the deformation gradient at the current Gauss
   *                           point
   */
  void CalculateButlerVolmerDispLinearizations(int kineticmodel, double alphaa, double alphac,
      double frt, double j0, double eta, double depd_ddetF, double& dj_dsqrtdetg, double& dj_ddetF);

  /*!
   * @brief Calculate the exchange mass flux density
   *
   * @param[in] kr                  charge transfer constant
   * @param[in] alpha_a             symmetry coefficient of anodic intercalation reaction
   * @param[in] alpha_c             symmetry coefficient of cathodic intercalation reaction
   * @param[in] c_max               saturation value of intercalated lithium concentration from
   *                                electrode material
   * @param[in] c_ed                electrode-side concentration
   * @param[in] c_el                electrolyte-side concentration
   * @param[in] kinetic_model       kinetic model of scatra-scatra interface condition
   * @param[in] s2i_condition_type  scatra-scatra interface condition type
   * @return exchange mass flux density
   */
  double CalculateButlerVolmerExchangeMassFluxDensity(double kr, double alpha_a, double alpha_c,
      double c_max, double c_ed, double c_el, int kinetic_model,
      const DRT::Condition::ConditionType& s2i_condition_type);

  /*!
   * \brief calculate modified Butler-Volmer mass flux density via Newton method
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
   * @return                     Butler-Volmer mass flux density
   */
  double CalculateModifiedButlerVolmerMassFluxDensity(double j0, double alphaa, double alphac,
      double frt, double pot_ed, double pot_el, double epd, double resistance, double itemax,
      double convtol, double faraday);

  //! Return, if kinetic model uses linearized Butler-Volmer equation
  bool IsButlerVolmerLinearized(int kineticmodel);

  //! Return, if kinetic model uses reduced prefactor in Butler-Volmer equation
  bool IsReducedButlerVolmer(int kineticmodel);

}  // namespace DRT::ELEMENTS

#endif
