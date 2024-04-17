/*----------------------------------------------------------------------*/
/*! \file
\brief inada myocard material model

\level 3

*/

/*----------------------------------------------------------------------*
 |  definitions                                              cbert 08/13 |
 *----------------------------------------------------------------------*/
#ifndef FOUR_C_MAT_MYOCARD_INADA_HPP
#define FOUR_C_MAT_MYOCARD_INADA_HPP

/*----------------------------------------------------------------------*
 |  headers                                                  cbert 08/13 |
 *----------------------------------------------------------------------*/
#include "baci_config.hpp"

#include "baci_comm_parobjectfactory.hpp"
#include "baci_linalg_serialdensematrix.hpp"
#include "baci_linalg_serialdensevector.hpp"
#include "baci_mat_material.hpp"
#include "baci_mat_myocard_general.hpp"
#include "baci_mat_myocard_tools.hpp"
#include "baci_mat_par_parameter.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/// Myocard material according to [1]
///
/// This is a reaction-diffusion law of anisotropic, instationary electric conductivity in cardiac
/// muscle tissue
///
/// <h3>References</h3>
/// <ul>
/// <li> [1] S Inada et. al., "One-Dinemsional Mathematical Model of the Atrioventricular Node
/// Including Atrio-Nodal, Nodal, and Nodal-His Cells", Biophysical Journal 97 (2009) 2117-2127
/// </ul>
///
/// \author ljag


/// \date 09/13

class Myocard_Inada : public Myocard_General

{
 public:
  /// construct empty material object
  Myocard_Inada();

  /// construct empty material object
  explicit Myocard_Inada(const double eps0_deriv_myocard, const std::string tissue);

  /// compute reaction coefficient
  double ReaCoeff(const double phi, const double dt) override;

  ///  returns number of internal state variables of the material
  int GetNumberOfInternalStateVariables() const override;

  ///  returns current internal state of the material
  double GetInternalState(const int k) const override;

  ///  set internal state of the material
  void SetInternalState(const int k, const double val) override;

  ///  return number of ionic currents
  int GetNumberOfIonicCurrents() const override;

  ///  return ionic currents
  double GetIonicCurrents(const int k) const override;

  /// time update for this material
  void Update(const double phi, const double dt) override;

 private:
  Myocard_Tools tools_;

  /// Global time
  double VOI_;

  /// perturbation for numerical approximation of the derivative
  double eps0_deriv_;

  /// gating variables Inada
  std::vector<double> s0_;
  std::vector<double> s_;
  std::vector<double> r_;
  std::vector<double> a_;
  std::vector<double> c_;

};  // Myocard_Inada


FOUR_C_NAMESPACE_CLOSE

#endif
