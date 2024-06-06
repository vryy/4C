/*! \file

\brief Helpers for solid elements

\level 1
*/

#ifndef FOUR_C_SOLID_3D_ELE_UTILS_HPP
#define FOUR_C_SOLID_3D_ELE_UTILS_HPP

#include "4C_config.hpp"

#include "4C_discretization_fem_general_element.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_solid_3D_ele_properties.hpp"

FOUR_C_NAMESPACE_OPEN

namespace STR::UTILS
{
  void NodalBlockInformationSolid(
      Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np);

  /*!
   * @brief Converts the 2nd Piola-Kirchhoff stress tensor in stress like Voigt notation to the
   * Cauchy stress tensor in stress like Voigt notation
   *
   * @param pk2 (in) : 2nd Piola-Kirchhoff stress tensor in stress like Voigt notation
   * @param defgrd (in) : Deformation gradient
   * @param cauchy (out) : Cauchy stress tensor in stress like Voigt notation
   */
  void Pk2ToCauchy(const Core::LinAlg::Matrix<6, 1>& pk2, const Core::LinAlg::Matrix<3, 3>& defgrd,
      Core::LinAlg::Matrix<6, 1>& cauchy);

  /*!
   * @brief Convert Green Lagrange strain tensor in strain like Voigt notation to Euler-Almansi
   * strain tensor in strain like Voigt notation.
   *
   * @param gl (in) : Green Lagrange strain tensor in strain like Voigt notation
   * @return Core::LinAlg::Matrix<6, 1> : Euler-Almansi strain tensor in strain like Voigt notation
   */
  Core::LinAlg::Matrix<6, 1> GreenLagrangeToEulerAlmansi(
      const Core::LinAlg::Matrix<6, 1>& gl, const Core::LinAlg::Matrix<3, 3>& defgrd);

  /*!
   * @brief Convert Green Lagrange strain tensor in strain like Voigt notation to Lograithmic strain
   * tensor in strain like Voigt notation.
   *
   * @param gl (in) : Green Lagrange strain tensor in strain like Voigt notation
   * @return Core::LinAlg::Matrix<6, 1> : Logarithmic strain tensor in strain like Voigt notation
   */
  Core::LinAlg::Matrix<6, 1> GreenLagrangeToLogStrain(const Core::LinAlg::Matrix<6, 1>& gl);

  namespace ReadElement
  {
    int ReadElementMaterial(Input::LineDefinition* linedef);

    Inpar::STR::KinemType ReadElementKinematicType(Input::LineDefinition* linedef);

    Discret::ELEMENTS::ElementTechnology ReadElementTechnology(Input::LineDefinition* linedef);

    Discret::ELEMENTS::PrestressTechnology ReadPrestressTechnology(Input::LineDefinition* linedef);

    Discret::ELEMENTS::SolidElementProperties ReadSolidElementProperties(
        Input::LineDefinition* linedef);

  }  // namespace ReadElement

}  // namespace STR::UTILS

FOUR_C_NAMESPACE_CLOSE

#endif
