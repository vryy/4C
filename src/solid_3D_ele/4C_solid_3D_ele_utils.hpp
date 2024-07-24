/*! \file

\brief Helpers for solid elements

\level 1
*/

#ifndef FOUR_C_SOLID_3D_ELE_UTILS_HPP
#define FOUR_C_SOLID_3D_ELE_UTILS_HPP

#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_solid_3D_ele_properties.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Solid::UTILS
{
  void nodal_block_information_solid(
      Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np);

  /*!
   * @brief Converts the 2nd Piola-Kirchhoff stress tensor in stress like Voigt notation to the
   * Cauchy stress tensor in stress like Voigt notation
   *
   * @param pk2 (in) : 2nd Piola-Kirchhoff stress tensor in stress like Voigt notation
   * @param defgrd (in) : Deformation gradient
   * @param cauchy (out) : Cauchy stress tensor in stress like Voigt notation
   */
  void pk2_to_cauchy(const Core::LinAlg::Matrix<6, 1>& pk2,
      const Core::LinAlg::Matrix<3, 3>& defgrd, Core::LinAlg::Matrix<6, 1>& cauchy);

  /*!
   * @brief Convert Green Lagrange strain tensor in strain like Voigt notation to Euler-Almansi
   * strain tensor in strain like Voigt notation.
   *
   * @param gl (in) : Green Lagrange strain tensor in strain like Voigt notation
   * @return Core::LinAlg::Matrix<6, 1> : Euler-Almansi strain tensor in strain like Voigt notation
   */
  Core::LinAlg::Matrix<6, 1> green_lagrange_to_euler_almansi(
      const Core::LinAlg::Matrix<6, 1>& gl, const Core::LinAlg::Matrix<3, 3>& defgrd);

  /*!
   * @brief Convert Green Lagrange strain tensor in strain like Voigt notation to Lograithmic strain
   * tensor in strain like Voigt notation.
   *
   * @param gl (in) : Green Lagrange strain tensor in strain like Voigt notation
   * @return Core::LinAlg::Matrix<6, 1> : Logarithmic strain tensor in strain like Voigt notation
   */
  Core::LinAlg::Matrix<6, 1> green_lagrange_to_log_strain(const Core::LinAlg::Matrix<6, 1>& gl);

  namespace read_element
  {
    int read_element_material(const Core::IO::InputParameterContainer& container);

    Inpar::Solid::KinemType read_element_kinematic_type(
        const Core::IO::InputParameterContainer& container);

    Discret::ELEMENTS::ElementTechnology read_element_technology(
        const Core::IO::InputParameterContainer& container);

    Discret::ELEMENTS::PrestressTechnology read_prestress_technology(
        const Core::IO::InputParameterContainer& container);

    Discret::ELEMENTS::SolidElementProperties read_solid_element_properties(
        const Core::IO::InputParameterContainer& container);

  }  // namespace read_element

}  // namespace Solid::UTILS

FOUR_C_NAMESPACE_CLOSE

#endif
