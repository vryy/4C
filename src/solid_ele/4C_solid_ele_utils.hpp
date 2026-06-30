// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SOLID_ELE_UTILS_HPP
#define FOUR_C_SOLID_ELE_UTILS_HPP

#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_io_input_spec.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_symmetric_tensor.hpp"
#include "4C_linalg_symmetric_tensor_eigen.hpp"
#include "4C_linalg_tensor_generators.hpp"
#include "4C_solid_ele_calc_lib.hpp"
#include "4C_solid_ele_properties.hpp"
#include "4C_structure_new_input.hpp"

#include <cmath>
#include <ranges>

FOUR_C_NAMESPACE_OPEN

namespace Solid::Utils
{
  void nodal_block_information_solid(Core::Elements::Element* dwele, int& numdf, int& dimns);

  /*!
   * @brief Converts the 2nd Piola-Kirchhoff stress tensor to the
   * Cauchy stress tensor
   *
   * @param pk2 (in) : 2nd Piola-Kirchhoff stress tensor
   * @param defgrd (in) : Deformation gradient
   * @return Cauchy stress tensor
   */
  template <std::size_t dim, Core::FE::CellType celltype>
  Core::LinAlg::SymmetricTensor<double, dim, dim> pk2_to_cauchy(
      const Discret::Elements::ElementProperties<celltype>& element_properties,
      const Core::LinAlg::SymmetricTensor<double, dim, dim>& pk2,
      const Core::LinAlg::Tensor<double, dim, dim>& defgrd)
  {
    if constexpr (dim == 2)
    {
      FOUR_C_ASSERT_ALWAYS(
          element_properties.plane_assumption == Discret::Elements::PlaneAssumption::plane_strain,
          "Cauchy stress output for 2D continua is only available for plane strain elements!");
    }
    return Core::LinAlg::assume_symmetry(defgrd * pk2 * Core::LinAlg::transpose(defgrd)) /
           Core::LinAlg::det(defgrd);
  }

  /*!
   * @brief Convert Green Lagrange strain tensor to Euler-Almansi
   * strain tensor.
   *
   * @param gl (in) : Green Lagrange strain tensor
   * @return Core::LinAlg::Matrix<6, 1> : Euler-Almansi strain tensor
   */
  template <std::size_t dim, Core::FE::CellType celltype>
  Core::LinAlg::SymmetricTensor<double, dim, dim> green_lagrange_to_euler_almansi(
      const Discret::Elements::ElementProperties<celltype>& element_properties,
      const Core::LinAlg::SymmetricTensor<double, dim, dim>& gl,
      const Core::LinAlg::Tensor<double, dim, dim>& defgrd)
  {
    if constexpr (dim == 2)
    {
      FOUR_C_ASSERT_ALWAYS(
          element_properties.plane_assumption == Discret::Elements::PlaneAssumption::plane_strain,
          "Euler Almansi strain output for 2D continua is only available for plane strain "
          "elements!");
    }
    Core::LinAlg::Tensor<double, dim, dim> invdefgrd = Core::LinAlg::inv(defgrd);

    return Core::LinAlg::assume_symmetry(Core::LinAlg::transpose(invdefgrd) * gl * invdefgrd);
  }

  /*!
   * @brief Convert Green Lagrange strain tensor in strain to Lograithmic strain
   * tensor.
   *
   * @param gl (in) : Green Lagrange strain tensor
   * @return Core::LinAlg::Matrix<6, 1> : Logarithmic strain tensor
   */
  template <std::size_t dim, Core::FE::CellType celltype>
  Core::LinAlg::SymmetricTensor<double, dim, dim> green_lagrange_to_log_strain(
      const Discret::Elements::ElementProperties<celltype>& element_properties,
      const Core::LinAlg::SymmetricTensor<double, dim, dim>& gl)
  {
    if constexpr (dim == 2)
    {
      FOUR_C_ASSERT_ALWAYS(
          element_properties.plane_assumption == Discret::Elements::PlaneAssumption::plane_strain,
          "Logarithmic strain output for 2D continua is only available for plane strain "
          "elements!");
    }
    auto [eigenvalues, eigenvectors] = Core::LinAlg::eig(gl);

    // compute principal logarithmic strains
    std::ranges::for_each(
        eigenvalues, [](double& value) { value = std::log(std::sqrt(2 * value + 1.0)); });

    const auto eig = Core::LinAlg::TensorGenerators::diagonal(eigenvalues);
    return Core::LinAlg::assume_symmetry(
        eigenvectors * eig * Core::LinAlg::transpose(eigenvectors));
  }

  namespace ReadElement
  {
    int read_element_material(const Core::IO::InputParameterContainer& container);

    template <unsigned dim>
    Discret::Elements::SolidElementProperties<dim> read_solid_element_properties(
        const Core::IO::InputParameterContainer& container);

  }  // namespace ReadElement

}  // namespace Solid::Utils

FOUR_C_NAMESPACE_CLOSE

#endif
