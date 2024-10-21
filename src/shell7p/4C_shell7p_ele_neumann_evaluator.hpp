// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SHELL7P_ELE_NEUMANN_EVALUATOR_HPP
#define FOUR_C_SHELL7P_ELE_NEUMANN_EVALUATOR_HPP

#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE


namespace Discret::ELEMENTS::Shell
{
  /*!
   * @brief Evaluates a Neumann condition @p condition for the element @p element.
   *
   * @note This function determines the shape of the element at runtime and calls the respective
   * templated version of @p evaluate_neumann. If you already know the Core::FE::CellType
   * of the element at compile-time, you could directly call @evaluate_neumann.
   *
   * @param element (in) : The element where we integrate
   * @param discretization (in) : discretization
   * @param condition (in) : The Neumann condition to be evaluated within the element.
   * @param dof_index_array (in) : The index array of the dofs of the element
   * @param element_force_vector (in/out) : The element force vector for the evaluated Neumann
   * condition
   * @param element_stiffness_matrix (in/out) : The element stiffness matrix for the evaluated
   * Neumann condition
   * @param total_time (in) : The total time for time dependent Neumann conditions
   */
  void evaluate_neumann_by_element(Core::Elements::Element& element,
      const Core::FE::Discretization& discretization, Core::Conditions::Condition& condition,
      const std::vector<int>& dof_index_array,
      Core::LinAlg::SerialDenseVector& element_force_vector,
      Core::LinAlg::SerialDenseMatrix* element_stiffness_matrix, double total_time);

  /*!
   * @brief Evaluates a Neumann condition @p condition for the element @p element with the
   * discretization type known at compile time.
   *
   *
   * @tparam distype discretization type known at compile time
   *
   * @param element (in) : The element where we integrate
   * @param discretization (in) : discretization
   * @param condition (in) : The Neumann condition to be evaluated within the element.
   * @param dof_index_array (in) : The index array of the dofs of the element
   * @param element_force_vector (in/out) : The element force vector for the evaluated Neumann
   * condition
   * @param element_stiffness_matrix (in/out) : The element stiffness matrix for the evaluated
   * Neumann condition
   * @param total_time (in) : The total time for time dependent Neumann conditions
   */
  template <Core::FE::CellType distype>
  void evaluate_neumann(Core::Elements::Element& element,
      const Core::FE::Discretization& discretization, Core::Conditions::Condition& condition,
      const std::vector<int>& dof_index_array,
      Core::LinAlg::SerialDenseVector& element_force_vector,
      Core::LinAlg::SerialDenseMatrix* element_stiffness_matrix, double total_time);

}  // namespace Discret::ELEMENTS::Shell

FOUR_C_NAMESPACE_CLOSE

#endif
