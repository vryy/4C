/*! \file

\brief Evaluation of neumann loads for shell7p element

\level 3
 */

#ifndef FOUR_C_SHELL7P_ELE_NEUMANN_EVALUATOR_HPP
#define FOUR_C_SHELL7P_ELE_NEUMANN_EVALUATOR_HPP

#include "4C_config.hpp"

#include "4C_lib_element.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  class Discretization;
  class Condition;

}  // namespace DRT


namespace DRT::ELEMENTS::SHELL
{
  /*!
   * @brief Evaluates a Neumann condition @p condition for the element @p element.
   *
   * @note This function determines the shape of the element at runtime and calls the respective
   * templated version of @p EvaluateNeumann. If you already know the CORE::FE::CellType
   * of the element at compile-time, you could directly call @EvaluateNeumann.
   *
   * @param element (in) : The element where we integrate
   * @param discretization (in) : Discretization
   * @param condition (in) : The Neumann condition to be evaluated within the element.
   * @param dof_index_array (in) : The index array of the dofs of the element
   * @param element_force_vector (in/out) : The element force vector for the evaluated Neumann
   * condition
   * @param element_stiffness_matrix (in/out) : The element stiffness matrix for the evaluated
   * Neumann condition
   * @param total_time (in) : The total time for time dependent Neumann conditions
   */
  void EvaluateNeumannByElement(DRT::Element& element, const DRT::Discretization& discretization,
      DRT::Condition& condition, const std::vector<int>& dof_index_array,
      CORE::LINALG::SerialDenseVector& element_force_vector,
      CORE::LINALG::SerialDenseMatrix* element_stiffness_matrix, double total_time);

  /*!
   * @brief Evaluates a Neumann condition @p condition for the element @p element with the
   * discretization type known at compile time.
   *
   *
   * @tparam distype Discretization type known at compile time
   *
   * @param element (in) : The element where we integrate
   * @param discretization (in) : Discretization
   * @param condition (in) : The Neumann condition to be evaluated within the element.
   * @param dof_index_array (in) : The index array of the dofs of the element
   * @param element_force_vector (in/out) : The element force vector for the evaluated Neumann
   * condition
   * @param element_stiffness_matrix (in/out) : The element stiffness matrix for the evaluated
   * Neumann condition
   * @param total_time (in) : The total time for time dependent Neumann conditions
   */
  template <CORE::FE::CellType distype>
  void EvaluateNeumann(DRT::Element& element, const DRT::Discretization& discretization,
      DRT::Condition& condition, const std::vector<int>& dof_index_array,
      CORE::LINALG::SerialDenseVector& element_force_vector,
      CORE::LINALG::SerialDenseMatrix* element_stiffness_matrix, double total_time);

}  // namespace DRT::ELEMENTS::SHELL

FOUR_C_NAMESPACE_CLOSE

#endif
