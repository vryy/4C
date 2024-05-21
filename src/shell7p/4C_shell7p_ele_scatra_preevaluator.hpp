/*! \file

\brief Preevalution methods for Shell-ScaTra elements

\level 1
*/

#ifndef FOUR_C_SHELL7P_ELE_SCATRA_PREEVALUATOR_HPP
#define FOUR_C_SHELL7P_ELE_SCATRA_PREEVALUATOR_HPP

#include "4C_config.hpp"

#include "4C_inpar_scatra.hpp"
#include "4C_lib_element.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  class Discretization;

}  // namespace DRT


namespace DRT::ELEMENTS::SHELL
{
  /*!
   * \brief Preevaluate anything that needs to be done before the standard Evaluate
   *
   * For the scatra coupling we need for example the concentrations. Empty function in the
   * base class that may be overloaded in derived elements.
   *
   * @param ele  (in) : Reference to the element
   * @param discretization (in) : Reference to the discretization
   * @param dof_index_array (in) : The location array of the owned dofs
   */
  void PreEvaluateScatraByElement(DRT::Element& ele, Teuchos::ParameterList& params,
      DRT::Discretization& discretization, DRT::Element::LocationArray& dof_index_array);

  /*!
   * @brief Preevaluate anything that needs to be done before the standard Evaluate of the element
   * @p element with the discretization type known at compile time.
   *
   *
   * @tparam distype : Discretization type known at compile time
   *
   * @param ele  (in) : Reference to the element
   * @param discretization (in) : Reference to the discretization
   * @param dof_index_array (in) : The location array of the owned dofs
   */
  template <CORE::FE::CellType distype>
  void PreEvaluateScatra(DRT::Element& ele, Teuchos::ParameterList& params,
      DRT::Discretization& discretization, DRT::Element::LocationArray& dof_index_array);

}  // namespace DRT::ELEMENTS::SHELL

FOUR_C_NAMESPACE_CLOSE

#endif
