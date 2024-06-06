/*! \file

\brief Preevalution methods for Shell-ScaTra elements

\level 1
*/

#ifndef FOUR_C_SHELL7P_ELE_SCATRA_PREEVALUATOR_HPP
#define FOUR_C_SHELL7P_ELE_SCATRA_PREEVALUATOR_HPP

#include "4C_config.hpp"

#include "4C_discretization_fem_general_element.hpp"
#include "4C_inpar_scatra.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  class Discretization;

}  // namespace Discret


namespace Discret::ELEMENTS::Shell
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
  void PreEvaluateScatraByElement(Core::Elements::Element& ele, Teuchos::ParameterList& params,
      Discret::Discretization& discretization,
      Core::Elements::Element::LocationArray& dof_index_array);

  /*!
   * @brief Preevaluate anything that needs to be done before the standard Evaluate of the element
   * @p element with the discretization type known at compile time.
   *
   *
   * @tparam distype : discretization type known at compile time
   *
   * @param ele  (in) : Reference to the element
   * @param discretization (in) : Reference to the discretization
   * @param dof_index_array (in) : The location array of the owned dofs
   */
  template <Core::FE::CellType distype>
  void PreEvaluateScatra(Core::Elements::Element& ele, Teuchos::ParameterList& params,
      Discret::Discretization& discretization,
      Core::Elements::Element::LocationArray& dof_index_array);

}  // namespace Discret::ELEMENTS::Shell

FOUR_C_NAMESPACE_CLOSE

#endif
