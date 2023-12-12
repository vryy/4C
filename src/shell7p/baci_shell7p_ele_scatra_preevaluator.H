/*! \file

\brief Preevalution methods for Shell-ScaTra elements

\level 1
*/

#ifndef BACI_SHELL7P_ELE_SCATRA_PREEVALUATOR_H
#define BACI_SHELL7P_ELE_SCATRA_PREEVALUATOR_H

#include "baci_inpar_scatra.H"
#include "baci_lib_element.H"

#include <vector>

BACI_NAMESPACE_OPEN

namespace DRT
{
  class Discretization;
  class Condition;

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

BACI_NAMESPACE_CLOSE

#endif
