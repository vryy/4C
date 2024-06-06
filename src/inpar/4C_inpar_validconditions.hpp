/*----------------------------------------------------------------------*/
/*! \file

\brief Setup of the list of valid conditions for input

\level 1

*/
/*----------------------------------------------------------------------*/



#ifndef FOUR_C_INPAR_VALIDCONDITIONS_HPP
#define FOUR_C_INPAR_VALIDCONDITIONS_HPP

#include "4C_config.hpp"

#include "4C_discretization_condition_definition.hpp"

#include <Teuchos_Array.hpp>
#include <Teuchos_RCP.hpp>

#include <iostream>
#include <string>
#include <vector>

FOUR_C_NAMESPACE_OPEN


namespace Input
{
  /// construct list with all conditions and documentation
  Teuchos::RCP<std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>> ValidConditions();

  /// print all known condition sections without contents
  void PrintEmptyConditionDefinitions(std::ostream& stream,
      std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>& condlist);

}  // namespace Input

/// print empty condition sections
void PrintConditionDatHeader();

FOUR_C_NAMESPACE_CLOSE

#endif
