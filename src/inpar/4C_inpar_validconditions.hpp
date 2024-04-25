/*----------------------------------------------------------------------*/
/*! \file

\brief Setup of the list of valid conditions for input

\level 1

*/
/*----------------------------------------------------------------------*/



#ifndef FOUR_C_INPAR_VALIDCONDITIONS_HPP
#define FOUR_C_INPAR_VALIDCONDITIONS_HPP

#include "4C_config.hpp"

#include <Teuchos_Array.hpp>
#include <Teuchos_RCP.hpp>

#include <iostream>
#include <string>
#include <vector>

FOUR_C_NAMESPACE_OPEN


namespace INPUT
{
  class ConditionDefinition;

  /// construct list with all conditions and documentation
  Teuchos::RCP<std::vector<Teuchos::RCP<INPUT::ConditionDefinition>>> ValidConditions();

  /// print all known condition sections without contents
  void PrintEmptyConditionDefinitions(
      std::ostream& stream, std::vector<Teuchos::RCP<INPUT::ConditionDefinition>>& condlist);

}  // namespace INPUT

/// print empty condition sections
void PrintConditionDatHeader();

FOUR_C_NAMESPACE_CLOSE

#endif
