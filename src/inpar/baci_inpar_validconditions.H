/*----------------------------------------------------------------------*/
/*! \file

\brief Setup of the list of valid conditions for input

\level 1

*/
/*----------------------------------------------------------------------*/



#ifndef BACI_INPAR_VALIDCONDITIONS_H
#define BACI_INPAR_VALIDCONDITIONS_H

#include "baci_config.H"

#include <Teuchos_Array.hpp>
#include <Teuchos_RCP.hpp>

#include <iostream>
#include <string>
#include <vector>

BACI_NAMESPACE_OPEN


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

BACI_NAMESPACE_CLOSE

#endif
