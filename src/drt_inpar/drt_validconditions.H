/*----------------------------------------------------------------------*/
/*!

\brief Setup of the list of valid conditions for input

\level 1

\maintainer Martin Kronbichler
*/
/*----------------------------------------------------------------------*/



#ifndef DRT_VALIDCONDITIONS_H
#define DRT_VALIDCONDITIONS_H

#include <vector>
#include <string>
#include <iostream>

#include <Teuchos_Array.hpp>
#include <Teuchos_RCP.hpp>


namespace DRT
{
  namespace INPUT
  {
    class ConditionDefinition;

    /// construct list with all conditions and documentation
    Teuchos::RCP<std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>> ValidConditions();

    /// print all known condition sections without contents
    void PrintEmptyConditionDefinitions(std::ostream& stream,
        std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist, bool color = true);

  }  // namespace INPUT
}  // namespace DRT

/// print empty condition sections
void PrintConditionDatHeader();

#endif
