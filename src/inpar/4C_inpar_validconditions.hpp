#ifndef FOUR_C_INPAR_VALIDCONDITIONS_HPP
#define FOUR_C_INPAR_VALIDCONDITIONS_HPP

#include "4C_config.hpp"

#include <Teuchos_Array.hpp>
#include <Teuchos_RCP.hpp>

#include <iostream>
#include <string>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::Conditions
{
  class ConditionDefinition;
}


namespace Input
{
  /// construct list with all conditions and documentation
  Teuchos::RCP<std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>> valid_conditions();

  /// print all known condition sections without contents
  void print_empty_condition_definitions(std::ostream& stream,
      std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>& condlist);

}  // namespace Input

/// print empty condition sections
void print_condition_dat_header();

FOUR_C_NAMESPACE_CLOSE

#endif
