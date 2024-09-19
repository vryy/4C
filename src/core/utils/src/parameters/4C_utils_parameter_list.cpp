/*---------------------------------------------------------------------*/
/*! \file
\brief A collection of helper functions for Teuchos::ParameterLists

\level 0

*/
/*---------------------------------------------------------------------*/
#include "4C_utils_parameter_list.hpp"

#include <Teuchos_Array.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_Tuple.hpp>

#include <array>

FOUR_C_NAMESPACE_OPEN

namespace Core::UTILS
{
  void bool_parameter(std::string const& paramName, std::string const& value,
      std::string const& docString, Teuchos::ParameterList* paramList)
  {
    Teuchos::Array<std::string> yesnotuple =
        Teuchos::tuple<std::string>("Yes", "No", "yes", "no", "YES", "NO");

    std::array<const bool, 6> yesnoarray = {true, false, true, false, true, false};
    Teuchos::ArrayView<const bool> yesnovalue =
        Teuchos::arrayView(yesnoarray.data(), yesnoarray.size());

    Teuchos::setStringToIntegralParameter<bool>(
        paramName, value, docString, yesnotuple, yesnovalue, paramList);
  }


  void int_parameter(std::string const& paramName, int const value, std::string const& docString,
      Teuchos::ParameterList* paramList)
  {
    const bool allow_all_types_per_default = false;
    Teuchos::AnyNumberParameterEntryValidator::AcceptedTypes validator(allow_all_types_per_default);
    validator.allowInt(true);
    Teuchos::setIntParameter(paramName, value, docString, paramList, validator);
  }


  void double_parameter(std::string const& paramName, double const& value,
      std::string const& docString, Teuchos::ParameterList* paramList)
  {
    // Create a validator that does not allow all types by default
    const bool allow_all_types_per_default = false;
    Teuchos::AnyNumberParameterEntryValidator::AcceptedTypes validator(allow_all_types_per_default);

    // Explicitly allow only double type
    validator.allowDouble(true);
    validator.allowInt(true);

    // Set the double parameter in the parameter list with the validator
    Teuchos::setDoubleParameter(paramName, value, docString, paramList, validator);
  }


  void string_parameter(std::string const& paramName, std::string const& value,
      std::string const& docString, Teuchos::ParameterList* paramList,
      std::vector<std::string> const& validParams)
  {
    Teuchos::RCP<Teuchos::StringValidator> validator = Teuchos::rcp(new Teuchos::StringValidator());
    // Set valid strings only if validParams is not empty
    if (!validParams.empty())
    {
      validator->setValidStrings(validParams);
    }
    // Set the parameter in the parameter list
    paramList->set(paramName, value, docString, validator);
  }

}  // namespace Core::UTILS

FOUR_C_NAMESPACE_CLOSE
