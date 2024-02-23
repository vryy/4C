/*---------------------------------------------------------------------*/
/*! \file
\brief A collection of helper functions for Teuchos::ParameterLists

\level 0

*/
/*---------------------------------------------------------------------*/
#include "baci_utils_parameter_list.hpp"

BACI_NAMESPACE_OPEN

namespace CORE::UTILS
{
  void BoolParameter(std::string const& paramName, std::string const& value,
      std::string const& docString, Teuchos::ParameterList* paramList)
  {
    Teuchos::Array<std::string> yesnotuple =
        Teuchos::tuple<std::string>("Yes", "No", "yes", "no", "YES", "NO");
    Teuchos::Array<int> yesnovalue = Teuchos::tuple<int>(true, false, true, false, true, false);
    Teuchos::setStringToIntegralParameter<int>(
        paramName, value, docString, yesnotuple, yesnovalue, paramList);
  }


  void IntParameter(std::string const& paramName, int const value, std::string const& docString,
      Teuchos::ParameterList* paramList)
  {
    Teuchos::AnyNumberParameterEntryValidator::AcceptedTypes validator(false);
    validator.allowInt(true);
    Teuchos::setIntParameter(paramName, value, docString, paramList, validator);
  }


  void DoubleParameter(std::string const& paramName, double const& value,
      std::string const& docString, Teuchos::ParameterList* paramList)
  {
    Teuchos::AnyNumberParameterEntryValidator::AcceptedTypes validator(false);
    validator.allowDouble(true);
    validator.allowInt(true);
    Teuchos::setDoubleParameter(paramName, value, docString, paramList, validator);
  }


  void StringParameter(std::string const& paramName, std::string const& value,
      std::string const& docString, Teuchos::ParameterList* paramList)
  {
    Teuchos::RCP<Teuchos::StringValidator> validator = Teuchos::rcp(new Teuchos::StringValidator());
    paramList->set(paramName, value, docString, validator);
  }

}  // namespace CORE::UTILS

BACI_NAMESPACE_CLOSE
