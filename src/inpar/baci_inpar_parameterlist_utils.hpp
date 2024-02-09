/*----------------------------------------------------------------------------*/
/*! \file
\brief Utility routines for the Teuchos::ParameterList

\level 1
 */
/*----------------------------------------------------------------------------*/

#ifndef BACI_INPAR_PARAMETERLIST_UTILS_HPP
#define BACI_INPAR_PARAMETERLIST_UTILS_HPP

#include "baci_config.hpp"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

BACI_NAMESPACE_OPEN


namespace INPUT
{
  /// cast internal integer to enum type
  template <class T>
  T IntegralValue(const Teuchos::ParameterList& params, const std::string& name)
  {
    int value = Teuchos::getIntegralValue<int>(params, name);
    return static_cast<T>(value);
  }

  template <class T>
  T get(const Teuchos::ParameterList& params, const std::string& name)
  {
    int value = params.get<int>(name);
    return static_cast<T>(value);
  }

  template <class T>
  T get(Teuchos::ParameterList& params, const std::string& name, T default_value)
  {
    int value = params.get<int>(name, default_value);
    return static_cast<T>(value);
  }
}  // namespace INPUT

BACI_NAMESPACE_CLOSE

#endif
