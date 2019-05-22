/*----------------------------------------------------------------------------*/
/*!
\brief Utility routines for the Teuchos::ParameterList

\level 1
\maintainer Martin Kronbichler
 */
/*----------------------------------------------------------------------------*/

#ifndef INPAR_PARAMETERLIST_UTILS_H
#define INPAR_PARAMETERLIST_UTILS_H

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

namespace DRT
{
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
}  // namespace DRT

#endif
