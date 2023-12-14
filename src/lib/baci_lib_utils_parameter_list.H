/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for namespace DRT

\level 0


*/
/*---------------------------------------------------------------------*/

#ifndef BACI_LIB_UTILS_PARAMETER_LIST_H
#define BACI_LIB_UTILS_PARAMETER_LIST_H

#include <Teuchos_StandardParameterEntryValidators.hpp>

BACI_NAMESPACE_OPEN

namespace DRT
{
  namespace UTILS
  {
    //! add entry as item of enum class @p value to @p list with name @p parameter_name
    template <class enumtype>
    void AddEnumClassToParameterList(
        const std::string& parameter_name, const enumtype value, Teuchos::ParameterList& list)
    {
      const std::string docu = "";
      const std::string value_name = "val";
      Teuchos::setStringToIntegralParameter<enumtype>(parameter_name, value_name, docu,
          Teuchos::tuple<std::string>(value_name), Teuchos::tuple<enumtype>(value), &list);
    }
  }  // namespace UTILS
}  // namespace DRT


BACI_NAMESPACE_CLOSE

#endif  // LIB_UTILS_PARAMETER_LIST_H
