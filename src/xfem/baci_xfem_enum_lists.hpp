/*----------------------------------------------------------------------------*/
/** \file

\brief list of important enumerators ( used for the necessary header inclusions
       till C++11 )


\level 3

*/
/*----------------------------------------------------------------------------*/


#ifndef BACI_XFEM_ENUM_LISTS_HPP
#define BACI_XFEM_ENUM_LISTS_HPP

#include "baci_config.hpp"

#include "baci_utils_exceptions.hpp"

#include <string>

BACI_NAMESPACE_OPEN

namespace XFEM
{
  /// supported field names
  enum FieldName
  {
    unknown = -1,
    structure = 0,
    xstructure = 1
  };

  /// map field name enumerator to string
  static inline std::string FieldName2String(const enum FieldName& field)
  {
    switch (field)
    {
      case structure:
        return "structure";
        break;
      case xstructure:
        return "xstructure";
        break;
      default:
        dserror("Unknown FieldName enumerator (enum = %d)!", field);
        break;
    }
    return "";
  }

  /// map field name enumerator to string
  static inline enum FieldName String2FieldName(const std::string& name)
  {
    enum FieldName field = unknown;
    if (name == "structure")
      field = structure;
    else if (name == "xstructure")
      field = xstructure;
    else
      dserror(
          "No known conversion for the given discretization "
          "name \"%s\"!",
          name.c_str());

    return field;
  }

  /// map types
  enum MapType
  {
    map_dofs = 0,  ///< extract/insert DoF's
    map_nodes = 1  ///< extract/insert nodes
  };

  namespace MULTIFIELD
  {
    /// block type enumerator
    enum BlockType
    {
      block_interface = 0,     ///< interface block
      block_non_interface = 1  ///< non-interface block
    };
  }  // namespace MULTIFIELD
}  // namespace XFEM


BACI_NAMESPACE_CLOSE

#endif  // XFEM_ENUM_LISTS_H
