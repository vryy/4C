/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for namespace DRT

\level 1


*/
/*---------------------------------------------------------------------*/

#ifndef BACI_DISCRETIZATION_GEOMETRY_UPDATE_MATERIAL_CONFIG_HPP
#define BACI_DISCRETIZATION_GEOMETRY_UPDATE_MATERIAL_CONFIG_HPP


#include "baci_config.hpp"

#include "baci_lib_element.hpp"
#include "baci_lib_node.hpp"

#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN

namespace CORE::GEO
{
  //! Update material configuration of @p dis with @p disp
  void UpdateMaterialConfigWithDispVector(
      Teuchos::RCP<const DRT::Discretization> dis, Teuchos::RCP<const Epetra_Vector> disp);
}  // namespace CORE::GEO

BACI_NAMESPACE_CLOSE

#endif  // LIB_UTILS_MATERIALS_H
