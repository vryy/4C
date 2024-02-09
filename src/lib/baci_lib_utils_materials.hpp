/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for namespace DRT

\level 1


*/
/*---------------------------------------------------------------------*/

#ifndef BACI_LIB_UTILS_MATERIALS_HPP
#define BACI_LIB_UTILS_MATERIALS_HPP


#include "baci_config.hpp"

#include "baci_lib_element.hpp"
#include "baci_lib_node.hpp"

#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN

namespace DRT
{
  namespace UTILS
  {
    /*! brief Update material configuration with displacements stored in dofrowmap format
     *
     *
     *
     */
    void UpdateMaterialConfigWithDispVector(
        const Teuchos::RCP<DRT::Discretization>&
            dis,  ///< discretization of which material configuration shall be updated
        const Teuchos::RCP<const Epetra_Vector>&
            disp  ///< displacement vector which is used for moving the nodes
    );

  }  // namespace UTILS
}  // namespace DRT


BACI_NAMESPACE_CLOSE

#endif  // LIB_UTILS_MATERIALS_H
