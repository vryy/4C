/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for namespace DRT

\level 1


*/
/*---------------------------------------------------------------------*/

#ifndef BACI_LIB_UTILS_MATERIALS_H
#define BACI_LIB_UTILS_MATERIALS_H


#include "baci_lib_node.H"
#include "baci_lib_element.H"

#include <Teuchos_RCP.hpp>

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


#endif  // LIB_UTILS_MATERIALS_H
