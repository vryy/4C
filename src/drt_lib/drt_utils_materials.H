/*---------------------------------------------------------------------*/
/*!

\brief A collection of helper methods for namespace DRT

\level 1

\maintainer Martin Kronbichler

*/
/*---------------------------------------------------------------------*/

#ifndef UTILS_MATERIALS_H
#define UTILS_MATERIALS_H


#include "drt_node.H"
#include "drt_element.H"

#include "Teuchos_RCP.hpp"

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
        Teuchos::RCP<DRT::Discretization>
            dis,  ///< discretization of which material configuration shall be updated
        Teuchos::RCP<const Epetra_Vector>
            disp  ///< displacement vector which is used for moving the nodes
    );

  }  // namespace UTILS
}  // namespace DRT


#endif  // #ifndef UTILS_MATERIALS_H
