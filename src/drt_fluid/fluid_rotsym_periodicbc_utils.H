/*-----------------------------------------------------------*/
/*!

\brief Methods needed to apply rotationally symmetric periodic boundary
       conditions for fluid problems

\maintainer Martin Kronbichler

\level 1

*/
/*-----------------------------------------------------------*/

#ifndef FLUID_ROTSYM_PERIODICBC_UTILS_H
#define FLUID_ROTSYM_PERIODICBC_UTILS_H


#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>


// forward declarations
namespace DRT
{
  class Discretization;
  class Node;
  class Condition;
}  // namespace DRT

namespace FLD
{
  //! return specific component of rotated global vector result
  double GetComponentOfRotatedVectorField(const int idf,  ///< the component id 0,1 (and 2)
      const Teuchos::RCP<const Epetra_Vector> proc0data,  ///< the data vector
      const int lid,         ///< local dof id (w.r.t proc0data) of component idf
      const double rotangle  ///< angle of rotation (RAD)
  );

  //! Is given node a slave node of rotationally symmetric periodic boundary conditions?
  bool IsSlaveNodeOfRotSymPBC(const DRT::Node* node,  ///< the node
      double& rotangle                                ///< the angle of slave plane rotation (RAD)
  );

  //! Access angle of rotation and convert it to RAD
  inline double GetRotAngleFromCondition(
      const DRT::Condition* cond  ///< pointer to desired periodic boundary condition
  );

  //! Get all relevant slave nodes of rotationally symmetric periodic bc's
  void GetRelevantSlaveNodesOfRotSymPBC(
      std::map<int, double>&
          pbcslavenodemap,  ///< map to be filled with node gids and rotation angles
      Teuchos::RCP<DRT::Discretization> dis);  ///< discretization

}  // namespace FLD

#endif  // #ifndef FLUID_ROTSYM_PERIODICBC_UTILS_H
