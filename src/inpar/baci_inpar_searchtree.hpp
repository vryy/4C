/*----------------------------------------------------------------------*/
/*! \file
\brief search tree input parameters
\level 2
*/

/*----------------------------------------------------------------------*/
#ifndef FOUR_C_INPAR_SEARCHTREE_HPP
#define FOUR_C_INPAR_SEARCHTREE_HPP

#include "baci_config.hpp"

#include "baci_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
namespace INPAR
{
  namespace GEO
  {
    /// specify tree type
    enum TreeType
    {
      Notree,
      Octree3D,
      Quadtree3D,
      Quadtree2D
    };

    /// set the searchtree parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace GEO

}  // namespace INPAR

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
