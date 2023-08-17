/*----------------------------------------------------------------------*/
/*! \file
\brief search tree input parameters
\level 2
*/

/*----------------------------------------------------------------------*/
#ifndef BACI_INPAR_SEARCHTREE_H
#define BACI_INPAR_SEARCHTREE_H

#include "baci_inpar_parameterlist_utils.H"

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
#endif  // INPAR_SEARCHTREE_H
