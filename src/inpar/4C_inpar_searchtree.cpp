/*----------------------------------------------------------------------*/
/*! \file
\brief search tree input parameters
\level 2
*/

/*----------------------------------------------------------------------*/



#include "4C_inpar_searchtree.hpp"

#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN



void Inpar::Geo::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& search_tree = list->sublist("SEARCH TREE", false, "");

  setStringToIntegralParameter<int>("TREE_TYPE", "notree", "set tree type",
      tuple<std::string>("notree", "octree3d", "quadtree3d", "quadtree2d"),
      tuple<int>(
          Inpar::Geo::Notree, Inpar::Geo::Octree3D, Inpar::Geo::Quadtree3D, Inpar::Geo::Quadtree2D),
      &search_tree);
}

FOUR_C_NAMESPACE_CLOSE
