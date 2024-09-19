/*----------------------------------------------------------------------*/
/*! \file
\brief search tree input parameters
\level 2
*/

/*----------------------------------------------------------------------*/



#include "4C_inpar_searchtree.hpp"

#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN



void Inpar::Geo::set_valid_parameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& search_tree = list->sublist("SEARCH TREE", false, "");

  setStringToIntegralParameter<Inpar::Geo::TreeType>("TREE_TYPE", "notree", "set tree type",
      tuple<std::string>("notree", "octree3d", "quadtree3d", "quadtree2d"),
      tuple<Inpar::Geo::TreeType>(
          Inpar::Geo::Notree, Inpar::Geo::Octree3D, Inpar::Geo::Quadtree3D, Inpar::Geo::Quadtree2D),
      &search_tree);
}

FOUR_C_NAMESPACE_CLOSE
