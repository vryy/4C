/*----------------------------------------------------------------------*/
/*!
\file inpar_searchtree.cpp
\brief search tree input parameters
\level 2
\maintainer Martin Kronbichler
*/

/*----------------------------------------------------------------------*/



#include "drt_validparameters.H"
#include "inpar_searchtree.H"



void INPAR::GEO::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;

  Teuchos::ParameterList& search_tree = list->sublist("SEARCH TREE",false,"");

  setStringToIntegralParameter<int>("TREE_TYPE","notree","set tree type",
                                   tuple<std::string>("notree","octree3d","quadtree3d","quadtree2d"),
                                   tuple<int>(
                                     INPAR::GEO::Notree,
                                     INPAR::GEO::Octree3D,
                                     INPAR::GEO::Quadtree3D,
                                     INPAR::GEO::Quadtree2D),
                                     &search_tree);
}
