/*!
 \file searchtree_geometry_service.cpp

 \brief provides geometry methods for a search tree

 \level 1

\maintainer Martin Kronbichler
 */

#include "element_coordtrafo.H"
#include "element_normals.H"
#include "intersection_interfacepoint.H"
#include "intersection_service_templates.H"
#include "intersection_service.H"
#include "position_array.H"
#include "searchtree_geometry_service.H"
#include "searchtree_nearestobject.H"
#include "../drt_lib/drt_discret.H"

/*----------------------------------------------------------------------*
 | delivers a axis-aligned bounding box for a given          peder 07/08|
 | discretization                                                       |
 *----------------------------------------------------------------------*/
LINALG::Matrix<3, 2> GEO::getXAABBofDis(
    const DRT::Discretization& dis, const std::map<int, LINALG::Matrix<3, 1>>& currentpositions)
{
  LINALG::Matrix<3, 2> XAABB(true);
  if (dis.NumGlobalElements() == 0) return XAABB;

  if (dis.NumMyColElements() == 0) return XAABB;

  // initialize XAABB as rectangle around the first point of dis
  const int nodeid = dis.lColElement(0)->Nodes()[0]->Id();
  const LINALG::Matrix<3, 1> pos = currentpositions.find(nodeid)->second;
  for (int dim = 0; dim < 3; ++dim)
  {
    XAABB(dim, 0) = pos(dim) - GEO::TOL7;
    XAABB(dim, 1) = pos(dim) + GEO::TOL7;
  }

  // TODO make more efficient by cahcking nodes
  // loop over elements and merge XAABB with their eXtendedAxisAlignedBoundingBox
  for (int j = 0; j < dis.NumMyColElements(); ++j)
  {
    const DRT::Element* element = dis.lColElement(j);
    const LINALG::SerialDenseMatrix xyze_element(
        GEO::getCurrentNodalPositions(element, currentpositions));
    GEO::EleGeoType eleGeoType(GEO::HIGHERORDER);
    GEO::checkRoughGeoType(element, xyze_element, eleGeoType);

    const LINALG::Matrix<3, 2> xaabbEle =
        GEO::computeFastXAABB(element->Shape(), xyze_element, eleGeoType);
    XAABB = mergeAABB(XAABB, xaabbEle);
  }
  return XAABB;
}

/*----------------------------------------------------------------------*
 | delivers a slightly enlarged axis-aligned bounding box for u.may09/09|
 | given elements with their current postions for sliding ALE           |
 *----------------------------------------------------------------------*/
LINALG::Matrix<3, 2> GEO::getXAABBofEles(std::map<int, Teuchos::RCP<DRT::Element>>& elements,
    const std::map<int, LINALG::Matrix<3, 1>>& currentpositions)
{
  LINALG::Matrix<3, 2> XAABB(true);
  if (elements.begin() == elements.end()) return XAABB;

  // initialize XAABB as rectangle around the first point of the first element
  const int nodeid = elements.begin()->second->Nodes()[0]->Id();
  const LINALG::Matrix<3, 1> pos = currentpositions.find(nodeid)->second;
  for (int dim = 0; dim < 3; ++dim)
  {
    XAABB(dim, 0) = pos(dim) - GEO::TOL7;
    XAABB(dim, 1) = pos(dim) + GEO::TOL7;
  }

  std::map<int, Teuchos::RCP<DRT::Element>>::const_iterator elemiter;
  for (elemiter = elements.begin(); elemiter != elements.end(); ++elemiter)
  {
    Teuchos::RCP<DRT::Element> currelement = elemiter->second;
    const LINALG::SerialDenseMatrix xyze_element(
        GEO::getCurrentNodalPositions(currelement, currentpositions));
    GEO::EleGeoType eleGeoType(GEO::HIGHERORDER);
    GEO::checkRoughGeoType(currelement, xyze_element, eleGeoType);
    const LINALG::Matrix<3, 2> xaabbEle =
        GEO::computeFastXAABB(currelement->Shape(), xyze_element, eleGeoType);
    XAABB = mergeAABB(XAABB, xaabbEle);
  }

  return XAABB;
}

/*----------------------------------------------------------------------*
 | delivers a axis-aligned bounding box for a given          u.may 12/09|
 | discretization                                                       |
 *----------------------------------------------------------------------*/
LINALG::Matrix<3, 2> GEO::getXAABBofDisPar(const DRT::Discretization& dis,
    const std::map<int, LINALG::Matrix<3, 1>>& currentpositions, double cutoff)
{
  LINALG::Matrix<3, 2> XAABB(true);
  if (dis.NumGlobalElements() == 0) return XAABB;

  if (dis.NumMyColElements() == 0)
  {
    dserror(
        "boundary should be ghosted and xfem discretization balanced, such that there are at least "
        "some xfem elements!");
  }

  // initialize XAABB as rectangle around the first point of dis
  const int nodeid = dis.lColElement(0)->Nodes()[0]->Id();
  const LINALG::Matrix<3, 1> pos = currentpositions.find(nodeid)->second;
  for (int dim = 0; dim < 3; ++dim)
  {
    XAABB(dim, 0) = pos(dim) - GEO::TOL7;
    XAABB(dim, 1) = pos(dim) + GEO::TOL7;
  }

  // TODO make more efficient by cahcking nodes
  // loop over elements and merge XAABB with their eXtendedAxisAlignedBoundingBox
  for (int j = 0; j < dis.NumMyColElements(); ++j)
  {
    const DRT::Element* element = dis.lColElement(j);
    const LINALG::SerialDenseMatrix xyze_element(
        GEO::getCurrentNodalPositions(element, currentpositions));
    GEO::EleGeoType eleGeoType(GEO::HIGHERORDER);
    GEO::checkRoughGeoType(element, xyze_element, eleGeoType);
    const LINALG::Matrix<3, 2> xaabbEle =
        GEO::computeFastXAABB(element->Shape(), xyze_element, eleGeoType);
    XAABB = mergeAABB(XAABB, xaabbEle);
  }
  // enlarge box by 5 times the cutoff-radius
  for (int dim = 0; dim < 3; ++dim)
  {
    XAABB(dim, 0) = XAABB(dim, 0) - 5 * cutoff;
    XAABB(dim, 1) = XAABB(dim, 1) + 5 * cutoff;
  }
  return XAABB;
}

/*----------------------------------------------------------------------*
 | delivers a axis-aligned bounding box for a given                     |
 | discretization in reference configuration                            |
 *----------------------------------------------------------------------*/
LINALG::Matrix<3, 2> GEO::getXAABBofDisPar(const DRT::Discretization& dis, double cutoff)
{
  std::map<int, LINALG::Matrix<3, 1>> currentpositions;

  for (int lid = 0; lid < dis.NumMyColNodes(); ++lid)
  {
    const DRT::Node* node = dis.lColNode(lid);
    LINALG::Matrix<3, 1> currpos;
    currpos(0) = node->X()[0];
    currpos(1) = node->X()[1];
    currpos(2) = node->X()[2];
    currentpositions[node->Id()] = currpos;
  }

  return getXAABBofDisPar(dis, currentpositions, cutoff);
}

/*----------------------------------------------------------------------*
 | delivers a axis-aligned bounding box for a given          peder 07/08|
 | discretization in reference configuration                            |
 *----------------------------------------------------------------------*/
LINALG::Matrix<3, 2> GEO::getXAABBofDis(const DRT::Discretization& dis)
{
  std::map<int, LINALG::Matrix<3, 1>> currentpositions;

  for (int lid = 0; lid < dis.NumMyColNodes(); ++lid)
  {
    const DRT::Node* node = dis.lColNode(lid);
    LINALG::Matrix<3, 1> currpos;
    currpos(0) = node->X()[0];
    currpos(1) = node->X()[1];
    currpos(2) = node->X()[2];
    currentpositions[node->Id()] = currpos;
  }

  return getXAABBofDis(dis, currentpositions);
}

/*----------------------------------------------------------------------*
 | delivers a axis-aligned bounding box for the nodes in     ghamm 09/12|
 | a given discretization in reference configuration                    |
 *----------------------------------------------------------------------*/
LINALG::Matrix<3, 2> GEO::getXAABBofNodes(const DRT::Discretization& dis)
{
  LINALG::Matrix<3, 2> XAABB(true);
  if (dis.NumMyRowNodes() == 0) dserror("no row node found on this proc");  // return XAABB;

  {
    // initialize XAABB as rectangle around the first node of dis
    const DRT::Node* node = dis.lRowNode(0);
    for (int dim = 0; dim < 3; ++dim)
    {
      XAABB(dim, 0) = node->X()[dim] - GEO::TOL7;
      XAABB(dim, 1) = node->X()[dim] + GEO::TOL7;
    }
  }

  // loop over remaining nodes and merge XAABB with their eXtendedAxisAlignedBoundingBox
  for (int lid = 1; lid < dis.NumMyRowNodes(); ++lid)
  {
    const DRT::Node* node = dis.lRowNode(lid);

    for (int dim = 0; dim < 3; dim++)
    {
      XAABB(dim, 0) = std::min(XAABB(dim, 0), node->X()[dim] - TOL7);
      XAABB(dim, 1) = std::max(XAABB(dim, 1), node->X()[dim] + TOL7);
    }
  }
  return XAABB;
}

/*----------------------------------------------------------------------*
 | delivers an axis-aligned bounding box for coords          ghamm 11/12|
 *----------------------------------------------------------------------*/
LINALG::Matrix<3, 2> GEO::getXAABBofPositions(
    const std::map<int, LINALG::Matrix<3, 1>>& currentpositions)
{
  LINALG::Matrix<3, 2> XAABB(true);

  if (currentpositions.size() == 0) dserror("map with current positions is emtpy");

  // initialize XAABB as rectangle around the first node
  const LINALG::Matrix<3, 1> initcurrentpos = currentpositions.begin()->second;
  for (int dim = 0; dim < 3; ++dim)
  {
    XAABB(dim, 0) = initcurrentpos(dim) - GEO::TOL7;
    XAABB(dim, 1) = initcurrentpos(dim) + GEO::TOL7;
  }

  // loop over remaining entries and merge XAABB with their eXtendedAxisAlignedBoundingBox
  std::map<int, LINALG::Matrix<3, 1>>::const_iterator iter;
  for (iter = currentpositions.begin(); iter != currentpositions.end(); ++iter)
  {
    const LINALG::Matrix<3, 1> currentpos = iter->second;
    for (int dim = 0; dim < 3; dim++)
    {
      XAABB(dim, 0) = std::min(XAABB(dim, 0), currentpos(dim) - TOL7);
      XAABB(dim, 1) = std::max(XAABB(dim, 1), currentpos(dim) + TOL7);
    }
  }
  return XAABB;
}

/*----------------------------------------------------------------------*
 | compute AABB s for all labeled strutcures in the          u.may 08/08|
 | element list                                                         |
 *----------------------------------------------------------------------*/
std::vector<LINALG::Matrix<3, 2>> GEO::computeXAABBForLabeledStructures(
    const DRT::Discretization& dis, const std::map<int, LINALG::Matrix<3, 1>>& currentpositions,
    const std::map<int, std::set<int>>& elementList)
{
  std::vector<LINALG::Matrix<3, 2>> XAABBs;

  for (std::map<int, std::set<int>>::const_iterator labelIter = elementList.begin();
       labelIter != elementList.end(); labelIter++)
  {
    LINALG::Matrix<3, 2> xaabb_label;
    // initialize xaabb_label with box around first point
    const int eleId = *((labelIter->second).begin());
    const int nodeId = dis.gElement(eleId)->Nodes()[0]->Id();
    const LINALG::Matrix<3, 1> pos = currentpositions.find(nodeId)->second;
    for (int dim = 0; dim < 3; ++dim)
    {
      xaabb_label(dim, 0) = pos(dim) - GEO::TOL7;
      xaabb_label(dim, 1) = pos(dim) + GEO::TOL7;
    }
    // run over set elements
    for (std::set<int>::const_iterator eleIter = (labelIter->second).begin();
         eleIter != (labelIter->second).end(); eleIter++)
    {
      const DRT::Element* element = dis.gElement(*eleIter);
      const LINALG::SerialDenseMatrix xyze_element(
          GEO::getCurrentNodalPositions(element, currentpositions));
      GEO::EleGeoType eleGeoType(GEO::HIGHERORDER);
      GEO::checkRoughGeoType(element, xyze_element, eleGeoType);
      LINALG::Matrix<3, 2> xaabbEle =
          GEO::computeFastXAABB(element->Shape(), xyze_element, eleGeoType);
      xaabb_label = mergeAABB(xaabb_label, xaabbEle);
    }
    XAABBs.push_back(xaabb_label);
  }
  return XAABBs;
}

/*----------------------------------------------------------------------*
 | returns a label for a given point                         u.may 07/08|
 | and element list                                                     |
 *----------------------------------------------------------------------*/
int GEO::getXFEMLabel(const DRT::Discretization& dis,
    const std::map<int, LINALG::Matrix<3, 1>>& currentpositions,
    const LINALG::Matrix<3, 1>& querypoint, const std::map<int, std::set<int>>& elementList)
{
  LINALG::Matrix<3, 1> minDistanceVec(true);

  GEO::NearestObject nearestObject;

  // compute the distance to the nearest object (surface, line, node) return label of nearest object
  // returns the label of the surface element structure the projection of the query point is lying
  // on
  int label = nearestObjectInNode(
      dis, currentpositions, elementList, querypoint, minDistanceVec, nearestObject);

  // compute normal in the point found on or in the object
  LINALG::Matrix<3, 1> normal = getNormalAtSurfacePoint(dis, currentpositions, nearestObject);

  // compare normals and set label
  const double scalarproduct =
      minDistanceVec(0) * normal(0) + minDistanceVec(1) * normal(1) + minDistanceVec(2) * normal(2);

  // if fluid
  if (scalarproduct > (-1) * GEO::TOL13) label = 0;
  if (minDistanceVec.Norm2() < GEO::TOL7) label = 0;

  return label;
}

/*----------------------------------------------------------------------*
 | returns a label for a given point                         u.may 07/08|
 | and element list                                                     |
 *----------------------------------------------------------------------*/
void GEO::moveNodeOutOfStructure(const std::vector<std::vector<int>>& triangleList,
    std::vector<GEO::InterfacePoint>& pointList, const int querypointId,
    const GEO::NearestObject& nearestObject, LINALG::Matrix<3, 1>& minDistanceVec)
{
  LINALG::Matrix<3, 1> querypoint = pointList[querypointId].getCoord();

  // std::cout << "query point" << querypoint << std::endl;

  // compute the distance to the nearest object (surface, line, node) return label of nearest object
  // returns the label of the surface element structure the projection of the query point is lying
  // on nearestObjectInNode(triangleList, pointList, elementList, querypoint, pointLabel,
  // minDistanceVec, nearestObject);

  std::cout << "minDistanceVec" << minDistanceVec << std::endl;

  // compute normal in the point found on or in the object
  LINALG::Matrix<3, 1> normal = getNormalAtSurfacePoint(triangleList, pointList, nearestObject);

  std::cout << "normal" << normal << std::endl;

  // compare normals and set label
  const double scalarproduct =
      minDistanceVec(0) * normal(0) + minDistanceVec(1) * normal(1) + minDistanceVec(2) * normal(2);

  const double length = minDistanceVec.Norm2();
  minDistanceVec.Scale(1.0 / length);

  std::cout << "length =" << length << std::endl;
  std::cout << "scalarproduct =" << scalarproduct << std::endl;

  // minDistanceVec.Norm2() < GEO::TOL7  &&
  // !(scalarproduct > (-1)*GEO::TOL13)
  if (scalarproduct < 0.0)
  {
    std::cout << "MOVE NODE = " << querypointId << std::endl;
    std::cout << "length =" << length << std::endl;

    std::cout << "scalarproduct =" << scalarproduct << std::endl;
    minDistanceVec.Scale((1.0 + 2 * GEO::TOL7) * length);
    querypoint += minDistanceVec;
    pointList[querypointId].setCoord(querypoint);
    // check if point now inside or outside the box adjust facemarker ...
    // and interface point
  }
  // if point lies within another structure move
  // it our of that structure to enable contact
  return;
}

/*----------------------------------------------------------------------*
 | returns a label for a given point                         u.may 07/08|
 | and element list                                                     |
 *----------------------------------------------------------------------*/
int GEO::getXFEMLabelAndNearestObject(const DRT::Discretization& dis,
    const std::map<int, LINALG::Matrix<3, 1>>& currentpositions,
    const LINALG::Matrix<3, 1>& querypoint, const std::map<int, std::set<int>>& elementList,
    GEO::NearestObject& nearestObject)
{
  LINALG::Matrix<3, 1> minDistanceVec(true);

  // compute the distance to the nearest object (surface, line, node) return label of nearest object
  // returns the label of the surface element structure the projection of the query point is lying
  // on
  int label = nearestObjectInNode(
      dis, currentpositions, elementList, querypoint, minDistanceVec, nearestObject);

  // compute normal in the point found on or in the object
  LINALG::Matrix<3, 1> normal = getNormalAtSurfacePoint(dis, currentpositions, nearestObject);

  // compare normals and set label
  const double scalarproduct =
      minDistanceVec(0) * normal(0) + minDistanceVec(1) * normal(1) + minDistanceVec(2) * normal(2);

  // if fluid
  if (scalarproduct > (-1) * GEO::TOL13) label = 0;
  if (minDistanceVec.Norm2() < GEO::TOL7) label = 0;

  return label;
}

/*----------------------------------------------------------------------*
 | returns a label for a given point                         u.may 07/08|
 | and element list                                                     |
 *----------------------------------------------------------------------*/
void GEO::getPotentialObjects(const DRT::Discretization& dis,
    const std::map<int, LINALG::Matrix<3, 1>>& currentpositions,
    const std::map<int, std::set<int>>& elementList,
    const std::vector<LINALG::Matrix<3, 1>>& gaussPoints,
    std::map<int, std::map<int, GEO::NearestObject>>& potObjects, const double cutoff_radius,
    const int label, const int projectiontype)
{
  // find nearest object for each structure
  fillPotObjectsInNode(dis, currentpositions, elementList, gaussPoints, potObjects, cutoff_radius,
      label, projectiontype);
  return;
}

/*----------------------------------------------------------------------*
 | a set of nodes in a given radius                          u.may 07/08|
 | from a query point                                                   |
 *----------------------------------------------------------------------*/
std::map<int, std::set<int>> GEO::getElementsInRadius(const DRT::Discretization& dis,
    const std::map<int, LINALG::Matrix<3, 1>>& currentpositions,
    const LINALG::Matrix<3, 1>& querypoint, const double radius, const int label,
    std::map<int, std::set<int>>& elementList)
{
  std::map<int, std::set<int>> nodeList;
  std::map<int, std::set<int>> elementMap;

  // collect all nodes with different label
  for (std::map<int, std::set<int>>::const_iterator labelIter = elementList.begin();
       labelIter != elementList.end(); labelIter++)
  {
    if (label != labelIter->first)  // don't collect nodes which belong to the same label
    {
      for (std::set<int>::const_iterator eleIter = (labelIter->second).begin();
           eleIter != (labelIter->second).end(); eleIter++)
      {
        DRT::Element* element = dis.gElement(*eleIter);
        for (int i = 0; i < DRT::UTILS::getNumberOfElementCornerNodes(element->Shape()); i++)
          nodeList[labelIter->first].insert(element->NodeIds()[i]);
      }
    }
  }

  for (std::map<int, std::set<int>>::const_iterator labelIter = nodeList.begin();
       labelIter != nodeList.end(); labelIter++)
    for (std::set<int>::const_iterator nodeIter = (labelIter->second).begin();
         nodeIter != (labelIter->second).end(); nodeIter++)
    {
      double distance = GEO::LARGENUMBER;
      const DRT::Node* node = dis.gNode(*nodeIter);
      GEO::getDistanceToPoint(node, currentpositions, querypoint, distance);

      if (distance < (radius + GEO::TOL7))
      {
        for (int i = 0; i < dis.gNode(*nodeIter)->NumElement(); i++)
          elementMap[labelIter->first].insert(dis.gNode(*nodeIter)->Elements()[i]->Id());
      }
    }

  return elementMap;
}

/*----------------------------------------------------------------------*
 | a set of nodes in a given radius from a element.          u.may 07/08|
 | The radius is approximated by a box with edge                        |
 | length 2xRadius + box width                                          |
 *----------------------------------------------------------------------*/
std::map<int, std::set<int>> GEO::getElementsInRadius(
    const DRT::Discretization& dis,  // von potential::discretRCP_
    const std::map<int, LINALG::Matrix<3, 1>>& currentpositions,
    const std::map<int, LINALG::Matrix<3, 2>>& elemXAABBList, LINALG::Matrix<3, 2> elemXAABB,
    const int label, std::map<int, std::set<int>>& elementList)
{
  std::map<int, std::set<int>> nodeList;
  std::map<int, std::set<int>> elementMap;

  // collect all nodes with different label
  for (std::map<int, std::set<int>>::const_iterator labelIter = elementList.begin();
       labelIter != elementList.end(); labelIter++)
  {
    if (label != labelIter->first)  // don't collect nodes which belong to the same label
    {
      for (std::set<int>::const_iterator eleIter = (labelIter->second).begin();
           eleIter != (labelIter->second).end(); eleIter++)
      {
        DRT::Element* element = dis.gElement(*eleIter);
        for (int i = 0; i < DRT::UTILS::getNumberOfElementCornerNodes(element->Shape()); i++)
          nodeList[labelIter->first].insert(element->NodeIds()[i]);
      }
    }
  }

  for (std::map<int, std::set<int>>::const_iterator labelIter = nodeList.begin();
       labelIter != nodeList.end(); labelIter++)
    for (std::set<int>::const_iterator nodeIter = (labelIter->second).begin();
         nodeIter != (labelIter->second).end(); nodeIter++)
    {
      const DRT::Node* node = dis.gNode(*nodeIter);
      const LINALG::Matrix<3, 1> x_node = currentpositions.find(node->Id())->second;

      if (GEO::pointInTreeNode(x_node, elemXAABB))
      {
        for (int i = 0; i < dis.gNode(*nodeIter)->NumElement(); i++)
          elementMap[labelIter->first].insert(dis.gNode(*nodeIter)->Elements()[i]->Id());
      }
    }

  return elementMap;
}

/*------------------------------------------------------------------------*
 | returns vector with IDs of nodes within given radius around querypoint;|
 | the kind of ID is the same one as the one on which currentpositions is |
 | based, thus typically node LIDs or node GIDs; however, also any other  |
 | kind of ID may be used.                                                |
 |                                                             cyron 02/09|
 *------------------------------------------------------------------------*/
std::vector<int> GEO::getPointsInRadius(const std::map<int, LINALG::Matrix<3, 1>>& currentpositions,
    const LINALG::Matrix<3, 1>& querypoint, const double radius,
    std::map<int, std::set<int>>& elementList)
{
  std::vector<int> nodes;

  // get set of (unlabled) elements in elementList
  std::set<int> pointlist = elementList[-1];

  /* looping through the set of points assigned to current tree node and checking for each
   * whether it's closer to querypoint than radius*/
  for (std::set<int>::const_iterator pointIter = pointlist.begin(); pointIter != pointlist.end();
       pointIter++)
  {
    // difference vector between querypoint and current point listed in pointList
    LINALG::Matrix<3, 1> difference = querypoint;
    difference -= (currentpositions.find(*pointIter))->second;

    /* if absolute value of difference is smaller than radius we add the GID of the node represented
     * by the i-th element of currentpositions to the vector with the GIDs of nodes within a given
     * radius*/
    if (difference.Norm2() < radius) nodes.push_back(*pointIter);
  }
  return nodes;
}

/*----------------------------------------------------------------------*
 | a vector of intersection elements                         u.may 02/09|
 | for a given query element (CONTACT)                                  | |
 *----------------------------------------------------------------------*/
void GEO::searchCollisions(const std::map<int, LINALG::Matrix<3, 2>>& currentBVs,
    const LINALG::Matrix<3, 2>& queryBV, const int label,
    const std::map<int, std::set<int>>& elementList, std::set<int>& collisions)
{
  // loop over all entries of elementList (= intersection candidates) with different label
  // run over global ids
  for (std::map<int, std::set<int>>::const_iterator labelIter = elementList.begin();
       labelIter != elementList.end(); labelIter++)
  {
    if (label != labelIter->first)  // don t collect nodes which belong to the same label
    {
      for (std::set<int>::const_iterator eleIter = (labelIter->second).begin();
           eleIter != (labelIter->second).end(); eleIter++)
      {
        if (intersectionOfBVs(queryBV, currentBVs.find(*eleIter)->second))
        {
          collisions.insert(*eleIter);
        }
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | a vector of intersection elements                         u.may 02/09|
 | for a given query element (CONTACT)                                  | |
 *----------------------------------------------------------------------*/
void GEO::searchCollisions(const std::map<int, LINALG::Matrix<9, 2>>& currentKDOPs,
    const LINALG::Matrix<9, 2>& queryKDOP, const int label,
    const std::map<int, std::set<int>>& elementList, std::set<int>& contactEleIds)
{
  // loop over all entries of elementList (= intersection candidates) with different label
  // run over global ids
  for (std::map<int, std::set<int>>::const_iterator labelIter = elementList.begin();
       labelIter != elementList.end(); labelIter++)
  {
    if (label != labelIter->first)  // don t collect nodes which belong to the same label
    {
      for (std::set<int>::const_iterator eleIter = (labelIter->second).begin();
           eleIter != (labelIter->second).end(); eleIter++)
      {
        if (intersectionOfKDOPs(queryKDOP, currentKDOPs.find(*eleIter)->second))
          contactEleIds.insert(*eleIter);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | searches a nearest object in tree node                    u.may 07/08|
 | object is either a node, line or surface element                     |
 *----------------------------------------------------------------------*/
int GEO::nearestObjectInNode(const DRT::Discretization& dis,
    const std::map<int, LINALG::Matrix<3, 1>>& currentpositions,
    const std::map<int, std::set<int>>& elementList, const LINALG::Matrix<3, 1>& point,
    LINALG::Matrix<3, 1>& minDistanceVec, GEO::NearestObject& nearestObject)
{
  bool pointFound = false;
  double min_distance = GEO::LARGENUMBER;
  double distance = GEO::LARGENUMBER;
  LINALG::Matrix<3, 1> normal(true);
  LINALG::Matrix<3, 1> x_surface(true);
  std::map<int, std::set<int>> nodeList;

  // run over all surface elements
  for (std::map<int, std::set<int>>::const_iterator labelIter = elementList.begin();
       labelIter != elementList.end(); labelIter++)
    for (std::set<int>::const_iterator eleIter = (labelIter->second).begin();
         eleIter != (labelIter->second).end(); eleIter++)
    {
      // not const because otherwise no lines can be obtained
      DRT::Element* element = dis.gElement(*eleIter);
      pointFound = GEO::getDistanceToSurface(element, currentpositions, point, x_surface, distance);
      if (pointFound && distance < min_distance)
      {
        pointFound = false;
        min_distance = distance;
        nearestObject.setSurfaceObjectType(*eleIter, labelIter->first, x_surface);
      }

      // run over all line elements
      const std::vector<Teuchos::RCP<DRT::Element>> eleLines = element->Lines();
      for (int i = 0; i < element->NumLine(); i++)
      {
        pointFound =
            GEO::getDistanceToLine(eleLines[i].get(), currentpositions, point, x_surface, distance);
        if (pointFound && distance < min_distance)
        {
          pointFound = false;
          min_distance = distance;
          nearestObject.setLineObjectType(i, *eleIter, labelIter->first, x_surface);
        }
      }
      // collect nodes
      for (int i = 0; i < DRT::UTILS::getNumberOfElementCornerNodes(element->Shape()); i++)
        nodeList[labelIter->first].insert(element->NodeIds()[i]);
    }

  // run over all nodes
  for (std::map<int, std::set<int>>::const_iterator labelIter = nodeList.begin();
       labelIter != nodeList.end(); labelIter++)
    for (std::set<int>::const_iterator nodeIter = (labelIter->second).begin();
         nodeIter != (labelIter->second).end(); nodeIter++)
    {
      const DRT::Node* node = dis.gNode(*nodeIter);
      GEO::getDistanceToPoint(node, currentpositions, point, distance);
      if (distance < min_distance)
      {
        min_distance = distance;
        nearestObject.setNodeObjectType(
            *nodeIter, labelIter->first, currentpositions.find(node->Id())->second);
      }
    }

  if (nearestObject.getObjectType() == GEO::NOTYPE_OBJECT) dserror("no nearest object obtained");

  // compute distance vector pointing away from the surface element
  minDistanceVec.Update(1.0, point, -1.0, nearestObject.getPhysCoord());

  return nearestObject.getLabel();
}

/*----------------------------------------------------------------------*
 | gives the coords of the nearest point on or in an object in  tk 01/10|
 | tree node; object is either a node or a line                         |
 *----------------------------------------------------------------------*/
void GEO::nearest2DObjectInNode(const Teuchos::RCP<DRT::Discretization> dis,
    std::map<int, Teuchos::RCP<DRT::Element>>& elements,
    const std::map<int, LINALG::Matrix<3, 1>>& currentpositions,
    const std::map<int, std::set<int>>& elementList, const LINALG::Matrix<3, 1>& point,
    LINALG::Matrix<3, 1>& minDistCoords)
{
  GEO::NearestObject nearestObject;
  bool pointFound = false;
  double min_distance = GEO::LARGENUMBER;
  double distance = GEO::LARGENUMBER;
  LINALG::Matrix<3, 1> normal(true);
  LINALG::Matrix<3, 1> x_surface(true);
  std::map<int, std::set<int>> nodeList;

  // run over all line elements
  for (std::map<int, std::set<int>>::const_iterator labelIter = elementList.begin();
       labelIter != elementList.end(); labelIter++)
    for (std::set<int>::const_iterator eleIter = (labelIter->second).begin();
         eleIter != (labelIter->second).end(); eleIter++)
    {
      // not const because otherwise no lines can be obtained
      DRT::Element* element = elements[*eleIter].get();
      pointFound = GEO::getDistanceToLine(element, currentpositions, point, x_surface, distance);
      if (pointFound && distance < min_distance)
      {
        pointFound = false;
        nearestObject.setLineObjectType(1, *eleIter, labelIter->first, x_surface);
        min_distance = distance;
      }

      // collect nodes
      for (int i = 0; i < DRT::UTILS::getNumberOfElementCornerNodes(element->Shape()); i++)
        nodeList[labelIter->first].insert(element->NodeIds()[i]);
    }

  // run over all nodes
  for (std::map<int, std::set<int>>::const_iterator labelIter = nodeList.begin();
       labelIter != nodeList.end(); labelIter++)
    for (std::set<int>::const_iterator nodeIter = (labelIter->second).begin();
         nodeIter != (labelIter->second).end(); nodeIter++)
    {
      const DRT::Node* node = dis->gNode(*nodeIter);
      GEO::getDistanceToPoint(node, currentpositions, point, distance);
      if (distance < min_distance)
      {
        min_distance = distance;
        nearestObject.setNodeObjectType(
            *nodeIter, labelIter->first, currentpositions.find(node->Id())->second);
      }
    }

  if (nearestObject.getObjectType() == GEO::NOTYPE_OBJECT) dserror("no nearest object obtained");

  // save projection point
  minDistCoords = nearestObject.getPhysCoord();

  return;
}

/*----------------------------------------------------------------------*
 | gives the coords of the nearest point on or in an object in  tk 01/10|
 | tree node; object is either a node, line or surface element;         |
 | also surface id of nearest object is returned, in case of a          |
 | line or node, a random adjacent surface id is returned               |
 *----------------------------------------------------------------------*/
int GEO::nearest3DObjectInNode(const Teuchos::RCP<DRT::Discretization> dis,
    std::map<int, Teuchos::RCP<DRT::Element>>& elements,
    const std::map<int, LINALG::Matrix<3, 1>>& currentpositions,
    const std::map<int, std::set<int>>& elementList, const LINALG::Matrix<3, 1>& point,
    LINALG::Matrix<3, 1>& minDistCoords)
{
  GEO::NearestObject nearestObject;
  bool pointFound = false;
  double min_distance = GEO::LARGENUMBER;
  double distance = GEO::LARGENUMBER;
  LINALG::Matrix<3, 1> normal(true);
  LINALG::Matrix<3, 1> x_surface(true);
  std::map<int, std::set<int>> nodeList;
  int surfid = -1;

  // run over all surface elements
  for (std::map<int, std::set<int>>::const_iterator labelIter = elementList.begin();
       labelIter != elementList.end(); labelIter++)
    for (std::set<int>::const_iterator eleIter = (labelIter->second).begin();
         eleIter != (labelIter->second).end(); eleIter++)
    {
      // not const because otherwise no lines can be obtained
      DRT::Element* element = elements[*eleIter].get();
      pointFound = GEO::getDistanceToSurface(element, currentpositions, point, x_surface, distance);
      if (pointFound && distance < min_distance)
      {
        pointFound = false;
        min_distance = distance;
        nearestObject.setSurfaceObjectType(*eleIter, labelIter->first, x_surface);
        surfid = element->Id();
      }

      // run over all line elements
      const std::vector<Teuchos::RCP<DRT::Element>> eleLines = element->Lines();
      for (int i = 0; i < element->NumLine(); i++)
      {
        pointFound =
            GEO::getDistanceToLine(eleLines[i].get(), currentpositions, point, x_surface, distance);
        if (pointFound && distance < min_distance)
        {
          pointFound = false;
          min_distance = distance;
          nearestObject.setLineObjectType(i, *eleIter, labelIter->first, x_surface);
          surfid = element->Id();
        }
      }
      // collect nodes
      for (int i = 0; i < DRT::UTILS::getNumberOfElementCornerNodes(element->Shape()); i++)
        nodeList[labelIter->first].insert(element->NodeIds()[i]);
    }

  // run over all nodes
  for (std::map<int, std::set<int>>::const_iterator labelIter = nodeList.begin();
       labelIter != nodeList.end(); labelIter++)
    for (std::set<int>::const_iterator nodeIter = (labelIter->second).begin();
         nodeIter != (labelIter->second).end(); nodeIter++)
    {
      const DRT::Node* node = dis->gNode(*nodeIter);
      GEO::getDistanceToPoint(node, currentpositions, point, distance);
      if (distance < min_distance)
      {
        min_distance = distance;
        nearestObject.setNodeObjectType(
            *nodeIter, labelIter->first, currentpositions.find(node->Id())->second);
        surfid = node->Elements()[0]->Id();  // surf id of any of the adjacent elements
      }
    }

  if (nearestObject.getObjectType() == GEO::NOTYPE_OBJECT) dserror("no nearest object obtained");

  // save projection point
  minDistCoords = nearestObject.getPhysCoord();

  return surfid;
}

/*----------------------------------------------------------------------*
 | gives the coords of the nearest point on a surface        ghamm 09/13|
 | element and return type of nearest object                            |
 *----------------------------------------------------------------------*/
GEO::ObjectType GEO::nearest3DObjectOnElement(DRT::Element* surfaceelement,
    const std::map<int, LINALG::Matrix<3, 1>>& currentpositions, const LINALG::Matrix<3, 1>& point,
    LINALG::Matrix<3, 1>& minDistCoords)
{
  GEO::NearestObject nearestObject;
  bool pointFound = false;
  double min_distance = GEO::LARGENUMBER;
  double distance = GEO::LARGENUMBER;
  LINALG::Matrix<3, 1> x_surface(true);

  pointFound =
      GEO::getDistanceToSurface(surfaceelement, currentpositions, point, x_surface, distance);
  if (pointFound && distance < min_distance)
  {
    pointFound = false;
    min_distance = distance;
    nearestObject.setSurfaceObjectType(surfaceelement->Id(), -1, x_surface);
  }

  // run over all line elements
  const std::vector<Teuchos::RCP<DRT::Element>> eleLines = surfaceelement->Lines();
  for (int i = 0; i < surfaceelement->NumLine(); i++)
  {
    pointFound =
        GEO::getDistanceToLine(eleLines[i].get(), currentpositions, point, x_surface, distance);
    if (pointFound && distance < min_distance)
    {
      pointFound = false;
      min_distance = distance;
      nearestObject.setLineObjectType(i, surfaceelement->Id(), -1, x_surface);
    }
  }

  // run over all nodes
  for (std::map<int, LINALG::Matrix<3, 1>>::const_iterator nodeIter = currentpositions.begin();
       nodeIter != currentpositions.end(); nodeIter++)
  {
    LINALG::Matrix<3, 1> distance_vector;
    // vector pointing away from the node towards physCoord
    distance_vector.Update(1.0, point, -1.0, nodeIter->second);

    // absolute distance between point and node
    distance = distance_vector.Norm2();

    if (distance < min_distance)
    {
      min_distance = distance;
      nearestObject.setNodeObjectType(nodeIter->first, -1, nodeIter->second);
    }
  }

  if (nearestObject.getObjectType() == GEO::NOTYPE_OBJECT) dserror("no nearest object obtained");

  // save projection point
  minDistCoords = nearestObject.getPhysCoord();

  return nearestObject.getObjectType();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double GEO::nearestNodeInNode(const Teuchos::RCP<DRT::Discretization> dis,
    std::map<int, Teuchos::RCP<DRT::Element>>& elements,
    const std::map<int, LINALG::Matrix<3, 1>>& currentpositions, const LINALG::Matrix<3, 1>& point,
    DRT::Node& nearnode)
{
  double min_distance = GEO::LARGENUMBER;
  double distance = GEO::LARGENUMBER;
  std::set<int> nodeList;

  // run over all surface elements to collect nodes
  std::map<int, Teuchos::RCP<DRT::Element>>::const_iterator eleIter;
  for (eleIter = elements.begin(); eleIter != elements.end(); eleIter++)
  {
    // not const because otherwise no lines can be obtained
    DRT::Element* element = (eleIter->second).get();

    // collect nodes
    for (int i = 0; i < DRT::UTILS::getNumberOfElementCornerNodes(element->Shape()); i++)
      nodeList.insert(element->NodeIds()[i]);
  }

  // run over all nodes collected above
  for (std::set<int>::const_iterator nodeIter = nodeList.begin(); nodeIter != nodeList.end();
       nodeIter++)
  {
    DRT::Node* node = dis->gNode(*nodeIter);
    GEO::getDistanceToPoint(node, currentpositions, point, distance);
    if (distance < min_distance)
    {
      min_distance = distance;
      nearnode = *node;
    }
  }

  // in case no node element was given return -1.0 as distance
  if (min_distance > 1.0e12) min_distance = -1.0;

  return min_distance;
}

/*----------------------------------------------------------------------*
 | searches a nearest object in tree node                    u.may 07/08|
 | object is either a node, line or surface element                     |
 *----------------------------------------------------------------------*/
void GEO::fillPotObjectsInNode(const DRT::Discretization& dis,
    const std::map<int, LINALG::Matrix<3, 1>>& currentpositions,
    const std::map<int, std::set<int>>& elementList,
    const std::vector<LINALG::Matrix<3, 1>>& gaussPoints,
    std::map<int, std::map<int, GEO::NearestObject>>& potObjects, const double cutoff_radius,
    const int label, const int projectiontype)
{
  potObjects.clear();
  // run over all surface elements

  for (unsigned int i_gp = 0; i_gp < gaussPoints.size(); i_gp++)
  {
    std::map<int, GEO::NearestObject> potObjectsAtGaussPoint;

    for (std::map<int, std::set<int>>::const_iterator labelIter = elementList.begin();
         labelIter != elementList.end(); labelIter++)
    {
      if (labelIter->first != label)
      {
        bool pointFound = false;
        double distance = GEO::LARGENUMBER;
        double min_distance = GEO::LARGENUMBER;
        LINALG::Matrix<3, 1> normal(true);
        LINALG::Matrix<3, 1> x_surface(true);
        std::set<int> nodeList;
        GEO::NearestObject nearestObject;

        for (std::set<int>::const_iterator eleIter = (labelIter->second).begin();
             eleIter != (labelIter->second).end(); eleIter++)
        {
          // not const because otherwise no lines can be obtained
          DRT::Element* element = dis.gElement(*eleIter);
          pointFound = GEO::getDistanceToSurface(
              element, currentpositions, gaussPoints[i_gp], x_surface, distance);
          if (pointFound && distance < cutoff_radius && distance < min_distance)
          {
            pointFound = false;
            min_distance = distance;
            nearestObject.setSurfaceObjectType(*eleIter, labelIter->first, x_surface);
          }

          // run over all line elements
          const std::vector<Teuchos::RCP<DRT::Element>> eleLines = element->Lines();
          for (int i_lines = 0; i_lines < element->NumLine(); i_lines++)
          {
            pointFound = GEO::getDistanceToLine(
                eleLines[i_lines].get(), currentpositions, gaussPoints[i_gp], x_surface, distance);
            if (pointFound && distance < cutoff_radius && distance < min_distance)
            {
              pointFound = false;
              min_distance = distance;
              nearestObject.setLineObjectType(i_lines, *eleIter, labelIter->first, x_surface);
            }
          }
          // collect nodes
          for (int i_node = 0; i_node < DRT::UTILS::getNumberOfElementCornerNodes(element->Shape());
               i_node++)
            nodeList.insert(element->NodeIds()[i_node]);
        }  // loop over structure

        // run over all nodes
        for (std::set<int>::const_iterator nodeIter = nodeList.begin(); nodeIter != nodeList.end();
             nodeIter++)
        {
          const DRT::Node* node = dis.gNode(*nodeIter);
          GEO::getDistanceToPoint(node, currentpositions, gaussPoints[i_gp], distance);
          if (distance < cutoff_radius && distance < min_distance)
          {
            min_distance = distance;
            nearestObject.setNodeObjectType(
                *nodeIter, labelIter->first, currentpositions.find(node->Id())->second);
          }
        }

        if (nearestObject.getObjectType() != GEO::NOTYPE_OBJECT)
          potObjectsAtGaussPoint[labelIter->first] = nearestObject;

      }  // labelIter->first != label
    }    // loop over structures
    if (potObjectsAtGaussPoint.size() != 0) potObjects[i_gp] = potObjectsAtGaussPoint;
  }  // loop over gaussian points

  return;
}

/*----------------------------------------------------------------------*
 | searches a nearest object in tree node                    u.may 02/09|
 | object is either a node, line or tri     element                     |
 *----------------------------------------------------------------------*/
int GEO::nearestObjectInNode(const std::vector<std::vector<int>>& triangleList,
    const std::vector<GEO::InterfacePoint>& pointList,
    const std::map<int, std::set<int>>& elementList, const LINALG::Matrix<3, 1>& point,
    const int pointLabel, LINALG::Matrix<3, 1>& minDistanceVec, GEO::NearestObject& nearestObject)
{
  bool pointFound = false;
  double min_distance = GEO::LARGENUMBER;
  double distance = GEO::LARGENUMBER;
  LINALG::Matrix<3, 1> normal(true);
  LINALG::Matrix<3, 1> x_surface(true);
  std::map<int, std::set<int>> nodeList;

  if (elementList.empty() || (elementList.size() == 1 && elementList.begin()->first == pointLabel))
    dserror("elelist empty");

  // clear nearest object
  nearestObject.clear();
  minDistanceVec.Clear();

  // run over all surface elements
  for (std::map<int, std::set<int>>::const_iterator labelIter = elementList.begin();
       labelIter != elementList.end(); labelIter++)
  {
    if (labelIter->first != pointLabel)
    {
      for (std::set<int>::const_iterator eleIter = (labelIter->second).begin();
           eleIter != (labelIter->second).end(); eleIter++)
      {
        // not const because otherwise no lines can be obtained
        pointFound = GEO::getDistanceToSurface(
            triangleList[*eleIter], pointList, point, x_surface, distance);
        if (pointFound && distance < min_distance)
        {
          pointFound = false;
          min_distance = distance;
          nearestObject.setSurfaceObjectType(*eleIter, labelIter->first, x_surface);
        }

        // run over all line elements
        for (int i = 0; i < 3; i++)
        {
          std::vector<int> linePoints(2, 0);
          if (i < 2)
          {
            linePoints[0] = triangleList[*eleIter][i];
            linePoints[1] = triangleList[*eleIter][i + 1];
          }
          else
          {
            linePoints[0] = triangleList[*eleIter][i];
            linePoints[1] = triangleList[*eleIter][0];
          }
          pointFound = GEO::getDistanceToLine(linePoints, pointList, point, x_surface, distance);
          if (pointFound && distance < min_distance)
          {
            pointFound = false;
            min_distance = distance;
            nearestObject.setLineObjectType(i, *eleIter, labelIter->first, x_surface);
          }
        }
        // collect nodes
        for (int i = 0; i < DRT::UTILS::getNumberOfElementCornerNodes(DRT::Element::tri3); i++)
          nodeList[labelIter->first].insert(triangleList[*eleIter][i]);
      }  // if(labelIter->first!=pointLabel)
    }    // for element iter
  }      // for label iter
  // run over all nodes

  for (std::map<int, std::set<int>>::const_iterator labelIter = nodeList.begin();
       labelIter != nodeList.end(); labelIter++)
  {
    if (labelIter->first != pointLabel)
    {
      for (std::set<int>::const_iterator nodeIter = (labelIter->second).begin();
           nodeIter != (labelIter->second).end(); nodeIter++)
      {
        GEO::getDistanceToPoint(*nodeIter, pointList, point, distance);
        if (distance < min_distance)
        {
          min_distance = distance;
          nearestObject.setNodeObjectType(
              *nodeIter, labelIter->first, pointList[*nodeIter].getCoord());
        }
      }  // if(labelIter->first!=pointLabel)
    }
  }

  if (nearestObject.getObjectType() == GEO::NOTYPE_OBJECT) dserror("no nearest object obtained");

  std::cout << "point = " << point << std::endl;
  std::cout << "nearestObject.getPhysCoord() = " << nearestObject.getPhysCoord() << std::endl;

  // compute distance vector pointing away from the surface element
  minDistanceVec.Update(1.0, point, -1.0, nearestObject.getPhysCoord());
  std::cout << "minDistanceVec = " << minDistanceVec << std::endl;

  return nearestObject.getLabel();
}

/*----------------------------------------------------------------------*
 |  computes the normal distance from a point to a           u.may 07/08|
 |  surface element, if it exits                                        |
 *----------------------------------------------------------------------*/
bool GEO::getDistanceToSurface(const DRT::Element* surfaceElement,
    const std::map<int, LINALG::Matrix<3, 1>>& currentpositions, const LINALG::Matrix<3, 1>& point,
    LINALG::Matrix<3, 1>& x_surface_phys, double& distance)
{
  bool pointFound = false;
  double min_distance = GEO::LARGENUMBER;
  LINALG::Matrix<3, 1> distance_vector(true);
  LINALG::Matrix<2, 1> elecoord(true);  // starting value at element center

  const LINALG::SerialDenseMatrix xyze_surfaceElement(
      GEO::getCurrentNodalPositions(surfaceElement, currentpositions));
  GEO::CurrentToSurfaceElementCoordinates(
      surfaceElement->Shape(), xyze_surfaceElement, point, elecoord);

  if (GEO::checkPositionWithinElementParameterSpace(elecoord, surfaceElement->Shape()))
  {
    GEO::elementToCurrentCoordinates(
        surfaceElement->Shape(), xyze_surfaceElement, elecoord, x_surface_phys);
    // normal pointing away from the surface towards point
    distance_vector.Update(1.0, point, -1.0, x_surface_phys);
    distance = distance_vector.Norm2();
    min_distance = distance;
    pointFound = true;
  }

  // if the element is curved
  GEO::EleGeoType eleGeoType = GEO::HIGHERORDER;
  checkRoughGeoType(surfaceElement, xyze_surfaceElement, eleGeoType);

  // TODO fix check if deformed in the linear case
  if (eleGeoType == GEO::HIGHERORDER)
  {
    LINALG::SerialDenseMatrix eleCoordMatrix =
        DRT::UTILS::getEleNodeNumbering_nodes_paramspace(surfaceElement->Shape());
    for (int i = 0; i < surfaceElement->NumNode(); i++)
    {
      // use nodes as starting values
      for (int j = 0; j < 2; j++) elecoord(j) = eleCoordMatrix(j, i);

      GEO::CurrentToSurfaceElementCoordinates(
          surfaceElement->Shape(), xyze_surfaceElement, point, elecoord);

      if (GEO::checkPositionWithinElementParameterSpace(elecoord, surfaceElement->Shape()))
      {
        LINALG::Matrix<3, 1> physcoord(true);
        GEO::elementToCurrentCoordinates(
            surfaceElement->Shape(), xyze_surfaceElement, elecoord, physcoord);
        // normal pointing away from the surface towards point
        distance_vector.Update(1.0, point, -1.0, physcoord);
        distance = distance_vector.Norm2();
        if (distance < min_distance)
        {
          x_surface_phys = physcoord;
          min_distance = distance;
          pointFound = true;
        }
      }
    }
  }
  return pointFound;
}

/*----------------------------------------------------------------------*
 |  computes the normal distance from a point to a           u.may 02/09|
 |  linear triangular element, if it exits                              |
 *----------------------------------------------------------------------*/
bool GEO::getDistanceToSurface(const std::vector<int>& triElement,
    const std::vector<GEO::InterfacePoint>& pointList, const LINALG::Matrix<3, 1>& point,
    LINALG::Matrix<3, 1>& x_surface_phys, double& distance)
{
  bool pointFound = false;
  LINALG::Matrix<3, 1> distance_vector(true);
  LINALG::Matrix<2, 1> elecoord(true);  // starting value at element center

  LINALG::SerialDenseMatrix xyze_triElement(3, 3);
  for (int i = 0; i < 3; i++)
  {
    LINALG::Matrix<3, 1> node = pointList[triElement[i]].getCoord();
    for (int j = 0; j < 3; j++) xyze_triElement(i, j) = node(j);
  }
  GEO::CurrentToSurfaceElementCoordinates(DRT::Element::tri3, xyze_triElement, point, elecoord);

  if (GEO::checkPositionWithinElementParameterSpace(elecoord, DRT::Element::tri3))
  {
    GEO::elementToCurrentCoordinates(DRT::Element::tri3, xyze_triElement, elecoord, x_surface_phys);
    // normal pointing away from the surface towards point
    distance_vector.Update(1.0, point, -1.0, x_surface_phys);
    distance = distance_vector.Norm2();
    pointFound = true;
  }

  return pointFound;
}

/*----------------------------------------------------------------------*
 |  computes the normal distance from a point to a           u.may 07/08|
 |  line element, if it exits                                           |
 *----------------------------------------------------------------------*/
bool GEO::getDistanceToLine(const DRT::Element* lineElement,
    const std::map<int, LINALG::Matrix<3, 1>>& currentpositions, const LINALG::Matrix<3, 1>& point,
    LINALG::Matrix<3, 1>& x_line_phys, double& distance)
{
  bool pointFound = false;
  double min_distance = GEO::LARGENUMBER;
  LINALG::Matrix<3, 1> distance_vector(true);
  LINALG::Matrix<1, 1> elecoord(true);  // starting value at element center

  const LINALG::SerialDenseMatrix xyze_lineElement(
      GEO::getCurrentNodalPositions(lineElement, currentpositions));
  GEO::CurrentToLineElementCoordinates(lineElement->Shape(), xyze_lineElement, point, elecoord);

  if (GEO::checkPositionWithinElementParameterSpace(elecoord, lineElement->Shape()))
  {
    GEO::elementToCurrentCoordinates(lineElement->Shape(), xyze_lineElement, elecoord, x_line_phys);
    // normal pointing away from the line towards point
    distance_vector.Update(1.0, point, -1.0, x_line_phys);
    distance = distance_vector.Norm2();
    min_distance = distance;
    pointFound = true;
  }

  // if the element is curved
  GEO::EleGeoType eleGeoType = GEO::HIGHERORDER;
  checkRoughGeoType(lineElement, xyze_lineElement, eleGeoType);

  if (eleGeoType == GEO::HIGHERORDER)
  {
    LINALG::SerialDenseMatrix eleCoordMatrix =
        DRT::UTILS::getEleNodeNumbering_nodes_paramspace(lineElement->Shape());
    // use end nodes as starting values in addition
    for (int i = 0; i < 2; i++)
    {
      // use end nodes as starting values in addition
      elecoord(0) = eleCoordMatrix(0, i);

      GEO::CurrentToLineElementCoordinates(lineElement->Shape(), xyze_lineElement, point, elecoord);

      if (GEO::checkPositionWithinElementParameterSpace(elecoord, lineElement->Shape()))
      {
        LINALG::Matrix<3, 1> physcoord(true);
        GEO::elementToCurrentCoordinates(
            lineElement->Shape(), xyze_lineElement, elecoord, physcoord);
        // normal pointing away from the line towards point
        distance_vector.Update(1.0, point, -1.0, physcoord);
        distance = distance_vector.Norm2();
        if (distance < min_distance)
        {
          x_line_phys = physcoord;
          min_distance = distance;
          pointFound = true;
        }
      }
    }
  }
  return pointFound;
}

/*----------------------------------------------------------------------*
 |  computes the normal distance from a point to a           u.may 07/08|
 |  line element, if it exits                                           |
 *----------------------------------------------------------------------*/
bool GEO::getDistanceToLine(const std::vector<int>& lineElement,
    const std::vector<GEO::InterfacePoint>& pointList, const LINALG::Matrix<3, 1>& point,
    LINALG::Matrix<3, 1>& x_line_phys, double& distance)
{
  bool pointFound = false;
  LINALG::Matrix<3, 1> distance_vector(true);
  LINALG::Matrix<1, 1> elecoord(true);  // starting value at element center

  LINALG::SerialDenseMatrix xyze_lineElement(2, 3);
  for (int i = 0; i < 2; i++)
  {
    LINALG::Matrix<3, 1> node = pointList[lineElement[i]].getCoord();
    for (int j = 0; j < 3; j++) xyze_lineElement(i, j) = node(j);
  }
  GEO::CurrentToLineElementCoordinates(DRT::Element::line2, xyze_lineElement, point, elecoord);

  if (GEO::checkPositionWithinElementParameterSpace(elecoord, DRT::Element::line2))
  {
    GEO::elementToCurrentCoordinates(DRT::Element::line2, xyze_lineElement, elecoord, x_line_phys);
    // normal pointing away from the line towards point
    distance_vector.Update(1.0, point, -1.0, x_line_phys);
    distance = distance_vector.Norm2();
    pointFound = true;
  }

  return pointFound;
}

/*----------------------------------------------------------------------*
 |  computes the distance from a point to a node             u.may 07/08|
 |  of an element                                                       |
 *----------------------------------------------------------------------*/
void GEO::getDistanceToPoint(const int nodeId, const std::vector<GEO::InterfacePoint>& pointList,
    const LINALG::Matrix<3, 1>& point, double& distance)
{
  // node position in physical coordinates
  const LINALG::Matrix<3, 1> x_node = pointList[nodeId].getCoord();

  LINALG::Matrix<3, 1> distance_vector;
  // vector pointing away from the node towards physCoord
  distance_vector.Update(1.0, point, -1.0, x_node);

  // absolute distance between point and node
  distance = distance_vector.Norm2();
}

/*----------------------------------------------------------------------*
 |  computes the distance from a point to a node             u.may 07/08|
 |  of an element                                                       |
 *----------------------------------------------------------------------*/
void GEO::getDistanceToPoint(const DRT::Node* node,
    const std::map<int, LINALG::Matrix<3, 1>>& currentpositions, const LINALG::Matrix<3, 1>& point,
    double& distance)
{
  // node position in physical coordinates
  const LINALG::Matrix<3, 1> x_node = currentpositions.find(node->Id())->second;

  LINALG::Matrix<3, 1> distance_vector;
  // vector pointing away from the node towards physCoord
  distance_vector.Update(1.0, point, -1.0, x_node);

  // absolute distance between point and node
  distance = distance_vector.Norm2();
}

/*----------------------------------------------------------------------*
 |  find adjacent elements                                   peder 07/08|
 |  for a line given by two end nodes                                   |
 *----------------------------------------------------------------------*/
std::vector<int> GEO::getAdjacentSurfaceElementsToLine(
    const DRT::Node* node1, const DRT::Node* node2)
{
  std::vector<int> adjacentElements;
  const DRT::Element* const* eleSet_node1 = node1->Elements();
  const DRT::Element* const* eleSet_node2 = node2->Elements();

  for (int i = 0; i < node1->NumElement(); i++)
    for (int j = 0; j < node2->NumElement(); j++)
    {
      if (eleSet_node1[i]->Id() == eleSet_node2[j]->Id())
        adjacentElements.push_back(eleSet_node1[i]->Id());
    }

  if (adjacentElements.size() > 2)
    dserror(
        "more than two surfaces adjacent to a line - xfem coupling conditions might be not "
        "correct");

  if (adjacentElements.size() == 0) dserror("no adjacent surface elements found");

  return adjacentElements;
}

/*----------------------------------------------------------------------*
 |  find adjacent elements                                   peder 07/08|
 |  for a line given by two end nodes                                   |
 *----------------------------------------------------------------------*/
std::vector<int> GEO::getAdjacentTriElementsToLine(
    const std::vector<std::vector<int>>& triangleList, const int triangleId, const int node1,
    const int node2)
{
  std::vector<int> adjacentElements;

  const int trinode1 = triangleList[triangleId][node1];
  const int trinode2 = triangleList[triangleId][node2];

  if (triangleList.empty()) dserror("triangle list is empty");

  for (int i = 0; i < (int)triangleList.size(); i++)
  {
    int count = 0;
    for (int j = 0; j < 3; j++)
    {
      if (triangleList[i][j] == trinode1 || triangleList[i][j] == trinode2) count++;
    }
    if (count == 2) adjacentElements.push_back(i);
  }

  if (adjacentElements.size() > 2)
    dserror(
        "more than two surfaces adjacent to a line - xfem coupling conditions might be not "
        "correct");

  if (adjacentElements.size() == 0) dserror("no adjacent surface elements found");

  return adjacentElements;
}

/*----------------------------------------------------------------------*
 |  find adjacent elements                                   u.may 05/09|
 |  for a given node                                                    |
 *----------------------------------------------------------------------*/
std::vector<int> GEO::getAdjacentTriElementsToNode(
    const std::vector<std::vector<int>>& triangleList, const int node)
{
  std::vector<int> adjacentElements;

  for (int i = 0; i < (int)triangleList.size(); i++)
    for (int j = 0; j < 3; j++)
      if (triangleList[i][j] == node)
      {
        adjacentElements.push_back(i);
        break;
      }

  if (adjacentElements.size() == 0) dserror("no adjacent surface elements found");

  return adjacentElements;
}

/*----------------------------------------------------------------------*
 |  computes the normal in a given point xsi in a            u.may 07/08|
 |  surface element (or elements if point is a on node                  |
 |  or on the line )                                                    |
 *----------------------------------------------------------------------*/
LINALG::Matrix<3, 1> GEO::getNormalAtSurfacePoint(const DRT::Discretization& dis,
    const std::map<int, LINALG::Matrix<3, 1>>& currentpositions, GEO::NearestObject& nearestObject)
{
  LINALG::Matrix<3, 1> normal(true);

  switch (nearestObject.getObjectType())
  {
    case GEO::SURFACE_OBJECT:
    {
      LINALG::Matrix<2, 1> elecoord(true);
      const DRT::Element* surfaceElement = dis.gElement(nearestObject.getSurfaceId());
      const LINALG::SerialDenseMatrix xyze_surfaceElement =
          GEO::getCurrentNodalPositions(surfaceElement, currentpositions);
      GEO::CurrentToSurfaceElementCoordinates(
          surfaceElement->Shape(), xyze_surfaceElement, nearestObject.getPhysCoord(), elecoord);
      GEO::computeNormalToSurfaceElement(
          surfaceElement->Shape(), xyze_surfaceElement, elecoord, normal);
      break;
    }
    case GEO::LINE_OBJECT:
    {
      DRT::Element* surfaceElement = dis.gElement(nearestObject.getSurfaceId());
      const std::vector<Teuchos::RCP<DRT::Element>> eleLines = surfaceElement->Lines();
      const DRT::Element* lineElement = eleLines[nearestObject.getLineId()].get();
      std::vector<int> adjacentElements =
          getAdjacentSurfaceElementsToLine(lineElement->Nodes()[0], lineElement->Nodes()[1]);

      // run over all elements adjacent ot a line
      for (std::vector<int>::const_iterator eleIter = adjacentElements.begin();
           eleIter != adjacentElements.end(); ++eleIter)
      {
        const DRT::Element* surfaceElement = dis.gElement(*eleIter);
        const LINALG::SerialDenseMatrix xyze_surfaceElement(
            GEO::getCurrentNodalPositions(surfaceElement, currentpositions));
        LINALG::Matrix<2, 1> eleCoord(true);
        LINALG::Matrix<3, 1> surface_normal(true);

        GEO::CurrentToSurfaceElementCoordinates(
            surfaceElement->Shape(), xyze_surfaceElement, nearestObject.getPhysCoord(), eleCoord);
        GEO::computeNormalToSurfaceElement(
            surfaceElement->Shape(), xyze_surfaceElement, eleCoord, surface_normal);
        normal += surface_normal;
      }
      normal.Scale(1.0 / ((double)adjacentElements.size()));
      break;
    }
    case GEO::NODE_OBJECT:
    {
      const DRT::Node* node = dis.gNode(nearestObject.getNodeId());
      // run over all elements adjacent ot a node
      for (int j = 0; j < node->NumElement(); j++)
      {
        const DRT::Element* surfaceElement = node->Elements()[j];
        const LINALG::SerialDenseMatrix xyze_surfaceElement(
            GEO::getCurrentNodalPositions(surfaceElement, currentpositions));
        LINALG::Matrix<2, 1> elecoord(true);
        LINALG::Matrix<3, 1> surface_normal(true);
        GEO::CurrentToSurfaceElementCoordinates(
            surfaceElement->Shape(), xyze_surfaceElement, nearestObject.getPhysCoord(), elecoord);
        GEO::computeNormalToSurfaceElement(
            surfaceElement->Shape(), xyze_surfaceElement, elecoord, surface_normal);
        normal += surface_normal;
      }
      normal.Scale(1.0 / ((double)node->NumElement()));
      break;
    }
    default:
      dserror("object type does not exist");
  }
  return normal;
}

/*----------------------------------------------------------------------*
 |  computes the normal in a given point xsi in a            u.may 02/09|
 |  surface element (or elements if point is a on node                  |
 |  or on the line )                                                    |
 *----------------------------------------------------------------------*/
LINALG::Matrix<3, 1> GEO::getNormalAtSurfacePoint(const std::vector<std::vector<int>>& triangleList,
    const std::vector<GEO::InterfacePoint>& pointList, const GEO::NearestObject& nearestObject)
{
  LINALG::Matrix<3, 1> normal(true);

  switch (nearestObject.getObjectType())
  {
    case GEO::SURFACE_OBJECT:
    {
      LINALG::Matrix<2, 1> elecoord(true);
      LINALG::SerialDenseMatrix xyze_triElement(3, 3);
      for (int i = 0; i < 3; i++)
      {
        LINALG::Matrix<3, 1> node =
            pointList[triangleList[nearestObject.getSurfaceId()][i]].getCoord();
        for (int j = 0; j < 3; j++) xyze_triElement(i, j) = node(j);
      }
      GEO::CurrentToSurfaceElementCoordinates(
          DRT::Element::tri3, xyze_triElement, nearestObject.getPhysCoord(), elecoord);
      GEO::computeNormalToSurfaceElement(DRT::Element::tri3, xyze_triElement, elecoord, normal);
      break;
    }
    case GEO::LINE_OBJECT:
    {
      int node1 = nearestObject.getLineId();
      int node2 = nearestObject.getLineId() + 1;
      if (nearestObject.getLineId() == 2) node2 = 0;

      std::vector<int> adjacentElements =
          getAdjacentTriElementsToLine(triangleList, nearestObject.getSurfaceId(), node1, node2);
      // run over all elements adjacent ot a line
      for (std::vector<int>::const_iterator eleIter = adjacentElements.begin();
           eleIter != adjacentElements.end(); ++eleIter)
      {
        LINALG::Matrix<2, 1> elecoord(true);
        LINALG::Matrix<3, 1> tri_normal(true);
        LINALG::SerialDenseMatrix xyze_triElement(3, 3);
        for (int i = 0; i < 3; i++)
        {
          LINALG::Matrix<3, 1> node = pointList[triangleList[*eleIter][i]].getCoord();
          for (int j = 0; j < 3; j++) xyze_triElement(i, j) = node(j);
        }
        GEO::CurrentToSurfaceElementCoordinates(
            DRT::Element::tri3, xyze_triElement, nearestObject.getPhysCoord(), elecoord);
        GEO::computeNormalToSurfaceElement(
            DRT::Element::tri3, xyze_triElement, elecoord, tri_normal);

        normal += tri_normal;
      }
      normal.Scale(1.0 / ((double)adjacentElements.size()));
      break;
    }
    case GEO::NODE_OBJECT:
    {
      // run over all elements adjacent ot a node
      std::vector<int> adjacentElements =
          getAdjacentTriElementsToNode(triangleList, nearestObject.getNodeId());

      // run over all elements adjacent ot a line
      for (std::vector<int>::const_iterator eleIter = adjacentElements.begin();
           eleIter != adjacentElements.end(); ++eleIter)
      {
        LINALG::Matrix<2, 1> elecoord(true);
        LINALG::Matrix<3, 1> tri_normal(true);
        LINALG::SerialDenseMatrix xyze_triElement(3, 3);
        for (int i = 0; i < 3; i++)
        {
          LINALG::Matrix<3, 1> node = pointList[triangleList[*eleIter][i]].getCoord();
          for (int j = 0; j < 3; j++) xyze_triElement(i, j) = node(j);
        }
        GEO::CurrentToSurfaceElementCoordinates(
            DRT::Element::tri3, xyze_triElement, nearestObject.getPhysCoord(), elecoord);
        GEO::computeNormalToSurfaceElement(
            DRT::Element::tri3, xyze_triElement, elecoord, tri_normal);

        normal += tri_normal;
      }
      normal.Scale(1.0 / ((double)adjacentElements.size()));
      break;
    }
    default:
      dserror("object type does not exist");
  }
  return normal;
}

/*----------------------------------------------------------------------*
 | check s if nearest point found on nearest object          u.may 09/08|
 | lies in a tree node specified by its box                             |
 *----------------------------------------------------------------------*/
bool GEO::pointInTreeNode(const LINALG::Matrix<3, 1>& point, const LINALG::Matrix<3, 2>& nodeBox)
{
  for (int dim = 0; dim < 3; dim++)
    if ((point(dim) < (nodeBox(dim, 0) - GEO::TOL6)) ||
        (point(dim) > (nodeBox(dim, 1) + GEO::TOL6)))
      return false;

  return true;
}

/*----------------------------------------------------------------------*
 | check s if nearest point found on nearest object          u.may 09/08|
 | lies in the minimum circle that fits inside the triangle             |
 *----------------------------------------------------------------------*/
bool GEO::pointInMinCircleInTreeNode(const LINALG::Matrix<3, 1>& nearestpoint,
    const LINALG::Matrix<3, 1>& querypoint, const LINALG::Matrix<3, 2>& nodeBox,
    const bool rootNode)
{
  double minRadius = GEO::LARGENUMBER;

  if (rootNode)
    return pointInTreeNode(nearestpoint, nodeBox);
  else
  {
    // determine minimum distance querypoint to node box wall
    for (int dim = 0; dim < 3; dim++)
      for (int i = 0; i < 2; i++)
      {
        double actRadius = fabs(querypoint(dim) - nodeBox(dim, i));
        if (actRadius < minRadius) minRadius = actRadius;
      }
  }

  // distance querypoint - nearest point
  LINALG::Matrix<3, 1> distance_vector;
  distance_vector.Update(1.0, querypoint, -1.0, nearestpoint);

  // absolute distance between point and node
  if (distance_vector.Norm2() < minRadius) return true;

  return false;
}

/*----------------------------------------------------------------------*
 | merge two AABB and deliver the resulting AABB's           peder 07/08|
 *----------------------------------------------------------------------*/
LINALG::Matrix<3, 2> GEO::mergeAABB(
    const LINALG::Matrix<3, 2>& AABB1, const LINALG::Matrix<3, 2>& AABB2)
{
  LINALG::Matrix<3, 2> mergedAABB;

  for (int dim = 0; dim < 3; dim++)
  {
    mergedAABB(dim, 0) = std::min(AABB1(dim, 0), AABB2(dim, 0));
    mergedAABB(dim, 1) = std::max(AABB1(dim, 1), AABB2(dim, 1));
  }
  return mergedAABB;
}

/*----------------------------------------------------------------------*
 | check if two AABBs are in the same node box             u.may   08/08|
 *----------------------------------------------------------------------*/
bool GEO::inSameNodeBox(const LINALG::Matrix<3, 2>& AABB_old, const LINALG::Matrix<3, 2>& AABB_new,
    const LINALG::Matrix<3, 2>& nodeBox)
{
  bool inSameNode = true;
  for (int i = 0; i < 3; i++)
  {
    if (AABB_old(i, 0) < (nodeBox(i, 0) - GEO::TOL7) ||
        AABB_old(i, 1) > (nodeBox(i, 1) + GEO::TOL7) ||
        AABB_new(i, 0) < (nodeBox(i, 0) - GEO::TOL7) ||
        AABB_new(i, 1) > (nodeBox(i, 1) + GEO::TOL7))
    {
      inSameNode = false;
      break;
    }
  }
  return inSameNode;
}

/*----------------------------------------------------------------------*
 | determines the geometry type of an element                peder 07/08|
 | not for Cartesian elements cjecked because no improvements are       |
 | expected                                                             |
 *----------------------------------------------------------------------*/
void GEO::checkRoughGeoType(const DRT::Element* element,
    const LINALG::SerialDenseMatrix xyze_element, GEO::EleGeoType& eleGeoType)
{
  const int order = DRT::UTILS::getOrder(element->Shape());

  if (order == 1)
    eleGeoType = GEO::LINEAR;  // TODO check for bilinear elements in the tree they count as
                               // higerorder fix it
  else if (order == 2)
    eleGeoType = GEO::HIGHERORDER;
  else
    dserror("order of element is not correct");

  // TODO check if a higher order element is linear
  // by checking if the higher order node is on the line
}

/*----------------------------------------------------------------------*
 | determines the geometry type of an element                 u.may 09/09|
 |  -->needed for Teuchos::RCP on element                                |
 *----------------------------------------------------------------------*/
void GEO::checkRoughGeoType(const Teuchos::RCP<DRT::Element> element,
    const LINALG::SerialDenseMatrix xyze_element, GEO::EleGeoType& eleGeoType)
{
  const int order = DRT::UTILS::getOrder(element->Shape());

  if (order == 1)
    eleGeoType = GEO::LINEAR;  // TODO check for bilinear elements in the tree they count as
                               // higerorder fix it
  else if (order == 2)
    eleGeoType = GEO::HIGHERORDER;
  else
    dserror("order of element is not correct");

  // TODO check if a higher order element is linear
  // by checking if the higher order node is on the line
}

/*----------------------------------------------------------------------*
 | returns a set of intersection candidate ids             u.may   09/08|
 *----------------------------------------------------------------------*/
void GEO::getIntersectionCandidates(const std::map<int, LINALG::Matrix<3, 2>>& currentXAABBs,
    const LINALG::Matrix<3, 2>& xfemXAABB, const std::map<int, std::set<int>>& elementList,
    std::set<int>& intersectionCandidateIds)
{
  // loop over all entries of elementList (= intersection candidates)
  // run over global ids

  for (std::set<int>::const_iterator elementIter = (elementList.begin()->second).begin();
       elementIter != (elementList.begin()->second).end(); elementIter++)
  {
    if (intersectionOfXAABB<3>(currentXAABBs.find(*elementIter)->second, xfemXAABB))
      intersectionCandidateIds.insert(*elementIter);
  }

  return;
}

/*----------------------------------------------------------------------*
 | returns a set of potential elements                     u.may   10/09|
 *----------------------------------------------------------------------*/
void GEO::getPotentialElements(const std::map<int, LINALG::Matrix<3, 2>>& currentXAABBs,
    const LINALG::Matrix<3, 2>& xfemXAABB, const std::map<int, std::set<int>>& elementList,
    std::map<int, std::set<int>>& intersectionCandidateIds, const int label)
{
  // loop over all entries of elementList (= intersection candidates)
  // run over global ids

  // collect all nodes with different label
  for (std::map<int, std::set<int>>::const_iterator labelIter = elementList.begin();
       labelIter != elementList.end(); labelIter++)
  {
    if (label != labelIter->first)  // don't collect nodes which belong to the same label
    {
      for (std::set<int>::const_iterator elementIter = (labelIter->second).begin();
           elementIter != (labelIter->second).end(); elementIter++)
      {
        if (intersectionOfXAABB<3>(currentXAABBs.find(*elementIter)->second, xfemXAABB))
          intersectionCandidateIds[labelIter->first].insert(*elementIter);
      }
    }
  }
  return;
}
