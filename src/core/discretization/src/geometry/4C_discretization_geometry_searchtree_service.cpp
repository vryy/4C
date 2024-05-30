/*----------------------------------------------------------------------*/
/*! \file

 \brief provides geometry methods for a search tree

 \level 1

 */

#include "4C_discretization_geometry_searchtree_service.hpp"

#include "4C_discretization_geometry_element_coordtrafo.hpp"
#include "4C_discretization_geometry_intersection_service.hpp"
#include "4C_discretization_geometry_intersection_service_templates.hpp"
#include "4C_discretization_geometry_position_array.hpp"
#include "4C_discretization_geometry_searchtree_nearestobject.hpp"
#include "4C_lib_discret.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 | delivers a axis-aligned bounding box for a given          peder 07/08|
 | discretization                                                       |
 *----------------------------------------------------------------------*/
CORE::LINALG::Matrix<3, 2> CORE::GEO::getXAABBofDis(const DRT::Discretization& dis,
    const std::map<int, CORE::LINALG::Matrix<3, 1>>& currentpositions)
{
  CORE::LINALG::Matrix<3, 2> XAABB(true);
  if (dis.NumGlobalElements() == 0) return XAABB;

  if (dis.NumMyColElements() == 0) return XAABB;

  // initialize XAABB as rectangle around the first point of dis
  const int nodeid = dis.lColElement(0)->Nodes()[0]->Id();
  const CORE::LINALG::Matrix<3, 1> pos = currentpositions.find(nodeid)->second;
  for (int dim = 0; dim < 3; ++dim)
  {
    XAABB(dim, 0) = pos(dim) - CORE::GEO::TOL7;
    XAABB(dim, 1) = pos(dim) + CORE::GEO::TOL7;
  }

  // TODO make more efficient by cahcking nodes
  // loop over elements and merge XAABB with their eXtendedAxisAlignedBoundingBox
  for (int j = 0; j < dis.NumMyColElements(); ++j)
  {
    const CORE::Elements::Element* element = dis.lColElement(j);
    const CORE::LINALG::SerialDenseMatrix xyze_element(
        CORE::GEO::getCurrentNodalPositions(element, currentpositions));
    CORE::GEO::EleGeoType eleGeoType(CORE::GEO::HIGHERORDER);
    CORE::GEO::checkRoughGeoType(element, xyze_element, eleGeoType);

    const CORE::LINALG::Matrix<3, 2> xaabbEle =
        CORE::GEO::computeFastXAABB(element->Shape(), xyze_element, eleGeoType);
    XAABB = mergeAABB(XAABB, xaabbEle);
  }
  return XAABB;
}

/*----------------------------------------------------------------------*
 | delivers a slightly enlarged axis-aligned bounding box for u.may09/09|
 | given elements with their current postions for sliding ALE           |
 *----------------------------------------------------------------------*/
CORE::LINALG::Matrix<3, 2> CORE::GEO::getXAABBofEles(
    std::map<int, Teuchos::RCP<CORE::Elements::Element>>& elements,
    const std::map<int, CORE::LINALG::Matrix<3, 1>>& currentpositions)
{
  CORE::LINALG::Matrix<3, 2> XAABB(true);
  if (elements.begin() == elements.end()) return XAABB;

  // initialize XAABB as rectangle around the first point of the first element
  const int nodeid = elements.begin()->second->Nodes()[0]->Id();
  const CORE::LINALG::Matrix<3, 1> pos = currentpositions.find(nodeid)->second;
  for (int dim = 0; dim < 3; ++dim)
  {
    XAABB(dim, 0) = pos(dim) - CORE::GEO::TOL7;
    XAABB(dim, 1) = pos(dim) + CORE::GEO::TOL7;
  }

  std::map<int, Teuchos::RCP<CORE::Elements::Element>>::const_iterator elemiter;
  for (elemiter = elements.begin(); elemiter != elements.end(); ++elemiter)
  {
    Teuchos::RCP<CORE::Elements::Element> currelement = elemiter->second;
    const CORE::LINALG::SerialDenseMatrix xyze_element(
        CORE::GEO::getCurrentNodalPositions(currelement, currentpositions));
    CORE::GEO::EleGeoType eleGeoType(CORE::GEO::HIGHERORDER);
    CORE::GEO::checkRoughGeoType(currelement, xyze_element, eleGeoType);
    const CORE::LINALG::Matrix<3, 2> xaabbEle =
        CORE::GEO::computeFastXAABB(currelement->Shape(), xyze_element, eleGeoType);
    XAABB = mergeAABB(XAABB, xaabbEle);
  }

  return XAABB;
}

/*----------------------------------------------------------------------*
 | delivers a axis-aligned bounding box for a given          peder 07/08|
 | discretization in reference configuration                            |
 *----------------------------------------------------------------------*/
CORE::LINALG::Matrix<3, 2> CORE::GEO::getXAABBofDis(const DRT::Discretization& dis)
{
  std::map<int, CORE::LINALG::Matrix<3, 1>> currentpositions;

  for (int lid = 0; lid < dis.NumMyColNodes(); ++lid)
  {
    const CORE::Nodes::Node* node = dis.lColNode(lid);
    CORE::LINALG::Matrix<3, 1> currpos;
    currpos(0) = node->X()[0];
    currpos(1) = node->X()[1];
    currpos(2) = node->X()[2];
    currentpositions[node->Id()] = currpos;
  }

  return getXAABBofDis(dis, currentpositions);
}


/*----------------------------------------------------------------------*
 | delivers an axis-aligned bounding box for coords          ghamm 11/12|
 *----------------------------------------------------------------------*/
CORE::LINALG::Matrix<3, 2> CORE::GEO::getXAABBofPositions(
    const std::map<int, CORE::LINALG::Matrix<3, 1>>& currentpositions)
{
  CORE::LINALG::Matrix<3, 2> XAABB(true);

  if (currentpositions.size() == 0) FOUR_C_THROW("map with current positions is emtpy");

  // initialize XAABB as rectangle around the first node
  const CORE::LINALG::Matrix<3, 1> initcurrentpos = currentpositions.begin()->second;
  for (int dim = 0; dim < 3; ++dim)
  {
    XAABB(dim, 0) = initcurrentpos(dim) - CORE::GEO::TOL7;
    XAABB(dim, 1) = initcurrentpos(dim) + CORE::GEO::TOL7;
  }

  // loop over remaining entries and merge XAABB with their eXtendedAxisAlignedBoundingBox
  std::map<int, CORE::LINALG::Matrix<3, 1>>::const_iterator iter;
  for (iter = currentpositions.begin(); iter != currentpositions.end(); ++iter)
  {
    const CORE::LINALG::Matrix<3, 1> currentpos = iter->second;
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
std::vector<CORE::LINALG::Matrix<3, 2>> CORE::GEO::computeXAABBForLabeledStructures(
    const DRT::Discretization& dis,
    const std::map<int, CORE::LINALG::Matrix<3, 1>>& currentpositions,
    const std::map<int, std::set<int>>& elementList)
{
  std::vector<CORE::LINALG::Matrix<3, 2>> XAABBs;

  for (std::map<int, std::set<int>>::const_iterator labelIter = elementList.begin();
       labelIter != elementList.end(); labelIter++)
  {
    CORE::LINALG::Matrix<3, 2> xaabb_label;
    // initialize xaabb_label with box around first point
    const int eleId = *((labelIter->second).begin());
    const int nodeId = dis.gElement(eleId)->Nodes()[0]->Id();
    const CORE::LINALG::Matrix<3, 1> pos = currentpositions.find(nodeId)->second;
    for (int dim = 0; dim < 3; ++dim)
    {
      xaabb_label(dim, 0) = pos(dim) - CORE::GEO::TOL7;
      xaabb_label(dim, 1) = pos(dim) + CORE::GEO::TOL7;
    }
    // run over set elements
    for (std::set<int>::const_iterator eleIter = (labelIter->second).begin();
         eleIter != (labelIter->second).end(); eleIter++)
    {
      const CORE::Elements::Element* element = dis.gElement(*eleIter);
      const CORE::LINALG::SerialDenseMatrix xyze_element(
          CORE::GEO::getCurrentNodalPositions(element, currentpositions));
      CORE::GEO::EleGeoType eleGeoType(CORE::GEO::HIGHERORDER);
      CORE::GEO::checkRoughGeoType(element, xyze_element, eleGeoType);
      CORE::LINALG::Matrix<3, 2> xaabbEle =
          CORE::GEO::computeFastXAABB(element->Shape(), xyze_element, eleGeoType);
      xaabb_label = mergeAABB(xaabb_label, xaabbEle);
    }
    XAABBs.push_back(xaabb_label);
  }
  return XAABBs;
}

/*----------------------------------------------------------------------*
 | a set of nodes in a given radius                          u.may 07/08|
 | from a query point                                                   |
 *----------------------------------------------------------------------*/
std::map<int, std::set<int>> CORE::GEO::getElementsInRadius(const DRT::Discretization& dis,
    const std::map<int, CORE::LINALG::Matrix<3, 1>>& currentpositions,
    const CORE::LINALG::Matrix<3, 1>& querypoint, const double radius, const int label,
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
        CORE::Elements::Element* element = dis.gElement(*eleIter);
        for (int i = 0; i < CORE::FE::getNumberOfElementCornerNodes(element->Shape()); i++)
          nodeList[labelIter->first].insert(element->NodeIds()[i]);
      }
    }
  }

  for (std::map<int, std::set<int>>::const_iterator labelIter = nodeList.begin();
       labelIter != nodeList.end(); labelIter++)
    for (std::set<int>::const_iterator nodeIter = (labelIter->second).begin();
         nodeIter != (labelIter->second).end(); nodeIter++)
    {
      double distance = CORE::GEO::LARGENUMBER;
      const CORE::Nodes::Node* node = dis.gNode(*nodeIter);
      CORE::GEO::getDistanceToPoint(node, currentpositions, querypoint, distance);

      if (distance < (radius + CORE::GEO::TOL7))
      {
        for (int i = 0; i < dis.gNode(*nodeIter)->NumElement(); i++)
          elementMap[labelIter->first].insert(dis.gNode(*nodeIter)->Elements()[i]->Id());
      }
    }

  return elementMap;
}

/*----------------------------------------------------------------------*
 | a vector of intersection elements                         u.may 02/09|
 | for a given query element (CONTACT)                                  | |
 *----------------------------------------------------------------------*/
void CORE::GEO::searchCollisions(const std::map<int, CORE::LINALG::Matrix<3, 2>>& currentBVs,
    const CORE::LINALG::Matrix<3, 2>& queryBV, const int label,
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
void CORE::GEO::searchCollisions(const std::map<int, CORE::LINALG::Matrix<9, 2>>& currentKDOPs,
    const CORE::LINALG::Matrix<9, 2>& queryKDOP, const int label,
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
 | gives the coords of the nearest point on or in an object in  tk 01/10|
 | tree node; object is either a node or a line                         |
 *----------------------------------------------------------------------*/
void CORE::GEO::nearest2DObjectInNode(const Teuchos::RCP<DRT::Discretization> dis,
    std::map<int, Teuchos::RCP<CORE::Elements::Element>>& elements,
    const std::map<int, CORE::LINALG::Matrix<3, 1>>& currentpositions,
    const std::map<int, std::set<int>>& elementList, const CORE::LINALG::Matrix<3, 1>& point,
    CORE::LINALG::Matrix<3, 1>& minDistCoords)
{
  CORE::GEO::NearestObject nearestObject;
  bool pointFound = false;
  double min_distance = CORE::GEO::LARGENUMBER;
  double distance = CORE::GEO::LARGENUMBER;
  CORE::LINALG::Matrix<3, 1> normal(true);
  CORE::LINALG::Matrix<3, 1> x_surface(true);
  std::map<int, std::set<int>> nodeList;

  // run over all line elements
  for (std::map<int, std::set<int>>::const_iterator labelIter = elementList.begin();
       labelIter != elementList.end(); labelIter++)
    for (std::set<int>::const_iterator eleIter = (labelIter->second).begin();
         eleIter != (labelIter->second).end(); eleIter++)
    {
      // not const because otherwise no lines can be obtained
      CORE::Elements::Element* element = elements[*eleIter].get();
      pointFound =
          CORE::GEO::getDistanceToLine(element, currentpositions, point, x_surface, distance);
      if (pointFound && distance < min_distance)
      {
        pointFound = false;
        nearestObject.setLineObjectType(1, *eleIter, labelIter->first, x_surface);
        min_distance = distance;
      }

      // collect nodes
      for (int i = 0; i < CORE::FE::getNumberOfElementCornerNodes(element->Shape()); i++)
        nodeList[labelIter->first].insert(element->NodeIds()[i]);
    }

  // run over all nodes
  for (std::map<int, std::set<int>>::const_iterator labelIter = nodeList.begin();
       labelIter != nodeList.end(); labelIter++)
    for (std::set<int>::const_iterator nodeIter = (labelIter->second).begin();
         nodeIter != (labelIter->second).end(); nodeIter++)
    {
      const CORE::Nodes::Node* node = dis->gNode(*nodeIter);
      CORE::GEO::getDistanceToPoint(node, currentpositions, point, distance);
      if (distance < min_distance)
      {
        min_distance = distance;
        nearestObject.setNodeObjectType(
            *nodeIter, labelIter->first, currentpositions.find(node->Id())->second);
      }
    }

  if (nearestObject.getObjectType() == CORE::GEO::NOTYPE_OBJECT)
    FOUR_C_THROW("no nearest object obtained");

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
int CORE::GEO::nearest3DObjectInNode(const Teuchos::RCP<DRT::Discretization> dis,
    std::map<int, Teuchos::RCP<CORE::Elements::Element>>& elements,
    const std::map<int, CORE::LINALG::Matrix<3, 1>>& currentpositions,
    const std::map<int, std::set<int>>& elementList, const CORE::LINALG::Matrix<3, 1>& point,
    CORE::LINALG::Matrix<3, 1>& minDistCoords)
{
  CORE::GEO::NearestObject nearestObject;
  bool pointFound = false;
  double min_distance = CORE::GEO::LARGENUMBER;
  double distance = CORE::GEO::LARGENUMBER;
  CORE::LINALG::Matrix<3, 1> normal(true);
  CORE::LINALG::Matrix<3, 1> x_surface(true);
  std::map<int, std::set<int>> nodeList;
  int surfid = -1;

  // run over all surface elements
  for (std::map<int, std::set<int>>::const_iterator labelIter = elementList.begin();
       labelIter != elementList.end(); labelIter++)
    for (std::set<int>::const_iterator eleIter = (labelIter->second).begin();
         eleIter != (labelIter->second).end(); eleIter++)
    {
      // not const because otherwise no lines can be obtained
      CORE::Elements::Element* element = elements[*eleIter].get();
      pointFound =
          CORE::GEO::getDistanceToSurface(element, currentpositions, point, x_surface, distance);
      if (pointFound && distance < min_distance)
      {
        pointFound = false;
        min_distance = distance;
        nearestObject.set_surface_object_type(*eleIter, labelIter->first, x_surface);
        surfid = element->Id();
      }

      // run over all line elements
      const std::vector<Teuchos::RCP<CORE::Elements::Element>> eleLines = element->Lines();
      for (int i = 0; i < element->NumLine(); i++)
      {
        pointFound = CORE::GEO::getDistanceToLine(
            eleLines[i].get(), currentpositions, point, x_surface, distance);
        if (pointFound && distance < min_distance)
        {
          pointFound = false;
          min_distance = distance;
          nearestObject.setLineObjectType(i, *eleIter, labelIter->first, x_surface);
          surfid = element->Id();
        }
      }
      // collect nodes
      for (int i = 0; i < CORE::FE::getNumberOfElementCornerNodes(element->Shape()); i++)
        nodeList[labelIter->first].insert(element->NodeIds()[i]);
    }

  // run over all nodes
  for (std::map<int, std::set<int>>::const_iterator labelIter = nodeList.begin();
       labelIter != nodeList.end(); labelIter++)
    for (std::set<int>::const_iterator nodeIter = (labelIter->second).begin();
         nodeIter != (labelIter->second).end(); nodeIter++)
    {
      const CORE::Nodes::Node* node = dis->gNode(*nodeIter);
      CORE::GEO::getDistanceToPoint(node, currentpositions, point, distance);
      if (distance < min_distance)
      {
        min_distance = distance;
        nearestObject.setNodeObjectType(
            *nodeIter, labelIter->first, currentpositions.find(node->Id())->second);
        surfid = node->Elements()[0]->Id();  // surf id of any of the adjacent elements
      }
    }

  if (nearestObject.getObjectType() == CORE::GEO::NOTYPE_OBJECT)
    FOUR_C_THROW("no nearest object obtained");

  // save projection point
  minDistCoords = nearestObject.getPhysCoord();

  return surfid;
}

/*----------------------------------------------------------------------*
 | gives the coords of the nearest point on a surface        ghamm 09/13|
 | element and return type of nearest object                            |
 *----------------------------------------------------------------------*/
CORE::GEO::ObjectType CORE::GEO::nearest3DObjectOnElement(CORE::Elements::Element* surfaceelement,
    const std::map<int, CORE::LINALG::Matrix<3, 1>>& currentpositions,
    const CORE::LINALG::Matrix<3, 1>& point, CORE::LINALG::Matrix<3, 1>& minDistCoords)
{
  CORE::GEO::NearestObject nearestObject;
  bool pointFound = false;
  double min_distance = CORE::GEO::LARGENUMBER;
  double distance = CORE::GEO::LARGENUMBER;
  CORE::LINALG::Matrix<3, 1> x_surface(true);

  pointFound =
      CORE::GEO::getDistanceToSurface(surfaceelement, currentpositions, point, x_surface, distance);
  if (pointFound && distance <= min_distance)
  {
    pointFound = false;
    min_distance = distance;
    nearestObject.set_surface_object_type(surfaceelement->Id(), -1, x_surface);
  }

  // run over all line elements
  const std::vector<Teuchos::RCP<CORE::Elements::Element>> eleLines = surfaceelement->Lines();
  for (int i = 0; i < surfaceelement->NumLine(); i++)
  {
    pointFound = CORE::GEO::getDistanceToLine(
        eleLines[i].get(), currentpositions, point, x_surface, distance);
    if (pointFound && distance <= min_distance)
    {
      pointFound = false;
      min_distance = distance;
      nearestObject.setLineObjectType(i, surfaceelement->Id(), -1, x_surface);
    }
  }

  // run over all nodes
  for (std::map<int, CORE::LINALG::Matrix<3, 1>>::const_iterator nodeIter =
           currentpositions.begin();
       nodeIter != currentpositions.end(); nodeIter++)
  {
    CORE::LINALG::Matrix<3, 1> distance_vector;
    // vector pointing away from the node towards physCoord
    distance_vector.Update(1.0, point, -1.0, nodeIter->second);

    // absolute distance between point and node
    distance = distance_vector.Norm2();

    if (distance <= min_distance)
    {
      min_distance = distance;
      nearestObject.setNodeObjectType(nodeIter->first, -1, nodeIter->second);
    }
  }

  if (nearestObject.getObjectType() == CORE::GEO::NOTYPE_OBJECT)
    FOUR_C_THROW("no nearest object obtained");

  // save projection point
  minDistCoords = nearestObject.getPhysCoord();

  return nearestObject.getObjectType();
}

/*----------------------------------------------------------------------*
 |  computes the normal distance from a point to a           u.may 07/08|
 |  surface element, if it exits                                        |
 *----------------------------------------------------------------------*/
bool CORE::GEO::getDistanceToSurface(const CORE::Elements::Element* surfaceElement,
    const std::map<int, CORE::LINALG::Matrix<3, 1>>& currentpositions,
    const CORE::LINALG::Matrix<3, 1>& point, CORE::LINALG::Matrix<3, 1>& x_surface_phys,
    double& distance)
{
  bool pointFound = false;
  double min_distance = CORE::GEO::LARGENUMBER;
  CORE::LINALG::Matrix<3, 1> distance_vector(true);
  CORE::LINALG::Matrix<2, 1> elecoord(true);  // starting value at element center

  const CORE::LINALG::SerialDenseMatrix xyze_surfaceElement(
      CORE::GEO::getCurrentNodalPositions(surfaceElement, currentpositions));
  CORE::GEO::CurrentToSurfaceElementCoordinates(
      surfaceElement->Shape(), xyze_surfaceElement, point, elecoord);

  if (CORE::GEO::checkPositionWithinElementParameterSpace(elecoord, surfaceElement->Shape()))
  {
    CORE::GEO::elementToCurrentCoordinates(
        surfaceElement->Shape(), xyze_surfaceElement, elecoord, x_surface_phys);
    // normal pointing away from the surface towards point
    distance_vector.Update(1.0, point, -1.0, x_surface_phys);
    distance = distance_vector.Norm2();
    min_distance = distance;
    pointFound = true;
  }

  // if the element is curved
  CORE::GEO::EleGeoType eleGeoType = CORE::GEO::HIGHERORDER;
  checkRoughGeoType(surfaceElement, xyze_surfaceElement, eleGeoType);

  // TODO fix check if deformed in the linear case
  if (eleGeoType == CORE::GEO::HIGHERORDER)
  {
    CORE::LINALG::SerialDenseMatrix eleCoordMatrix =
        CORE::FE::getEleNodeNumbering_nodes_paramspace(surfaceElement->Shape());
    for (int i = 0; i < surfaceElement->num_node(); i++)
    {
      // use nodes as starting values
      for (int j = 0; j < 2; j++) elecoord(j) = eleCoordMatrix(j, i);

      CORE::GEO::CurrentToSurfaceElementCoordinates(
          surfaceElement->Shape(), xyze_surfaceElement, point, elecoord);

      if (CORE::GEO::checkPositionWithinElementParameterSpace(elecoord, surfaceElement->Shape()))
      {
        CORE::LINALG::Matrix<3, 1> physcoord(true);
        CORE::GEO::elementToCurrentCoordinates(
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
 |  computes the normal distance from a point to a           u.may 07/08|
 |  line element, if it exits                                           |
 *----------------------------------------------------------------------*/
bool CORE::GEO::getDistanceToLine(const CORE::Elements::Element* lineElement,
    const std::map<int, CORE::LINALG::Matrix<3, 1>>& currentpositions,
    const CORE::LINALG::Matrix<3, 1>& point, CORE::LINALG::Matrix<3, 1>& x_line_phys,
    double& distance)
{
  bool pointFound = false;
  double min_distance = CORE::GEO::LARGENUMBER;
  CORE::LINALG::Matrix<3, 1> distance_vector(true);
  CORE::LINALG::Matrix<1, 1> elecoord(true);  // starting value at element center

  const CORE::LINALG::SerialDenseMatrix xyze_lineElement(
      CORE::GEO::getCurrentNodalPositions(lineElement, currentpositions));
  CORE::GEO::CurrentToLineElementCoordinates(
      lineElement->Shape(), xyze_lineElement, point, elecoord);

  if (CORE::GEO::checkPositionWithinElementParameterSpace(elecoord, lineElement->Shape()))
  {
    CORE::GEO::elementToCurrentCoordinates(
        lineElement->Shape(), xyze_lineElement, elecoord, x_line_phys);
    // normal pointing away from the line towards point
    distance_vector.Update(1.0, point, -1.0, x_line_phys);
    distance = distance_vector.Norm2();
    min_distance = distance;
    pointFound = true;
  }

  // if the element is curved
  CORE::GEO::EleGeoType eleGeoType = CORE::GEO::HIGHERORDER;
  checkRoughGeoType(lineElement, xyze_lineElement, eleGeoType);

  if (eleGeoType == CORE::GEO::HIGHERORDER)
  {
    CORE::LINALG::SerialDenseMatrix eleCoordMatrix =
        CORE::FE::getEleNodeNumbering_nodes_paramspace(lineElement->Shape());
    // use end nodes as starting values in addition
    for (int i = 0; i < 2; i++)
    {
      // use end nodes as starting values in addition
      elecoord(0) = eleCoordMatrix(0, i);

      CORE::GEO::CurrentToLineElementCoordinates(
          lineElement->Shape(), xyze_lineElement, point, elecoord);

      if (CORE::GEO::checkPositionWithinElementParameterSpace(elecoord, lineElement->Shape()))
      {
        CORE::LINALG::Matrix<3, 1> physcoord(true);
        CORE::GEO::elementToCurrentCoordinates(
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
 |  computes the distance from a point to a node             u.may 07/08|
 |  of an element                                                       |
 *----------------------------------------------------------------------*/
void CORE::GEO::getDistanceToPoint(const CORE::Nodes::Node* node,
    const std::map<int, CORE::LINALG::Matrix<3, 1>>& currentpositions,
    const CORE::LINALG::Matrix<3, 1>& point, double& distance)
{
  // node position in physical coordinates
  const CORE::LINALG::Matrix<3, 1> x_node = currentpositions.find(node->Id())->second;

  CORE::LINALG::Matrix<3, 1> distance_vector;
  // vector pointing away from the node towards physCoord
  distance_vector.Update(1.0, point, -1.0, x_node);

  // absolute distance between point and node
  distance = distance_vector.Norm2();
}

/*----------------------------------------------------------------------*
 | check s if nearest point found on nearest object          u.may 09/08|
 | lies in a tree node specified by its box                             |
 *----------------------------------------------------------------------*/
bool CORE::GEO::pointInTreeNode(
    const CORE::LINALG::Matrix<3, 1>& point, const CORE::LINALG::Matrix<3, 2>& nodeBox)
{
  for (int dim = 0; dim < 3; dim++)
    if ((point(dim) < (nodeBox(dim, 0) - CORE::GEO::TOL6)) ||
        (point(dim) > (nodeBox(dim, 1) + CORE::GEO::TOL6)))
      return false;

  return true;
}

/*----------------------------------------------------------------------*
 | merge two AABB and deliver the resulting AABB's           peder 07/08|
 *----------------------------------------------------------------------*/
CORE::LINALG::Matrix<3, 2> CORE::GEO::mergeAABB(
    const CORE::LINALG::Matrix<3, 2>& AABB1, const CORE::LINALG::Matrix<3, 2>& AABB2)
{
  CORE::LINALG::Matrix<3, 2> mergedAABB;

  for (int dim = 0; dim < 3; dim++)
  {
    mergedAABB(dim, 0) = std::min(AABB1(dim, 0), AABB2(dim, 0));
    mergedAABB(dim, 1) = std::max(AABB1(dim, 1), AABB2(dim, 1));
  }
  return mergedAABB;
}

/*----------------------------------------------------------------------*
 | determines the geometry type of an element                peder 07/08|
 | not for Cartesian elements cjecked because no improvements are       |
 | expected                                                             |
 *----------------------------------------------------------------------*/
void CORE::GEO::checkRoughGeoType(const CORE::Elements::Element* element,
    const CORE::LINALG::SerialDenseMatrix xyze_element, CORE::GEO::EleGeoType& eleGeoType)
{
  const int order = CORE::FE::getOrder(element->Shape());

  if (order == 1)
    eleGeoType = CORE::GEO::LINEAR;  // TODO check for bilinear elements in the tree they count
                                     // as higerorder fix it
  else if (order == 2)
    eleGeoType = CORE::GEO::HIGHERORDER;
  else
    FOUR_C_THROW("order of element is not correct");

  // TODO check if a higher order element is linear
  // by checking if the higher order node is on the line
}

/*----------------------------------------------------------------------*
 | determines the geometry type of an element                 u.may 09/09|
 |  -->needed for Teuchos::RCP on element                                |
 *----------------------------------------------------------------------*/
void CORE::GEO::checkRoughGeoType(const Teuchos::RCP<CORE::Elements::Element> element,
    const CORE::LINALG::SerialDenseMatrix xyze_element, CORE::GEO::EleGeoType& eleGeoType)
{
  const int order = CORE::FE::getOrder(element->Shape());

  if (order == 1)
    eleGeoType = CORE::GEO::LINEAR;  // TODO check for bilinear elements in the tree they count
                                     // as higerorder fix it
  else if (order == 2)
    eleGeoType = CORE::GEO::HIGHERORDER;
  else
    FOUR_C_THROW("order of element is not correct");

  // TODO check if a higher order element is linear
  // by checking if the higher order node is on the line
}

FOUR_C_NAMESPACE_CLOSE
