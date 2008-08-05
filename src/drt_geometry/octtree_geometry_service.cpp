/*!
\file octtree_geometry_service.cpp

\brief provides geometry methods for oct tree

<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
 */
#ifdef CCADISCRET
#include "../drt_geometry/octtree_geometry_service.H"
#include "../drt_xfem/intersection_service.H"
#include "../drt_lib/drt_utils.H"


/*----------------------------------------------------------------------*
 | delivers a axis-aligned bounding box for several given    peder 07/08|
 | discretizations                                                      |
 *----------------------------------------------------------------------*/
BlitzMat3x2 GEO::getXAABBofDis(
    const DRT::Discretization& dis,
    const std::map<int,BlitzVec3>& currentpositions)
{
  BlitzMat3x2 XAABB;

  // initialize XAABB as rectangle around the first point of dis
  const int nodeid = dis.lColElement(0)->Nodes()[0]->Id();
  const BlitzVec3 pos = currentpositions.find(nodeid)->second;
  for(int dim=0; dim<3; ++dim)
  {
    XAABB(dim, 0) = pos[dim] - XFEM::TOL7;
    XAABB(dim, 1) = pos[dim] + XFEM::TOL7;
  }

  // loop over elements and merge XAABB with their eXtendedAxisAlignedBoundingBox
  for (int j=0; j< dis.NumMyColElements(); ++j) 
  {
    const DRT::Element* element = dis.lColElement(j);
    const BlitzMat xyze_element(DRT::UTILS::getCurrentNodalPositions(element,currentpositions));
    XFEM::EleGeoType eleGeoType(XFEM::HIGHERORDER);
    GEO::checkRoughGeoType(element, xyze_element, eleGeoType);
    BlitzMat3x2 xaabbEle = XFEM::computeFastXAABB(element, xyze_element, eleGeoType);
    XAABB = mergeAABB(XAABB, xaabbEle);
  }

#ifdef DEBUG5
  cout << "_XAABB=" << XAABB(0,0) << ","<< XAABB(0,1) << "," << XAABB(1,0) << "," << XAABB(1,1) << "," << XAABB(2,0) << "," << XAABB(2,1) << endl;
#endif
  return XAABB;
}



/*----------------------------------------------------------------------*
 | compute AABB s for all labeled strutcures in the          u.may 08/08|
 | element list                                                         |
 *----------------------------------------------------------------------*/
std::vector< BlitzMat3x2 > GEO::computeXAABBForLabeledStructures(
    const DRT::Discretization&              dis,
    const std::map<int,BlitzVec3>&          currentpositions,
    const std::map<int, std::set<int> >&    elementList)
{
  std::vector< BlitzMat3x2 > XAABBs;

  for (std::map<int, std::set<int> >::const_iterator labelIter = elementList.begin(); labelIter != elementList.end(); labelIter++)
  {
    BlitzMat3x2 xaabb_label;
    // initialize xaabb_label with box around first point
    const int eleId = *((labelIter->second).begin());
    const int nodeId = dis.gElement(eleId)->Nodes()[0]->Id();
    const BlitzVec3 pos = currentpositions.find(nodeId)->second;
    for(int dim=0; dim<3; ++dim)
    {
      xaabb_label(dim, 0) = pos[dim] - XFEM::TOL7;
      xaabb_label(dim, 1) = pos[dim] + XFEM::TOL7;
    }
    // run over set elements
    for (std::set<int>::const_iterator eleIter = (labelIter->second).begin(); eleIter != (labelIter->second).end(); eleIter++)
    {   
      const DRT::Element* element = dis.gElement(*eleIter);
      const BlitzMat xyze_element(DRT::UTILS::getCurrentNodalPositions(element,currentpositions));
      XFEM::EleGeoType eleGeoType(XFEM::HIGHERORDER);
      GEO::checkRoughGeoType(element, xyze_element, eleGeoType);
      BlitzMat3x2 xaabbEle = XFEM::computeFastXAABB(element, xyze_element, eleGeoType);
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
int GEO::getXFEMLabel(
    const DRT::Discretization&              dis, 		
    const std::map<int,BlitzVec3>&          currentpositions, 	
    const BlitzVec3&                        point,
    std::map<int, std::set<int> >&          elementList)  		
{

  BlitzVec2 xsi = 0.0;
  BlitzVec3 minDistanceVec = 0.0;
 
  GEO::NearestObject nearestObject;

  // compute the distance to the nearest object (surface, line, node) return label of nearest object
  int label = nearestObjectInNode(dis, currentpositions, elementList, point, xsi, minDistanceVec, nearestObject);  	
  // compute normal in the point found on or in the object 
  BlitzVec3 normal= getNormalAtXsi(dis, currentpositions, point, xsi, nearestObject);  

  // compare normals and set label 
  const double scalarproduct = minDistanceVec(0)*normal(0) + minDistanceVec(1)*normal(1) + minDistanceVec(2)*normal(2);
  // if fluid
  if(scalarproduct > (-1)*XFEM::TOL7)
    label = 0;

  return label;
}



/*----------------------------------------------------------------------*
 | searches a nearest object in tree node                    u.may 07/08|
 | object is either a node, line or surface element                     |
 *----------------------------------------------------------------------*/
int GEO::nearestObjectInNode(
    const DRT::Discretization&                  dis,  
    const std::map<int,BlitzVec3>&              currentpositions,
    const std::map<int, std::set<int> >&        elementList,
    const BlitzVec3&                            point,
    BlitzVec2&                                  xsi,
    BlitzVec3&                                  minDistanceVec,
    GEO::NearestObject&                         nearestObject)
{
  bool pointFound = false;
  double min_distance = XFEM::LARGENUMBER;
  double distance = XFEM::LARGENUMBER;
  BlitzVec3 normal = 0.0;
  BlitzVec3 x_surface = 0.0;

  BlitzVec2 xsi_ob = 0.0;
  BlitzVec3 x_ob = 0.0;

  std::map< int, std::set<int> > nodeList;

  // run over all surface elements
  for (std::map<int, set<int> >::const_iterator labelIter = elementList.begin(); labelIter != elementList.end(); labelIter++)
    for (std::set<int>::const_iterator eleIter = (labelIter->second).begin(); eleIter != (labelIter->second).end(); eleIter++)
    {    
      // not const because other wise no lines can be obtained
      DRT::Element* element = dis.gElement(*eleIter);
      pointFound = GEO::getDistanceToSurface(element, currentpositions, point, xsi_ob, x_ob, distance);

      if(pointFound && distance < min_distance)
      {
        pointFound = false;
        xsi =  xsi_ob;
        x_surface =  x_ob;
        min_distance = distance;
        nearestObject.setSurfaceObjectType(*eleIter, labelIter->first);
      }

      // run over all line elements
      const std::vector<Teuchos::RCP< DRT::Element> > eleLines = element->Lines();
      for(int i = 0; i < element->NumLine(); i++)
      {
        pointFound = GEO::getDistanceToLine(eleLines[i].get(), currentpositions, point, xsi_ob, x_ob, distance);
        if(pointFound && distance < min_distance)
        {
          pointFound = false;
          xsi =  xsi_ob;
          x_surface =  x_ob;
          min_distance = distance;
          nearestObject.setLineObjectType(*eleIter, i, labelIter->first);
        }
      }
      // collect nodes
      for(int i = 0; i < element->NumNode(); i++)
        nodeList[labelIter->first].insert(element->NodeIds()[i]);
    }

  // run over all nodes
  for (std::map<int, std::set<int> >::const_iterator labelIter = nodeList.begin(); labelIter != nodeList.end(); labelIter++)
    for (std::set<int>::const_iterator nodeIter = (labelIter->second).begin(); nodeIter != (labelIter->second).end(); nodeIter++)
    {
      const DRT::Node* node = dis.gNode(*nodeIter);
      GEO::getDistanceToPoint(node, currentpositions, point, xsi_ob, x_ob, distance);
      if (distance < min_distance)
      {
        xsi =  xsi_ob;
        x_surface =  x_ob;
        min_distance = distance;
        nearestObject.setNodeObjectType(*nodeIter, labelIter->first);
      }
    }

  // compute vector
  minDistanceVec(0) = point(0) - x_surface(0);
  minDistanceVec(1) = point(1) - x_surface(1);
  minDistanceVec(2) = point(2) - x_surface(2);

  return nearestObject.getLabel();
}



/*----------------------------------------------------------------------*
 |  computes the normal distance from a point to a           u.may 07/08|
 |  surface element, if it exits                                        |
 *----------------------------------------------------------------------*/
bool GEO::getDistanceToSurface(
    const DRT::Element*                         surfaceElement,
    const std::map<int,BlitzVec3>&              currentpositions,
    const BlitzVec3&                            point,
    BlitzVec2&                                  xsi,
    BlitzVec3&                                  x_surface_phys,
    double&                                     distance)
{
  bool pointFound = false;
  double min_distance = XFEM::LARGENUMBER;
  BlitzVec3 distance_vector = 0.0;
  BlitzVec2 xsi_start = 0.0; // starting value at element center
  xsi = 0.0;   

  const BlitzMat xyze_surfaceElement(DRT::UTILS::getCurrentNodalPositions(surfaceElement, currentpositions));
  XFEM::CurrentToSurfaceElementCoordinates(surfaceElement, xyze_surfaceElement, x_surface_phys, xsi_start);

  if(XFEM::checkPositionWithinElementParameterSpace(xsi_start, surfaceElement->Shape()))
  { 
    // normal pointing away from the surface towards point
    distance_vector(0) = point(0) - x_surface_phys(0);
    distance_vector(1) = point(1) - x_surface_phys(1);
    distance_vector(2) = point(2) - x_surface_phys(2);
    distance = XFEM::Norm2(distance_vector);
    min_distance = distance;
    xsi = xsi_start;
    pointFound = true;
  }

  // if the element is curved
  XFEM::EleGeoType eleGeoType = XFEM::HIGHERORDER;
  checkRoughGeoType(surfaceElement, xyze_surfaceElement, eleGeoType);

  // T0DO fix check if deformed in the linear case
  if(eleGeoType == XFEM::HIGHERORDER)
  {
    for(int i = 0; i < surfaceElement->NumNode(); i++)
    { 
      // use nodes as starting values
      BlitzVec3 x_start = 0.0;
      for(int j = 0; j < 3; j++)
        x_start(j) = DRT::UTILS::getEleNodeNumbering_nodes_reference(surfaceElement->Shape())[i][j];

      XFEM::CurrentToSurfaceElementCoordinates(surfaceElement, xyze_surfaceElement, x_surface_phys, xsi_start);

      if(XFEM::checkPositionWithinElementParameterSpace(xsi_start, surfaceElement->Shape()))
      { 
        // normal pointing away from the surface towards point
        distance_vector(0) = point(0) - x_surface_phys(0);
        distance_vector(1) = point(1) - x_surface_phys(1);
        distance_vector(2) = point(2) - x_surface_phys(2);
        distance = XFEM::Norm2(distance_vector);
        min_distance = distance;
        xsi = xsi_start;
        pointFound = true;
      }
    }
  }
  return pointFound;
}



/*----------------------------------------------------------------------*
 |  computes the normal distance from a point to a           u.may 07/08|
 |  line element, if it exits                                           |
 *----------------------------------------------------------------------*/
bool GEO::getDistanceToLine(
    const DRT::Element*                         lineElement,
    const std::map<int,BlitzVec3>&              currentpositions,
    const BlitzVec3&                            point,
    BlitzVec2&                                  xsi,
    BlitzVec3&                                  x_line_phys,
    double&                                     distance)
{
  bool pointFound = false;
  double min_distance = XFEM::LARGENUMBER;
  BlitzVec3 distance_vector = 0.0;
  BlitzVec2 xsi_start = 0.0; // starting value at element center
  xsi = 0.0;   

  const BlitzMat xyze_lineElement(DRT::UTILS::getCurrentNodalPositions(lineElement, currentpositions));

  // TODO implement
  //XFEM::CurrentToLineElementCoordinates(lineElement, xyze_lineElement, x_line_phys, xsi_start);

  if(XFEM::checkPositionWithinElementParameterSpace(xsi_start, lineElement->Shape()))
  { 
    // normal pointing away from the surface towards point
    distance_vector(0) = point(0) - x_line_phys(0);
    distance_vector(1) = point(1) - x_line_phys(1);
    distance_vector(2) = point(2) - x_line_phys(2);
    distance = XFEM::Norm2(distance_vector);
    min_distance = distance;
    xsi = xsi_start;
    pointFound = true;
  }

  // if the element is curved
  XFEM::EleGeoType eleGeoType = XFEM::HIGHERORDER;
  checkRoughGeoType(lineElement, xyze_lineElement, eleGeoType);

  // T0DO fix check if deformed in the linear case
  if(eleGeoType == XFEM::HIGHERORDER)
  {
    for(int i = 0; i < lineElement->NumNode(); i++)
    { 
      // use nodes as starting values
      BlitzVec3 x_start = 0.0;
      for(int j = 0; j < 3; j++)
        x_start(j) = DRT::UTILS::getEleNodeNumbering_nodes_reference(lineElement->Shape())[i][j];

      // TODO implement
      //XFEM::CurrentToLineElementCoordinates(lineElement, xyze_lineElement, x_line_phys, xsi_start);

      if(XFEM::checkPositionWithinElementParameterSpace(xsi_start, lineElement->Shape()))
      { 
        // normal pointing away from the surface towards point
        distance_vector(0) = point(0) - x_line_phys(0);
        distance_vector(1) = point(1) - x_line_phys(1);
        distance_vector(2) = point(2) - x_line_phys(2);
        distance = XFEM::Norm2(distance_vector);
        min_distance = distance;
        xsi = xsi_start;
        pointFound = true;
      }
    }
  }
  return pointFound;
}



/*----------------------------------------------------------------------*
 |  computes the distance from a point to a node             u.may 07/08|
 |  of an element                                                       |
 *----------------------------------------------------------------------*/
void GEO::getDistanceToPoint(
    const DRT::Node*                            node,
    const std::map<int,BlitzVec3>&              currentpositions,
    const BlitzVec3&                            point,
    BlitzVec2&                                  xsi,
    BlitzVec3&                                  x_node,
    double&                                     distance)
{

  BlitzVec3 distance_vector = 0.0;
  // node position in physical coordinates

  x_node = currentpositions.find(node->Id())->second;

  // no defined element coordinates available, set to -2 -2 
  xsi(0) = -2;
  xsi(1) = -2;

  // vector pointing away from the node towards physCoord
  distance_vector(0) = point(0) - x_node(0);
  distance_vector(1) = point(1) - x_node(1);
  distance_vector(2) = point(2) - x_node(2);

  // absolute distance between point and node
  distance = XFEM::Norm2(distance_vector);
}


/*----------------------------------------------------------------------*
 |  find adjacent elements                                   peder 07/08|
 |  for a line given by two end nodes                                   |
 *----------------------------------------------------------------------*/
std::vector<int> GEO::getAdjacentSurfaceElementsToLine(
    const DRT::Node* node1,
    const DRT::Node* node2)
{
  std::vector<int> adjacentElements;
  const DRT::Element*const* eleSet_node1 = node1->Elements();
  const DRT::Element*const* eleSet_node2 = node2->Elements();

  for(int i=0; i<node1->NumElement();i++)
    for(int j=0; j<node2->NumElement();j++)
    {
      if (eleSet_node1[i]->Id() == eleSet_node2[j]->Id())
        adjacentElements.push_back((node1 -> Elements()[i])->Id());
    }

  if(adjacentElements.size() > 2)
    dserror("a surface line cannot have more than 2 adjacent surface elements");
  
  if(adjacentElements.size() == 0)
      dserror("no adjacent surface elements found");

  return adjacentElements;
}



/*----------------------------------------------------------------------*
 |  computes the normal in a given point xsi in a            u.may 07/08|
 |  surface element (or elements if point is a node                     |
 |  or on the line )                                                    |
 *----------------------------------------------------------------------*/
BlitzVec3 GEO::getNormalAtXsi(
    const DRT::Discretization&                  dis,
    const std::map<int,BlitzVec3>&              currentpositions,
    const BlitzVec3&                            point,
    const BlitzVec2&                            xsi,
    GEO::NearestObject&                         nearestObject)
{
  BlitzVec3 normal = 0.0;
  GEO::ObjectType objectType = nearestObject.getObjectType();

  // TODO not yet correct
  switch (objectType) 
  {
  case GEO::SURFACE:
  { 
    const int eleId = nearestObject.getSurfaceId();
    const DRT::Element* surfaceElement = dis.gElement(eleId);
    const BlitzMat xyze_surfaceElement = DRT::UTILS::getCurrentNodalPositions(surfaceElement, currentpositions);
    XFEM::computeNormalToSurfaceElement(surfaceElement, xyze_surfaceElement, xsi, normal);
    break;
  }
  case GEO::LINE:
  {
    const int eleId = nearestObject.getSurfaceId();
    const int lineId = nearestObject.getLineId();
    DRT::Element* surfaceElement = dis.gElement(eleId);
    const std::vector<Teuchos::RCP<DRT::Element> > eleLines = surfaceElement->Lines();
    const DRT::Element* lineElement = eleLines[lineId].get();
    const BlitzMat xyze_lineElement(DRT::UTILS::getCurrentNodalPositions(lineElement, currentpositions));
    std::vector<int> adjacentElements = getAdjacentSurfaceElementsToLine(lineElement->Nodes()[0], lineElement->Nodes()[1]);

    for(std::vector<int>::const_iterator eleIter = adjacentElements.begin(); eleIter != adjacentElements.end(); ++eleIter)
    {
      BlitzVec1 xsi_start = 0.0;
      BlitzVec2 xsi_normal = 0.0;
      BlitzVec3 tmp_normal = 0.0;
      const DRT::Element* surfaceElement = dis.gElement(*eleIter);
      const BlitzMat xyze_surfaceElement(DRT::UTILS::getCurrentNodalPositions(surfaceElement, currentpositions));

      // T0ODO finish line
      BlitzVec3 x_surface_phys = 0.0;
      XFEM::CurrentToLineElementCoordinates(lineElement, xyze_lineElement, x_surface_phys, xsi_start);
      XFEM::computeNormalToSurfaceElement(surfaceElement, xyze_surfaceElement, xsi_normal, tmp_normal);
      //TODO might be safer to scale normal ?
      normal += tmp_normal;
    }
    normal /= ((double) adjacentElements.size()); 
    break;
  }
  case GEO::NODE:
  {
    const int nodeId            = nearestObject.getNodeId();
    const DRT::Node* node       = dis.gNode(nodeId);
    const BlitzVec3 x_node      = currentpositions.find(nodeId)->second;

    for(int j=0; j<node->NumElement();j++)
    {
      BlitzVec2 xsi_start = 0.0;
      BlitzVec3 node_normal = 0.0;
      const DRT::Element* surfaceElement = node->Elements()[j];
      const BlitzMat xyze_surfaceElement(DRT::UTILS::getCurrentNodalPositions(surfaceElement, currentpositions));

      XFEM::CurrentToSurfaceElementCoordinates(surfaceElement, xyze_surfaceElement, x_node, xsi_start);
      XFEM::computeNormalToSurfaceElement(surfaceElement, xyze_surfaceElement, xsi_start, node_normal);
      normal += node_normal;
    }
    normal /= ((double) node->NumElement()); 
    break;
  }
  default:
    dserror("distance type does not exists");
  }
  return normal;
}



/*----------------------------------------------------------------------*
 | merge two AABB and deliver the resulting AABB's           peder 07/08|
 *----------------------------------------------------------------------*/
BlitzMat3x2 GEO::mergeAABB(
    const BlitzMat3x2& AABB1, 
    const BlitzMat3x2& AABB2)
{
  BlitzMat3x2 mergedAABB;

  for(int dim = 0; dim < 3; dim++)
  {
    mergedAABB(dim, 0) = std::min( AABB1(dim, 0),  AABB2(dim, 0));
    mergedAABB(dim, 1) = std::max( AABB1(dim, 1),  AABB2(dim, 1));
  }
  return mergedAABB;
}



/*----------------------------------------------------------------------*
 | check if two AABBs are in the same node box             u.may   08/08|
 *----------------------------------------------------------------------*/
bool GEO::inSameNodeBox(
    const BlitzMat3x2&  AABB_old, 
    const BlitzMat3x2&  AABB_new,
    const BlitzMat3x2&  nodeBox)
{
  bool inSameNode = true;
  for(int i = 0; i < 3; i++)
  {
    for(int j = 0; j < 2; j++)
    {
      if( AABB_old(i,j) < (nodeBox(i,0)-XFEM::TOL7) || AABB_old(i,j) > (nodeBox(i,1)+XFEM::TOL7) ||
          AABB_new(i,j) < (nodeBox(i,0)-XFEM::TOL7) || AABB_new(i,j) > (nodeBox(i,1)+XFEM::TOL7)  )
      {
        inSameNode = false;
        break;
      }
    }
    if(!inSameNode)
      break;
  }
  return inSameNode;
} 



/*----------------------------------------------------------------------*
 | determines the geometry type of an element                peder 07/08|
 | not for Cartesian elements cjecked because no improvements are       |
 | expected                                                             |
 *----------------------------------------------------------------------*/
void GEO::checkRoughGeoType(
    const DRT::Element*          element,
    const BlitzMat               xyze_element,
    XFEM::EleGeoType&            eleGeoType)
{
  const DRT::Element::DiscretizationType distype = element->Shape();

  if(DRT::UTILS::getOrder(distype) ==1)
    eleGeoType = XFEM::LINEAR;
  else if(DRT::UTILS::getOrder(distype)==2)
    eleGeoType = XFEM::HIGHERORDER;
  else
    dserror("order of element is not correct");

  // TODO maybe check if a higher order element is linear
  // by checking if the higher order node is on the line
}                             


#endif  // #ifdef CCADISCRET
