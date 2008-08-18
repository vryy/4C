/*!
\file searchtree_geometry_service.cpp

\brief provides geometry methods for a search tree

<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
 */
#ifdef CCADISCRET
#include "searchtree_geometry_service.H"
#include "intersection_service.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_contact/drt_celement.H"


/*----------------------------------------------------------------------*
 | delivers a axis-aligned bounding box for a given          peder 07/08|
 | discretization                                                       |
 *----------------------------------------------------------------------*/
BlitzMat3x2 GEO::getXAABBofDis(
    const DRT::Discretization& dis,
    const std::map<int,BlitzVec3>& currentpositions)
{
  BlitzMat3x2 XAABB;
  if (dis.NumGlobalElements() == 0){
    XAABB(0,0)=0;XAABB(0,1)=0;
    XAABB(1,0)=0;XAABB(1,1)=0;
    XAABB(2,0)=0;XAABB(2,1)=0;
    return XAABB;
    }

  // initialize XAABB as rectangle around the first point of dis
  const int nodeid = dis.lColElement(0)->Nodes()[0]->Id();
  const BlitzVec3 pos = currentpositions.find(nodeid)->second;
  for(int dim=0; dim<3; ++dim)
  {
    XAABB(dim, 0) = pos[dim] - GEO::TOL7;
    XAABB(dim, 1) = pos[dim] + GEO::TOL7;
  }

  // loop over elements and merge XAABB with their eXtendedAxisAlignedBoundingBox
  for (int j=0; j< dis.NumMyColElements(); ++j) 
  {
    const DRT::Element* element = dis.lColElement(j);
    const BlitzMat xyze_element(DRT::UTILS::getCurrentNodalPositions(element,currentpositions));
    GEO::EleGeoType eleGeoType(GEO::HIGHERORDER);
    GEO::checkRoughGeoType(element, xyze_element, eleGeoType);
    BlitzMat3x2 xaabbEle = GEO::computeFastXAABB(element, xyze_element, eleGeoType);
    XAABB = mergeAABB(XAABB, xaabbEle);
  }
  
  return XAABB;
}



/*----------------------------------------------------------------------*
 | delivers a axis-aligned bounding box for a given          peder 07/08|
 | discretization in reference configuration                            |
 *----------------------------------------------------------------------*/
BlitzMat3x2 GEO::getXAABBofDis(
    const DRT::Discretization& dis)
{
  std::map<int,BlitzVec3> currentpositions;

  for (int lid = 0; lid < dis.NumMyColNodes(); ++lid)
  {
    const DRT::Node* node = dis.lColNode(lid);
    BlitzVec3 currpos;
    currpos(0) = node->X()[0];
    currpos(1) = node->X()[1];
    currpos(2) = node->X()[2];
    currentpositions[node->Id()] = currpos;
  }
  
  return getXAABBofDis(dis, currentpositions);
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
      xaabb_label(dim, 0) = pos[dim] - GEO::TOL7;
      xaabb_label(dim, 1) = pos[dim] + GEO::TOL7;
    }
    // run over set elements
    for (std::set<int>::const_iterator eleIter = (labelIter->second).begin(); eleIter != (labelIter->second).end(); eleIter++)
    {   
      const DRT::Element* element = dis.gElement(*eleIter);
      const BlitzMat xyze_element(DRT::UTILS::getCurrentNodalPositions(element,currentpositions));
      GEO::EleGeoType eleGeoType(GEO::HIGHERORDER);
      GEO::checkRoughGeoType(element, xyze_element, eleGeoType);
      BlitzMat3x2 xaabbEle = GEO::computeFastXAABB(element, xyze_element, eleGeoType);
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
    const BlitzVec3&                        querypoint,
    std::map<int, std::set<int> >&          elementList)  		
{
  BlitzVec3 minDistanceVec = 0.0;
 
  GEO::NearestObject nearestObject;

  // compute the distance to the nearest object (surface, line, node) return label of nearest object
  // returns the label of the surface element structure the projection of the query point is lying on
  int label = nearestObjectInNode(dis, currentpositions, elementList, querypoint, minDistanceVec, nearestObject);  	
 
  // compute normal in the point found on or in the object 
  BlitzVec3 normal= getNormalAtSurfacePoint(dis, currentpositions, nearestObject);  

  // compare normals and set label 
  const double scalarproduct = minDistanceVec(0)*normal(0) + minDistanceVec(1)*normal(1) + minDistanceVec(2)*normal(2);
  
  // if fluid
  if(scalarproduct > (-1)*GEO::TOL7)
    label = 0;

  return label;
}



/*----------------------------------------------------------------------*
 | a set of nodes in a given radius                          u.may 07/08|
 | from a query point                                                   |
 *----------------------------------------------------------------------*/
std::set<int> GEO::getNodeSetInRadius(
    const DRT::Discretization&              dis, 		
    const std::map<int,BlitzVec3>&          currentpositions, 	
    const BlitzVec3&                        querypoint,
    const double		                        radius,
    const int                               label,
    std::map<int, std::set<int> >&          elementList)  
{
  std::map< int, std::set<int> > nodeList;
  std::set<int> nodeSet;
  
  // collect all nodes with different label
  for(std::map<int, std::set<int> >::const_iterator labelIter = elementList.begin(); labelIter != elementList.end(); labelIter++)
  {
    if(label != labelIter->first) // don t collect nodes which belong to the same label
    {
      for(std::set<int>::const_iterator eleIter = (labelIter->second).begin(); eleIter != (labelIter->second).end(); eleIter++)
      {
        DRT::Element* element = dis.gElement(*eleIter);
        // run over all corner points of the element
        for(int i = 0; i < 4; i++)
          nodeList[labelIter->first].insert(element->NodeIds()[i]);
      }
    }
  }
  
  for(std::map<int, std::set<int> >::const_iterator labelIter = nodeList.begin(); labelIter != nodeList.end(); labelIter++)
    for(std::set<int>::const_iterator nodeIter = (labelIter->second).begin(); nodeIter != (labelIter->second).end(); nodeIter++)
    {
      double distance = GEO::LARGENUMBER;
      const DRT::Node* node = dis.gNode(*nodeIter);
      GEO::getDistanceToPoint(node, currentpositions, querypoint, distance);
      if(distance < (radius + GEO::TOL7) )
        nodeSet.insert(node->Id());
    }
  
  return nodeSet;
}


/*----------------------------------------------------------------------*
 | a vector of intersection elements                          popp 07/08|
 | for a given query element (CONTACT)                                  |                                 |
 *----------------------------------------------------------------------*/
std::vector<int> GEO::getIntersectionElements(
    const DRT::Discretization&              dis,    
    const std::map<int,BlitzVec3>&          currentpositions,   
    DRT::Element*                           element,
    std::map<int, std::set<int> >&          elementList)  
{
  // vector containing sought-after element ids
  vector<int> intersectionids;
  
  // create XAABB for query element
  const BlitzMat xyze(DRT::UTILS::getCurrentNodalPositions(element,currentpositions));  
  BlitzMat3x2 slaveXAABB = GEO::computeContactXAABB(element,xyze);
  
  // loop over all entries of elementList (= intersection candidates)
  for(std::set<int>::const_iterator elementIter = (elementList.begin()->second).begin();
      elementIter != (elementList.begin()->second).end(); elementIter++)
  {
    DRT::Element* masterelement = dis.gElement(*elementIter);
    const BlitzMat masterxyze(DRT::UTILS::getCurrentNodalPositions(masterelement,currentpositions));
    BlitzMat3x2 masterXAABB = GEO::computeContactXAABB(masterelement,masterxyze);
    
    // compare slaveXAABB and masterXAABB (2D only)
    if(intersectionOfXAABB2D(slaveXAABB,masterXAABB))
      intersectionids.push_back(*elementIter);
  }

  return intersectionids;
}

/*----------------------------------------------------------------------*
 | checks if two XAABB's intersect (2D CONTACT)               popp 08/08|
 *----------------------------------------------------------------------*/
bool GEO::intersectionOfXAABB2D(  
    const BlitzMat3x2&     slaveXAABB, 
    const BlitzMat3x2&     masterXAABB)
{
    
  /*====================================================================*/
  /* bounding box topology*/
  /*--------------------------------------------------------------------*/
  /* parameter coordinates (x,y,z) of nodes
   * node 0: (minX, minY, minZ)
   * node 1: (maxX, minY, minZ)
   * node 2: (maxX, maxY, minZ)
   * node 3: (minX, maxY, minZ)
   * node 4: (minX, minY, maxZ)
   * node 5: (maxX, minY, maxZ)
   * node 6: (maxX, maxY, maxZ)
   * node 7: (minX, maxY, maxZ)
   * 
   *                      z
   *                      |           
   *             4========|================7
   *           //|        |               /||
   *          // |        |              //||
   *         //  |        |             // ||
   *        //   |        |            //  ||
   *       //    |        |           //   ||
   *      //     |        |          //    ||
   *     //      |        |         //     ||
   *     5=========================6       ||
   *    ||       |        |        ||      ||
   *    ||       |        o--------||---------y
   *    ||       |       /         ||      ||
   *    ||       0------/----------||------3
   *    ||      /      /           ||     //
   *    ||     /      /            ||    //
   *    ||    /      /             ||   //
   *    ||   /      /              ||  //
   *    ||  /      /               || //
   *    || /      x                ||//
   *    ||/                        ||/
   *     1=========================2
   *
   */
  /*====================================================================*/
    
    bool intersection =  false;
    static std::vector<BlitzVec3> nodes(8, BlitzVec3());
    
    nodes[0](0) = slaveXAABB(0,0); nodes[0](1) = slaveXAABB(1,0); nodes[0](2) = slaveXAABB(2,0); // node 0   
    nodes[1](0) = slaveXAABB(0,1); nodes[1](1) = slaveXAABB(1,0); nodes[1](2) = slaveXAABB(2,0); // node 1
    nodes[2](0) = slaveXAABB(0,1); nodes[2](1) = slaveXAABB(1,1); nodes[2](2) = slaveXAABB(2,0); // node 2
    nodes[3](0) = slaveXAABB(0,0); nodes[3](1) = slaveXAABB(1,1); nodes[3](2) = slaveXAABB(2,0); // node 3
    nodes[4](0) = slaveXAABB(0,0); nodes[4](1) = slaveXAABB(1,0); nodes[4](2) = slaveXAABB(2,1); // node 4
    nodes[5](0) = slaveXAABB(0,1); nodes[5](1) = slaveXAABB(1,0); nodes[5](2) = slaveXAABB(2,1); // node 5
    nodes[6](0) = slaveXAABB(0,1); nodes[6](1) = slaveXAABB(1,1); nodes[6](2) = slaveXAABB(2,1); // node 6
    nodes[7](0) = slaveXAABB(0,0); nodes[7](1) = slaveXAABB(1,1); nodes[7](2) = slaveXAABB(2,1); // node 7
    
    for (int i = 0; i < 8; i++)
        if(isPositionWithinXAABB(nodes[i], masterXAABB))
        {
            intersection = true;
            break;
        }
    
    if(!intersection)
    {
        nodes[0](0) = masterXAABB(0,0);   nodes[0](1) = masterXAABB(1,0);   nodes[0](2) = masterXAABB(2,0);   // node 0   
        nodes[1](0) = masterXAABB(0,1);   nodes[1](1) = masterXAABB(1,0);   nodes[1](2) = masterXAABB(2,0);   // node 1
        nodes[2](0) = masterXAABB(0,1);   nodes[2](1) = masterXAABB(1,1);   nodes[2](2) = masterXAABB(2,0);   // node 2
        nodes[3](0) = masterXAABB(0,0);   nodes[3](1) = masterXAABB(1,1);   nodes[3](2) = masterXAABB(2,0);   // node 3
        nodes[4](0) = masterXAABB(0,0);   nodes[4](1) = masterXAABB(1,0);   nodes[4](2) = masterXAABB(2,1);   // node 4
        nodes[5](0) = masterXAABB(0,1);   nodes[5](1) = masterXAABB(1,0);   nodes[5](2) = masterXAABB(2,1);   // node 5
        nodes[6](0) = masterXAABB(0,1);   nodes[6](1) = masterXAABB(1,1);   nodes[6](2) = masterXAABB(2,1);   // node 6
        nodes[7](0) = masterXAABB(0,0);   nodes[7](1) = masterXAABB(1,1);   nodes[7](2) = masterXAABB(2,1);   // node 7
    
        for (int i = 0; i < 8; i++)
            if(isPositionWithinXAABB(nodes[i], slaveXAABB))
            {
                intersection = true;
                break;
            }
    }   
    
    if(!intersection)
    {
        for (int i = 0; i < 12; i++)
        {
            const int index1 = DRT::UTILS::eleNodeNumbering_hex27_lines[i][0];
            const int index2 = DRT::UTILS::eleNodeNumbering_hex27_lines[i][1];
            if(isLineWithinXAABB(nodes[index1], nodes[index2], slaveXAABB))
            {
                intersection = true;
                break;
            }
        }
    }
    return intersection;
}


/*----------------------------------------------------------------------*
 | compute XAABB for contact search                           popp 08/08|
 | of a given query element (CONTACT)                                   |                             |
 *----------------------------------------------------------------------*/
BlitzMat3x2 GEO::computeContactXAABB( 
    DRT::Element*       element,
    const BlitzMat&     xyze)
{
  
    const int nsd = 2;
    BlitzMat3x2 XAABB;
    
    CONTACT::CElement* selement = static_cast<CONTACT::CElement*>(element);
    double area = selement->Area();
    
    // first node
    for(int dim=0; dim<nsd; ++dim)
    {
        XAABB(dim, 0) = xyze(dim, 0) - 0.1 * area;
        XAABB(dim, 1) = xyze(dim, 0) + 0.1 * area;
    }
    XAABB(2, 0) = 0.0;
    XAABB(2, 1) = 0.0;
    
    // remaining nodes
    for(int i=1; i<element->NumNode(); ++i)
    {
        for(int dim=0; dim<nsd; dim++)
        {
            XAABB(dim, 0) = std::min( XAABB(dim, 0), xyze(dim,i) - 0.1 * area);
            XAABB(dim, 1) = std::max( XAABB(dim, 1), xyze(dim,i) + 0.1 * area);
        }
    }
    
    /*
    printf("\n");
    printf("axis-aligned bounding box:\n minX = %f\n minY = %f\n minZ = %f\n maxX = %f\n maxY = %f\n maxZ = %f\n", 
              XAABB(0,0), XAABB(1,0), XAABB(2,0), XAABB(0,1), XAABB(1,1), XAABB(2,1));
    printf("\n");
    */
    
    return XAABB;
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
    BlitzVec3&                                  minDistanceVec,
    GEO::NearestObject&                         nearestObject)
{
  bool pointFound = false;
  double min_distance = GEO::LARGENUMBER;
  double distance = GEO::LARGENUMBER;
  BlitzVec3 normal = 0.0;
  BlitzVec3 x_surface = 0.0;
  std::map< int, std::set<int> > nodeList;

  // run over all surface elements
  for(std::map<int, std::set<int> >::const_iterator labelIter = elementList.begin(); labelIter != elementList.end(); labelIter++)
    for(std::set<int>::const_iterator eleIter = (labelIter->second).begin(); eleIter != (labelIter->second).end(); eleIter++)
    {    
      // not const because otherwise no lines can be obtained
      DRT::Element* element = dis.gElement(*eleIter);
      pointFound = GEO::getDistanceToSurface(element, currentpositions, point, x_surface, distance);
      if(pointFound && distance < min_distance)
      {
        pointFound = false;
        min_distance = distance;
        nearestObject.setSurfaceObjectType(*eleIter, labelIter->first, x_surface);
      }

      // run over all line elements
      const std::vector<Teuchos::RCP< DRT::Element> > eleLines = element->Lines();
      for(int i = 0; i < element->NumLine(); i++)
      {
        pointFound = GEO::getDistanceToLine(eleLines[i].get(), currentpositions, point, x_surface, distance);     
        if(pointFound && distance < min_distance)
        {
          pointFound = false;
          min_distance = distance;
          nearestObject.setLineObjectType(i, *eleIter, labelIter->first, x_surface);
        }
      }
      // collect nodes
      // 4 is correct because only the 4 corner nodes have to checked
      for(int i = 0; i < 4; i++)
        nodeList[labelIter->first].insert(element->NodeIds()[i]);
    }

  // run over all nodes
  for (std::map<int, std::set<int> >::const_iterator labelIter = nodeList.begin(); labelIter != nodeList.end(); labelIter++)
    for (std::set<int>::const_iterator nodeIter = (labelIter->second).begin(); nodeIter != (labelIter->second).end(); nodeIter++)
    {
      const DRT::Node* node = dis.gNode(*nodeIter);
      GEO::getDistanceToPoint(node, currentpositions, point, distance);
      if (distance < min_distance)
      {
        min_distance = distance;
        nearestObject.setNodeObjectType(*nodeIter, labelIter->first, currentpositions.find(node->Id())->second);
      }
    }

  // compute distance vector pointing away from the surface element
  // scalarproduct with element normal > 0 means fluid
  minDistanceVec(0) = point(0) - nearestObject.getPhysCoord()(0);
  minDistanceVec(1) = point(1) - nearestObject.getPhysCoord()(1);
  minDistanceVec(2) = point(2) - nearestObject.getPhysCoord()(2);

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
    BlitzVec3&                                  x_surface_phys,
    double&                                     distance)
{
  bool pointFound = false;
  double min_distance = GEO::LARGENUMBER;
  BlitzVec3 distance_vector = 0.0;
  BlitzVec2 elecoord = 0.0; // starting value at element center

  const BlitzMat xyze_surfaceElement(DRT::UTILS::getCurrentNodalPositions(surfaceElement, currentpositions));
  GEO::CurrentToSurfaceElementCoordinates(surfaceElement, xyze_surfaceElement, point, elecoord);

  if(GEO::checkPositionWithinElementParameterSpace(elecoord, surfaceElement->Shape()))
  { 
    GEO::elementToCurrentCoordinates(surfaceElement, xyze_surfaceElement, elecoord, x_surface_phys);
    // normal pointing away from the surface towards point
    distance_vector(0) = point(0) - x_surface_phys(0);
    distance_vector(1) = point(1) - x_surface_phys(1);
    distance_vector(2) = point(2) - x_surface_phys(2);
    distance = GEO::Norm2(distance_vector);
    min_distance = distance;
    pointFound = true;
  }

  // if the element is curved
  GEO::EleGeoType eleGeoType = GEO::HIGHERORDER;
  checkRoughGeoType(surfaceElement, xyze_surfaceElement, eleGeoType);

  // TODO fix check if deformed in the linear case
  if(eleGeoType == GEO::HIGHERORDER)
  {
    for(int i = 0; i < surfaceElement->NumNode(); i++)
    { 
      // use nodes as starting values
      for(int j = 0; j < 2; j++)
        elecoord(j) = DRT::UTILS::getEleNodeNumbering_nodes_reference(surfaceElement->Shape())[i][j];

      GEO::CurrentToSurfaceElementCoordinates(surfaceElement, xyze_surfaceElement, point, elecoord);

      if( GEO::checkPositionWithinElementParameterSpace(elecoord, surfaceElement->Shape()) )
      { 
        BlitzVec3 physcoord = 0.0;
        GEO::elementToCurrentCoordinates(surfaceElement, xyze_surfaceElement, elecoord, physcoord);
        // normal pointing away from the surface towards point
        distance_vector(0) = point(0) - physcoord(0);
        distance_vector(1) = point(1) - physcoord(1);
        distance_vector(2) = point(2) - physcoord(2);
        distance = GEO::Norm2(distance_vector);
        if(distance < min_distance)
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
bool GEO::getDistanceToLine(
    const DRT::Element*                         lineElement,
    const std::map<int,BlitzVec3>&              currentpositions,
    const BlitzVec3&                            point,
    BlitzVec3&                                  x_line_phys,
    double&                                     distance)
{
  bool pointFound = false;
  double min_distance = GEO::LARGENUMBER;
  BlitzVec3 distance_vector = 0.0;
  BlitzVec1 elecoord = 0.0; // starting value at element center  

  const BlitzMat xyze_lineElement(DRT::UTILS::getCurrentNodalPositions(lineElement, currentpositions));
  
  GEO::CurrentToLineElementCoordinates(lineElement, xyze_lineElement, point, elecoord);
  
  if(GEO::checkPositionWithinElementParameterSpace(elecoord, lineElement->Shape()))
  { 
    GEO::elementToCurrentCoordinates(lineElement, xyze_lineElement, elecoord, x_line_phys);
    // normal pointing away from the line towards point
    distance_vector(0) = point(0) - x_line_phys(0);
    distance_vector(1) = point(1) - x_line_phys(1);
    distance_vector(2) = point(2) - x_line_phys(2);
    distance = GEO::Norm2(distance_vector);
    min_distance = distance;
    pointFound = true;
  }

  // if the element is curved
  GEO::EleGeoType eleGeoType = GEO::HIGHERORDER;
  checkRoughGeoType(lineElement, xyze_lineElement, eleGeoType);

  if(eleGeoType == GEO::HIGHERORDER)
  {
    // use end nodes as starting values in addition
    for(int i = 0; i < 2; i++)
    { 
      // use end nodes as starting values in addition
      elecoord(0) = DRT::UTILS::getEleNodeNumbering_nodes_reference(lineElement->Shape())[i][0];

      GEO::CurrentToLineElementCoordinates(lineElement, xyze_lineElement, point, elecoord);

      if(GEO::checkPositionWithinElementParameterSpace(elecoord, lineElement->Shape()) )
      { 
        BlitzVec3 physcoord = 0.0;
        GEO::elementToCurrentCoordinates(lineElement, xyze_lineElement, elecoord, physcoord);
        // normal pointing away from the line towards point
        distance_vector(0) = point(0) - physcoord(0);
        distance_vector(1) = point(1) - physcoord(1);
        distance_vector(2) = point(2) - physcoord(2);
        distance = GEO::Norm2(distance_vector);
        if(distance < min_distance)
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
void GEO::getDistanceToPoint(
    const DRT::Node*                            node,
    const std::map<int,BlitzVec3>&              currentpositions,
    const BlitzVec3&                            point,
    double&                                     distance)
{

  BlitzVec3 distance_vector = 0.0;
  // node position in physical coordinates

  const BlitzVec3 x_node = currentpositions.find(node->Id())->second;

  // vector pointing away from the node towards physCoord
  distance_vector(0) = point(0) - x_node(0);
  distance_vector(1) = point(1) - x_node(1);
  distance_vector(2) = point(2) - x_node(2);

  // absolute distance between point and node
  distance = GEO::Norm2(distance_vector);
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
        adjacentElements.push_back( eleSet_node1[i]->Id() );
    }

  if(adjacentElements.size() > 2)
    dserror("more than two surfaces adjacent to a line - xfem coupling conditions might be not correct");
  
  if(adjacentElements.size() == 0)
      dserror("no adjacent surface elements found");

  return adjacentElements;
}



/*----------------------------------------------------------------------*
 |  computes the normal in a given point xsi in a            u.may 07/08|
 |  surface element (or elements if point is a on node                  |
 |  or on the line )                                                    |
 *----------------------------------------------------------------------*/
BlitzVec3 GEO::getNormalAtSurfacePoint(
    const DRT::Discretization&                  dis,
    const std::map<int,BlitzVec3>&              currentpositions,
    GEO::NearestObject&                         nearestObject)
{
  BlitzVec3 normal = 0.0;

  switch (nearestObject.getObjectType()) 
  {
  case GEO::SURFACE_OBJECT:
  { 
    BlitzVec2 elecoord = 0.0;
    const DRT::Element* surfaceElement = dis.gElement(nearestObject.getSurfaceId());
    const BlitzMat xyze_surfaceElement = DRT::UTILS::getCurrentNodalPositions(surfaceElement, currentpositions);
    GEO::CurrentToSurfaceElementCoordinates(surfaceElement, xyze_surfaceElement, nearestObject.getPhysCoord(), elecoord);
    GEO::computeNormalToSurfaceElement(surfaceElement, xyze_surfaceElement, elecoord, normal);
    break;
  }
  case GEO::LINE_OBJECT:
  {
    DRT::Element* surfaceElement = dis.gElement(nearestObject.getSurfaceId());
    const std::vector<Teuchos::RCP<DRT::Element> > eleLines = surfaceElement->Lines();
    const DRT::Element* lineElement = eleLines[nearestObject.getLineId()].get();
    std::vector<int> adjacentElements = getAdjacentSurfaceElementsToLine(lineElement->Nodes()[0], lineElement->Nodes()[1]);

    // run over all elements adjacent ot a line
    for(std::vector<int>::const_iterator eleIter = adjacentElements.begin(); eleIter != adjacentElements.end(); ++eleIter)
    {
      const DRT::Element* surfaceElement = dis.gElement(*eleIter);
      const BlitzMat xyze_surfaceElement(DRT::UTILS::getCurrentNodalPositions(surfaceElement, currentpositions));
      BlitzVec2 eleCoord = 0.0;
      BlitzVec3 surface_normal = 0.0;
      
      GEO::CurrentToSurfaceElementCoordinates(surfaceElement, xyze_surfaceElement, nearestObject.getPhysCoord(), eleCoord);
      GEO::computeNormalToSurfaceElement(surfaceElement, xyze_surfaceElement, eleCoord, surface_normal);
      normal += surface_normal;
    }
    normal /= ((double) adjacentElements.size()); 
    break;
  }
  case GEO::NODE_OBJECT:
  {
    const DRT::Node* node       = dis.gNode(nearestObject.getNodeId());
    // run over all elements adjacent ot a node
    for(int j=0; j<node->NumElement();j++)
    {      
      const DRT::Element* surfaceElement = node->Elements()[j];
      const BlitzMat xyze_surfaceElement(DRT::UTILS::getCurrentNodalPositions(surfaceElement, currentpositions));
      BlitzVec2 elecoord = 0.0;
      BlitzVec3 surface_normal = 0.0;
      GEO::CurrentToSurfaceElementCoordinates(surfaceElement, xyze_surfaceElement, nearestObject.getPhysCoord(), elecoord);
      GEO::computeNormalToSurfaceElement(surfaceElement, xyze_surfaceElement, elecoord, surface_normal);
      normal += surface_normal;
    }
    normal /= ((double) node->NumElement()); 
    break;
  }
  default:
    dserror("distance type does not exist");
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
    if( AABB_old(i,0) < (nodeBox(i,0)-GEO::TOL7) || AABB_old(i,1) > (nodeBox(i,1)+GEO::TOL7) ||
        AABB_new(i,0) < (nodeBox(i,0)-GEO::TOL7) || AABB_new(i,1) > (nodeBox(i,1)+GEO::TOL7)  )
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
void GEO::checkRoughGeoType(
    const DRT::Element*           element,
    const BlitzMat                xyze_element,
    GEO::EleGeoType&              eleGeoType)
{
  const DRT::Element::DiscretizationType distype = element->Shape();

  if(DRT::UTILS::getOrder(distype) ==1)
    eleGeoType = GEO::LINEAR;  //TODO check for bilinear elements in the tree they count as higerorder fix it
  else if(DRT::UTILS::getOrder(distype)==2)
    eleGeoType = GEO::HIGHERORDER;
  else
    dserror("order of element is not correct");

  // TODO check if a higher order element is linear
  // by checking if the higher order node is on the line
}

#endif  // #ifdef CCADISCRET
