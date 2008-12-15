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
#include "../drt_contact/drt_celement.H"


/*----------------------------------------------------------------------*
 | delivers a axis-aligned bounding box for a given          peder 07/08|
 | discretization                                                       |
 *----------------------------------------------------------------------*/
LINALG::Matrix<3,2> GEO::getXAABBofDis(
    const DRT::Discretization&                dis,
    const std::map<int,LINALG::Matrix<3,1> >& currentpositions)
{
  LINALG::Matrix<3,2> XAABB(true);
  if (dis.NumGlobalElements() == 0)
    return XAABB;
    
  if (dis.NumMyColElements() == 0){
    dserror("boundary should be ghosted and xfem discretization balanced, such that there are at least some xfem elements!");
    }

  // initialize XAABB as rectangle around the first point of dis
  const int nodeid = dis.lColElement(0)->Nodes()[0]->Id();
  const LINALG::Matrix<3,1> pos = currentpositions.find(nodeid)->second;
  for(int dim=0; dim<3; ++dim)
  {
    XAABB(dim, 0) = pos(dim) - GEO::TOL7;
    XAABB(dim, 1) = pos(dim) + GEO::TOL7;
  }

  // loop over elements and merge XAABB with their eXtendedAxisAlignedBoundingBox
  for (int j=0; j< dis.NumMyColElements(); ++j) 
  {
    const DRT::Element* element = dis.lColElement(j);
    const LINALG::SerialDenseMatrix xyze_element(GEO::getCurrentNodalPositions(element,currentpositions));
    GEO::EleGeoType eleGeoType(GEO::HIGHERORDER);
    GEO::checkRoughGeoType(element, xyze_element, eleGeoType);
    LINALG::Matrix<3,2> xaabbEle = GEO::computeFastXAABB(element, xyze_element, eleGeoType);
    XAABB = mergeAABB(XAABB, xaabbEle);
  }
  
  return XAABB;
}



/*----------------------------------------------------------------------*
 | delivers a axis-aligned bounding box for a given          peder 07/08|
 | discretization in reference configuration                            |
 *----------------------------------------------------------------------*/
LINALG::Matrix<3,2> GEO::getXAABBofDis(
    const DRT::Discretization&    dis)
{
  std::map<int,LINALG::Matrix<3,1> > currentpositions;

  for (int lid = 0; lid < dis.NumMyColNodes(); ++lid)
  {
    const DRT::Node* node = dis.lColNode(lid);
    LINALG::Matrix<3,1> currpos;
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
std::vector< LINALG::Matrix<3,2> > GEO::computeXAABBForLabeledStructures(
    const DRT::Discretization&                  dis,
    const std::map<int,LINALG::Matrix<3,1> >&   currentpositions,
    const std::map<int, std::set<int> >&        elementList)
{
  std::vector< LINALG::Matrix<3,2> > XAABBs;

  for (std::map<int, std::set<int> >::const_iterator labelIter = elementList.begin(); labelIter != elementList.end(); labelIter++)
  {
    LINALG::Matrix<3,2> xaabb_label;
    // initialize xaabb_label with box around first point
    const int eleId = *((labelIter->second).begin());
    const int nodeId = dis.gElement(eleId)->Nodes()[0]->Id();
    const LINALG::Matrix<3,1> pos = currentpositions.find(nodeId)->second;
    for(int dim=0; dim<3; ++dim)
    {
      xaabb_label(dim, 0) = pos(dim) - GEO::TOL7;
      xaabb_label(dim, 1) = pos(dim) + GEO::TOL7;
    }
    // run over set elements
    for (std::set<int>::const_iterator eleIter = (labelIter->second).begin(); eleIter != (labelIter->second).end(); eleIter++)
    {   
      const DRT::Element* element = dis.gElement(*eleIter);
      const LINALG::SerialDenseMatrix xyze_element(GEO::getCurrentNodalPositions(element,currentpositions));
      GEO::EleGeoType eleGeoType(GEO::HIGHERORDER);
      GEO::checkRoughGeoType(element, xyze_element, eleGeoType);
      LINALG::Matrix<3,2> xaabbEle = GEO::computeFastXAABB(element, xyze_element, eleGeoType);
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
    const DRT::Discretization&                    dis, 		
    const std::map<int,LINALG::Matrix<3,1> >&     currentpositions, 	
    const LINALG::Matrix<3,1>&                    querypoint,
    const std::map<int, std::set<int> >&          elementList)  		
{
  LINALG::Matrix<3,1> minDistanceVec(true);
 
  GEO::NearestObject nearestObject;

  // compute the distance to the nearest object (surface, line, node) return label of nearest object
  // returns the label of the surface element structure the projection of the query point is lying on
  int label = nearestObjectInNode(dis, currentpositions, elementList, querypoint, minDistanceVec, nearestObject);  	
 
  // compute normal in the point found on or in the object 
  LINALG::Matrix<3,1> normal = getNormalAtSurfacePoint(dis, currentpositions, nearestObject);  

  // compare normals and set label 
  const double scalarproduct = minDistanceVec(0)*normal(0) + minDistanceVec(1)*normal(1) + minDistanceVec(2)*normal(2);
  
  // if fluid
  if(scalarproduct > (-1)*GEO::TOL13)
    label = 0;

  return label;
}



/*----------------------------------------------------------------------*
 | returns a label for a given point                         u.may 07/08|
 | and element list                                                     |
 *----------------------------------------------------------------------*/
int GEO::getXFEMLabelAndNearestObject(
    const DRT::Discretization&                    dis,              
    const std::map<int,LINALG::Matrix<3,1> >&     currentpositions, 
    const LINALG::Matrix<3,1>&                    querypoint,      
    const std::map<int, std::set<int> >&          elementList,     
    GEO::NearestObject&                           nearestObject)
{
  LINALG::Matrix<3,1> minDistanceVec(true);

  // compute the distance to the nearest object (surface, line, node) return label of nearest object
  // returns the label of the surface element structure the projection of the query point is lying on
  int label = nearestObjectInNode(dis, currentpositions, elementList, querypoint, minDistanceVec, nearestObject);   
 
  // compute normal in the point found on or in the object 
  LINALG::Matrix<3,1> normal= getNormalAtSurfacePoint(dis, currentpositions, nearestObject);  

  // compare normals and set label 
  const double scalarproduct = minDistanceVec(0)*normal(0) + minDistanceVec(1)*normal(1) + minDistanceVec(2)*normal(2);
  
  // if fluid
  if(scalarproduct > (-1)*GEO::TOL13)
    label = 0;

  return label;
}



/*----------------------------------------------------------------------*
 | a set of nodes in a given radius                          u.may 07/08|
 | from a query point                                                   |
 *----------------------------------------------------------------------*/
std::map<int,std::set<int> > GEO::getElementsInRadius(
    const DRT::Discretization&                  dis, 		
    const std::map<int,LINALG::Matrix<3,1> >&   currentpositions, 	
    const LINALG::Matrix<3,1>&                  querypoint,
    const double		                            radius,
    const int                                   label,
    std::map<int, std::set<int> >&              elementList)  
{
  std::map< int, std::set<int> >  nodeList;
  std::map<int,std::set<int> >    elementMap;
  
  // collect all nodes with different label
  for(std::map<int, std::set<int> >::const_iterator labelIter = elementList.begin(); labelIter != elementList.end(); labelIter++)
  {
    if(label != labelIter->first) // don t collect nodes which belong to the same label
    {
      for(std::set<int>::const_iterator eleIter = (labelIter->second).begin(); eleIter != (labelIter->second).end(); eleIter++)
      {
        DRT::Element* element = dis.gElement(*eleIter);
        for(int i = 0; i < DRT::UTILS::getNumberOfElementCornerNodes(element->Shape()); i++)
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
      {
        for(int i=0; i<dis.gNode(*nodeIter)->NumElement();i++)  
          elementMap[labelIter->first].insert(dis.gNode(*nodeIter)->Elements()[i]->Id()); 
      }
    }
  
  return elementMap;
}


/*----------------------------------------------------------------------*
 | a vector of intersection elements                          popp 07/08|
 | for a given query element (CONTACT)                                  |                                 |
 *----------------------------------------------------------------------*/
std::vector<int> GEO::getIntersectionElements(
    const DRT::Discretization&                  dis,    
    const std::map<int,LINALG::Matrix<3,1> >&   currentpositions,   
    DRT::Element*                               element,
    std::map<int, std::set<int> >&              elementList)  
{
  // vector containing sought-after element ids
  vector<int> intersectionids;
  
  // create XAABB for query element
  const LINALG::SerialDenseMatrix xyze(GEO::getCurrentNodalPositions(element,currentpositions));  
  LINALG::Matrix<3,2> slaveXAABB = GEO::computeContactXAABB(element,xyze);
  
  // loop over all entries of elementList (= intersection candidates)
  for(std::set<int>::const_iterator elementIter = (elementList.begin()->second).begin();
      elementIter != (elementList.begin()->second).end(); elementIter++)
  {
    DRT::Element* masterelement = dis.gElement(*elementIter);
    const LINALG::SerialDenseMatrix masterxyze(GEO::getCurrentNodalPositions(masterelement,currentpositions));
    LINALG::Matrix<3,2> masterXAABB = GEO::computeContactXAABB(masterelement,masterxyze);
    
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
    const LINALG::Matrix<3,2>&     slaveXAABB, 
    const LINALG::Matrix<3,2>&     masterXAABB)
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
    static std::vector< LINALG::Matrix<3,1> > nodes(8, LINALG::Matrix<3,1>());
    
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
LINALG::Matrix<3,2> GEO::computeContactXAABB( 
    DRT::Element*                         element,
    const LINALG::SerialDenseMatrix&      xyze)
{
  
    const int nsd = 2;
    LINALG::Matrix<3,2> XAABB;
    
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
    const std::map<int,LINALG::Matrix<3,1> >&   currentpositions,
    const std::map<int, std::set<int> >&        elementList,
    const LINALG::Matrix<3,1>&                  point,
    LINALG::Matrix<3,1>&                        minDistanceVec,
    GEO::NearestObject&                         nearestObject)
{
  bool pointFound = false;
  double min_distance = GEO::LARGENUMBER;
  double distance = GEO::LARGENUMBER;
  LINALG::Matrix<3,1> normal(true);
  LINALG::Matrix<3,1> x_surface(true);
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
      // TODO make properly
      for(int i = 0; i < DRT::UTILS::getNumberOfElementCornerNodes(element->Shape()); i++)
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
  minDistanceVec.Update(1.0, point, -1.0, nearestObject.getPhysCoord());


  return nearestObject.getLabel();
}



/*----------------------------------------------------------------------*
 |  computes the normal distance from a point to a           u.may 07/08|
 |  surface element, if it exits                                        |
 *----------------------------------------------------------------------*/
bool GEO::getDistanceToSurface(
    const DRT::Element*                         surfaceElement,
    const std::map<int,LINALG::Matrix<3,1> >&   currentpositions,
    const LINALG::Matrix<3,1>&                  point,
    LINALG::Matrix<3,1>&                        x_surface_phys,
    double&                                     distance)
{
  bool pointFound = false;
  double min_distance = GEO::LARGENUMBER;
  LINALG::Matrix<3,1> distance_vector(true);
  LINALG::Matrix<2,1> elecoord(true); // starting value at element center

  const LINALG::SerialDenseMatrix xyze_surfaceElement(GEO::getCurrentNodalPositions(surfaceElement, currentpositions));
  GEO::CurrentToSurfaceElementCoordinates(surfaceElement->Shape(), xyze_surfaceElement, point, elecoord);

  if(GEO::checkPositionWithinElementParameterSpace(elecoord, surfaceElement->Shape()))
  { 
    GEO::elementToCurrentCoordinates(surfaceElement->Shape(), xyze_surfaceElement, elecoord, x_surface_phys);
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
  if(eleGeoType == GEO::HIGHERORDER)
  {
    for(int i = 0; i < surfaceElement->NumNode(); i++)
    { 
      // use nodes as starting values
      for(int j = 0; j < 2; j++)
        elecoord(j) = DRT::UTILS::getEleNodeNumbering_nodes_paramspace(surfaceElement->Shape())(j, i);

      GEO::CurrentToSurfaceElementCoordinates(surfaceElement->Shape(), xyze_surfaceElement, point, elecoord);

      if( GEO::checkPositionWithinElementParameterSpace(elecoord, surfaceElement->Shape()) )
      { 
        LINALG::Matrix<3,1> physcoord(true);
        GEO::elementToCurrentCoordinates(surfaceElement->Shape(), xyze_surfaceElement, elecoord, physcoord);
        // normal pointing away from the surface towards point
        distance_vector.Update(1.0, point, -1.0, physcoord);
        distance = distance_vector.Norm2();
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
    const std::map<int,LINALG::Matrix<3,1> >&   currentpositions,
    const LINALG::Matrix<3,1>&                  point,
    LINALG::Matrix<3,1>&                        x_line_phys,
    double&                                     distance)
{
  bool pointFound = false;
  double min_distance = GEO::LARGENUMBER;
  LINALG::Matrix<3,1> distance_vector(true);
  LINALG::Matrix<1,1> elecoord(true); // starting value at element center  

  const LINALG::SerialDenseMatrix xyze_lineElement(GEO::getCurrentNodalPositions(lineElement, currentpositions));
  
  GEO::CurrentToLineElementCoordinates(lineElement->Shape(), xyze_lineElement, point, elecoord);
  
  if(GEO::checkPositionWithinElementParameterSpace(elecoord, lineElement->Shape()))
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

  if(eleGeoType == GEO::HIGHERORDER)
  {
    // use end nodes as starting values in addition
    for(int i = 0; i < 2; i++)
    { 
      // use end nodes as starting values in addition
      elecoord(0) = DRT::UTILS::getEleNodeNumbering_nodes_paramspace(lineElement->Shape())(0,i);

      GEO::CurrentToLineElementCoordinates(lineElement->Shape(), xyze_lineElement, point, elecoord);

      if(GEO::checkPositionWithinElementParameterSpace(elecoord, lineElement->Shape()) )
      { 
        LINALG::Matrix<3,1> physcoord(true);
        GEO::elementToCurrentCoordinates(lineElement->Shape(), xyze_lineElement, elecoord, physcoord);
        // normal pointing away from the line towards point
        distance_vector.Update( 1.0, point, -1.0, physcoord);
        distance = distance_vector.Norm2();
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
    const std::map<int,LINALG::Matrix<3,1> >&   currentpositions,
    const LINALG::Matrix<3,1>&                  point,
    double&                                     distance)
{

  
  // node position in physical coordinates
  const LINALG::Matrix<3,1> x_node = currentpositions.find(node->Id())->second;

  LINALG::Matrix<3,1> distance_vector;
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
LINALG::Matrix<3,1> GEO::getNormalAtSurfacePoint(
    const DRT::Discretization&                  dis,
    const std::map<int,LINALG::Matrix<3,1> >&   currentpositions,
    GEO::NearestObject&                         nearestObject)
{
  LINALG::Matrix<3,1> normal(true);

  switch (nearestObject.getObjectType()) 
  {
  case GEO::SURFACE_OBJECT:
  { 
    LINALG::Matrix<2,1> elecoord(true);
    const DRT::Element* surfaceElement = dis.gElement(nearestObject.getSurfaceId());
    const LINALG::SerialDenseMatrix xyze_surfaceElement = GEO::getCurrentNodalPositions(surfaceElement, currentpositions);
    GEO::CurrentToSurfaceElementCoordinates(surfaceElement->Shape(), xyze_surfaceElement, nearestObject.getPhysCoord(), elecoord);
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
      const LINALG::SerialDenseMatrix xyze_surfaceElement(GEO::getCurrentNodalPositions(surfaceElement, currentpositions));
      LINALG::Matrix<2,1> eleCoord(true);
      LINALG::Matrix<3,1> surface_normal(true);
      
      GEO::CurrentToSurfaceElementCoordinates(surfaceElement->Shape(), xyze_surfaceElement, nearestObject.getPhysCoord(), eleCoord);
      GEO::computeNormalToSurfaceElement(surfaceElement, xyze_surfaceElement, eleCoord, surface_normal);
      normal += surface_normal;
    }
    normal.Scale(1.0/((double) adjacentElements.size())); 
    break;
  }
  case GEO::NODE_OBJECT:
  {
    const DRT::Node* node       = dis.gNode(nearestObject.getNodeId());
    // run over all elements adjacent ot a node
    for(int j=0; j<node->NumElement();j++)
    {      
      const DRT::Element* surfaceElement = node->Elements()[j];
      const LINALG::SerialDenseMatrix xyze_surfaceElement(GEO::getCurrentNodalPositions(surfaceElement, currentpositions));
      LINALG::Matrix<2,1> elecoord(true);
      LINALG::Matrix<3,1> surface_normal(true);
      GEO::CurrentToSurfaceElementCoordinates(surfaceElement->Shape(), xyze_surfaceElement, nearestObject.getPhysCoord(), elecoord);
      GEO::computeNormalToSurfaceElement(surfaceElement, xyze_surfaceElement, elecoord, surface_normal);
      normal += surface_normal;
    }
    normal.Scale(1.0/((double) node->NumElement())); 
    break;
  }
  default:
    dserror("distance type does not exist");
  }
  return normal;
}



/*----------------------------------------------------------------------*
 | check s if nearest point found on nearest object          u.may 09/08|
 | lies in a tree node specified by its box                             |
 *----------------------------------------------------------------------*/
bool  GEO::pointInTreeNode(
    const LINALG::Matrix<3,1>&    point,
    const LINALG::Matrix<3,2>&    nodeBox) 
{
  for(int dim=0; dim<3; dim++)
    if ( (point(dim) < nodeBox(dim,0)) || (point(dim) > nodeBox(dim,1)) ) // no tolerances here !!!
      return false;
  
  return true;
}



/*----------------------------------------------------------------------*
 | check s if nearest point found on nearest object          u.may 09/08|
 | lies in the minimum circle that fits inside the triangle             |
 *----------------------------------------------------------------------*/
bool  GEO::pointInMinCircleInTreeNode(
    const LINALG::Matrix<3,1>&          nearestpoint,
    const LINALG::Matrix<3,1>&          querypoint,
    const LINALG::Matrix<3,2>&          nodeBox,
    const bool                          rootNode) 
{
  double minRadius = GEO::LARGENUMBER;
  
  if(rootNode)
    return pointInTreeNode(nearestpoint,nodeBox); 
  else
  {
    // determine minimum distance querypoint to node box wall
    for(int dim=0; dim<3; dim++)
      for(int i=0; i<2; i++)
      {
        double actRadius = fabs(querypoint(dim) - nodeBox(dim,i));
        if ( actRadius < minRadius ) 
          minRadius = actRadius;
      }
  }
  
  // distance querypoint - nearest point
  LINALG::Matrix<3,1> distance_vector;
  distance_vector.Update(1.0, querypoint, -1.0, nearestpoint);
 
  // absolute distance between point and node
  if(distance_vector.Norm2() < minRadius)
    return true;
  
  return false;
}



/*----------------------------------------------------------------------*
 | merge two AABB and deliver the resulting AABB's           peder 07/08|
 *----------------------------------------------------------------------*/
LINALG::Matrix<3,2> GEO::mergeAABB(
    const LINALG::Matrix<3,2>& AABB1, 
    const LINALG::Matrix<3,2>& AABB2)
{
  LINALG::Matrix<3,2> mergedAABB;

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
    const LINALG::Matrix<3,2>&  AABB_old, 
    const LINALG::Matrix<3,2>&  AABB_new,
    const LINALG::Matrix<3,2>&  nodeBox)
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
    const DRT::Element*               element,
    const LINALG::SerialDenseMatrix   xyze_element,
    GEO::EleGeoType&                  eleGeoType)
{
  const int order = DRT::UTILS::getOrder(element->Shape());

  if(order ==1)
    eleGeoType = GEO::LINEAR;  //TODO check for bilinear elements in the tree they count as higerorder fix it
  else if(order==2)
    eleGeoType = GEO::HIGHERORDER;
  else
    dserror("order of element is not correct");

  // TODO check if a higher order element is linear
  // by checking if the higher order node is on the line
}



/*----------------------------------------------------------------------*
 | returns a set of intersection candidate ids             u.may   09/08|
 *----------------------------------------------------------------------*/
std::vector<int> GEO::getIntersectionCandidates(    
    const DRT::Discretization&                  dis,     
    const std::map<int,LINALG::Matrix<3,1> >&   currentpositions,   
    DRT::Element*                               xfemElement,
    std::map<int, std::set<int> >&              elementList)
{ 
  std::vector<int> intersectionCandidateIds;
  
  // create XAABB for query xfem  element
  GEO::EleGeoType xfemGeoType = HIGHERORDER;
  const LINALG::SerialDenseMatrix xyze_xfemElement(GEO::InitialPositionArray(xfemElement));
  GEO::checkGeoType(xfemElement, xyze_xfemElement, xfemGeoType);
  const LINALG::Matrix<3,2> xfemXAABB(GEO::computeFastXAABB(xfemElement, xyze_xfemElement, xfemGeoType));
  
  // loop over all entries of elementList (= intersection candidates)
  // run over global ids
  for(std::set<int>::const_iterator elementIter = (elementList.begin()->second).begin();
      elementIter != (elementList.begin()->second).end(); elementIter++)
  {
    DRT::Element*  cutterElement = dis.gElement(*elementIter);
    const LINALG::SerialDenseMatrix xyze_cutterElement(GEO::getCurrentNodalPositions(cutterElement, currentpositions));
    GEO::EleGeoType cutterGeoType = HIGHERORDER;
    GEO::checkGeoType(cutterElement, xyze_cutterElement, cutterGeoType);
    const LINALG::Matrix<3,2> cutterXAABB(GEO::computeFastXAABB(cutterElement, xyze_cutterElement, cutterGeoType));
    // compare 
    if(intersectionOfXAABB(cutterXAABB, xfemXAABB))
      intersectionCandidateIds.push_back(cutterElement->Id());
  }

  return intersectionCandidateIds;
}

#endif  // #ifdef CCADISCRET
