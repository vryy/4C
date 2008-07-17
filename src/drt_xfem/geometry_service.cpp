/*!
\file geometry_service.cpp

\brief provides a class with search tree

<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
 */
#ifdef CCADISCRET
#include "geometry_service.H"
#include "intersection_service.H"

using namespace std;
using namespace XFEM;

  
BlitzMat3x2 XFEM::getXAABBofDis(const DRT::Discretization& dis,const std::map<int,BlitzVec3>& currentpositions){
  const int nsd = 3;
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
  for(int dim=0; dim<nsd; ++dim)
  {
    XAABB(dim, 0) = pos[dim] - TOL7;
    XAABB(dim, 1) = pos[dim] + TOL7;
  }
  
  // loop over elements and merge XAABB with their eXtendedAxisAlignedBoundingBox
  //BlitzMat curr = BlitzMat(3, currentpositions.size());
  
  
  for (int j=0; j< dis.NumMyColElements(); ++j) {
    const DRT::Element* ele = dis.lColElement(j);
    BlitzMat curr(3, ele->NumNode());
    for (int inode = 0; inode < ele->NumNode(); ++inode)
    {
      const int nodeid = ele->Nodes()[inode]->Id();
      const BlitzVec3 pos = currentpositions.find(nodeid)->second;
      if (currentpositions.find(nodeid) == currentpositions.end())
      {
        dserror("global node id not found. bug!");
      }
      for (int k=0; k<3; k++)
      {
        curr(k, inode) = pos(k);
      }
    }

    BlitzMat3x2 xaabbEle = XFEM::computeFastXAABB(ele, curr, HIGHERORDER);
    XAABB = mergeAABB(XAABB, xaabbEle);
  }
  
  #ifdef DEBUG5  
    cout << "_XAABB=" << XAABB(0,0) << ","<< XAABB(0,1) << "," << XAABB(1,0) << "," << XAABB(1,1) << "," << XAABB(2,0) << "," << XAABB(2,1) << endl;
  #endif
  return XAABB;
}

BlitzMat3x2 XFEM::getXAABBofDis(const DRT::Discretization& dis){
  std::map<int,BlitzVec3> currentpositions;
  currentpositions.clear();
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



const DRT::Element* XFEM::nearestXAABBNeighbourInList(
    const DRT::Discretization& dis,
    const std::map<int,BlitzVec3>& currentpositions, 
    const list<const DRT::Element* >& ElementList, 
    const BlitzVec3& x_in, 
    double& dist)
{
  dist=1.0e12;
  const DRT::Element* closest_element;
  for (list< const DRT::Element* >::const_iterator myIt = ElementList.begin(); myIt != ElementList.end(); myIt++)
  {
    double distance = 1.0e12;
    const BlitzMat xyze(getCurrentNodalPositions(*myIt, currentpositions));
    distance = (getMaxDistanceFromAABB(x_in, computeFastXAABB(*myIt, xyze, HIGHERORDER)));
    if (distance < dist){
      closest_element = *myIt;   
      dist = distance;      
    }
  }
  return closest_element;
}

const DRT::Element* XFEM::nearestNeighbourInList(const DRT::Discretization& dis,const std::map<int,BlitzVec3>& currentpositions, const list<const DRT::Element* >& ElementList, const BlitzVec3& x_in, double& dist)
{
  return nearestNeighbourInListNew(dis,currentpositions, ElementList, x_in, dist);
}

const DRT::Element* XFEM::nearestNeighbourInListOld(const DRT::Discretization& dis,const std::map<int,BlitzVec3>& currentpositions, const list<const DRT::Element* >& ElementList, const BlitzVec3& x_in, double& dist)
{
  dist=1.0e12;
  bool in_element = false;
  double min_ele_distance = 1.0e12;
  double min_node_distance =1.0e12;
  BlitzVec3 vectorX2minNode;
  const DRT::Element* closest_element;
  const DRT::Element* closest_element_from_node;
  const DRT::Node* closest_node;
  bool foundNearSurfaceElement =false;
  
  for (list< const DRT::Element* >::const_iterator myIt = ElementList.begin(); myIt != ElementList.end(); myIt++)
  {
    double distance = 1.0e12;
    const DRT::Element* cutterele = &**myIt;
    const BlitzMat xyze_cutter(getCurrentNodalPositions(cutterele, currentpositions));
    BlitzVec2 eleCoord;
    BlitzVec3 normal;
    in_element = XFEM::searchForNearestPointOnSurface(cutterele,xyze_cutter,x_in,eleCoord,normal,distance);
    if (in_element && (fabs(distance) < fabs(min_ele_distance)))
    {
      closest_element = cutterele;
      min_ele_distance = distance;
    }
  }
  { 
    min_node_distance = fabs(min_ele_distance);
    for (list< const DRT::Element* >::const_iterator myIt2 = ElementList.begin(); myIt2 != ElementList.end(); myIt2++) {
      const DRT::Element* cutterele = &**myIt2;       
      const int numnode = cutterele->NumNode();
      const DRT::Node*const* nodes = cutterele->Nodes();          
      for (int inode = 0; inode < numnode; ++inode)
      {
        BlitzVec3 vector;
        const DRT::Node* cutternode = nodes[inode]; 
        // node position in physical coordinates
        const BlitzVec3 x_node = currentpositions.find(cutternode->Id())->second;
        
        // vector pointing away from the node towards physCoord
        vector(0) = x_in(0) - x_node(0);
        vector(1) = x_in(1) - x_node(1);
        vector(2) = x_in(2) - x_node(2);
        // absolute distance between point and node
        const double distance = sqrt(vector(0)*vector(0) + vector(1)*vector(1) + vector(2)*vector(2));
        
        if (distance < min_node_distance) {
          closest_node = cutternode;
          min_node_distance = distance;
          vectorX2minNode=vector;
        }
      }          
    }
    if (fabs(min_node_distance)< fabs(min_ele_distance)) {
    BlitzVec3 normal;
    normal(0)=0;normal(1)=0;normal(2)=0;
    DRT::Node*  node = dis.gNode(closest_node->Id());
    const BlitzVec3 closest_node_pos = currentpositions.find(closest_node->Id())->second;
    for(int j=0; j<node->NumElement();j++)
    {
      DRT::Element* surfaceElement = node->Elements()[j];
      BlitzMat xyze_surfaceElement(getCurrentNodalPositions(surfaceElement, currentpositions));
      BlitzVec3 eleNormalAtXsi;
      BlitzVec2 xsi;
      CurrentToSurfaceElementCoordinates(surfaceElement, xyze_surfaceElement, closest_node_pos, xsi);
      computeNormalToBoundaryElement(surfaceElement, xyze_surfaceElement, xsi, eleNormalAtXsi);
      normal(0) += eleNormalAtXsi(0);
      normal(1) += eleNormalAtXsi(1);
      normal(2) += eleNormalAtXsi(2);
    }
    closest_element_from_node=node->Elements()[0]; 
    const double scalarproduct = vectorX2minNode(0)*normal(0) + vectorX2minNode(1)*normal(1) + vectorX2minNode(2)*normal(2);
    const double vorzeichen = scalarproduct/fabs(scalarproduct);
    min_node_distance *= vorzeichen;
    }
  }
  
  if (fabs(min_ele_distance) <= fabs(min_node_distance)){
    dist = min_ele_distance;
    return closest_element;
  }
  else {
    dist = min_node_distance;
    return closest_element_from_node;
  }
  
  return closest_element;
}


const DRT::Element* XFEM::nearestNeighbourInListNew(
    const DRT::Discretization& dis,
    const std::map<int,BlitzVec3>& currentpositions, 
    const list<const DRT::Element* >& ElementList, 
    const BlitzVec3& X, 
    double& distance)
{
  double min_ele_distance = 1.0e12;
  const DRT::Element* closest_element;
  int closestElement_Type=0;
  BlitzVec3 vectorX2minDistPoint;
  BlitzVec2 xsi;

  bool printDEBUG = false;

  BlitzVec3 testpoint(0.407,0.407,0.0933);
  double eps = 0.01;
  if (X(0)>testpoint(0)-eps && X(0)<testpoint(0)+eps && X(1)>testpoint(1)-eps && X(1)<testpoint(1)+eps && X(2)>testpoint(2)-eps && X(2)<testpoint(2)+eps)
    printDEBUG = true;
  if (printDEBUG)  
  {
    cout << "nbr of candidates in list : "<< ElementList.size()<< endl;
  }
  
  int ElementType=0;
  for (list< const DRT::Element* >::const_iterator myIt = ElementList.begin(); myIt != ElementList.end(); myIt++)
  {
    double distance = 1.0e12;
   
    const DRT::Element* cutterele = &**myIt;
    if (printDEBUG)
      cout << X(0) << "/" << X(1)<< "/" << X(2) <<" -> candID: " << cutterele->Id() << endl;
    const BlitzMat xyze_cutter(getCurrentNodalPositions(cutterele, currentpositions));
    BlitzVec2 eleCoord(0.0,0.0);
    BlitzVec3 vector2minDistPoint(0.0,0.0,0.0);
    distance = XFEM::getSquaredElementDistance(cutterele,currentpositions,X,eleCoord,vector2minDistPoint,ElementType);
    if (distance < min_ele_distance)
    {
      closest_element = cutterele;
      closestElement_Type = ElementType;
      xsi = eleCoord;
      min_ele_distance = distance;
      vectorX2minDistPoint = vector2minDistPoint;
    }
  }
  
  BlitzVec3 normal= getNormalAtXsi(dis, closest_element, currentpositions, xsi,X, closestElement_Type);  
  const double scalarproduct = vectorX2minDistPoint(0)*normal(0) + vectorX2minDistPoint(1)*normal(1) + vectorX2minDistPoint(2)*normal(2);
  min_ele_distance = sqrt(min_ele_distance);
  if (scalarproduct<0.0)
    min_ele_distance *= (-1);
  if (scalarproduct==0.0){
    cout << "scalarprod is 0, just for testing ";
  }
    
  distance =  min_ele_distance;
 return closest_element;
}
 
BlitzVec3 XFEM::getNormalAtXsi(
    const DRT::Discretization& dis,
    const DRT::Element* surfaceElement,
    const std::map<int,BlitzVec3>& currentpositions, 
    const BlitzVec2& xsi, 
    const BlitzVec3& X,
    const int& ElementType)
{ 
  BlitzVec3 normal(0.0,0.0,0.0);

  switch (ElementType) {
  case DISTANCE_TO_ELEMENT_SURFACE:{     
    const BlitzMat xyze_surfaceElement = getCurrentNodalPositions(surfaceElement, currentpositions);
    computeNormalToBoundaryElement(surfaceElement, xyze_surfaceElement, xsi, normal);
    return normal;
    break;
  }
  case DISTANCE_TO_ELEMENT_LINE:
  {
    if (!( (xsi(0)==(-1)) || (xsi(1)==(-1)) || (xsi(0)==1) || (xsi(1)==1) )){
      cout << "xsi: " << xsi(0) <<", "<<xsi(1)<<endl;
      dserror("no line problem wrong path");
    }
    BlitzVec2 xsiA(0.0,0.0);
    BlitzVec2 xsiB(0.0,0.0);
    if (xsi(0)==(-1) || xsi(0)==1){
      xsiA(0) =xsi(0);
      xsiA(1) =(-1.0);
      xsiB(0) =xsi(0);
      xsiB(1) =1.0;
    }        
    else {
      xsiA(0)=(-1.0);
      xsiA(1)=xsi(1);
      xsiB(0)=1.0;
      xsiB(1)=xsi(1);
    }
    const BlitzMat xyze_surfaceElement = getCurrentNodalPositions(surfaceElement, currentpositions);
    const DRT::Node* NodeA = getNodeAtXsi(surfaceElement, currentpositions, xsiA);
    const DRT::Node* NodeB = getNodeAtXsi(surfaceElement, currentpositions, xsiB);
    list< const DRT::Element* > commonElements = getCommonElements(NodeA, NodeB);
    list< BlitzVec3 > normalVectors;
    normalVectors.clear();
    if (commonElements.size() > 2 or commonElements.size() == 0){
      // note that 1 line is ok if "2d" problems are calculated, where in z-direction we have no closed surface
      // it would be also ok, for adaptivity with hanging nodes and hence hanging lines
      dserror("a line should be bounded by 2 surface elements");
      }
    for(list< const DRT::Element* >::const_iterator myIt = commonElements.begin(); myIt != commonElements.end(); ++myIt ){
      BlitzMat xyze_surfaceElement(getCurrentNodalPositions(*myIt, currentpositions));
      BlitzVec3 tmpNormal;
      BlitzVec2 tmpXsi;
      CurrentToSurfaceElementCoordinates(*myIt, xyze_surfaceElement, X, tmpXsi);
      computeNormalToBoundaryElement(*myIt, xyze_surfaceElement, tmpXsi, tmpNormal);
      normalVectors.push_back(tmpNormal);
    }
    return addVectors(normalVectors);
    break;
  } 
  case DISTANCE_TO_ELEMENT_POINT:
  {
    bool printDEBUG = false;
    const BlitzMat xyze_surfaceElement = getCurrentNodalPositions(surfaceElement, currentpositions);
    const DRT::Node*  node = getNodeAtXsi(surfaceElement, currentpositions, xsi);
    const BlitzVec3 closest_node_pos = currentpositions.find(node->Id())->second;
    for(int j=0; j<node->NumElement();j++)
    {
      const DRT::Element* tmpSurfaceElement = node->Elements()[j];
      if (printDEBUG)
        cout << X(0) << "/" << X(1)<< "/" << X(2) <<" -> candID2: " << tmpSurfaceElement->Id() << endl;
      BlitzMat xyze_surfaceElement(getCurrentNodalPositions(tmpSurfaceElement, currentpositions));
      BlitzVec3 eleNormalAtXsi;
      BlitzVec2 tmpXsi;
      CurrentToSurfaceElementCoordinates(tmpSurfaceElement, xyze_surfaceElement, closest_node_pos, tmpXsi);
      computeNormalToBoundaryElement(tmpSurfaceElement, xyze_surfaceElement, tmpXsi, eleNormalAtXsi);
      normal(0) = normal(0) +  eleNormalAtXsi(0);
      normal(1) = normal(1) +  eleNormalAtXsi(1);
      normal(2) = normal(2) +  eleNormalAtXsi(2);
    }
    return normal;
    break;
  }
  }
  cout << "ElementType: " << ElementType << endl;
  return normal;
  dserror("should not get here");
  return normal;
}

list < const DRT::Element* > XFEM::getCommonElements(const DRT::Node* A, const DRT::Node* B){
  list < const DRT::Element* > commonEles;
  commonEles.clear();
  const DRT::Element*const* tmp1 = A->Elements();
  const DRT::Element*const* tmp2 = B->Elements();
  for(int i=0; i<A->NumElement();i++)
   {
    for(int j=0; j<B->NumElement();j++)
     {      
      if (tmp1[i]->Id()==tmp2[j]->Id())
        commonEles.push_back(A->Elements()[i]);
     }
   }
  return commonEles;
}
    
BlitzVec3 XFEM::addVectors(const list< BlitzVec3 > vectors){
    BlitzVec3 vecSum(0.0,0.0,0.0); 
    for(list< BlitzVec3 >::const_iterator myIt = vectors.begin(); myIt!=vectors.end();++myIt)
    {
      vecSum(0) = vecSum(0) +  (*myIt)(0);
      vecSum(1) = vecSum(1) +  (*myIt)(1);
      vecSum(2) = vecSum(2) +  (*myIt)(2);
    }
    return vecSum;
  }
  
/*----------------------------------------------------------------------*
 |  RQI:    searches the nearest point on a surface          peder 07/08|
 |          element for a given point in physical coordinates           |
 *----------------------------------------------------------------------*/
double XFEM::getSquaredElementDistance(
    const DRT::Element*                     surfaceElement,
    const std::map<int,BlitzVec3>&          currentpositions,
    const BlitzVec3&                        physCoord,
    BlitzVec2&                              xsi,
    BlitzVec3&                              vector2minDistPoint,
    int&                                    distanceType)
{  
  /* 1.) get nearest point on element support in xsi-coordinates
   *  2.) if ((xsi0<xsi0_min||xsi0>xsi0_max) XOR (xsi1<xsi1_min || xsi1>xsi1_max)) 
   *      nearest point is on element surface (which is a line because the element itself is 2D) 
   *  3.) if xsi is totaly out of bounds the nearest point on this surface is one of the node points
   *  4.) if its perfectly in xsi-limits nearest point is on the element surface
   */ 
 
  double distance = 1.0e12;
  BlitzVec3 normal;
  BlitzMat2x2 xsiBoundingBox = XFEM::getXsiBoundingBox(surfaceElement);
  BlitzMat xyze_surfaceElement(getCurrentNodalPositions(surfaceElement, currentpositions));
  CurrentToSurfaceElementCoordinates(surfaceElement, xyze_surfaceElement, physCoord, xsi);

  if ( (xsi(0) < xsiBoundingBox(0,0) || xsi(0) > xsiBoundingBox(0,1) ) && 
      (xsi(1) < xsiBoundingBox(1,0) || xsi(1) > xsiBoundingBox(1,1) )    )   {
    distanceType = DISTANCE_TO_ELEMENT_POINT;
//    cout<<" point type xsi:" <<  xsi(0)<<", "<<xsi(1)<< endl;
  }
  else {
    if ( xsi(0)>=xsiBoundingBox(0,0) && xsi(0)<=xsiBoundingBox(0,1) && xsi(1)>=xsiBoundingBox(1,0) && xsi(1)<=xsiBoundingBox(1,1) ){
      distanceType = DISTANCE_TO_ELEMENT_SURFACE;      
//      cout<<" surface type xsi:" <<  xsi(0)<<", "<<xsi(1)<< endl;
    }
    else {
      distanceType = DISTANCE_TO_ELEMENT_LINE;   
//      cout<<" line type xsi:" <<  xsi(0)<<", "<<xsi(1)<< endl;
    }
  }

  switch (distanceType) {
  case DISTANCE_TO_ELEMENT_SURFACE:    
    distance = getSquaredElementDistance_Surface(surfaceElement, currentpositions,physCoord, xsi, vector2minDistPoint);
    break;
  case DISTANCE_TO_ELEMENT_LINE: 
    distance = getSquaredElementDistance_Line(surfaceElement, currentpositions,physCoord, xsi, vector2minDistPoint);
    break;
  case DISTANCE_TO_ELEMENT_POINT:        
    distance = getSquaredElementDistance_Point(surfaceElement, currentpositions,physCoord, xsi, vector2minDistPoint);
    break;
  }
  return distance;
}  


double XFEM::getSquaredElementDistance_Surface(
    const DRT::Element*                     surfaceElement,
    const std::map<int,BlitzVec3>&          currentpositions,
    const BlitzVec3&                        X,
    BlitzVec2&                              xsi,
    BlitzVec3&                              vector2minDistPoint)
{
  double distance = 1.0e12;
  BlitzVec3 normal(0.0,0.0,0.0);

  // normal vector at position xsi
  BlitzVec3 eleNormalAtXsi;
  BlitzVec3 x_surface_phys;
  BlitzMat xyze_surfaceElement(getCurrentNodalPositions(surfaceElement, currentpositions));
  elementToCurrentCoordinates(surfaceElement, xyze_surfaceElement, xsi, x_surface_phys);
  // normal pointing away from the surface towards physCoord
  vector2minDistPoint(0) = X(0) - x_surface_phys(0);
  vector2minDistPoint(1) = X(1) - x_surface_phys(1);
  vector2minDistPoint(2) = X(2) - x_surface_phys(2);
  // absolute distance between point and surface
  distance = (vector2minDistPoint(0)*vector2minDistPoint(0) + vector2minDistPoint(1)*vector2minDistPoint(1) + vector2minDistPoint(2)*vector2minDistPoint(2));
  return distance; 
}

double XFEM::getSquaredElementDistance_Line(
    const DRT::Element*                     surfaceElement,
    const std::map<int,BlitzVec3>&          currentpositions,
    const BlitzVec3&                        X,
    BlitzVec2&                              xsi,
    BlitzVec3&                              vector2minDistPoint)
{
  double distance = 1.0e12;
  BlitzVec3 normal(0.0,0.0,0.0);

  // normal vector at position xsi
  BlitzVec3 eleNormalAtXsi;
  BlitzVec3 x_surface_phys;
  
  if (xsi(0)>=1) xsi(0)=1;
  if (xsi(0)<=(-1)) xsi(0)=(-1);
  if (xsi(1)>=1) xsi(1)=1;
  if (xsi(1)<=(-1)) xsi(1)=(-1);
  
  BlitzMat xyze_surfaceElement(getCurrentNodalPositions(surfaceElement, currentpositions));
  elementToCurrentCoordinates(surfaceElement, xyze_surfaceElement, xsi, x_surface_phys);
  // normal pointing away from the surface towards physCoord
  vector2minDistPoint(0) = X(0) - x_surface_phys(0);
  vector2minDistPoint(1) = X(1) - x_surface_phys(1);
  vector2minDistPoint(2) = X(2) - x_surface_phys(2);
  // absolute distance between point and surface
  distance = (vector2minDistPoint(0)*vector2minDistPoint(0) + vector2minDistPoint(1)*vector2minDistPoint(1) + vector2minDistPoint(2)*vector2minDistPoint(2));
  return distance; 
}

double XFEM::getSquaredElementDistance_Point(
    const DRT::Element*                     surfaceElement,
    const std::map<int,BlitzVec3>&          currentpositions,
    const BlitzVec3&                        X,
    BlitzVec2&                              xsi,
    BlitzVec3&                              vector2minDistPoint)
{
  double min_node_distance = 1.0e12;
  const DRT::Node* closestNode;
  const DRT::Node*const* nodes = surfaceElement->Nodes();  
  BlitzVec3 xNodePos;
  int numnode = surfaceElement->NumNode();
  for (int inode = 0; inode < numnode; ++inode)
  {
    BlitzVec3 vector;
    const DRT::Node* surfaceElementNode = nodes[inode]; 
    // node position in physical coordinates
    const BlitzVec3 x_node = currentpositions.find(surfaceElementNode->Id())->second;

    // vector pointing away from the node towards physCoord
    vector(0) = X(0) - x_node(0);
    vector(1) = X(1) - x_node(1);
    vector(2) = X(2) - x_node(2);
    // absolute distance between point and node
    const double distance = vector(0)*vector(0) + vector(1)*vector(1) + vector(2)*vector(2);

    if (distance < min_node_distance) {
      closestNode = surfaceElementNode;
      min_node_distance = distance;
      vector2minDistPoint=vector;
      xNodePos = x_node;
    }
  }     
  BlitzMat xyze_surfaceElement(getCurrentNodalPositions(surfaceElement, currentpositions));
  CurrentToSurfaceElementCoordinates(surfaceElement, xyze_surfaceElement, xNodePos, xsi);

  return min_node_distance; 
}


  
BlitzMat2x2 XFEM::getXsiBoundingBox(const DRT::Element* surfaceElement){
  BlitzMat2x2 xsiBB;
  xsiBB(0,0)=(-1);
  xsiBB(0,1)=1;
  xsiBB(1,0)=(-1);
  xsiBB(1,1)=1;
  return xsiBB;
}


BlitzMat3x2 XFEM::mergeAABB(const BlitzMat3x2& A, const BlitzMat3x2& B){
  BlitzMat3x2 MergedAABB;
  const int nsd=3;
  for(int dim=0; dim<nsd; dim++)
  {
    MergedAABB(dim, 0) = std::min( A(dim, 0),  B(dim, 0));
    MergedAABB(dim, 1) = std::max( A(dim, 1),  B(dim, 1));
  }
  return MergedAABB;
}

double XFEM::getOverlapArea(const BlitzMat3x2& A, const BlitzMat3x2& B, const BlitzMat3x2& C){
  const int nsd=3;
  double Area = 1;
  double edge = 0;
  for(int dim=0; dim<nsd; dim++)
  {
    edge = std::min( A(dim, 1), std::min( B(dim, 1), C(dim,1) ) )-std::max( A(dim, 0), std::max( B(dim, 0),C(dim,0) ) );
    if (edge<=0)
      return 0;
    Area = Area * edge;
  }
  return Area;
}

bool XFEM::isContainedAinB(const BlitzMat3x2& A, const BlitzMat3x2& B){
  const int nsd=3;
  for(int dim=0; dim<nsd; dim++){
    if ( (B(dim,0)<A(dim,0)) || (B(dim,1)>A(dim,1)) )
      return false;
  }
  return true;
}

double XFEM::getOverlapArea(const list<BlitzMat3x2 > AABBs){
  const int nsd=3;
  double Area = 1;
  double edge = 0;
  for(int dim=0; dim<nsd; dim++)
  {
    list<double> maxCoords;
    maxCoords.clear();
    list<double> minCoords;
    minCoords.clear();
    for(list<BlitzMat3x2 >::const_iterator myIt= AABBs.begin(); myIt != AABBs.end(); myIt++){
      maxCoords.push_back( (*myIt)(dim, 1) );
      minCoords.push_back( (*myIt)(dim, 0) );
    }
    edge = *min_element(maxCoords.begin(),maxCoords.end()) -*max_element(minCoords.begin(),minCoords.end());
    if (edge<=0)
      return 0;
    Area = Area * edge;
  }
  return Area;  
}

double XFEM::getVolume(const BlitzMat3x2& AABB){
  const int nsd=3;
  double A=1;
  for(int dim=0; dim<nsd; dim++)
  {
    A = A*(AABB(dim, 1)-AABB(dim,0));
  }
  return A;
}



void XFEM::checkRoughGeoType(
           DRT::Element*                element,
           const BlitzMat               xyze_element,
           EleGeoType&                  eleGeoType)
{
  const DRT::Element::DiscretizationType distype = element->Shape();
  
  if(DRT::UTILS::getOrder(distype) ==1)
    eleGeoType = LINEAR;
  else if(DRT::UTILS::getOrder(distype)==2)
    eleGeoType = HIGHERORDER;
  else
    dserror("order of element shapefuntion is not correct");
}

// gives maximum distance of a point from an AABB
double XFEM::getMaxDistanceFromAABB(const BlitzVec3& X, const BlitzMat3x2 AABB){
  double maxDistance2=0.0;
  double AABB_X2=0.0;
  double AABB_Y2=0.0;
  double AABB_Z2=0.0; 

  for(int x=0; x<2; x++)
  {
    AABB_X2 = X(0) - AABB(0,x);
    AABB_X2=AABB_X2*AABB_X2;
    for(int y=0; y<2; y++)
    {
      AABB_Y2=(X(1) - AABB(1,y))*(X(1) - AABB(0,y));
      for(int z=0; z<2; z++)
      {
        AABB_Z2=(X(2) - AABB(2,z))*(X(2) - AABB(0,z));
        maxDistance2 = max(maxDistance2, (AABB_X2+AABB_Y2+AABB_Z2) );
      }
    }
  } 
  return sqrt(maxDistance2); 
}

// gives AABB of an circle
BlitzMat3x2 XFEM::getAABBofSphere(const BlitzVec3& X, const double radius){
  const int nsd=3;
  BlitzMat3x2 AABB;
  for(int dim=0; dim<nsd; dim++)
  {
    AABB(dim, 1)= X(dim)+radius;
    AABB(dim, 0)= X(dim)-radius;
  }
  return AABB; 
}

const DRT::Node* XFEM::getNodeAtXsi(
    const DRT::Element*             surfaceElement, 
    const std::map<int,BlitzVec3>&  currentpositions,
    const BlitzVec2&                xsi)
{
  const double TOL = XFEM::TOL7;
  BlitzVec3 X;
  BlitzMat xyze_surfaceElement(getCurrentNodalPositions(surfaceElement, currentpositions));
  elementToCurrentCoordinates(surfaceElement, xyze_surfaceElement, xsi, X);
//  cout << "X:" << X(0)<<", "<<X(1)<<", "<<X(2)<< endl;
//  cout << "xsi:" << xsi(0)<<", "<<xsi(1)<< endl;
  const DRT::Node*const* nodes = surfaceElement->Nodes();  
  int numnode = surfaceElement->NumNode();
  for (int inode = 0; inode < numnode; ++inode)
  {
    BlitzVec3 vector;
    const DRT::Node* surfaceElementNode = nodes[inode]; 
    // node position in physical coordinates
    const BlitzVec3 x_node = currentpositions.find(surfaceElementNode->Id())->second;
//    cout << "X_node:" << x_node(0)<<", "<<x_node(1)<<", "<<x_node(2)<< endl;
    
    if ( (fabs(X(0)-x_node(0))<=TOL) && 
         (fabs(X(1)-x_node(1))<=TOL) &&
         (fabs(X(2)-x_node(2))<=TOL) ){
      return surfaceElementNode;
    }
  } 
  dserror ("there is no node at xsi ==> wrong tolerance?, bug??");
  return NULL; 
}

#endif  // #ifdef CCADISCRET
