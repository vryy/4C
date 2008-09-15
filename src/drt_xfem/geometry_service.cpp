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
#include "../drt_geometry/intersection_service.H"

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
    XAABB(dim, 0) = pos[dim] - GEO::TOL7;
    XAABB(dim, 1) = pos[dim] + GEO::TOL7;
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
    GEO::EleGeoType eleGeoType(GEO::HIGHERORDER);
    XFEM::checkRoughGeoType(ele, curr, eleGeoType);
    BlitzMat3x2 xaabbEle = GEO::computeFastXAABB(ele, curr, eleGeoType);
    XAABB = mergeAABB(XAABB, xaabbEle);
  }

  #ifdef DEBUG5
    cout << "_XAABB=" << XAABB(0,0) << ","<< XAABB(0,1) << "," << XAABB(1,0) << "," << XAABB(1,1) << "," << XAABB(2,0) << "," << XAABB(2,1) << endl;
  #endif
  return XAABB;
}

BlitzMat3x2 XFEM::getXAABBofDis(
    const DRT::Discretization& dis
    )
{
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
    const DRT::Discretization&        dis,
    const std::map<int,BlitzVec3>&    currentpositions,
    const list<const DRT::Element* >& ElementList,
    const BlitzVec3&                  X,
    double&                           dist)
{
  dist=1.0e12;
  const DRT::Element* closest_element;
  for (list< const DRT::Element* >::const_iterator myIt = ElementList.begin(); myIt != ElementList.end(); myIt++)
  {
    double distance = 1.0e12;
    const BlitzMat xyze(GEO::getCurrentNodalPositions(*myIt, currentpositions));
    GEO::EleGeoType eleGeoType(GEO::HIGHERORDER);
    checkRoughGeoType(*myIt, xyze, eleGeoType);
    distance = (getMaxDistanceFromAABB(X, computeFastXAABB(*myIt, xyze, eleGeoType)));
    if (distance < dist){
      closest_element = *myIt;
      dist = distance;
    }
  }
  return closest_element;
}

const DRT::Element* XFEM::nearestNeighbourInList(
    const DRT::Discretization&              dis, ///< discretization containing elements
    const std::map<int,BlitzVec3>&          currentpositions, ///< current position of elements in dis
    const list<const DRT::Element* >&       ElementList, ///<list of elements to be examined
    const BlitzVec3&                        X,  ///<coords of point which will be examined
    map<const DRT::Element*, double >&      squaredDistanceMap, ///<> map of already known distances (may be empty)
    BlitzVec3&                              vector2minDistPoint,  ///>vector from querypoint X to point on surface element with minimal distance
    double&                                 eleDistance, ///>distance btw. x and nearest element
    DistanceType&                           distanceType ///>distance Typ
    )
{
  double min_ele_distance = 1.0e12;
  const DRT::Element* closest_element=NULL;
  DistanceType closestElement_distanceType=ELEMENT_SURFACE;
  BlitzVec2 xsi;
  DistanceType tmpdistanceType=ELEMENT_SURFACE;
  bool nearestDistanceIsFromList=false;

  for (list< const DRT::Element* >::const_iterator myIt = ElementList.begin(); myIt != ElementList.end(); myIt++)
  {
    double distance = 1.0e12;
    const DRT::Element* cutterele = &**myIt;
    map<const DRT::Element*, double >::iterator distIterator=squaredDistanceMap.find(cutterele);
    BlitzVec2 eleCoord(0.0,0.0);
    BlitzVec3 tmpvector2minDistPoint(0.0,0.0,0.0);
    if (distIterator!=squaredDistanceMap.end())  // means that element is saved in list
    {
      distance = fabs(distIterator->second);
    }
    else{
      const BlitzMat xyze_cutter(GEO::getCurrentNodalPositions(cutterele, currentpositions));
      distance = XFEM::getSquaredElementDistance(cutterele,currentpositions,X,eleCoord,tmpvector2minDistPoint,tmpdistanceType);
      squaredDistanceMap[cutterele]=distance;
    }
    if (distance < min_ele_distance)
    {
      if (distIterator!=squaredDistanceMap.end())
        nearestDistanceIsFromList = true;
      else
        nearestDistanceIsFromList = false;
      closest_element = cutterele;
      closestElement_distanceType = tmpdistanceType;
      xsi = eleCoord;
      min_ele_distance = distance;
      vector2minDistPoint = tmpvector2minDistPoint;
    }

  }

  if (nearestDistanceIsFromList){
    const BlitzMat xyze_cutter(GEO::getCurrentNodalPositions(closest_element, currentpositions));
    min_ele_distance = XFEM::getSquaredElementDistance(closest_element,currentpositions,X,xsi,vector2minDistPoint,closestElement_distanceType);
  }
  distanceType = closestElement_distanceType;
  BlitzVec3 normal= getNormalAtXsi(dis, closest_element, currentpositions, xsi,X, closestElement_distanceType);
  const BlitzVec3 vectorX2minDistPoint((X(0)-vector2minDistPoint(0)),(X(1)-vector2minDistPoint(1)),(X(2)-vector2minDistPoint(2)));
  const double scalarproduct = vectorX2minDistPoint(0)*normal(0) + vectorX2minDistPoint(1)*normal(1) + vectorX2minDistPoint(2)*normal(2);
  min_ele_distance = sqrt(min_ele_distance);
  if (scalarproduct<0.0)
    min_ele_distance *= (-1);

  #ifdef DEBUG
  if (scalarproduct==0.0){
    cout << "scalarprod is 0" << endl;
  }
  #endif

  eleDistance =  min_ele_distance;
  return closest_element;
}

BlitzVec3 XFEM::getNormalAtXsi(
    const DRT::Discretization&        dis,
    const DRT::Element*               surfaceElement,
    const std::map<int,BlitzVec3>&    currentpositions,
    const BlitzVec2&                  xsi,
    const BlitzVec3&                  X,
    const DistanceType&               distanceType)
{
  BlitzVec3 normal(0.0,0.0,0.0);

  switch (distanceType) {
  case ELEMENT_SURFACE:{
    const BlitzMat xyze_surfaceElement = GEO::getCurrentNodalPositions(surfaceElement, currentpositions);
    GEO::computeNormalToSurfaceElement(surfaceElement, xyze_surfaceElement, xsi, normal);
    return normal;
    break;
  }
  case ELEMENT_LINE:
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
    const BlitzMat xyze_surfaceElement = GEO::getCurrentNodalPositions(surfaceElement, currentpositions);
    const DRT::Node* NodeA = getNodeAtXsi(surfaceElement, currentpositions, xsiA);
    const DRT::Node* NodeB = getNodeAtXsi(surfaceElement, currentpositions, xsiB);
    list< const DRT::Element* > commonElements = getCommonElements(NodeA, NodeB);
    list< BlitzVec3 > normalVectors;
    normalVectors.clear();
    if (commonElements.size() > 2 or commonElements.size() == 0){
      // note that 1 is ok if quasi-2d(thin 3d) problems are calculated, where in z-direction we have no closed surface
      // it would be also ok, for adaptivity with hanging nodes and hence hanging lines
      //dserror("a line should be bounded by 2 surface elements");
      }
    for(list< const DRT::Element* >::const_iterator myIt = commonElements.begin(); myIt != commonElements.end(); ++myIt ){
      const BlitzMat xyze_surfaceElement(GEO::getCurrentNodalPositions(*myIt, currentpositions));
      BlitzVec3 tmpNormal;
      BlitzVec2 tmpXsi;
      GEO::CurrentToSurfaceElementCoordinates((*myIt)->Shape(), xyze_surfaceElement, X, tmpXsi);
      GEO::computeNormalToSurfaceElement(*myIt, xyze_surfaceElement, tmpXsi, tmpNormal);
      normalVectors.push_back(tmpNormal);
    }
    return addVectors(normalVectors);
    break;
  }
  case ELEMENT_POINT:
  {
    bool printDEBUG = false;
    const BlitzMat xyze_surfaceElement = GEO::getCurrentNodalPositions(surfaceElement, currentpositions);
    const DRT::Node*  node = getNodeAtXsi(surfaceElement, currentpositions, xsi);
    const BlitzVec3 closest_node_pos = currentpositions.find(node->Id())->second;
    for(int j=0; j<node->NumElement();j++)
    {
      const DRT::Element* tmpSurfaceElement = node->Elements()[j];
      if (printDEBUG)
        cout << X(0) << "/" << X(1)<< "/" << X(2) <<" -> candID2: " << tmpSurfaceElement->Id() << endl;
      const BlitzMat xyze_surfaceElement(GEO::getCurrentNodalPositions(tmpSurfaceElement, currentpositions));
      BlitzVec3 eleNormalAtXsi;
      BlitzVec2 tmpXsi;
      GEO::CurrentToSurfaceElementCoordinates(tmpSurfaceElement->Shape(), xyze_surfaceElement, closest_node_pos, tmpXsi);
      GEO::computeNormalToSurfaceElement(tmpSurfaceElement, xyze_surfaceElement, tmpXsi, eleNormalAtXsi);
      normal(0) = normal(0) +  eleNormalAtXsi(0);
      normal(1) = normal(1) +  eleNormalAtXsi(1);
      normal(2) = normal(2) +  eleNormalAtXsi(2);
    }
    return normal;
    break;
  }
  }
  cout << "DistanceType: " << distanceType << endl;
  return normal;
  dserror("should not get here");
  return normal;
}

list < const DRT::Element* > XFEM::getCommonElements(
    const DRT::Node* A,
    const DRT::Node* B)
{
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

BlitzVec3 XFEM::addVectors(
    const list< BlitzVec3 > vectors
    )
{
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
    BlitzVec3&                              x_surface_phys,
    DistanceType&                           distanceType)
{
  double distance = 1.0e12;
  BlitzVec3 normal;
  const BlitzMat2x2 xsiBoundingBox = XFEM::getXsiBoundingBox(surfaceElement);
  const BlitzMat xyze_surfaceElement(GEO::getCurrentNodalPositions(surfaceElement, currentpositions));
  GEO::CurrentToSurfaceElementCoordinates(surfaceElement->Shape(), xyze_surfaceElement, physCoord, xsi);

  if ( (xsi(0) < xsiBoundingBox(0,0) || xsi(0) > xsiBoundingBox(0,1) ) &&
      (xsi(1) < xsiBoundingBox(1,0) || xsi(1) > xsiBoundingBox(1,1) )    )   {
    distanceType = ELEMENT_POINT;
//    cout<<" point type xsi:" <<  xsi(0)<<", "<<xsi(1)<< endl;
  }
  else {
    if ( xsi(0)>=xsiBoundingBox(0,0) && xsi(0)<=xsiBoundingBox(0,1) && xsi(1)>=xsiBoundingBox(1,0) && xsi(1)<=xsiBoundingBox(1,1) ){
      distanceType = ELEMENT_SURFACE;
//      cout<<" surface type xsi:" <<  xsi(0)<<", "<<xsi(1)<< endl;
    }
    else {
      distanceType = ELEMENT_LINE;
//      cout<<" line type xsi:" <<  xsi(0)<<", "<<xsi(1)<< endl;
    }
  }

  switch (distanceType) {
  case ELEMENT_SURFACE:
    distance = getSquaredElementDistance_Surface(surfaceElement, currentpositions,physCoord, xsi, x_surface_phys);
    break;
  case ELEMENT_LINE:
    distance = getSquaredElementDistance_Line(surfaceElement, currentpositions,physCoord, xsi, x_surface_phys);
    break;
  case ELEMENT_POINT:
    distance = getSquaredElementDistance_Point(surfaceElement, currentpositions,physCoord, xsi, x_surface_phys);
    break;
  }
  return distance;
}


double XFEM::getSquaredElementDistance_Surface(
    const DRT::Element*                     surfaceElement,
    const std::map<int,BlitzVec3>&          currentpositions,
    const BlitzVec3&                        X,
    BlitzVec2&                              xsi,
    BlitzVec3&                              x_surface_phys)
{
  double distance = 1.0e12;
  BlitzVec3 normal(0.0,0.0,0.0);
  BlitzVec3 vector2minDistPoint(0.0,0.0,0.0);

  // normal vector at position xsi
  BlitzVec3 eleNormalAtXsi;
  const BlitzMat xyze_surfaceElement(GEO::getCurrentNodalPositions(surfaceElement, currentpositions));
  GEO::elementToCurrentCoordinates(surfaceElement, xyze_surfaceElement, xsi, x_surface_phys);
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
    BlitzVec3&                              x_surface_phys)
{
  double distance = 1.0e12;
  BlitzVec3 normal(0.0,0.0,0.0);

  // normal vector at position xsi
  BlitzVec3 eleNormalAtXsi;
  BlitzVec3 vector2minDistPoint(0.0,0.0,0.0);

  if (xsi(0)>=1) xsi(0)=1;
  if (xsi(0)<=(-1)) xsi(0)=(-1);
  if (xsi(1)>=1) xsi(1)=1;
  if (xsi(1)<=(-1)) xsi(1)=(-1);

  const BlitzMat xyze_surfaceElement(GEO::getCurrentNodalPositions(surfaceElement, currentpositions));
  GEO::elementToCurrentCoordinates(surfaceElement, xyze_surfaceElement, xsi, x_surface_phys);
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
    BlitzVec3&                              x_surface_phys)
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
      xNodePos = x_node;
    }
  }
  const BlitzMat xyze_surfaceElement(GEO::getCurrentNodalPositions(surfaceElement, currentpositions));
  GEO::CurrentToSurfaceElementCoordinates(surfaceElement->Shape(), xyze_surfaceElement, xNodePos, xsi);
  x_surface_phys = xNodePos;
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

bool XFEM::isContainedXinAABB(
    const BlitzMat3x2& AABB,
    const BlitzVec3& X)
{
  const int nsd=3;
  for(int dim=0; dim<nsd; dim++){
    if ( (X(dim)<AABB(dim,0)) || (X(dim)>AABB(dim,1)) )
      return false;
  }
  return true;
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

double XFEM::getVolume(
    const BlitzMat3x2& AABB
    ){
  const int nsd=3;
  double A=1;
  for(int dim=0; dim<nsd; dim++)
  {
    A = A*(AABB(dim, 1)-AABB(dim,0));
  }
  return A;
}



void XFEM::checkRoughGeoType(
           const DRT::Element*          element,
           const BlitzMat               xyze_element,
           GEO::EleGeoType&             eleGeoType)
{
  const DRT::Element::DiscretizationType distype = element->Shape();

  if(DRT::UTILS::getOrder(distype) ==1)
    eleGeoType = GEO::LINEAR;
  else if(DRT::UTILS::getOrder(distype)==2)
    eleGeoType = GEO::HIGHERORDER;
  else
    dserror("order of element shapefuntion is not correct");
}

//! gives maximum distance of a point from an AABB
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

//! gives AABB of an circle
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
  const double TOL = GEO::TOL7;
  BlitzVec3 X;
  const BlitzMat xyze_surfaceElement(GEO::getCurrentNodalPositions(surfaceElement, currentpositions));
  GEO::elementToCurrentCoordinates(surfaceElement, xyze_surfaceElement, xsi, X);
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

double XFEM::biggestRadiusInAABB(
    const BlitzMat3x2&                      AABB, ///<bounding box
    const BlitzVec3&                        X //point coords
    )
{
  double radius=0.0;
  for (int i=0; i<3; i++){
    radius = min(radius,fabs(X(i)-AABB(i,0)));
    radius = min(radius,fabs(X(i)-AABB(i,1)));
  }
  return radius;
}


#endif  // #ifdef CCADISCRET
