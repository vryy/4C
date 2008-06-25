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


const DRT::Element* XFEM::nearestNeighbourInList(const DRT::Discretization& dis,const std::map<int,BlitzVec3>& currentpositions, const list<const DRT::Element* >& ElementList, const BlitzVec3& x_in, double& dist)
{
  bool in_element = false;
  double min_ele_distance = 1.0e12;
  BlitzVec3 vectorX2minNode;
  const DRT::Element* closest_element;
  const DRT::Node* closest_node;
  bool foundNearSurfaceElement =false;
  
  for (list< const DRT::Element* >::const_iterator myIt = ElementList.begin(); myIt != ElementList.end(); myIt++)
  {
    double distance = 1.0e12;
    const DRT::Element* cutterele = &**myIt;
    const BlitzMat xyze_cutter(getCurrentNodalPositions(cutterele, currentpositions));
    static BlitzVec2 eleCoord;
    static BlitzVec3 normal;
    in_element = XFEM::searchForNearestPointOnSurface(cutterele,xyze_cutter,x_in,eleCoord,normal,distance);
    if (in_element && (abs(distance) < abs(min_ele_distance)))
    {
      closest_element = cutterele;
      min_ele_distance = distance;
      foundNearSurfaceElement =true;
    }
  }
  if (foundNearSurfaceElement)  
  {
    dist = min_ele_distance;
  }
  else 
  { 
    min_ele_distance = 1.0e12;
    for (list< const DRT::Element* >::const_iterator myIt2 = ElementList.begin(); myIt2 != ElementList.end(); myIt2++) {
      double distance = 1.0e12;
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
        distance = sqrt(vector(0)*vector(0) + vector(1)*vector(1) + vector(2)*vector(2));
        
        if (distance < min_ele_distance) {
          closest_node = cutternode;
          min_ele_distance = distance;
          vectorX2minNode=vector;
        }
      }          
    }
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
      normal(0) = normal(0) +  eleNormalAtXsi(0);
      normal(1) = normal(1) +  eleNormalAtXsi(1);
      normal(2) = normal(2) +  eleNormalAtXsi(2);
    }
    closest_element=node->Elements()[0];
    const double scalarproduct = vectorX2minNode(0)*normal(0) + vectorX2minNode(1)*normal(1) + vectorX2minNode(2)*normal(2);
    const double vorzeichen = scalarproduct/abs(scalarproduct);
    min_ele_distance *= vorzeichen;
    dist = min_ele_distance;
  }
  return closest_element;
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


double XFEM::getArea(const BlitzMat3x2& AABB){
  const int nsd=3;
  double A=1;
  for(int dim=0; dim<nsd; dim++)
  {
    A = A*(AABB(dim, 1)-AABB(dim,0));
  }
  return A;
}

#endif  // #ifdef CCADISCRET
