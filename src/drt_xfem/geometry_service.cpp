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



BlitzMat3x2 XFEM::getXAABBofDis(const DRT::Discretization& dis){
  const int nsd = 3;
  BlitzMat3x2 XAABB;
  
  // XAABB of cutterElements
  // first node
  const double* pos = dis.lRowElement(0)->Nodes()[0]->X();
  for(int dim=0; dim<nsd; ++dim)
  {
    XAABB(dim, 0) = pos[dim] - TOL7;
    XAABB(dim, 1) = pos[dim] + TOL7;
  }
  for (int j=0; j< dis.NumMyRowElements(); ++j) {
    // remaining node
    for(int i=0; i< dis.lRowElement(j)->NumNode(); ++i)
    {
      const double* posEle = dis.lRowElement(j)->Nodes()[i]->X();
      for(int dim=0; dim<nsd; dim++)
      {
        XAABB(dim, 0) = std::min( XAABB(dim, 0), posEle[dim] - TOL7);
        XAABB(dim, 1) = std::max( XAABB(dim, 1), posEle[dim] + TOL7);
      }
    }  
  }
   XAABB(0,0) = XAABB(0,0) - TOL7;
   XAABB(0,1) = XAABB(0,1) + TOL7;
   XAABB(1,0) = XAABB(1,0) - TOL7;
   XAABB(1,1) = XAABB(1,1) + TOL7;
   XAABB(2,0) = XAABB(2,0) - TOL7;
   XAABB(2,1) = XAABB(2,1) + TOL7;
   
   //  cout << "_XAABB=" << XAABB(0,0) << ","<< XAABB(0,1) << "," << XAABB(1,0) << "," << XAABB(1,1) << "," << XAABB(2,0) << "," << XAABB(2,1) << endl;
   return XAABB;
}
  
BlitzMat3x2 XFEM::getXAABBofDis(const DRT::Discretization& dis,const std::map<int,BlitzVec3>& currentpositions){
  const int nsd = 3;
  BlitzMat3x2 XAABB;
  if (dis.NumGlobalElements() == 0){ 
    XAABB(0,0)=0;XAABB(0,1)=0;
    XAABB(1,0)=0;XAABB(1,1)=0;
    XAABB(2,0)=0;XAABB(2,1)=0;
    return XAABB;
    }
  // extend XAABB to xfem elements
  const double* pos = dis.lRowElement(0)->Nodes()[0]->X();
  for(int dim=0; dim<nsd; ++dim)
  {
    XAABB(dim, 0) = pos[dim] - TOL7;
    XAABB(dim, 1) = pos[dim] + TOL7;
  }
  for (int j=0; j< dis.NumMyRowElements(); ++j) {
    // remaining node
    for(int i=0; i< dis.lRowElement(j)->NumNode(); ++i)
    {
      const double* posEle = dis.lRowElement(j)->Nodes()[i]->X();
      for(int dim=0; dim<nsd; dim++)
      {
        XAABB(dim, 0) = std::min( XAABB(dim, 0), posEle[dim] - TOL7);
        XAABB(dim, 1) = std::max( XAABB(dim, 1), posEle[dim] + TOL7);
      }
    }  
  }
  XAABB(0,0) = XAABB(0,0) - TOL7;
  XAABB(0,1) = XAABB(0,1) + TOL7;
  XAABB(1,0) = XAABB(1,0) - TOL7;
  XAABB(1,1) = XAABB(1,1) + TOL7;
  XAABB(2,0) = XAABB(2,0) - TOL7;
  XAABB(2,1) = XAABB(2,1) + TOL7;
  
  //  cout << "_XAABB=" << XAABB(0,0) << ","<< XAABB(0,1) << "," << XAABB(1,0) << "," << XAABB(1,1) << "," << XAABB(2,0) << "," << XAABB(2,1) << endl;
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
    //      printf("BRANCHING to node-distance ************************************************** \n");
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
    //        if (tmp!=closest_node->Id())
    //          printf("found nearest node with id %d (%f,%f,%f) for x_in (%f,%f,%f)", closest_node->Id(), closest_node->X()[0], closest_node->X()[1], closest_node->X()[2], x_in(0),x_in(1),x_in(2));
    // getElements corseponding to node
    // calculate normal at node by (vectorial) addition of element normals
    BlitzVec3 normal;
    normal(0)=0;normal(1)=0;normal(2)=0;
    DRT::Node*  node = dis.gNode(closest_node->Id());
    for(int j=0; j<node->NumElement();j++)
    {
      DRT::Element* surfaceElement = node->Elements()[j];
      BlitzMat xyze_surfaceElement(getCurrentNodalPositions(surfaceElement, currentpositions));
      BlitzVec3 eleNormalAtXsi;
      BlitzVec2 xsi;
      CurrentToSurfaceElementCoordinates(surfaceElement, xyze_surfaceElement, node->X(), xsi);
      //          printf("xsi (%f,%f)\n", xsi(0), xsi(1));
      // normal vector at position xsi
      computeNormalToBoundaryElement(surfaceElement, xyze_surfaceElement, xsi, eleNormalAtXsi);
      normal(0) = normal(0) +  eleNormalAtXsi(0);
      normal(1) = normal(1) +  eleNormalAtXsi(1);
      normal(2) = normal(2) +  eleNormalAtXsi(2);
    }
    //      double norm = sqrt(normal(0)*normal(0)+normal(1)*normal(1)+normal(2)*normal(2));
    //        if (tmp!=closest_node->Id())
    //          printf(", normal(%f,%f,%f)", normal(0),normal(1),normal(2));
    
    //        normal(0) = normal(0)/norm;
    //        normal(1) = normal(1)/norm;
    //        normal(2) = normal(2)/norm;
    //        printf("normed normal at node (%f,%f,%f)\n", normal(0),normal(1),normal(2));
    closest_element=node->Elements()[0];
    // compute distance with sign
    //        if (tmp!=closest_node->Id())
    //              printf(", x2N: %f,%f,%f ", vectorX2minNode(0),vectorX2minNode(1),vectorX2minNode(2));
    const double scalarproduct = vectorX2minNode(0)*normal(0) + vectorX2minNode(1)*normal(1) + vectorX2minNode(2)*normal(2);
    const double vorzeichen = scalarproduct/abs(scalarproduct);
    //        if (tmp!=closest_node->Id())
    //              printf(", vz: %f\n", vorzeichen);
    min_ele_distance *= vorzeichen;
    dist = min_ele_distance;
    //        printf("end BRANCHING to node-distance ********************************************** \n");
  }
  return closest_element;
}


#endif  // #ifdef CCADISCRET
