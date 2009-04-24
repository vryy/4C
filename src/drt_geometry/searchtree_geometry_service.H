/*!
\file searchtree_geometry_service.H

\brief provides geometry methods for a search tree

<pre>
-------------------------------------------------------------------------
                 BACI finite element library subsystem
            Copyright (2008) Technical University of Munich
              
Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed, 
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions and of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library may solemnly used in conjunction with the BACI contact library
for purposes described in the above mentioned contract.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de) 
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de                   

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
 */

#ifdef CCADISCRET
#ifndef SEARCHTREE_GEOMETRY_SERVICE_H_
#define SEARCHTREE_GEOMETRY_SERVICE_H_

#include "../drt_geometry/element_coordtrafo.H"
#include "../drt_geometry/intersection_service_templates.H"
#include "../drt_geometry/searchtree_nearestobject.H"
#include "../drt_geometry/intersection_interfacepoint.H"

namespace GEO
{

  /*!
  \brief Returns the eXtendedAxisAlignedBoundingBox of a discretization, 
         also takes into account current displacements
  \param dis                  discretization 
  \param currentpositions     current nodal positions in discretization
  \return  XAxisAlignedBoundingBox 
   */
  LINALG::Matrix<3,2> getXAABBofDis(
    const DRT::Discretization&                  dis,               
    const std::map<int,LINALG::Matrix<3,1> >&   currentpositions  
    );

  
  /*!
  \brief Returns the eXtendedAxisAlignedBoundingBox of a discretization, 
         in reference configuration
  \param dis                  discretization 
  \return  XAxisAlignedBoundingBox 
   */
  LINALG::Matrix<3,2> getXAABBofDis(
      const DRT::Discretization& dis
      );
  
  
  /*!
  \brief Returns a vector of eXtendedAxisAlignedBoundingBox for labeled structures 
  \param dis                  discretization 
  \param currentpositions     current nodal positions in discretization
  \param elementList          map of int = label and set of element ids
  \return  AxisAlignedBoundingBox as BlitzMat3x2
   */
  std::vector< LINALG::Matrix<3,2> > computeXAABBForLabeledStructures(
      const DRT::Discretization&                  dis,
      const std::map<int,LINALG::Matrix<3,1> >&   currentpositions,
      const std::map< int, std::set<int> >&       elementList
      );
  
  
  /*!
  \brief Returns the label for a point and a given map of elements
  \param dis                  discretization
  \param currentpositions     current nodal positions in discretization
  \param querypoint           point to be examined
  \param elementList          map of elements label gid
  \return                     label
   */
  int getXFEMLabel(
      const DRT::Discretization&                  dis, 		
      const std::map<int,LINALG::Matrix<3,1> >&   currentpositions, 	
      const LINALG::Matrix<3,1>&                  querypoint,
      const std::map<int, std::set<int> >&        elementList
      ); 
      
    
    /*!
  \brief Returns the label for a point and a given map of elements
  \param triangleList         triangleList
  \param currentpositions     current nodal positions in discretization
  \param querypoint           point to be examined
  \param elementList          map of elements label gid
  \return                     label
   */  
  int getXFEMLabel(
    const std::vector<std::vector<int> >&         triangleList,    
    const std::vector<GEO::InterfacePoint>&       pointList,
    const LINALG::Matrix<3,1>&                    querypoint,
    const std::map<int, std::set<int> >&          elementList
    );       	
  
  
  /*!
  \brief Returns the label for a point and a given map of elements
  \param triangleList         list of triangles
  \param interfacePointList   list of interface points
  \param currentpositions     current nodal positions in discretization
  \param querypoint           point to be examined
  \param elementList          map of elements label gid
  \return                     label
   */
  void moveNodeOutOfStructure(
    const vector<vector<int> >&                   triangleList,
    vector< GEO::InterfacePoint >&                pointList,
    const int                                     querypointId,
    const std::map<int, std::set<int> >&          elementList
    ); 
  
  
  /*!
  \brief Returns the label and the nearest object for a point and a given map of elements
  \param dis                  discretization
  \param currentpositions     current nodal positions in discretization
  \param querypoint           point to be examined
  \param elementList          map of elements label gid
  \param nearestObject        nearestObject  
  \return                     label
   */
  int getXFEMLabelAndNearestObject(
      const DRT::Discretization&                dis,              
      const std::map<int,LINALG::Matrix<3,1> >& currentpositions, 
      const LINALG::Matrix<3,1>&                querypoint,      
      const std::map<int, std::set<int> >&      elementList,     
      GEO::NearestObject&                       nearestObject
      );
  
  
  /*!
  \brief Returns a map of element ids lying in a circle of a given querypoint and a radius
  \param dis                  discretization
  \param currentpositions     current nodal positions in discretization
  \param querypoint           point to be examined
  \param radius               radius
  \param label                label
  \param elementList          map of elements label gid
  \return                     label
   */
  std::map<int,std::set<int> > getElementsInRadius(
      const DRT::Discretization&                  dis,    
      const std::map<int,LINALG::Matrix<3,1> >&   currentpositions,   
      const LINALG::Matrix<3,1>&                  querypoint,
      const double                                radius,
      const int                                   label,
      std::map<int, std::set<int> >&              elementList
      );
  
  /*!
  \brief returns vector with GIDs of nodes within given radius around querypoint
  \param dis                  discretization
  \param currentpositions     current nodal positions in discretization
  \param querypoint           point to be examined
  \param radius               radius
  \param elementList          map of elements label gid
  \return                     label
   */
  std::vector<int> getPointsInRadius(   
      const std::map<int,LINALG::Matrix<3,1> >&   currentpositions,   
      const LINALG::Matrix<3,1>&                  querypoint,
      const double                                radius,
      std::map<int, std::set<int> >&              elementList
      );
  
  
  /*!
  \brief Returns a vector of element ids intersecting with a given query element (XAABB)
         This is for CONTACT search!
  \param dis                  discretization
  \param currentpositions     current nodal positions in discretization
  \param element              element to be examined
  \param elementList          map of elements label gid
  \return                     label
   */
  void searchMultibodyContactElements(
    const std::map<int,LINALG::Matrix<9,2> >&   currentKDOPs, 
    const LINALG::Matrix<9,2>&                  queryKDOP,
    const int                                   label,
    const std::map<int, std::set<int> >&        elementList,
    std::set<int>&                              contactEleIds
    );  
  
  
  /*!
  \brief Computes a 2D extended axis-aligned bounding box for an element (XAABB)
         This is for CONTACT search!
  \param element        (in)  element
  \param xyze           (in)  nodal position array (3,numnode)
  \return extended axis-aligned bounding box (XAABB) as BlitzMat3x2
   */
  LINALG::Matrix<3,2> computeContactXAABB( 
      DRT::Element*       element,
      const LINALG::SerialDenseMatrix&     xyze);
  
  
  /*!
  \brief Returns the nearest object to a point 
  \param dis                  discretization
  \param currentpositions     current nodal positions in discretization
  \param elementList          element list
  \param point                point to be examined
  \param x_surface_phys       basis point in current coordinates
  \param nearestObject        info about nearest object
  \return label of nearest object
   */
  int nearestObjectInNode(
      const DRT::Discretization&                  dis,  
      const std::map<int,LINALG::Matrix<3,1> >&   currentpositions,
      const std::map<int, std::set<int> >&        elementList,
      const LINALG::Matrix<3,1>&                  point,
      LINALG::Matrix<3,1>&                        minDistanceVec,
      GEO::NearestObject&                         nearestObject
      );
  
  
  /*!
  \brief Returns the nearest object to a point 
  \param triangleList         list of linear triangular elements
  \param currentpositions     current nodal positions in discretization
  \param elementList          element list
  \param point                point to be examined
  \param x_surface_phys       basis point in current coordinates
  \param nearestObject        info about nearest object
  \return label of nearest object
   */
  int nearestObjectInNode(
    const vector<vector<int> >&                 triangleList,
    const std::vector<GEO::InterfacePoint>&     pointList,
    const std::map<int, std::set<int> >&        elementList,
    const LINALG::Matrix<3,1>&                  point,
    LINALG::Matrix<3,1>&                        minDistanceVec,
    GEO::NearestObject&                         nearestObject
    );


  /*!
  \brief Searches for the nearest surface element to a given point
  \param surfaceElement       surface element
  \param currentpositions     current nodal positions in discretization
  \param point                point to be examined
  \param x_surface_phys       physical coordinates 
  \param distance             distance to nearest element 
  \return true if nearest surface element found
   */
  bool getDistanceToSurface(
      const DRT::Element*                         surfaceElement,
      const std::map<int,LINALG::Matrix<3,1> >&   currentpositions,
      const LINALG::Matrix<3,1>&                  point,
      LINALG::Matrix<3,1>&                        x_surface_phys,
      double&                                     distance
      );


  /*!
  \brief Searches for the nearest linear triangular element to a given point
  \param triElement           triangular element
  \param currentpositions     current nodal positions in discretization
  \param point                point to be examined
  \param x_surface_phys       physical coordinates 
  \param distance             distance to nearest element 
  \return true if nearest surface element found
   */
  bool getDistanceToSurface(
      const std::vector<int>&                     triElement,
      const std::vector<GEO::InterfacePoint>&     pointList,
      const LINALG::Matrix<3,1>&                  point,
      LINALG::Matrix<3,1>&                        x_surface_phys,
      double&                                     distance
      );
  
  
  /*!
  \brief Searches for the nearest line element to a given point
  \param lineElement          lineElement
  \param currentpositions     current nodal positions in discretization
  \param point                point to be examined
  \param x_line_phys          physical coordinates 
  \param distance             distance to nearest point
  \return true if nearest line element found
  */
  bool getDistanceToLine(
      const DRT::Element*                         lineElement,
      const std::map<int,LINALG::Matrix<3,1> >&   currentpositions,
      const LINALG::Matrix<3,1>&                  point,
      LINALG::Matrix<3,1>&                        x_line_phys,
      double&                                     distance
      );
  
  
  /*!
  \brief Searches for the nearest line element to a given point
  \param lineElement          lineElement
  \param currentpositions     current nodal positions in discretization
  \param point                point to be examined
  \param x_line_phys          physical coordinates 
  \param distance             distance to nearest point
  \return true if nearest line element found
  */
  bool getDistanceToLine(
      const vector<int>&                          lineElement,
      const std::vector<GEO::InterfacePoint>&     pointList,
      const LINALG::Matrix<3,1>&                  point,
      LINALG::Matrix<3,1>&                        x_line_phys,
      double&                                     distance
      );
  
  
  /*!
  \brief Searches for the nearest node to a given point
  \param nodeId               nodeId
  \param currentpositions     current nodal positions in discretization
  \param point                point to be examined
  \param distance             distance to nearest point
   */
  void getDistanceToPoint(
      const int                                   nodeId,
      const std::vector<GEO::InterfacePoint>&     pointList,
      const LINALG::Matrix<3,1>&                  point,
      double&                                     distance
      );
       
  
  /*!
  \brief Searches for the nearest node to a given point
  \param node                 node
  \param currentpositions     current nodal positions in discretization
  \param point                point to be examined
  \param distance             distance to nearest point
   */    
  void getDistanceToPoint(
    const DRT::Node*                            node,
    const std::map<int,LINALG::Matrix<3,1> >&   currentpositions,
    const LINALG::Matrix<3,1>&                  point,
    double&                                     distance
    );


  /*!
  \brief  Returns the normal to a surface element or a an average normal 
          of all surface elements adjacent to a line or node
  \param dis                  discretization
  \param currentpositions     current nodal positions in discretization
  \param nearestObject        info about nearest object
  \return normal vector in physical coordinates 
   */
  LINALG::Matrix<3,1> getNormalAtSurfacePoint(
      const DRT::Discretization&                  dis,
      const std::map<int,LINALG::Matrix<3,1> >&   currentpositions,
      GEO::NearestObject&                         nearestObject
      );
      
      
    /*!
  \brief  Returns the normal to a surface element or a an average normal 
          of all surface elements adjacent to a line or node
  \param triangleList         triangleList 
  \param currentpositions     current nodal positions in discretization
  \param nearestObject        info about nearest object
  \return normal vector in physical coordinates 
   */
  LINALG::Matrix<3,1> getNormalAtSurfacePoint(
    const std::vector< std::vector<int> >&      triangleList,
    const std::vector<GEO::InterfacePoint>&     pointList,
    GEO::NearestObject&                         nearestObject
    );

  
  /*!
  \brief  Finds all elements adjacent to a line described by two end points
  \param node1         node 1
  \param node2         node 2
  \return list of common elements
   */
  std::vector <int> getAdjacentSurfaceElementsToLine(
      const DRT::Node*                        node1, 
      const DRT::Node*                        node2  
      );
      
  
  /*!
  \brief  Finds all elements adjacent to a line described by two end points
  \param node1         node 1
  \param node2         node 2
  \return list of common elements
   */    
  std::vector<int> getAdjacentTriElementsToLine(
    const std::vector< vector<int> >& triangleList,
    const int                         triangleId,
    const int                         node1,
    const int                         node2
    );
    
    
  /*!
  \brief  Finds all elements adjacent to a node
  \param triangleList       list of trianglea
  \param triangleId         triangle id
  \param node               node id
  \return list of common elements
   */  
  std::vector<int> getAdjacentTriElementsToNode(
    const std::vector< vector<int> >& triangleList,
    const int                         triangleId,
    const int                         node
    );

  
  /*!
  \brief  Checks if a point lies with in a node box
  \param point         query point
  \param nodeBox       box of tree node
  \return true, if a point lies with in a node box, false otherwise 
   */
  bool  pointInTreeNode(
      const LINALG::Matrix<3,1>&              point,
      const LINALG::Matrix<3,2>&            nodeBox
      );

  
  /*!
  \brief  Checks if a point , which classifies a nearest object, lies within a circle (radius is equal 
          minimum distance between the querypoint and a node box wall)
  \param nearestpoint         nearestpoint
  \param querypoint           query point
  \param nodeBox              nodeBox
  \param rootNode             rootNode
  \return true if nearest point lies within a circle around the query point to a wall of the tree node box with minimum distance
   */
  bool  pointInMinCircleInTreeNode(
    const LINALG::Matrix<3,1>&                nearestpoint,
    const LINALG::Matrix<3,1>&                querypoint,
    const LINALG::Matrix<3,2>&                nodeBox,
    const bool                                rootNode
    ); 
  
  
  /*!
  \brief Merges two axis-aligned bounding boxes
  \param AABB1          axis - aligned bounding box 1
  \param AABB2          axis - aligned bounding box 2
  \return merged axis-aligned bounding box
   */
  LINALG::Matrix<3,2> mergeAABB(
      const LINALG::Matrix<3,2>&                      AABB1, 
      const LINALG::Matrix<3,2>&                      AABB2  
      );

  
  /*!
  \brief Tests if  AABB_old AND AABB_new are completely contained in node box
  \param AABB_old             AABB_old
  \param AABB_new             AABB_new
  \param nodeBox              AABB (e.g. node box)
  \return true if both AABB's are in node box
   */
  bool inSameNodeBox(
      const LINALG::Matrix<3,2>&  AABB_old, 
      const LINALG::Matrix<3,2>&  AABB_new,
      const LINALG::Matrix<3,2>&  nodeBox
      );

  
  /*!
  \brief Check the type of geometry of an element, 
         rough because only linear or nonlinear is checked
  \param element              element
  \param xyze_element         nodal coordinates 
  \param eleGeoType           geometry type
   */ 
  void checkRoughGeoType(
      const DRT::Element*                     element,      
      const LINALG::SerialDenseMatrix         xyze_element, 
      GEO::EleGeoType&                        eleGeoType  
      );
  
  
  /*!
  \brief Returns a vector of intersection candidate ids
  \param dis                  discretization
  \param currentpositions     current nodal positions in discretization
  \param xfemElement          element to be examined
  \param elementList          map of elements label gid
  \return                     label
   */
  void getIntersectionCandidates(      
      const std::map<int,LINALG::Matrix<3,2> >&   currentXAABBs, 
      const LINALG::Matrix<3,2>&                  xfemXAABB,
      const std::map<int, std::set<int> >&        elementList,
      std::set<int>&                              intersectionCandidateIds
      );


  
}

#endif
#endif


