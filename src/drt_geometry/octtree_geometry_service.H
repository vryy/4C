/*!
\file geometry_service.H

\brief provides service functions for algorithmic geometry

<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
 */

#ifdef CCADISCRET
#ifndef OCTTREE_GEOMETRY_SERVICE_H_
#define OCTTREE_GEOMETRY_SERVICE_H_

#include <blitz/array.h>
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_element.H"
#include "../drt_xfem/intersection_service.H"
#include "../drt_xfem/blitz_tiny_operation.H"
#include "../drt_xfem/xfem.H"
#include "Teuchos_RefCountPtr.hpp"

namespace GEO
{


//! this proporty decides which distance algorithm is used
enum DistanceType {ELEMENT_SURFACE, ELEMENT_LINE, ELEMENT_POINT};

  /*!
  \brief Returns the eXtendedAxisAlignedBoundingBox of a complete Discretization and adds TOL07, 
         also takes into account current displacements
  \param dis                  discretization 
  \param currentpositions     current nodal positions in discretization
  \return  AxisAlignedBoundingBox as BlitzMat3x2
   */
  BlitzMat3x2 getXAABBofDis(
    const DRT::Discretization&            dis,               
    const std::map<int,BlitzVec3>&        currentpositions  
    );

  
  /*!
  \brief Returns a vector of eXtendedAxisAlignedBoundingBox for labeled structures 
  \param dis                  discretization 
  \param currentpositions     current nodal positions in discretization
  \param elementList          map of int = label and set of element ids
  \return  AxisAlignedBoundingBox as BlitzMat3x2
   */
  std::vector< BlitzMat3x2 > computeXAABBForLabeledStructures(
      const DRT::Discretization&              dis,
      const std::map<int,BlitzVec3>&          currentpositions,
      const std::map< int, std::set<int> >&   elementList
      );
  
  
  /*!
  \brief get a label for a point and a given map of elements
  \param dis                  discretization
  \param currentpositions     current nodal positions in discretization
  \param point                point to be examined
  \param elementList          map of elements label gid
  \return                     label
   */
  int getXFEMLabel(
      const DRT::Discretization&              dis, 		
      const std::map<int,BlitzVec3>&          currentpositions, 	
      const BlitzVec3&                        point,
      std::map<int, std::set<int> >&          elementList
      );  	


  /*!
  \brief  returns the nearest object to a point 
  \param dis                  discretization
  \param currentpositions     current nodal positions in discretization
  \param elementList          element list
  \param point                point to be examined
  \param xsi                  basis point in element coordinates
  \param x_surface_phys       basis point in current coordinates
  \param closestObject        info about closest object
  \return label of nearest object
   */
  int nearestObjectInNode(
      const DRT::Discretization&                  dis,  
      const std::map<int,BlitzVec3>&              currentpositions,
      const std::map<int, std::set<int> >&        elementList,
      const BlitzVec3&                            point,
      BlitzVec2&                                  xsi,
      BlitzVec3&                                  minDistanceVec,
      std::map<DistanceType, std::vector<int> >&  closestObject
      );


  /*!
  \brief  returns the distance of a point to an surface element 
  \param surfaceElement       surface element
  \param currentpositions     current nodal positions in discretization
  \param point                point to be examined
  \param xsi                  basis point in element coordinates
  \param x_surface_phys       basis point in current coordinates
  \param closestObject        info about closest object
  \param distance             distance to closest element if found other wise -1
  \return true if closest surface element found
   */
  bool getDistanceToSurface(
      const DRT::Element*                         surfaceElement,
      const std::map<int,BlitzVec3>&              currentpositions,
      const BlitzVec3&                            point,
      BlitzVec2&                                  xsi,
      BlitzVec3&                                  x_surface_phys,
      std::map<DistanceType, std::vector<int> >&  closestObject,
      double&                                     distance
      );


  /*!
  \brief finds the closest line belonging to an element of the element list to a given point
  \param dis                  discretization
  \param currentpositions     current nodal positions in discretization
  \param point                point to be examined
  \param xsi                  basis point in element coordinates
  \param x_line_phys          basis point in element coordinates
  \param closestObject        info about closest object
  \param distance             distance to closest point
  */
  bool getDistanceToLine(
      const DRT::Element*                         lineElement,
      const std::map<int,BlitzVec3>&              currentpositions,
      const BlitzVec3&                            point,
      BlitzVec2&                                  xsi,
      BlitzVec3&                                  x_line_phys,
      std::map<DistanceType, std::vector<int> >&  closestObject,
      double&                                     distance
      );

  
  /*!
  \brief finds the closest node belonging to an element of the element list to a given point
  \param dis                  discretization
  \param currentpositions     current nodal positions in discretization
  \param point                point to be examined
  \param xsi                  basis point in element coordinates
  \param x_node               physical coordinates of tree node
  \param closestObject        info about closest object
  \param distance             distance to closest point
   */
  void getDistanceToPoint(
      const DRT::Node*                            node,
      const std::map<int,BlitzVec3>&              currentpositions,
      const BlitzVec3&                            point,
      BlitzVec2&                                  xsi,
      BlitzVec3&                                  x_node,
      std::map<DistanceType, std::vector<int> >&  closestObject,
      double&                                     distance
      );


  /*!
  \brief  gives back the normal to an surface element at position xsi 
  \param dis                  discretization
  \param currentpositions     current nodal positions in discretization
  \param point                point to be examined
  \param xsi                  basis point in element coordinates
  \param closestObject        info about closest object
  \return normal vector in physical coordinates 
   */
  BlitzVec3 getNormalAtXsi(
      const DRT::Discretization&                  dis,
      const std::map<int,BlitzVec3>&              currentpositions,
      const BlitzVec3&                            point,
      const BlitzVec2&                            xsi,
      std::map<DistanceType, std::vector<int> >&  closestObject
      );

  
  /*!
  \brief  find all elements adjacent to a line described by two end points
  \param node1         node 1
  \param node2         node 2
  \return list of common elements
   */
  std::vector <int> getAdjacentSurfaceElementsToLine(
      const DRT::Node*                        node1, 
      const DRT::Node*                        node2  
      );

  
  /*!
  \brief merges 2 bounding boxes
  \param AABB1          axis - aligned bounding box 1
  \param AABB2          axis - aligned bounding box 2
  \return merged axis-aligned bounding box
   */
  BlitzMat3x2 mergeAABB(
      const BlitzMat3x2&                      AABB1, 
      const BlitzMat3x2&                      AABB2  
      );

  
  /*!
  \brief tests if  AABB_old AND AABB_new are completely contained in node box
  \param AABB_old             AABB_old
  \param AABB_new             AABB_new
  \param nodeBox              AABB (e.g. node box)
   */
  bool inSameNodeBox(
      const BlitzMat3x2&  AABB_old, 
      const BlitzMat3x2&  AABB_new,
      const BlitzMat3x2&  nodeBox
      );
  
  /*!
  \brief  checks if point X is within AABB
  \param AABB          axis - aligned bounding box
  \param point         point coords
  \return true if point in AABB, false otherwise
  */
  bool isPointContainedInAABB(
      const BlitzMat3x2&  AABB,   
      const BlitzVec3&    point       
      );

  
  /*!
  \brief check the type of geometry of an element, rough because only linear or nonlinear is checked
  \param element              element
  \param xyze_element         nodal coordinates 
  \param eleGeoType           geometry type
   */ 
  void checkRoughGeoType(
      const DRT::Element*                     element,      
      const BlitzMat                          xyze_element, 
      XFEM::EleGeoType&                       eleGeoType  
      );


  
}

#endif
#endif


