/*!
\file octtree_geometry_service.H

\brief provides geometry methods for oct tree

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
#include "../drt_geometry/vector_definitions.H"
#include "../drt_geometry/octtree_nearestobject.H"
#include "Teuchos_RefCountPtr.hpp"

namespace GEO
{

  /*!
  \brief Returns the eXtendedAxisAlignedBoundingBox of a discretization, 
         also takes into account current displacements
  \param dis                  discretization 
  \param currentpositions     current nodal positions in discretization
  \return  XAxisAlignedBoundingBox 
   */
  BlitzMat3x2 getXAABBofDis(
    const DRT::Discretization&            dis,               
    const std::map<int,BlitzVec3>&        currentpositions  
    );

  
  /*!
  \brief Returns the eXtendedAxisAlignedBoundingBox of a discretization, 
         in reference configuration
  \param dis                  discretization 
  \return  XAxisAlignedBoundingBox 
   */
  BlitzMat3x2 getXAABBofDis(
      const DRT::Discretization& dis
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
  \brief Returns the label for a point and a given map of elements
  \param dis                  discretization
  \param currentpositions     current nodal positions in discretization
  \param querypoint           point to be examined
  \param elementList          map of elements label gid
  \return                     label
   */
  int getXFEMLabel(
      const DRT::Discretization&              dis, 		
      const std::map<int,BlitzVec3>&          currentpositions, 	
      const BlitzVec3&                        querypoint,
      std::map<int, std::set<int> >&          elementList
      );  	


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
      const std::map<int,BlitzVec3>&              currentpositions,
      const std::map<int, std::set<int> >&        elementList,
      const BlitzVec3&                            point,
      BlitzVec3&                                  minDistanceVec,
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
      const std::map<int,BlitzVec3>&              currentpositions,
      const BlitzVec3&                            point,
      BlitzVec3&                                  x_surface_phys,
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
      const std::map<int,BlitzVec3>&              currentpositions,
      const BlitzVec3&                            point,
      BlitzVec3&                                  x_line_phys,
      double&                                     distance
      );

  
  /*!
  \brief Searches for the nearest node to a given point
  \param dis                  discretization
  \param currentpositions     current nodal positions in discretization
  \param point                point to be examined
  \param x_node               physical coordinates of tree node
  \param distance             distance to nearest point
   */
  void getDistanceToPoint(
      const DRT::Node*                            node,
      const std::map<int,BlitzVec3>&              currentpositions,
      const BlitzVec3&                            point,
      BlitzVec3&                                  x_node,
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
  BlitzVec3 getNormalAtSurfacePoint(
      const DRT::Discretization&                  dis,
      const std::map<int,BlitzVec3>&              currentpositions,
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
  \brief Merges two axis-aligned bounding boxes
  \param AABB1          axis - aligned bounding box 1
  \param AABB2          axis - aligned bounding box 2
  \return merged axis-aligned bounding box
   */
  BlitzMat3x2 mergeAABB(
      const BlitzMat3x2&                      AABB1, 
      const BlitzMat3x2&                      AABB2  
      );

  
  /*!
  \brief Tests if  AABB_old AND AABB_new are completely contained in node box
  \param AABB_old             AABB_old
  \param AABB_new             AABB_new
  \param nodeBox              AABB (e.g. node box)
  \return true if both AABB's are in node box
   */
  bool inSameNodeBox(
      const BlitzMat3x2&  AABB_old, 
      const BlitzMat3x2&  AABB_new,
      const BlitzMat3x2&  nodeBox
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
      const BlitzMat                          xyze_element, 
      XFEM::EleGeoType&                       eleGeoType  
      );


  
}

#endif
#endif


