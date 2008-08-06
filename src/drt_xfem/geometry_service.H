/*!
\file geometry_service.H

\brief provides service functions for algorithmic geomerty

<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
 */

#ifdef CCADISCRET
#ifndef GEOMETRY_SERVICE_H_
#define GEOMETRY_SERVICE_H_

#include "blitz_tiny_operation.H"
#include "../drt_geometry/vector_definitions.H"
#include "../drt_lib/drt_discret.H"
#include "intersection_service.H"

namespace XFEM
{

  //! this proporty decides which distance algorithm is used
  enum DistanceType {ELEMENT_SURFACE, ELEMENT_LINE, ELEMENT_POINT};

  /*!
  \brief Gives the eXtendedAxisAlignedBoundingBox of a complete Discretization and adds TOL07, 
         also takes into account current displacements
  \return  AxisAlignedBoundingBox as BlitzMat3x2
  */
  BlitzMat3x2 getXAABBofDis(
      const DRT::Discretization&            dis, ///<  discretization containing elements
      const std::map<int,BlitzVec3>&        currentpositions ///<current position of elements in dis
      );  

  /*!
  \brief Gives the eXtendedAxisAlignedBoundingBox of a complete Discretization and adds TOL07, 
         also takes into account current displacements
  \return                             AxisAlignedBoundingBox as BlitzMat3x2
  */
  BlitzMat3x2 getXAABBofDis(
      const DRT::Discretization&              dis ///<discretization containing elements
      );  

  /*!
  \brief gives back the nearest Element in an List of Elements (wrapper)
  \return                             pointer to nearest element
  */
  const DRT::Element* nearestNeighbourInList(
      const DRT::Discretization&              dis, ///< discretization containing elements 
      const std::map<int,BlitzVec3>&          currentpositions, ///< current position of elements in dis
      const list<const DRT::Element* >&       ElementList, ///<list of elements to be examined
      const BlitzVec3&                        X,  ///<coords of point which will be examined
      map<const DRT::Element*, double >&      squaredDistanceMap, ///<> map of already known distances (may be empty)
      BlitzVec3&                              vector2minDistPoint,  ///>vector to point on surface element with minimal distance
      double&                                 eleDistance, ///>distance btw. x and nearest element
      DistanceType&                           distanceType ///>distance Typ                            
      );

  /*!
   \brief  gives the Element from a list whose XAABB has the minimum maxDistance to the Point
   \return                             pointer to nearest element
   */
  const DRT::Element* nearestXAABBNeighbourInList(
      const DRT::Discretization&              dis,///< discretization containing elements
      const std::map<int,BlitzVec3>&          currentpositions,///< current position of elements in dis 
      const list<const DRT::Element* >&       ElementList, ///<list of elements to be examined
      const BlitzVec3&                        X, ///<coords of point which will be examined
      double&                                 dist ///>distance btw. x and nearest element
      );

  /*!
  \brief merges 2 bounding boxes
  \return                             merged bounding box
  */
  BlitzMat3x2 mergeAABB(
      const BlitzMat3x2&                      A, ///<bounding box A 
      const BlitzMat3x2&                      B  ///<bounding box B
      );

  /*!
  \brief gives overlap-percentage of 2 AABBs
  \return true if A is contained in B, else false
  */
  bool isContainedAinB(
      const BlitzMat3x2&                      A, ///<bounding box A 
      const BlitzMat3x2&                      B  ///<bounding box B
      );
    
  /*!
  \brief gives overlap-area of 3 AABBs
  \return Volume of the overlap
  */
  double getOverlapArea(
      const BlitzMat3x2&                      A, ///<bounding box A 
      const BlitzMat3x2&                      B, ///<bounding box B 
      const BlitzMat3x2&                      C  ///<bounding box C 
      );
 
  /*!
  \brief gives overlap-area of n AABBs#
  \return Volume of the overlap
  */
  double getOverlapArea(
      const list<BlitzMat3x2 >                AABBs ///< list of bounding boxes
      );
  
  /*!
  \brief  gives of an AABB
  \return Volume of the bounding box 
  */
  double getVolume(
      const BlitzMat3x2&                      AABB ///<bounding box
      );
  
  /*!
  \brief  gives maximum distance of a point from an AABB (distance from farest point)
  \return distance (unsigned)
  */
  double getMaxDistanceFromAABB(
      const BlitzVec3&                        X, ///< point coords
      const BlitzMat3x2                       AABB ///<bounding box
      );

  /*!
  \brief  gives AABB of an circle
  \return bounding box
  */
  BlitzMat3x2 getAABBofSphere(
      const BlitzVec3&                        X, ///< point coords
      const double                            radius ///<radius
      );

  /*!
  \brief check the type of geometry of an element, rough because only linear or nonlinear is checked
  */ 
  void checkRoughGeoType(
      const DRT::Element*                     element, ///< element
      const BlitzMat                          xyze_element,///< actual point positions in element
      XFEM::EleGeoType&                       eleGeoType ///> geometry type
      );
  
  
  /*!
  \brief  gives back the squared (euclidian) distance of a point to an surface element 
          1.) get nearest point on element support in xsi-coordinates
          4.) if is truly(< not <=) in xsi-limits, nearest point is on the element surface
          2.) if ((xsi0<xsi0_min||xsi0>xsi0_max) XOR (xsi1<xsi1_min || xsi1>xsi1_max)) 
              nearest point is on element surface (which is a line because the element itself is 2D) 
          3.) if xsi is totally out of bounds the nearest point on this element is one of the node points
  \return squared distance
  */
  double getSquaredElementDistance(
      const DRT::Element*                     surfaceElement, //< surface element to examine
      const std::map<int,BlitzVec3>&          currentpositions,///< actual point positions in all elements
      const BlitzVec3&                        physCoord, ///< coordinates of point (physical)
      BlitzVec2&                              xsi, ///> coordinates of nearest point (parameter space)
      BlitzVec3&                              x_surface_phys,///> vector 2 minimal distance point on surface
      DistanceType&                           distanceTyp///> type of distance (nearest point in element, nearest point on boundary: line or node)
      );

  /*!
  \brief  gives back the squared (euclidian) distance of a point to an surface element in case
          the nearest point is on the surface
  \return squared distance
  */
  double getSquaredElementDistance_Surface(
      const DRT::Element*                     surfaceElement,//< surface element to examine
      const std::map<int,BlitzVec3>&          currentpositions,///< actual point positions in all element
      const BlitzVec3&                        physCoord,///<  coordinates of point (physical)
      BlitzVec2&                              xsi, ///<> coordinates of nearest point (parameter space)
      BlitzVec3&                              x_surface_phys///> vector 2 minimal distance point on surface
      );

  /*!
  \brief  gives back the squared (euclidian) distance of a point to an surface element in case
          the nearest point is on a boundary line
  \return squared distance
  */
  double getSquaredElementDistance_Line(
      const DRT::Element*                     surfaceElement,//< surface element to examine
      const std::map<int,BlitzVec3>&          currentpositions,///< actual point positions in all elements
      const BlitzVec3&                        physCoord,///<  coordinates of point (physical)
      BlitzVec2&                              xsi,///<> coordinates of nearest point (parameter space)
      BlitzVec3&                              x_surface_phys);///> vector 2 minimal distance point on surface

  
  /*!
  \brief  gives back the squared (euclidian) distance of a point to an surface element in case
          the nearest point is a node
  \return squared distance
  */
  double getSquaredElementDistance_Point(    
      const DRT::Element*                     surfaceElement,//< surface element to examine
      const std::map<int,BlitzVec3>&          currentpositions,///< actual point positions in all elements
      const BlitzVec3&                        physCoord,///<  coordinates of point (physical)
      BlitzVec2&                              xsi,///<> coordinates of nearest point (parameter space)
      BlitzVec3&                              x_surface_phys);///<> coordinates of nearest point (physical space)
  
  
  /*!
  \brief  gives back the normal to an surface element at position xsi 
  \return normal vector in physical coordinates
  */
  BlitzVec3 getNormalAtXsi(
      const DRT::Discretization&              dis, ///<  discretization containing elements
      const DRT::Element*                     surfaceElement,//< surface element to examine
      const std::map<int,BlitzVec3>&          currentpositions, ///< actual point positions in all elements
      const BlitzVec2&                        xsi, //< coordinates of point to get normal on (parameter space)
      const BlitzVec3&                        X, //< coordinates of point we at least want to know the distance from (physical)
      const DistanceType&                     DistanceType///< type of distance (nearest point in element, nearest point on boundary: line or node)
      );

  
  /*!
  \brief  find all elements which A and B have in common
  \return list of common elements
  */
  list < const DRT::Element* > getCommonElements(
      const DRT::Node*                        A, ///< node A
      const DRT::Node*                        B  ///< node B
      );
  
  /*!
  \brief  calculate vector sum of a list of vectors
  \return vector sum
  */
  BlitzVec3 addVectors(
      const list< BlitzVec3 >                 vectors ///< list of vectors
      );
  
  /*!
  \brief  get the minimal and maximal xsi-values of a 2D Element in parameter-space
  \return bouding box for xsi values
  */
  BlitzMat2x2 getXsiBoundingBox(
      const DRT::Element*                     surfaceElement ///< element
      );
  
  
  /*!
  \brief  get Node of the Element which has coordinates xsi, within TOL7
          ATTENTION: throws dserror() if Node is not found, because this indicates a bug
  \return Node
  */
  const DRT::Node* getNodeAtXsi(
      const DRT::Element*                     surfaceElement, ///< element to search on
      const std::map<int,BlitzVec3>&          currentpositions, ///< actual point positions in all elements
      const BlitzVec2&                        xsi ///< xsi value of searched Node
      );
  
  
  /*!
  \brief  get the maximum radius of a circle segment(which is bounded by 
          AABB and centered at C) pointing away from X
  \return Node
  */
  double biggestRadiusInAABB(
      const BlitzMat3x2&                      AABB, ///<bounding box
      const BlitzVec3&                        X ///<point coords                  
      );
  
  /*!
  \brief  checks if point X is within AABB
  \return true/false
  */
  bool isContainedXinAABB(
      const BlitzMat3x2& AABB, ///<bounding box
      const BlitzVec3& X //<point coords    
      );

  
  
}

#endif
#endif


