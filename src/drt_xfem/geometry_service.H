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

//#include "intersection_service.H"
#include <blitz/array.h>
#include "blitz_tiny_operation.H"
#include "xfem.H"
#include "../drt_lib/drt_discret.H"
#include "intersection_service.H"

namespace XFEM
{
  const static int DISTANCE_TO_ELEMENT_SURFACE=1;
  const static int DISTANCE_TO_ELEMENT_LINE=2;
  const static int DISTANCE_TO_ELEMENT_POINT=3;

  struct compareByEleGid
  {
    inline bool operator()(const DRT::Element* A, const DRT::Element* B) const
    {
      return (A->Id())<(B->Id()) ? true : false;
    }
  };
  
  /*!
   * An set with element pointers, which is ordered by global Element id
   */
  typedef std::set< const DRT::Element*, compareByEleGid > EleSet;
  
  /*!
  \brief Gives the eXtendedAxisAlignedBoundingBox of a complete Discretization and adds TOL07, 
         also takes into account current displacements

  \param  dis               (in)    : discretization containing elements 
  \param  currentpositions  (in)    : current position of elements in dis
  \return                             AxisAlignedBoundingBox as BlitzMat3x2
  */
  BlitzMat3x2 getXAABBofDis(const DRT::Discretization& dis,const std::map<int,BlitzVec3>& currentpositions);  

  /*!
  \brief Gives the eXtendedAxisAlignedBoundingBox of a complete Discretization and adds TOL07, 
         also takes into account current displacements

  \param  dis               (in)    : discretization containing elements 
  \return                             AxisAlignedBoundingBox as BlitzMat3x2
  */
  BlitzMat3x2 getXAABBofDis(const DRT::Discretization& dis);  
  
  
  /*!
  \brief gives back the nearest Element in an List of Elements

  \param  dis               (in)    : discretization containing elements 
  \param  currentpositions  (in)    : current position of elements in dis
  \param  ElementList       (in)    : list of elements to be examined
  \param  x                 (in)    : coords of point which will be examined
  \param  dist              (out)   : distance btw. x and nearest element
  \return                             pointer to nearest element
  */
  const DRT::Element* nearestNeighbourInList(
      const DRT::Discretization& dis,
      const std::map<int,BlitzVec3>& currentpositions, 
      const list<const DRT::Element* >& ElementList, 
      const BlitzVec3& x, 
      double& dist);
  
  const DRT::Element* nearestNeighbourInListOld(
      const DRT::Discretization& dis,
      const std::map<int,BlitzVec3>& currentpositions, 
      const list<const DRT::Element* >& ElementList, 
      const BlitzVec3& x, 
      double& dist);
  
  const DRT::Element* nearestNeighbourInListNew(
      const DRT::Discretization& dis,
      const std::map<int,BlitzVec3>& currentpositions, 
      const list<const DRT::Element* >& ElementList, 
      const BlitzVec3& x, 
      double& distance);
  
  //! gives the Element of which the XAABB has the minimum maxDistance to the Point
  const DRT::Element* nearestXAABBNeighbourInList(
      const DRT::Discretization& dis,
      const std::map<int,BlitzVec3>& currentpositions, 
      const list<const DRT::Element* >& ElementList, 
      const BlitzVec3& x, 
      double& dist);

  /*!
  \brief merges 2 bounding boxes

  \param  A                 (in)    : bounding box A 
  \param  B                 (in)    : bounding box B
  \return                             merged bounding box
  */
  BlitzMat3x2 mergeAABB(const BlitzMat3x2& A, const BlitzMat3x2& B);

  //! gives overlap-percentage of 2 AABBs
  bool isContainedAinB(const BlitzMat3x2& A, const BlitzMat3x2& B);
    
  //! gives overlap-area of 3 AABBs
  double getOverlapArea(const BlitzMat3x2& A, const BlitzMat3x2& B, const BlitzMat3x2& C);
 
  //! gives overlap-area of n AABBs
  double getOverlapArea(const list<BlitzMat3x2 > AABBs);
  
  //! gives of an AABB
  double getVolume(const BlitzMat3x2& AABB);
  
  //! gives maximum distance of a point from an AABB
  double getMaxDistanceFromAABB(const BlitzVec3& X, const BlitzMat3x2 AABB);

  //! gives AABB of an circle
  BlitzMat3x2 getAABBofSphere(const BlitzVec3& X, const double radius);

  
  void checkRoughGeoType(
           DRT::Element*                element,
           const BlitzMat               xyze_element,
           EleGeoType&                  eleGeoType
           );
  
  double getSquaredElementDistance(
      const DRT::Element*                     surfaceElement,
      const std::map<int,BlitzVec3>&                currentpositions,
      const BlitzVec3&                        physCoord,
      BlitzVec2&                              xsi,
      BlitzVec3&                              normal,
      int&                                    distanceTyp);
  
  double getSquaredElementDistance_Surface(
      const DRT::Element*                     surfaceElement,
      const std::map<int,BlitzVec3>&                currentpositions,
      const BlitzVec3&                        physCoord,
      BlitzVec2&                              xsi,
      BlitzVec3&                              vector2minDistPoint);

  double getSquaredElementDistance_Line(
      const DRT::Element*                     surfaceElement,
      const std::map<int,BlitzVec3>&                currentpositions,
      const BlitzVec3&                        physCoord,
      BlitzVec2&                              xsi,
      BlitzVec3&                              vector2minDistPoint);

  double getSquaredElementDistance_Point(    
      const DRT::Element*                     surfaceElement,
      const std::map<int,BlitzVec3>&                currentpositions,
      const BlitzVec3&                        physCoord,
      BlitzVec2&                              xsi,
      BlitzVec3&                              vector2minDistPoint);
  
  BlitzVec3 getNormalAtXsi(
      const DRT::Discretization& dis,
      const DRT::Element* surfaceElement,
      const std::map<int,BlitzVec3>& currentpositions, 
      const BlitzVec2& xsi, 
      const BlitzVec3& X,
      const int& ElementType);
  
  list < const DRT::Element* > getCommonElements(const DRT::Node* A,const DRT::Node* B);
  
  BlitzVec3 addVectors(const list< BlitzVec3 > vectors);
  BlitzMat2x2 getXsiBoundingBox(const DRT::Element* surfaceElement);
  const DRT::Node* getNodeAtXsi(
      const DRT::Element*       surfaceElement, 
      const std::map<int,BlitzVec3>&  currentpositions,
      const BlitzVec2&          xsi);
  
}

#endif
#endif


