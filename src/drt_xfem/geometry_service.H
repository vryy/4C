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

namespace XFEM
{
  /*!
  \brief Gives the eXtendedAxisAlignedBoundingBox of a complete Discretization and adds TOL07, 
         also takes into account current displacements

  \param  dis               (in)    : discretization containing elements 
  \param  currentpositions  (in)    : current position of elements in dis
  \return                             AxisAlignedBoundingBox as BlitzMat3x2
  */
  BlitzMat3x2 getXAABBofDis(const DRT::Discretization& dis,const std::map<int,BlitzVec3>& currentpositions);  

  /*!
  \brief gives back the nearest Element in an List of Elements

  \param  dis               (in)    : discretization containing elements 
  \param  currentpositions  (in)    : current position of elements in dis
  \param  ElementList       (in)    : list of elements to be examined
  \param  x                 (in)    : coords of point which will be examined
  \param  dist              (out)   : distance btw. x and nearest element
  \return                             pointer to nearest element
  */
  const DRT::Element* nearestNeighbourInList(const DRT::Discretization& dis,const std::map<int,BlitzVec3>& currentpositions, const list<const DRT::Element* >& ElementList, const BlitzVec3& x, double& dist);

  /*!
  \brief merges 2 bounding boxes

  \param  A                 (in)    : bounding box A 
  \param  B                 (in)    : bounding box B
  \return                             merged bounding box
  */
  BlitzMat3x2 mergeAABB(const BlitzMat3x2& A, const BlitzMat3x2& B);
  
  //! gives overlap-area of 2 AABB
  double getOverlapArea(const BlitzMat3x2& A, const BlitzMat3x2& B, const BlitzMat3x2& C);
 
  //! gives overlap-area of 2 AABB
  double getOverlapArea(const list<BlitzMat3x2 > AABBs);
  
  //! gives of an AABB
  double getArea(const BlitzMat3x2& AABB);
  
  
}

#endif
#endif


