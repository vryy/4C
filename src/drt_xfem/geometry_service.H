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
  //! gives the AABB of a complete Discretization and adds TOL07, doesn't care about displacements
  BlitzMat3x2         getXAABBofDis         (const DRT::Discretization& dis);  
  //! gives the AABB of a complete Discretization and adds TOL07, also takes into account current displacements
  BlitzMat3x2         getXAABBofDis         (const DRT::Discretization& dis,const std::map<int,BlitzVec3>& currentpositions);  
  //! gives back the nearest Element in an List of Elements
  const DRT::Element* nearestNeighbourInList(const DRT::Discretization& dis,const std::map<int,BlitzVec3>& currentpositions, const list<const DRT::Element* >& ElementList, const BlitzVec3& x, double& dist);
}

#endif
#endif


