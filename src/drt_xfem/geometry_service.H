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

#include "intersection_service.H"
#include <blitz/array.h>
#include "blitz_tiny_operation.H"
#include "xfem.H"
#include "interface.H"

namespace XFEM
{
  BlitzMat3x2 getXAABBofDis(const XFEM::InterfaceHandle& ih);
  const DRT::Element* nearestNeighbourInList(const XFEM::InterfaceHandle& ih, const list<const DRT::Element* > ElementList, const BlitzVec3& x, double& dist);
}
  
#endif
#endif


