/*!-----------------------------------------------------------------------------------------------*
 \file levelset_intersection_utils.H

 \brief

 - collection of helper functions for level-set problems
 - allow for computing intersection of zero level-set iso-contour with discretization
   and related quantities, i.g., volume of subdomains, interface discretization, ...

\level 2

<pre>
\maintainer Magnus Winter
            winter@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
 *------------------------------------------------------------------------------------------------*/


#ifndef LEVELSET_INTERSECTION_UTILS_H
#define LEVELSET_INTERSECTION_UTILS_H

#include <Teuchos_RCP.hpp>
#include <Epetra_MultiVector.h>
#include <Epetra_MpiComm.h>

#include "../drt_geometry/geo_utils.H"

namespace DRT
{
  class Discretization;
  class PackBuffer;
}

namespace SCATRA
{
  /// construct zero iso-contour of level-set field
  void CaptureZeroLevelSet(
    const Teuchos::RCP<const Epetra_Vector>& phi,
    const Teuchos::RCP<const DRT::Discretization> & scatradis,
    double& volumedomainminus,
    double& volumedomainplus,
    double& zerosurface,
    std::map<int,GEO::BoundaryIntCells >& elementBoundaryIntCells);

  /// export boundary integration cells from this proc to parallel distribution
  void ExportInterface(
    std::map<int, GEO::BoundaryIntCells>& myinterface,
    const Epetra_Comm& comm);

  /// pack boundary integration cells from set into char array
  void packBoundaryIntCells(
      const std::map<int, GEO::BoundaryIntCells>& intcellmap,
      DRT::PackBuffer&                            dataSend
  );

  /// unpack boundary integration cells from char array
  void unpackBoundaryIntCells(
      const std::vector<char>&              dataRecv,
      std::map<int, GEO::BoundaryIntCells>& intcellmap
  );
}


#endif //LEVELSET_INTERSECTION_UTILS_H
