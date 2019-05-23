/*!------------------------------------------------------------------------------------------------*

\brief Within this file all geometry typedefs shall be defined. Thus they can
be used from outside in header files without including the according geometry
header file
--> THIS FUNCTIONALITY IS JUST USED IN COMBUST AND WILL LEAVE BACI SOON

\level 3

\maintainer Martin Kronbichler
 *------------------------------------------------------------------------------------------------*/


#ifndef GEO_UTILS_H_
#define GEO_UTILS_H_


#include <vector>
#include "../drt_geometry/integrationcell.H"


// forward declarations
namespace GEO
{
  //! shortcut for a vector of BoundaryIntCells
  typedef std::vector<GEO::BoundaryIntCell> BoundaryIntCells;

  //! shortcut for a vector of BoundaryCell pointers
  typedef std::vector<Teuchos::RCP<GEO::BoundaryIntCell>> BoundaryIntCellPtrs;

  //! shortcut for a vector of DomainIntCells
  typedef std::vector<GEO::DomainIntCell> DomainIntCells;


  //! based on this element property, one can speed up geometry algorithms
  enum EleGeoType
  {
    CARTESIAN,
    LINEAR,
    HIGHERORDER
  };
}  // namespace GEO

#endif /* GEO_UTILS_H_ */
