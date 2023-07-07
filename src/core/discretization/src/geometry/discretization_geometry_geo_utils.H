/*----------------------------------------------------------------------*/
/*! \file

\brief Within this file all geometry typedefs shall be defined. Thus they can
be used from outside in header files without including the according geometry
header file
--> THIS FUNCTIONALITY IS JUST USED IN COMBUST AND WILL LEAVE BACI SOON

\level 3

 *------------------------------------------------------------------------------------------------*/


#ifndef DISCRETIZATION_GEOMETRY_GEO_UTILS_H
#define DISCRETIZATION_GEOMETRY_GEO_UTILS_H


#include <vector>
#include "discretization_geometry_integrationcell.H"


// forward declarations
namespace CORE::GEO
{
  //! shortcut for a vector of BoundaryIntCells
  typedef std::vector<CORE::GEO::BoundaryIntCell> BoundaryIntCells;

  //! shortcut for a vector of BoundaryCell pointers
  typedef std::vector<Teuchos::RCP<CORE::GEO::BoundaryIntCell>> BoundaryIntCellPtrs;

  //! shortcut for a vector of DomainIntCells
  typedef std::vector<CORE::GEO::DomainIntCell> DomainIntCells;


  //! based on this element property, one can speed up geometry algorithms
  enum EleGeoType
  {
    CARTESIAN,
    LINEAR,
    HIGHERORDER
  };
}  // namespace CORE::GEO

#endif  // DISCRETIZATION_GEOMETRY_GEO_UTILS_H
