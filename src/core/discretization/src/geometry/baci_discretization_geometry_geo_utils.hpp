/*----------------------------------------------------------------------*/
/*! \file

\brief Within this file all geometry typedefs shall be defined. Thus they can
be used from outside in header files without including the according geometry
header file
--> THIS FUNCTIONALITY IS JUST USED IN COMBUST AND WILL LEAVE BACI SOON

\level 3

 *------------------------------------------------------------------------------------------------*/


#ifndef BACI_DISCRETIZATION_GEOMETRY_GEO_UTILS_HPP
#define BACI_DISCRETIZATION_GEOMETRY_GEO_UTILS_HPP


#include "baci_config.hpp"

#include "baci_discretization_geometry_integrationcell.hpp"

#include <vector>

BACI_NAMESPACE_OPEN


// forward declarations
namespace CORE::GEO
{
  //! shortcut for a vector of BoundaryIntCells
  typedef std::vector<CORE::GEO::BoundaryIntCell> BoundaryIntCells;

  //! shortcut for a vector of BoundaryCell pointers
  typedef std::vector<Teuchos::RCP<CORE::GEO::BoundaryIntCell>> BoundaryIntCellPtrs;


  //! based on this element property, one can speed up geometry algorithms
  enum EleGeoType
  {
    CARTESIAN,
    LINEAR,
    HIGHERORDER
  };
}  // namespace CORE::GEO

BACI_NAMESPACE_CLOSE

#endif  // DISCRETIZATION_GEOMETRY_GEO_UTILS_H
