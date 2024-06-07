/*----------------------------------------------------------------------*/
/*! \file

\brief Within this file all geometry typedefs shall be defined. Thus they can
be used from outside in header files without including the according geometry
header file
--> THIS FUNCTIONALITY IS JUST USED IN COMBUST AND WILL LEAVE 4C SOON

\level 3

 *------------------------------------------------------------------------------------------------*/


#ifndef FOUR_C_FEM_GEOMETRY_GEO_UTILS_HPP
#define FOUR_C_FEM_GEOMETRY_GEO_UTILS_HPP


#include "4C_config.hpp"

#include "4C_fem_geometry_integrationcell.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN


// forward declarations
namespace Core::Geo
{
  //! shortcut for a vector of BoundaryIntCells
  typedef std::vector<Core::Geo::BoundaryIntCell> BoundaryIntCells;

  //! shortcut for a vector of BoundaryCell pointers
  typedef std::vector<Teuchos::RCP<Core::Geo::BoundaryIntCell>> BoundaryIntCellPtrs;


  //! based on this element property, one can speed up geometry algorithms
  enum EleGeoType
  {
    CARTESIAN,
    LINEAR,
    HIGHERORDER
  };
}  // namespace Core::Geo

FOUR_C_NAMESPACE_CLOSE

#endif
