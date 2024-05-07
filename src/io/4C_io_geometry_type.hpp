/*----------------------------------------------------------------------*/
/*! \file
\brief Parameters for geometry input
\level 1
*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_IO_GEOMETRY_TYPE_HPP
#define FOUR_C_IO_GEOMETRY_TYPE_HPP

#include "4C_config.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
namespace IO
{
  //! Geometry reading specification
  enum GeometryType
  {
    geometry_full,
    geometry_box,
    geometry_file
  };

}  // namespace IO

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
