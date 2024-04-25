/*----------------------------------------------------------------------*/
/*! \file

\brief General input parameters


\level 1
*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_INPAR_HPP
#define FOUR_C_INPAR_HPP

#include "4C_config.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
namespace INPAR
{
  //! Geometry reading specification
  enum GeometryType
  {
    geometry_full,
    geometry_box,
    geometry_file
  };

}  // namespace INPAR

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
