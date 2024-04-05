/*----------------------------------------------------------------------*/
/*! \file

\brief General input parameters


\level 1
*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_INPAR_HPP
#define FOUR_C_INPAR_HPP

#include "baci_config.hpp"

BACI_NAMESPACE_OPEN

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
BACI_NAMESPACE_CLOSE

#endif
