/*----------------------------------------------------------------------*/
/*! \file

\brief All enumerators used in cut libraries

\level 3

*/

/*----------------------------------------------------------------------*/

#ifndef FOUR_C_CUT_ENUM_HPP
#define FOUR_C_CUT_ENUM_HPP

#include "4C_config.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
namespace Core::Geo
{
  namespace Cut
  {
    enum FacetShape
    {
      Convex,           // convex facet
      SinglePtConcave,  // only one concave point --> special triangulation can be used
      Concave           // concave facet
    };

    enum ProjectionDirection
    {
      proj_x,
      proj_y,
      proj_z
    };

  }  // namespace Cut

}  // namespace Core::Geo

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
