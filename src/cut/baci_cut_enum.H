/*----------------------------------------------------------------------*/
/*! \file

\brief All enumerators used in cut libraries

\level 3

*/

/*----------------------------------------------------------------------*/

#ifndef BACI_CUT_ENUM_H
#define BACI_CUT_ENUM_H

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
namespace CORE::GEO
{
  namespace CUT
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

  }  // namespace CUT

}  // namespace CORE::GEO

/*----------------------------------------------------------------------*/
BACI_NAMESPACE_CLOSE

#endif
