/*----------------------------------------------------------------------*/
/*! \file

\brief All enumerators used in cut libraries

\level 3

\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
*/

/*----------------------------------------------------------------------*/

#ifndef CUT_ENUM_H
#define CUT_ENUM_H

/*----------------------------------------------------------------------*/
namespace GEO
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

}  // namespace GEO

/*----------------------------------------------------------------------*/
#endif
