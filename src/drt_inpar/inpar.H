/*----------------------------------------------------------------------*/
/*!

\brief General input parameters

\maintainer Martin Kronbichler

\level 1
*/
/*----------------------------------------------------------------------*/
#ifndef INPAR_H
#define INPAR_H

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
#endif  // #ifndef INPAR_H
