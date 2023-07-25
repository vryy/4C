/*----------------------------------------------------------------------*/
/*! \file

\brief General input parameters


\level 1
*/
/*----------------------------------------------------------------------*/
#ifndef BACI_INPAR_H
#define BACI_INPAR_H

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
#endif  // INPAR_H
