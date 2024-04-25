/*----------------------------------------------------------------------*/
/*! \file
\brief Definitions of cell types
\level 1
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_DISCRETIZATION_FEM_GENERAL_CELL_TYPE_HPP
#define FOUR_C_DISCRETIZATION_FEM_GENERAL_CELL_TYPE_HPP

#include "4C_config.hpp"

FOUR_C_NAMESPACE_OPEN

namespace CORE::FE
{
  enum class CellType
  {
    dis_none,    ///< unknown dis type
    quad4,       ///< 4 noded quadrilateral
    quad6,       ///< 6 noded quadrilateral (linear in one direction, quadratic in the other)
    quad8,       ///< 8 noded quadrilateral
    quad9,       ///< 9 noded quadrilateral
    tri3,        ///< 3 noded triangle
    tri6,        ///< 6 noded triangle
    hex8,        ///< 8 noded hexahedra
    hex16,       ///< 16 noded hexahedra
    hex18,       ///< 18 noded hexahedra
    hex20,       ///< 20 noded hexahedra
    hex27,       ///< 27 noded hexahedra
    tet4,        ///< 4 noded tetrahedra
    tet10,       ///< 10 noded tetrahedra
    wedge6,      ///< 6 noded wedge
    wedge15,     ///< 15 noded wedge
    pyramid5,    ///< 5 noded pyramid
    line2,       ///< 2 noded line
    line3,       ///< 3 noded line
    line4,       ///< 4 noded line
    line5,       ///< 5 noded line
    line6,       ///< 6 noded line
    point1,      ///< 1 noded point
    nurbs2,      ///< 2 control point first order nurbs line element
    nurbs3,      ///< 3 control point second order nurbs line element
    nurbs4,      ///< 4 control point first order nurbs surface element
    nurbs9,      ///< 9 control point second order nurbs surface element
    nurbs8,      ///< 8 control point first order nurbs volume element
    nurbs27,     ///< 27 control point second order nurbs volume element
    max_distype  ///< end marker. must be the last entry
  };

  /*!
   * @brief A compile time sequence of celltypes.
   *
   * @tparam celltypes Celltypes in the sequence
   */
  template <CellType... celltypes>
  struct CelltypeSequence
  {
  };
}  // namespace CORE::FE

FOUR_C_NAMESPACE_CLOSE

#endif