/*----------------------------------------------------------------------*/
/*! \file
\brief Definitions of cell types traits.

In this file, the helper defining the mapping from celltype to their properties are defined

\level 1
*/
/*----------------------------------------------------------------------*/

#ifndef BACI_DISCRETIZATION_FEM_GENERAL_CELL_TYPE_TRAITS_DEF_H
#define BACI_DISCRETIZATION_FEM_GENERAL_CELL_TYPE_TRAITS_DEF_H

#include "baci_discretization_fem_general_cell_type.H"

#include <tuple>
#include <type_traits>

namespace CORE::FE::DETAILS
{
  //! @name Definition of the mapping celltype to string
  /// @{
  template <CellType celltype>
  struct CellTypeToString
  {
  };

  template <>
  struct CellTypeToString<CORE::FE::CellType::dis_none>
  {
    static constexpr auto value = "DIS_NONE";
  };

  template <>
  struct CellTypeToString<CORE::FE::CellType::quad4>
  {
    static constexpr auto value = "QUAD4";
  };

  template <>
  struct CellTypeToString<CORE::FE::CellType::quad6>
  {
    static constexpr auto value = "QUAD6";
  };

  template <>
  struct CellTypeToString<CORE::FE::CellType::quad8>
  {
    static constexpr auto value = "QUAD8";
  };

  template <>
  struct CellTypeToString<CORE::FE::CellType::quad9>
  {
    static constexpr auto value = "QUAD9";
  };

  template <>
  struct CellTypeToString<CORE::FE::CellType::tri3>
  {
    static constexpr auto value = "TRI3";
  };

  template <>
  struct CellTypeToString<CORE::FE::CellType::tri6>
  {
    static constexpr auto value = "TRI6";
  };

  template <>
  struct CellTypeToString<CORE::FE::CellType::hex8>
  {
    static constexpr auto value = "HEX8";
  };

  template <>
  struct CellTypeToString<CORE::FE::CellType::hex16>
  {
    static constexpr auto value = "HEX16";
  };

  template <>
  struct CellTypeToString<CORE::FE::CellType::hex18>
  {
    static constexpr auto value = "HEX18";
  };

  template <>
  struct CellTypeToString<CORE::FE::CellType::hex20>
  {
    static constexpr auto value = "HEX20";
  };

  template <>
  struct CellTypeToString<CORE::FE::CellType::hex27>
  {
    static constexpr auto value = "HEX27";
  };

  template <>
  struct CellTypeToString<CORE::FE::CellType::tet4>
  {
    static constexpr auto value = "TET4";
  };

  template <>
  struct CellTypeToString<CORE::FE::CellType::tet10>
  {
    static constexpr auto value = "TET10";
  };

  template <>
  struct CellTypeToString<CORE::FE::CellType::wedge6>
  {
    static constexpr auto value = "WEDGE6";
  };

  template <>
  struct CellTypeToString<CORE::FE::CellType::wedge15>
  {
    static constexpr auto value = "WEDGE15";
  };

  template <>
  struct CellTypeToString<CORE::FE::CellType::pyramid5>
  {
    static constexpr auto value = "PYRAMID5";
  };

  template <>
  struct CellTypeToString<CORE::FE::CellType::line2>
  {
    static constexpr auto value = "LINE2";
  };

  template <>
  struct CellTypeToString<CORE::FE::CellType::line3>
  {
    static constexpr auto value = "LINE3";
  };

  template <>
  struct CellTypeToString<CORE::FE::CellType::line4>
  {
    static constexpr auto value = "LINE4";
  };

  template <>
  struct CellTypeToString<CORE::FE::CellType::line5>
  {
    static constexpr auto value = "LINE5";
  };

  template <>
  struct CellTypeToString<CORE::FE::CellType::line6>
  {
    static constexpr auto value = "LINE6";
  };

  template <>
  struct CellTypeToString<CORE::FE::CellType::point1>
  {
    static constexpr auto value = "POINT1";
  };

  template <>
  struct CellTypeToString<CORE::FE::CellType::nurbs2>
  {
    static constexpr auto value = "NURBS2";
  };

  template <>
  struct CellTypeToString<CORE::FE::CellType::nurbs3>
  {
    static constexpr auto value = "NURBS3";
  };

  template <>
  struct CellTypeToString<CORE::FE::CellType::nurbs4>
  {
    static constexpr auto value = "NURBS4";
  };

  template <>
  struct CellTypeToString<CORE::FE::CellType::nurbs9>
  {
    static constexpr auto value = "NURBS9";
  };

  template <>
  struct CellTypeToString<CORE::FE::CellType::nurbs8>
  {
    static constexpr auto value = "NURBS8";
  };

  template <>
  struct CellTypeToString<CORE::FE::CellType::nurbs27>
  {
    static constexpr auto value = "NURBS27";
  };

  template <>
  struct CellTypeToString<CORE::FE::CellType::max_distype>
  {
    static constexpr auto value = "MAX_DISTYPE";
  };
  /// @}

  //! @name Helper functions to convert cell types to the shape of the cell
  /// @{
  template <CORE::FE::CellType celltype>
  using is_tet = std::bool_constant<celltype == CORE::FE::CellType::tet4 ||
                                    celltype == CORE::FE::CellType::tet10>;

  template <CORE::FE::CellType celltype>
  using is_hex = std::bool_constant<
      celltype == CORE::FE::CellType::hex8 || celltype == CORE::FE::CellType::hex16 ||
      celltype == CORE::FE::CellType::hex18 || celltype == CORE::FE::CellType::hex20 ||
      celltype == CORE::FE::CellType::hex27>;

  template <CORE::FE::CellType celltype>
  using is_wedge = std::bool_constant<celltype == CORE::FE::CellType::wedge6 ||
                                      celltype == CORE::FE::CellType::wedge15>;

  template <CORE::FE::CellType celltype>
  using is_pyramid = std::bool_constant<celltype == CORE::FE::CellType::pyramid5>;

  template <CORE::FE::CellType celltype>
  using is_quad = std::bool_constant<
      celltype == CORE::FE::CellType::quad4 || celltype == CORE::FE::CellType::quad6 ||
      celltype == CORE::FE::CellType::quad8 || celltype == CORE::FE::CellType::quad9>;

  template <CORE::FE::CellType celltype>
  using is_tri = std::bool_constant<celltype == CORE::FE::CellType::tri3 ||
                                    celltype == CORE::FE::CellType::tri6>;

  template <CORE::FE::CellType celltype>
  using is_line = std::bool_constant<
      celltype == CORE::FE::CellType::line2 || celltype == CORE::FE::CellType::line3 ||
      celltype == CORE::FE::CellType::line4 || celltype == CORE::FE::CellType::line5 ||
      celltype == CORE::FE::CellType::line6>;
  /// @}
}  // namespace CORE::FE::DETAILS

#endif