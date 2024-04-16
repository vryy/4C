/*----------------------------------------------------------------------*/
/*! \file
\brief Definitions of cell types traits.

In this file, the helper defining the mapping from celltype to their properties are defined

\level 1
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_DISCRETIZATION_FEM_GENERAL_CELL_TYPE_TRAITS_DEF_HPP
#define FOUR_C_DISCRETIZATION_FEM_GENERAL_CELL_TYPE_TRAITS_DEF_HPP

#include "baci_config.hpp"

#include "baci_discretization_fem_general_cell_type.hpp"

#include <tuple>
#include <type_traits>

BACI_NAMESPACE_OPEN

namespace CORE::FE::DETAILS
{
  //! @name Definition of the mapping celltype to string
  /// @{
  template <CellType celltype>
  struct CellTypeInformation
  {
  };

  template <>
  struct CellTypeInformation<CORE::FE::CellType::dis_none>
  {
    static constexpr auto name = "DIS_NONE";
  };

  template <>
  struct CellTypeInformation<CORE::FE::CellType::quad4>
  {
    static constexpr auto name = "QUAD4";
    static constexpr int dim = 2;
    static constexpr int num_nodes = 4;
    static constexpr int num_faces = 4;
  };

  template <>
  struct CellTypeInformation<CORE::FE::CellType::quad6>
  {
    static constexpr auto name = "QUAD6";
    static constexpr int dim = 2;
    static constexpr int num_nodes = 6;
    static constexpr int num_faces = 4;
  };

  template <>
  struct CellTypeInformation<CORE::FE::CellType::quad8>
  {
    static constexpr auto name = "QUAD8";
    static constexpr int dim = 2;
    static constexpr int num_nodes = 8;
    static constexpr int num_faces = 4;
  };

  template <>
  struct CellTypeInformation<CORE::FE::CellType::quad9>
  {
    static constexpr auto name = "QUAD9";
    static constexpr int dim = 2;
    static constexpr int num_nodes = 9;
    static constexpr int num_faces = 4;
  };

  template <>
  struct CellTypeInformation<CORE::FE::CellType::tri3>
  {
    static constexpr auto name = "TRI3";
    static constexpr int dim = 2;
    static constexpr int num_nodes = 3;
    static constexpr int num_faces = 3;
  };

  template <>
  struct CellTypeInformation<CORE::FE::CellType::tri6>
  {
    static constexpr auto name = "TRI6";
    static constexpr int dim = 2;
    static constexpr int num_nodes = 6;
    static constexpr int num_faces = 3;
  };

  template <>
  struct CellTypeInformation<CORE::FE::CellType::hex8>
  {
    static constexpr auto name = "HEX8";
    static constexpr int dim = 3;
    static constexpr int num_nodes = 8;
    static constexpr int num_faces = 6;
  };

  template <>
  struct CellTypeInformation<CORE::FE::CellType::hex16>
  {
    static constexpr auto name = "HEX16";
    static constexpr int dim = 3;
    static constexpr int num_nodes = 16;
    static constexpr int num_faces = 6;
  };

  template <>
  struct CellTypeInformation<CORE::FE::CellType::hex18>
  {
    static constexpr auto name = "HEX18";
    static constexpr int dim = 3;
    static constexpr int num_nodes = 18;
    static constexpr int num_faces = 6;
  };

  template <>
  struct CellTypeInformation<CORE::FE::CellType::hex20>
  {
    static constexpr auto name = "HEX20";
    static constexpr int dim = 3;
    static constexpr int num_nodes = 20;
    static constexpr int num_faces = 6;
  };

  template <>
  struct CellTypeInformation<CORE::FE::CellType::hex27>
  {
    static constexpr auto name = "HEX27";
    static constexpr int dim = 3;
    static constexpr int num_nodes = 27;
    static constexpr int num_faces = 6;
  };

  template <>
  struct CellTypeInformation<CORE::FE::CellType::tet4>
  {
    static constexpr auto name = "TET4";
    static constexpr int dim = 3;
    static constexpr int num_nodes = 4;
    static constexpr int num_faces = 4;
  };

  template <>
  struct CellTypeInformation<CORE::FE::CellType::tet10>
  {
    static constexpr auto name = "TET10";
    static constexpr int dim = 3;
    static constexpr int num_nodes = 10;
    static constexpr int num_faces = 4;
  };

  template <>
  struct CellTypeInformation<CORE::FE::CellType::wedge6>
  {
    static constexpr auto name = "WEDGE6";
    static constexpr int dim = 3;
    static constexpr int num_nodes = 6;
    static constexpr int num_faces = 5;
  };

  template <>
  struct CellTypeInformation<CORE::FE::CellType::wedge15>
  {
    static constexpr auto name = "WEDGE15";
    static constexpr int dim = 3;
    static constexpr int num_nodes = 15;
    static constexpr int num_faces = 5;
  };

  template <>
  struct CellTypeInformation<CORE::FE::CellType::pyramid5>
  {
    static constexpr auto name = "PYRAMID5";
    static constexpr int dim = 3;
    static constexpr int num_nodes = 5;
    static constexpr int num_faces = 5;
  };

  template <>
  struct CellTypeInformation<CORE::FE::CellType::line2>
  {
    static constexpr auto name = "LINE2";
    static constexpr int dim = 1;
    static constexpr int num_nodes = 2;
    static constexpr int num_faces = 2;
    static constexpr int order = 1;
  };

  template <>
  struct CellTypeInformation<CORE::FE::CellType::line3>
  {
    static constexpr auto name = "LINE3";
    static constexpr int dim = 1;
    static constexpr int num_nodes = 3;
    static constexpr int num_faces = 2;
    static constexpr int order = 2;
  };

  template <>
  struct CellTypeInformation<CORE::FE::CellType::line4>
  {
    static constexpr auto name = "LINE4";
    static constexpr int dim = 1;
    static constexpr int num_nodes = 4;
    static constexpr int num_faces = 2;
    static constexpr int order = 3;
  };

  template <>
  struct CellTypeInformation<CORE::FE::CellType::line5>
  {
    static constexpr auto name = "LINE5";
    static constexpr int dim = 1;
    static constexpr int num_nodes = 5;
    static constexpr int num_faces = 2;
    static constexpr int order = 4;
  };

  template <>
  struct CellTypeInformation<CORE::FE::CellType::line6>
  {
    static constexpr auto name = "LINE6";
    static constexpr int dim = 1;
    static constexpr int num_nodes = 6;
    static constexpr int num_faces = 2;
    static constexpr int order = 5;
  };

  template <>
  struct CellTypeInformation<CORE::FE::CellType::point1>
  {
    static constexpr auto name = "POINT1";
    static constexpr int dim = 0;
    static constexpr int num_nodes = 1;
    static constexpr int num_faces = 0;
  };

  template <>
  struct CellTypeInformation<CORE::FE::CellType::nurbs2>
  {
    static constexpr auto name = "NURBS2";
    static constexpr int dim = 1;
    static constexpr int num_nodes = 2;
    static constexpr int num_faces = 2;
  };

  template <>
  struct CellTypeInformation<CORE::FE::CellType::nurbs3>
  {
    static constexpr auto name = "NURBS3";
    static constexpr int dim = 1;
    static constexpr int num_nodes = 3;
    static constexpr int num_faces = 2;
  };

  template <>
  struct CellTypeInformation<CORE::FE::CellType::nurbs4>
  {
    static constexpr auto name = "NURBS4";
    static constexpr int dim = 2;
    static constexpr int num_nodes = 4;
    static constexpr int num_faces = 4;
  };

  template <>
  struct CellTypeInformation<CORE::FE::CellType::nurbs9>
  {
    static constexpr auto name = "NURBS9";
    static constexpr int dim = 2;
    static constexpr int num_nodes = 9;
    static constexpr int num_faces = 4;
  };

  template <>
  struct CellTypeInformation<CORE::FE::CellType::nurbs8>
  {
    static constexpr auto name = "NURBS8";
    static constexpr int dim = 3;
    static constexpr int num_nodes = 8;
    static constexpr int num_faces = 6;
  };

  template <>
  struct CellTypeInformation<CORE::FE::CellType::nurbs27>
  {
    static constexpr auto name = "NURBS27";
    static constexpr int dim = 3;
    static constexpr int num_nodes = 27;
    static constexpr int num_faces = 6;
  };

  template <>
  struct CellTypeInformation<CORE::FE::CellType::max_distype>
  {
    static constexpr auto name = "MAX_DISTYPE";
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
  using is_nurbs = std::bool_constant<celltype == CORE::FE::CellType::nurbs27>;

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

  template <CORE::FE::CellType celltype>
  using use_lagrange_shapefnct =
      std::bool_constant<is_tet<celltype>::value || is_hex<celltype>::value ||
                         is_wedge<celltype>::value || is_pyramid<celltype>::value ||
                         is_quad<celltype>::value || is_tri<celltype>::value ||
                         is_line<celltype>::value>;
  /// @}
}  // namespace CORE::FE::DETAILS

BACI_NAMESPACE_CLOSE

#endif