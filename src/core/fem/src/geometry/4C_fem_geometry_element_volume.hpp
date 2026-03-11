// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FEM_GEOMETRY_ELEMENT_VOLUME_HPP
#define FOUR_C_FEM_GEOMETRY_ELEMENT_VOLUME_HPP


#include "4C_config.hpp"

#include "4C_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_fem_general_utils_gausspoints.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::Geo
{
  //! calculates the length of an element in given configuration
  template <Core::FE::CellType distype, class Matrixtype>
  double element_length(const Matrixtype& xyze)
  {
    static constexpr int numnode = Core::FE::num_nodes(distype);
    const Core::FE::GaussIntegration intpoints(distype);

    Core::LinAlg::Matrix<1, 1> eleCoord;
    Core::LinAlg::Matrix<1, numnode> deriv;
    Core::LinAlg::Matrix<1, 3> xjm;

    double length = 0.0;

    for (int iquad = 0; iquad < intpoints.num_points(); ++iquad)
    {
      // coordinates of the current integration point in element coordinates \xi
      eleCoord(0) = intpoints.point(iquad)[0];

      // shape functions and their first derivatives
      Core::FE::shape_function_1d_deriv1(deriv, eleCoord(0), distype);

      // get transposed of the jacobian matrix d x / d \xi
      xjm = 0;

      for (int inode = 0; inode < numnode; ++inode)
        for (int i = 0; i < 1; ++i)
          for (int j = 0; j < 3; ++j) xjm(i, j) += deriv(i, inode) * xyze(j, inode);

      const double fac = intpoints.weight(iquad) * xjm.norm2();

      length += fac;
    }

    return length;
  }

  //! calculates the area of an element in given configuration
  template <Core::FE::CellType distype, class Matrixtype>
  double element_area(const Matrixtype& xyze)
  {
    static constexpr int numnode = Core::FE::num_nodes(distype);
    const Core::FE::GaussIntegration intpoints(distype);

    Core::LinAlg::Matrix<2, 1> eleCoord;
    Core::LinAlg::Matrix<2, numnode> deriv;
    Core::LinAlg::Matrix<2, 3> xjm;
    Core::LinAlg::Matrix<2, 2> xjm_xjmt;

    double area = 0.0;

    for (int iquad = 0; iquad < intpoints.num_points(); ++iquad)
    {
      // coordinates of the current integration point in element coordinates \xi
      eleCoord(0) = intpoints.point(iquad)[0];
      eleCoord(1) = intpoints.point(iquad)[1];

      // shape functions and their first derivatives
      Core::FE::shape_function_2d_deriv1(deriv, eleCoord(0), eleCoord(1), distype);

      // get transposed of the jacobian matrix d x / d \xi
      xjm = 0;

      for (int inode = 0; inode < numnode; ++inode)
        for (int i = 0; i < 2; ++i)
          for (int j = 0; j < 3; ++j) xjm(i, j) += deriv(i, inode) * xyze(j, inode);

      xjm_xjmt.multiply_nt<3>(xjm, xjm);

      const double det = xjm_xjmt.determinant();
      const double fac = intpoints.weight(iquad) * std::sqrt(det);

      area += fac;
    }

    return area;
  }

  //! calculates the volume of an element in given configuration
  template <Core::FE::CellType distype, class Matrixtype>
  double element_volume(const Matrixtype& xyze)
  {
    const int numnode = Core::FE::num_nodes(distype);
    const Core::FE::GaussIntegration intpoints(distype);

    Core::LinAlg::Matrix<3, 1> eleCoord;
    Core::LinAlg::Matrix<3, numnode> deriv;
    Core::LinAlg::Matrix<3, 3> xjm;

    double vol = 0.0;

    for (int iquad = 0; iquad < intpoints.num_points(); ++iquad)
    {
      // coordinates of the current integration point in element coordinates \xi
      eleCoord(0) = intpoints.point(iquad)[0];
      eleCoord(1) = intpoints.point(iquad)[1];
      eleCoord(2) = intpoints.point(iquad)[2];

      // shape functions and their first derivatives
      Core::FE::shape_function_3d_deriv1(deriv, eleCoord(0), eleCoord(1), eleCoord(2), distype);

      // get transposed of the jacobian matrix d x / d \xi
      xjm = 0;

      for (int inode = 0; inode < numnode; ++inode)
        for (int i = 0; i < 3; ++i)
          for (int j = 0; j < 3; ++j) xjm(i, j) += deriv(i, inode) * xyze(j, inode);

      const double det = xjm.determinant();
      const double fac = intpoints.weight(iquad) * det;

      if (det <= 0.0) FOUR_C_THROW("NEGATIVE JACOBIAN DETERMINANT: {}", det);

      vol += fac;
    }

    return vol;
  }


  //! calculates the volume of an element in given configuration
  template <class Matrixtype>
  double element_volume(const Core::FE::CellType distype, const Matrixtype& xyze)
  {
    switch (distype)
    {
      case Core::FE::CellType::line2:
        return element_length<Core::FE::CellType::line2>(xyze);
      case Core::FE::CellType::line3:
        return element_length<Core::FE::CellType::line3>(xyze);
      case Core::FE::CellType::tri3:
        return element_area<Core::FE::CellType::tri3>(xyze);
      case Core::FE::CellType::tri6:
        return element_area<Core::FE::CellType::tri6>(xyze);
      case Core::FE::CellType::quad4:
        return element_area<Core::FE::CellType::quad4>(xyze);
      case Core::FE::CellType::quad8:
        return element_area<Core::FE::CellType::quad8>(xyze);
      case Core::FE::CellType::quad9:
        return element_area<Core::FE::CellType::quad9>(xyze);
      case Core::FE::CellType::hex8:
        return element_volume<Core::FE::CellType::hex8>(xyze);
      case Core::FE::CellType::hex20:
        return element_volume<Core::FE::CellType::hex20>(xyze);
      case Core::FE::CellType::hex27:
        return element_volume<Core::FE::CellType::hex27>(xyze);
      case Core::FE::CellType::tet4:
        return element_volume<Core::FE::CellType::tet4>(xyze);
      case Core::FE::CellType::tet10:
        return element_volume<Core::FE::CellType::tet10>(xyze);
      case Core::FE::CellType::wedge6:
        return element_volume<Core::FE::CellType::wedge6>(xyze);
      case Core::FE::CellType::wedge15:
        return element_volume<Core::FE::CellType::wedge15>(xyze);
      case Core::FE::CellType::pyramid5:
        return element_volume<Core::FE::CellType::pyramid5>(xyze);
      default:
        FOUR_C_THROW(
            "Element volume calculation is current not implemented for the given cell type.");
    }
  }

}  // namespace Core::Geo


FOUR_C_NAMESPACE_CLOSE

#endif
