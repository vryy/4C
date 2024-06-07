/*----------------------------------------------------------------------*/
/*! \file

\brief computes element volume

\level 2

*----------------------------------------------------------------------*/


#ifndef FOUR_C_FEM_GEOMETRY_ELEMENT_VOLUME_HPP
#define FOUR_C_FEM_GEOMETRY_ELEMENT_VOLUME_HPP


#include "4C_config.hpp"

#include "4C_cut_kernel.hpp"
#include "4C_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_fem_general_utils_gausspoints.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::Geo
{
  //! calculates the length of an element in given configuration
  template <Core::FE::CellType distype, class matrixtype>
  double ElementLengthT(const matrixtype& xyze)  ///> xyze nsd = 3 coords, number of nodes
  {
    // gaussian points
    static constexpr int numnode = Core::FE::num_nodes<distype>;
    const Core::FE::GaussIntegration intpoints(distype);

    Core::LinAlg::Matrix<1, 1> eleCoord;
    Core::LinAlg::Matrix<1, numnode> deriv;
    Core::LinAlg::Matrix<1, 3> xjm;

    double length = 0.0;

    // integration loop
    for (int iquad = 0; iquad < intpoints.NumPoints(); ++iquad)
    {
      // coordinates of the current integration point in element coordinates \xi
      eleCoord(0) = intpoints.Point(iquad)[0];

      // shape functions and their first derivatives
      Core::FE::shape_function_1D_deriv1(deriv, eleCoord(0), distype);

      // get transposed of the jacobian matrix d x / d \xi
      xjm = 0;

      for (int inode = 0; inode < numnode; ++inode)
        for (int i = 0; i < 1; ++i)
          for (int j = 0; j < 3; ++j) xjm(i, j) += deriv(i, inode) * xyze(j, inode);

      const double fac = intpoints.Weight(iquad) * xjm.Norm2();

      length += fac;
    }  // end loop over gauss points

    return length;
  }

  //! calculates the area of an element in given configuration
  template <Core::FE::CellType distype, class matrixtype>
  double ElementAreaT(const matrixtype& xyze)  ///> xyze nsd = 3 coords, number of nodes
  {
    // gaussian points
    static constexpr int numnode = Core::FE::num_nodes<distype>;
    const Core::FE::GaussIntegration intpoints(distype);

    Core::LinAlg::Matrix<2, 1> eleCoord;
    Core::LinAlg::Matrix<2, numnode> deriv;
    Core::LinAlg::Matrix<2, 3> xjm;
    Core::LinAlg::Matrix<2, 2> xjm_xjmt;

    double area = 0.0;

    // integration loop
    for (int iquad = 0; iquad < intpoints.NumPoints(); ++iquad)
    {
      // coordinates of the current integration point in element coordinates \xi
      eleCoord(0) = intpoints.Point(iquad)[0];
      eleCoord(1) = intpoints.Point(iquad)[1];

      // shape functions and their first derivatives
      Core::FE::shape_function_2D_deriv1(deriv, eleCoord(0), eleCoord(1), distype);

      // get transposed of the jacobian matrix d x / d \xi
      xjm = 0;

      for (int inode = 0; inode < numnode; ++inode)
        for (int i = 0; i < 2; ++i)
          for (int j = 0; j < 3; ++j) xjm(i, j) += deriv(i, inode) * xyze(j, inode);

      xjm_xjmt.MultiplyNT<3>(xjm, xjm);

      const double det = xjm_xjmt.Determinant();
      const double fac = intpoints.Weight(iquad) * std::sqrt(det);

      area += fac;
    }  // end loop over gauss points

    return area;
  }

  //! calculates the volume of an element in given configuration          u.may
  template <Core::FE::CellType distype, class matrixtype>
  double ElementVolumeT(const matrixtype& xyze  ///> xyze nsd = 3 coords, number of nodes)
  )
  {
    // number of nodes for element
    const int numnode = Core::FE::num_nodes<distype>;
    // gaussian points
    const Core::FE::GaussIntegration intpoints(distype);

    Core::LinAlg::Matrix<3, 1> eleCoord;
    Core::LinAlg::Matrix<3, numnode> deriv;
    Core::LinAlg::Matrix<3, 3> xjm;

    double vol = 0.0;

    // integration loop
    for (int iquad = 0; iquad < intpoints.NumPoints(); ++iquad)
    {
      // coordinates of the current integration point in element coordinates \xi
      eleCoord(0) = intpoints.Point(iquad)[0];
      eleCoord(1) = intpoints.Point(iquad)[1];
      eleCoord(2) = intpoints.Point(iquad)[2];

      // shape functions and their first derivatives
      Core::FE::shape_function_3D_deriv1(deriv, eleCoord(0), eleCoord(1), eleCoord(2), distype);

      // get transposed of the jacobian matrix d x / d \xi
      xjm = 0;

      for (int inode = 0; inode < numnode; ++inode)
        for (int i = 0; i < 3; ++i)
          for (int j = 0; j < 3; ++j) xjm(i, j) += deriv(i, inode) * xyze(j, inode);

      const double det = xjm.Determinant();
      const double fac = intpoints.Weight(iquad) * det;

      if (det <= 0.0) FOUR_C_THROW("NEGATIVE JACOBIAN DETERMINANT: %g", det);

      vol += fac;
    }  // end loop over gauss points
    return vol;
  }

  /** \brief calculates the length of a edge element in given configuration
   *
   *  \params distype (in) : discretization type of the given element
   *  \params xyze    (in) : spatial coordinates of the elememnt nodes
   *                         (row = dim, col = number of nodes)
   */
  template <class matrixtype>
  double ElementLength(const Core::FE::CellType& distype, const matrixtype& xyze)
  {
    if (distype != Core::FE::CellType::line2 or xyze.numCols() != 2)
      FOUR_C_THROW("Currently only line2 elements are supported!");

    // calculate the distance between the two given nodes and return
    // the value
    Core::LinAlg::Matrix<3, 1> d(&xyze(0, 0));
    const Core::LinAlg::Matrix<3, 1> x1(&xyze(0, 1));

    d.Update(1.0, x1, -1.0);
    return d.Norm2();
  }


  /** \brief calculates the area of a surface element in given configuration
   */
  template <class matrixtype>
  double ElementArea(const Core::FE::CellType distype, const matrixtype& xyze)
  {
    switch (distype)
    {
      // --- 2-D boundary elements
      case Core::FE::CellType::line2:
        return ElementLength(distype, xyze);
      case Core::FE::CellType::line3:
        return ElementLengthT<Core::FE::CellType::line3>(xyze);
      // --- 3-D boundary elements
      case Core::FE::CellType::tri3:
        return ElementAreaT<Core::FE::CellType::tri3>(xyze);
      case Core::FE::CellType::tri6:
        return ElementAreaT<Core::FE::CellType::tri6>(xyze);
      case Core::FE::CellType::quad4:
        return ElementAreaT<Core::FE::CellType::quad4>(xyze);
      case Core::FE::CellType::quad8:
        return ElementAreaT<Core::FE::CellType::quad8>(xyze);
      case Core::FE::CellType::quad9:
        return ElementAreaT<Core::FE::CellType::quad9>(xyze);
      default:
        FOUR_C_THROW("Unsupported surface element type!");
        exit(EXIT_FAILURE);
    }

    return -1.0;
  }

  //! calculates the volume of an element in given configuration          u.may
  template <class matrixtype>
  double ElementVolume(const Core::FE::CellType distype,
      const matrixtype& xyze  ///> xyze nsd = 3 coords, number of nodes
  )
  {
    switch (distype)
    {
      // --- 1-D elements -----------------------------------------------------
      case Core::FE::CellType::line2:
        return ElementLength(distype, xyze);
      // --- 2-D elements -----------------------------------------------------
      case Core::FE::CellType::quad4:
      case Core::FE::CellType::tri3:
        return ElementArea(distype, xyze);
      // --- 3-D elements -----------------------------------------------------
      case Core::FE::CellType::hex8:
        return ElementVolumeT<Core::FE::CellType::hex8>(xyze);
      case Core::FE::CellType::hex20:
        return ElementVolumeT<Core::FE::CellType::hex20>(xyze);
      case Core::FE::CellType::hex27:
        return ElementVolumeT<Core::FE::CellType::hex27>(xyze);
      case Core::FE::CellType::tet4:
        return ElementVolumeT<Core::FE::CellType::tet4>(xyze);
      case Core::FE::CellType::tet10:
        return ElementVolumeT<Core::FE::CellType::tet10>(xyze);
      case Core::FE::CellType::wedge6:
        return ElementVolumeT<Core::FE::CellType::wedge6>(xyze);
      case Core::FE::CellType::wedge15:
        return ElementVolumeT<Core::FE::CellType::wedge15>(xyze);
      case Core::FE::CellType::pyramid5:
        return ElementVolumeT<Core::FE::CellType::pyramid5>(xyze);
      default:
        FOUR_C_THROW("add you distype here");
    }
    return -1.0;
  }

}  // namespace Core::Geo


FOUR_C_NAMESPACE_CLOSE

#endif
