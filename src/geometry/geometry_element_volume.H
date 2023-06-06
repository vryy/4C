/*----------------------------------------------------------------------*/
/*! \file

\brief computes element volume

\level 2

*----------------------------------------------------------------------*/


#ifndef GEOMETRY_ELEMENT_VOLUME_H
#define GEOMETRY_ELEMENT_VOLUME_H


#include "discretization_fem_general_utils_gausspoints.H"
#include "discretization_fem_general_utils_fem_shapefunctions.H"
#include "cut_kernel.H"

namespace GEO
{
  //! calculates the length of an element in given configuration
  template <DRT::Element::DiscretizationType distype, class matrixtype>
  double ElementLengthT(const matrixtype& xyze)  ///> xyze nsd = 3 coords, number of nodes
  {
    // gaussian points
    static const int numnode = CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
    const CORE::DRT::UTILS::GaussIntegration intpoints(distype);

    LINALG::Matrix<1, 1> eleCoord;
    LINALG::Matrix<1, numnode> deriv;
    LINALG::Matrix<1, 3> xjm;

    double length = 0.0;

    // integration loop
    for (int iquad = 0; iquad < intpoints.NumPoints(); ++iquad)
    {
      // coordinates of the current integration point in element coordinates \xi
      eleCoord(0) = intpoints.Point(iquad)[0];

      // shape functions and their first derivatives
      CORE::DRT::UTILS::shape_function_1D_deriv1(deriv, eleCoord(0), distype);

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
  template <DRT::Element::DiscretizationType distype, class matrixtype>
  double ElementAreaT(const matrixtype& xyze)  ///> xyze nsd = 3 coords, number of nodes
  {
    // gaussian points
    static const int numnode = CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
    const CORE::DRT::UTILS::GaussIntegration intpoints(distype);

    LINALG::Matrix<2, 1> eleCoord;
    LINALG::Matrix<2, numnode> deriv;
    LINALG::Matrix<2, 3> xjm;
    LINALG::Matrix<2, 2> xjm_xjmt;

    double area = 0.0;

    // integration loop
    for (int iquad = 0; iquad < intpoints.NumPoints(); ++iquad)
    {
      // coordinates of the current integration point in element coordinates \xi
      eleCoord(0) = intpoints.Point(iquad)[0];
      eleCoord(1) = intpoints.Point(iquad)[1];

      // shape functions and their first derivatives
      CORE::DRT::UTILS::shape_function_2D_deriv1(deriv, eleCoord(0), eleCoord(1), distype);

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
  template <DRT::Element::DiscretizationType distype, class matrixtype>
  double ElementVolumeT(const matrixtype& xyze  ///> xyze nsd = 3 coords, number of nodes)
  )
  {
    // number of nodes for element
    const int numnode = CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
    // gaussian points
    const CORE::DRT::UTILS::GaussIntegration intpoints(distype);

    LINALG::Matrix<3, 1> eleCoord;
    LINALG::Matrix<3, numnode> deriv;
    LINALG::Matrix<3, 3> xjm;

    double vol = 0.0;

    // integration loop
    for (int iquad = 0; iquad < intpoints.NumPoints(); ++iquad)
    {
      // coordinates of the current integration point in element coordinates \xi
      eleCoord(0) = intpoints.Point(iquad)[0];
      eleCoord(1) = intpoints.Point(iquad)[1];
      eleCoord(2) = intpoints.Point(iquad)[2];

      // shape functions and their first derivatives
      CORE::DRT::UTILS::shape_function_3D_deriv1(
          deriv, eleCoord(0), eleCoord(1), eleCoord(2), distype);

      // get transposed of the jacobian matrix d x / d \xi
      xjm = 0;

      for (int inode = 0; inode < numnode; ++inode)
        for (int i = 0; i < 3; ++i)
          for (int j = 0; j < 3; ++j) xjm(i, j) += deriv(i, inode) * xyze(j, inode);

      const double det = xjm.Determinant();
      const double fac = intpoints.Weight(iquad) * det;

      if (det <= 0.0) dserror("NEGATIVE JACOBIAN DETERMINANT: %g", det);

      vol += fac;
    }  // end loop over gauss points
    return vol;
  }

  /** \brief calculates the length of a edge element in given configuration
   *
   *  \params distype (in) : discretization type of the given element
   *  \params xyze    (in) : spatial coordinates of the elememnt nodes
   *                         (row = dim, col = number of nodes)
   *
   *  \author hiermeier \date 11/16 */
  template <class matrixtype>
  double ElementLength(const DRT::Element::DiscretizationType& distype, const matrixtype& xyze)
  {
    if (distype != DRT::Element::line2 or xyze.N() != 2)
      dserror("Currently only line2 elements are supported!");

    // calculate the distance between the two given nodes and return
    // the value
    LINALG::Matrix<3, 1> d(&xyze(0, 0));
    const LINALG::Matrix<3, 1> x1(&xyze(0, 1));

    d.Update(1.0, x1, -1.0);
    return d.Norm2();
  }


  /** \brief calculates the area of a surface element in given configuration
   *
   *  \author hiermeier \date 11/16 */
  template <class matrixtype>
  double ElementArea(const DRT::Element::DiscretizationType distype, const matrixtype& xyze)
  {
    switch (distype)
    {
      // --- 2-D boundary elements
      case DRT::Element::line2:
        return ElementLength(distype, xyze);
      case DRT::Element::line3:
        return ElementLengthT<DRT::Element::line3>(xyze);
      // --- 3-D boundary elements
      case DRT::Element::tri3:
        return ElementAreaT<DRT::Element::tri3>(xyze);
      case DRT::Element::tri6:
        return ElementAreaT<DRT::Element::tri6>(xyze);
      case DRT::Element::quad4:
        return ElementAreaT<DRT::Element::quad4>(xyze);
      case DRT::Element::quad8:
        return ElementAreaT<DRT::Element::quad8>(xyze);
      case DRT::Element::quad9:
        return ElementAreaT<DRT::Element::quad9>(xyze);
      default:
        dserror("Unsupported surface element type!");
        exit(EXIT_FAILURE);
    }

    return -1.0;
  }

  //! calculates the volume of an element in given configuration          u.may
  template <class matrixtype>
  double ElementVolume(const DRT::Element::DiscretizationType distype,
      const matrixtype& xyze  ///> xyze nsd = 3 coords, number of nodes
  )
  {
    switch (distype)
    {
      // --- 1-D elements -----------------------------------------------------
      case DRT::Element::line2:
        return ElementLength(distype, xyze);
      // --- 2-D elements -----------------------------------------------------
      case DRT::Element::quad4:
      case DRT::Element::tri3:
        return ElementArea(distype, xyze);
      // --- 3-D elements -----------------------------------------------------
      case DRT::Element::hex8:
        return ElementVolumeT<DRT::Element::hex8>(xyze);
      case DRT::Element::hex20:
        return ElementVolumeT<DRT::Element::hex20>(xyze);
      case DRT::Element::hex27:
        return ElementVolumeT<DRT::Element::hex27>(xyze);
      case DRT::Element::tet4:
        return ElementVolumeT<DRT::Element::tet4>(xyze);
      case DRT::Element::tet10:
        return ElementVolumeT<DRT::Element::tet10>(xyze);
      case DRT::Element::wedge6:
        return ElementVolumeT<DRT::Element::wedge6>(xyze);
      case DRT::Element::wedge15:
        return ElementVolumeT<DRT::Element::wedge15>(xyze);
      case DRT::Element::pyramid5:
        return ElementVolumeT<DRT::Element::pyramid5>(xyze);
      default:
        dserror("add you distype here");
    }
    return -1.0;
  }

}  // namespace GEO


#endif  // GEOMETRY_ELEMENT_VOLUME_H
