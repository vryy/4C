/*----------------------------------------------------------------------*/
/*! \file

\brief computes element volume

\level 2

\maintainer Martin Kronbichler
*----------------------------------------------------------------------*/


#ifndef ELEMENT_VOLUME_H
#define ELEMENT_VOLUME_H


#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_cut/cut_kernel.H"

namespace GEO
{
  //! calculates the volume of an element in given configuration          u.may
  template <DRT::Element::DiscretizationType DISTYPE, class M1>
  double ElementVolumeT(const M1& xyze  ///> xyze nsd = 3 coords, number of nodes)
  )
  {
    // number of nodes for element
    const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

    DRT::UTILS::GaussRule3D gaussrule = DRT::UTILS::intrule3D_undefined;

    switch (DISTYPE)
    {
      case DRT::Element::hex8:
      case DRT::Element::hex20:
      case DRT::Element::hex27:
      {
        gaussrule = DRT::UTILS::intrule_hex_8point;
        break;
      }
      case DRT::Element::tet4:
      case DRT::Element::tet10:
      {
        gaussrule = DRT::UTILS::intrule_tet_4point;
        break;
      }
      case DRT::Element::wedge6:
      case DRT::Element::wedge15:
      {
        gaussrule = DRT::UTILS::intrule_wedge_6point;
        break;
      }
      case DRT::Element::pyramid5:
      {
        gaussrule = DRT::UTILS::intrule_pyramid_8point;
        break;
      }
      default:
        dserror("add your element type here");
    }

    // gaussian points
    const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);

    LINALG::Matrix<3, 1> eleCoord;
    LINALG::Matrix<3, numnode> deriv;
    LINALG::Matrix<3, 3> xjm;

    double vol = 0.0;

    // integration loop
    for (int iquad = 0; iquad < intpoints.nquad; ++iquad)
    {
      // coordinates of the current integration point in element coordinates \xi
      eleCoord(0) = intpoints.qxg[iquad][0];
      eleCoord(1) = intpoints.qxg[iquad][1];
      eleCoord(2) = intpoints.qxg[iquad][2];

      // shape functions and their first derivatives
      DRT::UTILS::shape_function_3D_deriv1(deriv, eleCoord(0), eleCoord(1), eleCoord(2), DISTYPE);

      // get transposed of the jacobian matrix d x / d \xi
      xjm = 0;

      for (int inode = 0; inode < numnode; ++inode)
        for (int i = 0; i < 3; ++i)
          for (int j = 0; j < 3; ++j) xjm(i, j) += deriv(i, inode) * xyze(j, inode);

      const double det = xjm.Determinant();
      const double fac = intpoints.qwgt[iquad] * det;

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
  template <class M1>
  double ElementLength(const DRT::Element::DiscretizationType& distype, const M1& xyze)
  {
    if (distype != DRT::Element::line2 or xyze.N() != 2)
      run_time_error("Currently only line2 elements are supported!");

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
  template <class M1>
  double ElementArea(const DRT::Element::DiscretizationType distype, const M1& xyze)
  {
    switch (distype)
    {
      // --- 2-D boundary elements
      case DRT::Element::line2:
        return ElementLength(distype, xyze);
      // --- 3-D boundary elements
      case DRT::Element::quad4:
      {
        LINALG::Matrix<3, 3> xyze_tri3;
        double area = 0.0;
        for (unsigned id = 0; id < 2; ++id)
        {
          CUT::KERNEL::SplitQuad4IntoTri3(xyze, id, xyze_tri3);
          area += CUT::KERNEL::getAreaTri(xyze_tri3);
        }
        return area;
      }
      case DRT::Element::tri3:
      {
        return CUT::KERNEL::getAreaTri(xyze);
      }
      default:
        dserror("Unsupported surface element type!");
        exit(EXIT_FAILURE);
    }

    return -1.0;
  }

  //! calculates the volume of an element in given configuration          u.may
  template <class M1>
  double ElementVolume(const DRT::Element::DiscretizationType distype,
      const M1& xyze  ///> xyze nsd = 3 coords, number of nodes
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


#endif /*ELEMENT_VOLUME_H*/
