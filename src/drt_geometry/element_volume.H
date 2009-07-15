/*!----------------------------------------------------------------------
\file element_volume.H

\brief computes element volume

<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef ELEMENT_VOLUME_H
#define ELEMENT_VOLUME_H

#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_geometry/position_array.H"

namespace GEO
{

  //! calculates the volume of an element in given configuration          u.may
  template <DRT::Element::DiscretizationType DISTYPE, class M1>
    static inline double ElementVolumeT(
          const M1&                     xyze   ///> xyze nsd = 3 coords, number of nodes)
          )
  {
    // number of nodes for element
    const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

    DRT::UTILS::GaussRule3D gaussrule = DRT::UTILS::intrule3D_undefined;

    switch (DISTYPE)
    {
      case DRT::Element::hex8: case DRT::Element::hex20: case DRT::Element::hex27:
      {
        gaussrule = DRT::UTILS::intrule_hex_8point;
        break;
      }
      case DRT::Element::tet4: case DRT::Element::tet10:
      {
        gaussrule = DRT::UTILS::intrule_tet_4point;
        break;
      }
      default:
        dserror("add your element type here...");
    }

    // gaussian points
    const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);

    double vol = 0.0;
    // integration loop
    for (int iquad=0; iquad<intpoints.nquad; ++iquad)
    {
      // coordinates of the current integration point in element coordinates \xi
      static LINALG::Matrix<3,1> eleCoord;
      eleCoord(0) = intpoints.qxg[iquad][0];
      eleCoord(1) = intpoints.qxg[iquad][1];
      eleCoord(2) = intpoints.qxg[iquad][2];

      // shape functions and their first derivatives
      static LINALG::Matrix<3,numnode> deriv;
      DRT::UTILS::shape_function_3D_deriv1(deriv,eleCoord(0),eleCoord(1),eleCoord(2),DISTYPE);

      // get transposed of the jacobian matrix d x / d \xi
      static LINALG::Matrix<3,3> xjm;
      xjm.Clear();

      for (int inode = 0; inode < numnode; ++inode)
        for (int i = 0; i < 3; ++i)
          for (int j = 0; j < 3; ++j)
            xjm(i,j) += deriv(i,inode) * xyze(j,inode);

      const double det = xjm.Determinant();
      const double fac = intpoints.qwgt[iquad]*det;

      if (det <= 0.0)
        dserror("NEGATIVE JACOBIAN DETERMINANT: %f", det);

      vol += fac;
    } // end loop over gauss points
    return vol;
  }



  //! calculates the volume of an element in given configuration          u.may
  template <class M1>
    double ElementVolume(
        const DRT::Element::DiscretizationType    distype,
        const M1&                                 xyze       ///> xyze nsd = 3 coords, number of nodes
        )
  {
    switch(distype)
    {
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
      default:
        dserror("add you distype here...");

    }
    return -1.0;
  }
}  // namespace GEO


#endif /*ELEMENT_VOLUME_H*/
#endif /*CCADISCRET*/
