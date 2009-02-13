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

#ifndef ELEMENT_VOLUME_H_
#define ELEMENT_VOLUME_H_

#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_geometry/position_array.H"

namespace GEO
{


  //! calculates the volume of an element in initial configuration      a.ger
  template <DRT::Element::DiscretizationType DISTYPE>
    static inline double ElementVolumeT(
          const DRT::Element&           ele)
  {
    if (ele.Shape() != DISTYPE)
      dserror("mismatch between element shape and template parameter DISTYPE! This is a bug!");
  
    // number of nodes for element
    const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

    // get node coordinates of the current element
    static LINALG::Matrix<3,numnode> xyze;
    GEO::fillInitialPositionArray<DISTYPE>(&ele, xyze);
    
    // physical element volume 
    double vol = 0.0;
    
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
      default:
        dserror("add your element type here...");
    }
    
    // integration points
    const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);
  
    // integration loop
    for (int iquad=0; iquad<intpoints.nquad; ++iquad)
    {
      // coordinates of the current integration point in element coordinates \xi
      static LINALG::Matrix<3,1> posXiDomain;
      posXiDomain(0) = intpoints.qxg[iquad][0];
      posXiDomain(1) = intpoints.qxg[iquad][1];
      posXiDomain(2) = intpoints.qxg[iquad][2];
      
      // shape functions and their first derivatives
      static LINALG::Matrix<3,numnode> deriv;
      DRT::UTILS::shape_function_3D_deriv1(deriv,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
  
      // get transposed of the jacobian matrix d x / d \xi
      static LINALG::Matrix<3,3> xjm;
      xjm.MultiplyNT(deriv,xyze);
  
      const double det = xjm.Determinant();
      const double fac = intpoints.qwgt[iquad]*det;
  
      if (det <= 0.0)
        dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", ele.Id(), det);
     
      vol += fac;
        
    }
    return vol;
  }

  
  //! calculates the volume of an element in initial configuration          a.ger
  double ElementVolume(
          const DRT::Element&           ele)
  {
    switch (ele.Shape())
    {
      case DRT::Element::hex8:
        return ElementVolumeT<DRT::Element::hex8>(ele);
      case DRT::Element::hex20:
        return ElementVolumeT<DRT::Element::hex20>(ele);
      case DRT::Element::hex27:
        return ElementVolumeT<DRT::Element::hex27>(ele);
      case DRT::Element::tet4:
        return ElementVolumeT<DRT::Element::tet4>(ele);
      case DRT::Element::tet10:
        return ElementVolumeT<DRT::Element::tet10>(ele);
      default:
        dserror("add you distype here...");
      exit(1);
    }
  }

}  // namespace GEO


#endif /*ELEMENT_VOLUME_H_*/
#endif /*CCADISCRET*/
