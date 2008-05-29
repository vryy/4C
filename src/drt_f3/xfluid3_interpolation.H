/*----------------------------------------------------------------------*/
/*!
\file xfluid3_interpolation.H

\brief rotuines to assemble element stiffness matrix and vector

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
 */
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET


#ifndef XFLUID3_INTERPOLATION_H
#define XFLUID3_INTERPOLATION_H

#include "../drt_lib/drt_element.H"


//! support routines for the xfluid3 element
namespace XFLUID
{
  
  /*!
   * \brief Returns the distype for stress interpolation
   *        for each parent element discretization type in 3d
   * 
   * \return distype for stress interpolation
   */
  template<DRT::Element::DiscretizationType veldistype>
  struct StressInterpolation3D {};
  template<> struct StressInterpolation3D<DRT::Element::tet4>  {static const DRT::Element::DiscretizationType distype = DRT::Element::point1;};
  template<> struct StressInterpolation3D<DRT::Element::hex8>  {static const DRT::Element::DiscretizationType distype = DRT::Element::hex8;};
  template<> struct StressInterpolation3D<DRT::Element::hex20> {static const DRT::Element::DiscretizationType distype = DRT::Element::hex20;};
  template<> struct StressInterpolation3D<DRT::Element::hex27> {static const DRT::Element::DiscretizationType distype = DRT::Element::hex8;};
  
  
  /*!
   * \brief Returns the distype for discont. pressure interpolation
   *        for each parent element discretization type in 3d
   * 
   * \return distype for discont. pressure interpolation
   */
  template<DRT::Element::DiscretizationType veldistype>
  struct DiscPressureInterpolation3D {};
  template<> struct DiscPressureInterpolation3D<DRT::Element::tet4>  {static const DRT::Element::DiscretizationType distype = DRT::Element::point1;};
  template<> struct DiscPressureInterpolation3D<DRT::Element::hex8>  {static const DRT::Element::DiscretizationType distype = DRT::Element::hex8;};
  template<> struct DiscPressureInterpolation3D<DRT::Element::hex20> {static const DRT::Element::DiscretizationType distype = DRT::Element::hex20;};
  template<> struct DiscPressureInterpolation3D<DRT::Element::hex27> {static const DRT::Element::DiscretizationType distype = DRT::Element::hex8;};
  
  
  /*!
   * \brief Returns the distype for stress interpolation
   *        for each parent element discretization type in 3d
   * 
   * \return distype for stress interpolation
   */  
  inline DRT::Element::DiscretizationType getStressInterpolationType3D(DRT::Element::DiscretizationType veldistype)
  {
    DRT::Element::DiscretizationType distype = DRT::Element::dis_none;
    switch (veldistype)
    {
      case DRT::Element::tet4:
        distype = StressInterpolation3D<DRT::Element::tet4>::distype;
        break;
      case DRT::Element::hex8:
        distype = StressInterpolation3D<DRT::Element::hex8>::distype;
        break;
      case DRT::Element::hex20:
        distype = StressInterpolation3D<DRT::Element::hex20>::distype;
        break;
      case DRT::Element::hex27:
        distype = StressInterpolation3D<DRT::Element::hex27>::distype;
        break;
      default:
        dserror("no stressdistype defined for velocity distype");
    }
    return distype;
  }
  
  /*!
   * \brief Returns the distype for discont. pressure interpolation
   *        for each parent element discretization type in 3d
   * 
   * \return distype for discont. pressure interpolation
   */
  template <DRT::Element::DiscretizationType veldistype>
  inline DRT::Element::DiscretizationType getDiscPressureInterpolationType3D()
  {
    DRT::Element::DiscretizationType distype = DRT::Element::dis_none;
    switch (veldistype)
    {
      case DRT::Element::tet4:
        distype = DiscPressureInterpolation3D<DRT::Element::tet4>::distype;
        break;
      case DRT::Element::hex8:
        distype = DiscPressureInterpolation3D<DRT::Element::hex8>::distype;
        break;
      case DRT::Element::hex20:
        distype = DiscPressureInterpolation3D<DRT::Element::hex20>::distype;
        break;
      case DRT::Element::hex27:
        distype = DiscPressureInterpolation3D<DRT::Element::hex27>::distype;
        break;
      default:
        dserror("no discpressuredistype defined for velocity distype");
    }
    return distype;
  }
  
}


#endif /*XFLUID3_INTERPOLATION_H*/
#endif
