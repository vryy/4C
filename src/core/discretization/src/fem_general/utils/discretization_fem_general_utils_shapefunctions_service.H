/*----------------------------------------------------------------------*/
/*! \file

 \brief provide service routines involving shape functions

\level 2
 */
/*----------------------------------------------------------------------*/

#ifndef DISCRETIZATION_FEM_GENERAL_UTILS_SHAPEFUNCTIONS_SERVICE_H
#define DISCRETIZATION_FEM_GENERAL_UTILS_SHAPEFUNCTIONS_SERVICE_H


#include "lib_element.H"


namespace CORE::DRT::UTILS
{
  /*!

    \brief decide whether second derivatives are zero

    \detailed

    For some element types second derivatives are always zero. They do not have
    to be computed.

    - For elements with linear shape functions second derivatives are 0.
    - For elements with bi- or trilinear shape functions second derivatives
      are not 0, because mixed components (e.g. d^2 N/dxdy) are not zero.
    - For higher order elements (hoele) (P > 1), second derivatives are not 0.

    \author henke
  */
  template <::DRT::Element::DiscretizationType DISTYPE>
  bool secondDerivativesZero()
  {
    bool hoele = true;
    switch (DISTYPE)
    {
      case ::DRT::Element::quad4:
      case ::DRT::Element::quad9:
      case ::DRT::Element::hex8:
      case ::DRT::Element::hex20:
      case ::DRT::Element::hex27:
      case ::DRT::Element::tet10:
      case ::DRT::Element::wedge6:
      case ::DRT::Element::wedge15:
      case ::DRT::Element::pyramid5:
      case ::DRT::Element::nurbs9:
      case ::DRT::Element::nurbs27:
        hoele = true;
        break;
      case ::DRT::Element::tri3:
      case ::DRT::Element::tet4:
        hoele = false;
        break;
      default:
        dserror("DISTYPE unknown!");
        break;
    }
    return hoele;
  }


}  // namespace CORE::DRT::UTILS

#endif  // FEM_GENERAL_UTILS_SHAPEFUNCTIONS_SERVICE_H
