/*!
\file condif3_utils.H

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef CONDIF3_UTILS_H
#define CONDIF3_UTILS_H


#include "condif3.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_lib/drt_discret.H"


namespace SCATRA
{

  //! check for higher order derivatives for shape functions
  bool is3DHigherOrderElement
  (
   const DRT::Element::DiscretizationType& distype ///< distype of the element
  );


  //! get optimal gaussrule for discretization type
  DRT::UTILS::GaussRule3D get3DOptimalGaussrule
  (
   const DRT::Element::DiscretizationType& distype ///< distype of the element
  );


  //! decide, whether second derivatives are needed  (template version)
  template <DRT::Element::DiscretizationType  DISTYPE>
  bool is3DHigherOrderElement()
  {
    switch (DISTYPE)
    {
    case DRT::Element::hex8:
    case DRT::Element::tet4:
      return false;
      break;
    case DRT::Element::hex20:
    case DRT::Element::hex27:
    case DRT::Element::tet10:
      return true;
      break;
    case DRT::Element::wedge6:
    case DRT::Element::pyramid5:
    case DRT::Element::wedge15:
      //!!!TODO:  wedge und pyramid have 2nd derivatives!!!!!!
      dserror("wedges and pyramids have second derivatives!");
      break;
    default:
      dserror("distype unknown!");
    }
    return true;
  }


  //! get optimal gaussrule for discretization type (template version)
  template <DRT::Element::DiscretizationType  DISTYPE>
  DRT::UTILS::GaussRule3D get3DOptimalGaussrule()
  {
    DRT::UTILS::GaussRule3D rule = DRT::UTILS::intrule3D_undefined;
    switch (DISTYPE)
    {
      case DRT::Element::hex8:
        rule = DRT::UTILS::intrule_hex_8point;
        break;
      case DRT::Element::hex20: case DRT::Element::hex27:
        rule = DRT::UTILS::intrule_hex_27point;
        break;
      case DRT::Element::tet4:
        rule = DRT::UTILS::intrule_tet_4point;
        break;
      case DRT::Element::tet10:
        rule = DRT::UTILS::intrule_tet_5point;
        break;
      default:
        dserror("unknown number of nodes for gaussrule initialization");
    }
    return rule;
  }


} // namespace SCATRA

#endif
#endif  // #ifdef CCADISCRET
