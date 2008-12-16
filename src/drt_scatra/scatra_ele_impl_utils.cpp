/*!
\file scatra_ele_impl_utils.cpp

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "scatra_ele_impl_utils.H"


/*-------------------------------------------------------------------------*
 | check for higher order derivatives for shape functions         gjb 08/08|
 *-------------------------------------------------------------------------*/
bool SCATRA::is3DHigherOrderElement(const DRT::Element::DiscretizationType& distype)
{
  bool hoel = true;
  switch (distype)
  {
    case DRT::Element::hex8:
    case DRT::Element::tet4:
      hoel = false;
      break;
    case DRT::Element::hex20:
    case DRT::Element::hex27:
    case DRT::Element::tet10:
      hoel = true;
      break;
    case DRT::Element::wedge6:
    case DRT::Element::pyramid5:
    case DRT::Element::wedge15:
      //!!!TODO:  wedge und pyramid have 2nd derivatives!!!!!!!!!!!!!!!!!!!!!!!!
      dserror("wedges and pyramids have second derivatives!");
      break;
    default:
      dserror("distype unknown!");
  }
  return hoel;
}


/*----------------------------------------------------------------------*
 |  get optimal gaussrule for discretization type              gjb 08/08|
 *----------------------------------------------------------------------*/
DRT::UTILS::GaussRule3D SCATRA::get3DOptimalGaussrule
(const DRT::Element::DiscretizationType& distype)
{
  DRT::UTILS::GaussRule3D rule = DRT::UTILS::intrule3D_undefined;
  switch (distype)
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


#endif  // #ifdef CCADISCRET
