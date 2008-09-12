/*!
\file condif2_utils.cpp

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "condif2_utils.H"


/*-------------------------------------------------------------------------*
 | check for higher order derivatives for shape functions          vg 08/08|
 *-------------------------------------------------------------------------*/
bool SCATRA::is2DHigherOrderElement(const DRT::Element::DiscretizationType& distype)
{
  bool hoel = true;
  switch (distype)
  {
    case DRT::Element::quad4:
    case DRT::Element::tri3:
      hoel = false;
      break;
    case DRT::Element::quad8:
    case DRT::Element::quad9:
    case DRT::Element::tri6:
      hoel = true;
      break;
    default:
      dserror("distype unknown!");
  }
  return hoel;
}


/*----------------------------------------------------------------------*
 |  get optimal gaussrule for discretization type               vg 08/08|
 *----------------------------------------------------------------------*/
DRT::UTILS::GaussRule2D SCATRA::get2DOptimalGaussrule
(const DRT::Element::DiscretizationType& distype)
{
  DRT::UTILS::GaussRule2D rule = DRT::UTILS::intrule2D_undefined;
  switch (distype)
  {
    case DRT::Element::quad4:
      rule = DRT::UTILS::intrule_quad_4point;
      break;
    case DRT::Element::quad8: case DRT::Element::quad9:
      rule = DRT::UTILS::intrule_quad_9point;
      break;
    case DRT::Element::tri3:
      rule = DRT::UTILS::intrule_tri_3point;
      break;
    case DRT::Element::tri6:
      rule = DRT::UTILS::intrule_tri_6point;
      break;
    default:
      dserror("unknown number of nodes for gaussrule initialization");
  }
  return rule;
}


#endif  // #ifdef CCADISCRET
