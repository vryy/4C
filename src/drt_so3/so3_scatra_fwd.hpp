/*!----------------------------------------------------------------------
\file so3_scatra_fwd.hpp

<pre>
   Maintainer: Cristobal Bertoglio
               bertoglio@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
</pre>

*----------------------------------------------------------------------*/

// template classes
template class DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex8,DRT::Element::hex8>;
template class DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_tet4,DRT::Element::tet4>;
template class DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_tet10,DRT::Element::tet10>;
template class DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_weg6,DRT::Element::wedge6>;
// template class DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::NStet5,DRT::Element::tet4>; // Require minor modifications in nstet5 element
