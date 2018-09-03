/*
 * so3_thermo_fwd.hpp
 *
 *  forward declarations of discretisation types
 *
   Maintainer: Alexander Seitz
               seitz@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15271
 *
 */

// template classes
template class DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>;
template class DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex8fbar, DRT::Element::hex8>;
template class DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex27, DRT::Element::hex27>;
template class DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex20, DRT::Element::hex20>;
template class DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>;
template class DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_tet10, DRT::Element::tet10>;
template class DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::NURBS::So_nurbs27, DRT::Element::nurbs27>;
