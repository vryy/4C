#ifndef BACI_SO3_THERMO_FWD_HPP
#define BACI_SO3_THERMO_FWD_HPP

/*----------------------------------------------------------------------*/
/*! \file

\level 1

\brief forward declarations of discretisation types

 *
 */

// template classes
template class DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex8,
    DRT::Element::DiscretizationType::hex8>;
template class DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex8fbar,
    DRT::Element::DiscretizationType::hex8>;
template class DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex27,
    DRT::Element::DiscretizationType::hex27>;
template class DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex20,
    DRT::Element::DiscretizationType::hex20>;
template class DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_tet4,
    DRT::Element::DiscretizationType::tet4>;
template class DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_tet10,
    DRT::Element::DiscretizationType::tet10>;
template class DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::NURBS::So_nurbs27,
    DRT::Element::DiscretizationType::nurbs27>;

#endif
