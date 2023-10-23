#ifndef BACI_SO3_PORO_FWD_HPP
#define BACI_SO3_PORO_FWD_HPP

/*----------------------------------------------------------------------*/
/*! \file
\brief forward declarations (explicit instantiation) of the templated poro elements

\level 2
*----------------------------------------------------------------------*/
// template classes
template class DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_hex8,
    DRT::Element::DiscretizationType::hex8>;
template class DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_tet4,
    DRT::Element::DiscretizationType::tet4>;
template class DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_hex27,
    DRT::Element::DiscretizationType::hex27>;
template class DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_tet10,
    DRT::Element::DiscretizationType::tet10>;
template class DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::NURBS::So_nurbs27,
    DRT::Element::DiscretizationType::nurbs27>;

#endif
