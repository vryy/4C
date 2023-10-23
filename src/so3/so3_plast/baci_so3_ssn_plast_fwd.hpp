/*----------------------------------------------------------------------*/
/*! \file
\brief template combinations
\level 2
*/

#ifndef BACI_SO3_SSN_PLAST_FWD_HPP
#define BACI_SO3_SSN_PLAST_FWD_HPP

// template classes
template class DRT::ELEMENTS::So3_Plast<DRT::Element::DiscretizationType::tet4>;
template class DRT::ELEMENTS::So3_Plast<DRT::Element::DiscretizationType::hex8>;
template class DRT::ELEMENTS::So3_Plast<DRT::Element::DiscretizationType::hex18>;
template class DRT::ELEMENTS::So3_Plast<DRT::Element::DiscretizationType::hex27>;
template class DRT::ELEMENTS::So3_Plast<DRT::Element::DiscretizationType::nurbs27>;

#endif