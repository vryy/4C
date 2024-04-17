/*----------------------------------------------------------------------*/
/*! \file
\brief template combinations
\level 2
*/

#ifndef FOUR_C_SO3_SSN_PLAST_FWD_HPP
#define FOUR_C_SO3_SSN_PLAST_FWD_HPP

FOUR_C_NAMESPACE_OPEN

// template classes
template class DRT::ELEMENTS::So3_Plast<CORE::FE::CellType::tet4>;
template class DRT::ELEMENTS::So3_Plast<CORE::FE::CellType::hex8>;
template class DRT::ELEMENTS::So3_Plast<CORE::FE::CellType::hex18>;
template class DRT::ELEMENTS::So3_Plast<CORE::FE::CellType::hex27>;
template class DRT::ELEMENTS::So3_Plast<CORE::FE::CellType::nurbs27>;

FOUR_C_NAMESPACE_CLOSE

#endif