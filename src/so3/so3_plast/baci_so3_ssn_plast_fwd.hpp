/*----------------------------------------------------------------------*/
/*! \file
\brief template combinations
\level 2
*/

#ifndef BACI_SO3_SSN_PLAST_FWD_HPP
#define BACI_SO3_SSN_PLAST_FWD_HPP

// template classes
template class DRT::ELEMENTS::So3_Plast<CORE::FE::CellType::tet4>;
template class DRT::ELEMENTS::So3_Plast<CORE::FE::CellType::hex8>;
template class DRT::ELEMENTS::So3_Plast<CORE::FE::CellType::hex18>;
template class DRT::ELEMENTS::So3_Plast<CORE::FE::CellType::hex27>;
template class DRT::ELEMENTS::So3_Plast<CORE::FE::CellType::nurbs27>;

#endif