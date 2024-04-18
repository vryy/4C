#ifndef FOUR_C_SO3_PORO_FWD_HPP
#define FOUR_C_SO3_PORO_FWD_HPP

/*----------------------------------------------------------------------*/
/*! \file
\brief forward declarations (explicit instantiation) of the templated poro elements

\level 2
*----------------------------------------------------------------------*/

FOUR_C_NAMESPACE_OPEN

// template classes
template class DRT::ELEMENTS::So3Poro<DRT::ELEMENTS::SoHex8, CORE::FE::CellType::hex8>;
template class DRT::ELEMENTS::So3Poro<DRT::ELEMENTS::SoTet4, CORE::FE::CellType::tet4>;
template class DRT::ELEMENTS::So3Poro<DRT::ELEMENTS::SoHex27, CORE::FE::CellType::hex27>;
template class DRT::ELEMENTS::So3Poro<DRT::ELEMENTS::SoTet10, CORE::FE::CellType::tet10>;
template class DRT::ELEMENTS::So3Poro<DRT::ELEMENTS::NURBS::SoNurbs27, CORE::FE::CellType::nurbs27>;

FOUR_C_NAMESPACE_CLOSE

#endif
