#ifndef FOUR_C_SO3_THERMO_FWD_HPP
#define FOUR_C_SO3_THERMO_FWD_HPP

/*----------------------------------------------------------------------*/
/*! \file

\level 1

\brief forward declarations of discretisation types

 *
 */

FOUR_C_NAMESPACE_OPEN

// template classes
template class DRT::ELEMENTS::So3Thermo<DRT::ELEMENTS::SoHex8, CORE::FE::CellType::hex8>;
template class DRT::ELEMENTS::So3Thermo<DRT::ELEMENTS::SoHex8fbar, CORE::FE::CellType::hex8>;
template class DRT::ELEMENTS::So3Thermo<DRT::ELEMENTS::SoHex27, CORE::FE::CellType::hex27>;
template class DRT::ELEMENTS::So3Thermo<DRT::ELEMENTS::SoHex20, CORE::FE::CellType::hex20>;
template class DRT::ELEMENTS::So3Thermo<DRT::ELEMENTS::SoTet4, CORE::FE::CellType::tet4>;
template class DRT::ELEMENTS::So3Thermo<DRT::ELEMENTS::SoTet10, CORE::FE::CellType::tet10>;
template class DRT::ELEMENTS::So3Thermo<DRT::ELEMENTS::NURBS::SoNurbs27,
    CORE::FE::CellType::nurbs27>;

FOUR_C_NAMESPACE_CLOSE

#endif
