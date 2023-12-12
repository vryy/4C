#ifndef BACI_SO3_THERMO_FWD_HPP
#define BACI_SO3_THERMO_FWD_HPP

/*----------------------------------------------------------------------*/
/*! \file

\level 1

\brief forward declarations of discretisation types

 *
 */

BACI_NAMESPACE_OPEN

// template classes
template class DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex8, CORE::FE::CellType::hex8>;
template class DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex8fbar, CORE::FE::CellType::hex8>;
template class DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex27, CORE::FE::CellType::hex27>;
template class DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex20, CORE::FE::CellType::hex20>;
template class DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_tet4, CORE::FE::CellType::tet4>;
template class DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_tet10, CORE::FE::CellType::tet10>;
template class DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::NURBS::So_nurbs27,
    CORE::FE::CellType::nurbs27>;

BACI_NAMESPACE_CLOSE

#endif
