#ifndef FOUR_C_SO3_PORO_FWD_HPP
#define FOUR_C_SO3_PORO_FWD_HPP

/*----------------------------------------------------------------------*/
/*! \file
\brief forward declarations (explicit instantiation) of the templated poro elements

\level 2
*----------------------------------------------------------------------*/

FOUR_C_NAMESPACE_OPEN

// template classes
template class DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_hex8, CORE::FE::CellType::hex8>;
template class DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_tet4, CORE::FE::CellType::tet4>;
template class DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_hex27, CORE::FE::CellType::hex27>;
template class DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_tet10, CORE::FE::CellType::tet10>;
template class DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::NURBS::So_nurbs27,
    CORE::FE::CellType::nurbs27>;

FOUR_C_NAMESPACE_CLOSE

#endif
