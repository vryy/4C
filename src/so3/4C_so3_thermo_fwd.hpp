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
template class Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::SoHex8, Core::FE::CellType::hex8>;
template class Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::SoHex8fbar,
    Core::FE::CellType::hex8>;
template class Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::SoHex27, Core::FE::CellType::hex27>;
template class Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::SoHex20, Core::FE::CellType::hex20>;
template class Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::SoTet4, Core::FE::CellType::tet4>;
template class Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::SoTet10, Core::FE::CellType::tet10>;
template class Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::Nurbs::SoNurbs27,
    Core::FE::CellType::nurbs27>;

FOUR_C_NAMESPACE_CLOSE

#endif
