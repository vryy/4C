#ifndef FOUR_C_SO3_PORO_P1_FWD_HPP
#define FOUR_C_SO3_PORO_P1_FWD_HPP

/*----------------------------------------------------------------------*/
/*! \file

 \brief forward declarations (explicit instantiation) of the templated p1 (mixed )poro-elements

 \level 2

 *----------------------------------------------------------------------*/

FOUR_C_NAMESPACE_OPEN

template class Discret::ELEMENTS::So3PoroP1<Discret::ELEMENTS::SoHex8, Core::FE::CellType::hex8>;
template class Discret::ELEMENTS::So3PoroP1<Discret::ELEMENTS::SoTet4, Core::FE::CellType::tet4>;

FOUR_C_NAMESPACE_CLOSE

#endif
