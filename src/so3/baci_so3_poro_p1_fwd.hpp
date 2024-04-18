#ifndef FOUR_C_SO3_PORO_P1_FWD_HPP
#define FOUR_C_SO3_PORO_P1_FWD_HPP

/*----------------------------------------------------------------------*/
/*! \file

 \brief forward declarations (explicit instantiation) of the templated p1 (mixed )poro-elements

 \level 2

 *----------------------------------------------------------------------*/

FOUR_C_NAMESPACE_OPEN

template class DRT::ELEMENTS::So3PoroP1<DRT::ELEMENTS::SoHex8, CORE::FE::CellType::hex8>;
template class DRT::ELEMENTS::So3PoroP1<DRT::ELEMENTS::SoTet4, CORE::FE::CellType::tet4>;

FOUR_C_NAMESPACE_CLOSE

#endif
