#ifndef FOUR_C_SCATRA_ELE_CALC_FWD_HPP
#define FOUR_C_SCATRA_ELE_CALC_FWD_HPP

/*----------------------------------------------------------------------*/
/*! \file

 \brief forward declarations for scatra_ele_calc classes

\level 1

 *----------------------------------------------------------------------*/

FOUR_C_NAMESPACE_OPEN

// 1D elements
template class Discret::ELEMENTS::ScaTraEleCalc<Core::FE::CellType::line2, 1>;
template class Discret::ELEMENTS::ScaTraEleCalc<Core::FE::CellType::line2, 2>;
template class Discret::ELEMENTS::ScaTraEleCalc<Core::FE::CellType::line2, 3>;
template class Discret::ELEMENTS::ScaTraEleCalc<Core::FE::CellType::line3, 1>;
// template class Discret::ELEMENTS::ScaTraEleCalc<Core::FE::CellType::line3,2>;
// template class Discret::ELEMENTS::ScaTraEleCalc<Core::FE::CellType::line3,3>;

// 2D elements
template class Discret::ELEMENTS::ScaTraEleCalc<Core::FE::CellType::tri3, 2>;
template class Discret::ELEMENTS::ScaTraEleCalc<Core::FE::CellType::tri3, 3>;
template class Discret::ELEMENTS::ScaTraEleCalc<Core::FE::CellType::tri6, 2>;
template class Discret::ELEMENTS::ScaTraEleCalc<Core::FE::CellType::quad4, 2>;
template class Discret::ELEMENTS::ScaTraEleCalc<Core::FE::CellType::quad4, 3>;
// template class Discret::ELEMENTS::ScaTraEleCalc<Core::FE::CellType::quad8>;
template class Discret::ELEMENTS::ScaTraEleCalc<Core::FE::CellType::quad9, 2>;
// template class Discret::ELEMENTS::ScaTraEleCalc<Core::FE::CellType::quad9,3>;
template class Discret::ELEMENTS::ScaTraEleCalc<Core::FE::CellType::nurbs9, 2>;
// template class Discret::ELEMENTS::ScaTraEleCalc<Core::FE::CellType::nurbs9,3>;

// 3D elements
template class Discret::ELEMENTS::ScaTraEleCalc<Core::FE::CellType::hex8, 3>;
// template class Discret::ELEMENTS::ScaTraEleCalc<Core::FE::CellType::hex20>;
template class Discret::ELEMENTS::ScaTraEleCalc<Core::FE::CellType::hex27, 3>;
template class Discret::ELEMENTS::ScaTraEleCalc<Core::FE::CellType::tet4, 3>;
template class Discret::ELEMENTS::ScaTraEleCalc<Core::FE::CellType::tet10, 3>;
// template class Discret::ELEMENTS::ScaTraEleCalc<Core::FE::CellType::wedge6>;
template class Discret::ELEMENTS::ScaTraEleCalc<Core::FE::CellType::pyramid5, 3>;
// template class Discret::ELEMENTS::ScaTraEleCalc<Core::FE::CellType::nurbs27>;

FOUR_C_NAMESPACE_CLOSE

#endif
