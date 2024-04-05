#ifndef FOUR_C_SCATRA_ELE_CALC_FWD_HPP
#define FOUR_C_SCATRA_ELE_CALC_FWD_HPP

/*----------------------------------------------------------------------*/
/*! \file

 \brief forward declarations for scatra_ele_calc classes

\level 1

 *----------------------------------------------------------------------*/

BACI_NAMESPACE_OPEN

// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalc<CORE::FE::CellType::line2, 1>;
template class DRT::ELEMENTS::ScaTraEleCalc<CORE::FE::CellType::line2, 2>;
template class DRT::ELEMENTS::ScaTraEleCalc<CORE::FE::CellType::line2, 3>;
template class DRT::ELEMENTS::ScaTraEleCalc<CORE::FE::CellType::line3, 1>;
// template class DRT::ELEMENTS::ScaTraEleCalc<CORE::FE::CellType::line3,2>;
// template class DRT::ELEMENTS::ScaTraEleCalc<CORE::FE::CellType::line3,3>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalc<CORE::FE::CellType::tri3, 2>;
template class DRT::ELEMENTS::ScaTraEleCalc<CORE::FE::CellType::tri3, 3>;
template class DRT::ELEMENTS::ScaTraEleCalc<CORE::FE::CellType::tri6, 2>;
template class DRT::ELEMENTS::ScaTraEleCalc<CORE::FE::CellType::quad4, 2>;
template class DRT::ELEMENTS::ScaTraEleCalc<CORE::FE::CellType::quad4, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalc<CORE::FE::CellType::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalc<CORE::FE::CellType::quad9, 2>;
// template class DRT::ELEMENTS::ScaTraEleCalc<CORE::FE::CellType::quad9,3>;
template class DRT::ELEMENTS::ScaTraEleCalc<CORE::FE::CellType::nurbs9, 2>;
// template class DRT::ELEMENTS::ScaTraEleCalc<CORE::FE::CellType::nurbs9,3>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalc<CORE::FE::CellType::hex8, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalc<CORE::FE::CellType::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalc<CORE::FE::CellType::hex27, 3>;
template class DRT::ELEMENTS::ScaTraEleCalc<CORE::FE::CellType::tet4, 3>;
template class DRT::ELEMENTS::ScaTraEleCalc<CORE::FE::CellType::tet10, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalc<CORE::FE::CellType::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalc<CORE::FE::CellType::pyramid5, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalc<CORE::FE::CellType::nurbs27>;

BACI_NAMESPACE_CLOSE

#endif
