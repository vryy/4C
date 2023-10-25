#ifndef BACI_SCATRA_ELE_CALC_NO_PHYSICS_FWD_HPP
#define BACI_SCATRA_ELE_CALC_NO_PHYSICS_FWD_HPP

/*----------------------------------------------------------------------*/
/*! \file

\brief forward declarations for scatra_ele_calc_no_physics classes

\level 2


*/

// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<CORE::FE::CellType::line2, 1>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<CORE::FE::CellType::line2, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<CORE::FE::CellType::line2, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<CORE::FE::CellType::line3, 1>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<CORE::FE::CellType::tri3, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<CORE::FE::CellType::tri3, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<CORE::FE::CellType::tri6, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<CORE::FE::CellType::quad4, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<CORE::FE::CellType::quad4, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<CORE::FE::CellType::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<CORE::FE::CellType::quad9, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<CORE::FE::CellType::nurbs9, 2>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<CORE::FE::CellType::hex8, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<CORE::FE::CellType::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<CORE::FE::CellType::hex27, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<CORE::FE::CellType::tet4, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<CORE::FE::CellType::tet10, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<CORE::FE::CellType::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<CORE::FE::CellType::pyramid5, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<CORE::FE::CellType::nurbs27>;
#endif
