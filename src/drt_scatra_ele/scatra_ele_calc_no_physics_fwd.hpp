/*----------------------------------------------------------------------*/
/*!

\brief forward declarations for scatra_ele_calc_no_physics classes

\level 2

\maintainer Amadeus Gebauer

*/

// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::line2, 1>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::line2, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::line2, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::line3, 1>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::tri3, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::tri3, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::tri6, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::quad4, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::quad4, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::quad9, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::nurbs9, 2>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::hex8, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::hex27, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::tet4, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::tet10, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::pyramid5, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::nurbs27>;