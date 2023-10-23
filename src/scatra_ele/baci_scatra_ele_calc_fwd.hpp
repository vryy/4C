#ifndef BACI_SCATRA_ELE_CALC_FWD_HPP
#define BACI_SCATRA_ELE_CALC_FWD_HPP

/*----------------------------------------------------------------------*/
/*! \file

 \brief forward declarations for scatra_ele_calc classes

\level 1

 *----------------------------------------------------------------------*/


// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::DiscretizationType::line2, 1>;
template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::DiscretizationType::line2, 2>;
template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::DiscretizationType::line2, 3>;
template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::DiscretizationType::line3, 1>;
// template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::DiscretizationType::line3,2>;
// template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::DiscretizationType::line3,3>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::DiscretizationType::tri3, 2>;
template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::DiscretizationType::tri3, 3>;
template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::DiscretizationType::tri6, 2>;
template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::DiscretizationType::quad4, 2>;
template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::DiscretizationType::quad4, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::DiscretizationType::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::DiscretizationType::quad9, 2>;
// template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::DiscretizationType::quad9,3>;
template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::DiscretizationType::nurbs9, 2>;
// template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::DiscretizationType::nurbs9,3>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::DiscretizationType::hex8, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::DiscretizationType::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::DiscretizationType::hex27, 3>;
template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::DiscretizationType::tet4, 3>;
template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::DiscretizationType::tet10, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::DiscretizationType::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::DiscretizationType::pyramid5, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::DiscretizationType::nurbs27>;

#endif
