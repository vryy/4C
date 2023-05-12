/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluation of scalar transport elements for standard scalar transport problems

\level 1

*/
/*--------------------------------------------------------------------------*/
#include "scatra_ele_calc_std.H"
#include "utils_singleton_owner.H"


/*----------------------------------------------------------------------*
 | singleton access method                                   fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
DRT::ELEMENTS::ScaTraEleCalcStd<distype, probdim>*
DRT::ELEMENTS::ScaTraEleCalcStd<distype, probdim>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = ::UTILS::MakeSingletonMap<std::pair<std::string, int>>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleCalcStd<distype, probdim>>(
            new ScaTraEleCalcStd<distype, probdim>(numdofpernode, numscal, disname));
      });

  return singleton_map[std::make_pair(disname, numdofpernode)].Instance(
      ::UTILS::SingletonAction::create, numdofpernode, numscal, disname);
}


/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
DRT::ELEMENTS::ScaTraEleCalcStd<distype, probdim>::ScaTraEleCalcStd(
    const int numdofpernode, const int numscal, const std::string& disname)
    : DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::ScaTraEleCalc(numdofpernode, numscal, disname)
{
  return;
}


// template classes

// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::line2, 1>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::line2, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::line2, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::line3, 1>;
// template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::line3,2>;
// template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::line3,3>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::tri3, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::tri3, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::tri6, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::quad4, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::quad4, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::quad9, 2>;
// template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::quad9,3>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::nurbs9, 2>;
// template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::nurbs9,3>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::hex8, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::hex27, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::tet4, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::tet10, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::pyramid5, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::nurbs27>;
