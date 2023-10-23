/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluation of scalar transport elements for standard scalar transport problems

\level 1

*/
/*--------------------------------------------------------------------------*/
#include "baci_scatra_ele_calc_std.H"

#include "baci_utils_singleton_owner.H"


/*----------------------------------------------------------------------*
 | singleton access method                                   fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
DRT::ELEMENTS::ScaTraEleCalcStd<distype, probdim>*
DRT::ELEMENTS::ScaTraEleCalcStd<distype, probdim>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = CORE::UTILS::MakeSingletonMap<std::pair<std::string, int>>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleCalcStd<distype, probdim>>(
            new ScaTraEleCalcStd<distype, probdim>(numdofpernode, numscal, disname));
      });

  return singleton_map[std::make_pair(disname, numdofpernode)].Instance(
      CORE::UTILS::SingletonAction::create, numdofpernode, numscal, disname);
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
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::DiscretizationType::line2, 1>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::DiscretizationType::line2, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::DiscretizationType::line2, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::DiscretizationType::line3, 1>;
// template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::DiscretizationType::line3,2>;
// template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::DiscretizationType::line3,3>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::DiscretizationType::tri3, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::DiscretizationType::tri3, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::DiscretizationType::tri6, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::DiscretizationType::quad4, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::DiscretizationType::quad4, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::DiscretizationType::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::DiscretizationType::quad9, 2>;
// template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::DiscretizationType::quad9,3>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::DiscretizationType::nurbs9, 2>;
// template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::DiscretizationType::nurbs9,3>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::DiscretizationType::hex8, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::DiscretizationType::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::DiscretizationType::hex27, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::DiscretizationType::tet4, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::DiscretizationType::tet10, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::DiscretizationType::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::DiscretizationType::pyramid5, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::DiscretizationType::nurbs27>;
