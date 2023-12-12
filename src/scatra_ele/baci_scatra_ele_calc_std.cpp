/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluation of scalar transport elements for standard scalar transport problems

\level 1

*/
/*--------------------------------------------------------------------------*/
#include "baci_scatra_ele_calc_std.H"

#include "baci_utils_singleton_owner.H"

BACI_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 | singleton access method                                   fang 02/15 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
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
template <CORE::FE::CellType distype, int probdim>
DRT::ELEMENTS::ScaTraEleCalcStd<distype, probdim>::ScaTraEleCalcStd(
    const int numdofpernode, const int numscal, const std::string& disname)
    : DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::ScaTraEleCalc(numdofpernode, numscal, disname)
{
  return;
}


// template classes

// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcStd<CORE::FE::CellType::line2, 1>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<CORE::FE::CellType::line2, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<CORE::FE::CellType::line2, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<CORE::FE::CellType::line3, 1>;
// template class DRT::ELEMENTS::ScaTraEleCalcStd<CORE::FE::CellType::line3,2>;
// template class DRT::ELEMENTS::ScaTraEleCalcStd<CORE::FE::CellType::line3,3>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcStd<CORE::FE::CellType::tri3, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<CORE::FE::CellType::tri3, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<CORE::FE::CellType::tri6, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<CORE::FE::CellType::quad4, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<CORE::FE::CellType::quad4, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcStd<CORE::FE::CellType::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<CORE::FE::CellType::quad9, 2>;
// template class DRT::ELEMENTS::ScaTraEleCalcStd<CORE::FE::CellType::quad9,3>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<CORE::FE::CellType::nurbs9, 2>;
// template class DRT::ELEMENTS::ScaTraEleCalcStd<CORE::FE::CellType::nurbs9,3>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcStd<CORE::FE::CellType::hex8, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcStd<CORE::FE::CellType::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<CORE::FE::CellType::hex27, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<CORE::FE::CellType::tet4, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<CORE::FE::CellType::tet10, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcStd<CORE::FE::CellType::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcStd<CORE::FE::CellType::pyramid5, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcStd<CORE::FE::CellType::nurbs27>;

BACI_NAMESPACE_CLOSE
