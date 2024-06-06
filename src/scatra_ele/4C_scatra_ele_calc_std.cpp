/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluation of scalar transport elements for standard scalar transport problems

\level 1

*/
/*--------------------------------------------------------------------------*/
#include "4C_scatra_ele_calc_std.hpp"

#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 | singleton access method                                   fang 02/15 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
Discret::ELEMENTS::ScaTraEleCalcStd<distype, probdim>*
Discret::ELEMENTS::ScaTraEleCalcStd<distype, probdim>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = Core::UTILS::MakeSingletonMap<std::pair<std::string, int>>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleCalcStd<distype, probdim>>(
            new ScaTraEleCalcStd<distype, probdim>(numdofpernode, numscal, disname));
      });

  return singleton_map[std::make_pair(disname, numdofpernode)].Instance(
      Core::UTILS::SingletonAction::create, numdofpernode, numscal, disname);
}


/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 02/15 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
Discret::ELEMENTS::ScaTraEleCalcStd<distype, probdim>::ScaTraEleCalcStd(
    const int numdofpernode, const int numscal, const std::string& disname)
    : Discret::ELEMENTS::ScaTraEleCalc<distype, probdim>::ScaTraEleCalc(
          numdofpernode, numscal, disname)
{
  return;
}


// template classes

// 1D elements
template class Discret::ELEMENTS::ScaTraEleCalcStd<Core::FE::CellType::line2, 1>;
template class Discret::ELEMENTS::ScaTraEleCalcStd<Core::FE::CellType::line2, 2>;
template class Discret::ELEMENTS::ScaTraEleCalcStd<Core::FE::CellType::line2, 3>;
template class Discret::ELEMENTS::ScaTraEleCalcStd<Core::FE::CellType::line3, 1>;
// template class Discret::ELEMENTS::ScaTraEleCalcStd<Core::FE::CellType::line3,2>;
// template class Discret::ELEMENTS::ScaTraEleCalcStd<Core::FE::CellType::line3,3>;

// 2D elements
template class Discret::ELEMENTS::ScaTraEleCalcStd<Core::FE::CellType::tri3, 2>;
template class Discret::ELEMENTS::ScaTraEleCalcStd<Core::FE::CellType::tri3, 3>;
template class Discret::ELEMENTS::ScaTraEleCalcStd<Core::FE::CellType::tri6, 2>;
template class Discret::ELEMENTS::ScaTraEleCalcStd<Core::FE::CellType::quad4, 2>;
template class Discret::ELEMENTS::ScaTraEleCalcStd<Core::FE::CellType::quad4, 3>;
// template class Discret::ELEMENTS::ScaTraEleCalcStd<Core::FE::CellType::quad8>;
template class Discret::ELEMENTS::ScaTraEleCalcStd<Core::FE::CellType::quad9, 2>;
// template class Discret::ELEMENTS::ScaTraEleCalcStd<Core::FE::CellType::quad9,3>;
template class Discret::ELEMENTS::ScaTraEleCalcStd<Core::FE::CellType::nurbs9, 2>;
// template class Discret::ELEMENTS::ScaTraEleCalcStd<Core::FE::CellType::nurbs9,3>;

// 3D elements
template class Discret::ELEMENTS::ScaTraEleCalcStd<Core::FE::CellType::hex8, 3>;
// template class Discret::ELEMENTS::ScaTraEleCalcStd<Core::FE::CellType::hex20>;
template class Discret::ELEMENTS::ScaTraEleCalcStd<Core::FE::CellType::hex27, 3>;
template class Discret::ELEMENTS::ScaTraEleCalcStd<Core::FE::CellType::tet4, 3>;
template class Discret::ELEMENTS::ScaTraEleCalcStd<Core::FE::CellType::tet10, 3>;
// template class Discret::ELEMENTS::ScaTraEleCalcStd<Core::FE::CellType::wedge6>;
template class Discret::ELEMENTS::ScaTraEleCalcStd<Core::FE::CellType::pyramid5, 3>;
// template class Discret::ELEMENTS::ScaTraEleCalcStd<Core::FE::CellType::nurbs27>;

FOUR_C_NAMESPACE_CLOSE
