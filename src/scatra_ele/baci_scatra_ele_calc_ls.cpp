/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluations for level sets

\level 2

*/
/*--------------------------------------------------------------------------*/
#include "baci_scatra_ele_calc_ls.hpp"

#include "baci_scatra_ele_parameter_std.hpp"
#include "baci_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
DRT::ELEMENTS::ScaTraEleCalcLS<distype>* DRT::ELEMENTS::ScaTraEleCalcLS<distype>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = CORE::UTILS::MakeSingletonMap<std::string>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleCalcLS<distype>>(
            new ScaTraEleCalcLS<distype>(numdofpernode, numscal, disname));
      });

  return singleton_map[disname].Instance(
      CORE::UTILS::SingletonAction::create, numdofpernode, numscal, disname);
}


/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 02/15 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
DRT::ELEMENTS::ScaTraEleCalcLS<distype>::ScaTraEleCalcLS(
    const int numdofpernode, const int numscal, const std::string& disname)
    : DRT::ELEMENTS::ScaTraEleCalc<distype>::ScaTraEleCalc(numdofpernode, numscal, disname)
{
  // safety check
  if (my::scatrapara_->RBSubGrVel())
    FOUR_C_THROW("CalcSubgrVelocityLevelSet not available anymore");

  return;
}


// template classes

// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcLS<CORE::FE::CellType::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<CORE::FE::CellType::line3>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcLS<CORE::FE::CellType::tri3>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<CORE::FE::CellType::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<CORE::FE::CellType::quad4>;
// template class DRT::ELEMENTS::ScaTraEleCalcLS<CORE::FE::CellType::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<CORE::FE::CellType::quad9>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<CORE::FE::CellType::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcLS<CORE::FE::CellType::hex8>;
// template class DRT::ELEMENTS::ScaTraEleCalcLS<CORE::FE::CellType::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<CORE::FE::CellType::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<CORE::FE::CellType::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<CORE::FE::CellType::tet10>;
// template class DRT::ELEMENTS::ScaTraEleCalcLS<CORE::FE::CellType::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<CORE::FE::CellType::pyramid5>;
// template class DRT::ELEMENTS::ScaTraEleCalcLS<CORE::FE::CellType::nurbs27>;

FOUR_C_NAMESPACE_CLOSE
