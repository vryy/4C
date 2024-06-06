/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluations for level sets

\level 2

*/
/*--------------------------------------------------------------------------*/
#include "4C_scatra_ele_calc_ls.hpp"

#include "4C_scatra_ele_parameter_std.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::ELEMENTS::ScaTraEleCalcLS<distype>* Discret::ELEMENTS::ScaTraEleCalcLS<distype>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = Core::UTILS::MakeSingletonMap<std::string>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleCalcLS<distype>>(
            new ScaTraEleCalcLS<distype>(numdofpernode, numscal, disname));
      });

  return singleton_map[disname].Instance(
      Core::UTILS::SingletonAction::create, numdofpernode, numscal, disname);
}


/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 02/15 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::ELEMENTS::ScaTraEleCalcLS<distype>::ScaTraEleCalcLS(
    const int numdofpernode, const int numscal, const std::string& disname)
    : Discret::ELEMENTS::ScaTraEleCalc<distype>::ScaTraEleCalc(numdofpernode, numscal, disname)
{
  // safety check
  if (my::scatrapara_->RBSubGrVel())
    FOUR_C_THROW("CalcSubgrVelocityLevelSet not available anymore");

  return;
}


// template classes

// 1D elements
template class Discret::ELEMENTS::ScaTraEleCalcLS<Core::FE::CellType::line2>;
template class Discret::ELEMENTS::ScaTraEleCalcLS<Core::FE::CellType::line3>;

// 2D elements
template class Discret::ELEMENTS::ScaTraEleCalcLS<Core::FE::CellType::tri3>;
template class Discret::ELEMENTS::ScaTraEleCalcLS<Core::FE::CellType::tri6>;
template class Discret::ELEMENTS::ScaTraEleCalcLS<Core::FE::CellType::quad4>;
// template class Discret::ELEMENTS::ScaTraEleCalcLS<Core::FE::CellType::quad8>;
template class Discret::ELEMENTS::ScaTraEleCalcLS<Core::FE::CellType::quad9>;
template class Discret::ELEMENTS::ScaTraEleCalcLS<Core::FE::CellType::nurbs9>;

// 3D elements
template class Discret::ELEMENTS::ScaTraEleCalcLS<Core::FE::CellType::hex8>;
// template class Discret::ELEMENTS::ScaTraEleCalcLS<Core::FE::CellType::hex20>;
template class Discret::ELEMENTS::ScaTraEleCalcLS<Core::FE::CellType::hex27>;
template class Discret::ELEMENTS::ScaTraEleCalcLS<Core::FE::CellType::tet4>;
template class Discret::ELEMENTS::ScaTraEleCalcLS<Core::FE::CellType::tet10>;
// template class Discret::ELEMENTS::ScaTraEleCalcLS<Core::FE::CellType::wedge6>;
template class Discret::ELEMENTS::ScaTraEleCalcLS<Core::FE::CellType::pyramid5>;
// template class Discret::ELEMENTS::ScaTraEleCalcLS<Core::FE::CellType::nurbs27>;

FOUR_C_NAMESPACE_CLOSE
