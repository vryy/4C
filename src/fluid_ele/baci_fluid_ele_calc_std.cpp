/*----------------------------------------------------------------------*/
/*! \file

\brief standard routines for calculation of fluid element

\level 1


*/
/*----------------------------------------------------------------------*/

#include "baci_fluid_ele_calc_std.hpp"

#include "baci_fluid_ele_parameter_std.hpp"

FOUR_C_NAMESPACE_OPEN


template <CORE::FE::CellType distype>
DRT::ELEMENTS::FluidEleCalcStd<distype>* DRT::ELEMENTS::FluidEleCalcStd<distype>::Instance(
    CORE::UTILS::SingletonAction action)
{
  static auto singleton_owner = CORE::UTILS::MakeSingletonOwner(
      []()
      {
        return std::unique_ptr<DRT::ELEMENTS::FluidEleCalcStd<distype>>(
            new DRT::ELEMENTS::FluidEleCalcStd<distype>());
      });

  return singleton_owner.Instance(action);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
DRT::ELEMENTS::FluidEleCalcStd<distype>::FluidEleCalcStd()
    : DRT::ELEMENTS::FluidEleCalc<distype>::FluidEleCalc()
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/

// template classes
template class DRT::ELEMENTS::FluidEleCalcStd<CORE::FE::CellType::hex8>;
template class DRT::ELEMENTS::FluidEleCalcStd<CORE::FE::CellType::hex20>;
template class DRT::ELEMENTS::FluidEleCalcStd<CORE::FE::CellType::hex27>;
template class DRT::ELEMENTS::FluidEleCalcStd<CORE::FE::CellType::tet4>;
template class DRT::ELEMENTS::FluidEleCalcStd<CORE::FE::CellType::tet10>;
template class DRT::ELEMENTS::FluidEleCalcStd<CORE::FE::CellType::wedge6>;
template class DRT::ELEMENTS::FluidEleCalcStd<CORE::FE::CellType::wedge15>;
template class DRT::ELEMENTS::FluidEleCalcStd<CORE::FE::CellType::pyramid5>;
template class DRT::ELEMENTS::FluidEleCalcStd<CORE::FE::CellType::quad4>;
template class DRT::ELEMENTS::FluidEleCalcStd<CORE::FE::CellType::quad8>;
template class DRT::ELEMENTS::FluidEleCalcStd<CORE::FE::CellType::quad9>;
template class DRT::ELEMENTS::FluidEleCalcStd<CORE::FE::CellType::tri3>;
template class DRT::ELEMENTS::FluidEleCalcStd<CORE::FE::CellType::tri6>;
template class DRT::ELEMENTS::FluidEleCalcStd<CORE::FE::CellType::nurbs9>;
template class DRT::ELEMENTS::FluidEleCalcStd<CORE::FE::CellType::nurbs27>;

FOUR_C_NAMESPACE_CLOSE
