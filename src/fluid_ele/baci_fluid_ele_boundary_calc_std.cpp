/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of fluid terms at integration points of boundaries

\level 1


*/
/*----------------------------------------------------------------------*/

#include "baci_fluid_ele_boundary_calc_std.H"

#include "baci_fluid_ele_parameter_std.H"
#include "baci_lib_elementtype.H"

template <CORE::FE::CellType distype>
DRT::ELEMENTS::FluidEleBoundaryCalcStd<distype>*
DRT::ELEMENTS::FluidEleBoundaryCalcStd<distype>::Instance(CORE::UTILS::SingletonAction action)
{
  static auto singleton_owner = CORE::UTILS::MakeSingletonOwner(
      []()
      {
        return std::unique_ptr<DRT::ELEMENTS::FluidEleBoundaryCalcStd<distype>>(
            new DRT::ELEMENTS::FluidEleBoundaryCalcStd<distype>());
      });

  return singleton_owner.Instance(action);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
DRT::ELEMENTS::FluidEleBoundaryCalcStd<distype>::FluidEleBoundaryCalcStd()
    : DRT::ELEMENTS::FluidBoundaryImpl<distype>::FluidBoundaryImpl()
{
  // pointer to class FluidImplParameter
  my::fldpara_ = DRT::ELEMENTS::FluidEleParameterStd::Instance();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
// template classes
template class DRT::ELEMENTS::FluidEleBoundaryCalcStd<CORE::FE::CellType::quad4>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcStd<CORE::FE::CellType::quad8>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcStd<CORE::FE::CellType::quad9>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcStd<CORE::FE::CellType::tri3>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcStd<CORE::FE::CellType::tri6>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcStd<CORE::FE::CellType::line2>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcStd<CORE::FE::CellType::line3>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcStd<CORE::FE::CellType::nurbs2>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcStd<CORE::FE::CellType::nurbs3>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcStd<CORE::FE::CellType::nurbs4>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcStd<CORE::FE::CellType::nurbs9>;
