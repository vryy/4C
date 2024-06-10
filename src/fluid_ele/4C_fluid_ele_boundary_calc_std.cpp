/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of fluid terms at integration points of boundaries

\level 1


*/
/*----------------------------------------------------------------------*/

#include "4C_fluid_ele_boundary_calc_std.hpp"

#include "4C_fem_general_elementtype.hpp"
#include "4C_fluid_ele_parameter_std.hpp"

FOUR_C_NAMESPACE_OPEN

template <Core::FE::CellType distype>
Discret::ELEMENTS::FluidEleBoundaryCalcStd<distype>*
Discret::ELEMENTS::FluidEleBoundaryCalcStd<distype>::Instance(Core::UTILS::SingletonAction action)
{
  static auto singleton_owner = Core::UTILS::MakeSingletonOwner(
      []()
      {
        return std::unique_ptr<Discret::ELEMENTS::FluidEleBoundaryCalcStd<distype>>(
            new Discret::ELEMENTS::FluidEleBoundaryCalcStd<distype>());
      });

  return singleton_owner.Instance(action);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::ELEMENTS::FluidEleBoundaryCalcStd<distype>::FluidEleBoundaryCalcStd()
    : Discret::ELEMENTS::FluidBoundaryImpl<distype>::FluidBoundaryImpl()
{
  // pointer to class FluidImplParameter
  my::fldpara_ = Discret::ELEMENTS::FluidEleParameterStd::Instance();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
// template classes
template class Discret::ELEMENTS::FluidEleBoundaryCalcStd<Core::FE::CellType::quad4>;
template class Discret::ELEMENTS::FluidEleBoundaryCalcStd<Core::FE::CellType::quad8>;
template class Discret::ELEMENTS::FluidEleBoundaryCalcStd<Core::FE::CellType::quad9>;
template class Discret::ELEMENTS::FluidEleBoundaryCalcStd<Core::FE::CellType::tri3>;
template class Discret::ELEMENTS::FluidEleBoundaryCalcStd<Core::FE::CellType::tri6>;
template class Discret::ELEMENTS::FluidEleBoundaryCalcStd<Core::FE::CellType::line2>;
template class Discret::ELEMENTS::FluidEleBoundaryCalcStd<Core::FE::CellType::line3>;
template class Discret::ELEMENTS::FluidEleBoundaryCalcStd<Core::FE::CellType::nurbs2>;
template class Discret::ELEMENTS::FluidEleBoundaryCalcStd<Core::FE::CellType::nurbs3>;
template class Discret::ELEMENTS::FluidEleBoundaryCalcStd<Core::FE::CellType::nurbs4>;
template class Discret::ELEMENTS::FluidEleBoundaryCalcStd<Core::FE::CellType::nurbs9>;

FOUR_C_NAMESPACE_CLOSE
