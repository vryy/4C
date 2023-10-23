/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of fluid terms at integration points of boundaries

\level 1


*/
/*----------------------------------------------------------------------*/

#include "baci_fluid_ele_boundary_calc_std.H"

#include "baci_fluid_ele_parameter_std.H"
#include "baci_lib_elementtype.H"

template <DRT::Element::DiscretizationType distype>
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
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidEleBoundaryCalcStd<distype>::FluidEleBoundaryCalcStd()
    : DRT::ELEMENTS::FluidBoundaryImpl<distype>::FluidBoundaryImpl()
{
  // pointer to class FluidImplParameter
  my::fldpara_ = DRT::ELEMENTS::FluidEleParameterStd::Instance();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
// template classes
template class DRT::ELEMENTS::FluidEleBoundaryCalcStd<DRT::Element::DiscretizationType::quad4>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcStd<DRT::Element::DiscretizationType::quad8>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcStd<DRT::Element::DiscretizationType::quad9>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcStd<DRT::Element::DiscretizationType::tri3>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcStd<DRT::Element::DiscretizationType::tri6>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcStd<DRT::Element::DiscretizationType::line2>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcStd<DRT::Element::DiscretizationType::line3>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcStd<DRT::Element::DiscretizationType::nurbs2>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcStd<DRT::Element::DiscretizationType::nurbs3>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcStd<DRT::Element::DiscretizationType::nurbs4>;
template class DRT::ELEMENTS::FluidEleBoundaryCalcStd<DRT::Element::DiscretizationType::nurbs9>;
