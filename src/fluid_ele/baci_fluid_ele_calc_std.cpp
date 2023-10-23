/*----------------------------------------------------------------------*/
/*! \file

\brief standard routines for calculation of fluid element

\level 1


*/
/*----------------------------------------------------------------------*/

#include "baci_fluid_ele_calc_std.H"

#include "baci_fluid_ele_parameter_std.H"


template <DRT::Element::DiscretizationType distype>
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
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidEleCalcStd<distype>::FluidEleCalcStd()
    : DRT::ELEMENTS::FluidEleCalc<distype>::FluidEleCalc()
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/

// template classes
template class DRT::ELEMENTS::FluidEleCalcStd<DRT::Element::DiscretizationType::hex8>;
template class DRT::ELEMENTS::FluidEleCalcStd<DRT::Element::DiscretizationType::hex20>;
template class DRT::ELEMENTS::FluidEleCalcStd<DRT::Element::DiscretizationType::hex27>;
template class DRT::ELEMENTS::FluidEleCalcStd<DRT::Element::DiscretizationType::tet4>;
template class DRT::ELEMENTS::FluidEleCalcStd<DRT::Element::DiscretizationType::tet10>;
template class DRT::ELEMENTS::FluidEleCalcStd<DRT::Element::DiscretizationType::wedge6>;
template class DRT::ELEMENTS::FluidEleCalcStd<DRT::Element::DiscretizationType::wedge15>;
template class DRT::ELEMENTS::FluidEleCalcStd<DRT::Element::DiscretizationType::pyramid5>;
template class DRT::ELEMENTS::FluidEleCalcStd<DRT::Element::DiscretizationType::quad4>;
template class DRT::ELEMENTS::FluidEleCalcStd<DRT::Element::DiscretizationType::quad8>;
template class DRT::ELEMENTS::FluidEleCalcStd<DRT::Element::DiscretizationType::quad9>;
template class DRT::ELEMENTS::FluidEleCalcStd<DRT::Element::DiscretizationType::tri3>;
template class DRT::ELEMENTS::FluidEleCalcStd<DRT::Element::DiscretizationType::tri6>;
template class DRT::ELEMENTS::FluidEleCalcStd<DRT::Element::DiscretizationType::nurbs9>;
template class DRT::ELEMENTS::FluidEleCalcStd<DRT::Element::DiscretizationType::nurbs27>;
