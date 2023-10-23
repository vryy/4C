/*----------------------------------------------------------------------*/
/*! \file

\brief factory class into templated evaluators for fluid boundary integration

\level 1


*/
/*----------------------------------------------------------------------*/

#include "baci_fluid_ele_boundary_factory.H"

#include "baci_fluid_ele_boundary_calc_poro.H"
#include "baci_fluid_ele_boundary_calc_std.H"
#include "baci_fluid_ele_boundary_interface.H"
#include "baci_fluid_ele_calc.H"

/*--------------------------------------------------------------------------*
 |                                                 (public) rasthofer 11/13 |
 *--------------------------------------------------------------------------*/
DRT::ELEMENTS::FluidBoundaryInterface* DRT::ELEMENTS::FluidBoundaryFactory::ProvideImpl(
    DRT::Element::DiscretizationType distype, std::string problem)
{
  switch (distype)
  {
    case DRT::Element::DiscretizationType::quad4:
    {
      return DefineProblemType<DRT::Element::DiscretizationType::quad4>(problem);
    }
    case DRT::Element::DiscretizationType::quad8:
    {
      return DefineProblemType<DRT::Element::DiscretizationType::quad8>(problem);
    }
    case DRT::Element::DiscretizationType::quad9:
    {
      return DefineProblemType<DRT::Element::DiscretizationType::quad9>(problem);
    }
    case DRT::Element::DiscretizationType::tri3:
    {
      return DefineProblemType<DRT::Element::DiscretizationType::tri3>(problem);
    }
    case DRT::Element::DiscretizationType::tri6:
    {
      return DefineProblemType<DRT::Element::DiscretizationType::tri6>(problem);
    }
    case DRT::Element::DiscretizationType::line2:
    {
      return DefineProblemType<DRT::Element::DiscretizationType::line2>(problem);
    }
    case DRT::Element::DiscretizationType::line3:
    {
      return DefineProblemType<DRT::Element::DiscretizationType::line3>(problem);
    }
    case DRT::Element::DiscretizationType::nurbs2:
    {
      return DefineProblemType<DRT::Element::DiscretizationType::nurbs2>(problem);
    }
    case DRT::Element::DiscretizationType::nurbs3:
    {
      return DefineProblemType<DRT::Element::DiscretizationType::nurbs3>(problem);
    }
    case DRT::Element::DiscretizationType::nurbs4:
    {
      return DefineProblemType<DRT::Element::DiscretizationType::nurbs4>(problem);
    }
    case DRT::Element::DiscretizationType::nurbs9:
    {
      return DefineProblemType<DRT::Element::DiscretizationType::nurbs9>(problem);
    }
    default:
      dserror("Element shape %s not activated. Just do it.", DRT::DistypeToString(distype).c_str());
      break;
  }
  return nullptr;
}

/*--------------------------------------------------------------------------*
 |                                                 (public) rasthofer 11/13 |
 *--------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidBoundaryInterface* DRT::ELEMENTS::FluidBoundaryFactory::DefineProblemType(
    std::string problem)
{
  if (problem == "std")
    return DRT::ELEMENTS::FluidEleBoundaryCalcStd<distype>::Instance();
  else if (problem == "poro")
    return DRT::ELEMENTS::FluidEleBoundaryCalcPoro<distype>::Instance();
  else if (problem == "poro_p1")
    return DRT::ELEMENTS::FluidEleBoundaryCalcPoroP1<distype>::Instance();
  else
    dserror("Defined problem type does not exist!!");

  return nullptr;
}
