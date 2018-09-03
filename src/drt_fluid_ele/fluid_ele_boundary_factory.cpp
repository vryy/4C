/*!----------------------------------------------------------------------
\file fluid_ele_boundary_factory.cpp

\brief factory for fluid boundary evaluation

<pre>
\maintainer  Ursula Rasthofer & Volker Gravemeier
             {rasthofer,vgravem}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15236/-245
</pre>

\level 1

*----------------------------------------------------------------------*/

#include "fluid_ele_boundary_factory.H"
#include "fluid_ele_boundary_interface.H"
#include "fluid_ele_boundary_calc_std.H"
#include "fluid_ele_boundary_calc_poro.H"
#include "fluid_ele_calc.H"

#include "../drt_lib/drt_globalproblem.H"

/*--------------------------------------------------------------------------*
 |                                                 (public) rasthofer 11/13 |
 *--------------------------------------------------------------------------*/
DRT::ELEMENTS::FluidBoundaryInterface* DRT::ELEMENTS::FluidBoundaryFactory::ProvideImpl(
    DRT::Element::DiscretizationType distype, std::string problem)
{
  switch (distype)
  {
    case DRT::Element::quad4:
    {
      return DefineProblemType<DRT::Element::quad4>(problem);
    }
    case DRT::Element::quad8:
    {
      return DefineProblemType<DRT::Element::quad8>(problem);
    }
    case DRT::Element::quad9:
    {
      return DefineProblemType<DRT::Element::quad9>(problem);
    }
    case DRT::Element::tri3:
    {
      return DefineProblemType<DRT::Element::tri3>(problem);
    }
    case DRT::Element::tri6:
    {
      return DefineProblemType<DRT::Element::tri6>(problem);
    }
    case DRT::Element::line2:
    {
      return DefineProblemType<DRT::Element::line2>(problem);
    }
    case DRT::Element::line3:
    {
      return DefineProblemType<DRT::Element::line3>(problem);
    }
    case DRT::Element::nurbs2:
    {
      return DefineProblemType<DRT::Element::nurbs2>(problem);
    }
    case DRT::Element::nurbs3:
    {
      return DefineProblemType<DRT::Element::nurbs3>(problem);
    }
    case DRT::Element::nurbs4:
    {
      return DefineProblemType<DRT::Element::nurbs4>(problem);
    }
    case DRT::Element::nurbs9:
    {
      return DefineProblemType<DRT::Element::nurbs9>(problem);
    }
    default:
      dserror("Element shape %s not activated. Just do it.", DRT::DistypeToString(distype).c_str());
      break;
  }
  return NULL;
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

  return NULL;
}
