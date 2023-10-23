/*--------------------------------------------------------------------------*/
/*! \file

\brief Factory of electromagnetic elements

\level 2

*/
/*--------------------------------------------------------------------------*/

#include "baci_elemag_ele_factory.H"

#include "baci_elemag_diff_ele_calc.H"
#include "baci_elemag_ele_calc.H"
#include "baci_elemag_ele_interface.H"

/*--------------------------------------------------------------------------*
 |                                                (public) berardocco 02/18 |
 *--------------------------------------------------------------------------*/
DRT::ELEMENTS::ElemagEleInterface* DRT::ELEMENTS::ElemagFactory::ProvideImpl(
    DRT::Element::DiscretizationType distype, std::string problem)
{
  switch (distype)
  {
    case DRT::Element::DiscretizationType::hex8:
    {
      return DefineProblemType<DRT::Element::DiscretizationType::hex8>(problem);
    }
    case DRT::Element::DiscretizationType::hex20:
    {
      return DefineProblemType<DRT::Element::DiscretizationType::hex20>(problem);
    }
    case DRT::Element::DiscretizationType::hex27:
    {
      return DefineProblemType<DRT::Element::DiscretizationType::hex27>(problem);
    }
    case DRT::Element::DiscretizationType::tet4:
    {
      return DefineProblemType<DRT::Element::DiscretizationType::tet4>(problem);
    }
    case DRT::Element::DiscretizationType::tet10:
    {
      return DefineProblemType<DRT::Element::DiscretizationType::tet10>(problem);
    }
    case DRT::Element::DiscretizationType::wedge6:
    {
      return DefineProblemType<DRT::Element::DiscretizationType::wedge6>(problem);
    }
    /* wedge15 cannot be used since no mesh generator exists
    case DRT::Element::DiscretizationType::wedge15:
    {
      return DefineProblemType<DRT::Element::DiscretizationType::wedge15>(problem);
    }
    */
    case DRT::Element::DiscretizationType::pyramid5:
    {
      return DefineProblemType<DRT::Element::DiscretizationType::pyramid5>(problem);
    }
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
    // Nurbs support
    case DRT::Element::DiscretizationType::nurbs9:
    {
      return DefineProblemType<DRT::Element::DiscretizationType::nurbs9>(problem);
    }
    case DRT::Element::DiscretizationType::nurbs27:
    {
      return DefineProblemType<DRT::Element::DiscretizationType::nurbs27>(problem);
    }
    // no 1D elements
    default:
      dserror("Element shape %s not activated. Just do it.", DRT::DistypeToString(distype).c_str());
      break;
  }
  return nullptr;
}

/*--------------------------------------------------------------------------*
 |                                                (public) berardocco 02/18 |
 *--------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ElemagEleInterface* DRT::ELEMENTS::ElemagFactory::DefineProblemType(
    std::string problem)
{
  if (problem == "std")
    return DRT::ELEMENTS::ElemagEleCalc<distype>::Instance();
  else if (problem == "diff")
    return DRT::ELEMENTS::ElemagDiffEleCalc<distype>::Instance();
  else
    dserror("Defined problem type does not exist!!");

  return nullptr;
}
