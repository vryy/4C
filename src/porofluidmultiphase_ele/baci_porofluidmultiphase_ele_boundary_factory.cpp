/*----------------------------------------------------------------------*/
/*! \file
 \brief factory class providing the implementations of the porofluidmultiphase
        boundary element evaluation routines

   \level 3

 *----------------------------------------------------------------------*/


#include "baci_porofluidmultiphase_ele_boundary_factory.H"

#include "baci_lib_element.H"
#include "baci_lib_globalproblem.H"
#include "baci_porofluidmultiphase_ele_boundary_calc.H"
#include "baci_porofluidmultiphase_ele_interface.H"


/*--------------------------------------------------------------------------*
 | provide the implementation of evaluation class      (public) vuong 08/16 |
 *--------------------------------------------------------------------------*/
DRT::ELEMENTS::PoroFluidMultiPhaseEleInterface*
DRT::ELEMENTS::PoroFluidMultiPhaseBoundaryFactory::ProvideImpl(
    const DRT::Element* ele, const int numdofpernode, const std::string& disname)
{
  switch (ele->Shape())
  {
    case DRT::Element::DiscretizationType::quad4:
    {
      return DefineProblemType<DRT::Element::DiscretizationType::quad4>(numdofpernode, disname);
    }
    case DRT::Element::DiscretizationType::quad8:
    {
      return DefineProblemType<DRT::Element::DiscretizationType::quad8>(numdofpernode, disname);
    }
    case DRT::Element::DiscretizationType::quad9:
    {
      return DefineProblemType<DRT::Element::DiscretizationType::quad9>(numdofpernode, disname);
    }
    case DRT::Element::DiscretizationType::tri3:
    {
      return DefineProblemType<DRT::Element::DiscretizationType::tri3>(numdofpernode, disname);
    }
    case DRT::Element::DiscretizationType::tri6:
    {
      return DefineProblemType<DRT::Element::DiscretizationType::tri6>(numdofpernode, disname);
    }
    case DRT::Element::DiscretizationType::line2:
    {
      return DefineProblemType<DRT::Element::DiscretizationType::line2>(numdofpernode, disname);
    }
    case DRT::Element::DiscretizationType::line3:
    {
      return DefineProblemType<DRT::Element::DiscretizationType::line3>(numdofpernode, disname);
    }
    default:
    {
      dserror(
          "Element shape %d (%d nodes) not activated. Just do it.", ele->Shape(), ele->NumNode());
      break;
    }
  }

  return nullptr;
}


/*--------------------------------------------------------------------------*
 | provide the implementation of evaluation class      (public) vuong 08/16 |
 *--------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::PoroFluidMultiPhaseEleInterface*
DRT::ELEMENTS::PoroFluidMultiPhaseBoundaryFactory::DefineProblemType(
    const int numdofpernode, const std::string& disname)
{
  return DRT::ELEMENTS::PoroFluidMultiPhaseEleBoundaryCalc<distype>::Instance(
      numdofpernode, disname);
}
