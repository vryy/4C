/*----------------------------------------------------------------------*/
/*! \file

\brief Factory of solid elements

\level 1
*/
/*--------------------------------------------------------------------------*/

#include "solid_ele_factory.H"
#include "solid_ele_calc.H"
#include "solid_ele_calc_eas.H"
#include "../drt_inpar/inpar_structure.H"
#include "solid_ele.H"


DRT::ELEMENTS::SolidEleInterface* DRT::ELEMENTS::SolidFactory::ProvideImpl(
    DRT::ELEMENTS::Solid* ele)
{
  switch (ele->Shape())
  {
    case DRT::Element::hex8:
      return ProvideImpl<DRT::Element::hex8>(ele);
      break;
    case DRT::Element::hex27:
      return ProvideImpl<DRT::Element::hex27>(ele);
      break;
    default:
      dserror("unknown distype provided");
      break;
  }
  return nullptr;
}

template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::SolidEleInterface* DRT::ELEMENTS::SolidFactory::ProvideImpl(
    DRT::ELEMENTS::Solid* ele)
{
  // here we go into the different cases for element technology
  switch (ele->GetEleTech().size())
  {
    // no element technology
    case 0:
      //    return DRT::ELEMENTS::SolidEleCalc<distype>::Instance();
      return DRT::ELEMENTS::SolidEleCalc<DRT::Element::hex8>::Instance();
      break;
    // simple: just one element technology
    case 1:
      switch (*ele->GetEleTech().begin())
      {
        case INPAR::STR::eletech_eas:
          //      dserror("not implemented");
          switch (ele->GetEAStype())
          {
            case ::STR::ELEMENTS::EASType::eastype_h8_9:
              if constexpr (distype != DRT::Element::hex8)
              {
                dserror("EAS type h8_9 is only for hex8 elements (you are using %s)", distype);
              }
              return DRT::ELEMENTS::SolidEleCalcEas<DRT::Element::hex8, 9>::Instance();
              break;
            case ::STR::ELEMENTS::EASType::eastype_h8_21:
              if constexpr (distype != DRT::Element::hex8)
              {
                dserror("EAS type h8_21 is only for hex8 elements (you are using %s)", distype);
              }
              return DRT::ELEMENTS::SolidEleCalcEas<DRT::Element::hex8, 21>::Instance();
              break;
            default:
              dserror("eastype not implemented %d", (int)ele->GetEAStype());
              break;
          }
          break;
        default:
          dserror("unknown element technology");
          break;
      }
    // complicated: combination of element technologies
    default:
    {
      if (ele->GetEleTech().size() == 2 and
          ele->GetEleTech().find(INPAR::STR::eletech_eas) != ele->GetEleTech().end() and
          ele->GetEleTech().find(INPAR::STR::eletech_plasticity) != ele->GetEleTech().end())
        //      return DRT::ELEMENTS::SolidEleCalcEasPlast<distype>::Instance();
        dserror("this combination of element technology is not implemented yet");
      else
      {
        dserror("unknown combination of element technologies");
      }
    }
    break;
  }
  return nullptr;
}
