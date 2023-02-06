/*----------------------------------------------------------------------*/
/*! \file

\brief Factory of solid elements

\level 1
*/
/*--------------------------------------------------------------------------*/

#include "solid_ele_factory.H"
#include <memory>
// #include "solid_ele_calc_eas.H"
#include "solid_ele_calc_fbar.H"
#include "solid_ele_calc.H"
#include "inpar_structure.H"
#include "solid_ele.H"


std::unique_ptr<DRT::ELEMENTS::SolidEleCalcInterface> DRT::ELEMENTS::SolidFactory::ProvideImpl(
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
    case DRT::Element::hex20:
      return ProvideImpl<DRT::Element::hex20>(ele);
      break;
    case DRT::Element::hex18:
      return ProvideImpl<DRT::Element::hex18>(ele);
      break;
    case DRT::Element::pyramid5:
      return ProvideImpl<DRT::Element::pyramid5>(ele);
      break;
    case DRT::Element::wedge6:
      return ProvideImpl<DRT::Element::wedge6>(ele);
      break;
    case DRT::Element::tet4:
      return ProvideImpl<DRT::Element::tet4>(ele);
      break;
    case DRT::Element::tet10:
      return ProvideImpl<DRT::Element::tet10>(ele);
      break;
    default:
      dserror("unknown distype provided");
      break;
  }
  return nullptr;
}

template <DRT::Element::DiscretizationType distype>
std::unique_ptr<DRT::ELEMENTS::SolidEleCalcInterface> DRT::ELEMENTS::SolidFactory::ProvideImpl(
    DRT::ELEMENTS::Solid* ele)
{
  // here we go into the different cases for element technology
  switch (ele->GetEleTech().size())
  {
    // no element technology
    case 0:
      return std::make_unique<DRT::ELEMENTS::SolidEleCalc<distype>>();
      break;
    // simple: just one element technology
    case 1:
      switch (*ele->GetEleTech().begin())
      {
        case INPAR::STR::EleTech::eas:
          /* switch (ele->GetEAStype())
          {
            case ::STR::ELEMENTS::EASType::eastype_h8_9:
              if constexpr (distype != DRT::Element::hex8)
              {
                dserror("EAS type h8_9 is only for hex8 elements (you are using %s)", distype);
              }
              return std::make_unique<DRT::ELEMENTS::SolidEleCalcEas<DRT::Element::hex8, 9>>();
              break;
            case ::STR::ELEMENTS::EASType::eastype_h8_21:
              if constexpr (distype != DRT::Element::hex8)
              {
                dserror("EAS type h8_21 is only for hex8 elements (you are using %s)", distype);
              }
              return std::make_unique<DRT::ELEMENTS::SolidEleCalcEas<DRT::Element::hex8, 21>>();
              break;
            default:
              dserror("eastype not implemented %d", (int)ele->GetEAStype());
              break;
          } */
        case INPAR::STR::EleTech::fbar:
          if constexpr (distype != DRT::Element::hex8)
          {
            dserror("FBAR is only implemented for hex8 elements.");
          }
          if (ele->GetKinemType() != INPAR::STR::KinemType::kinem_nonlinearTotLag)
          {
            dserror("FBAR only usable for KINEM nonlinear (you are using %s)", ele->GetKinemType());
          }
          // template only instantiated for distype == DRT::Element::hex8,
          // hence this check is necessary
          if constexpr (distype == DRT::Element::hex8)
            return std::make_unique<DRT::ELEMENTS::SolidEleCalcFbar<distype>>();
          break;
        default:
          dserror("unknown element technology");
          break;
      }
    // complicated: combination of element technologies
    default:
    {
      if (ele->GetEleTech().size() == 2 and
          ele->GetEleTech().find(INPAR::STR::EleTech::eas) != ele->GetEleTech().end() and
          ele->GetEleTech().find(INPAR::STR::EleTech::plasticity) != ele->GetEleTech().end())
        //      return DRT::ELEMENTS::SolidEleCalcEasPlast<distype>::Instance();
        dserror("this combination of element technology is not implemented yet");
      else
      {
        dserror("unknown combination of element technologies");
      }
    }
  }
  return nullptr;
}
