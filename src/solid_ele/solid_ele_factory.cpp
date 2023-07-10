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
    DRT::Element* ele, std::set<INPAR::STR::EleTech> eletech, ::INPAR::STR::KinemType kinemtype,
    ::STR::ELEMENTS::EASType eastype)
{
  switch (ele->Shape())
  {
    case DRT::Element::hex8:
      return ProvideImpl<DRT::Element::hex8>(eletech, kinemtype, eastype);
      break;
    case DRT::Element::hex27:
      return ProvideImpl<DRT::Element::hex27>(eletech, kinemtype, eastype);
      break;
    case DRT::Element::hex20:
      return ProvideImpl<DRT::Element::hex20>(eletech, kinemtype, eastype);
      break;
    case DRT::Element::hex18:
      return ProvideImpl<DRT::Element::hex18>(eletech, kinemtype, eastype);
      break;
    case DRT::Element::pyramid5:
      return ProvideImpl<DRT::Element::pyramid5>(eletech, kinemtype, eastype);
      break;
    case DRT::Element::wedge6:
      return ProvideImpl<DRT::Element::wedge6>(eletech, kinemtype, eastype);
      break;
    case DRT::Element::tet4:
      return ProvideImpl<DRT::Element::tet4>(eletech, kinemtype, eastype);
      break;
    case DRT::Element::tet10:
      return ProvideImpl<DRT::Element::tet10>(eletech, kinemtype, eastype);
      break;
    default:
      dserror("unknown distype provided");
      break;
  }
  return nullptr;
}

template <DRT::Element::DiscretizationType distype>
std::unique_ptr<DRT::ELEMENTS::SolidEleCalcInterface> DRT::ELEMENTS::SolidFactory::ProvideImpl(
    std::set<INPAR::STR::EleTech> eletech, ::INPAR::STR::KinemType kinemtype,
    ::STR::ELEMENTS::EASType eastype)
{
  // here we go into the different cases for element technology
  switch (eletech.size())
  {
    // no element technology
    case 0:
      return std::make_unique<DRT::ELEMENTS::SolidEleCalc<distype>>();
      break;
    // simple: just one element technology
    case 1:
      switch (*eletech.begin())
      {
        case INPAR::STR::EleTech::eas:
          /* switch (eastype)
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
              dserror("eastype not implemented %d", (int)eastype);
              break;
          } */
        case INPAR::STR::EleTech::fbar:
          if constexpr (distype != DRT::Element::hex8)
          {
            dserror("FBAR is only implemented for hex8 elements.");
          }
          if (kinemtype != INPAR::STR::KinemType::kinem_nonlinearTotLag)
          {
            dserror("FBAR only usable for KINEM nonlinear (you are using %s)", kinemtype);
          }
          // template only instantiated for distype == DRT::Element::hex8,
          // hence this check is necessary
          if constexpr (distype == DRT::Element::hex8)
          {
            return std::make_unique<DRT::ELEMENTS::SolidEleCalcFbar<distype>>();
          }
          break;
        default:
          dserror("unknown element technology");
          break;
      }
    // complicated: combination of element technologies
    default:
    {
      if (eletech.size() == 2 and eletech.find(INPAR::STR::EleTech::eas) != eletech.end() and
          eletech.find(INPAR::STR::EleTech::plasticity) != eletech.end())
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
