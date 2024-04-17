/*! \file

\brief Properties of solid elements

\level 1
*/

#ifndef FOUR_C_SOLID_3D_ELE_PROPERTIES_HPP
#define FOUR_C_SOLID_3D_ELE_PROPERTIES_HPP

#include "baci_config.hpp"

#include "baci_inpar_structure.hpp"
#include "baci_utils_exceptions.hpp"

#include <string>
FOUR_C_NAMESPACE_OPEN

namespace CORE::COMM
{
  class PackBuffer;
}
namespace STR::ELEMENTS
{
  enum class EasType;
}
namespace DRT::ELEMENTS
{
  enum class ElementTechnology
  {
    none,
    fbar,
    eas_mild,
    eas_full
  };

  static inline std::string ElementTechnologyString(const ElementTechnology ele_tech)
  {
    switch (ele_tech)
    {
      case ElementTechnology::none:
        return "none";
      case ElementTechnology::fbar:
        return "fbar";
      case ElementTechnology::eas_mild:
        return "eas_mild";
      case ElementTechnology::eas_full:
        return "eas_full";
    }

    dserror("Unknown element technology %d", ele_tech);
  }

  template <typename Function>
  auto ElementTechnologySwitch(ElementTechnology eletech, Function fct)
  {
    switch (eletech)
    {
      case ElementTechnology::none:
        return fct(std::integral_constant<ElementTechnology, ElementTechnology::none>{});
      case ElementTechnology::fbar:
        return fct(std::integral_constant<ElementTechnology, ElementTechnology::fbar>{});
      case ElementTechnology::eas_mild:
        return fct(std::integral_constant<ElementTechnology, ElementTechnology::eas_mild>{});
      case ElementTechnology::eas_full:
        return fct(std::integral_constant<ElementTechnology, ElementTechnology::eas_full>{});
    }

    dserror("Your element technology is unknown: %d", eletech);
  }

  enum class PrestressTechnology
  {
    none,
    mulf
  };

  static inline std::string PrestressTechnologyString(const PrestressTechnology prestress_tech)
  {
    switch (prestress_tech)
    {
      case PrestressTechnology::none:
        return "none";
      case PrestressTechnology::mulf:
        return "mulf";
    }

    dserror("Unknown prestress technology %d", prestress_tech);
  }

  template <typename Function>
  auto PrestressTechnologySwitch(PrestressTechnology prestress_technology, Function fct)
  {
    switch (prestress_technology)
    {
      case PrestressTechnology::none:
        return fct(std::integral_constant<PrestressTechnology, PrestressTechnology::none>{});
      case PrestressTechnology::mulf:
        return fct(std::integral_constant<PrestressTechnology, PrestressTechnology::mulf>{});
    }

    dserror("Your prestress technology is unknown: %d", prestress_technology);
  }


  /*!
   *  @brief struct for managing solid element properties
   */
  struct SolidElementProperties
  {
    //! kinematic type
    INPAR::STR::KinemType kintype{INPAR::STR::KinemType::vague};

    //! element technology (none, F-Bar, EAS full, EAS mild)
    ElementTechnology element_technology{ElementTechnology::none};

    //! specify prestress technology (none, MULF)
    PrestressTechnology prestress_technology{PrestressTechnology::none};
  };

  void AddToPack(
      CORE::COMM::PackBuffer& data, const DRT::ELEMENTS::SolidElementProperties& properties);

  void ExtractFromPack(std::size_t& position, const std::vector<char>& data,
      DRT::ELEMENTS::SolidElementProperties& properties);

}  // namespace DRT::ELEMENTS


FOUR_C_NAMESPACE_CLOSE

#endif
