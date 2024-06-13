/*! \file

\brief Properties of solid elements

\level 1
*/

#ifndef FOUR_C_SOLID_3D_ELE_PROPERTIES_HPP
#define FOUR_C_SOLID_3D_ELE_PROPERTIES_HPP

#include "4C_config.hpp"

#include "4C_inpar_structure.hpp"
#include "4C_utils_exceptions.hpp"

#include <string>
FOUR_C_NAMESPACE_OPEN

namespace Core::Communication
{
  class PackBuffer;
}
namespace STR::ELEMENTS
{
  enum class EasType;
}
namespace Discret::ELEMENTS
{
  enum class ElementTechnology
  {
    none,
    fbar,
    eas_mild,
    eas_full
  };

  static inline std::string element_technology_string(const ElementTechnology ele_tech)
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

    FOUR_C_THROW("Unknown element technology %d", ele_tech);
  }

  template <typename Function>
  auto element_technology_switch(ElementTechnology eletech, Function fct)
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

    FOUR_C_THROW("Your element technology is unknown: %d", eletech);
  }

  enum class PrestressTechnology
  {
    none,
    mulf
  };

  static inline std::string prestress_technology_string(const PrestressTechnology prestress_tech)
  {
    switch (prestress_tech)
    {
      case PrestressTechnology::none:
        return "none";
      case PrestressTechnology::mulf:
        return "mulf";
    }

    FOUR_C_THROW("Unknown prestress technology %d", prestress_tech);
  }

  template <typename Function>
  auto prestress_technology_switch(PrestressTechnology prestress_technology, Function fct)
  {
    switch (prestress_technology)
    {
      case PrestressTechnology::none:
        return fct(std::integral_constant<PrestressTechnology, PrestressTechnology::none>{});
      case PrestressTechnology::mulf:
        return fct(std::integral_constant<PrestressTechnology, PrestressTechnology::mulf>{});
    }

    FOUR_C_THROW("Your prestress technology is unknown: %d", prestress_technology);
  }


  /*!
   *  @brief struct for managing solid element properties
   */
  struct SolidElementProperties
  {
    //! kinematic type
    Inpar::STR::KinemType kintype{Inpar::STR::KinemType::vague};

    //! element technology (none, F-Bar, EAS full, EAS mild)
    ElementTechnology element_technology{ElementTechnology::none};

    //! specify prestress technology (none, MULF)
    PrestressTechnology prestress_technology{PrestressTechnology::none};
  };

  void add_to_pack(Core::Communication::PackBuffer& data,
      const Discret::ELEMENTS::SolidElementProperties& properties);

  void ExtractFromPack(std::size_t& position, const std::vector<char>& data,
      Discret::ELEMENTS::SolidElementProperties& properties);

}  // namespace Discret::ELEMENTS


FOUR_C_NAMESPACE_CLOSE

#endif
