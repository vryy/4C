/*! \file

\brief Helpers for solid-poro elements

\level 1
*/

#ifndef FOUR_C_SOLID_PORO_3D_ELE_UTILS_HPP
#define FOUR_C_SOLID_PORO_3D_ELE_UTILS_HPP

#include "4C_config.hpp"

#include "4C_inpar_poro.hpp"
#include "4C_inpar_scatra.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Solid::Utils::read_element
{
  inline Inpar::ScaTra::ImplType read_type(const Core::IO::InputParameterContainer& container)
  {
    auto impltype = container.get_or<std::string>("TYPE", "Undefined");

    if (impltype == "Undefined")
      return Inpar::ScaTra::impltype_undefined;
    else if (impltype == "AdvReac")
      return Inpar::ScaTra::impltype_advreac;
    else if (impltype == "CardMono")
      return Inpar::ScaTra::impltype_cardiac_monodomain;
    else if (impltype == "Chemo")
      return Inpar::ScaTra::impltype_chemo;
    else if (impltype == "ChemoReac")
      return Inpar::ScaTra::impltype_chemoreac;
    else if (impltype == "Loma")
      return Inpar::ScaTra::impltype_loma;
    else if (impltype == "Poro")
      return Inpar::ScaTra::impltype_poro;
    else if (impltype == "PoroReac")
      return Inpar::ScaTra::impltype_pororeac;
    else if (impltype == "PoroReacECM")
      return Inpar::ScaTra::impltype_pororeacECM;
    else if (impltype == "PoroMultiReac")
      return Inpar::ScaTra::impltype_multipororeac;
    else if (impltype == "RefConcReac")
      return Inpar::ScaTra::impltype_refconcreac;
    else if (impltype == "Std")
      return Inpar::ScaTra::impltype_std;
    else
    {
      FOUR_C_THROW("Invalid TYPE for SOLIDPORO_PRESSURE_BASED elements!");
      return Inpar::ScaTra::impltype_undefined;
    }
  }

}  // namespace Solid::Utils::read_element

FOUR_C_NAMESPACE_CLOSE

#endif