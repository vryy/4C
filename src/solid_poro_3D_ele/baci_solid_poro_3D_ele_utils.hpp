/*! \file

\brief Helpers for solid-poro elements

\level 1
*/

#ifndef FOUR_C_SOLID_PORO_3D_ELE_UTILS_HPP
#define FOUR_C_SOLID_PORO_3D_ELE_UTILS_HPP

#include "baci_config.hpp"

#include "baci_inpar_poro.hpp"
#include "baci_inpar_scatra.hpp"
#include "baci_io_linedefinition.hpp"
#include "baci_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

namespace STR::UTILS::READELEMENT
{
  INPAR::SCATRA::ImplType ReadType(INPUT::LineDefinition* linedef)
  {
    std::string impltype;
    linedef->ExtractString("TYPE", impltype);

    if (impltype == "Undefined")
      return INPAR::SCATRA::impltype_undefined;
    else if (impltype == "AdvReac")
      return INPAR::SCATRA::impltype_advreac;
    else if (impltype == "CardMono")
      return INPAR::SCATRA::impltype_cardiac_monodomain;
    else if (impltype == "Chemo")
      return INPAR::SCATRA::impltype_chemo;
    else if (impltype == "ChemoReac")
      return INPAR::SCATRA::impltype_chemoreac;
    else if (impltype == "Loma")
      return INPAR::SCATRA::impltype_loma;
    else if (impltype == "Poro")
      return INPAR::SCATRA::impltype_poro;
    else if (impltype == "PoroReac")
      return INPAR::SCATRA::impltype_pororeac;
    else if (impltype == "PoroReacECM")
      return INPAR::SCATRA::impltype_pororeacECM;
    else if (impltype == "PoroMultiReac")
      return INPAR::SCATRA::impltype_multipororeac;
    else if (impltype == "RefConcReac")
      return INPAR::SCATRA::impltype_refconcreac;
    else if (impltype == "Std")
      return INPAR::SCATRA::impltype_std;
    else
    {
      FOUR_C_THROW("Invalid TYPE for SOLIDPORO elements!");
      return INPAR::SCATRA::impltype_undefined;
    }
  }

  INPAR::PORO::PoroType ReadPoroType(INPUT::LineDefinition* linedef)
  {
    std::string impltype;
    linedef->ExtractString("POROTYPE", impltype);

    if (impltype == "PressureVelocityBased")
      return INPAR::PORO::PoroType::pressure_velocity_based;
    else if (impltype == "PressureBased")
      return INPAR::PORO::PoroType::pressure_based;
    else
    {
      FOUR_C_THROW("Invalid POROTYPE for SolidPoro elements!");
      return INPAR::PORO::PoroType::undefined;
    }
  }



}  // namespace STR::UTILS::READELEMENT

FOUR_C_NAMESPACE_CLOSE

#endif