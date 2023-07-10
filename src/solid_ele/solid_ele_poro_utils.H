/*----------------------------------------------------------------------*/
/*! \file

\brief little helpers for solid elements

\level 1

*----------------------------------------------------------------------*/

#ifndef SOLID_ELE_PORO_UTILS_H
#define SOLID_ELE_PORO_UTILS_H

#include "discretization_fem_general_utils_integration.H"
#include "linalg_utils_densematrix_eigen.H"
#include "lib_element.H"
#include "lib_linedefinition.H"
#include "solid_ele_poro.H"

namespace STR::UTILS
{
  namespace READELEMENT
  {
    INPAR::SCATRA::ImplType ReadType(DRT::INPUT::LineDefinition* linedef)
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
        dserror("Invalid TYPE for SOLIDPORO elements!");
        return INPAR::SCATRA::impltype_undefined;
      }
    }

    INPAR::PORO::PoroType ReadPoroType(DRT::INPUT::LineDefinition* linedef)
    {
      std::string impltype;
      linedef->ExtractString("POROTYPE", impltype);

      if (impltype == "PressureVelocityBased")
        return INPAR::PORO::PoroType::porotype_pressure_velocity_based;
      else if (impltype == "PressureBased")
        return INPAR::PORO::PoroType::porotype_pressure_based;
      else
      {
        dserror("Invalid POROTYPE for SolidPoro elements!");
        return INPAR::PORO::PoroType::porotype_undefined;
      }
    }
  }  // namespace READELEMENT



}  // namespace STR::UTILS

#endif  // SOLID_ELE_PORO_UTILS_H