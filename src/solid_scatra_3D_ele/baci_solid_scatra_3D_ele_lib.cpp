/*! \file

\brief Implementation of the solid-scatra element

\level 1
*/
#include "baci_solid_scatra_3D_ele_lib.hpp"

#include "baci_utils_exceptions.hpp"

BACI_NAMESPACE_OPEN

INPAR::SCATRA::ImplType DRT::ELEMENTS::ReadScatraImplType(
    const INPUT::LineDefinition& line_definition)
{
  std::string impltype;
  line_definition.ExtractString("TYPE", impltype);

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
  else if (impltype == "ElchDiffCond")
    return INPAR::SCATRA::impltype_elch_diffcond;
  else if (impltype == "ElchElectrode")
    return INPAR::SCATRA::impltype_elch_electrode;
  else if (impltype == "Loma")
    return INPAR::SCATRA::impltype_loma;
  else if (impltype == "RefConcReac")
    return INPAR::SCATRA::impltype_refconcreac;
  else if (impltype == "Std")
    return INPAR::SCATRA::impltype_std;

  dserror("The input type %s is not valud for SOLIDSCATRA elements!", impltype.c_str());
}

BACI_NAMESPACE_CLOSE
