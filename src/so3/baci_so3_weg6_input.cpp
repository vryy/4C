/*----------------------------------------------------------------------*/
/*! \file
\brief Solid Wedge6 Element
\level 1

*----------------------------------------------------------------------*/

#include "baci_io_linedefinition.H"
#include "baci_mat_so3_material.H"
#include "baci_so3_weg6.H"

BACI_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_weg6::ReadElement(
    const std::string& eletype, const std::string& distype, INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(material);

  SolidMaterial()->Setup(NUMGPT_WEG6, linedef);

  std::string buffer;
  linedef->ExtractString("KINEM", buffer);
  if (buffer == "linear")
  {
    kintype_ = INPAR::STR::KinemType::linear;
    dserror("Reading of SO_WEG6 element failed only nonlinear kinematics implemented");
  }
  else if (buffer == "nonlinear")
  {
    kintype_ = INPAR::STR::KinemType::nonlinearTotLag;
  }
  else
    dserror("Reading SO_WEG6 element failed KINEM unknwon");

  // check if material kinematics is compatible to element kinematics
  SolidMaterial()->ValidKinematics(kintype_);

  // Validate that materials doesn't use extended update call.
  if (SolidMaterial()->UsesExtendedUpdate())
    dserror("This element currently does not support the extended update call.");

  return true;
}

BACI_NAMESPACE_CLOSE
