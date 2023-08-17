/*----------------------------------------------------------------------*/
/*! \file
\brief 3D quadratic serendipity element
\level 1

*----------------------------------------------------------------------*/

#include "baci_lib_linedefinition.H"
#include "baci_mat_so3_material.H"
#include "baci_so3_hex20.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_hex20::ReadElement(
    const std::string& eletype, const std::string& distype, DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(material);

  SolidMaterial()->Setup(NUMGPT_SOH20, linedef);

  std::string buffer;
  linedef->ExtractString("KINEM", buffer);

  if (buffer == "linear")
  {
    kintype_ = INPAR::STR::kinem_linear;
  }
  else if (buffer == "nonlinear")
  {
    kintype_ = INPAR::STR::kinem_nonlinearTotLag;
  }
  else
    dserror("Reading SO_HEX20 element failed KINEM unknown");

  // check if material kinematics is compatible to element kinematics
  SolidMaterial()->ValidKinematics(kintype_);

  // Validate that materials doesn't use extended update call.
  if (SolidMaterial()->UsesExtendedUpdate())
    dserror("This element currently does not support the extended update call.");

  return true;
}
