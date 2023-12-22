/*----------------------------------------------------------------------*/
/*! \file
\brief pyramid shaped solid element
\level 2

*----------------------------------------------------------------------*/

#include "baci_io_linedefinition.H"
#include "baci_mat_so3_material.H"
#include "baci_so3_pyramid5fbar.H"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_pyramid5fbar::ReadElement(
    const std::string& eletype, const std::string& distype, INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(material);

  // set up of materials with GP data (e.g., history variables)
  SolidMaterial()->Setup(NUMGPT_SOP5, linedef);

  // temporary variable for read-in
  std::string buffer;

  // read kinematic flag
  linedef->ExtractString("KINEM", buffer);
  if (buffer == "linear")
  {
    dserror("Only nonlinear kinematics for SO_PYRAMID5FBAR implemented!");
  }
  else if (buffer == "nonlinear")
  {
    kintype_ = INPAR::STR::KinemType::nonlinearTotLag;
  }
  else
    dserror("Reading SO_PYRAMID5FBAR element failed KINEM unknown");

  // check if material kinematics is compatible to element kinematics
  SolidMaterial()->ValidKinematics(kintype_);

  // Validate that materials doesn't use extended update call.
  if (SolidMaterial()->UsesExtendedUpdate())
    dserror("This element currently does not support the extended update call.");

  return true;
}

BACI_NAMESPACE_CLOSE
