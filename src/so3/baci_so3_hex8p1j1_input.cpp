/*----------------------------------------------------------------------*/
/*! \file
\brief 'Q1P0' element in 8-node hexahedron shape

\level 2

*/
/*----------------------------------------------------------------------*/


#include "baci_io_linedefinition.H"
#include "baci_mat_so3_material.H"
#include "baci_so3_hex8p1j1.H"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_Hex8P1J1::ReadElement(
    const std::string& eletype, const std::string& distype, INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(material);

  // set up of materials with GP data (e.g., history variables)
  SolidMaterial()->Setup(NUMGPT_SOH8, linedef);

  // temporary variable for read-in
  std::string buffer;
  // read kinematic flag
  linedef->ExtractString("KINEM", buffer);
  // no linear case implemented so far, hence just a dummy check
  if (buffer == "linear")
  {
    // kintype_ = soh8_linear;
    dserror("Only nonlinear kinematics for SO_HEX8p1j1 implemented!");
  }
  else if (buffer == "nonlinear")
  {
    kintype_ = INPAR::STR::KinemType::nonlinearTotLag;
  }
  else
    dserror("Reading SO_HEX8p1j1 element failed KINEM unknown");

  // check if material kinematics is compatible to element kinematics
  SolidMaterial()->ValidKinematics(kintype_);

  // Validate that materials doesn't use extended update call.
  if (SolidMaterial()->UsesExtendedUpdate())
    dserror("This element currently does not support the extended update call.");

  return true;
}

BACI_NAMESPACE_CLOSE
