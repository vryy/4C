/*----------------------------------------------------------------------*/
/*! \file
\brief Solid Hex8 element with F-bar modification

\level 1

*----------------------------------------------------------------------*/

#include "4C_io_linedefinition.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_so3_hex8fbar.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool Discret::ELEMENTS::SoHex8fbar::ReadElement(
    const std::string& eletype, const std::string& distype, Input::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->extract_int("MAT", material);
  SetMaterial(0, Mat::Factory(material));

  // set up of materials with GP data (e.g., history variables)
  SolidMaterial()->setup(NUMGPT_SOH8, linedef);

  // temporary variable for read-in
  std::string buffer;

  // read kinematic flag
  linedef->extract_string("KINEM", buffer);
  if (buffer == "linear")
  {
    FOUR_C_THROW("Only nonlinear kinematics for SO_HEX8FBAR implemented!");
  }
  else if (buffer == "nonlinear")
  {
    kintype_ = Inpar::Solid::KinemType::nonlinearTotLag;
  }
  else
    FOUR_C_THROW("Reading SO_HEX8FBAR element failed KINEM unknown");

  // check if material kinematics is compatible to element kinematics
  SolidMaterial()->ValidKinematics(kintype_);

  return true;
}

FOUR_C_NAMESPACE_CLOSE
