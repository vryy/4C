/*----------------------------------------------------------------------*/
/*! \file
\brief Solid Hex8 element

\level 1

*----------------------------------------------------------------------*/

#include "4C_io_linedefinition.hpp"
#include "4C_mat_robinson.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_so3_hex8.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool Discret::ELEMENTS::SoHex8::read_element(
    const std::string& eletype, const std::string& distype, Input::LineDefinition* linedef)
{
  // read number of material model
  int material_id = 0;
  linedef->extract_int("MAT", material_id);

  set_material(0, Mat::Factory(material_id));

  // set up of materials with GP data (e.g., history variables)
  solid_material()->setup(NUMGPT_SOH8, linedef);

  // temporary variable for read-in
  std::string buffer;

  // read kinematic flag
  linedef->extract_string("KINEM", buffer);
  if (buffer == "linear")
  {
    kintype_ = Inpar::Solid::KinemType::linear;
  }
  else if (buffer == "nonlinear")
  {
    kintype_ = Inpar::Solid::KinemType::nonlinearTotLag;
  }
  else
    FOUR_C_THROW("Reading SO_HEX8 element failed");

  // check if material kinematics is compatible to element kinematics
  solid_material()->valid_kinematics(kintype_);

  // read EAS technology flag
  linedef->extract_string("EAS", buffer);

  // full EAS technology
  if (buffer == "full")
  {
    eastype_ = soh8_easfull;
    neas_ = 21;  // number of eas parameters for full EAS
    soh8_easinit();
  }
  // mild EAS technology
  else if (buffer == "mild")
  {
    eastype_ = soh8_easmild;
    neas_ = 9;  // number of eas parameters for mild EAS
    soh8_easinit();
  }
  // no EAS technology
  else if (buffer == "none")
  {
    eastype_ = soh8_easnone;
    neas_ = 0;
  }
  else
    FOUR_C_THROW("Reading of SO_HEX8 EAS technology failed");

  return true;
}

FOUR_C_NAMESPACE_CLOSE
