/*----------------------------------------------------------------------*/
/*! \file
\brief pyramid shaped solid element
\level 2

*----------------------------------------------------------------------*/

#include "4C_io_linedefinition.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_so3_pyramid5fbar.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool Discret::ELEMENTS::SoPyramid5fbar::ReadElement(
    const std::string& eletype, const std::string& distype, Input::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->extract_int("MAT", material);
  SetMaterial(0, Mat::Factory(material));

  // set up of materials with GP data (e.g., history variables)
  SolidMaterial()->setup(NUMGPT_SOP5, linedef);

  // temporary variable for read-in
  std::string buffer;

  // read kinematic flag
  linedef->extract_string("KINEM", buffer);
  if (buffer == "linear")
  {
    FOUR_C_THROW("Only nonlinear kinematics for SO_PYRAMID5FBAR implemented!");
  }
  else if (buffer == "nonlinear")
  {
    kintype_ = Inpar::STR::KinemType::nonlinearTotLag;
  }
  else
    FOUR_C_THROW("Reading SO_PYRAMID5FBAR element failed KINEM unknown");

  // check if material kinematics is compatible to element kinematics
  SolidMaterial()->ValidKinematics(kintype_);

  // Validate that materials doesn't use extended update call.
  if (SolidMaterial()->UsesExtendedUpdate())
    FOUR_C_THROW("This element currently does not support the extended update call.");

  return true;
}

FOUR_C_NAMESPACE_CLOSE
