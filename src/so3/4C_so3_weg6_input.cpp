/*----------------------------------------------------------------------*/
/*! \file
\brief Solid Wedge6 Element
\level 1

*----------------------------------------------------------------------*/

#include "4C_io_linedefinition.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_so3_weg6.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool Discret::ELEMENTS::SoWeg6::ReadElement(
    const std::string& eletype, const std::string& distype, Input::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->extract_int("MAT", material);
  SetMaterial(0, Mat::Factory(material));

  SolidMaterial()->setup(NUMGPT_WEG6, linedef);

  std::string buffer;
  linedef->extract_string("KINEM", buffer);
  if (buffer == "linear")
  {
    kintype_ = Inpar::Solid::KinemType::linear;
    FOUR_C_THROW("Reading of SO_WEG6 element failed only nonlinear kinematics implemented");
  }
  else if (buffer == "nonlinear")
  {
    kintype_ = Inpar::Solid::KinemType::nonlinearTotLag;
  }
  else
    FOUR_C_THROW("Reading SO_WEG6 element failed KINEM unknwon");

  // check if material kinematics is compatible to element kinematics
  SolidMaterial()->ValidKinematics(kintype_);

  // Validate that materials doesn't use extended update call.
  if (SolidMaterial()->UsesExtendedUpdate())
    FOUR_C_THROW("This element currently does not support the extended update call.");

  return true;
}

FOUR_C_NAMESPACE_CLOSE
