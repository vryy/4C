/*----------------------------------------------------------------------*/
/*! \file
\brief Input for solid shell 18
\level 3

*----------------------------------------------------------------------*/

#include "4C_io_linedefinition.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_so3_sh18.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool Discret::ELEMENTS::SoSh18::ReadElement(
    const std::string& eletype, const std::string& distype, Input::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->extract_int("MAT", material);

  SetMaterial(0, Mat::Factory(material));

  // set up of materials with GP data (e.g., history variables)

  Teuchos::RCP<Core::Mat::Material> mat = Material();

  SolidMaterial()->setup(NUMGPT_SOH18, linedef);

  // temporary variable for read-in
  std::string buffer;

  // read kinematic flag
  linedef->extract_string("KINEM", buffer);
  if (buffer == "linear")
  {
    // kintype_ = soh8_linear;
    FOUR_C_THROW("Only nonlinear kinematics for SO_SH8 implemented!");
  }
  else if (buffer == "nonlinear")
  {
    kintype_ = Inpar::STR::KinemType::nonlinearTotLag;
  }
  else
    FOUR_C_THROW("Reading SO_HEX18 element failed KINEM unknown");

  // check if material kinematics is compatible to element kinematics
  SolidMaterial()->ValidKinematics(Inpar::STR::KinemType::nonlinearTotLag);

  // transverse shear locking
  linedef->extract_string("TSL", buffer);
  if (buffer == "dsg")
    dsg_shear_ = true;
  else if (buffer == "none")
    dsg_shear_ = false;
  else
    FOUR_C_THROW("unknown transverse shear locking method");

  // membrane locking
  linedef->extract_string("MEL", buffer);
  if (buffer == "dsg")
    dsg_membrane_ = true;
  else if (buffer == "none")
    dsg_membrane_ = false;
  else
    FOUR_C_THROW("unknown membrane locking method");

  // curvature thickness locking
  linedef->extract_string("CTL", buffer);
  if (buffer == "dsg")
    dsg_ctl_ = true;
  else if (buffer == "none")
    dsg_ctl_ = false;
  else
    FOUR_C_THROW("unknown curvature thickness locking method");

  // volumetric locking
  linedef->extract_string("VOL", buffer);
  if (buffer == "eas9")
    eas_ = true;
  else if (buffer == "none")
    eas_ = false;
  else
    FOUR_C_THROW("unknown volumetric locking method");

  setup_dsg();

  return true;
}

FOUR_C_NAMESPACE_CLOSE
