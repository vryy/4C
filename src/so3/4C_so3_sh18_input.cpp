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
bool Discret::ELEMENTS::SoSh18::read_element(const std::string& eletype, const std::string& distype,
    const Core::IO::InputParameterContainer& container)
{
  // read number of material model
  int material_id = container.get<int>("MAT");

  set_material(0, Mat::Factory(material_id));

  // set up of materials with GP data (e.g., history variables)

  Teuchos::RCP<Core::Mat::Material> mat = material();

  solid_material()->setup(NUMGPT_SOH18, container);

  // temporary variable for read-in
  std::string buffer;

  // read kinematic flag
  buffer = container.get<std::string>("KINEM");
  if (buffer == "linear")
  {
    // kintype_ = soh8_linear;
    FOUR_C_THROW("Only nonlinear kinematics for SO_SH8 implemented!");
  }
  else if (buffer == "nonlinear")
  {
    kintype_ = Inpar::Solid::KinemType::nonlinearTotLag;
  }
  else
    FOUR_C_THROW("Reading SO_HEX18 element failed KINEM unknown");

  // check if material kinematics is compatible to element kinematics
  solid_material()->valid_kinematics(Inpar::Solid::KinemType::nonlinearTotLag);

  // transverse shear locking
  buffer = container.get<std::string>("TSL");

  if (buffer == "dsg")
    dsg_shear_ = true;
  else if (buffer == "none")
    dsg_shear_ = false;
  else
    FOUR_C_THROW("unknown transverse shear locking method");

  // membrane locking
  buffer = container.get<std::string>("MEL");
  if (buffer == "dsg")
    dsg_membrane_ = true;
  else if (buffer == "none")
    dsg_membrane_ = false;
  else
    FOUR_C_THROW("unknown membrane locking method");

  // curvature thickness locking
  buffer = container.get<std::string>("CTL");
  if (buffer == "dsg")
    dsg_ctl_ = true;
  else if (buffer == "none")
    dsg_ctl_ = false;
  else
    FOUR_C_THROW("unknown curvature thickness locking method");

  // volumetric locking
  buffer = container.get<std::string>("VOL");

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
