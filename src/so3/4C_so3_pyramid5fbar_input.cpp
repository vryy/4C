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
bool Discret::ELEMENTS::SoPyramid5fbar::read_element(const std::string& eletype,
    const std::string& distype, const Core::IO::InputParameterContainer& container)
{
  // read number of material model
  int material_id = container.get<int>("MAT");
  set_material(0, Mat::Factory(material_id));

  // set up of materials with GP data (e.g., history variables)
  solid_material()->setup(NUMGPT_SOP5, container);

  // read kinematic flag
  std::string kinem = container.get<std::string>("KINEM");

  if (kinem == "linear")
  {
    FOUR_C_THROW("Only nonlinear kinematics for SO_PYRAMID5FBAR implemented!");
  }
  else if (kinem == "nonlinear")
  {
    kintype_ = Inpar::Solid::KinemType::nonlinearTotLag;
  }
  else
    FOUR_C_THROW("Reading SO_PYRAMID5FBAR element failed KINEM unknown");

  // check if material kinematics is compatible to element kinematics
  solid_material()->valid_kinematics(kintype_);

  // Validate that materials doesn't use extended update call.
  if (solid_material()->uses_extended_update())
    FOUR_C_THROW("This element currently does not support the extended update call.");

  return true;
}

FOUR_C_NAMESPACE_CLOSE
