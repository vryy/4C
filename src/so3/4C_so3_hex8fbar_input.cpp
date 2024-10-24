// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_io_linedefinition.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_so3_hex8fbar.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool Discret::Elements::SoHex8fbar::read_element(const std::string& eletype,
    const std::string& distype, const Core::IO::InputParameterContainer& container)
{
  // read number of material model
  int material_id = container.get<int>("MAT");
  set_material(0, Mat::factory(material_id));

  // set up of materials with GP data (e.g., history variables)
  solid_material()->setup(NUMGPT_SOH8, container);

  // read kinematic flag
  std::string kinem = container.get<std::string>("KINEM");
  if (kinem == "linear")
  {
    FOUR_C_THROW("Only nonlinear kinematics for SO_HEX8FBAR implemented!");
  }
  else if (kinem == "nonlinear")
  {
    kintype_ = Inpar::Solid::KinemType::nonlinearTotLag;
  }
  else
    FOUR_C_THROW("Reading SO_HEX8FBAR element failed KINEM unknown");

  // check if material kinematics is compatible to element kinematics
  solid_material()->valid_kinematics(kintype_);

  return true;
}

FOUR_C_NAMESPACE_CLOSE
