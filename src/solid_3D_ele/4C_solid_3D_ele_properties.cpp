// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_config.hpp"

#include "4C_solid_3D_ele_properties.hpp"

#include "4C_solid_scatra_3D_ele_factory.hpp"
#include "4C_utils_enum.hpp"

FOUR_C_NAMESPACE_OPEN

std::string Discret::Elements::element_technology_string(const ElementTechnology ele_tech)
{
  switch (ele_tech)
  {
    case ElementTechnology::none:
      return "none";
    case ElementTechnology::fbar:
      return "fbar";
    case ElementTechnology::eas_mild:
      return "eas_mild";
    case ElementTechnology::eas_full:
      return "eas_full";
    case ElementTechnology::shell_ans:
      return "shell_ans";
    case ElementTechnology::shell_eas:
      return "shell_eas";
    case ElementTechnology::shell_eas_ans:
      return "shell_eas_ans";
  }

  FOUR_C_THROW("Unknown element technology {}", ele_tech);
}

FOUR_C_NAMESPACE_CLOSE
