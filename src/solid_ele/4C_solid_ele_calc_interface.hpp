// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SOLID_ELE_CALC_INTERFACE_HPP
#define FOUR_C_SOLID_ELE_CALC_INTERFACE_HPP


#include "4C_config.hpp"

#include "4C_structure_new_input.hpp"

FOUR_C_NAMESPACE_OPEN


namespace Discret::Elements
{
  struct StressIO
  {
    FourC::Solid::StressType type;
    std::vector<char>& mutable_data;
  };

  struct StrainIO
  {
    FourC::Solid::StrainType type;
    std::vector<char>& mutable_data;
  };

}  // namespace Discret::Elements
FOUR_C_NAMESPACE_CLOSE

#endif
