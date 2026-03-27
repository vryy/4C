// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FS3I_BIOFILM_FSI_INPUT_HPP
#define FOUR_C_FS3I_BIOFILM_FSI_INPUT_HPP

#include "4C_config.hpp"

#include "4C_io_input_spec.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

// forward declaration
namespace Core::Conditions
{
  class ConditionDefinition;
}

namespace BioFilm
{
  /// biofilm parameters
  Core::IO::InputSpec valid_parameters();

  /// set specific biofilm conditions
  void set_valid_conditions(std::vector<Core::Conditions::ConditionDefinition>& condlist);
}  // namespace BioFilm


/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
