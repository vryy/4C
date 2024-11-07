// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_INPAR_VALIDPARAMETERS_HPP
#define FOUR_C_INPAR_VALIDPARAMETERS_HPP

#include "4C_config.hpp"

#include "4C_utils_parameter_list.fwd.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

#include <iostream>
#include <memory>
#include <string>

FOUR_C_NAMESPACE_OPEN

namespace Core::IO
{
  class Pstream;
}


namespace Input
{
  /**
   * Construct a `Teuchos::ParameterList` with all parameters and their documentation.
   *
   * @return A `std::shared_ptr` to a constant `Teuchos::ParameterList` containing all valid
   * parameters.
   */
  std::shared_ptr<const Teuchos::ParameterList> valid_parameters();

}  // namespace Input


/*! print list of valid parameters with documentation */
void print_valid_parameters();

/*! print help message */
void print_help_message();

/*! print flag sections of dat file with default flags */
void print_default_dat_header();


FOUR_C_NAMESPACE_CLOSE

#endif
