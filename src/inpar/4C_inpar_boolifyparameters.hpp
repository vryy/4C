// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_INPAR_BOOLIFYPARAMETERS_HPP
#define FOUR_C_INPAR_BOOLIFYPARAMETERS_HPP

#include "4C_config.hpp"

#include "4C_utils_parameter_list.fwd.hpp"

FOUR_C_NAMESPACE_OPEN


namespace Input
{
  /// Auxiliar routine to boolify integral Yes/No data
  ///
  /// Parameters consisting integral values for Yes/No tuples
  /// are removed and replaced by a bool parameter holding true/false
  /// respectively.
  void boolify_valid_input_parameters(
      Teuchos::ParameterList& list  ///< the valid input parameter list
  );

}  // namespace Input

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
