// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_INPAR_SOLVER_HPP
#define FOUR_C_INPAR_SOLVER_HPP

#include "4C_config.hpp"

#include "4C_utils_parameter_list.fwd.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Inpar::SOLVER
{
  //! set the parameters for the linear solver
  void set_valid_parameters(Teuchos::ParameterList& list);

}  // namespace Inpar::SOLVER

FOUR_C_NAMESPACE_CLOSE

#endif
