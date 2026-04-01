// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SSI_PROBLEM_ACCESS_HPP
#define FOUR_C_SSI_PROBLEM_ACCESS_HPP

#include "4C_config.hpp"

#include "4C_global_data.hpp"

FOUR_C_NAMESPACE_OPEN

namespace SSI::Utils
{
  inline Global::Problem* problem_from_instance() { return Global::Problem::instance(); }

  inline const Teuchos::ParameterList& ssi_control_params_from_problem()
  {
    return problem_from_instance()->ssi_control_params();
  }
}  // namespace SSI::Utils

FOUR_C_NAMESPACE_CLOSE

#endif
