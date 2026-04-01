// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_TSI_PROBLEM_ACCESS_HPP
#define FOUR_C_TSI_PROBLEM_ACCESS_HPP

#include "4C_config.hpp"

#include "4C_global_data.hpp"

FOUR_C_NAMESPACE_OPEN

namespace TSI::Utils
{
  inline Global::Problem* problem_from_instance() { return Global::Problem::instance(); }

  inline const Teuchos::ParameterList& tsi_dynamic_params_from_problem()
  {
    return problem_from_instance()->tsi_dynamic_params();
  }

  inline const Teuchos::ParameterList& tsi_monolithic_dynamic_params_from_problem()
  {
    return tsi_dynamic_params_from_problem().sublist("MONOLITHIC");
  }
}  // namespace TSI::Utils

FOUR_C_NAMESPACE_CLOSE

#endif
