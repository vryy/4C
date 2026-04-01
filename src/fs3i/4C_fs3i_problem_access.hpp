// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FS3I_PROBLEM_ACCESS_HPP
#define FOUR_C_FS3I_PROBLEM_ACCESS_HPP

#include "4C_config.hpp"

#include "4C_global_data.hpp"

FOUR_C_NAMESPACE_OPEN

namespace FS3I::Utils
{
  inline Global::Problem* problem_from_instance() { return Global::Problem::instance(); }

  inline const Teuchos::ParameterList& fs3i_dynamic_params_from_problem()
  {
    return problem_from_instance()->f_s3_i_dynamic_params();
  }
}  // namespace FS3I::Utils

FOUR_C_NAMESPACE_CLOSE

#endif
