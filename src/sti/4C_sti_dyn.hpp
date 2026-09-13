// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_STI_DYN_HPP
#define FOUR_C_STI_DYN_HPP

#include "4C_config.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Global
{
  class Problem;
}  // namespace Global

/*!
 * \brief Entry point for simulations of scalar-thermo interaction problems
 *
 * \param problem      global problem providing parameters, discretizations and solvers
 * \param restartstep  time step for restart (0 if no restart is performed)
 */
FOUR_C_API(FOUR_C_CORE) void sti_dyn(Global::Problem& problem, const int& restartstep);

FOUR_C_NAMESPACE_CLOSE

#endif
