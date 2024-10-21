// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_STRU_MULTI_MICROSTATIC_NPSUPPORT_HPP
#define FOUR_C_STRU_MULTI_MICROSTATIC_NPSUPPORT_HPP

#include "4C_config.hpp"

FOUR_C_NAMESPACE_OPEN

namespace MultiScale
{
  // std::endless loop for supporting procs in multi scale problems
  void np_support_drt();


}  // namespace MultiScale

FOUR_C_NAMESPACE_CLOSE

#endif
