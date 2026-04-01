// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later


#ifndef FOUR_C_REDUCED_LUNG_MAIN_HPP
#define FOUR_C_REDUCED_LUNG_MAIN_HPP

#include "4C_config.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Global
{
  class Problem;
}

namespace ReducedLung
{
  void reduced_lung_main(Global::Problem& problem);
  void reduced_lung_main();
}  // namespace ReducedLung

FOUR_C_NAMESPACE_CLOSE

#endif
