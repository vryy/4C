// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_INPAR_PORO_HPP
#define FOUR_C_INPAR_PORO_HPP

#include "4C_config.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Inpar::Poro
{
  //! poro element implementation type
  enum class PoroType
  {
    undefined,
    pressure_based,
    pressure_velocity_based
  };

}  // namespace Inpar::Poro

FOUR_C_NAMESPACE_CLOSE

#endif