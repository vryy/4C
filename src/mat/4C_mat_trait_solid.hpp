// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_MAT_TRAIT_SOLID_HPP
#define FOUR_C_MAT_TRAIT_SOLID_HPP

#include "4C_config.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  namespace Trait
  {
    /*!
     * This interface should define all core functionality for solid materials. Currently this is
     * done in So3Material.
     */
    class Solid
    {
     public:
      virtual ~Solid() = default;
    };
  }  // namespace Trait
}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif
