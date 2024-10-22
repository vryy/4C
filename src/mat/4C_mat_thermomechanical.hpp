// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_MAT_THERMOMECHANICAL_HPP
#define FOUR_C_MAT_THERMOMECHANICAL_HPP

#include "4C_config.hpp"

#include "4C_mat_so3_material.hpp"
#include "4C_mat_trait_thermo_solid.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  /*!
   * Interface for all thermo-mechanical materials
   *
   * @note this interface  is inspired by the way materials are supposed to work in 4C, not all
   * materials are moved to this new interface
   */
  class ThermoMechanicalMaterial : public So3Material, public Trait::ThermoSolid
  {
  };

}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif
