// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fluid_ele_parameter_std.hpp"

FOUR_C_NAMESPACE_OPEN

Discret::Elements::FluidEleParameterStd* Discret::Elements::FluidEleParameterStd::instance(
    Core::Utils::SingletonAction action)
{
  static auto singleton_owner = Core::Utils::make_singleton_owner(
      []()
      {
        return std::unique_ptr<Discret::Elements::FluidEleParameterStd>(
            new Discret::Elements::FluidEleParameterStd());
      });

  return singleton_owner.instance(action);
}

//----------------------------------------------------------------------*/
//    constructor
//----------------------------------------------------------------------*/
Discret::Elements::FluidEleParameterStd::FluidEleParameterStd()
    : Discret::Elements::FluidEleParameter::FluidEleParameter()
{
}

FOUR_C_NAMESPACE_CLOSE
