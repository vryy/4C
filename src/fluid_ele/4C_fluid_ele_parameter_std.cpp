// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fluid_ele_parameter_std.hpp"

FOUR_C_NAMESPACE_OPEN

Discret::ELEMENTS::FluidEleParameterStd* Discret::ELEMENTS::FluidEleParameterStd::instance(
    Core::Utils::SingletonAction action)
{
  static auto singleton_owner = Core::Utils::make_singleton_owner(
      []()
      {
        return std::unique_ptr<Discret::ELEMENTS::FluidEleParameterStd>(
            new Discret::ELEMENTS::FluidEleParameterStd());
      });

  return singleton_owner.instance(action);
}

//----------------------------------------------------------------------*/
//    constructor
//----------------------------------------------------------------------*/
Discret::ELEMENTS::FluidEleParameterStd::FluidEleParameterStd()
    : Discret::ELEMENTS::FluidEleParameter::FluidEleParameter()
{
}

FOUR_C_NAMESPACE_CLOSE
