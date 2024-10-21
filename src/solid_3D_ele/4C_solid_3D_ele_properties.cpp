// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_config.hpp"

#include "4C_solid_3D_ele_properties.hpp"

#include "4C_comm_parobject.hpp"
#include "4C_solid_scatra_3D_ele_factory.hpp"

FOUR_C_NAMESPACE_OPEN

void Discret::ELEMENTS::add_to_pack(Core::Communication::PackBuffer& data,
    const Discret::ELEMENTS::SolidElementProperties& properties)
{
  add_to_pack(data, static_cast<int>(properties.kintype));
  add_to_pack(data, static_cast<int>(properties.element_technology));
  add_to_pack(data, static_cast<int>(properties.prestress_technology));
}

void Discret::ELEMENTS::extract_from_pack(Core::Communication::UnpackBuffer& buffer,
    Discret::ELEMENTS::SolidElementProperties& properties)
{
  properties.kintype = static_cast<Inpar::Solid::KinemType>(extract_int(buffer));
  properties.element_technology = static_cast<ElementTechnology>(extract_int(buffer));
  properties.prestress_technology = static_cast<PrestressTechnology>(extract_int(buffer));
}

FOUR_C_NAMESPACE_CLOSE
