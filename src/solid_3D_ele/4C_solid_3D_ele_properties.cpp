/*! \file

\brief Properties of solid elements

\level 1
*/


#include "4C_config.hpp"

#include "4C_solid_3D_ele_properties.hpp"

#include "4C_comm_parobject.hpp"
#include "4C_solid_scatra_3D_ele_factory.hpp"

FOUR_C_NAMESPACE_OPEN

void Discret::ELEMENTS::add_to_pack(Core::Communication::PackBuffer& data,
    const Discret::ELEMENTS::SolidElementProperties& properties)
{
  Core::Communication::ParObject::add_to_pack(data, static_cast<int>(properties.kintype));
  Core::Communication::ParObject::add_to_pack(
      data, static_cast<int>(properties.element_technology));
  Core::Communication::ParObject::add_to_pack(
      data, static_cast<int>(properties.prestress_technology));
}

void Discret::ELEMENTS::ExtractFromPack(std::size_t& position, const std::vector<char>& data,
    Discret::ELEMENTS::SolidElementProperties& properties)
{
  properties.kintype = static_cast<Inpar::STR::KinemType>(
      Core::Communication::ParObject::ExtractInt(position, data));
  properties.element_technology =
      static_cast<ElementTechnology>(Core::Communication::ParObject::ExtractInt(position, data));
  properties.prestress_technology =
      static_cast<PrestressTechnology>(Core::Communication::ParObject::ExtractInt(position, data));
}

FOUR_C_NAMESPACE_CLOSE
