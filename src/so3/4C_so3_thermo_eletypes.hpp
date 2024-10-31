// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SO3_THERMO_ELETYPES_HPP
#define FOUR_C_SO3_THERMO_ELETYPES_HPP

#include "4C_config.hpp"

#include "4C_so3_hex8.hpp"


FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Discret
{
  namespace Elements
  {
    /*----------------------------------------------------------------------*
     * HEX8 element
     *----------------------------------------------------------------------*/
    class SoHex8ThermoType : public SoHex8Type
    {
     public:
      std::string name() const override { return "So_hex8ThermoType"; }

      static SoHex8ThermoType& instance();

      Core::Communication::ParObject* create(Core::Communication::UnpackBuffer& buffer) override;

      Teuchos::RCP<Core::Elements::Element> create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> create(const int id, const int owner) override;

      int initialize(Core::FE::Discretization& dis) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static SoHex8ThermoType instance_;

      std::string get_element_type_string() const { return "SOLIDH8THERMO"; }
    };  // class So_hex8ThermoType

  }  // namespace Elements

}  // namespace Discret


/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
