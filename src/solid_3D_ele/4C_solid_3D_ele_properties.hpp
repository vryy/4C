// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SOLID_3D_ELE_PROPERTIES_HPP
#define FOUR_C_SOLID_3D_ELE_PROPERTIES_HPP

#include "4C_config.hpp"

#include "4C_inpar_structure.hpp"
#include "4C_utils_exceptions.hpp"

#include <string>

FOUR_C_NAMESPACE_OPEN

namespace Core::Communication
{
  class PackBuffer;
  class UnpackBuffer;
}  // namespace Core::Communication

namespace Discret::Elements
{
  enum class ElementTechnology
  {
    none,
    fbar,
    eas_mild,
    eas_full,
    shell_ans,
    shell_eas,
    shell_eas_ans
  };

  std::string element_technology_string(const ElementTechnology ele_tech);

  enum class PrestressTechnology
  {
    none,
    mulf
  };

  /*!
   *  @brief struct for managing solid element properties
   */
  template <unsigned dim>
  struct SolidElementProperties;

  template <>
  struct SolidElementProperties<3>
  {
    //! kinematic type
    Inpar::Solid::KinemType kintype{Inpar::Solid::KinemType::vague};

    //! element technology (none, F-Bar, EAS full, EAS mild)
    ElementTechnology element_technology{ElementTechnology::none};

    //! specify prestress technology (none, MULF)
    PrestressTechnology prestress_technology{PrestressTechnology::none};
  };

}  // namespace Discret::Elements


FOUR_C_NAMESPACE_CLOSE

#endif
