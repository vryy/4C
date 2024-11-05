// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_INPAR_VALIDCONTACTCONSTITUTIVELAW_HPP
#define FOUR_C_INPAR_VALIDCONTACTCONSTITUTIVELAW_HPP

/*----------------------------------------------------------------------*/
#include "4C_config.hpp"

#include <iostream>
#include <memory>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace CONTACT
{
  namespace CONSTITUTIVELAW
  {
    class LawDefinition;
  }
}  // namespace CONTACT

namespace Input
{
  /// construct list with all contact constitutive laws and documentation
  std::shared_ptr<std::vector<std::shared_ptr<CONTACT::CONSTITUTIVELAW::LawDefinition>>>
  valid_contact_constitutive_laws();

  /** \brief print all known contact constitutive law sections without contents
   *
   * \param[in] contactconstitutivlawlist list of contact constitutive law definitions
   */
  void print_empty_contact_constitutive_law_definitions(std::ostream& stream,
      std::vector<std::shared_ptr<CONTACT::CONSTITUTIVELAW::LawDefinition>>&
          contactconstitutivlawlist);

}  // namespace Input

/// print empty contact constitutivelaw sections
void print_contact_constitutive_law_dat_header();


FOUR_C_NAMESPACE_CLOSE

#endif
