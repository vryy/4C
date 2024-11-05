// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_INPAR_VALIDMATERIALS_HPP
#define FOUR_C_INPAR_VALIDMATERIALS_HPP

/*----------------------------------------------------------------------*/
/* headers */
#include "4C_config.hpp"

#include <Teuchos_Array.hpp>

#include <iostream>
#include <memory>
#include <string>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  class MaterialDefinition;
}

namespace Input
{
  /// construct list with all materials and documentation
  std::shared_ptr<std::vector<std::shared_ptr<Mat::MaterialDefinition>>> valid_materials();

  /// print all known material sections without contents
  void print_empty_material_definitions(
      std::ostream& stream, std::vector<std::shared_ptr<Mat::MaterialDefinition>>& matlist);
}  // namespace Input

/// print empty material sections
void print_material_dat_header();


FOUR_C_NAMESPACE_CLOSE

#endif
