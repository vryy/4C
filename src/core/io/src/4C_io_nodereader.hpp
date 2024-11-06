// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_IO_NODEREADER_HPP
#define FOUR_C_IO_NODEREADER_HPP

#include "4C_config.hpp"

#include "4C_io_elementreader.hpp"

#include <string>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::IO
{
  /**
   * Read all nodes that are defined in section @p node_section_name in the @p input.
   * Nodes are added to the discretization objects associated with the @p element_readers.
   * The @p max_node_id is tracked for consistency checks.
   */
  void read_nodes(Core::IO::InputFile& input, const std::string& node_section_name,
      std::vector<ElementReader>& element_readers, int& max_node_id);

}  // namespace Core::IO

FOUR_C_NAMESPACE_CLOSE

#endif
