/*----------------------------------------------------------------------*/
/*! \file

\brief Read node sections of dat files.

\level 0


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_IO_NODEREADER_HPP
#define FOUR_C_IO_NODEREADER_HPP

#include "baci_config.hpp"

#include "baci_io_elementreader.hpp"
#include "baci_io_inputreader.hpp"

#include <string>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace IO
{
  /**
   * Read all nodes that are defined in section @p node_section_name in the @p reader.
   * Nodes are added to the Discretization objects associated with the @p element_readers.
   * The @p max_node_id is tracked for consistency checks.
   */
  void ReadNodes(const INPUT::DatFileReader& reader, const std::string& node_section_name,
      std::vector<ElementReader>& element_readers, int& max_node_id);

}  // namespace IO

FOUR_C_NAMESPACE_CLOSE

#endif
