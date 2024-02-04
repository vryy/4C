/*----------------------------------------------------------------------*/
/*! \file

\brief Read node sections of dat files.

\level 0


*/
/*----------------------------------------------------------------------*/

#ifndef BACI_LIB_NODEREADER_H
#define BACI_LIB_NODEREADER_H

#include "baci_config.H"

#include "baci_io_inputreader.H"
#include "baci_lib_elementreader.H"

#include <string>
#include <vector>

BACI_NAMESPACE_OPEN

namespace INPUT
{
  /**
   * Read all nodes that are defined in section @p node_section_name in the @p reader.
   * Nodes are added to the Discretization objects associated with the @p element_readers.
   * The @p max_node_id is tracked for consistency checks.
   */
  void ReadNodes(const DatFileReader& reader, const std::string& node_section_name,
      std::vector<ElementReader>& element_readers, int& max_node_id);

}  // namespace INPUT

BACI_NAMESPACE_CLOSE

#endif
