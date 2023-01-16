/*----------------------------------------------------------------------*/
/*! \file

\brief Read node sections of dat files.

\level 0

*/
/*----------------------------------------------------------------------*/

#include "lib_nodereader.H"
#include "lib_standardtypes_cpp.H"
#include "lib_elementdefinition.H"
#include "lib_globalproblem.H"
#include "rebalance_utils.H"
#include "lib_utils_factory.H"
#include "lib_utils_parallel.H"

#include <Epetra_Time.h>

namespace DRT::INPUT
{
  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  NodeReader::NodeReader(const DRT::INPUT::DatFileReader& reader, std::string sectionname)
      : reader_(reader), comm_(reader.Comm()), sectionname_(sectionname)
  {
  }

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void NodeReader::Read()
  {
    // TODO: Implement node reading
  }

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void NodeReader::ThrowIfNotEnoughNodes(int max_node_id) const
  {
    int local_max_node_id = max_node_id;
    comm_->MaxAll(&local_max_node_id, &max_node_id, 1);

    if ((max_node_id < comm_->NumProc()) && (reader_.ExcludedSectionLength(sectionname_) != 0))
      dserror("Bad idea: Simulation with %d procs for problem with %d nodes", comm_->NumProc(),
          max_node_id);
  }
}  // namespace DRT::INPUT