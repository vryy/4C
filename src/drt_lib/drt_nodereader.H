

#ifndef DRT_NODEREADER_H
#define DRT_NODEREADER_H

#include "drt_inputreader.H"
#include "drt_elementreader.H"

namespace DRT
{
namespace INPUT
{

/*----------------------------------------------------------------------*/
/*!
  \brief helper class to read all nodes shared by many discretizations

  The interface class to the reading subsystem.
  \author u.kue
  \date 05/07
 */
/*----------------------------------------------------------------------*/
class NodeReader
{
public:

  /// construct a reader that reads a given section
  NodeReader(const DatFileReader& reader, string sectionname);

  /// add an element reader for each discretization
  /*!
    Each discretization needs its own ElementReader. These readers
    have to be registered at the NodeReader.

    \param er (i) a reader of one discretization that uses (a fraction of) our nodes
   */
  void AddElementReader(Teuchos::RCP<ElementReader> er) { ereader_.push_back(er); }

  /// do the actual reading
  /*!
    This method contains the whole machinery. The reading consists of
    three steps:

    - Reading and distributing elements using each registered
      ElementReader. This includes creating the connectivity graph,
      partitioning, building the node row and column maps.

    - Reading and distributing all nodes. Each node gets assigned to
      its discretization.

    - Finalizing the discretizations.

    Actually most of the work gets done by the ElementReader. The
    reading of both elements and nodes happens in blocks on processor
    0. After each block read the discretizations are redistributed.

   */
  void Read();

private:

  /// find the discretization that contains the given node
  vector<Teuchos::RCP<DRT::Discretization> > FindDisNode(int nodeid);

  /// the main dat file reader
  const DatFileReader& reader_;

  /// my comm
  Teuchos::RCP<Epetra_Comm> comm_;

  /// my section to read
  string sectionname_;

  /// my element readers
  std::vector<Teuchos::RCP<ElementReader> > ereader_;
};

}
}

#endif
