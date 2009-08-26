#ifdef CCADISCRET

#ifndef DRT_ELEMENTREADER_H
#define DRT_ELEMENTREADER_H

#include "drt_inputreader.H"

namespace DRT
{
namespace INPUT
{

/*----------------------------------------------------------------------*/
/*!
  \brief helper class to read the elements of a discretization

  Together with NodeReader this class constitutes a (almost) parallel
  and efficient reading mechanism for discretizations from dat files.

  We face the following problem:

  - There are elements and nodes. One set of elements per
    discretization. One set of nodes with the nodes from all
    discretizations.

  - Elements and nodes have ids. These are unique but otherwise
    arbitrary.

  - We cannot afford to read all elements or all nodes on one
    processor.

  - Only processor 0 can actually read the (ascii) input file

  - We do not want to setup (that is read) elements more than once.

  The idea is to read blocks of elements and nodes and distribute them
  to different processors at first. Afterwards a reasonable
  distribution can be calculated and the discretizations can be
  redistributed. How this work is done is a mere
  technicality. However, we need to be able to use the discretization
  in a partially constructed state. In particular we need to read, add
  and distribute elements even if nodes not yet are available.

  \author u.kue
  \date 05/07
 */
/*----------------------------------------------------------------------*/
class ElementReader
{
  friend class NodeReader;

public:

  /// construct element reader for a given field that reads a given section
  /*!
    Create empty discretization and append it to given field.

    \param dis (i) the new discretization
    \param comm (i) our communicator
    \param sectionname (i) the section that contains the element lines
   */
  ElementReader(RCP<Discretization> dis, const DatFileReader& reader, string sectionname);

  /// give the discretization this reader fills
  RCP<Discretization> MyDis() const { return dis_; }

private:

  /// read elements, partition node graph
  void Partition();

  /// tell whether the given node belongs to us
  /*!
    \note This is based on the redundant nodes_ set and only available
    on processor 0.
   */
  bool HasNode(int nodeid) { return nodes_.find(nodeid)!=nodes_.end(); }

  /// finialize reading. FillComplete(false,false,false), that is, do not
  /// initialize elements. This is done later after reading boundary conditions.
  void Complete();

  /// discretization name
  string name_;

  /// the main dat file reader
  const DatFileReader& reader_;

  /// my comm
  RCP<Epetra_Comm> comm_;

  /// my section to read
  string sectionname_;

  /// all global node ids of a discretization on processor 0
  /*!
    This is a redundant set of all node numbers. But it is only valid
    on processor 0. We need it to easily figure out to which
    discretization a node belongs.
   */
  set<int> nodes_;

  /// my discretization
  RCP<DRT::Discretization> dis_;

  /// node row map
  RCP<Epetra_Map> rownodes_;

  /// node col map
  RCP<Epetra_Map> colnodes_;

  /// element row map
  RCP<Epetra_Map> roweles_;

  /// element col map
  RCP<Epetra_Map> coleles_;
};


}
}

#endif
#endif
