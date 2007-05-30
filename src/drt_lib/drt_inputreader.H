/*----------------------------------------------------------------------*/
/*!
\file drt_inputreader.H

\brief Internal classes to read elements and nodes

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#ifndef DRT_INPUTREADER_H
#define DRT_INPUTREADER_H

#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <list>
#include <string>
#include <algorithm>

#include <Epetra_Comm.h>
#include <Teuchos_RefCountPtr.hpp>

#include "drt_discret.H"

using namespace std;
using namespace Teuchos;

namespace DRT
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
   */
  ElementReader(string name, int pos, RefCountPtr<Epetra_Comm> comm, string sectionname);

  /// give the discretization this reader fills
  RefCountPtr<DRT::Discretization> Discretization() const { return dis_; }

private:

  /// read elements, partition node graph
  void Partition();

  /// allreduce node row map to simplify node reading
  void CollapseRowMap();

  /// finialize reading. FillComplete.
  void Complete();

  /// discretization name
  string name_;

  /// my comm
  RefCountPtr<Epetra_Comm> comm_;

  /// my section to read
  string sectionname_;

  /// my discretization
  RefCountPtr<DRT::Discretization> dis_;

  /// node row map
  RefCountPtr<Epetra_Map> rownodes_;

  /// node col map
  RefCountPtr<Epetra_Map> colnodes_;

  /// allreduced node map
  RefCountPtr<Epetra_Map> collapsedrows_;

  /// element row map
  RefCountPtr<Epetra_Map> roweles_;

  /// element col map
  RefCountPtr<Epetra_Map> coleles_;
};


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
  NodeReader(RefCountPtr<Epetra_Comm> comm, string sectionname);

  /// add an element reader for each discretization
  void AddElementReader(RefCountPtr<ElementReader> er) { ereader_.push_back(er); }

  /// do the actual reading
  void Read();

private:

  /// find the discretization that contains the given node
  RefCountPtr<DRT::Discretization> FindDisNode(int nodeid);

  /// my comm
  RefCountPtr<Epetra_Comm> comm_;

  /// my section to read
  string sectionname_;

  /// my element readers
  vector<RefCountPtr<ElementReader> > ereader_;
};

}

#endif
#endif
#endif
