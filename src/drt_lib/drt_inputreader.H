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
/// reading, broadcasting and storing of dat file contents
/*!
  The normal baci job needs to read one dat file at the beginning of
  its execution. This file contains all input parameters as well as
  the meshes and so on and needs to be available on all
  processors. However, we cannot rely on a shared file system, so only
  on process (the first one) can do the actual reading. All others
  need to get the data via MPI communication.

  This class manages the reading and broadcasting of the input
  file. To do so processor 0 reads the input file line by line and
  stores all lines in an array. Afterwards these lines are
  communicated to all other processors, so each processor has an
  internal copy of the input file and subsequently parses it.  There
  is one instance of this class during the input phase of baci. When
  everything is set up properly this instance goes away and the
  internal copy of the input file vanishes.

  There are two more details to be considered: First we do not want to
  read all elements and nodes on each processor. So we skip reading
  these sections here. Elements and nodes are read in parallel by the
  ElementReader and NodeReader, respectively.

  The second point is that the DatFileReader has to interact with the
  old ccarat input system. That is once the input file is read you
  have to Activate() the reader, which publishes some internal data of
  the DatFileReader to the ccarat input system. After that the old fr*
  functions can be used to read as long as the DatFileReader object is
  alive. You must not read anymore once the active DatFileReader
  object has been destroyed!

  The point of having this reader is to allow for special problem
  types that need to read more that one dat file. (Scary!) This can be
  done be creating a new (additional) reader to read the additional
  file.

  \author u.kue
  \date 06/07
 */
/*----------------------------------------------------------------------*/
class DatFileReader
{
public:

  /// construct a reader for a given file
  DatFileReader(string filename, RefCountPtr<Epetra_Comm> comm, int outflag=0);

  /// Publish (some) internal data to the old ccarat reading system
  /*!
    Note that these links remain intact even when the reader goes
    away. So you have to be careful not to use the ccarat reading
    functions any more once the active reader is gone.
  */
  void Activate();

  /// return my inputfile name
  string MyInputfileName() const;

  /// return my output flag
  int MyOutputFlag() const;

  /// get file position of excluded section
  /*!
    The big sections (nodes and elements) are excluded from general
    read so that we can read them in parallel. This is done using the
    NodeReader and the ElementReader.
   */
  ifstream::pos_type ExcludedSectionPosition(string section) const { return excludepositions_.find(section)->second; }

  /// return my communicator
  Teuchos::RCP<Epetra_Comm> Comm() const { return comm_; }

  /// convert a parameter section in extended format in a parameter list
  bool ReadSection(string name, Teuchos::ParameterList& list);

  /// convert a parameter section in GiD format in a parameter list
  bool ReadGidSection(string name, Teuchos::ParameterList& list);

private:

  /// parse name string an get sublist of global input list
  Teuchos::ParameterList& FindSublist(string name, Teuchos::ParameterList& list);

  /// parse the value and add the key,value pair to the list
  void AddEntry(string key, string value, Teuchos::ParameterList& list);

  /// actual read dat file, store and broadcast general sections
  void ReadDat();

  /// dump input file contents to error file
  void DumpInput();

  /// input file name
  string filename_;

  /// my communicator
  RefCountPtr<Epetra_Comm> comm_;

  /// number of lines
  int numrows_;

  /// flag for output (default: output should be written)
  int outflag_;

  /// the whole input file, lines separated by 0
  vector<char> inputfile_;

  /// pointers to line beginnings
  vector<char*> lines_;

  /// file positions of skipped sections
  map<string,ifstream::pos_type> excludepositions_;

  /// section positions of all sections
  map<string,unsigned> positions_;
};


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
  ElementReader(RefCountPtr<Discretization> dis, const DRT::DatFileReader& reader, string sectionname);

  /// give the discretization this reader fills
  RefCountPtr<Discretization> MyDis() const { return dis_; }

private:

  /// read elements, partition node graph
  void Partition();

  /// tell whether the given node belongs to us
  /*!
    \note This is based on the redundant nodes_ set and only available
    on processor 0.
   */
  bool HasNode(int nodeid) { return nodes_.find(nodeid)!=nodes_.end(); }

  /// finialize reading. FillComplete.
  void Complete();

  /// discretization name
  string name_;

  /// the main dat file reader
  const DRT::DatFileReader& reader_;

  /// my comm
  RefCountPtr<Epetra_Comm> comm_;

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
  RefCountPtr<DRT::Discretization> dis_;

  /// node row map
  RefCountPtr<Epetra_Map> rownodes_;

  /// node col map
  RefCountPtr<Epetra_Map> colnodes_;

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
  NodeReader(const DRT::DatFileReader& reader, string sectionname);

  /// add an element reader for each discretization
  /*!
    Each discretization needs its own ElementReader. These readers
    have to be registered at the NodeReader.

    \param er (i) a reader of one discretization that uses (a fraction of) our nodes
   */
  void AddElementReader(RefCountPtr<ElementReader> er) { ereader_.push_back(er); }

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
  RefCountPtr<DRT::Discretization> FindDisNode(int nodeid);

  /// the main dat file reader
  const DRT::DatFileReader& reader_;

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
