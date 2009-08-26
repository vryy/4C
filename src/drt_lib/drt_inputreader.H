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
#include "../drt_nurbs_discret/drt_knotvector.H"
#include "../drt_nurbs_discret/drt_control_point.H"

using namespace std;
using namespace Teuchos;

namespace DRT
{

/// Helper classes to read and validate input files
namespace INPUT
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
  DatFileReader(string filename, RefCountPtr<Epetra_Comm> comm, int outflag=0, bool dumpinput=true);

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

  /// return the content of a section
  vector<const char*> Section(string name) const;

  /// read a node-design topology section
  void ReadDesign(const std::string& name, std::vector<std::vector<int> >& dobj_fenode) const;

  /*!
    \brief read the knotvector section (for isogeometric analysis)

    \param  dim            (in ): dimension, i.e. number of knotvectors in one patch
    \param  name           (in ): Name/type of discretisation
    \param  disknots       (out): node vector coordinates

  */
  void ReadKnots(
    const int                              dim     ,
    const string                           name    ,
    Teuchos::RCP<DRT::NURBS::Knotvector> & disknots
    ) const;

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


}
}

#endif
#endif
