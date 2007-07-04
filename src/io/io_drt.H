/*----------------------------------------------------------------------*/
/*!
\file io_drt.H

\brief output context of one discretization

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef IO_DRT_H
#define IO_DRT_H

#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include <Teuchos_RefCountPtr.hpp>
#include <Epetra_BlockMap.h>

#include <map>
#include <string>

#include "../drt_lib/drt_discret.H"
#include "io.h"
#include "hdf_reader.H"

#ifdef BINIO
#include <H5LT.h>
#endif

extern "C" {
#include "../headers/standardtypes.h"
#include "../pss_full/pss_table.h"
}

#include "../drt_lib/dstrc.H" /* needs standardtypes.h */

/// IO: input/output facility
namespace IO
{

/*!
  \class DiscretizationReader

  \brief base class of baci restart

  \author m.kue
  \date 04/07
 */
class DiscretizationReader
{
public:

  /// construct reader for a given discretization to read a particular time step
  explicit DiscretizationReader(RefCountPtr<DRT::Discretization> dis, int step);

  /// read into given vector
  void ReadVector(RefCountPtr<Epetra_Vector> vec, string name);

  /// read an integer value from the control file
  int ReadInt(string name);

  /// read a double value from the control file
  double ReadDouble(string name);

  /// read into the discretization given in the constructor
  void ReadMesh(int step);

private:

  /// find control file entry to given time step
  MAP* FindResultGroup(int step);

  /// find control file entry to given time step
  MAP* FindMeshGroup(int step);

  /// Open data files. You guess it.
  void OpenDataFiles(MAP* result_step);

  /// Open mesh files. You guess it.
  void OpenMeshFiles(MAP* result_step);

  //! my discretization
  RefCountPtr<DRT::Discretization> dis_;

  /// control file entry of this step
  MAP *restart_step_;

#ifdef BINIO
  RefCountPtr<HDFReader> reader_;
  RefCountPtr<HDFReader> meshreader_;
#endif
};


/*!
  \class DiscretizationWriter

  \brief The output context of a discretization

  Create an object of this class for every discretization those mesh
  and results you want to write. Data are written in parallel to
  processor local files. The first process additionally maintains the
  (plain text) control file that glues all result files together.

  \author m.kue
  \date 02/07
*/
class DiscretizationWriter
{
public:

  /// Writer constructor always needs a discretization
  explicit DiscretizationWriter(RefCountPtr<DRT::Discretization> dis, int probnum=0);

  /// cleanup, close hdf5 files
  ~DiscretizationWriter();


  //!@name Output methods
  //@{

  //! write result header to control file
  /*!
    You will want to call this once each time step _before_ the
    result data is written.
    \param step : current time step
    \param time : current absolute time
  */
  void NewStep(int step, double time);

  //! write a result vector
  /*!
    There will be an entry in the current result step in the control
    file that points to this vector

    \param name : control file entry name
    \param vec  : the result data vector
  */
  void WriteVector(string name, RefCountPtr<Epetra_Vector> vec);

  //! write new "field" group to control file including node and element chunks
  void WriteMesh(int step, double time);

  //@}

protected:
  string cfname_;
  FILE* cf_;
  int steps_per_file_;

private:

  //! open new mesh file
  void CreateMeshFile(int step);

  //! open new result file
  void CreateResultFile(int step);

  //! my discretization
  RefCountPtr<DRT::Discretization> dis_;

  int disnum_;
  int field_pos_;
  int step_;
  double time_;

#ifdef BINIO
  hid_t meshfile_;
  hid_t resultfile_;
  string meshfilename_;
  string resultfilename_;
  hid_t meshgroup_;
  hid_t resultgroup_;
#endif

  std::map<const Epetra_BlockMapData*, std::string> mapcache_;

  int resultfile_changed_;
  int meshfile_changed_;
};

}

#endif
#endif
#endif
