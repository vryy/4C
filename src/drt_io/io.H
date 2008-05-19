/*----------------------------------------------------------------------*/
/*!
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

#include <Teuchos_RCP.hpp>
#include <Epetra_Comm.h>
#include <Epetra_BlockMap.h>

#include <map>
#include <string>

#include "../drt_lib/drt_discret.H"
//#include "io.h"
#include "io_hdf.H"

#ifdef BINIO
#include <H5LT.h>
#endif

extern "C" {
#include "../headers/standardtypes.h"
#include "../pss_full/pss_table.h"
}

/// IO: input/output facility
namespace IO
{

  class InputControl;
  class OutputControl;

/*!
  \brief base class of baci restart

  \author m.kue
  \date 04/07
 */
class DiscretizationReader
{
public:

  /// construct reader for a given discretization to read a particular time step
  DiscretizationReader(Teuchos::RCP<DRT::Discretization> dis,
                       Teuchos::RCP<IO::InputControl> input,
                       int step);

  /// construct reader for a given discretization to read a particular time step
  /*!
    Use default control file from global problem.
   */
  DiscretizationReader(Teuchos::RCP<DRT::Discretization> dis,
                       int step);

  /// read into given vector
  void ReadVector(Teuchos::RCP<Epetra_Vector> vec, string name);

  /// read into given vector
  void ReadMultiVector(Teuchos::RCP<Epetra_MultiVector> vec, string name);

  /// read into given std::map<int, RCP<Epetra_SerialDenseMatrix> >
  void ReadSerialDenseMatrix(RefCountPtr<std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> > > mapdata,
                             string name);

  /// read an integer value from the control file
  int ReadInt(string name);

  /// read a double value from the control file
  double ReadDouble(string name);

  /// read into the discretization given in the constructor
  void ReadMesh(int step);

protected:

  /// find control file entry to given time step
  void FindResultGroup(int step, MAP* file);


private:

  /// find control file entry to given time step
  void FindMeshGroup(int step, MAP* file);

  /// find control file entry to given time step
  /*!
    The control file entry with the given caption those field and step match
    my discretization and step. From that we need a backward search to find
    the entry that links to the binary files that cover our entry.
   */
  void FindGroup(int step,
                 MAP* file,
                 const char* caption,
                 const char* filestring,
                 MAP*& result_info,
                 MAP*& file_info);

#ifdef BINIO

  /// Open data files.
  Teuchos::RCP<HDFReader> OpenFiles(const char* filestring, MAP* result_step);
#endif

  //! my discretization
  Teuchos::RCP<DRT::Discretization> dis_;

  /// my input control file
  Teuchos::RCP<IO::InputControl> input_;

  /// control file entry of this step
  MAP *restart_step_;

#ifdef BINIO
  Teuchos::RCP<HDFReader> reader_;
  Teuchos::RCP<HDFReader> meshreader_;
#endif
};


/*!
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

  enum VectorType {
    dofvector,
    nodevector,
    elementvector
  };

  /// Writer constructor always needs a discretization
  DiscretizationWriter(Teuchos::RCP<DRT::Discretization> dis,
                       Teuchos::RCP<OutputControl> output);

  /// Writer constructor always needs a discretization
  /*!
    Use default control file from global problem.
   */
  explicit DiscretizationWriter(Teuchos::RCP<DRT::Discretization> dis);

  /// cleanup, close hdf5 files
  virtual ~DiscretizationWriter();


  //!@name Output methods
  //@{

  //! write result header to control file
  /*!
    You will want to call this once each time step _before_ the
    result data is written.
    \param step : current time step
    \param time : current absolute time
  */
  virtual void NewStep(const int step, const double time);

  //! write a result double to constrol file
  /*!
    There will be an entry in the current result step in the control
    file that points to this vector

    \param name : control file entry name
    \param value  : the result data value

    \author tk
    \date 04/08
  */
  virtual void WriteDouble(const string name, const double value);

  //! write a result integer to constrol file
  /*!
    There will be an entry in the current result step in the control
    file that points to this vector

    \param name : control file entry name
    \param value  : the result data value

    \author tk
    \date 04/08
  */
  virtual void WriteInt(const string name, const int value);


  //! write a result vector
  /*!
    There will be an entry in the current result step in the control
    file that points to this vector

    \param name : control file entry name
    \param vec  : the result data vector
    \param nodevector : whether this vector has a row dof map (nodevector==true) or a row element map.
  */
  virtual void WriteVector(const string name, Teuchos::RCP<Epetra_MultiVector> vec, VectorType vt=dofvector);

  //! write a result vector
  /*!
    There will be an entry in the current result step in the control
    file that points to this vector

    \param name : control file entry name
    \param vec  : the result data vector
    \param elemap: element map of discretization
    \param nodevector : whether this vector has a row dof map (nodevector==true) or a row element map.
  */
  virtual void WriteVector(const string name, const std::vector<char>& vec, const Epetra_Map& elemap, VectorType vt=dofvector);

  //! write a condition identified by its name
  virtual void WriteCondition(const string condname) const;

  //! write new "field" group to control file including node and element chunks
  virtual void WriteMesh(const int step, const double time);

  //! write element data to file
  virtual void WriteElementData();

  //! reset resultfile_changed_ to -1 in case of inverse analysis
  void RevertResultfileChanged() { resultfile_changed_ = -1;}

  //@}

private:

  //! open new mesh file
  virtual void CreateMeshFile(const int step);

  //! open new result file
  virtual void CreateResultFile(const int step);

  /// test writer preconditions
  void Check();

  //! my discretization
  Teuchos::RCP<DRT::Discretization> dis_;

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

  /// cache to remember maps we have already written
  std::map<const Epetra_BlockMapData*, std::string> mapcache_;

  /// dummy stack to really save the maps we cache
  std::vector<Epetra_BlockMap> mapstack_;

  int resultfile_changed_;
  int meshfile_changed_;

  Teuchos::RCP<OutputControl> output_;
};

}

#endif
#endif
