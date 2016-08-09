/*!----------------------------------------------------------------------*

\file io.H

\brief output context of one discretization

\level 1

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235

*----------------------------------------------------------------------*/

#ifndef IO_DRT_H
#define IO_DRT_H


#include "io_hdf.H"

#include <Teuchos_RCP.hpp>

#include <map>
#include <string>
#include <vector>

extern "C"
{
#include "../pss_full/pss_types.h" // typedef struct MAP etc.
}


class Epetra_Vector;
class Epetra_MultiVector;
class Epetra_SerialDenseMatrix;
class Epetra_Map;
class Epetra_BlockMap;
class Epetra_BlockMapData;
namespace DRT
{
  class Discretization;
}

/// IO: input/output facility
namespace IO
{

  class InputControl;
  class OutputControl;
  class HDFReader;

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
  void ReadVector(Teuchos::RCP<Epetra_Vector> vec, std::string name);

  /// read into given vector
  void ReadMultiVector(Teuchos::RCP<Epetra_MultiVector> vec, std::string name);

  /// read into given std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix> >
  void ReadSerialDenseMatrix(Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix> > > mapdata,
                             std::string name);

  /// read an integer value from the control file
  int ReadInt(std::string name);

  /// read a double value from the control file
  double ReadDouble(std::string name);

  /// read into the discretization given in the constructor
  void ReadMesh(int step);

  /// read nodes into the discretization given in the constructor
  void ReadNodesOnly(int step);

  /// Read the history data of elements and nodes from restart files
  void ReadHistoryData(int step);

  //! read a non discretisation based vector of doubles
  /*!
    This vector should have been written only by proc0.
    It is assumed that this is a 'small' vector which has to be present on all procs.
    It is read from proc0 again and then communicated to all present procs.
   */
  void ReadRedundantDoubleVector(Teuchos::RCP<std::vector<double> >& doublevec, const std::string name);

  //! read a non discretisation based vector of integers
  /*!
    This vector should have been written only by proc0.
    It is assumed that this is a 'small' vector which has to be present on all procs.
    It is read from proc0 again and then communicated to all present procs.
   */
  void ReadRedundantIntVector( Teuchos::RCP<std::vector<int> >& intvec, const std::string name);

  /// return number of procs which were used for restart output (read from control file)
  int GetNumOutputProc(int step);

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



  /// Open data files.
  Teuchos::RCP<HDFReader> OpenFiles(const char* filestring, MAP* result_step);


  //! my discretization
  Teuchos::RCP<DRT::Discretization> dis_;

  /// my input control file
  Teuchos::RCP<IO::InputControl> input_;

  /// control file entry of this step
  MAP *restart_step_;


  Teuchos::RCP<HDFReader> reader_;
  Teuchos::RCP<HDFReader> meshreader_;

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

  /*!
    Writer constructor always needs a discretization
    Use default control file from global problem.
   */
  DiscretizationWriter(
    Teuchos::RCP<DRT::Discretization> dis);

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

  //! write a result double to control file
  /*!
    There will be an entry in the current result step in the control
    file that points to this vector

    \param name : control file entry name
    \param value  : the result data value

    \author tk
    \date 04/08
  */
  virtual void WriteDouble(const std::string name, const double value);

  //! write a result integer to constrol file
  /*!
    There will be an entry in the current result step in the control
    file that points to this vector

    \param name : control file entry name
    \param value  : the result data value

    \author tk
    \date 04/08
  */
  virtual void WriteInt(const std::string name, const int value);


  //! write a result vector
  /*!
    There will be an entry in the current result step in the control
    file that points to this vector

    \param name : control file entry name
    \param vec  : the result data vector
    \param vt   : vector type
  */
  virtual void WriteVector(const std::string name, Teuchos::RCP<const Epetra_MultiVector> vec, VectorType vt=dofvector);

  //! write a result vector
  /*!
    There will be an entry in the current result step in the control
    file that points to this vector

    \param name : control file entry name
    \param vec  : the result data vector
    \param elemap: element map of discretization
    \param vt   : vector type
  */
  virtual void WriteVector(const std::string name, const std::vector<char>& vec, const Epetra_Map& elemap, VectorType vt=dofvector);

  //! write a condition identified by its name
  virtual void WriteCondition(const std::string condname) const;


  /*!
    \brief write a knotvector for a nurbs discretisation

    \author gammi
    \date 06/08
  */
  virtual void WriteKnotvector() const;

  //! write new mesh and result file next time it is possible
  void CreateNewResultAndMeshFile(){resultfile_changed_=-1;meshfile_changed_=-1;};
  //! write new "field" group to control file including node and element chunks
  virtual void WriteMesh(const int step, const double time);
  // for MLMC purposes do not write new meshfile but write name of base mesh file to controlfile
  virtual void WriteMesh(const int step, const double time, std::string name_base_file );
  // for particle simulations: write only nodes in new "field" group to control file
  virtual void ParticleOutput(const int step, const double time, const bool writerestart);

  //! write element data to file
  virtual void WriteElementData(bool writeowner);

  //! write node data to file
  virtual void WriteNodeData(bool writeowner);

  //! write a non discretisation based vector of doubles
  /*!
    Write this vector only from proc0. It is assumed that this is a 'small' vector
    which is present on all procs. It shall be read from proc0 again and then
    communicated to all present procs.
   */
  virtual void WriteRedundantDoubleVector(const std::string name,
               Teuchos::RCP<std::vector<double> > doublevec);

  //! write a non discretisation based vector of integers
  /*!
    Write this vector only from proc0. It is assumed that this is a 'small' vector
    which is present on all procs. It shall be read from proc0 again and then
    communicated to all present procs.
   */
  virtual void WriteRedundantIntVector(const std::string name,
               Teuchos::RCP<std::vector<int> > setint );


  /// overwrite result files for the inverse analysis
  void OverwriteResultFile();

  /// creating new result files for the inverse analysis
  void NewResultFile(int numb_run);

  /// creating new result files for the mlmc
  void NewResultFile(std::string name_appendix, int numb_run);

  //@}

  //!@name Data management
  //@{

  /// clear all stored map data
  void ClearMapCache();

  //@}

  /// get output control
  Teuchos::RCP<OutputControl> Output() const { return output_;}

  /// set output control
  void SetOutput(Teuchos::RCP<OutputControl> output);


protected:

  //! open new mesh file
  virtual void CreateMeshFile(const int step);

  //! open new result file
  virtual void CreateResultFile(const int step);

  //! my discretization
  Teuchos::RCP<DRT::Discretization> dis_;

  int step_;
  double time_;


  hid_t meshfile_;
  hid_t resultfile_;
  std::string meshfilename_;
  std::string resultfilename_;
  hid_t meshgroup_;
  hid_t resultgroup_;


  /// cache to remember maps we have already written
  std::map<const Epetra_BlockMapData*, std::string> mapcache_;

  /// dummy stack to really save the maps we cache
  std::vector<Epetra_BlockMap> mapstack_;

  int resultfile_changed_;
  int meshfile_changed_;

  Teuchos::RCP<OutputControl> output_;

  // do we want binary output
  bool binio_;

};

/*!
  \brief The output context of a meshfree discretization

  \author Keijo Nissen
  \date Jan14
*/
class MeshfreeDiscretizationWriter : public DiscretizationWriter
{
public:

  /*!
    Writer constructor always needs a discretization
    Use default control file from global problem.
   */
  MeshfreeDiscretizationWriter(Teuchos::RCP<DRT::Discretization> dis)
    : DiscretizationWriter(dis)
  {};

  //! write new "field" group to control file including node and element chunks
  virtual void WriteMesh(const int step, const double time);
};


}

#endif
