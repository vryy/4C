/*----------------------------------------------------------------------*/
/*!
\file hdf_reader.H

\brief Helpers to read HDF5 based output.

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef HDF_READER_H
#define HDF_READER_H

#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifdef BINIO

#include <H5LT.h>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "../drt_lib/drt_discret.H"
#include "Teuchos_RefCountPtr.hpp"

extern "C" {
#include "../headers/standardtypes.h"
}
#include "../drt_lib/dstrc.H"

#ifdef PARALLEL
#include <Epetra_MpiComm.h>
#endif
#include <Epetra_SerialComm.h>


#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"

namespace IO
{

/*!
  \class HDFReader
  \brief Helper class that handles the HDF5 files while reading.

  After parallel output we need parallel input. Each processor writes
  its own file. On input therefore each processor needs to read some
  files (zero to many, depending on the number of processors that
  wrote and read.) Additionally each processor might use several files
  (for several time steps, as it happens with restart.) This class
  handles the basic HDF5 file access.

  \author m.kue
  \date 02/07
*/
class HDFReader
{
public:

  //! Construct with the directory where the files can be found.
  explicit HDFReader(string dir);
  ~HDFReader();

  //! Open a new set of input files.
  /*!
    With num_output_procs==1 this function opens the result data
    file with name basename. If num_output_procs>1 it opens the result
    files of all processors, by appending .p<proc_num> to the
    basename.
  */
  void Open(string basename,int num_output_procs);
  //!
  void Close();

  //! read the packed element data from the mesh files
  /*!
    Note: this function should only be called while the HDFReader reads
    mesh files
  */
  RefCountPtr<std::vector<char> > ReadElementData(int step, int new_proc_num=0, int my_id=0);

  //! read the packed nodal data from the mesh files
  /*!
    Note: this function should only be called while the HDFReader reads
    mesh files
  */
  RefCountPtr<std::vector<char> > ReadNodeData(int step, int new_proc_num=0, int my_id=0);

#if 1
/*!reads the packed periodic boundary condition data from the mesh files
 *!
  Note: this function should only be called when the HDFReader opened
  the mesh files
 *----------------------------------------------------------------------*/
RefCountPtr<std::vector<char> > ReadPeriodicBoundaryConditions (
  int step, int new_proc_num, int my_id);
#endif


  //! read an Epetra_Vector from the result files
  /*!
    Right now an Epetra_Vector has to be read along with its map. Thus
    we read an integer and a double array here.

    \param id_path      (in): hdf5 path to map array (from control file)
    \param value_path   (in): hdf5 path to value array (from control file)
    \param Comm         (in): the communicator
   */
  RefCountPtr<Epetra_Vector> ReadResultData(string id_path, string value_path, const Epetra_Comm& Comm);

private:

  //! reads the dataset 'path' in all the files in the range [start,end)
  /*!
    Here we finally loop all the files the local processor has to read.
    returns all the data in one vector. The data is assumed to by of type char.
   */
  RefCountPtr<std::vector<char> > ReadCharData(std::string path,
                                               int start=0, int end=-1);

  //! reads the dataset 'path' in all the files in the range [start,end)
  /*!
    Here we finally loop all the files the local processor has to read.
    returns all the data in one vector<int>
  */
  RefCountPtr<std::vector<int> > ReadIntData(std::string path,
                                             int start=0, int end=-1);

  //! reads the dataset 'path' in all the files in the range [start,end)
  /*!
    Here we finally loop all the files the local processor has to
    read.
  */
  RefCountPtr<std::vector<double> > ReadDoubleData(std::string path,
                                                   int start=0, int end=-1);

  //! Figure out which subset of files this process needs to read
  /*!
    In a parallel run each processor writes one file. If we are to
    read this set of files on with a different number of processors,
    we need to know how many files each one has to read. One file will
    always be read by one processor. The reader already knows how many
    files there are (how many processors wrote this set.)

    We assign a consecutive range of files to each processor. Linearly
    distributed. After the read we can redistribute. :)

    \param new_proc_num (in): numbers of processors that read
    \param my_id        (in): my rank. id. processor number.
    \param start       (out): first file I have to read
    \param end         (out): first file I do not read
   */
  void CalculateRange(int new_proc_num, int my_id, int& start, int& end);

  //! the names of the files opened here
  std::vector<std::string> filenames_;

  //! HDF5 ids of the open files
  std::vector<hid_t> files_;

  //! directory where to find the files
  string input_dir_;

  //! number of processors that wrote this set of files
  int num_output_proc_;
};
}

#endif
#endif
#endif
#endif // HDF_READER_H
