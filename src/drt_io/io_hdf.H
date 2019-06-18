/*----------------------------------------------------------------------*/
/*!
 * \brief Helpers to read HDF5 based output.

\level 1

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235

*----------------------------------------------------------------------*/

#ifndef IO_HDF_H
#define IO_HDF_H


#ifndef H5GCreate_vers
#define H5Gcreate_vers 1
#endif
#ifndef H5Dopen_vers
#define H5Dopen_vers 1
#endif

#include <hdf5.h>

#if H5_VERS_MINOR == 6
#include <H5LT.h>
#elif H5_VERS_MINOR == 8 || H5_VERS_MINOR == 10
#include <hdf5_hl.h>
#else
#error Unsupported hdf5 version
#endif

#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

#include "Teuchos_RCP.hpp"

#ifdef PARALLEL
#include <Epetra_MpiComm.h>
#endif
#include <Epetra_SerialComm.h>


#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>

namespace IO
{
  /*!
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
    explicit HDFReader(std::string dir);
    ~HDFReader();

    //! Open a new set of input files.
    /*!
      With num_output_procs==1 this function opens the result data
      file with name basename. If num_output_procs>1 it opens the result
      files of all processors, by appending .p<proc_num> to the
      basename.
    */
    void Open(std::string basename, int num_output_procs, int new_proc_num, int my_id);
    //!
    void Close();

    /*!
     * \brief read the packed element data from the mesh files
     * \note  this function should only be called while the HDFReader reads
     *        mesh files
     */
    Teuchos::RCP<std::vector<char>> ReadElementData(int step, int new_proc_num, int my_id) const;

    //!
    /*!
     * \brief read the packed nodal data from the mesh files
     * \note: this function should only be called while the HDFReader reads
     *         mesh files
     */
    Teuchos::RCP<std::vector<char>> ReadNodeData(int step, int new_proc_num, int my_id) const;

    /*!
     * \brief reads the packed periodic boundary condition data from the mesh files
     *
     *
     * \note this function should only be called when the HDFReader opened
     *       the mesh files
     */
    Teuchos::RCP<std::vector<char>> ReadCondition(
        const int step, const int new_proc_num, const int my_id, const std::string condname) const;

    /*!
    //  \brief reads the packed knotvector data from the mesh files
    //
    //    \param   step (i)
    //
    //      \return  The whole knotvector data in a char vector
    //      */
    Teuchos::RCP<std::vector<char>> ReadKnotvector(const int step) const;


    //! read an Epetra_MultiVector from the result files
    /*!
      Right now an Epetra_MultiVector has to be read along with its map. Thus
      we read an integer and a double array here.

      \note If columns==1, we create an Epetra_Vector.

      \param id_path      (in): hdf5 path to map array (from control file)
      \param value_path   (in): hdf5 path to value array (from control file)
      \param columns      (in): number of vector columns
      \param Comm         (in): the communicator
     */
    Teuchos::RCP<Epetra_MultiVector> ReadResultData(
        std::string id_path, std::string value_path, int columns, const Epetra_Comm& Comm) const;

    //! read a std::vector<char> from the result files
    /*!
      Right now a std::vector<char> has to be read along with its map. Thus
      we read an integer and a double array here.

      \param id_path      (in): hdf5 path to map array (from control file)
      \param value_path   (in): hdf5 path to value array (from control file)
      \param columns      (in): number of vector columns
      \param Comm         (in): the communicator
      \param elemap      (out): element map
     */
    Teuchos::RCP<std::vector<char>> ReadResultDataVecChar(std::string id_path,
        std::string value_path, int columns, const Epetra_Comm& Comm,
        Teuchos::RCP<Epetra_Map>& elemap) const;

    Teuchos::RCP<std::vector<char>> ReadCharVector(
        std::string value_path, const Epetra_Comm& Comm) const;

    Teuchos::RCP<std::vector<double>> ReadDoubleVector(std::string path) const
    {
      std::vector<int> length;
      Teuchos::RCP<std::vector<double>> values = ReadDoubleData(path, 0, 1, length);
      return values;
    };

    Teuchos::RCP<std::vector<int>> ReadIntVector(std::string path) const
    {
      std::vector<int> length;
      Teuchos::RCP<std::vector<int>> values = ReadIntData(path, 0, 1);
      return values;
    };

   private:
    /// reads the dataset 'path' in all the files in the range [start,end)
    /*!
      Here we finally loop all the files the local processor has to read.
      returns all the data in one vector. The data is assumed to by of type char.
     */
    Teuchos::RCP<std::vector<char>> ReadCharData(std::string path, int start, int end) const;

    /// reads the dataset 'path' in all the files in the range [start,end)
    /*!
      Here we finally loop all the files the local processor has to read.
      returns all the data in one vector<int>
    */
    Teuchos::RCP<std::vector<int>> ReadIntData(std::string path, int start, int end) const;

    /// reads the dataset 'path' in all the files in the range [start,end)
    /*!
      Here we finally loop all the files the local processor has to
      read.
    */
    Teuchos::RCP<std::vector<double>> ReadDoubleData(
        std::string path, int start, int end, std::vector<int>& lengths) const;

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
    void CalculateRange(int new_proc_num, int my_id, int& start, int& end) const;

    //! the names of the files opened here
    std::vector<std::string> filenames_;

    //! HDF5 ids of the open files
    std::vector<hid_t> files_;

    //! directory where to find the files
    std::string input_dir_;

    //! number of processors that wrote this set of files
    int num_output_proc_;

    //! file access property list for HDF5 files
    hid_t H5Plist_;
  };
}  // namespace IO


#endif
