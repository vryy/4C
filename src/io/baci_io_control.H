/*----------------------------------------------------------------------*/
/*! \file
 * \brief output control
\level 0
*/
/*----------------------------------------------------------------------*/

#ifndef BACI_IO_CONTROL_H
#define BACI_IO_CONTROL_H

#include "baci_discretization_fem_general_shape_function_type.H"

#include <Epetra_Comm.h>

#include <fstream>
#include <string>

extern "C"
{
#include "baci_io_legacy_types.h"
}

namespace IO
{
  /// control class to manage a control file for output
  class OutputControl
  {
   public:
    OutputControl(const Epetra_Comm& comm, std::string problemtype,
        CORE::FE::ShapeFunctionType type_of_spatial_approx, std::string inputfile,
        const std::string& outputname, int ndim, int restart, int filesteps,
        int create_controlfile);

    OutputControl(const Epetra_Comm& comm, std::string problemtype,
        CORE::FE::ShapeFunctionType type_of_spatial_approx, std::string inputfile,
        const std::string& restartname, std::string outputname, int ndim, int restart,
        int filesteps, int create_controlfile, bool adaptname = true);

    /// \brief copy constructor
    /** \param[in] ocontrol   Copy this object of same type
     *  \param[in] new_prefix Add the prefix to restart and file names (optional)
     *
     *  \author hiermeier \date 08/17 */
    OutputControl(const OutputControl& ocontrol, const char* new_prefix = nullptr);

    /// output prefix we write to
    /*!
      In case of restart this will be different from the read prefix.

      \note might contain path
     */
    std::string FileName() const { return filename_; }

    /**
     * @brief Return the file name prefix, i.e., the file name that is given to the baci call
     */
    std::string FileNameOnlyPrefix() const;

    /**
     * @brief Base output directory
     */
    std::string DirectoryName() const;

    /// original prefix as given
    /*!
      In case of restart this prefix specifies the control file we read.

      \note might contain path
     */
    std::string RestartName() const { return restartname_; }

    std::string NewOutputFileName() const { return filename_; }

    /// open control file
    std::fstream& ControlFile() { return controlfile_; }

    /// number of output steps per binary file
    int FileSteps() const { return filesteps_; }

    // modify the number of output steps per binary file
    // (necessary for the structural debugging option "OUTPUTEVERYITER")
    void SetFileSteps(int filesteps) { filesteps_ = filesteps; }

    /// input filename
    std::string InputFileName() const { return inputfile_; }

    int BinIO() const { return create_controlfile_; }

    /// overwrites result files
    void OverwriteResultFile(const CORE::FE::ShapeFunctionType& spatial_approx);
    /// creates new result files
    void NewResultFile(int numb_run, const CORE::FE::ShapeFunctionType& spatial_approx);
    /// creates new result files for the mlmc
    void NewResultFile(const std::string& name_appendix, int numb_run,
        const CORE::FE::ShapeFunctionType& spatial_approx);

    /// creates new result files
    void NewResultFile(std::string name, const CORE::FE::ShapeFunctionType& spatial_approx);

    /// return my processor ID
    inline int MyRank() const { return myrank_; };

   private:
    void WriteHeader(
        const std::string& control_file_name, const CORE::FE::ShapeFunctionType& spatial_approx);

    void InsertRestartBackReference(int restart, const std::string& outputname);

   private:
    std::string problemtype_;
    std::string inputfile_;  ///< input file name
    int ndim_;
    std::string filename_;  ///< prefix of outputfiles (might contain path)
    std::string restartname_;
    std::fstream controlfile_;
    int filesteps_;
    int create_controlfile_;
    int myrank_;
  };


  /// control class to manage a control file for input
  class InputControl
  {
   public:
    InputControl(const std::string& filename, const bool serial = false);
    InputControl(const std::string& filename, const Epetra_Comm& comm);
    ~InputControl();

    MAP* ControlFile() { return &table_; }

    std::string FileName() const { return filename_; }


   private:
    InputControl(const InputControl&);
    InputControl& operator=(const InputControl&);

    std::string filename_;
    MAP table_;
  };


  /// Controller for error files
  ///
  /// This object determines the names of error files #errname_, e.g. 'xxx2.err',
  /// and attaches file handles #errfile_ to them.
  ///
  /// \author bborn
  /// \date 11/08
  class ErrorFileControl
  {
   public:
    /// constructor
    ErrorFileControl(const Epetra_Comm& comm,  ///< communicator
        std::string outputname,                ///< output file name prefix
        int restart,                           ///< restart step
        int create_errorfiles                  // do we need to create errorfiles
    );

    /// destructor
    ~ErrorFileControl();

    /// the error file name
    std::string Name() { return errname_; }

    /// the error file handle
    FILE* Handle() { return errfile_; }

   private:
    /// output file name prefix
    std::string filename_;

    /// name of error file on processor
    std::string errname_;

    /// file handle to error file
    FILE* errfile_;
  };

  /// find position of restart number in filename (if existing):
  /// for "outname-5" will return position of the "-"
  /// returns std::string::npos if not found
  size_t RestartFinder(const std::string& filename);

}  // namespace IO

#endif
