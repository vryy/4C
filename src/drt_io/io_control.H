/*----------------------------------------------------------------------*/
/*!
 * \file io_control.H
\brief output control

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef IO_CONTROL_H
#define IO_CONTROL_H

#include <fstream>
#include <string>

#include <Epetra_Comm.h>

extern "C" {

#include "../pss_full/pss_types.h"

}

namespace IO
{

  /// control class to manage a control file for output
  class OutputControl
  {
  public:
    OutputControl(const Epetra_Comm& comm,
                  std::string problemtype,
                  std::string type_of_spatial_approx,
                  std::string inputfile,
                  std::string filename,
                  int ndim,
                  int restart,
                  int filesteps,
                  int create_controlfile);

    OutputControl(const Epetra_Comm& comm,
                  std::string problemtype,
                  std::string type_of_spatial_approx,
                  std::string inputfile,
                  std::string restartfilename,
                  std::string outputfilename,
                  int ndim,
                  int restart,
                  int filesteps,
                  int create_controlfile,
                  bool adaptname = true);

    /// output prefix we write to
    /*!
      In case of restart this will be different from the read prefix.

      \note might contain path
     */
    std::string FileName() const { return filename_; }

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

    /// input filename
    std::string InputFileName() const { return inputfile_;}

    int BinIO() const { return create_controlfile_; }

    /// overwrites result files for the inverse analysis
    void OverwriteResultFile();
    /// creates new result files for the inverse analysis
    void NewResultFile(int numb_run);
    /// creates new result files for the mlmc
    void NewResultFile(string name_appendix, int numb_run);

  private:
    std::string problemtype_;
    std::string inputfile_;  ///< input file name
    int ndim_;
    std::string filename_;  ///< prefix of outputfiles (might contain path)
    std::string restartname_;
    std::fstream controlfile_;
    int filesteps_;
    int create_controlfile_;
  };


  /// control class to manage a control file for input
  class InputControl
  {
  public:
    InputControl(std::string filename, const bool serial=false);
    InputControl(std::string filename, const Epetra_Comm& comm);
    ~InputControl();

    MAP* ControlFile() { return &table_; }

    std::string FileName() const { return filename_; }


  private:

    InputControl(const InputControl& );
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
    ErrorFileControl(
      const Epetra_Comm& comm,  ///< communicator
      const std::string filename, ///< output file name prefix
      int restart,   ///< restart step
      int create_errorfiles       // do we need to create errorfiles
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

}

#endif
