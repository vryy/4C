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
#ifdef CCADISCRET
#ifndef IO_CONTROL_H
#define IO_CONTROL_H

#include <fstream>
#include <string>

#include <Epetra_Comm.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

extern "C" {

#include "../headers/standardtypes.h"
#include "../pss_full/pss_table.h"

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
                  int filesteps);

    OutputControl(const Epetra_Comm& comm,
                  std::string problemtype,
                  std::string type_of_spatial_approx,
                  std::string inputfile,
                  std::string restartfilename,
                  std::string outputfilename,
                  int ndim,
                  int restart,
                  int filesteps);

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
    
    /// creates new result files for the inverse analysis
    void NewResultFile(int numb_run);

  private:
    std::string problemtype_;
    std::string inputfile_;
    int ndim_;
    std::string filename_;
    std::string restartname_;
    std::fstream controlfile_;
    int filesteps_;
  };


  /// control class to manage a control file for input
  class InputControl
  {
  public:
    InputControl(std::string filename, const bool serial=false);
    ~InputControl();

    MAP* ControlFile() { return &table_; }

    std::string FileName() const { return filename_; }

  private:

    InputControl(const InputControl& );
    InputControl& operator=(const InputControl&);

    std::string filename_;
    MAP table_;
  };

}

#endif
#endif
