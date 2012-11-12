/*----------------------------------------------------------------------*/
/*!
\file io_pstream.H

\brief A substitute for STL cout for parallel and complex output schemes.

<pre>
Maintainer: Karl-Robert Wichmann
            wichmann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef IO_PSTREAM_H
#define IO_PSTREAM_H

/*----------------------------------------------------------------------*/
/* headers */
#include "../drt_lib/drt_dserror.H"

#include <Teuchos_RCP.hpp>
#include <Epetra_Comm.h>

#include <iostream>
#include <fstream>
#include <sstream>


/*----------------------------------------------------------------------*/
/* namespace */
namespace IO
{
  /*====================================================================*/
  /*!
   * \brief This object allows to write to std::cout in a very fancy way
   *        The output can be regulated by various parameters, in particular
   *        there are:
   *        o writing to screen
   *        o writing to files
   *        o modifying the output for use with nested parallelism
   *        o limit output to certain procs
   *
   * \author wichmann & hammerl \date 11/12
   */
  class Pstream
  {

  public :
    /// This empty constructor is called when the global object is instantiated.
    Pstream();

    /// configure the output. Must be called before IO::cout can be used. Is
    /// currently called in the global problem using the params specified in the
    /// IO section of the .dat-file.
    void setup(
      const bool writetoscreen,       ///< bool whether output is written to screen
      const bool writetofile,         ///< bool whether output is written to file
      const bool prefixgroupID,       ///< bool whether group ID is prefixed in each line
      Teuchos::RCP<Epetra_Comm> comm, ///< MPI communicator
      const int targetpid,            ///< target processor ID from which to print
      const int groupID,              ///< the ID
      const std::string fileprefix    ///< prefix for the output file
    );

    /// must be called to close open file handles and resets the IO::cout object to a pristine state.
    /// This is currently called by the destructor of the global problem singleton.
    void close();

    /// writes the buffer to screen
    void flush();

    /// This handles the actual printing to screen/writing to a file.
    template <typename CharT>
    Pstream& stream(CharT s)  ///< text to be added
    {
      if (not is_initialized_)
        dserror("Setup the output before you use it!");

      // are we on a proc that is allowed to output
      if (OnPid())
      {
        // store formatting in this semi persistent buffer
        fmt_tmp_ << s;
        if (fmt_tmp_.str().size() <= 0)
          return *this;

        // actual content goes into the string str
        std::string str = fmt_tmp_.str();
        fmt_tmp_.str(std::string());

        // write to file
        if (writetofile_)
        {
    #ifdef DEBUG
          if (!outfile_)
            dserror("outputfile does not exist - file handle lost");
    #endif
          (*outfile_) << str;
        }

        // write to screen
        if(writetoscreen_)
        {
          size_t oldpos = 0;
          size_t pos    = 0;
          while (true)
          {
            pos = str.find('\n',oldpos); // do include the carriage return itself

            if(pos == std::string::npos)
            {
              buffer_ << str.substr(oldpos, str.size()-oldpos);
              break;
            }

            buffer_ << str.substr(oldpos ,pos-oldpos+1);

            oldpos = pos+1;

            // we are using a hard coded std::cout here, as writing to screen means std::cout by definition
            std::cout << buffer_.str();
            std::flush(std::cout);
            buffer_.str(std::string());

            if (prefixgroupID_)
              buffer_ << groupID_ << ": ";
          }
        }
      }

      return *this;
    }

    std::stringstream& cout_replacement()
    {
      return fmt_tmp_;
    }

  private :
    /// Return whether this is a target processor for output
    bool OnPid();

    /// Shelter copy constructor
    Pstream(const Pstream& old);

    /// signal whether setup has been called
    bool is_initialized_;

    /// MPI communicator
    Teuchos::RCP<Epetra_Comm> comm_;

    /// Target processor ID from which to print
    int targetpid_;

    /// bool whether output is written to screen
    bool writetoscreen_;

    /// bool whether output is written to file
    bool writetofile_;

    /// file descriptor
    std::ofstream* outfile_;

    /// bool whether group ID is prefixed in each line
    bool prefixgroupID_;

    /// bool whether group ID is prefixed in each line
    bool groupID_;

    /// buffer the output in a stringstream
    std::stringstream buffer_;

    /// format containing temporary buffer
    std::stringstream fmt_tmp_;

  };  // class Pstream

  /// Imitate the std::endl behavior w/out the flush
  Pstream& endl(Pstream& out);

  /// Imitate the std::flush behavior
  Pstream& flush(Pstream& out);

  /// this is the IO::cout that everyone can refer to
  extern Pstream cout;

}  // namspace IO

/// Handle streaming to Pstream objects
template<typename charT>
IO::Pstream& operator<<(IO::Pstream& out, charT s)
{
  return out.stream(s);
}

/// Handle special manipulators (currently only IO::endl) that are streamed to Pstream
IO::Pstream& operator<<(IO::Pstream& out, IO::Pstream& (*pf)(IO::Pstream&));

/*----------------------------------------------------------------------*/
#endif // IO_PSTREAM_H
