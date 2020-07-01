/*----------------------------------------------------------------------*/
/*! \file

\brief A substitute for STL cout for parallel and complex output schemes.

\level 0

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
  // forward declaration
  class Level;

  /// enum for verbosity level
  enum verbositylevel
  {
    undef = -1,  ///< do not use!
    minimal,     ///< one line per time step should be enough
    standard,    ///< a few lines per time step is okay
    verbose,     ///< detailed information, e.g. nonlinear convergence output, execution steps
    debug        ///< to be used for permanent debugging information (not for temporary hacks)
  };


  /*====================================================================*/
  /*!
   * \brief This object allows to write to std::cout in a very fancy way.
   *        The output can be regulated by various parameters, in particular
   *        there are:
   *        o writing to screen
   *        o writing to files
   *        o modifying the output for use with nested parallelism
   *        o limit output to certain procs
   *        o selecting an output level
   *
   * \author wichmann & hammerl \date 11/12
   */
  class Pstream
  {
   public:
    /// This empty constructor is called when the global object is instantiated.
    Pstream();

    /// Destructor
    virtual ~Pstream();

    /// configure the output. Must be called before IO::cout can be used. Is
    /// currently called in the global problem using the params specified in the
    /// IO section of the .dat-file.
    void setup(const bool writetoscreen,  ///< bool whether output is written to screen
        const bool writetofile,           ///< bool whether output is written to file
        const bool prefixgroupID,         ///< bool whether group ID is prefixed in each line
        const IO::verbositylevel level,   ///< verbosity level
        Teuchos::RCP<Epetra_Comm> comm,   ///< MPI communicator
        const int targetpid,              ///< target processor ID from which to print
        const int groupID,                ///< the ID
        const std::string fileprefix      ///< prefix for the output file
    );

    /// must be called to close open file handles and resets the IO::cout object to a pristine
    /// state. This is currently called by the destructor of the global problem singleton.
    void close();

    /// writes the buffer to screen
    void flush();

    /// This handles the actual printing to screen/writing to a file.
    template <typename CharT>
    Pstream& stream(CharT s)  ///< text to be added
    {
      if (not is_initialized_) dserror("Setup the output before you use it!");

      // are we on a proc that is allowed to output
      if (OnPid())
      {
        // store formatting in this semi persistent buffer
        fmt_tmp_ << s;
        if (fmt_tmp_.str().size() <= 0) return *this;

        // actual content goes into the string str
        std::string str = fmt_tmp_.str();
        fmt_tmp_.str(std::string());

        // write to file
        if (writetofile_)
        {
#ifdef DEBUG
          if (!outfile_) dserror("outputfile does not exist - file handle lost");
#endif
          (*outfile_) << str;
        }

        // write to screen
        if (writetoscreen_)
        {
          size_t oldpos = 0;
          size_t pos = 0;
          while (true)
          {
            pos = str.find('\n', oldpos);  // do include the carriage return itself

            if (pos == std::string::npos)
            {
              buffer_ << str.substr(oldpos, str.size() - oldpos);
              break;
            }

            buffer_ << str.substr(oldpos, pos - oldpos + 1);

            oldpos = pos + 1;

            // we are using a hard coded std::cout here, as writing to screen means std::cout by
            // definition
            std::cout << buffer_.str();
            std::flush(std::cout);
            buffer_.str(std::string());

            if (prefixgroupID_) buffer_ << groupID_ << ": ";
          }
        }
      }

      return *this;
    }

    std::stringstream& cout_replacement() { return fmt_tmp_; }

    /// Set verbosity level
    Level& operator()(const verbositylevel level);

    /// \brief Return a pure screen ostream (a.k.a. std::cout)
    /** \note This is mainly a workaround to use the verbosity levels as well
     *  as the processor restriction also for standard print() methods which
     *  typically expect a std::ostream object as input argument.
     *
     *  \author hiermeier 12/17 */
    std::ostream& os(const verbositylevel level = undef) const;

    /// Return verbosity level
    inline verbositylevel RequestedOutputLevel() const { return requestedoutputlevel_; }

   private:
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

    /// group ID that is prefixed in each line if desired
    int groupID_;

    /// buffer the output in a stringstream
    std::stringstream buffer_;

    /// format containing temporary buffer
    std::stringstream fmt_tmp_;

    verbositylevel requestedoutputlevel_;

    // internal variable of Pstream
    Level* level_;

    /// ostream for all non-target procs and inadmissible cases
    std::ostream* blackholestream_ = NULL;

    /// ostream if all prerequisites are fulfilled
    std::ostream* mystream_ = NULL;

  };  // class Pstream

  /// Imitate the std::endl behavior w/out the flush
  Pstream& endl(Pstream& out);

  /// Imitate the std::flush behavior
  Pstream& flush(Pstream& out);

  /// this is the IO::cout that everyone can refer to
  extern Pstream cout;


  /*====================================================================*/
  /*!
   * \brief This object handles output if a level was specified.
   *        Everything is then streamed to Pstream.
   *
   * \author wichmann & hammerl \date 09/16
   */
  class Level
  {
   public:
    Level(Pstream* pstream) : pstream_(pstream), level_(undef) {}

    /// set verbosity level for data to be written
    Level& SetLevel(const verbositylevel level)
    {
      level_ = level;
      return *this;
    }

    /// Handle streaming to Level objects
    template <typename CharT>
    Level& stream(CharT s)  ///< text to be added
    {
      if (level_ <= pstream_->RequestedOutputLevel()) pstream_->stream(s);
      return *this;
    }

    /// writes the buffer to screen
    void flush();

   private:
    /// Shelter copy constructor
    Level(const Level& old);

    /// wrapped pstream object
    Pstream* pstream_;

    /// current verbosity level of data to be printed
    verbositylevel level_;
  };

  /// Imitate the std::endl behavior w/out the flush
  Level& endl(Level& out);

  /// Imitate the std::flush behavior
  Level& flush(Level& out);

}  // namespace IO

/// Handle streaming to Pstream objects
template <typename charT>
IO::Pstream& operator<<(IO::Pstream& out, charT s)
{
  return out.stream(s);
}

/// Handle streaming to Level objects
template <typename charT>
IO::Level& operator<<(IO::Level& out, charT s)
{
  return out.stream(s);
}

/// Handle special manipulators (currently only IO::endl) that are streamed to Pstream
IO::Pstream& operator<<(IO::Pstream& out, IO::Pstream& (*pf)(IO::Pstream&));

/// Handle special manipulators (currently only IO::endl) that are streamed to Level
IO::Level& operator<<(IO::Level& out, IO::Level& (*pf)(IO::Level&));

/*----------------------------------------------------------------------*/
#endif  // IO_PSTREAM_H
