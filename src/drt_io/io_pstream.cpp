/*----------------------------------------------------------------------*/
/*!

\brief A substitute for STL cout for parallel and complex output schemes.

\level 0

\maintainer Martin Kronbichler
*/

/*----------------------------------------------------------------------*/
/* headers */
#include "io_pstream.H"
#include <Teuchos_oblackholestream.hpp>

namespace IO
{
  /// this is the IO::cout that everyone can refer to
  Pstream cout;
}  // namespace IO


/*----------------------------------------------------------------------*
 * empty constructor                                          wic 11/12 *
 *----------------------------------------------------------------------*/
IO::Pstream::Pstream()
    : is_initialized_(false),
      comm_(Teuchos::null),
      targetpid_(-2),
      writetoscreen_(false),
      writetofile_(false),
      outfile_(NULL),
      prefixgroupID_(false),
      groupID_(-2),
      buffer_(std::string()),
      requestedoutputlevel_(undef),
      level_(new Level(this))
{
}


/*----------------------------------------------------------------------*
 * destructor                                                 wic 09/16 *
 *----------------------------------------------------------------------*/
IO::Pstream::~Pstream()
{
  if (level_) delete level_;
  level_ = NULL;

  if (blackholestream_) delete blackholestream_;
  blackholestream_ = NULL;

  mystream_ = NULL;

  this->close();
}

/*----------------------------------------------------------------------*
 * configure the output                                       wic 11/12 *
 *----------------------------------------------------------------------*/
void IO::Pstream::setup(const bool writetoscreen, const bool writetofile, const bool prefixgroupID,
    const IO::verbositylevel level, Teuchos::RCP<Epetra_Comm> comm, const int targetpid,
    const int groupID, const std::string fileprefix)
{
  // make sure that setup is called only once or we get unpredictable behavior
  if (is_initialized_) dserror("Thou shalt not call setup on the output twice!");
  is_initialized_ = true;

  requestedoutputlevel_ = level;
  comm_ = comm;
  targetpid_ = targetpid;
  writetoscreen_ = writetoscreen;
  writetofile_ = writetofile;
  outfile_ = NULL;
  prefixgroupID_ = prefixgroupID;
  groupID_ = groupID;

  // make sure the target processor exists
  if (targetpid_ >= comm_->NumProc()) dserror("Chosen target processor does not exist.");

  // prepare the file handle
  if (OnPid() and writetofile_)
  {
    std::stringstream fname;
    fname << fileprefix << ".p" << std::setfill('0') << std::setw(2) << comm_->MyPID() << ".log";
    outfile_ = new std::ofstream(fname.str().c_str());
    if (!outfile_) dserror("could not open output file");
  }

  // prepare the very first line of output
  if (OnPid() and prefixgroupID_) buffer_ << groupID_ << ": ";

  // setup mystream
  blackholestream_ = new Teuchos::oblackholestream;
  if (writetoscreen_ and (comm_->MyPID() == targetpid_ or targetpid_ < 0))
    mystream_ = &std::cout;
  else
    mystream_ = blackholestream_;
}

/*-----------------------------------------------------------------------*
 * return a std::ostream following the restrictions      hiermeier 12/17 *
 *-----------------------------------------------------------------------*/
std::ostream& IO::Pstream::os(const verbositylevel level) const
{
  if (not is_initialized_) dserror("Setup the output before you use it!");

  if (level <= RequestedOutputLevel())
  {
    if (prefixgroupID_) *mystream_ << prefixgroupID_ << ": ";
    return *mystream_;
  }
  else
    return *blackholestream_;
}

/*----------------------------------------------------------------------*
 * close open file handles and reset pstream                  wic 11/12 *
 *----------------------------------------------------------------------*/
void IO::Pstream::close()
{
  if (not is_initialized_) return;

  is_initialized_ = false;
  comm_ = Teuchos::null;
  targetpid_ = -2;
  writetoscreen_ = false;
  writetofile_ = false;

  // close file
  if (outfile_)
  {
    outfile_->close();
    delete outfile_;
  }
  outfile_ = NULL;

  prefixgroupID_ = false;
  groupID_ = -2;

  // flush the buffer
  if (writetoscreen_ and OnPid() and buffer_.str().size() > 0)
    std::cout << buffer_.str() << std::flush;
  buffer_.str(std::string());
}


/*----------------------------------------------------------------------*
 * writes the buffer to screen                                wic 11/12 *
 *----------------------------------------------------------------------*/
void IO::Pstream::flush()
{
  if (not is_initialized_) dserror("Setup the output before you use it!");

  if (OnPid() and writetoscreen_ and buffer_.str().size() > 0)
  {
    std::cout << buffer_.str();
    std::flush(std::cout);
    buffer_.str(std::string());
  }

  return;
}


/*----------------------------------------------------------------------*
 * return whether this is a target processor                  wic 11/12 *
 *----------------------------------------------------------------------*/
bool IO::Pstream::OnPid()
{
  if (targetpid_ < 0) return true;
  return (comm_->MyPID() == targetpid_);
}


/*----------------------------------------------------------------------*
 * set output level                                           wic 09/16 *
 *----------------------------------------------------------------------*/
IO::Level& IO::Pstream::operator()(const verbositylevel level) { return level_->SetLevel(level); }


/*----------------------------------------------------------------------*
 * Imitate the std::endl behavior w/out the flush             wic 11/12 *
 *----------------------------------------------------------------------*/
IO::Pstream& IO::endl(IO::Pstream& out)
{
  out << "\n";
  return out;
}

/*----------------------------------------------------------------------*
 * Imitate the std::endl behavior w/out the flush             wic 09/16 *
 *----------------------------------------------------------------------*/
IO::Level& IO::endl(IO::Level& out)
{
  out << "\n";
  return out;
}

/*----------------------------------------------------------------------*
 * Imitate the std::flush behavior                            wic 11/12 *
 *----------------------------------------------------------------------*/
IO::Pstream& IO::flush(IO::Pstream& out)
{
  out.flush();
  return out;
}


/*----------------------------------------------------------------------*
 * Imitate the std::flush behavior                            wic 11/12 *
 *----------------------------------------------------------------------*/
IO::Level& IO::flush(IO::Level& out)
{
  out.flush();
  return out;
}


/*----------------------------------------------------------------------*
 * writes the buffer to screen                                wic 09/16 *
 *----------------------------------------------------------------------*/
void IO::Level::flush()
{
  if (level_ <= pstream_->RequestedOutputLevel()) pstream_->flush();

  return;
}


/*----------------------------------------------------------------------*
 * Handle special manipulators                                wic 11/12 *
 *----------------------------------------------------------------------*/
IO::Pstream& operator<<(IO::Pstream& out, IO::Pstream& (*pf)(IO::Pstream&)) { return pf(out); }

/*----------------------------------------------------------------------*
 * Handle special manipulators                                wic 09/16 *
 *----------------------------------------------------------------------*/
IO::Level& operator<<(IO::Level& out, IO::Level& (*pf)(IO::Level&)) { return pf(out); }
