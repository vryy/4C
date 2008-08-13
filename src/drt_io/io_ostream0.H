/*----------------------------------------------------------------------*/
/*!
\file io_ostream0.H

\brief An extension of STL cout for parallel, i.e. it prints only
       if called on processor 0, but does nothing on every other chip

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */
#ifndef IO_OSTREAM0_H
#define IO_OSTREAM0_H
#ifdef CCADISCRET

/*----------------------------------------------------------------------*/
/* headers */
#include <ios>
#include <iostream>
#include <sstream>
#include <Epetra_Comm.h>

#include "../drt_lib/drt_discret.H"

/*----------------------------------------------------------------------*/
/* namespace */
namespace IO
{
  
  /*====================================================================*/
  /*!
   * \brief This object allows to write to std::cout
   *        (if used at construction time) only on the target processor,
   *        which is ordinarily the master processor with rank/PID equal to zero.
   *        The remaining processors keep quiet.
   *
   * References:
   * - http://www.cplusplus.com/reference/iostream/ostream/
   * - http://www.cplusplus.com/reference/iostream/ostringstream/
   * - http://www.tacc.utexas.edu/services/userguides/pgi/pgC++_lib/stdlibug/cre_2288.htm
   *
   * \author bborn \date 08/08
   */
  class Ostream0 : public std::ostream
  {

  public :

    /// Constructor
    Ostream0(
      const Epetra_Comm& comm,  ///< MPI communicator
      std::ostream& ostr,  ///< output stream, this is going to be std::cout
      const int targetpid = 0  ///< target processor ID from which to print,
                               ///< default master
    )
    : std::ostream(ostr.rdbuf()),  // this grabs the std::streambuf of ostr
      comm_(comm),
      emptyoss_(),
      targetpid_(targetpid)
    {
      if (targetpid_ >= comm_.NumProc())
        dserror("Chosen target processor does not exist.");
    }

    /// Return wether on target processor
    bool OnPid() { return (comm_.MyPID() == targetpid_); }

    /// Return reference to empty out-stream
    std::ostringstream& EmptyOSS() { return emptyoss_; }

  private :

    /// Shelter copy constructor
    Ostream0(const Ostream0& old);

    /// MPI communicator
    const Epetra_Comm& comm_;

    /// empty std::ostringstream which <b>is a</b> std::ostream
    std::ostringstream emptyoss_;

    /// Target processor ID from which to print
    const int targetpid_;

  };  // class Ostream0

}  // namspace IO

/*======================================================================*/
/*! 
 * \brief The core function. This allows to write to std::cout
 *        (if out is a refence on this variable) only on
 *        the single processor with rank equal IO::Ostream0::targetpid_. 
 *        The remaining processors keep quiet.
 *
 * \author bborn \date 08/08
 */
template <typename CharT> inline std::ostream& operator<<(
  IO::Ostream0& out,  ///< incoming of output stream
  CharT s  ///< text to be added
)
{
  if (out.OnPid())
  {
    std::ostream& str = out;
    str << s;
    return str;
  }
  else
  {
    return out.EmptyOSS();
  }
}

/*----------------------------------------------------------------------*/
#endif // CCADISCRET
#endif // IO_OSTREAM0_H
