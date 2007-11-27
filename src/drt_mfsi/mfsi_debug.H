/*----------------------------------------------------------------------*/
/*!
\file mfsi_debug.H

\brief Very simple output facility

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef MFSI_DEBUG_H
#define MFSI_DEBUG_H

#include "../drt_lib/drt_discret.H"

namespace MFSI
{
  /// Very simple output facility for debugging only
  /*!

  I'm not sure whether we want to keep this around. But then, some
  debug output is needed during developement. And I keep writing code
  like that. So for the time being...

  This is single processor, 2d, quad4 only. Extend it if you have to.

  \author u.kue
  \date 11/07
  */
  class Debug
  {
  public:

    /// write a vector on this discretization to a gnuplot file
    void DumpVector(std::string name,
                    const DRT::Discretization& discret,
                    const Epetra_Vector& v) const;

  private:
    mutable std::map<std::string,int> file_;
  };
}

#endif
#endif
