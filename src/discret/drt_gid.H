/*!----------------------------------------------------------------------
\file drt_gid.H

\class DRT::Gidout

\brief output to gid

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef GID_H
#define GID_H


using namespace std;
#include "drt_dserror.H"

/*!
\brief DRT: namespace of the ccarat discretization module

*/
namespace DRT
{

/*!
\class Gid

\brief A virtual class all nodes that are used in DRT have to implement

\author gee (gee@lnm.mw.tum.de)
*/
class Gidout
{
public:

  // @{ \name Constructors and destructors

  /*!
  \brief Standard Constructor
  */
  Gidout();

  /*!
  \brief Copy Constructor
  */
  Gidout(const DRT::Gidout& old);

  /*!
  \brief Destructor

  */
  virtual ~Gidout();

  //@}

  // @{ \name Acess methods

  //@}

private:

}; // class Gidout
} // namespace DRT











#endif  // #ifndef GID_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
