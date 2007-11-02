/*!----------------------------------------------------------------------
\file drt_cdofset.H

\brief A set of degrees of freedom special for contact

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef DRT_CDOFSET_H
#define DRT_CDOFSET_H


#include "../drt_lib/drt_dofset.H"
#include "../drt_lib/drt_discret.H"
using namespace std;

/*!
\brief CONTACT: namespace of the baci contact algorithms

*/
namespace CONTACT
{
/*!
\class CONTACT::CDofSet

\brief A set of degrees of freedom special for contact

*/
class CDofSet : public DRT::DofSet
{
public:



  /*!
  \brief Standard Constructor

  */
  explicit CDofSet();

  /*!
  \brief Destructor

  */
  virtual ~CDofSet() { return; }


  //@}

  //! @name Construction

  /// Assign dof numbers to all elements and nodes of the discretization.
  virtual int AssignDegreesOfFreedom(const DRT::Discretization& dis, const int start);

  //@}

protected:


}; // class CDofSet
} // namespace CONTACT




#endif  // #ifndef DRT_CDOFSET_H
#endif  // #ifdef CCADISCRET
