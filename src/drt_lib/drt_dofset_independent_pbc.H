/*!----------------------------------------------------------------------
\file drt_dofset_independent_pbc.H

\brief A modified set of degrees of freedom. This class inherits from
       both the independent dofset as well as the pbc dofset.

<pre>
Maintainer: Ursula Rasthofer
            rasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifndef DOFSET_INDEPENDENT_PBC_H
#define DOFSET_INDEPENDENT_PBC_H

#include <vector>
#include <map>

#include <Epetra_IntVector.h>
#include <Teuchos_RefCountPtr.hpp>
#include <Epetra_Comm.h>
#include <Epetra_Map.h>

#include "../drt_lib/drt_dofset_independent.H"
#include "../drt_lib/drt_dofset_pbc.H"


// forward declarations
namespace DRT
{
class Discretization;

/*!
\brief A set of degrees of freedom

\author
*/
class IndependentPBCDofSet: public IndependentDofSet, public PBCDofSet
{
public:

  /*!
  \brief Standard Constructor


  <pre>

  create a dofset that is independent of the other dofsets

  </pre>

  \return void

  */
  IndependentPBCDofSet(RefCountPtr<std::map<int,vector<int> > >  couplednodes);

  /*!
  \brief Destructor

  */
  virtual ~IndependentPBCDofSet();

  /// create a copy of this object
  virtual Teuchos::RCP<DofSet> Clone()
  {
    return Teuchos::rcp(new IndependentPBCDofSet(*this));
  }

  /// Assign dof numbers using all elements and nodes of the discretization.
  virtual int AssignDegreesOfFreedom(const DRT::Discretization& dis, const unsigned dspos, const int start);

protected:

private:

}; // class IndependentPBCDofSet

} // namespace DRT

#endif  // #ifndef DOFSET_INDEPENDENT_H
