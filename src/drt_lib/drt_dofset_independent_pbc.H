/*!----------------------------------------------------------------------
\file drt_dofset_independent_pbc.H

\brief A modified set of degrees of freedom. This class inherits from
       both the independent dofset as well as the pbc dofset.

<pre>
\level 2

\maintainer Ursula Rasthofer
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifndef DOFSET_INDEPENDENT_PBC_H
#define DOFSET_INDEPENDENT_PBC_H

#include <vector>
#include <map>

#include <Epetra_IntVector.h>
#include <Teuchos_RCP.hpp>
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
  IndependentPBCDofSet(Teuchos::RCP<std::map<int,std::vector<int> > >  couplednodes);

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
