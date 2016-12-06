/*!----------------------------------------------------------------------
\file drt_dofset_independent.H

\brief This class is inherited from dofset and replaces the
       method AssignDegreesOfFreedom by a version that does not query
       the static_dofsets_ list for the max GID, but always starts from
       0. Also, it does not register with the static_dofsets_ list. This
       class is intended to be used for xfem approaches. It provides a
       dofset without xfem dofs for output routines.

\level 2

\maintainer Martin Kronbichler

*----------------------------------------------------------------------*/
#ifndef DOFSET_INDEPENDENT_H
#define DOFSET_INDEPENDENT_H

#include <vector>
#include <map>

#include <Epetra_IntVector.h>
#include <Teuchos_RCP.hpp>
#include <Epetra_Comm.h>
#include <Epetra_Map.h>

#include "../drt_lib/drt_dofset.H"


// forward declarations
namespace DRT
{
class Discretization;

/*!
\brief A set of degrees of freedom

\author
*/
class IndependentDofSet: virtual public DofSet
{
public:

  /*!
  \brief Standard Constructor


  <pre>

  create a dofset that is independent of the other dofsets

  </pre>

  \return void

  */
  IndependentDofSet(bool ignoreminnodegid = false);

  /*!
  \brief Copy constructor

  */
  IndependentDofSet(const IndependentDofSet& old);

  /*!
  \brief Destructor

  */
  virtual ~IndependentDofSet();

  /// create a copy of this object
  virtual Teuchos::RCP<DofSet> Clone()
  {
    return Teuchos::rcp(new IndependentDofSet(*this));
  }

  /// Assign dof numbers using all elements and nodes of the discretization.
  virtual int AssignDegreesOfFreedom(const DRT::Discretization& dis, const unsigned dspos, const int start);

protected:
  bool ignoreminnodegid_;   //< bool whether minnodegid is taken from the discretization or ignored

private:

}; // class StdDofSet

}

#endif  // #ifndef DOFSET_INDEPENDENT_H
