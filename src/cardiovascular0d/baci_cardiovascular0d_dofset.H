/*----------------------------------------------------------------------*/
/*! \file

\brief A set of degrees of freedom

\level 2


*----------------------------------------------------------------------*/
#ifndef BACI_CARDIOVASCULAR0D_DOFSET_H
#define BACI_CARDIOVASCULAR0D_DOFSET_H

#include "baci_config.H"

#include "baci_lib_discret.H"
#include "baci_lib_dofset.H"
#include "baci_mor_pod.H"

#include <Epetra_IntVector.h>
#include <Epetra_Map.h>
#include <Teuchos_RCP.hpp>

#include <list>
#include <vector>

BACI_NAMESPACE_OPEN

namespace UTILS
{
  /*!
  \brief A set of degrees of freedom

  \note This is an internal class of the Cardiovascular0D manager that one
  should not need to touch on an ordinary day. It is here to support the
  Cardiovascular0D manager class. And does all the degree of freedom assignmets
  for the Cardiovascular0Ds.

  <h3>Purpose</h3>

  This class represents one set of degrees of freedom for the
  Cardiovascular0Ds in the usual parallel fashion. That is there is a
  DofRowMap() and a DofColMap() that return the maps of the global FE
  system of equation in row and column setting respectively. These maps
  are used by the algorithm's Epetra_Vector classes amoung others.

  It is not connected to elements or nodes.
  <h3>Invariants</h3>

  There are two possible states in this class: Reset and setup. To
  change back and forth use AssignDegreesOfFreedom() and Reset().


  \author tk     */
  class Cardiovascular0DDofSet : public DRT::DofSet
  {
   public:
    /*!
    \brief Standard Constructor

    */
    Cardiovascular0DDofSet();



    //! @name Access methods

    virtual int FirstGID()
    {
      int lmin = dofrowmap_->MinMyGID();
      if (dofrowmap_->NumMyElements() == 0) lmin = INT_MAX;
      int gmin = INT_MAX;
      dofrowmap_->Comm().MinAll(&lmin, &gmin, 1);
      return gmin;
    };

    //@}

    //! @name Construction

    /// Assign dof numbers using all elements and nodes of the discretization.
    virtual int AssignDegreesOfFreedom(const Teuchos::RCP<DRT::Discretization> dis, const int ndofs,
        const int start, const Teuchos::RCP<MOR::ProperOrthogonalDecomposition> mor);

    /// reset all internal variables
    void Reset() override;

    //@}

   protected:
  };  // class Cardiovascular0DDofSet
}  // namespace UTILS


// << operator
std::ostream& operator<<(std::ostream& os, const UTILS::Cardiovascular0DDofSet& dofset);


BACI_NAMESPACE_CLOSE

#endif  // CARDIOVASCULAR0D_DOFSET_H
