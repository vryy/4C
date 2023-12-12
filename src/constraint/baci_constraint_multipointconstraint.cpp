/*----------------------------------------------------------------------*/
/*! \file
\brief Basic constraint class, dealing with multi point constraints
\level 2


*----------------------------------------------------------------------*/


#include "baci_constraint_multipointconstraint.H"

#include "baci_lib_discret.H"
#include "baci_linalg_utils_sparse_algebra_create.H"
#include "baci_linalg_utils_sparse_algebra_manipulation.H"

#include <iostream>

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  ctor (public)                                               tk 07/08|
 *----------------------------------------------------------------------*/
UTILS::MPConstraint::MPConstraint(Teuchos::RCP<DRT::Discretization> discr,
    const std::string& conditionname, int& minID, int& maxID)
    : UTILS::Constraint(discr, conditionname, minID, maxID)
{
  return;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                               tk 07/08|
 *----------------------------------------------------------------------*/
UTILS::MPConstraint::MPConstraint(
    Teuchos::RCP<DRT::Discretization> discr, const std::string& conditionname)
    : UTILS::Constraint(discr, conditionname)
{
  return;
}

/// Set state of the underlying constraint discretization
void UTILS::MPConstraint::SetConstrState(const std::string& state,  ///< name of state to set
    Teuchos::RCP<const Epetra_Vector> V                             ///< values to set
)
{
  if (constrtype_ != none)
  {
    std::map<int, Teuchos::RCP<DRT::Discretization>>::iterator discrit;
    for (discrit = constraintdis_.begin(); discrit != constraintdis_.end(); ++discrit)
    {
      Teuchos::RCP<Epetra_Vector> tmp =
          CORE::LINALG::CreateVector(*(discrit->second)->DofColMap(), false);
      CORE::LINALG::Export(*V, *tmp);
      (discrit->second)->ClearState();
      (discrit->second)->SetState(state, tmp);
    }
  }
}

BACI_NAMESPACE_CLOSE
