/*----------------------------------------------------------------------*/
/*! \file
\brief Basic constraint class, dealing with multi point constraints
\level 2


*----------------------------------------------------------------------*/


#include "baci_constraint_multipointconstraint.hpp"

#include "baci_lib_discret.hpp"
#include "baci_linalg_utils_sparse_algebra_create.hpp"
#include "baci_linalg_utils_sparse_algebra_manipulation.hpp"

#include <iostream>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  ctor (public)                                               tk 07/08|
 *----------------------------------------------------------------------*/
CONSTRAINTS::MPConstraint::MPConstraint(Teuchos::RCP<DRT::Discretization> discr,
    const std::string& conditionname, int& minID, int& maxID)
    : CONSTRAINTS::Constraint(discr, conditionname, minID, maxID)
{
  return;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                               tk 07/08|
 *----------------------------------------------------------------------*/
CONSTRAINTS::MPConstraint::MPConstraint(
    Teuchos::RCP<DRT::Discretization> discr, const std::string& conditionname)
    : CONSTRAINTS::Constraint(discr, conditionname)
{
  return;
}

/// Set state of the underlying constraint discretization
void CONSTRAINTS::MPConstraint::SetConstrState(const std::string& state,  ///< name of state to set
    Teuchos::RCP<const Epetra_Vector> V                                   ///< values to set
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

FOUR_C_NAMESPACE_CLOSE
