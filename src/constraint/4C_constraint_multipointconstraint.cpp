/*----------------------------------------------------------------------*/
/*! \file
\brief Basic constraint class, dealing with multi point constraints
\level 2


*----------------------------------------------------------------------*/


#include "4C_constraint_multipointconstraint.hpp"

#include "4C_lib_discret.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"

#include <iostream>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  ctor (public)                                               tk 07/08|
 *----------------------------------------------------------------------*/
CONSTRAINTS::MPConstraint::MPConstraint(Teuchos::RCP<Discret::Discretization> discr,
    const std::string& conditionname, int& minID, int& maxID)
    : CONSTRAINTS::Constraint(discr, conditionname, minID, maxID)
{
  return;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                               tk 07/08|
 *----------------------------------------------------------------------*/
CONSTRAINTS::MPConstraint::MPConstraint(
    Teuchos::RCP<Discret::Discretization> discr, const std::string& conditionname)
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
    std::map<int, Teuchos::RCP<Discret::Discretization>>::iterator discrit;
    for (discrit = constraintdis_.begin(); discrit != constraintdis_.end(); ++discrit)
    {
      Teuchos::RCP<Epetra_Vector> tmp =
          Core::LinAlg::CreateVector(*(discrit->second)->DofColMap(), false);
      Core::LinAlg::Export(*V, *tmp);
      (discrit->second)->ClearState();
      (discrit->second)->set_state(state, tmp);
    }
  }
}

FOUR_C_NAMESPACE_CLOSE
