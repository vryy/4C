/*----------------------------------------------------------------------*/
/*! \file
\brief Lagrange contact solving strategy with standard Lagrangian
       multipliers based on the augmented Lagrange formulation.

\level 3

*/
/*----------------------------------------------------------------------*/

#include "4C_contact_aug_lagrange_strategy.hpp"

#include "4C_contact_aug_lagrange_interface.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_utils_epetra_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::LAGRANGE::Strategy::Strategy(
    const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
    const Epetra_Map* dof_row_map, const Epetra_Map* NodeRowMap,
    const Teuchos::ParameterList& params, const plain_interface_set& interfaces, int dim,
    const Teuchos::RCP<const Epetra_Comm>& comm, int maxdof)
    : CONTACT::AUG::Strategy(
          data_ptr, dof_row_map, NodeRowMap, params, interfaces, dim, comm, maxdof)
{
  // cast to steepest ascent interfaces
  for (plain_interface_set::const_iterator cit = interfaces.begin(); cit != interfaces.end(); ++cit)
  {
    const Teuchos::RCP<CONTACT::Interface>& interface = *cit;
    // test interfaces for the correct type
    Teuchos::rcp_dynamic_cast<CONTACT::AUG::LAGRANGE::Interface>(interface, true);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::LAGRANGE::Strategy::eval_str_contact_rhs()
{
  if (!IsInContact() and !WasInContact() and !was_in_contact_last_time_step())
  {
    data().StrContactRhsPtr() = Teuchos::null;
    return;
  }
  data().StrContactRhsPtr() = Teuchos::rcp(new Epetra_Vector(*ProblemDofs(), true));


  // For self contact, slave and master sets may have changed,
  if (IsSelfContact())
    FOUR_C_THROW(
        "ERROR: Augmented Lagrange Formulation: Self contact is not yet "
        "considered!");

  // --- add contact force terms ----------------------------------------------
  // *** Slave side ***
  Epetra_Vector augfs_exp(*ProblemDofs());
  CORE::LINALG::Export(data().SlForceLm(), augfs_exp);
  data().StrContactRhs().Scale(-1.0, augfs_exp);

  // Master side
  Epetra_Vector augfm_exp(*ProblemDofs());
  CORE::LINALG::Export(data().MaForceLm(), augfm_exp);
  CATCH_EPETRA_ERROR(data().StrContactRhs().Update(-1.0, augfm_exp, 1.0));

  return;
}

FOUR_C_NAMESPACE_CLOSE
