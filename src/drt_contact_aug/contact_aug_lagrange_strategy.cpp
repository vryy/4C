/*----------------------------------------------------------------------*/
/*!
\file contact_aug_lagrange_strategy.cpp

\brief Lagrange contact solving strategy with standard Lagrangian
       multipliers based on the augmented Lagrange formulation.

\level 3

\maintainer Michael Hiermeier

\date Mar 28, 2017

*/
/*----------------------------------------------------------------------*/

#include "contact_aug_lagrange_strategy.H"
#include "contact_aug_lagrange_interface.H"

#include "../linalg/linalg_utils.H"
#include "../drt_lib/epetra_utils.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::LAGRANGE::Strategy::Strategy(
    const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr, const Epetra_Map* DofRowMap,
    const Epetra_Map* NodeRowMap, const Teuchos::ParameterList& params,
    const plain_interface_set& interfaces, int dim, const Teuchos::RCP<const Epetra_Comm>& comm,
    int maxdof)
    : CONTACT::AUG::Strategy(data_ptr, DofRowMap, NodeRowMap, params, interfaces, dim, comm, maxdof)
{
  // cast to steepest ascent interfaces
  for (plain_interface_set::const_iterator cit = interfaces.begin(); cit != interfaces.end(); ++cit)
  {
    const Teuchos::RCP<CONTACT::CoInterface>& interface = *cit;
    // test interfaces for the correct type
    Teuchos::rcp_dynamic_cast<CONTACT::AUG::LAGRANGE::Interface>(interface, true);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::LAGRANGE::Strategy::EvalStrContactRHS()
{
  if (!IsInContact() and !WasInContact() and !WasInContactLastTimeStep())
  {
    Data().StrContactRhsPtr() = Teuchos::null;
    return;
  }
  Data().StrContactRhsPtr() = Teuchos::rcp(new Epetra_Vector(*ProblemDofs(), true));


  // For self contact, slave and master sets may have changed,
  if (IsSelfContact())
    dserror(
        "ERROR: Augmented Lagrange Formulation: Self contact is not yet "
        "considered!");

  // --- add contact force terms ----------------------------------------------
  // *** Slave side ***
  Epetra_Vector augfs_exp(*ProblemDofs());
  LINALG::Export(Data().SlForceLm(), augfs_exp);
  Data().StrContactRhs().Scale(-1.0, augfs_exp);

  // Master side
  Epetra_Vector augfm_exp(*ProblemDofs());
  LINALG::Export(Data().MaForceLm(), augfm_exp);
  CATCH_EPETRA_ERROR(Data().StrContactRhs().Update(-1.0, augfm_exp, 1.0));

  return;
}
