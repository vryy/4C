/*----------------------------------------------------------------------*/
/*!
\brief Interface class for the Lagrange solving strategy of the augmented
       framework.

\level 3

\maintainer Matthias Mayr
*/
/*----------------------------------------------------------------------*/

#include "contact_aug_lagrange_interface.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::LAGRANGE::Interface::Interface(
    const Teuchos::RCP<CONTACT::AUG::IDataContainer>& idata_ptr)
    : ::CONTACT::AUG::Interface(idata_ptr)
{
  /* do nothing */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::LAGRANGE::Interface::Interface(
    const Teuchos::RCP<MORTAR::InterfaceDataContainer>& interfaceData_ptr, const int id,
    const Epetra_Comm& comm, const int dim, const Teuchos::ParameterList& icontact,
    const bool selfcontact, INPAR::MORTAR::RedundantStorage redundant)
    : ::CONTACT::AUG::Interface(interfaceData_ptr, id, comm, dim, icontact, selfcontact, redundant)
{
  /* left blank, nothing to do here */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::LAGRANGE::Interface::AssembleDGGLinMatrix(
    LINALG::SparseMatrix& dGGSlLinMatrix, const Epetra_Vector& cnVec) const
{
  /* do nothing */
}
