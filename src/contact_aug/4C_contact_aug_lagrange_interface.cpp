/*----------------------------------------------------------------------*/
/*! \file
\brief Interface class for the Lagrange solving strategy of the augmented
       framework.

\level 3

*/
/*----------------------------------------------------------------------*/

#include "4C_contact_aug_lagrange_interface.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::Aug::Lagrange::Interface::Interface(
    const Teuchos::RCP<CONTACT::Aug::InterfaceDataContainer>& idata_ptr)
    : CONTACT::Aug::Interface(idata_ptr)
{
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::Aug::Lagrange::Interface::Interface(
    const Teuchos::RCP<Mortar::InterfaceDataContainer>& interfaceData_ptr, const int id,
    const Epetra_Comm& comm, const int dim, const Teuchos::ParameterList& icontact,
    const bool selfcontact)
    : CONTACT::Aug::Interface(interfaceData_ptr, id, comm, dim, icontact, selfcontact)
{
}
FOUR_C_NAMESPACE_CLOSE
