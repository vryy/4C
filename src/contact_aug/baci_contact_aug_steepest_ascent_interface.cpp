/*---------------------------------------------------------------------*/
/*! \file
\brief Steepest ascent interface based on the augmented contact
       formulation.

\level 3

*/
/*---------------------------------------------------------------------*/


#include "baci_contact_aug_steepest_ascent_interface.hpp"

#include "baci_contact_node.hpp"
#include "baci_lib_discret.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::STEEPESTASCENT::Interface::Interface(
    const Teuchos::RCP<CONTACT::AUG::InterfaceDataContainer>& interfaceData_ptr)
    : CONTACT::AUG::Interface(interfaceData_ptr)
{
  /* do nothing */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::STEEPESTASCENT::Interface::Interface(
    const Teuchos::RCP<MORTAR::InterfaceDataContainer>& interfaceData_ptr, int id,
    const Epetra_Comm& comm, int dim, const Teuchos::ParameterList& icontact, bool selfcontact)
    : CONTACT::AUG::Interface(interfaceData_ptr, id, comm, dim, icontact, selfcontact)
{
  /* left blank, nothing to do here */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::AUG::INTERFACE::AssembleStrategy>
CONTACT::AUG::STEEPESTASCENT::Interface::CreateNodeBasedAssembleStrategy()
{
  const enum INPAR::CONTACT::VariationalApproach var_type = GetVariationalApproachType();

  switch (var_type)
  {
    case INPAR::CONTACT::var_complete:
    {
      typedef CONTACT::AUG::INTERFACE::CompleteAssemblePolicy complete_policy;

      return Teuchos::rcp(
          new STEEPESTASCENT::INTERFACE::NodeBasedAssembleStrategy<complete_policy>(this));
    }
    case INPAR::CONTACT::var_incomplete:
    {
      typedef CONTACT::AUG::INTERFACE::IncompleteAssemblePolicy incomplete_policy;

      return Teuchos::rcp(
          new STEEPESTASCENT::INTERFACE::NodeBasedAssembleStrategy<incomplete_policy>(this));
    }
    default:
    {
      dserror("Unknown variational approach! (var_type= \"%s\" | %d)",
          INPAR::CONTACT::VariationalApproach2String(var_type).c_str(), var_type);
      exit(EXIT_FAILURE);
    }
  }
}

FOUR_C_NAMESPACE_CLOSE
