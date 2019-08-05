/*---------------------------------------------------------------------*/
/*!
\brief Steepest ascent interface based on the augmented contact
       formulation.

\level 3

\maintainer Matthias Mayr
*/
/*---------------------------------------------------------------------*/


#include "contact_aug_steepest_ascent_interface.H"
#include "../drt_contact/contact_node.H"

#include "../drt_lib/drt_discret.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::STEEPESTASCENT::Interface::Interface(
    const Teuchos::RCP<CONTACT::AUG::InterfaceDataContainer>& interfaceData_ptr)
    : ::CONTACT::AUG::Interface(interfaceData_ptr)
{
  /* do nothing */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::STEEPESTASCENT::Interface::Interface(
    const Teuchos::RCP<MORTAR::InterfaceDataContainer>& interfaceData_ptr, int id,
    const Epetra_Comm& comm, int dim, const Teuchos::ParameterList& icontact, bool selfcontact,
    INPAR::MORTAR::RedundantStorage redundant)
    : ::CONTACT::AUG::Interface(interfaceData_ptr, id, comm, dim, icontact, selfcontact, redundant)
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
