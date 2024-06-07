/*---------------------------------------------------------------------*/
/*! \file
\brief Steepest ascent interface based on the augmented contact
       formulation.

\level 3

*/
/*---------------------------------------------------------------------*/


#include "4C_contact_aug_steepest_ascent_interface.hpp"

#include "4C_contact_node.hpp"
#include "4C_fem_discretization.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::Aug::SteepestAscent::Interface::Interface(
    const Teuchos::RCP<CONTACT::Aug::InterfaceDataContainer>& interfaceData_ptr)
    : CONTACT::Aug::Interface(interfaceData_ptr)
{
  /* do nothing */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::Aug::SteepestAscent::Interface::Interface(
    const Teuchos::RCP<Mortar::InterfaceDataContainer>& interfaceData_ptr, int id,
    const Epetra_Comm& comm, int dim, const Teuchos::ParameterList& icontact, bool selfcontact)
    : CONTACT::Aug::Interface(interfaceData_ptr, id, comm, dim, icontact, selfcontact)
{
  /* left blank, nothing to do here */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::Aug::INTERFACE::AssembleStrategy>
CONTACT::Aug::SteepestAscent::Interface::create_node_based_assemble_strategy()
{
  const enum Inpar::CONTACT::VariationalApproach var_type = get_variational_approach_type();

  switch (var_type)
  {
    case Inpar::CONTACT::var_complete:
    {
      typedef CONTACT::Aug::INTERFACE::CompleteAssemblePolicy complete_policy;

      return Teuchos::rcp(
          new SteepestAscent::INTERFACE::NodeBasedAssembleStrategy<complete_policy>(this));
    }
    case Inpar::CONTACT::var_incomplete:
    {
      typedef CONTACT::Aug::INTERFACE::IncompleteAssemblePolicy incomplete_policy;

      return Teuchos::rcp(
          new SteepestAscent::INTERFACE::NodeBasedAssembleStrategy<incomplete_policy>(this));
    }
    default:
    {
      FOUR_C_THROW("Unknown variational approach! (var_type= \"%s\" | %d)",
          Inpar::CONTACT::VariationalApproach2String(var_type).c_str(), var_type);
      exit(EXIT_FAILURE);
    }
  }
}

FOUR_C_NAMESPACE_CLOSE
