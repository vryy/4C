/*---------------------------------------------------------------------*/
/*! \file
\brief Deprecated wrapper for the mesh-tying and contact managers.

\level 2


*/
/*---------------------------------------------------------------------*/

#include "4C_contact_meshtying_contact_bridge.hpp"

#include "4C_contact_manager.hpp"
#include "4C_contact_meshtying_manager.hpp"
#include "4C_fem_condition.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_linalg_mapextractor.hpp"
#include "4C_mortar_manager_base.hpp"
#include "4C_mortar_strategy_base.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  ctor (public)                                            farah 06/14|
 *----------------------------------------------------------------------*/
CONTACT::MeshtyingContactBridge::MeshtyingContactBridge(Core::FE::Discretization& dis,
    std::vector<Core::Conditions::Condition*>& meshtyingConditions,
    std::vector<Core::Conditions::Condition*>& contactConditions, double timeIntegrationMidPoint)
    : cman_(Teuchos::null), mtman_(Teuchos::null)
{
  bool onlymeshtying = false;
  bool onlycontact = false;
  bool meshtyingandcontact = false;

  // check for case
  if (meshtyingConditions.size() != 0 and contactConditions.size() != 0) meshtyingandcontact = true;

  if (meshtyingConditions.size() != 0 and contactConditions.size() == 0) onlymeshtying = true;

  if (meshtyingConditions.size() == 0 and contactConditions.size() != 0) onlycontact = true;

  // create meshtying and contact manager
  if (onlymeshtying)
  {
    mtman_ = Teuchos::rcp(new CONTACT::MtManager(dis, timeIntegrationMidPoint));
  }
  else if (onlycontact)
  {
    cman_ = Teuchos::rcp(new CONTACT::Manager(dis, timeIntegrationMidPoint));
  }
  else if (meshtyingandcontact)
  {
    mtman_ = Teuchos::rcp(new CONTACT::MtManager(dis, timeIntegrationMidPoint));
    cman_ = Teuchos::rcp(new CONTACT::Manager(dis, timeIntegrationMidPoint));
  }

  // Sanity check for writing output for each interface
  {
    const bool writeInterfaceOutput =
        Core::UTILS::IntegralValue<bool>(get_strategy().params(), "OUTPUT_INTERFACES");

    if (writeInterfaceOutput && have_contact() && contact_manager()->get_strategy().is_friction())
      FOUR_C_THROW(
          "Output for each interface does not work yet, if friction is enabled. Switch off the "
          "interface-based output in the input file (or implement/fix it for frictional contact "
          "problems.");
  }

  return;
}

/*----------------------------------------------------------------------*
 |  store_dirichlet_status                                     farah 06/14|
 *----------------------------------------------------------------------*/
void CONTACT::MeshtyingContactBridge::store_dirichlet_status(
    Teuchos::RCP<Core::LinAlg::MapExtractor> dbcmaps)
{
  if (have_meshtying()) mt_manager()->get_strategy().store_dirichlet_status(dbcmaps);
  if (have_contact()) contact_manager()->get_strategy().store_dirichlet_status(dbcmaps);

  return;
}

/*----------------------------------------------------------------------*
 |  Set displacement state                                   farah 06/14|
 *----------------------------------------------------------------------*/
void CONTACT::MeshtyingContactBridge::set_state(Teuchos::RCP<Epetra_Vector> zeros)
{
  if (have_meshtying())
    mt_manager()->get_strategy().set_state(Mortar::state_new_displacement, *zeros);
  if (have_contact())
    contact_manager()->get_strategy().set_state(Mortar::state_new_displacement, *zeros);

  return;
}

/*----------------------------------------------------------------------*
 |  Get Strategy                                             farah 06/14|
 *----------------------------------------------------------------------*/
Mortar::StrategyBase& CONTACT::MeshtyingContactBridge::get_strategy() const
{
  // if contact is involved use contact strategy!
  // contact conditions/strategies are dominating the algorithm!
  if (have_meshtying() and !have_contact())
    return mt_manager()->get_strategy();
  else
    return contact_manager()->get_strategy();
}

/*----------------------------------------------------------------------*
 |  PostprocessTractions                                     farah 06/14|
 *----------------------------------------------------------------------*/
void CONTACT::MeshtyingContactBridge::postprocess_quantities(
    Teuchos::RCP<Core::IO::DiscretizationWriter>& output)
{
  // contact
  if (have_contact()) contact_manager()->postprocess_quantities(*output);

  // meshtying
  if (have_meshtying()) mt_manager()->postprocess_quantities(*output);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::MeshtyingContactBridge::postprocess_quantities_per_interface(
    Teuchos::RCP<Teuchos::ParameterList> outputParams)
{
  // This is an optional feature, so we check if it has been enabled in the input file
  const bool writeInterfaceOutput =
      Core::UTILS::IntegralValue<bool>(get_strategy().params(), "OUTPUT_INTERFACES");
  if (writeInterfaceOutput)
  {
    // contact
    if (have_contact()) contact_manager()->postprocess_quantities_per_interface(outputParams);

    // meshtying
    if (have_meshtying()) mt_manager()->postprocess_quantities_per_interface(outputParams);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Recover lagr. mult and slave displ                       farah 06/14|
 *----------------------------------------------------------------------*/
void CONTACT::MeshtyingContactBridge::recover(Teuchos::RCP<Epetra_Vector> disi)
{
  // meshtying
  if (have_meshtying()) mt_manager()->get_strategy().recover(disi);

  // contact
  if (have_contact()) contact_manager()->get_strategy().recover(disi);

  return;
}

/*----------------------------------------------------------------------*
 |  Recover lagr. mult and slave displ                       farah 06/14|
 *----------------------------------------------------------------------*/
void CONTACT::MeshtyingContactBridge::read_restart(Core::IO::DiscretizationReader& reader,
    Teuchos::RCP<Epetra_Vector> dis, Teuchos::RCP<Epetra_Vector> zero)
{
  // contact
  if (have_contact()) contact_manager()->read_restart(reader, dis, zero);

  // meshtying
  if (have_meshtying()) mt_manager()->read_restart(reader, dis, zero);

  return;
}

/*----------------------------------------------------------------------*
 |  Write restart                                            farah 06/14|
 *----------------------------------------------------------------------*/
void CONTACT::MeshtyingContactBridge::write_restart(
    Teuchos::RCP<Core::IO::DiscretizationWriter>& output, bool forcedrestart)
{
  // contact
  if (have_contact()) contact_manager()->write_restart(*output, forcedrestart);

  // meshtying
  if (have_meshtying()) mt_manager()->write_restart(*output, forcedrestart);

  return;
}

/*----------------------------------------------------------------------*
 |  Write restart                                            farah 06/14|
 *----------------------------------------------------------------------*/
void CONTACT::MeshtyingContactBridge::update(Teuchos::RCP<Epetra_Vector> dis)
{
  // contact
  if (have_contact()) contact_manager()->get_strategy().update(dis);

  // meshtying
  if (have_meshtying()) mt_manager()->get_strategy().update(dis);

  return;
}

/*----------------------------------------------------------------------*
 |  Write restart                                             popp 11/14|
 *----------------------------------------------------------------------*/
void CONTACT::MeshtyingContactBridge::visualize_gmsh(const int istep, const int iter)
{
  // contact
  if (have_contact()) contact_manager()->get_strategy().visualize_gmsh(istep, iter);

  // meshtying
  if (have_meshtying()) mt_manager()->get_strategy().visualize_gmsh(istep, iter);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Comm& CONTACT::MeshtyingContactBridge::get_comm() const
{
  if (cman_ != Teuchos::null) return cman_->get_comm();

  if (mtman_ != Teuchos::null) return mtman_->get_comm();

  FOUR_C_THROW("can't get get_comm()");
  return cman_->get_comm();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Mortar::ManagerBase> CONTACT::MeshtyingContactBridge::contact_manager() const
{
  return cman_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Mortar::ManagerBase> CONTACT::MeshtyingContactBridge::mt_manager() const
{
  return mtman_;
}

FOUR_C_NAMESPACE_CLOSE
