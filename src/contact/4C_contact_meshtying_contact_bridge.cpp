/*---------------------------------------------------------------------*/
/*! \file
\brief Deprecated wrapper for the mesh-tying and contact managers.

\level 2


*/
/*---------------------------------------------------------------------*/

#include "4C_contact_meshtying_contact_bridge.hpp"

#include "4C_contact_manager.hpp"
#include "4C_contact_meshtying_manager.hpp"
#include "4C_discretization_condition.hpp"
#include "4C_lib_discret.hpp"
#include "4C_linalg_mapextractor.hpp"
#include "4C_mortar_manager_base.hpp"
#include "4C_mortar_strategy_base.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  ctor (public)                                            farah 06/14|
 *----------------------------------------------------------------------*/
CONTACT::MeshtyingContactBridge::MeshtyingContactBridge(Discret::Discretization& dis,
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
        Core::UTILS::IntegralValue<bool>(GetStrategy().Params(), "OUTPUT_INTERFACES");

    if (writeInterfaceOutput && HaveContact() && ContactManager()->GetStrategy().Friction())
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
  if (HaveMeshtying()) MtManager()->GetStrategy().store_dirichlet_status(dbcmaps);
  if (HaveContact()) ContactManager()->GetStrategy().store_dirichlet_status(dbcmaps);

  return;
}

/*----------------------------------------------------------------------*
 |  Set displacement state                                   farah 06/14|
 *----------------------------------------------------------------------*/
void CONTACT::MeshtyingContactBridge::set_state(Teuchos::RCP<Epetra_Vector> zeros)
{
  if (HaveMeshtying()) MtManager()->GetStrategy().set_state(Mortar::state_new_displacement, *zeros);
  if (HaveContact())
    ContactManager()->GetStrategy().set_state(Mortar::state_new_displacement, *zeros);

  return;
}

/*----------------------------------------------------------------------*
 |  Get Strategy                                             farah 06/14|
 *----------------------------------------------------------------------*/
Mortar::StrategyBase& CONTACT::MeshtyingContactBridge::GetStrategy() const
{
  // if contact is involved use contact strategy!
  // contact conditions/strategies are dominating the algorithm!
  if (HaveMeshtying() and !HaveContact())
    return MtManager()->GetStrategy();
  else
    return ContactManager()->GetStrategy();
}

/*----------------------------------------------------------------------*
 |  PostprocessTractions                                     farah 06/14|
 *----------------------------------------------------------------------*/
void CONTACT::MeshtyingContactBridge::postprocess_quantities(
    Teuchos::RCP<Core::IO::DiscretizationWriter>& output)
{
  // contact
  if (HaveContact()) ContactManager()->postprocess_quantities(*output);

  // meshtying
  if (HaveMeshtying()) MtManager()->postprocess_quantities(*output);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::MeshtyingContactBridge::postprocess_quantities_per_interface(
    Teuchos::RCP<Teuchos::ParameterList> outputParams)
{
  // This is an optional feature, so we check if it has been enabled in the input file
  const bool writeInterfaceOutput =
      Core::UTILS::IntegralValue<bool>(GetStrategy().Params(), "OUTPUT_INTERFACES");
  if (writeInterfaceOutput)
  {
    // contact
    if (HaveContact()) ContactManager()->postprocess_quantities_per_interface(outputParams);

    // meshtying
    if (HaveMeshtying()) MtManager()->postprocess_quantities_per_interface(outputParams);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Recover lagr. mult and slave displ                       farah 06/14|
 *----------------------------------------------------------------------*/
void CONTACT::MeshtyingContactBridge::Recover(Teuchos::RCP<Epetra_Vector> disi)
{
  // meshtying
  if (HaveMeshtying()) MtManager()->GetStrategy().Recover(disi);

  // contact
  if (HaveContact()) ContactManager()->GetStrategy().Recover(disi);

  return;
}

/*----------------------------------------------------------------------*
 |  Recover lagr. mult and slave displ                       farah 06/14|
 *----------------------------------------------------------------------*/
void CONTACT::MeshtyingContactBridge::read_restart(Core::IO::DiscretizationReader& reader,
    Teuchos::RCP<Epetra_Vector> dis, Teuchos::RCP<Epetra_Vector> zero)
{
  // contact
  if (HaveContact()) ContactManager()->read_restart(reader, dis, zero);

  // meshtying
  if (HaveMeshtying()) MtManager()->read_restart(reader, dis, zero);

  return;
}

/*----------------------------------------------------------------------*
 |  Write restart                                            farah 06/14|
 *----------------------------------------------------------------------*/
void CONTACT::MeshtyingContactBridge::write_restart(
    Teuchos::RCP<Core::IO::DiscretizationWriter>& output, bool forcedrestart)
{
  // contact
  if (HaveContact()) ContactManager()->write_restart(*output, forcedrestart);

  // meshtying
  if (HaveMeshtying()) MtManager()->write_restart(*output, forcedrestart);

  return;
}

/*----------------------------------------------------------------------*
 |  Write restart                                            farah 06/14|
 *----------------------------------------------------------------------*/
void CONTACT::MeshtyingContactBridge::Update(Teuchos::RCP<Epetra_Vector> dis)
{
  // contact
  if (HaveContact()) ContactManager()->GetStrategy().Update(dis);

  // meshtying
  if (HaveMeshtying()) MtManager()->GetStrategy().Update(dis);

  return;
}

/*----------------------------------------------------------------------*
 |  Write restart                                             popp 11/14|
 *----------------------------------------------------------------------*/
void CONTACT::MeshtyingContactBridge::VisualizeGmsh(const int istep, const int iter)
{
  // contact
  if (HaveContact()) ContactManager()->GetStrategy().VisualizeGmsh(istep, iter);

  // meshtying
  if (HaveMeshtying()) MtManager()->GetStrategy().VisualizeGmsh(istep, iter);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Comm& CONTACT::MeshtyingContactBridge::Comm() const
{
  if (cman_ != Teuchos::null) return cman_->Comm();

  if (mtman_ != Teuchos::null) return mtman_->Comm();

  FOUR_C_THROW("can't get comm()");
  return cman_->Comm();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Mortar::ManagerBase> CONTACT::MeshtyingContactBridge::ContactManager() const
{
  return cman_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Mortar::ManagerBase> CONTACT::MeshtyingContactBridge::MtManager() const
{
  return mtman_;
}

FOUR_C_NAMESPACE_CLOSE
