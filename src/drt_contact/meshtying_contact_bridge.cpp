/*---------------------------------------------------------------------*/
/*!
\file meshtying_contact_bridge.cpp

\brief Deprecated wrapper for the mesh-tying and contact managers.

\level 2

\maintainer Alexander Popp

*/
/*---------------------------------------------------------------------*/

#include "meshtying_contact_bridge.H"
#include "../drt_contact/meshtying_manager.H"
#include "../drt_contact/contact_manager.H"
#include "../drt_contact/smoothing_manager.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_condition.H"

#include "../linalg/linalg_mapextractor.H"

#include "../drt_mortar/mortar_strategy_base.H"
#include "../drt_mortar/mortar_manager_base.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            farah 06/14|
 *----------------------------------------------------------------------*/
CONTACT::MeshtyingContactBridge::MeshtyingContactBridge(DRT::Discretization& dis,
    std::vector<DRT::Condition*>& mtcond, std::vector<DRT::Condition*>& ccond, double alphaf,
    bool smoothing)
    : cman_(Teuchos::null), mtman_(Teuchos::null)
{
  bool onlymeshtying = false;
  bool onlycontact = false;
  bool meshtyingandcontact = false;

  // check for case
  if (mtcond.size() != 0 and ccond.size() != 0) meshtyingandcontact = true;

  if (mtcond.size() != 0 and ccond.size() == 0) onlymeshtying = true;

  if (mtcond.size() == 0 and ccond.size() != 0) onlycontact = true;

  if (smoothing and !meshtyingandcontact)
    dserror(
        "ERROR: Interface smoothing with additional discr. only possible with meshtying and "
        "contact conditions!");

  // create meshtying and contact manager
  if (onlymeshtying)
  {
    mtman_ = Teuchos::rcp(new CONTACT::MtManager(dis, alphaf));
  }
  else if (onlycontact)
  {
    cman_ = Teuchos::rcp(new CONTACT::CoManager(dis, alphaf));
  }
  else if (meshtyingandcontact)
  {
    if (!smoothing)
    {
      mtman_ = Teuchos::rcp(new CONTACT::MtManager(dis, alphaf));
      cman_ = Teuchos::rcp(new CONTACT::CoManager(dis, alphaf));
    }
    else
    {
      sman_ = Teuchos::rcp(new CONTACT::SmoothingManager(dis, alphaf));
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  StoreDirichletStatus                                     farah 06/14|
 *----------------------------------------------------------------------*/
void CONTACT::MeshtyingContactBridge::StoreDirichletStatus(
    Teuchos::RCP<LINALG::MapExtractor> dbcmaps)
{
  if (HaveSmoothing())
  {
    SManager()->GetStrategy().StoreDirichletStatus(dbcmaps);
  }
  else
  {
    if (HaveMeshtying()) MtManager()->GetStrategy().StoreDirichletStatus(dbcmaps);
    if (HaveContact()) ContactManager()->GetStrategy().StoreDirichletStatus(dbcmaps);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Set displacement state                                   farah 06/14|
 *----------------------------------------------------------------------*/
void CONTACT::MeshtyingContactBridge::SetState(Teuchos::RCP<Epetra_Vector> zeros)
{
  if (HaveSmoothing())
  {
    SManager()->GetStrategy().SetState(MORTAR::state_new_displacement, *zeros);
  }
  else
  {
    if (HaveMeshtying())
      MtManager()->GetStrategy().SetState(MORTAR::state_new_displacement, *zeros);
    if (HaveContact())
      ContactManager()->GetStrategy().SetState(MORTAR::state_new_displacement, *zeros);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Get Strategy                                             farah 06/14|
 *----------------------------------------------------------------------*/
MORTAR::StrategyBase& CONTACT::MeshtyingContactBridge::GetStrategy()
{
  if (HaveSmoothing())
  {
    return SManager()->GetStrategy();
  }
  else
  {
    // if contact is involved use contact strategy!
    // contact conditions/strategies are dominating the algorithm!
    if (HaveMeshtying() and !HaveContact())
      return MtManager()->GetStrategy();
    else
      return ContactManager()->GetStrategy();
  }
}

/*----------------------------------------------------------------------*
 |  PostprocessTractions                                     farah 06/14|
 *----------------------------------------------------------------------*/
void CONTACT::MeshtyingContactBridge::PostprocessQuantities(
    Teuchos::RCP<IO::DiscretizationWriter>& output)
{
  // contact
  if (HaveContact()) ContactManager()->PostprocessQuantities(*output);

  // meshtying
  if (HaveMeshtying()) MtManager()->PostprocessQuantities(*output);

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
void CONTACT::MeshtyingContactBridge::ReadRestart(IO::DiscretizationReader& reader,
    Teuchos::RCP<Epetra_Vector> dis, Teuchos::RCP<Epetra_Vector> zero)
{
  // contact
  if (HaveContact()) ContactManager()->ReadRestart(reader, dis, zero);

  // meshtying
  if (HaveMeshtying()) MtManager()->ReadRestart(reader, dis, zero);

  return;
}

/*----------------------------------------------------------------------*
 |  Write restart                                            farah 06/14|
 *----------------------------------------------------------------------*/
void CONTACT::MeshtyingContactBridge::WriteRestart(
    Teuchos::RCP<IO::DiscretizationWriter>& output, bool forcedrestart)
{
  // contact
  if (HaveContact()) ContactManager()->WriteRestart(*output, forcedrestart);

  // meshtying
  if (HaveMeshtying()) MtManager()->WriteRestart(*output, forcedrestart);

  return;
}

/*----------------------------------------------------------------------*
 |  Write restart                                            farah 06/14|
 *----------------------------------------------------------------------*/
void CONTACT::MeshtyingContactBridge::Update(Teuchos::RCP<Epetra_Vector> dis)
{
  if (HaveSmoothing())
  {
    SManager()->GetStrategy().Update(dis);
  }
  else
  {
    // contact
    if (HaveContact()) ContactManager()->GetStrategy().Update(dis);

    // meshtying
    if (HaveMeshtying()) MtManager()->GetStrategy().Update(dis);
  }


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
