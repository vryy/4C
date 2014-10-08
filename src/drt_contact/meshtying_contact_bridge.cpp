/*!----------------------------------------------------------------------
\file meshtying_contact_bridge.cpp

<pre>
-------------------------------------------------------------------------
                        BACI Contact library
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Philipp Farah
            farah@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/
#include "meshtying_contact_bridge.H"
#include "../drt_contact/meshtying_manager.H"
#include "../drt_contact/contact_manager.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_condition.H"
#include "../linalg/linalg_mapextractor.H"
#include "../drt_mortar/mortar_strategy_base.H"
#include "../drt_mortar/mortar_manager_base.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            farah 06/14|
 *----------------------------------------------------------------------*/
CONTACT::MeshtyingContactBridge::MeshtyingContactBridge(
    DRT::Discretization& dis,
    std::vector<DRT::Condition*>& mtcond,
    std::vector<DRT::Condition*>& ccond,
    double alphaf)
 : cman_(Teuchos::null),
   mtman_(Teuchos::null)
{
  bool onlymeshtying       = false;
  bool onlycontact         = false;
  bool meshtyingandcontact = false;

  // check for case
  if(mtcond.size()!=0 and ccond.size()!=0)
    meshtyingandcontact = true;

  if(mtcond.size()!=0 and ccond.size()==0)
    onlymeshtying = true;

  if(mtcond.size()==0 and ccond.size()!=0)
    onlycontact = true;

  // create meshtying and contact manager
  if (onlymeshtying)
  {
    mtman_ = Teuchos::rcp(new CONTACT::MtManager(dis,alphaf));
  }
  else if (onlycontact)
  {
    cman_ = Teuchos::rcp(new CONTACT::CoManager(dis,alphaf));
  }
  else if (meshtyingandcontact)
  {
    mtman_ = Teuchos::rcp(new CONTACT::MtManager(dis,alphaf));
    cman_  = Teuchos::rcp(new CONTACT::CoManager(dis,alphaf));
  }

  return;
}

/*----------------------------------------------------------------------*
 |  StoreDirichletStatus                                     farah 06/14|
 *----------------------------------------------------------------------*/
void CONTACT::MeshtyingContactBridge::StoreDirichletStatus(Teuchos::RCP<LINALG::MapExtractor> dbcmaps)
{
  if(HaveMeshtying())
    MtManager()->GetStrategy().StoreDirichletStatus(dbcmaps);
  if(HaveContact())
    ContactManager()->GetStrategy().StoreDirichletStatus(dbcmaps);

  return;
}

/*----------------------------------------------------------------------*
 |  Set displacement state                                   farah 06/14|
 *----------------------------------------------------------------------*/
void CONTACT::MeshtyingContactBridge::SetState(Teuchos::RCP<Epetra_Vector> zeros)
{
  if(HaveMeshtying())
    MtManager()->GetStrategy().SetState("displacement",zeros);
  if(HaveContact())
    ContactManager()->GetStrategy().SetState("displacement",zeros);

  return;
}

/*----------------------------------------------------------------------*
 |  Get Strategy                                             farah 06/14|
 *----------------------------------------------------------------------*/
MORTAR::StrategyBase& CONTACT::MeshtyingContactBridge::GetStrategy()
{
  // if contact is involved use contact strategy!
  // contact conditions/strategies are dominating the algorithm!
  if(HaveMeshtying() and !HaveContact())
    return MtManager()->GetStrategy();
  else
    return ContactManager()->GetStrategy();
}

/*----------------------------------------------------------------------*
 |  PostprocessTractions                                     farah 06/14|
 *----------------------------------------------------------------------*/
void CONTACT::MeshtyingContactBridge::PostprocessTractions(Teuchos::RCP<IO::DiscretizationWriter>& output)
{
  // contact
  if (HaveContact())
    ContactManager()->PostprocessTractions(*output);

  // meshtying
  if (HaveMeshtying())
    MtManager()->PostprocessTractions(*output);

  return;
}

/*----------------------------------------------------------------------*
 |  Recover lagr. mult and slave displ                       farah 06/14|
 *----------------------------------------------------------------------*/
void CONTACT::MeshtyingContactBridge::Recover(Teuchos::RCP<Epetra_Vector> disi)
{
  // meshtying
  if (HaveMeshtying())
    MtManager()->GetStrategy().Recover(disi);

  // contact
  if (HaveContact())
    ContactManager()->GetStrategy().Recover(disi);

  return;
}

/*----------------------------------------------------------------------*
 |  Recover lagr. mult and slave displ                       farah 06/14|
 *----------------------------------------------------------------------*/
void CONTACT::MeshtyingContactBridge::ReadRestart(IO::DiscretizationReader& reader,
                                                  Teuchos::RCP<Epetra_Vector> dis,
                                                  Teuchos::RCP<Epetra_Vector> zero)
{
  // contact
  if (HaveContact())
    ContactManager()->ReadRestart(reader,dis,zero);

  // meshtying
  if (HaveMeshtying())
    MtManager()->ReadRestart(reader,dis,zero);

  return;
}

/*----------------------------------------------------------------------*
 |  Write restart                                            farah 06/14|
 *----------------------------------------------------------------------*/
void CONTACT::MeshtyingContactBridge::WriteRestart(Teuchos::RCP<IO::DiscretizationWriter>& output,
                                                   bool forcedrestart)
{
  // contact
  if (HaveContact())
    ContactManager()->WriteRestart(*output,forcedrestart);

  // meshtying
  if (HaveMeshtying())
    MtManager()->WriteRestart(*output,forcedrestart);

  return;
}

/*----------------------------------------------------------------------*
 |  Write restart                                            farah 06/14|
 *----------------------------------------------------------------------*/
void CONTACT::MeshtyingContactBridge::Update(int istep, Teuchos::RCP<Epetra_Vector> dis)
{
  // contact
  if (HaveContact())
    ContactManager()->GetStrategy().Update(istep,dis);

  // meshtying
  if (HaveMeshtying())
    MtManager()->GetStrategy().Update(istep,dis);

  return;
}

