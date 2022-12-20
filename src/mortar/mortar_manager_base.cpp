/*-----------------------------------------------------------------------*/
/*! \file
\brief Abstract base class to control all mortar coupling


\level 2
*/
/*---------------------------------------------------------------------*/

#include "Epetra_SerialComm.h"
#include "mortar_manager_base.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 01/10|
 *----------------------------------------------------------------------*/
MORTAR::ManagerBase::ManagerBase()
{
  //**********************************************************************
  // empty constructor (this is an abstract base class)
  //**********************************************************************
  // Setup of the mortar contact library is done by a derived class. This
  // derived class is specific to the FEM code into which the mortar contact
  // library is meant to be integrated. For BACI this is realized via the
  // CONTACT::ContactManager class! There the following actions are performed:
  //**********************************************************************
  // 1) get problem dimension (2D or 3D)
  // 2) read and check contact input parameters
  // 3) read and check contact boundary conditions
  // 4) build contact interfaces
  //**********************************************************************
  // A similar process also applies to mortar meshtying libraries. Again
  // a specific derived class is needed. For BACI this is realized via the
  // CONTACT::MeshtyingManager class!
  //**********************************************************************

  // create a simple serial communicator
  comm_ = Teuchos::rcp(new Epetra_SerialComm());

  return;
}
