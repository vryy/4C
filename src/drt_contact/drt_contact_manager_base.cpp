/*!----------------------------------------------------------------------
\file drt_contact_manager_base.cpp
\brief Main class to control all contact

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
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "Epetra_SerialComm.h"
#include "drt_contact_manager_base.H"
#include "drt_cnode.H"
#include "drt_celement.H"
#include "contactdefines.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_inpar/inpar_contact.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 03/08|
 *----------------------------------------------------------------------*/
CONTACT::ManagerBase::ManagerBase()
{
  //**********************************************************************
  // empty constructor
  //**********************************************************************
  // Setup of the contact library has to be done by a derived class. This
  // derived class is specific to the FEM code into which the contact
  // library is meant to be integrated. For BACI this is realized via the
  // CONTACT::Manager class! There the following actions are performed:
  //**********************************************************************
  // 1) get problem dimension (2D or 3D) and store into dim_
  // 2) read and check contact input parameters
  // 3) read and check contact boundary conditions
  // 4) build contact interfaces
  //**********************************************************************

  // create a simple serial communicator
  comm_ = rcp(new Epetra_SerialComm());

  return;
}

#endif  // #ifdef CCADISCRET
