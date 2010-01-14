/*!----------------------------------------------------------------------
\file meshtying_interface.cpp
\brief One meshtying interface

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

#ifndef PARALLEL
#include "Epetra_SerialComm.h"
#endif
#include "meshtying_interface.H"
#include "../drt_inpar/inpar_mortar.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 10/07|
 *----------------------------------------------------------------------*/
CONTACT::MtInterface::MtInterface(const int id, const Epetra_Comm& comm,
                                  const int dim,
                                  const Teuchos::ParameterList& icontact) :
MORTAR::MortarInterface(id,comm,dim,icontact)
{
  // empty constructor

  return;
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 10/07|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const CONTACT::MtInterface& interface)
{
  interface.Print(os);
  return os;
}


/*----------------------------------------------------------------------*
 |  print interface (public)                                 mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::MtInterface::Print(ostream& os) const
{
  if (Comm().MyPID()==0)
    os << "Meshyting ";
  MORTAR::MortarInterface::Print(os);

  return;
}

#endif  // #ifdef CCADISCRET
