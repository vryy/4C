/*----------------------------------------------------------------------*/
/*!
 \file drt_dofset_subproxy.cpp

 \brief

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *----------------------------------------------------------------------*/

#include "drt_dofset_subproxy.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::DofSetSubProxy::DofSetSubProxy(DofSet* dofset,
                               Teuchos::RCP<const Epetra_Map> subcolnodes,
                               Teuchos::RCP<const Epetra_Map> subcoleles)
  : DofSetProxy(dofset),
    subcolnodes_(subcolnodes),
    subcoleles_(subcoleles)
{
  if(subcolnodes_==Teuchos::null or subcoleles_==Teuchos::null)
    dserror("no node or element map provided for DofSetSubProxy");
}

