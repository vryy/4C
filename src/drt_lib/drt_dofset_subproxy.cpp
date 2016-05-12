/*----------------------------------------------------------------------*/
/*!
 \file drt_dofset_subproxy.cpp

 \brief Implementation of subproxy for dofset

 <pre>
   \level 1

   \maintainer Anh-Tu Vuong
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *----------------------------------------------------------------------*/

#include "drt_dofset_subproxy.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::DofSetSubProxy::DofSetSubProxy(DofSet* dofset,
                                    Teuchos::RCP<Epetra_Map>& subcolnodes,
                                    Teuchos::RCP<Epetra_Map>& subcoleles)
  : DofSetProxy(dofset),
    subcolnodes_(subcolnodes),
    subcoleles_(subcoleles)
{
  if(subcolnodes_==Teuchos::null or subcoleles_==Teuchos::null)
    dserror("no node or element map provided for DofSetSubProxy");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::DofSetSubProxy::NotifyReset()
{
//  subcolnodes_ = NULL;
//  subcoleles_ = NULL;

  // clear my Teuchos::rcps.
  DofSetProxy::NotifyReset();
}
