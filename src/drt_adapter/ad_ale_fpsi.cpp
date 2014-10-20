/*----------------------------------------------------------------------------*/
/*!
 \file ad_ale_fpsi.cpp

 <pre>
Maintainer: Andreas Rauch
            rauch@lnm.mw.tum.de
            089 - 289-15240
 </pre>
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "ad_ale_fpsi.H"
#include "../drt_lib/drt_dserror.H"

#include "../drt_ale_new/ale_utils_mapextractor.H"
#include "../drt_ale_new/ale.cpp"

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
ADAPTER::AleFpsiWrapper::AleFpsiWrapper(Teuchos::RCP<Ale> ale)
  : AleWrapper(ale)
{
  // create the FSI interface
  interface_ = Teuchos::rcp(new ALENEW::UTILS::MapExtractor);
  interface_->Setup(*Discretization(),true); //create overlapping maps for fpsi problem

  return;
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ADAPTER::AleFpsiWrapper::ApplyInterfaceDisplacements(Teuchos::RCP<Epetra_Vector> idisp)
{
  interface_->InsertFPSICondVector(idisp,WriteAccessDispnp());

  return;

}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ADAPTER::AleFpsiWrapper::ApplyFSIInterfaceDisplacements(Teuchos::RCP<Epetra_Vector> idisp)
{
  interface_->InsertFSICondVector(idisp,WriteAccessDispnp());

  return;

}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<const ALENEW::UTILS::MapExtractor>
ADAPTER::AleFpsiWrapper::Interface() const
{
  return interface_;
}


