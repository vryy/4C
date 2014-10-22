/*----------------------------------------------------------------------------*/
/*!
 \file ad_ale_fs.cpp
 <pre>
 Maintainer: Raffaela Kruse
            kruse@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
 </pre>
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "ad_ale_fs.H"

#include "../drt_ale_new/ale_utils_mapextractor.H"

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
ADAPTER::AleFSWrapper::AleFSWrapper(Teuchos::RCP<Ale> ale)
  : AleWrapper(ale)
{
  // create the FSI interface
  interface_ = Teuchos::rcp(new ALENEW::UTILS::MapExtractor);
  interface_->Setup(*Discretization());
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<const ALENEW::UTILS::MapExtractor>
ADAPTER::AleFSWrapper::Interface() const
{
  return interface_;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ADAPTER::AleFSWrapper::ApplyFreeSurfaceDisplacements(
    Teuchos::RCP<const Epetra_Vector> fsdisp)
{
  interface_->InsertFSCondVector(fsdisp, WriteAccessDispnp());
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ADAPTER::AleFSWrapper::ApplyAleUpdateDisplacements(
    Teuchos::RCP<const Epetra_Vector> audisp)
{
  interface_->InsertAUCondVector(audisp, WriteAccessDispnp());
}
