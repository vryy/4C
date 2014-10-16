/*----------------------------------------------------------------------------*/
/*!
 \file ad_ale_xffsi.cpp
 <pre>
       Maintainer: Matthias Mayr
       mayr@mhpc.mw.tum.de
       089 - 289-10362
 </pre>
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "ad_ale_xffsi.H"

#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_utils.H"

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
ADAPTER::AleXFFsiWrapper::AleXFFsiWrapper(Teuchos::RCP<Ale> ale)
  : AleFsiWrapper(ale)
{
  // setup the map extractor for both fluids
  xffinterface_ = Teuchos::rcp(new ALENEW::UTILS::XFluidFluidMapExtractor);
  xffinterface_->Setup(*Discretization());

  // mark the fluid-fluid interface dof as Dirichlet values
  if (xffinterface_->XFluidFluidCondRelevant())
  {
    // create the toggle vector for fluid-fluid-Coupling
    Teuchos::RCP<Epetra_Vector> dispnp_xff = LINALG::CreateVector(*xffinterface_->XFluidFluidCondMap(),true);
    dispnp_xff->PutScalar(1.0);
    xfftoggle_ = LINALG::CreateVector(*Discretization()->DofRowMap(),true);
    xffinterface_->InsertXFluidFluidCondVector(dispnp_xff,xfftoggle_);
  }
}

//! ToDo (mayr) move this to XFluidFluid adapter
void ADAPTER::AleXFFsiWrapper::SolveAleXFluidFluidFSI()
{
  // At the beginning of the fluid-fluid-fsi step the xfem-dofs are
  // dirichlet values so that they can not change in the next
  // iterations. After the fsi step we put the ALE FSI-dofs to
  // dirichlet and we solve the ALE again to find the real ALE
  // displacement.

  // turn the toggle vector off
  xfftoggle_->PutScalar(0.0);

  // new toggle vector which is on for the fsi-dofs_
  Teuchos::RCP<Epetra_Vector> dispnp_fsicond = LINALG::CreateVector(*Interface()->FSICondMap(),true);
  dispnp_fsicond->PutScalar(1.0);
  Interface()->InsertFSICondVector(dispnp_fsicond,xfftoggle_);

  CreateSystemMatrix(true);

  Evaluate(dispnp_fsicond);

  Solve();

  // for the next time step, set the xfem dofs to dirichlet values
  Teuchos::RCP<Epetra_Vector> dispnp_xff = LINALG::CreateVector(*xffinterface_->XFluidFluidCondMap(),true);
  dispnp_xff->PutScalar(1.0);
  xfftoggle_->PutScalar(0.0);
  xffinterface_->InsertXFluidFluidCondVector(dispnp_xff,xfftoggle_);
}
