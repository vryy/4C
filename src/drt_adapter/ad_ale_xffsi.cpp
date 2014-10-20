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
  // create the FSI interface
  xff_interface_ = Teuchos::rcp(new ALENEW::UTILS::XFluidFluidMapExtractor);
  xff_interface_->Setup(*Discretization());
  SetupDBCMapEx(ALENEW::UTILS::MapExtractor::dbc_set_x_ff,Interface(),xff_interface_);
  SetupDBCMapEx(ALENEW::UTILS::MapExtractor::dbc_set_x_fsi,Interface());
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
const Teuchos::RCP<const LINALG::MapExtractor> ADAPTER::AleXFFsiWrapper::GetDBCMapExtractor()
{
  return AleWrapper::GetDBCMapExtractor(ALENEW::UTILS::MapExtractor::dbc_set_x_ff);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ADAPTER::AleXFFsiWrapper::Evaluate(
  Teuchos::RCP<const Epetra_Vector> disiterinc ///< step increment such that \f$ x_{n+1}^{k+1} = x_{n}^{converged}+ stepinc \f$
)
{
  AleFsiWrapper::Evaluate(disiterinc,ALENEW::UTILS::MapExtractor::dbc_set_x_ff);
  // set dispnp_ of xfem dofs to dispn_
  xff_interface_->InsertXFluidFluidCondVector(xff_interface_->ExtractXFluidFluidCondVector(Dispn()), WriteAccessDispnp());
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ADAPTER::AleXFFsiWrapper::SolveAleXFluidFluidFSI()
{
  // At the beginning of the fluid-fluid-fsi step the xfem-dofs are
  // dirichlet values so that they can not change in the next
  // iterations. After the fsi step we put the ALE FSI-dofs to
  // dirichlet and we solve the ALE again to find the real ALE
  // displacement.

  AleFsiWrapper::CreateSystemMatrix();

  AleFsiWrapper::Evaluate(Teuchos::null,ALENEW::UTILS::MapExtractor::dbc_set_x_fsi);

  Solve();

  UpdateIter();
}
