/*----------------------------------------------------------------------------*/
/*! \file

\brief Wrapper for the ALE time integration

\level 2


*/
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "4C_adapter_ale_xffsi.hpp"

#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Adapter::AleXFFsiWrapper::AleXFFsiWrapper(Teuchos::RCP<Ale> ale) : AleFsiWrapper(ale)
{
  // create the FSI interface
  xff_interface_ = Teuchos::rcp(new ALE::UTILS::XFluidFluidMapExtractor);
  xff_interface_->Setup(*discretization());
  SetupDBCMapEx(ALE::UTILS::MapExtractor::dbc_set_x_ff, Interface(), xff_interface_);
  SetupDBCMapEx(ALE::UTILS::MapExtractor::dbc_set_x_fsi, Interface());
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<const Core::LinAlg::MapExtractor> Adapter::AleXFFsiWrapper::GetDBCMapExtractor()
{
  return AleWrapper::GetDBCMapExtractor(ALE::UTILS::MapExtractor::dbc_set_x_ff);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void Adapter::AleXFFsiWrapper::Evaluate(Teuchos::RCP<const Epetra_Vector> stepinc)
{
  AleFsiWrapper::Evaluate(stepinc, ALE::UTILS::MapExtractor::dbc_set_x_ff);

  // set dispnp_ of xfem dofs to dispn_
  xff_interface_->InsertXFluidFluidCondVector(
      xff_interface_->ExtractXFluidFluidCondVector(Dispn()), WriteAccessDispnp());
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
int Adapter::AleXFFsiWrapper::Solve()
{
  AleFsiWrapper::Evaluate(Teuchos::null, ALE::UTILS::MapExtractor::dbc_set_x_fsi);

  int err = AleFsiWrapper::Solve();

  UpdateIter();

  return err;
}

FOUR_C_NAMESPACE_CLOSE
