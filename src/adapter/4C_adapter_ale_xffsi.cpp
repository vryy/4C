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
  xff_interface_->setup(*discretization());
  setup_dbc_map_ex(ALE::UTILS::MapExtractor::dbc_set_x_ff, interface(), xff_interface_);
  setup_dbc_map_ex(ALE::UTILS::MapExtractor::dbc_set_x_fsi, interface());
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<const Core::LinAlg::MapExtractor> Adapter::AleXFFsiWrapper::get_dbc_map_extractor()
{
  return AleWrapper::get_dbc_map_extractor(ALE::UTILS::MapExtractor::dbc_set_x_ff);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void Adapter::AleXFFsiWrapper::evaluate(Teuchos::RCP<const Epetra_Vector> stepinc)
{
  AleFsiWrapper::evaluate(stepinc, ALE::UTILS::MapExtractor::dbc_set_x_ff);

  // set dispnp_ of xfem dofs to dispn_
  xff_interface_->insert_xfluid_fluid_cond_vector(
      xff_interface_->extract_xfluid_fluid_cond_vector(dispn()), write_access_dispnp());
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
int Adapter::AleXFFsiWrapper::solve()
{
  AleFsiWrapper::evaluate(Teuchos::null, ALE::UTILS::MapExtractor::dbc_set_x_fsi);

  int err = AleFsiWrapper::solve();

  update_iter();

  return err;
}

FOUR_C_NAMESPACE_CLOSE
