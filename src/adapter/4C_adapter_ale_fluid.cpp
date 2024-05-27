/*----------------------------------------------------------------------------*/
/*! \file

\brief Wrapper for the ALE time integration for fluid problems with moving boundaries

\level 2

 */
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "4C_adapter_ale_fluid.hpp"

#include "4C_ale_utils_mapextractor.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
ADAPTER::AleFluidWrapper::AleFluidWrapper(Teuchos::RCP<Ale> ale) : AleWrapper(ale)
{
  // create the FSI interface
  interface_ = Teuchos::rcp(new ALE::UTILS::MapExtractor);
  interface_->Setup(*discretization());
  // extend dirichlet map by the dof
  if (interface_->FSICondRelevant())
    SetupDBCMapEx(ALE::UTILS::MapExtractor::dbc_set_part_fsi, interface_);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<const ALE::UTILS::MapExtractor> ADAPTER::AleFluidWrapper::Interface() const
{
  return interface_;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
int ADAPTER::AleFluidWrapper::Solve()
{
  if (interface_->FSICondRelevant())
    Evaluate(Teuchos::null, ALE::UTILS::MapExtractor::dbc_set_part_fsi);
  else
    Evaluate();

  int err = AleWrapper::Solve();
  UpdateIter();

  return err;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ADAPTER::AleFluidWrapper::apply_free_surface_displacements(
    Teuchos::RCP<const Epetra_Vector> fsdisp)
{
  interface_->InsertFSCondVector(fsdisp, WriteAccessDispnp());
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ADAPTER::AleFluidWrapper::apply_ale_update_displacements(
    Teuchos::RCP<const Epetra_Vector> audisp)
{
  interface_->InsertAUCondVector(audisp, WriteAccessDispnp());
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ADAPTER::AleFluidWrapper::apply_interface_displacements(
    Teuchos::RCP<const Epetra_Vector> idisp)
{
  interface_->InsertFSICondVector(idisp, WriteAccessDispnp());
}

FOUR_C_NAMESPACE_CLOSE
