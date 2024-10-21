#include "4C_adapter_ale_fluid.hpp"

#include "4C_ale_utils_mapextractor.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Adapter::AleFluidWrapper::AleFluidWrapper(Teuchos::RCP<Ale> ale) : AleWrapper(ale)
{
  // create the FSI interface
  interface_ = Teuchos::make_rcp<ALE::Utils::MapExtractor>();
  interface_->setup(*discretization());
  // extend dirichlet map by the dof
  if (interface_->fsi_cond_relevant())
    setup_dbc_map_ex(ALE::Utils::MapExtractor::dbc_set_part_fsi, interface_);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<const ALE::Utils::MapExtractor> Adapter::AleFluidWrapper::interface() const
{
  return interface_;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
int Adapter::AleFluidWrapper::solve()
{
  if (interface_->fsi_cond_relevant())
    evaluate(Teuchos::null, ALE::Utils::MapExtractor::dbc_set_part_fsi);
  else
    evaluate();

  int err = AleWrapper::solve();
  update_iter();

  return err;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void Adapter::AleFluidWrapper::apply_free_surface_displacements(
    Teuchos::RCP<const Core::LinAlg::Vector<double>> fsdisp)
{
  interface_->insert_fs_cond_vector(*fsdisp, *write_access_dispnp());
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void Adapter::AleFluidWrapper::apply_ale_update_displacements(
    Teuchos::RCP<const Core::LinAlg::Vector<double>> audisp)
{
  interface_->insert_au_cond_vector(*audisp, *write_access_dispnp());
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void Adapter::AleFluidWrapper::apply_interface_displacements(
    Teuchos::RCP<const Core::LinAlg::Vector<double>> idisp)
{
  interface_->insert_fsi_cond_vector(*idisp, *write_access_dispnp());
}

FOUR_C_NAMESPACE_CLOSE
