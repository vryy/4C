/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid field adapter for embedded (ALE-)fluid-fluid problems using XFEM

\level 2


*/
/*----------------------------------------------------------------------*/

#include "4C_adapter_fld_fluid_fluid_fsi.hpp"

#include "4C_fluid_utils_mapextractor.hpp"
#include "4C_fluid_xfluid_fluid.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Adapter::FluidFluidFSI::FluidFluidFSI(Teuchos::RCP<Fluid> xfluidfluid, Teuchos::RCP<Fluid> embfluid,
    Teuchos::RCP<Core::LinAlg::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
    bool isale, bool dirichletcond)
    : FluidFSI(xfluidfluid, embfluid->discretization(), solver, params,
          embfluid->discretization()->writer(), isale, dirichletcond)
{
  // cast fluid to XFluidFluid
  xfluidfluid_ = Teuchos::rcp_dynamic_cast<FLD::XFluidFluid>(xfluidfluid);
  if (xfluidfluid_ == Teuchos::null)
    FOUR_C_THROW("Failed to cast Adapter::Fluid to FLD::XFluidFluid.");
  fluidimpl_ = Teuchos::rcp_dynamic_cast<FLD::FluidImplicitTimeInt>(embfluid);
  if (fluidimpl_ == Teuchos::null)
    FOUR_C_THROW("Failed to cast Adapter::Fluid to FLD::FluidImplicitTimInt.");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FluidFluidFSI::init()
{
  // determine the type of monolithic approach
  const Teuchos::ParameterList& xfluiddyn = params_->sublist("XFLUID DYNAMIC/GENERAL");
  monolithic_approach_ = Core::UTILS::IntegralValue<Inpar::XFEM::MonolithicXffsiApproach>(
      xfluiddyn, "MONOLITHIC_XFFSI_APPROACH");

  // should ALE-relaxation be carried out?
  relaxing_ale_ = (bool)Core::UTILS::IntegralValue<int>(xfluiddyn, "RELAXING_ALE");

  // get no. of timesteps, after which ALE-mesh should be relaxed
  relaxing_ale_every_ = xfluiddyn.get<int>("RELAXING_ALE_EVERY");

  if (!relaxing_ale_ && relaxing_ale_every_ != 0)
    FOUR_C_THROW("You don't want to relax the ALE but provide a relaxation interval != 0 ?!");

  if (relaxing_ale_every_ < 0)
    FOUR_C_THROW(
        "Please provide a reasonable relaxation interval. We can't travel back in time yet.");

  // create map extractor for combined fluid domains
  // (to distinguish between FSI interface DOF / merged inner embedded & background fluid DOF)
  mergedfluidinterface_ = Teuchos::rcp(new FLD::UTILS::MapExtractor());
  // call base class init
  FluidFSI::init();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FluidFluidFSI::prepare_time_step()
{
  if (interface()->fsi_cond_relevant() &&
      (monolithic_approach_ == Inpar::XFEM::XFFSI_FixedALE_Partitioned ||
          monolithic_approach_ == Inpar::XFEM::XFFSI_FixedALE_Interpolation))
  {
    xfluidfluid_->set_interface_fixed();
  }
  else
  {
    xfluidfluid_->set_interface_free();
  }
  xfluidfluid_->prepare_time_step();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> Adapter::FluidFluidFSI::dof_row_map()
{
  return xfluidfluid_->dof_row_map();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FluidFluidFSI::solve()
{
  // cut and do XFEM time integration, solve
  xfluidfluid_->solve();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FluidFluidFSI::update()
{
  if (interface()->fsi_cond_relevant() && is_ale_relaxation_step(step()) &&
      (monolithic_approach_ == Inpar::XFEM::XFFSI_FixedALE_Partitioned ||
          monolithic_approach_ == Inpar::XFEM::XFFSI_FixedALE_Interpolation))
  {
    // allow new interface position
    xfluidfluid_->set_interface_free();

    // cut with new interface location and do XFEM time integration
    xfluidfluid_->prepare_xfem_solve();

    // fix interface position again
    xfluidfluid_->set_interface_fixed();

    if (monolithic_approach_ == Inpar::XFEM::XFFSI_FixedALE_Partitioned)
      xfluidfluid_->update_monolithic_fluid_solution(FluidFSI::interface()->fsi_cond_map());

    if (monolithic_approach_ == Inpar::XFEM::XFFSI_FixedALE_Interpolation)
      xfluidfluid_->interpolate_embedded_state_vectors();

    // refresh the merged fluid map extractor
    setup_interface();

    // create new extended shape derivatives matrix
    prepare_shape_derivatives();
  }

  FluidWrapper::update();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<FLD::UTILS::XFluidFluidMapExtractor> const&
Adapter::FluidFluidFSI::x_fluid_fluid_map_extractor()
{
  return xfluidfluid_->x_fluid_fluid_map_extractor();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FluidFluidFSI::apply_mesh_displacement(Teuchos::RCP<const Epetra_Vector> fluiddisp)
{
  // store old state
  Teuchos::RCP<const Epetra_Vector> disp = meshmap_->extract_cond_vector(fluidimpl_->dispnp());
  meshmap_->insert_cond_vector(disp, xfluidfluid_->write_access_disp_old_state());
  // apply mesh displacement and update grid velocity
  FluidFSI::apply_mesh_displacement(fluiddisp);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> Adapter::FluidFluidFSI::block_system_matrix()
{
  if (mergedfluidinterface_ == Teuchos::null)
    FOUR_C_THROW(
        "Uninitialized map FSI/inner fluid map extractor! Failed to create fluid block matrix.");

  // Create a local copy of the inner & conditioned map
  // Reason: the matrix splitting method from Core::LINALG expects non-const maps
  Teuchos::RCP<Epetra_Map> innermap =
      Teuchos::rcp(new Epetra_Map(*mergedfluidinterface_->other_map()));
  Teuchos::RCP<Epetra_Map> condmap =
      Teuchos::rcp(new Epetra_Map(*mergedfluidinterface_->fsi_cond_map()));
  return xfluidfluid_->block_system_matrix(innermap, condmap);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Adapter::FluidFluidFSI::evaluate(
    Teuchos::RCP<const Epetra_Vector> stepinc  ///< solution increment between time step n and n+1
)
{
  if (monolithic_approach_ == Inpar::XFEM::XFFSI_Full_Newton)
    *xfluidfluid_->write_access_disp_old_state() = *fluidimpl_->dispnp();

  // call the usual routine
  xfluidfluid_->evaluate(stepinc);

  // for fixed ALE approach, we only refresh the global fluid map extractor in update()
  if (monolithic_approach_ != Inpar::XFEM::XFFSI_Full_Newton) return;

  // this is the case of a full Newton approach: update the map extractor, as fluid DOFs possibly
  // have changed!
  setup_interface();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> Adapter::FluidFluidFSI::grid_vel()
{
  return fluidimpl_->grid_vel();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Adapter::FluidFluidFSI::write_access_grid_vel()
{
  return fluidimpl_->write_access_grid_vel();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> Adapter::FluidFluidFSI::dispnp() { return fluidimpl_->dispnp(); }

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Adapter::FluidFluidFSI::write_access_dispnp()
{
  return fluidimpl_->write_access_dispnp();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> Adapter::FluidFluidFSI::dispn() { return fluidimpl_->dispn(); }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Teuchos::RCP<Core::FE::Discretization>& Adapter::FluidFluidFSI::discretization()
{
  return fluidimpl_->discretization();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> Adapter::FluidFluidFSI::velocity_row_map()
{
  // in case of fsi with fluidsplit, return the embedded velocity DOF
  // (to understand the motivation behind this, have a look at the recovery of the
  // Lagrange multiplier in standard ALE-FSI class (fluidsplit) in case of active
  // shape derivatives)
  std::vector<Teuchos::RCP<const Epetra_Map>> maps;
  maps.push_back(xfluidfluid_->x_fluid_fluid_map_extractor()->fluid_map());
  maps.push_back(xfluidfluid_->velocity_row_map());
  Teuchos::RCP<const Epetra_Map> innervelocitymap =
      Core::LinAlg::MultiMapExtractor::intersect_maps(maps);
  return innervelocitymap;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Adapter::FluidFluidFSI::use_block_matrix(bool split_fluidsysmat)
{
  prepare_shape_derivatives();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool Adapter::FluidFluidFSI::is_ale_relaxation_step(int step) const
{
  return relaxing_ale_ && step % relaxing_ale_every_ == 0;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> Adapter::FluidFluidFSI::shape_derivatives()
{
  return xfluidfluid_->extended_shape_derivatives();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Adapter::FluidFluidFSI::setup_interface(const int nds_master)
{
  // check nds_master
  if (nds_master != 0) FOUR_C_THROW("nds_master is supposed to be 0 here");

  if (mergedfluidinterface_ == Teuchos::null)
  {
    std::stringstream errmsg;
    errmsg
        << "Uninitialized map  map extractor for merged background & embedded inner/FSI fluid DOFs."
        << "\nFailed to perform map extractor setup.";
    FOUR_C_THROW(errmsg.str());
  }

  FluidFSI::setup_interface();

  // get background fluid map
  Teuchos::RCP<const Epetra_Map> xfluidmap =
      xfluidfluid_->x_fluid_fluid_map_extractor()->x_fluid_map();
  // do the setup
  mergedfluidinterface_->setup(xfluidmap, *FluidFSI::interface());

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FluidFluidFSI::prepare_shape_derivatives()
{
  // the dof-maps may have changed: create a new shape derivatives matrix
  Teuchos::RCP<std::set<int>> condelements =
      mergedfluidinterface_->conditioned_element_map(*fluidimpl_->discretization()());
  xfluidfluid_->prepare_shape_derivatives(mergedfluidinterface_, condelements);
}

FOUR_C_NAMESPACE_CLOSE
