/*----------------------------------------------------------------------*/
/*! \file

\brief Basis of all algorithms that perform a coupling between Navier-Stokes
       and scalar transport equations including deforming meshes
\level 2


*/
/*----------------------------------------------------------------------*/

#include "4C_adapter_scatra_fluid_ale_coupling_algo.hpp"

#include "4C_ale_utils_mapextractor.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fluid_utils_mapextractor.hpp"
#include "4C_global_data.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Adapter::ScaTraFluidAleCouplingAlgorithm::ScaTraFluidAleCouplingAlgorithm(const Epetra_Comm& comm,
    const Teuchos::ParameterList& prbdyn, const std::string condname,
    const Teuchos::ParameterList& solverparams)
    : ScaTraFluidCouplingAlgorithm(
          comm, prbdyn, true, "scatra", solverparams),  // yes, we need the ALE formulation
      AleBaseAlgorithm(prbdyn,
          Global::Problem::instance()->get_dis("ale")),  // construct ale base algorithm as well
      condname_(condname)
{
  // keep constructor empty
  return;
}


/*----------------------------------------------------------------------*
| Setup                                                     rauch 08/16 |
*----------------------------------------------------------------------*/
void Adapter::ScaTraFluidAleCouplingAlgorithm::init()
{
  // call init() in base class
  Adapter::ScaTraFluidCouplingAlgorithm::init();

  ale_ = Teuchos::rcp_dynamic_cast<AleFluidWrapper>(AleBaseAlgorithm::ale_field(), true);
}


/*----------------------------------------------------------------------*
| Init                                                      rauch 08/16 |
*----------------------------------------------------------------------*/
void Adapter::ScaTraFluidAleCouplingAlgorithm::setup()
{
  // call setup() in base class
  Adapter::ScaTraFluidCouplingAlgorithm::setup();

  const int ndim = Global::Problem::instance()->n_dim();

  // set up couplings
  icoupfa_ = Teuchos::rcp(new Coupling::Adapter::Coupling());
  icoupfa_->setup_condition_coupling(*fluid_field()->discretization(),
      fluid_field()->interface()->fsi_cond_map(), *ale_field()->discretization(),
      ale_field()->interface()->fsi_cond_map(), condname_, ndim);

  // the fluid-ale coupling always matches
  const Epetra_Map* fluidnodemap = fluid_field()->discretization()->node_row_map();
  const Epetra_Map* alenodemap = ale_field()->discretization()->node_row_map();

  coupfa_ = Teuchos::rcp(new Coupling::Adapter::Coupling());
  coupfa_->setup_coupling(*fluid_field()->discretization(), *ale_field()->discretization(),
      *fluidnodemap, *alenodemap, ndim);

  fluid_field()->set_mesh_map(coupfa_->master_dof_map());

  // the ale matrix might be build just once!
  ale_field()->create_system_matrix(ale_field()->interface());

  return;
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::ScaTraFluidAleCouplingAlgorithm::fluid_ale_nonlinear_solve(
    Teuchos::RCP<Core::LinAlg::Vector> idisp, Teuchos::RCP<Core::LinAlg::Vector> ivel,
    bool pseudotransient)
{
  if (idisp != Teuchos::null)
  {
    // if we have values at the interface we need to apply them
    ale_field()->apply_interface_displacements(fluid_to_ale(idisp));
    if (not pseudotransient)
    {
      fluid_field()->apply_interface_velocities(ivel);
    }
  }

  // Note: We do not look for moving ale boundaries (outside the coupling
  // interface) on the fluid side. Thus if you prescribe time variable ale
  // Dirichlet conditions the according fluid Dirichlet conditions will not
  // notice.

  ale_field()->solve();
  Teuchos::RCP<Core::LinAlg::Vector> fluiddisp =
      ale_to_fluid_field(ale_field()->write_access_dispnp());
  fluid_field()->apply_mesh_displacement(fluiddisp);

  // no computation of fluid velocities in case only ScaTra and ALE are to compute
  if (not pseudotransient)
  {
    fluid_field()->solve();
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::Vector> Adapter::ScaTraFluidAleCouplingAlgorithm::ale_to_fluid_field(
    Teuchos::RCP<Core::LinAlg::Vector> iv) const
{
  return coupfa_->slave_to_master(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::Vector> Adapter::ScaTraFluidAleCouplingAlgorithm::ale_to_fluid_field(
    Teuchos::RCP<const Core::LinAlg::Vector> iv) const
{
  return coupfa_->slave_to_master(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::Vector> Adapter::ScaTraFluidAleCouplingAlgorithm::fluid_to_ale(
    Teuchos::RCP<Core::LinAlg::Vector> iv) const
{
  return icoupfa_->master_to_slave(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::Vector> Adapter::ScaTraFluidAleCouplingAlgorithm::fluid_to_ale(
    Teuchos::RCP<const Core::LinAlg::Vector> iv) const
{
  return icoupfa_->master_to_slave(iv);
}

FOUR_C_NAMESPACE_CLOSE
