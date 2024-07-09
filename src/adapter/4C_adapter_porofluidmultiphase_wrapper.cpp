/*----------------------------------------------------------------------*/
/*! \file
 \brief a wrapper for porous multiphase flow algorithms

   \level 3

 *----------------------------------------------------------------------*/

#include "4C_adapter_porofluidmultiphase_wrapper.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_validparameters.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_porofluidmultiphase_timint_implicit.hpp"
#include "4C_porofluidmultiphase_timint_ost.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Adapter::PoroFluidMultiphaseWrapper::PoroFluidMultiphaseWrapper(
    Teuchos::RCP<PoroFluidMultiphase> porofluid)
    : porofluid_(porofluid)
{
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::PoroFluidMultiphaseWrapper::init(const bool isale,  ///< ALE flag
    const int nds_disp,           ///< number of dofset associated with displacements
    const int nds_vel,            ///< number of dofset associated with fluid velocities
    const int nds_solidpressure,  ///< number of dofset associated with solid pressure
    const int
        ndsporofluid_scatra,  ///< number of dofset associated with scalar on fluid discretization
    const std::map<int, std::set<int>>* nearbyelepairs  ///< possible interaction partners between
                                                        ///< porofluid and artery discretization
)
{
  // initialize algorithm for specific time-integration scheme
  porofluid_->init(
      isale, nds_disp, nds_vel, nds_solidpressure, ndsporofluid_scatra, nearbyelepairs);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Core::UTILS::ResultTest> Adapter::PoroFluidMultiphaseWrapper::create_field_test()
{
  return porofluid_->create_field_test();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> Adapter::PoroFluidMultiphaseWrapper::dof_row_map(unsigned nds) const
{
  return porofluid_->dof_row_map(nds);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> Adapter::PoroFluidMultiphaseWrapper::artery_dof_row_map() const
{
  return porofluid_->artery_dof_row_map();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase>
Adapter::PoroFluidMultiphaseWrapper::artery_porofluid_sysmat() const
{
  return porofluid_->artery_porofluid_sysmat();
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Core::FE::Discretization> Adapter::PoroFluidMultiphaseWrapper::discretization() const
{
  return porofluid_->discretization();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::PoroFluidMultiphaseWrapper::read_restart(int restart)
{
  return porofluid_->read_restart(restart);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::PoroFluidMultiphaseWrapper::apply_mesh_movement(
    Teuchos::RCP<const Epetra_Vector> dispnp  //!< displacement vector
)
{
  porofluid_->apply_mesh_movement(dispnp);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::PoroFluidMultiphaseWrapper::set_velocity_field(Teuchos::RCP<const Epetra_Vector> vel)
{
  porofluid_->set_velocity_field(vel);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::PoroFluidMultiphaseWrapper::set_state(
    unsigned nds, const std::string& name, Teuchos::RCP<const Epetra_Vector> state)
{
  porofluid_->set_state(nds, name, state);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::PoroFluidMultiphaseWrapper::set_scatra_solution(
    unsigned nds, Teuchos::RCP<const Epetra_Vector> scalars)
{
  set_state(nds, "scalars", scalars);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> Adapter::PoroFluidMultiphaseWrapper::phinp() const
{
  return porofluid_->phinp();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> Adapter::PoroFluidMultiphaseWrapper::phin() const
{
  return porofluid_->phin();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> Adapter::PoroFluidMultiphaseWrapper::solid_pressure() const
{
  return porofluid_->solid_pressure();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> Adapter::PoroFluidMultiphaseWrapper::pressure() const
{
  return porofluid_->pressure();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> Adapter::PoroFluidMultiphaseWrapper::saturation() const
{
  return porofluid_->saturation();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> Adapter::PoroFluidMultiphaseWrapper::valid_vol_frac_spec_dofs()
    const
{
  return porofluid_->valid_vol_frac_spec_dofs();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_MultiVector> Adapter::PoroFluidMultiphaseWrapper::flux() const
{
  return porofluid_->flux();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int Adapter::PoroFluidMultiphaseWrapper::get_dof_set_number_of_solid_pressure() const
{
  return porofluid_->get_dof_set_number_of_solid_pressure();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::PoroFluidMultiphaseWrapper::time_loop() { porofluid_->time_loop(); }

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::PoroFluidMultiphaseWrapper::prepare_time_step() { porofluid_->prepare_time_step(); }

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::PoroFluidMultiphaseWrapper::output() { porofluid_->output(); }

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::PoroFluidMultiphaseWrapper::update() { porofluid_->update(); }

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::PoroFluidMultiphaseWrapper::evaluate_error_compared_to_analytical_sol()
{
  porofluid_->evaluate_error_compared_to_analytical_sol();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::PoroFluidMultiphaseWrapper::solve() { porofluid_->solve(); }

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::PoroFluidMultiphaseWrapper::prepare_time_loop() { porofluid_->prepare_time_loop(); }
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Core::LinAlg::MapExtractor>
Adapter::PoroFluidMultiphaseWrapper::get_dbc_map_extractor() const
{
  return porofluid_->get_dbc_map_extractor();
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> Adapter::PoroFluidMultiphaseWrapper::rhs() const
{
  return porofluid_->rhs();
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> Adapter::PoroFluidMultiphaseWrapper::artery_porofluid_rhs() const
{
  return porofluid_->artery_porofluid_rhs();
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::PoroFluidMultiphaseWrapper::update_iter(const Teuchos::RCP<const Epetra_Vector> inc)
{
  porofluid_->update_iter(inc);
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::PoroFluidMultiphaseWrapper::reconstruct_pressures_and_saturations()
{
  porofluid_->reconstruct_pressures_and_saturations();
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::PoroFluidMultiphaseWrapper::reconstruct_flux() { porofluid_->reconstruct_flux(); }
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::PoroFluidMultiphaseWrapper::calculate_phase_velocities()
{
  porofluid_->calculate_phase_velocities();
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::PoroFluidMultiphaseWrapper::evaluate() { porofluid_->evaluate(); }
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::SparseMatrix> Adapter::PoroFluidMultiphaseWrapper::system_matrix()
{
  return porofluid_->system_matrix();
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::PoroFluidMultiphaseWrapper::assemble_fluid_struct_coupling_mat(
    Teuchos::RCP<Core::LinAlg::SparseOperator> k_fs)
{
  return porofluid_->assemble_fluid_struct_coupling_mat(k_fs);
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::PoroFluidMultiphaseWrapper::assemble_fluid_scatra_coupling_mat(
    Teuchos::RCP<Core::LinAlg::SparseOperator> k_pfs)
{
  return porofluid_->assemble_fluid_scatra_coupling_mat(k_pfs);
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Adapter::ArtNet> Adapter::PoroFluidMultiphaseWrapper::art_net_tim_int()
{
  return porofluid_->art_net_tim_int();
}

FOUR_C_NAMESPACE_CLOSE
