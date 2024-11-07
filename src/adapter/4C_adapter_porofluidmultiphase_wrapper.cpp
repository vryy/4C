// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

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
    std::shared_ptr<PoroFluidMultiphase> porofluid)
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
std::shared_ptr<Core::Utils::ResultTest> Adapter::PoroFluidMultiphaseWrapper::create_field_test()
{
  return porofluid_->create_field_test();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<const Epetra_Map> Adapter::PoroFluidMultiphaseWrapper::dof_row_map(
    unsigned nds) const
{
  return porofluid_->dof_row_map(nds);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<const Epetra_Map> Adapter::PoroFluidMultiphaseWrapper::artery_dof_row_map() const
{
  return porofluid_->artery_dof_row_map();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::LinAlg::BlockSparseMatrixBase>
Adapter::PoroFluidMultiphaseWrapper::artery_porofluid_sysmat() const
{
  return porofluid_->artery_porofluid_sysmat();
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::FE::Discretization> Adapter::PoroFluidMultiphaseWrapper::discretization()
    const
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
    std::shared_ptr<const Core::LinAlg::Vector<double>> dispnp  //!< displacement vector
)
{
  porofluid_->apply_mesh_movement(dispnp);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::PoroFluidMultiphaseWrapper::set_velocity_field(
    std::shared_ptr<const Core::LinAlg::Vector<double>> vel)
{
  porofluid_->set_velocity_field(vel);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::PoroFluidMultiphaseWrapper::set_state(unsigned nds, const std::string& name,
    std::shared_ptr<const Core::LinAlg::Vector<double>> state)
{
  porofluid_->set_state(nds, name, state);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::PoroFluidMultiphaseWrapper::set_scatra_solution(
    unsigned nds, std::shared_ptr<const Core::LinAlg::Vector<double>> scalars)
{
  set_state(nds, "scalars", scalars);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<const Core::LinAlg::Vector<double>> Adapter::PoroFluidMultiphaseWrapper::phinp()
    const
{
  return porofluid_->phinp();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<const Core::LinAlg::Vector<double>> Adapter::PoroFluidMultiphaseWrapper::phin()
    const
{
  return porofluid_->phin();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<const Core::LinAlg::Vector<double>>
Adapter::PoroFluidMultiphaseWrapper::solid_pressure() const
{
  return porofluid_->solid_pressure();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<const Core::LinAlg::Vector<double>> Adapter::PoroFluidMultiphaseWrapper::pressure()
    const
{
  return porofluid_->pressure();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<const Core::LinAlg::Vector<double>>
Adapter::PoroFluidMultiphaseWrapper::saturation() const
{
  return porofluid_->saturation();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<const Core::LinAlg::Vector<double>>
Adapter::PoroFluidMultiphaseWrapper::valid_vol_frac_spec_dofs() const
{
  return porofluid_->valid_vol_frac_spec_dofs();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<const Core::LinAlg::MultiVector<double>> Adapter::PoroFluidMultiphaseWrapper::flux()
    const
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
std::shared_ptr<const Core::LinAlg::MapExtractor>
Adapter::PoroFluidMultiphaseWrapper::get_dbc_map_extractor() const
{
  return porofluid_->get_dbc_map_extractor();
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<const Core::LinAlg::Vector<double>> Adapter::PoroFluidMultiphaseWrapper::rhs() const
{
  return porofluid_->rhs();
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<const Core::LinAlg::Vector<double>>
Adapter::PoroFluidMultiphaseWrapper::artery_porofluid_rhs() const
{
  return porofluid_->artery_porofluid_rhs();
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::PoroFluidMultiphaseWrapper::update_iter(
    const std::shared_ptr<const Core::LinAlg::Vector<double>> inc)
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
std::shared_ptr<Core::LinAlg::SparseMatrix> Adapter::PoroFluidMultiphaseWrapper::system_matrix()
{
  return porofluid_->system_matrix();
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::PoroFluidMultiphaseWrapper::assemble_fluid_struct_coupling_mat(
    std::shared_ptr<Core::LinAlg::SparseOperator> k_fs)
{
  return porofluid_->assemble_fluid_struct_coupling_mat(k_fs);
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::PoroFluidMultiphaseWrapper::assemble_fluid_scatra_coupling_mat(
    std::shared_ptr<Core::LinAlg::SparseOperator> k_pfs)
{
  return porofluid_->assemble_fluid_scatra_coupling_mat(k_pfs);
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Adapter::ArtNet> Adapter::PoroFluidMultiphaseWrapper::art_net_tim_int()
{
  return porofluid_->art_net_tim_int();
}

FOUR_C_NAMESPACE_CLOSE
