// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ADAPTER_POROFLUIDMULTIPHASE_WRAPPER_HPP
#define FOUR_C_ADAPTER_POROFLUIDMULTIPHASE_WRAPPER_HPP

#include "4C_config.hpp"

#include "4C_adapter_porofluidmultiphase.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace POROFLUIDMULTIPHASE
{
  class TimIntImpl;
}

namespace Adapter
{
  /// basic multiphase porous flow adapter
  class PoroFluidMultiphaseWrapper : public PoroFluidMultiphase
  {
   public:
    /// constructor
    explicit PoroFluidMultiphaseWrapper(std::shared_ptr<PoroFluidMultiphase> porofluid);

    /// initialization
    void init(const bool isale,         ///< ALE flag
        const int nds_disp,             ///< number of dofset associated with displacements
        const int nds_vel,              ///< number of dofset associated with fluid velocities
        const int nds_solidpressure,    ///< number of dofset associated with solid pressure
        const int ndsporofluid_scatra,  ///< number of dofset associated with scalar on fluid
                                        ///< discretization
        const std::map<int, std::set<int>>*
            nearbyelepairs  ///< possible interaction partners between porofluid and artery
                            ///< discretization
        ) override;

    /// create result test for multiphase porous fluid field
    std::shared_ptr<Core::Utils::ResultTest> create_field_test() override;

    /// read restart
    void read_restart(int restart) override;

    /// access dof row map
    std::shared_ptr<const Epetra_Map> dof_row_map(unsigned nds = 0) const override;

    /// access dof row map
    std::shared_ptr<const Epetra_Map> artery_dof_row_map() const override;

    /// access coupled system matrix
    std::shared_ptr<Core::LinAlg::BlockSparseMatrixBase> artery_porofluid_sysmat() const override;

    /// direct access to discretization
    std::shared_ptr<Core::FE::Discretization> discretization() const override;

    //! apply moving mesh data
    void apply_mesh_movement(
        std::shared_ptr<const Core::LinAlg::Vector<double>> dispnp  //!< displacement vector
        ) override;

    //! set state on discretization
    void set_state(unsigned nds, const std::string& name,
        std::shared_ptr<const Core::LinAlg::Vector<double>> state) override;

    //! set convective velocity field (+ pressure and acceleration field as
    //! well as fine-scale velocity field, if required)
    void set_velocity_field(
        std::shared_ptr<const Core::LinAlg::Vector<double>> vel  //!< velocity vector
        ) override;

    //! set solution of scatra problem
    void set_scatra_solution(
        unsigned nds, std::shared_ptr<const Core::LinAlg::Vector<double>> scalars);

    //! return primary field at time n+1
    std::shared_ptr<const Core::LinAlg::Vector<double>> phinp() const override;

    //! return primary field at time n
    std::shared_ptr<const Core::LinAlg::Vector<double>> phin() const override;

    //! return solid pressure field at time n+1
    std::shared_ptr<const Core::LinAlg::Vector<double>> solid_pressure() const override;

    //! return pressure field at time n+1
    std::shared_ptr<const Core::LinAlg::Vector<double>> pressure() const override;

    //! return saturation field at time n+1
    std::shared_ptr<const Core::LinAlg::Vector<double>> saturation() const override;

    //! return valid volume fraction species dof vector
    std::shared_ptr<const Core::LinAlg::Vector<double>> valid_vol_frac_spec_dofs() const override;

    //! return phase flux field at time n+1
    std::shared_ptr<const Core::LinAlg::MultiVector<double>> flux() const override;

    //! return number of dof set associated with solid pressure
    int get_dof_set_number_of_solid_pressure() const override;

    //! do time integration (time loop)
    void time_loop() override;

    //! initialization procedure prior to evaluation of a time step
    void prepare_time_step() override;

    //! output solution and restart data to file
    void output() override;

    //! update the solution after convergence of the nonlinear iteration.
    void update() override;

    //! calculate error compared to analytical solution
    void evaluate_error_compared_to_analytical_sol() override;

    //! general solver call for coupled algorithms
    void solve() override;

    /// prepare timeloop of coupled problem
    void prepare_time_loop() override;

    //! Return MapExtractor for Dirichlet boundary conditions
    std::shared_ptr<const Core::LinAlg::MapExtractor> get_dbc_map_extractor() const override;

    //! right-hand side alias the dynamic force residual
    std::shared_ptr<const Core::LinAlg::Vector<double>> rhs() const override;

    //! right-hand side alias the dynamic force residual for coupled system
    std::shared_ptr<const Core::LinAlg::Vector<double>> artery_porofluid_rhs() const override;

    //! iterative update of phinp
    void update_iter(const std::shared_ptr<const Core::LinAlg::Vector<double>> inc) override;

    //! reconstruct pressures and saturation from current solution
    void reconstruct_pressures_and_saturations() override;

    //! reconstruct flux from current solution
    void reconstruct_flux() override;

    //! calculate phase velocities from current solution
    void calculate_phase_velocities() override;

    //! build linear system tangent matrix, rhs/force residual
    void evaluate() override;

    // Assemble Off-Diagonal Fluid-Structure Coupling matrix
    void assemble_fluid_struct_coupling_mat(
        std::shared_ptr<Core::LinAlg::SparseOperator> k_fs) override;

    // Assemble Off-Diagonal Fluid-Scatra Coupling matrix
    void assemble_fluid_scatra_coupling_mat(
        std::shared_ptr<Core::LinAlg::SparseOperator> k_pfs) override;

    /// direct access to system matrix
    std::shared_ptr<Core::LinAlg::SparseMatrix> system_matrix() override;

    // return arterial network time integrator
    std::shared_ptr<Adapter::ArtNet> art_net_tim_int() override;

   private:
    /// multiphase porous flow time integrator
    std::shared_ptr<PoroFluidMultiphase> porofluid_;

  };  // class PoroFluidMultiphaseWrapper

}  // namespace Adapter



FOUR_C_NAMESPACE_CLOSE

#endif
