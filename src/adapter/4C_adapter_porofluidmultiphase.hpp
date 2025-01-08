// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ADAPTER_POROFLUIDMULTIPHASE_HPP
#define FOUR_C_ADAPTER_POROFLUIDMULTIPHASE_HPP


#include "4C_config.hpp"

#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_utils_result_test.hpp"

#include <Epetra_Map.h>

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
  // forward declaration
  class ArtNet;

  /// basic multiphase porous flow adapter
  class PoroFluidMultiphase
  {
   public:
    /// constructor
    PoroFluidMultiphase() {};

    /// virtual destructor to support polymorph destruction
    virtual ~PoroFluidMultiphase() = default;

    /// initialization
    virtual void init(const bool isale,  ///< ALE flag
        const int nds_disp,              ///< number of dofset associated with displacements
        const int nds_vel,               ///< number of dofset associated with fluid velocities
        const int nds_solidpressure,     ///< number of dofset associated with solid pressure
        const int ndsporofluid_scatra,   ///< number of dofset associated with scalar on fluid
                                         ///< discretization
        const std::map<int, std::set<int>>*
            nearbyelepairs  ///< possible interaction partners between porofluid and artery
                            ///< discretization
        ) = 0;

    /// create result test for multiphase porous fluid field
    virtual std::shared_ptr<Core::Utils::ResultTest> create_field_test() = 0;

    /// read restart
    virtual void read_restart(int restart) = 0;

    /// access dof row map
    virtual std::shared_ptr<const Epetra_Map> dof_row_map(unsigned nds = 0) const = 0;

    /// access dof row map of artery discretization
    virtual std::shared_ptr<const Epetra_Map> artery_dof_row_map() const = 0;

    /// direct access to discretization
    virtual std::shared_ptr<Core::FE::Discretization> discretization() const = 0;

    //! apply moving mesh data
    virtual void apply_mesh_movement(
        std::shared_ptr<const Core::LinAlg::Vector<double>> dispnp  //!< displacement vector
        ) = 0;

    //! set convective velocity field (+ pressure and acceleration field as
    //! well as fine-scale velocity field, if required)
    virtual void set_velocity_field(
        std::shared_ptr<const Core::LinAlg::Vector<double>> vel  //!< velocity vector
        ) = 0;

    //! set state on discretization
    virtual void set_state(unsigned nds, const std::string& name,
        std::shared_ptr<const Core::LinAlg::Vector<double>> state) = 0;

    //! return primary field at time n+1
    virtual std::shared_ptr<const Core::LinAlg::Vector<double>> phinp() const = 0;

    //! return primary field at time n
    virtual std::shared_ptr<const Core::LinAlg::Vector<double>> phin() const = 0;

    //! return solid pressure field at time n+1
    virtual std::shared_ptr<const Core::LinAlg::Vector<double>> solid_pressure() const = 0;

    //! return pressure field at time n+1
    virtual std::shared_ptr<const Core::LinAlg::Vector<double>> pressure() const = 0;

    //! return saturation field at time n+1
    virtual std::shared_ptr<const Core::LinAlg::Vector<double>> saturation() const = 0;

    //! return valid volume fraction species dof vector
    virtual std::shared_ptr<const Core::LinAlg::Vector<double>> valid_vol_frac_spec_dofs()
        const = 0;

    //! return phase flux field at time n+1
    virtual std::shared_ptr<const Core::LinAlg::MultiVector<double>> flux() const = 0;

    //! do time integration (time loop)
    virtual void time_loop() = 0;

    //! initialization procedure prior to evaluation of a time step
    virtual void prepare_time_step() = 0;

    //! output solution and restart data to file
    virtual void output() = 0;

    //! update the solution after convergence of the nonlinear iteration.
    virtual void update() = 0;

    //! calculate error compared to analytical solution
    virtual void evaluate_error_compared_to_analytical_sol() = 0;

    //! general solver call for coupled algorithms
    virtual void solve() = 0;

    /// prepare timeloop of coupled problem
    virtual void prepare_time_loop() = 0;

    //! return number of dof set associated with solid pressure
    virtual int get_dof_set_number_of_solid_pressure() const = 0;

    //! Return MapExtractor for Dirichlet boundary conditions
    virtual std::shared_ptr<const Core::LinAlg::MapExtractor> get_dbc_map_extractor() const = 0;

    //! right-hand side alias the dynamic force residual
    virtual std::shared_ptr<const Core::LinAlg::Vector<double>> rhs() const = 0;

    //! right-hand side alias the dynamic force residual for coupled system
    virtual std::shared_ptr<const Core::LinAlg::Vector<double>> artery_porofluid_rhs() const = 0;

    //! iterative update of phinp
    virtual void update_iter(const std::shared_ptr<const Core::LinAlg::Vector<double>> inc) = 0;

    //! reconstruct pressures and saturation from current solution
    virtual void reconstruct_pressures_and_saturations() = 0;

    //! reconstruct flux from current solution
    virtual void reconstruct_flux() = 0;

    //! calculate phase velocities from current solution
    virtual void calculate_phase_velocities() = 0;

    //! build linear system tangent matrix, rhs/force residual
    virtual void evaluate() = 0;

    // Assemble Off-Diagonal Fluid-Structure Coupling matrix
    virtual void assemble_fluid_struct_coupling_mat(
        std::shared_ptr<Core::LinAlg::SparseOperator> k_fs) = 0;

    // Assemble Off-Diagonal Fluid-scatra Coupling matrix
    virtual void assemble_fluid_scatra_coupling_mat(
        std::shared_ptr<Core::LinAlg::SparseOperator> k_pfs) = 0;

    //! direct access to system matrix
    virtual std::shared_ptr<Core::LinAlg::SparseMatrix> system_matrix() = 0;

    //! direct access to block system matrix of artery poro problem
    virtual std::shared_ptr<Core::LinAlg::BlockSparseMatrixBase> artery_porofluid_sysmat()
        const = 0;

    // return arterial network time integrator
    virtual std::shared_ptr<Adapter::ArtNet> art_net_tim_int() = 0;


  };  // class PoroFluidMultiphase

}  // namespace Adapter


FOUR_C_NAMESPACE_CLOSE

#endif
