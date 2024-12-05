// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FS3I_FPS3I_PARTITIONED_HPP
#define FOUR_C_FS3I_FPS3I_PARTITIONED_HPP


#include "4C_config.hpp"

#include "4C_fs3i.hpp"

FOUR_C_NAMESPACE_OPEN


// forward declarations
namespace FPSI
{
  class MonolithicPlain;
}

namespace Adapter
{
  class ScaTraBaseAlgorithm;
}  // namespace Adapter

namespace FSI
{
  class Monolithic;

  namespace Utils
  {
    class MatrixRowTransform;
    class MatrixColTransform;
    class MatrixRowColTransform;
  }  // namespace Utils
}  // namespace FSI

namespace Core::LinAlg
{
  class MultiMapExtractor;
  class BlockSparseMatrixBase;
  class SparseMatrix;
  class Solver;
}  // namespace Core::LinAlg


namespace FS3I
{
  class PartFPS3I : public FS3IBase
  {
   public:
    //! constructor of base class for partitioned FPS3I
    PartFPS3I(MPI_Comm comm);

    //! initialize this class
    void init() override;

    //! setup this class
    void setup() override;

    //! time loop to be defined in inherited classes (structure depends on
    //! considered coupling, i.e. one-way or two-way)
    void timeloop() override = 0;

    //! flag whether time loop should be finished
    bool not_finished() { return step_ < numstep_ and time_ <= timemax_; };

    //! read and set fields needed for restart
    void read_restart() override;

    /// redistribute FPS3I interface, if running on parallel
    void redistribute_interface() override;

    //! set-up of FPSI and ScaTra systems
    void setup_system() override;

    //! test results for individual fields
    void test_results(MPI_Comm comm) override;

    //! evaluate ScaTra fields
    void evaluate_scatra_fields() override;

    //! information transfer FPSI -> ScaTra
    void set_fpsi_solution();

    /// set scatra solution on structure field
    void set_struct_scatra_solution();

    //! return communicator
    MPI_Comm get_comm() const { return comm_; }


    /// extract fluid convective and structure convective velocities
    void extract_vel(std::vector<std::shared_ptr<const Core::LinAlg::Vector<double>>>& vel,
        std::vector<std::shared_ptr<const Core::LinAlg::Vector<double>>>& convel);

    /// extract Wall Shear Stresses at the interface
    void extract_wss(std::vector<std::shared_ptr<const Core::LinAlg::Vector<double>>>& wss);

    /// extracts pressures at the interface
    void extract_pressure(
        std::vector<std::shared_ptr<const Core::LinAlg::Vector<double>>>& pressure);

    /// provide velocities from FPSI subproblem for scatra subproblem
    void set_velocity_fields();

    /// provide wall shear stresses from FPSI subproblem for scatra subproblem
    void set_wall_shear_stresses();

    /// provide pressures from FPSI subproblem for scatra subproblem
    void set_pressure_fields();

    /// provide displacements from FPSI subproblem for scatra subproblem
    void set_mesh_disp();

   protected:
    /// fpsi algorithm
    std::shared_ptr<FPSI::MonolithicPlain> fpsi_;


   private:
    /// communication (mainly for screen output)
    MPI_Comm comm_;

    /// scatra field on fluid
    std::shared_ptr<Adapter::ScaTraBaseAlgorithm> fluidscatra_;

    /// scatra field on structure
    std::shared_ptr<Adapter::ScaTraBaseAlgorithm> structscatra_;
  };
}  // namespace FS3I

FOUR_C_NAMESPACE_CLOSE

#endif
