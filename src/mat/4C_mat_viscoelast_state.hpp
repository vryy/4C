// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_MAT_VISCOELAST_STATE_HPP
#define FOUR_C_MAT_VISCOELAST_STATE_HPP

#include "4C_config.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_mat_material_factory.hpp"

#include <cstddef>
#include <memory>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  struct ViscoElastState
  {
    using StressVector = Core::LinAlg::Matrix<NUM_STRESS_3D, 1>;
    using PointHistory = std::vector<StressVector>;
    using BranchHistory = std::vector<PointHistory>;
    using FslsHistory = std::vector<PointHistory>;

    struct IsoRateState
    {
      std::shared_ptr<PointHistory> histscgcurr_;
      std::shared_ptr<PointHistory> histscglast_;
      std::shared_ptr<PointHistory> histmodrcgcurr_;
      std::shared_ptr<PointHistory> histmodrcglast_;
    };

    struct GeneralizedMaxwellState
    {
      std::shared_ptr<BranchHistory> histbranchstresscurr_;
      std::shared_ptr<BranchHistory> histbranchstresslast_;
      std::shared_ptr<BranchHistory> histbranchelaststresscurr_;
      std::shared_ptr<BranchHistory> histbranchelaststresslast_;
    };

    struct FslsState
    {
      std::shared_ptr<PointHistory> histfslsartstresscurr_;
      std::shared_ptr<FslsHistory> histfslsartstresslastall_;
    };

    IsoRateState iso_rate_;
    GeneralizedMaxwellState generalized_maxwell_;
    FslsState fsls_;

    bool isinitvis_ = false;

    [[nodiscard]] bool initialized() const;
    void set_state_initialized(bool initialized);

    [[nodiscard]] int packed_history_size() const;
    void pack_kinematic_history(Core::Communication::PackBuffer& data, int histsize) const;
    void pack_generalized_maxwell_history(
        Core::Communication::PackBuffer& data, int histsize) const;
    void pack_fsls_history(Core::Communication::PackBuffer& data) const;
    void unpack_kinematic_history(Core::Communication::UnpackBuffer& buffer, int histsize);
    void unpack_generalized_maxwell_history(
        Core::Communication::UnpackBuffer& buffer, int histsize);
    void unpack_fsls_history(Core::Communication::UnpackBuffer& buffer, int histsize);

    [[nodiscard]] int gauss_point_count() const;
    void initialize_kinematic_history(int numgp);
    void initialize_generalized_maxwell_history(int numgp, std::size_t numbranch);
    void initialize_fsls_history(int numgp);
    void initialize_for_unpack(int histsize);
    void initialize_generalized_maxwell_for_unpack(int histsize);
    void initialize_fsls_for_unpack(int histsize, std::size_t history_size);

    void commit_kinematic_history();
    void append_fsls_history(unsigned int max_hist);
    void reset_current_iteration(int numgp);
    std::size_t rollover_generalized_maxwell_history(int numgp, int visco_mat_id);
    void reset_generalized_maxwell_current(int numgp, std::size_t numbranch);

    const StressVector& scg_last_at(int gp) const;
    const StressVector& modrcg_last_at(int gp) const;
    void set_scg_current_at(int gp, const StressVector& value);
    void set_modrcg_current_at(int gp, const StressVector& value);

    const PointHistory& branch_elastic_stress_last_at(int gp) const;
    const PointHistory& branch_stress_last_at(int gp) const;
    void set_branch_elastic_stress_current_at(int gp, const PointHistory& value);
    void set_branch_stress_current_at(int gp, const PointHistory& value);

    [[nodiscard]] bool has_fsls_history() const;
    [[nodiscard]] int fsls_num_gauss_points() const;
    [[nodiscard]] int fsls_history_size_at(int gp) const;
    const StressVector& fsls_history_at(int gp, int step) const;
    void set_fsls_current_at(int gp, const StressVector& value);

    void clear();
  };
}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif
