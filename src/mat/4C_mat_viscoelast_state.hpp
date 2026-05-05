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
      PointHistory scg_current_;
      PointHistory scg_previous_;
      PointHistory modrcg_current_;
      PointHistory modrcg_previous_;
    };

    struct GeneralizedMaxwellState
    {
      BranchHistory branch_stress_current_;
      BranchHistory branch_stress_previous_;
      BranchHistory branch_elastic_stress_current_;
      BranchHistory branch_elastic_stress_previous_;
    };

    struct FslsState
    {
      PointHistory artificial_stress_current_;
      FslsHistory artificial_stress_previous_history_;
    };

    struct IsoRatePrevPointState
    {
      const StressVector& scg;
      const StressVector& modrcg;
    };

    struct GeneralizedMaxwellPrevPointState
    {
      const PointHistory& branch_elastic_stress;
      const PointHistory& branch_stress;
    };

    IsoRateState iso_rate_;
    GeneralizedMaxwellState generalized_maxwell_;
    FslsState fsls_;

    bool isinitvis_ = false;
    int gp_count_ = 0;
    bool has_iso_rate_state_ = false;
    bool has_generalized_maxwell_state_ = false;
    bool has_fsls_state_ = false;

    [[nodiscard]] bool initialized() const;

    void serialize_state(Core::Communication::PackBuffer& data, bool has_iso_rate_state,
        bool has_generalized_maxwell_state, bool has_fsls_state) const;
    void deserialize_state(Core::Communication::UnpackBuffer& buffer, bool has_iso_rate_state,
        bool has_generalized_maxwell_state, bool has_fsls_state);

    void initialize_from_setup(int gp_count, bool has_iso_rate_state,
        bool has_generalized_maxwell_state, std::size_t generalized_maxwell_branch_count,
        bool has_fsls_state);

    void advance_time_step(bool has_iso_rate_state, bool has_generalized_maxwell_state,
        bool has_fsls_state, unsigned int fsls_max_history_size, int visco_mat_id);

    IsoRatePrevPointState iso_rate_prev_point(int gp) const;
    void set_iso_rate_current_point(int gp, const StressVector& scg, const StressVector& modrcg);

    GeneralizedMaxwellPrevPointState generalized_maxwell_prev_point(int gp) const;
    void set_generalized_maxwell_current_point(
        int gp, const PointHistory& branch_elastic_stress, const PointHistory& branch_stress);

    const FslsHistory& fsls_previous_history() const;
    void set_fsls_current_artificial_stress(int gp, const StressVector& value);

    void clear();

   private:
    void mark_initialized(bool initialized);

    [[nodiscard]] int serialized_gp_count() const;
    [[nodiscard]] int gp_count() const;

    void serialize_iso_rate_state(Core::Communication::PackBuffer& data, int gp_count) const;
    void serialize_generalized_maxwell_state(
        Core::Communication::PackBuffer& data, int gp_count) const;
    void serialize_fsls_state(Core::Communication::PackBuffer& data) const;
    void deserialize_iso_rate_state(Core::Communication::UnpackBuffer& buffer, int gp_count);
    void deserialize_generalized_maxwell_state(
        Core::Communication::UnpackBuffer& buffer, int gp_count);
    void deserialize_fsls_state(Core::Communication::UnpackBuffer& buffer, int gp_count);

    void initialize_iso_rate_state_for_deserialize(int gp_count);
    void initialize_generalized_maxwell_state_for_deserialize(int gp_count);
    void initialize_fsls_state_for_deserialize(int gp_count, std::size_t history_size);

    void initialize_iso_rate_state(int gp_count);
    void initialize_generalized_maxwell_state(int gp_count, std::size_t branch_count);
    void initialize_fsls_state(int gp_count);

    void commit_iso_rate_to_previous();
    void append_fsls_current_to_previous_history(unsigned int max_history_size);
    void reset_current_state(int gp_count, bool has_iso_rate_state, bool has_fsls_state);
    std::size_t rotate_generalized_maxwell_to_previous(int gp_count, int visco_mat_id);
    void reset_generalized_maxwell_current_state(int gp_count, std::size_t branch_count);
  };
}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif
