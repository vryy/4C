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
  /**
   * \brief Per-Gauss-point viscoelastic history state for MAT_ViscoElastHyper.
   *
   * This struct owns all time-dependent internal variables that must be stored,
   * advanced between previous/current time levels, and serialized for MPI/restart.
   */
  struct ViscoElastState
  {
    using StressVector = Core::LinAlg::Matrix<NUM_STRESS_3D, 1>;
    using PointHistory = std::vector<StressVector>;
    using BranchHistory = std::vector<PointHistory>;
    using FslsHistory = std::vector<PointHistory>;

    /// \brief Previous/current kinematic history for the iso-rate visco model.
    struct IsoRateState
    {
      /// Right Cauchy-Green tensor \f$C\f$ at current time level, per GP.
      PointHistory scg_current_;
      /// Right Cauchy-Green tensor \f$C\f$ at previous converged time level, per GP.
      PointHistory scg_previous_;
      /// Modified right Cauchy-Green tensor \f$\bar{C}\f$ at current time level, per GP.
      PointHistory modrcg_current_;
      /// Modified right Cauchy-Green tensor \f$\bar{C}\f$ at previous time level, per GP.
      PointHistory modrcg_previous_;
    };

    /// \brief Previous/current branch stress history for generalized Maxwell visco model.
    struct GeneralizedMaxwellState
    {
      /// Viscous branch stress at current time level, indexed by GP and branch.
      BranchHistory branch_stress_current_;
      /// Viscous branch stress at previous time level, indexed by GP and branch.
      BranchHistory branch_stress_previous_;
      /// Elastic branch stress at current time level, indexed by GP and branch.
      BranchHistory branch_elastic_stress_current_;
      /// Elastic branch stress at previous time level, indexed by GP and branch.
      BranchHistory branch_elastic_stress_previous_;
    };

    /// \brief Previous/current history for FSLS artificial stress state.
    struct FslsState
    {
      /// Artificial FSLS stress at current time level, per GP.
      PointHistory artificial_stress_current_;
      /// Full converged FSLS history per GP (time history of artificial stress).
      FslsHistory artificial_stress_previous_history_;
    };

    /// \brief Read-only per-GP view of previous iso-rate kinematic state.
    struct IsoRatePrevPointState
    {
      const StressVector& scg;
      const StressVector& modrcg;
    };

    /// \brief Read-only per-GP view of previous generalized Maxwell branch state.
    struct GeneralizedMaxwellPrevPointState
    {
      const PointHistory& branch_elastic_stress;
      const PointHistory& branch_stress;
    };

    /// Iso-rate model state containers.
    IsoRateState iso_rate_;
    /// Generalized Maxwell model state containers.
    GeneralizedMaxwellState generalized_maxwell_;
    /// FSLS model state containers.
    FslsState fsls_;

    /// True once state was initialized via setup or deserialization.
    bool isinitvis_ = false;
    /// Number of Gauss points represented by this state instance.
    int gp_count_ = 0;
    /// Activation flag for iso-rate visco model state.
    bool has_iso_rate_state_ = false;
    /// Activation flag for generalized Maxwell model state.
    bool has_generalized_maxwell_state_ = false;
    /// Activation flag for FSLS model state.
    bool has_fsls_state_ = false;

    /// \brief Returns whether the state has been initialized.
    [[nodiscard]] bool initialized() const;

    /// \brief Serialize active model histories to communication buffer.
    void serialize_state(Core::Communication::PackBuffer& data, bool has_iso_rate_state,
        bool has_generalized_maxwell_state, bool has_fsls_state) const;
    /// \brief Deserialize active model histories from communication buffer.
    void deserialize_state(Core::Communication::UnpackBuffer& buffer, bool has_iso_rate_state,
        bool has_generalized_maxwell_state, bool has_fsls_state);

    /// \brief Initialize state containers from setup metadata.
    void initialize_from_setup(int gp_count, bool has_iso_rate_state,
        bool has_generalized_maxwell_state, std::size_t generalized_maxwell_branch_count,
        bool has_fsls_state);

    /// \brief Advance active histories from current to previous and reset current state.
    void advance_time_step(bool has_iso_rate_state, bool has_generalized_maxwell_state,
        bool has_fsls_state, unsigned int fsls_max_history_size, int visco_mat_id);

    /// \brief Return read-only previous iso-rate state for one Gauss point.
    IsoRatePrevPointState iso_rate_prev_point(int gp) const;
    /// \brief Store current iso-rate state for one Gauss point.
    void set_iso_rate_current_point(int gp, const StressVector& scg, const StressVector& modrcg);

    /// \brief Return read-only previous generalized Maxwell state for one Gauss point.
    GeneralizedMaxwellPrevPointState generalized_maxwell_prev_point(int gp) const;
    /// \brief Store current generalized Maxwell state for one Gauss point.
    void set_generalized_maxwell_current_point(
        int gp, const PointHistory& branch_elastic_stress, const PointHistory& branch_stress);

    /// \brief Return full converged FSLS history (per GP and time step).
    const FslsHistory& fsls_previous_history() const;
    /// \brief Store current FSLS artificial stress for one Gauss point.
    void set_fsls_current_artificial_stress(int gp, const StressVector& value);

    /// \brief Reset all state containers and activation metadata.
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
