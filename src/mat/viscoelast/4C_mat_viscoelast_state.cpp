// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_mat_viscoelast_state.hpp"

#include <algorithm>

FOUR_C_NAMESPACE_OPEN

namespace
{
  Mat::ViscoElastState::StressVector make_zero_stress()
  {
    return Mat::ViscoElastState::StressVector(Core::LinAlg::Initialization::zero);
  }


  Mat::ViscoElastState::StressVector make_identity_strain_like()
  {
    Mat::ViscoElastState::StressVector identity(Core::LinAlg::Initialization::zero);
    for (int i = 0; i < 3; ++i) identity(i) = 1.;
    return identity;
  }


  void validate_iso_rate_state(
      const Mat::ViscoElastState::IsoRateState& iso_rate, const char* context)
  {
    const std::size_t gp_count = iso_rate.scg_current_.size();
    if (iso_rate.scg_previous_.size() != gp_count || iso_rate.modrcg_current_.size() != gp_count ||
        iso_rate.modrcg_previous_.size() != gp_count)
      FOUR_C_THROW(
          "Inconsistent iso-rate visco state sizes while {}: scg_current={}, scg_previous={}, "
          "modrcg_current={}, modrcg_previous={}",
          context, iso_rate.scg_current_.size(), iso_rate.scg_previous_.size(),
          iso_rate.modrcg_current_.size(), iso_rate.modrcg_previous_.size());
  }


  void validate_generalized_maxwell_state(
      const Mat::ViscoElastState::GeneralizedMaxwellState& generalized_maxwell,
      const std::size_t gp_count, const char* context)
  {
    if (generalized_maxwell.branch_stress_current_.size() != gp_count ||
        generalized_maxwell.branch_stress_previous_.size() != gp_count ||
        generalized_maxwell.branch_elastic_stress_current_.size() != gp_count ||
        generalized_maxwell.branch_elastic_stress_previous_.size() != gp_count)
      FOUR_C_THROW(
          "Inconsistent generalized Maxwell state sizes while {}: branch_stress_current={}, "
          "branch_stress_previous={}, branch_elastic_stress_current={}, "
          "branch_elastic_stress_previous={}, expected={}",
          context, generalized_maxwell.branch_stress_current_.size(),
          generalized_maxwell.branch_stress_previous_.size(),
          generalized_maxwell.branch_elastic_stress_current_.size(),
          generalized_maxwell.branch_elastic_stress_previous_.size(), gp_count);
  }
}  // namespace


bool Mat::ViscoElastState::initialized() const { return isinitvis_; }


void Mat::ViscoElastState::mark_initialized(const bool initialized) { isinitvis_ = initialized; }


int Mat::ViscoElastState::serialized_gp_count() const
{
  if (!initialized()) return 0;
  return gp_count_;
}


void Mat::ViscoElastState::serialize_state(Core::Communication::PackBuffer& data,
    const bool has_iso_rate_state, const bool has_generalized_maxwell_state,
    const bool has_fsls_state) const
{
  const int gp_count = serialized_gp_count();
  add_to_pack(data, gp_count);

  if (has_iso_rate_state) serialize_iso_rate_state(data, gp_count);
  if (has_generalized_maxwell_state) serialize_generalized_maxwell_state(data, gp_count);
  if (has_fsls_state) serialize_fsls_state(data);
}


void Mat::ViscoElastState::deserialize_state(Core::Communication::UnpackBuffer& buffer,
    const bool has_iso_rate_state, const bool has_generalized_maxwell_state,
    const bool has_fsls_state)
{
  has_iso_rate_state_ = has_iso_rate_state;
  has_generalized_maxwell_state_ = has_generalized_maxwell_state;
  has_fsls_state_ = has_fsls_state;

  mark_initialized(true);

  int gp_count = 0;
  extract_from_pack(buffer, gp_count);
  gp_count_ = gp_count;
  if (gp_count == 0) mark_initialized(false);

  if (has_iso_rate_state)
    deserialize_iso_rate_state(buffer, gp_count);
  else
  {
    iso_rate_.scg_previous_.clear();
    iso_rate_.modrcg_previous_.clear();
    iso_rate_.scg_current_.clear();
    iso_rate_.modrcg_current_.clear();
  }

  if (has_generalized_maxwell_state)
    deserialize_generalized_maxwell_state(buffer, gp_count);
  else
  {
    generalized_maxwell_.branch_stress_previous_.clear();
    generalized_maxwell_.branch_elastic_stress_previous_.clear();
    generalized_maxwell_.branch_stress_current_.clear();
    generalized_maxwell_.branch_elastic_stress_current_.clear();
  }

  if (has_fsls_state)
    deserialize_fsls_state(buffer, gp_count);
  else
  {
    fsls_.artificial_stress_current_.clear();
    fsls_.artificial_stress_previous_history_.clear();
  }
}


void Mat::ViscoElastState::initialize_from_setup(const int gp_count, const bool has_iso_rate_state,
    const bool has_generalized_maxwell_state, const std::size_t generalized_maxwell_branch_count,
    const bool has_fsls_state)
{
  gp_count_ = gp_count;
  has_iso_rate_state_ = has_iso_rate_state;
  has_generalized_maxwell_state_ = has_generalized_maxwell_state;
  has_fsls_state_ = has_fsls_state;

  if (has_iso_rate_state)
    initialize_iso_rate_state(gp_count);
  else
  {
    iso_rate_.scg_current_.clear();
    iso_rate_.scg_previous_.clear();
    iso_rate_.modrcg_current_.clear();
    iso_rate_.modrcg_previous_.clear();
  }

  if (has_generalized_maxwell_state)
  {
    if (generalized_maxwell_branch_count == 0)
      FOUR_C_THROW(
          "Invalid setup for generalized Maxwell visco state: branch count must be positive.");
    initialize_generalized_maxwell_state(gp_count, generalized_maxwell_branch_count);
  }
  else
  {
    generalized_maxwell_.branch_stress_current_.clear();
    generalized_maxwell_.branch_stress_previous_.clear();
    generalized_maxwell_.branch_elastic_stress_current_.clear();
    generalized_maxwell_.branch_elastic_stress_previous_.clear();
  }

  if (has_fsls_state)
    initialize_fsls_state(gp_count);
  else
  {
    fsls_.artificial_stress_current_.clear();
    fsls_.artificial_stress_previous_history_.clear();
  }

  mark_initialized(true);
}


void Mat::ViscoElastState::serialize_iso_rate_state(
    Core::Communication::PackBuffer& data, const int gp_count) const
{
  if (gp_count == 0) return;

  validate_iso_rate_state(iso_rate_, "serializing");
  if (gp_count != static_cast<int>(iso_rate_.scg_previous_.size()))
    FOUR_C_THROW("Invalid iso-rate state size for serialization: expected {}, got {}.",
        iso_rate_.scg_previous_.size(), gp_count);

  for (int gp = 0; gp < gp_count; ++gp)
  {
    add_to_pack(data, iso_rate_.scg_previous_.at(gp));
    add_to_pack(data, iso_rate_.modrcg_previous_.at(gp));
    add_to_pack(data, iso_rate_.scg_current_.at(gp));
    add_to_pack(data, iso_rate_.modrcg_current_.at(gp));
  }
}


void Mat::ViscoElastState::serialize_generalized_maxwell_state(
    Core::Communication::PackBuffer& data, const int gp_count) const
{
  if (gp_count == 0) return;

  validate_generalized_maxwell_state(generalized_maxwell_, gp_count, "serializing");
  for (int gp = 0; gp < gp_count; ++gp)
  {
    add_to_pack(data, generalized_maxwell_.branch_stress_previous_.at(gp));
    add_to_pack(data, generalized_maxwell_.branch_elastic_stress_previous_.at(gp));
    add_to_pack(data, generalized_maxwell_.branch_stress_current_.at(gp));
    add_to_pack(data, generalized_maxwell_.branch_elastic_stress_current_.at(gp));
  }
}


void Mat::ViscoElastState::serialize_fsls_state(Core::Communication::PackBuffer& data) const
{
  const bool have_previous_history = !fsls_.artificial_stress_previous_history_.empty();
  add_to_pack(data, have_previous_history);
  if (!have_previous_history) FOUR_C_THROW("Missing FSLS previous history for serialization.");

  add_to_pack(data, fsls_.artificial_stress_previous_history_.at(0).size());
  for (std::size_t gp = 0; gp < fsls_.artificial_stress_previous_history_.size(); ++gp)
    for (std::size_t step = 0; step < fsls_.artificial_stress_previous_history_.at(gp).size();
        ++step)
      add_to_pack(data, fsls_.artificial_stress_previous_history_.at(gp).at(step));
}


void Mat::ViscoElastState::deserialize_iso_rate_state(
    Core::Communication::UnpackBuffer& buffer, const int gp_count)
{
  initialize_iso_rate_state_for_deserialize(gp_count);
  for (int gp = 0; gp < gp_count; ++gp)
  {
    extract_from_pack(buffer, iso_rate_.scg_previous_.at(gp));
    extract_from_pack(buffer, iso_rate_.modrcg_previous_.at(gp));
    extract_from_pack(buffer, iso_rate_.scg_current_.at(gp));
    extract_from_pack(buffer, iso_rate_.modrcg_current_.at(gp));
  }
}


void Mat::ViscoElastState::deserialize_generalized_maxwell_state(
    Core::Communication::UnpackBuffer& buffer, const int gp_count)
{
  initialize_generalized_maxwell_state_for_deserialize(gp_count);
  for (int gp = 0; gp < gp_count; ++gp)
  {
    extract_from_pack(buffer, generalized_maxwell_.branch_stress_previous_.at(gp));
    extract_from_pack(buffer, generalized_maxwell_.branch_elastic_stress_previous_.at(gp));
    extract_from_pack(buffer, generalized_maxwell_.branch_stress_current_.at(gp));
    extract_from_pack(buffer, generalized_maxwell_.branch_elastic_stress_current_.at(gp));
  }
}


void Mat::ViscoElastState::deserialize_fsls_state(
    Core::Communication::UnpackBuffer& buffer, const int gp_count)
{
  bool have_previous_history;
  extract_from_pack(buffer, have_previous_history);
  if (!have_previous_history) FOUR_C_THROW("Missing FSLS previous history in serialized data.");

  std::size_t fsls_previous_history_size;
  extract_from_pack(buffer, fsls_previous_history_size);
  initialize_fsls_state_for_deserialize(gp_count, fsls_previous_history_size);
  for (std::size_t gp = 0; gp < static_cast<std::size_t>(gp_count); ++gp)
    for (std::size_t step = 0; step < fsls_previous_history_size; ++step)
      extract_from_pack(buffer, fsls_.artificial_stress_previous_history_.at(gp).at(step));
}


int Mat::ViscoElastState::gp_count() const { return gp_count_; }


void Mat::ViscoElastState::initialize_iso_rate_state(const int gp_count)
{
  const StressVector identity = make_identity_strain_like();

  iso_rate_.scg_current_.assign(gp_count, identity);
  iso_rate_.scg_previous_.assign(gp_count, identity);
  iso_rate_.modrcg_current_.assign(gp_count, identity);
  iso_rate_.modrcg_previous_.assign(gp_count, identity);
}


void Mat::ViscoElastState::initialize_generalized_maxwell_state(
    const int gp_count, const std::size_t branch_count)
{
  const PointHistory empty_branch_values(branch_count, make_zero_stress());
  generalized_maxwell_.branch_stress_current_.assign(gp_count, empty_branch_values);
  generalized_maxwell_.branch_stress_previous_.assign(gp_count, empty_branch_values);
  generalized_maxwell_.branch_elastic_stress_current_.assign(gp_count, empty_branch_values);
  generalized_maxwell_.branch_elastic_stress_previous_.assign(gp_count, empty_branch_values);
}


void Mat::ViscoElastState::initialize_fsls_state(const int gp_count)
{
  const StressVector zero = make_zero_stress();
  fsls_.artificial_stress_current_.assign(gp_count, zero);
  fsls_.artificial_stress_previous_history_.assign(gp_count, PointHistory(1, zero));
}


void Mat::ViscoElastState::initialize_iso_rate_state_for_deserialize(const int gp_count)
{
  const StressVector zero = make_zero_stress();
  iso_rate_.scg_previous_.assign(gp_count, zero);
  iso_rate_.modrcg_previous_.assign(gp_count, zero);
  iso_rate_.scg_current_.assign(gp_count, zero);
  iso_rate_.modrcg_current_.assign(gp_count, zero);
}


void Mat::ViscoElastState::initialize_generalized_maxwell_state_for_deserialize(const int gp_count)
{
  generalized_maxwell_.branch_stress_previous_.assign(gp_count, PointHistory{});
  generalized_maxwell_.branch_elastic_stress_previous_.assign(gp_count, PointHistory{});
  generalized_maxwell_.branch_stress_current_.assign(gp_count, PointHistory{});
  generalized_maxwell_.branch_elastic_stress_current_.assign(gp_count, PointHistory{});
}


void Mat::ViscoElastState::initialize_fsls_state_for_deserialize(
    const int gp_count, const std::size_t history_size)
{
  const StressVector zero = make_zero_stress();
  fsls_.artificial_stress_current_.assign(gp_count, zero);
  fsls_.artificial_stress_previous_history_.assign(gp_count, PointHistory(history_size, zero));
}


void Mat::ViscoElastState::advance_time_step(const bool has_iso_rate_state,
    const bool has_generalized_maxwell_state, const bool has_fsls_state,
    const unsigned int fsls_max_history_size, const int visco_mat_id)
{
  has_iso_rate_state_ = has_iso_rate_state;
  has_generalized_maxwell_state_ = has_generalized_maxwell_state;
  has_fsls_state_ = has_fsls_state;

  if (has_iso_rate_state) commit_iso_rate_to_previous();
  if (has_fsls_state) append_fsls_current_to_previous_history(fsls_max_history_size);

  const int gp_count = this->gp_count();
  reset_current_state(gp_count, has_iso_rate_state, has_fsls_state);

  if (has_generalized_maxwell_state)
  {
    const std::size_t branch_count = rotate_generalized_maxwell_to_previous(gp_count, visco_mat_id);
    reset_generalized_maxwell_current_state(gp_count, branch_count);
  }
}


void Mat::ViscoElastState::commit_iso_rate_to_previous()
{
  validate_iso_rate_state(iso_rate_, "advancing time step");
  std::swap(iso_rate_.scg_previous_, iso_rate_.scg_current_);
  std::swap(iso_rate_.modrcg_previous_, iso_rate_.modrcg_current_);
}


void Mat::ViscoElastState::append_fsls_current_to_previous_history(
    const unsigned int max_history_size)
{
  if (fsls_.artificial_stress_previous_history_.size() != fsls_.artificial_stress_current_.size())
    FOUR_C_THROW(
        "Inconsistent FSLS state sizes while advancing time step: previous history gauss points "
        "{}, current stress entries {}.",
        fsls_.artificial_stress_previous_history_.size(), fsls_.artificial_stress_current_.size());

  for (int gp = 0; gp < static_cast<int>(fsls_.artificial_stress_previous_history_.size()); ++gp)
  {
    fsls_.artificial_stress_previous_history_.at(gp).push_back(
        fsls_.artificial_stress_current_.at(gp));

    const std::size_t max_history_size_per_gp = max_history_size;
    if (fsls_.artificial_stress_previous_history_.at(gp).size() > max_history_size_per_gp)
    {
      auto& gp_history = fsls_.artificial_stress_previous_history_.at(gp);
      gp_history.erase(gp_history.begin(), gp_history.begin() + 1);
    }
  }
}


void Mat::ViscoElastState::reset_current_state(
    const int gp_count, const bool has_iso_rate_state, const bool has_fsls_state)
{
  const StressVector identity = make_identity_strain_like();
  const StressVector zero = make_zero_stress();

  if (has_iso_rate_state)
  {
    if (iso_rate_.scg_current_.size() != static_cast<std::size_t>(gp_count))
      iso_rate_.scg_current_.assign(gp_count, identity);
    else
      std::fill(iso_rate_.scg_current_.begin(), iso_rate_.scg_current_.end(), identity);

    if (iso_rate_.modrcg_current_.size() != static_cast<std::size_t>(gp_count))
      iso_rate_.modrcg_current_.assign(gp_count, identity);
    else
      std::fill(iso_rate_.modrcg_current_.begin(), iso_rate_.modrcg_current_.end(), identity);
  }

  if (has_fsls_state)
  {
    if (fsls_.artificial_stress_current_.size() != static_cast<std::size_t>(gp_count))
      fsls_.artificial_stress_current_.assign(gp_count, zero);
    else
      std::fill(
          fsls_.artificial_stress_current_.begin(), fsls_.artificial_stress_current_.end(), zero);
  }
}


std::size_t Mat::ViscoElastState::rotate_generalized_maxwell_to_previous(
    const int gp_count, const int visco_mat_id)
{
  std::swap(
      generalized_maxwell_.branch_stress_previous_, generalized_maxwell_.branch_stress_current_);
  std::swap(generalized_maxwell_.branch_elastic_stress_previous_,
      generalized_maxwell_.branch_elastic_stress_current_);

  validate_generalized_maxwell_state(generalized_maxwell_, gp_count, "advancing time step");

  const std::size_t branch_count =
      gp_count > 0 ? generalized_maxwell_.branch_stress_previous_.at(0).size() : 0;
  if (branch_count == 0)
    FOUR_C_THROW(
        "Invalid generalized Maxwell state in MAT_ViscoElastHyper (MAT {}) during time-step "
        "advance: branch count is zero.",
        visco_mat_id);

  for (int gp = 0; gp < gp_count; ++gp)
  {
    if (generalized_maxwell_.branch_stress_previous_.at(gp).size() != branch_count ||
        generalized_maxwell_.branch_elastic_stress_previous_.at(gp).size() != branch_count)
      FOUR_C_THROW(
          "Inconsistent generalized Maxwell branch sizes in MAT_ViscoElastHyper (MAT {}) "
          "during time-step advance at GP {}: expected {}, got branch_stress={} and "
          "branch_elastic_stress={}",
          visco_mat_id, gp, branch_count,
          generalized_maxwell_.branch_stress_previous_.at(gp).size(),
          generalized_maxwell_.branch_elastic_stress_previous_.at(gp).size());
  }

  return branch_count;
}


void Mat::ViscoElastState::reset_generalized_maxwell_current_state(
    const int gp_count, const std::size_t branch_count)
{
  const StressVector zero = make_zero_stress();

  if (generalized_maxwell_.branch_stress_current_.size() != static_cast<std::size_t>(gp_count))
    generalized_maxwell_.branch_stress_current_.assign(
        gp_count, PointHistory(branch_count, make_zero_stress()));
  if (generalized_maxwell_.branch_elastic_stress_current_.size() !=
      static_cast<std::size_t>(gp_count))
    generalized_maxwell_.branch_elastic_stress_current_.assign(
        gp_count, PointHistory(branch_count, make_zero_stress()));

  for (int gp = 0; gp < gp_count; ++gp)
  {
    auto& branch_stress_current = generalized_maxwell_.branch_stress_current_.at(gp);
    auto& branch_elastic_stress_current =
        generalized_maxwell_.branch_elastic_stress_current_.at(gp);

    if (branch_stress_current.size() != branch_count)
      branch_stress_current.assign(branch_count, zero);
    else
      std::fill(branch_stress_current.begin(), branch_stress_current.end(), zero);

    if (branch_elastic_stress_current.size() != branch_count)
      branch_elastic_stress_current.assign(branch_count, zero);
    else
      std::fill(branch_elastic_stress_current.begin(), branch_elastic_stress_current.end(), zero);
  }
}


Mat::ViscoElastState::IsoRatePrevPointState Mat::ViscoElastState::iso_rate_prev_point(
    const int gp) const
{
  if (!has_iso_rate_state_)
    FOUR_C_THROW("Requested iso-rate previous state, but iso-rate visco model is inactive.");

  return {iso_rate_.scg_previous_.at(gp), iso_rate_.modrcg_previous_.at(gp)};
}


void Mat::ViscoElastState::set_iso_rate_current_point(
    const int gp, const StressVector& scg, const StressVector& modrcg)
{
  if (!has_iso_rate_state_)
    FOUR_C_THROW("Attempted to write iso-rate current state, but iso-rate model is inactive.");

  iso_rate_.scg_current_.at(gp) = scg;
  iso_rate_.modrcg_current_.at(gp) = modrcg;
}


Mat::ViscoElastState::GeneralizedMaxwellPrevPointState
Mat::ViscoElastState::generalized_maxwell_prev_point(const int gp) const
{
  if (!has_generalized_maxwell_state_)
    FOUR_C_THROW("Requested generalized Maxwell previous state, but model is inactive.");

  return {generalized_maxwell_.branch_elastic_stress_previous_.at(gp),
      generalized_maxwell_.branch_stress_previous_.at(gp)};
}


void Mat::ViscoElastState::set_generalized_maxwell_current_point(
    const int gp, const PointHistory& branch_elastic_stress, const PointHistory& branch_stress)
{
  if (!has_generalized_maxwell_state_)
    FOUR_C_THROW("Attempted to write generalized Maxwell current state, but model is inactive.");

  generalized_maxwell_.branch_elastic_stress_current_.at(gp) = branch_elastic_stress;
  generalized_maxwell_.branch_stress_current_.at(gp) = branch_stress;
}


const Mat::ViscoElastState::FslsHistory& Mat::ViscoElastState::fsls_previous_history() const
{
  if (!has_fsls_state_)
    FOUR_C_THROW("Requested FSLS previous history, but FSLS model is inactive.");

  return fsls_.artificial_stress_previous_history_;
}


void Mat::ViscoElastState::set_fsls_current_artificial_stress(
    const int gp, const StressVector& value)
{
  if (!has_fsls_state_)
    FOUR_C_THROW("Attempted to write FSLS current state, but FSLS model is inactive.");

  fsls_.artificial_stress_current_.at(gp) = value;
}


void Mat::ViscoElastState::clear()
{
  iso_rate_.scg_current_.clear();
  iso_rate_.scg_previous_.clear();
  iso_rate_.modrcg_current_.clear();
  iso_rate_.modrcg_previous_.clear();

  generalized_maxwell_.branch_stress_current_.clear();
  generalized_maxwell_.branch_stress_previous_.clear();
  generalized_maxwell_.branch_elastic_stress_current_.clear();
  generalized_maxwell_.branch_elastic_stress_previous_.clear();

  fsls_.artificial_stress_current_.clear();
  fsls_.artificial_stress_previous_history_.clear();

  gp_count_ = 0;
  has_iso_rate_state_ = false;
  has_generalized_maxwell_state_ = false;
  has_fsls_state_ = false;
  isinitvis_ = false;
}

FOUR_C_NAMESPACE_CLOSE
