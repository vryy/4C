// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_mat_viscoelast_state.hpp"

FOUR_C_NAMESPACE_OPEN

namespace
{
  void validate_kinematic_state(
      const Mat::ViscoElastState::IsoRateState& iso_rate, const char* context)
  {
    if (iso_rate.histscgcurr_ == nullptr || iso_rate.histscglast_ == nullptr ||
        iso_rate.histmodrcgcurr_ == nullptr || iso_rate.histmodrcglast_ == nullptr)
      FOUR_C_THROW("Incomplete kinematic visco history state while {}.", context);

    const std::size_t numgp = iso_rate.histscgcurr_->size();
    if (iso_rate.histscglast_->size() != numgp || iso_rate.histmodrcgcurr_->size() != numgp ||
        iso_rate.histmodrcglast_->size() != numgp)
      FOUR_C_THROW(
          "Inconsistent kinematic visco history sizes while {}: curr={}, last={}, modcurr={}, "
          "modlast={}",
          context, iso_rate.histscgcurr_->size(), iso_rate.histscglast_->size(),
          iso_rate.histmodrcgcurr_->size(), iso_rate.histmodrcglast_->size());
  }
}  // namespace


bool Mat::ViscoElastState::initialized() const
{
  return isinitvis_ && (iso_rate_.histscgcurr_ != nullptr);
}


void Mat::ViscoElastState::set_state_initialized(const bool initialized)
{
  isinitvis_ = initialized;
}


int Mat::ViscoElastState::packed_history_size() const
{
  if (!initialized()) return 0;
  return gauss_point_count();
}


void Mat::ViscoElastState::pack_kinematic_history(
    Core::Communication::PackBuffer& data, const int histsize) const
{
  if (histsize == 0) return;

  validate_kinematic_state(iso_rate_, "packing");
  if (histsize != static_cast<int>(iso_rate_.histscglast_->size()))
    FOUR_C_THROW("Invalid kinematic history size for packing: expected {}, got {}.",
        iso_rate_.histscglast_->size(), histsize);

  for (int gp = 0; gp < histsize; ++gp)
  {
    add_to_pack(data, iso_rate_.histscglast_->at(gp));
    add_to_pack(data, iso_rate_.histmodrcglast_->at(gp));
    add_to_pack(data, iso_rate_.histscgcurr_->at(gp));
    add_to_pack(data, iso_rate_.histmodrcgcurr_->at(gp));
  }
}


void Mat::ViscoElastState::pack_generalized_maxwell_history(
    Core::Communication::PackBuffer& data, const int histsize) const
{
  for (int gp = 0; gp < histsize; ++gp)
  {
    add_to_pack(data, generalized_maxwell_.histbranchstresslast_->at(gp));
    add_to_pack(data, generalized_maxwell_.histbranchelaststresslast_->at(gp));
    add_to_pack(data, generalized_maxwell_.histbranchstresscurr_->at(gp));
    add_to_pack(data, generalized_maxwell_.histbranchelaststresscurr_->at(gp));
  }
}


void Mat::ViscoElastState::pack_fsls_history(Core::Communication::PackBuffer& data) const
{
  add_to_pack(data, (fsls_.histfslsartstresslastall_ != nullptr));
  if (!(int)(fsls_.histfslsartstresslastall_ != nullptr))
    FOUR_C_THROW("Something got wrong with your history data.");

  add_to_pack(data, fsls_.histfslsartstresslastall_->at(0).size());
  for (std::size_t gp = 0; gp < fsls_.histfslsartstresslastall_->size(); ++gp)
    for (std::size_t step = 0; step < fsls_.histfslsartstresslastall_->at(gp).size(); ++step)
      add_to_pack(data, fsls_.histfslsartstresslastall_->at(gp).at(step));
}


void Mat::ViscoElastState::unpack_kinematic_history(
    Core::Communication::UnpackBuffer& buffer, const int histsize)
{
  initialize_for_unpack(histsize);
  for (int gp = 0; gp < histsize; ++gp)
  {
    extract_from_pack(buffer, iso_rate_.histscglast_->at(gp));
    extract_from_pack(buffer, iso_rate_.histmodrcglast_->at(gp));
    extract_from_pack(buffer, iso_rate_.histscgcurr_->at(gp));
    extract_from_pack(buffer, iso_rate_.histmodrcgcurr_->at(gp));
  }
}


void Mat::ViscoElastState::unpack_generalized_maxwell_history(
    Core::Communication::UnpackBuffer& buffer, const int histsize)
{
  initialize_generalized_maxwell_for_unpack(histsize);
  for (int gp = 0; gp < histsize; ++gp)
  {
    extract_from_pack(buffer, generalized_maxwell_.histbranchstresslast_->at(gp));
    extract_from_pack(buffer, generalized_maxwell_.histbranchelaststresslast_->at(gp));
    extract_from_pack(buffer, generalized_maxwell_.histbranchstresscurr_->at(gp));
    extract_from_pack(buffer, generalized_maxwell_.histbranchelaststresscurr_->at(gp));
  }
}


void Mat::ViscoElastState::unpack_fsls_history(
    Core::Communication::UnpackBuffer& buffer, const int histsize)
{
  bool have_historyalldata;
  extract_from_pack(buffer, have_historyalldata);
  if (!have_historyalldata) FOUR_C_THROW("Something got wrong with your history data.");

  std::size_t histfslsartstressall_stepsize;
  extract_from_pack(buffer, histfslsartstressall_stepsize);
  initialize_fsls_for_unpack(histsize, histfslsartstressall_stepsize);
  for (std::size_t gp = 0; gp < static_cast<std::size_t>(histsize); ++gp)
    for (std::size_t step = 0; step < histfslsartstressall_stepsize; ++step)
      extract_from_pack(buffer, fsls_.histfslsartstresslastall_->at(gp).at(step));
}


int Mat::ViscoElastState::gauss_point_count() const
{
  if (iso_rate_.histscglast_ == nullptr) return 0;
  return static_cast<int>(iso_rate_.histscglast_->size());
}


void Mat::ViscoElastState::initialize_kinematic_history(const int numgp)
{
  StressVector idvec(Core::LinAlg::Initialization::zero);
  for (int i = 0; i < 3; ++i) idvec(i) = 1.;

  iso_rate_.histscgcurr_ = std::make_shared<PointHistory>(numgp, idvec);
  iso_rate_.histscglast_ = std::make_shared<PointHistory>(numgp, idvec);
  iso_rate_.histmodrcgcurr_ = std::make_shared<PointHistory>(numgp, idvec);
  iso_rate_.histmodrcglast_ = std::make_shared<PointHistory>(numgp, idvec);
}


void Mat::ViscoElastState::initialize_generalized_maxwell_history(
    const int numgp, const std::size_t numbranch)
{
  const PointHistory empty_branch_values(numbranch);
  generalized_maxwell_.histbranchstresscurr_ =
      std::make_shared<BranchHistory>(numgp, empty_branch_values);
  generalized_maxwell_.histbranchstresslast_ =
      std::make_shared<BranchHistory>(numgp, empty_branch_values);
  generalized_maxwell_.histbranchelaststresscurr_ =
      std::make_shared<BranchHistory>(numgp, empty_branch_values);
  generalized_maxwell_.histbranchelaststresslast_ =
      std::make_shared<BranchHistory>(numgp, empty_branch_values);
}


void Mat::ViscoElastState::initialize_fsls_history(const int numgp)
{
  const StressVector emptyvec(Core::LinAlg::Initialization::zero);
  fsls_.histfslsartstresscurr_ = std::make_shared<PointHistory>(numgp, emptyvec);
  fsls_.histfslsartstresslastall_ = std::make_shared<FslsHistory>(numgp, PointHistory(1, emptyvec));
}


void Mat::ViscoElastState::initialize_for_unpack(const int histsize)
{
  iso_rate_.histscglast_ = std::make_shared<PointHistory>(histsize);
  iso_rate_.histmodrcglast_ = std::make_shared<PointHistory>(histsize);
  iso_rate_.histscgcurr_ = std::make_shared<PointHistory>(histsize);
  iso_rate_.histmodrcgcurr_ = std::make_shared<PointHistory>(histsize);
}


void Mat::ViscoElastState::initialize_generalized_maxwell_for_unpack(const int histsize)
{
  generalized_maxwell_.histbranchstresslast_ = std::make_shared<BranchHistory>(histsize);
  generalized_maxwell_.histbranchelaststresslast_ = std::make_shared<BranchHistory>(histsize);
  generalized_maxwell_.histbranchstresscurr_ = std::make_shared<BranchHistory>(histsize);
  generalized_maxwell_.histbranchelaststresscurr_ = std::make_shared<BranchHistory>(histsize);
}


void Mat::ViscoElastState::initialize_fsls_for_unpack(
    const int histsize, const std::size_t history_size)
{
  fsls_.histfslsartstresscurr_ = std::make_shared<PointHistory>(histsize);
  fsls_.histfslsartstresslastall_ =
      std::make_shared<FslsHistory>(histsize, PointHistory(history_size));
}


void Mat::ViscoElastState::commit_kinematic_history()
{
  validate_kinematic_state(iso_rate_, "updating");
  iso_rate_.histscglast_ = iso_rate_.histscgcurr_;
  iso_rate_.histmodrcglast_ = iso_rate_.histmodrcgcurr_;
}


void Mat::ViscoElastState::append_fsls_history(const unsigned int max_hist)
{
  if (fsls_.histfslsartstresslastall_ == nullptr || fsls_.histfslsartstresscurr_ == nullptr)
    FOUR_C_THROW("Missing FSLS history state in ViscoElastState during update.");

  if (fsls_.histfslsartstresslastall_->size() != fsls_.histfslsartstresscurr_->size())
    FOUR_C_THROW(
        "Inconsistent FSLS history sizes in ViscoElastState during update: {} gauss-point "
        "histories but {} current entries.",
        fsls_.histfslsartstresslastall_->size(), fsls_.histfslsartstresscurr_->size());

  for (int gp = 0; gp < static_cast<int>(fsls_.histfslsartstresslastall_->size()); ++gp)
  {
    fsls_.histfslsartstresslastall_->at(gp).push_back(fsls_.histfslsartstresscurr_->at(gp));

    if (fsls_.histfslsartstresslastall_->at(gp).size() > max_hist)
    {
      PointHistory tmp_vec(++fsls_.histfslsartstresslastall_->at(gp).begin(),
          fsls_.histfslsartstresslastall_->at(gp).end());
      fsls_.histfslsartstresslastall_->at(gp) = tmp_vec;
    }
  }
}


void Mat::ViscoElastState::reset_current_iteration(const int numgp)
{
  StressVector idvec(Core::LinAlg::Initialization::zero);
  for (int i = 0; i < 3; ++i) idvec(i) = 1.;

  const StressVector emptyvec(Core::LinAlg::Initialization::zero);
  iso_rate_.histscgcurr_ = std::make_shared<PointHistory>(numgp, idvec);
  iso_rate_.histmodrcgcurr_ = std::make_shared<PointHistory>(numgp, idvec);
  fsls_.histfslsartstresscurr_ = std::make_shared<PointHistory>(numgp, emptyvec);
}


std::size_t Mat::ViscoElastState::rollover_generalized_maxwell_history(
    const int numgp, const int visco_mat_id)
{
  generalized_maxwell_.histbranchstresslast_ = generalized_maxwell_.histbranchstresscurr_;
  generalized_maxwell_.histbranchelaststresslast_ = generalized_maxwell_.histbranchelaststresscurr_;

  if (generalized_maxwell_.histbranchstresslast_ == nullptr ||
      generalized_maxwell_.histbranchelaststresslast_ == nullptr)
    FOUR_C_THROW(
        "Missing generalized Maxwell history state in MAT_ViscoElastHyper (MAT {}) during "
        "update.",
        visco_mat_id);

  if (generalized_maxwell_.histbranchstresslast_->size() != static_cast<unsigned int>(numgp) ||
      generalized_maxwell_.histbranchelaststresslast_->size() != static_cast<unsigned int>(numgp))
    FOUR_C_THROW(
        "Invalid generalized Maxwell history container size in MAT_ViscoElastHyper (MAT {}) "
        "during update: expected {} gauss points, got {} and {}.",
        visco_mat_id, numgp, generalized_maxwell_.histbranchstresslast_->size(),
        generalized_maxwell_.histbranchelaststresslast_->size());

  const std::size_t numbranch =
      numgp > 0 ? generalized_maxwell_.histbranchstresslast_->at(0).size() : 0;
  if (numbranch == 0)
    FOUR_C_THROW(
        "Invalid generalized Maxwell history in MAT_ViscoElastHyper (MAT {}) during update: "
        "branch history size is zero.",
        visco_mat_id);

  for (int gp = 0; gp < numgp; ++gp)
  {
    if (generalized_maxwell_.histbranchstresslast_->at(gp).size() != numbranch ||
        generalized_maxwell_.histbranchelaststresslast_->at(gp).size() != numbranch)
      FOUR_C_THROW(
          "Inconsistent generalized Maxwell history sizes in MAT_ViscoElastHyper (MAT {}) "
          "during update at GP {}: expected {} branch entries, got {} and {}.",
          visco_mat_id, gp, numbranch, generalized_maxwell_.histbranchstresslast_->at(gp).size(),
          generalized_maxwell_.histbranchelaststresslast_->at(gp).size());
  }

  return numbranch;
}


void Mat::ViscoElastState::reset_generalized_maxwell_current(
    const int numgp, const std::size_t numbranch)
{
  const PointHistory empty_branch_values(numbranch);
  generalized_maxwell_.histbranchstresscurr_ =
      std::make_shared<BranchHistory>(numgp, empty_branch_values);
  generalized_maxwell_.histbranchelaststresscurr_ =
      std::make_shared<BranchHistory>(numgp, empty_branch_values);
}


const Mat::ViscoElastState::StressVector& Mat::ViscoElastState::scg_last_at(const int gp) const
{
  return iso_rate_.histscglast_->at(gp);
}


const Mat::ViscoElastState::StressVector& Mat::ViscoElastState::modrcg_last_at(const int gp) const
{
  return iso_rate_.histmodrcglast_->at(gp);
}


void Mat::ViscoElastState::set_scg_current_at(const int gp, const StressVector& value)
{
  iso_rate_.histscgcurr_->at(gp) = value;
}


void Mat::ViscoElastState::set_modrcg_current_at(const int gp, const StressVector& value)
{
  iso_rate_.histmodrcgcurr_->at(gp) = value;
}


const Mat::ViscoElastState::PointHistory& Mat::ViscoElastState::branch_elastic_stress_last_at(
    const int gp) const
{
  return generalized_maxwell_.histbranchelaststresslast_->at(gp);
}


const Mat::ViscoElastState::PointHistory& Mat::ViscoElastState::branch_stress_last_at(
    const int gp) const
{
  return generalized_maxwell_.histbranchstresslast_->at(gp);
}


void Mat::ViscoElastState::set_branch_elastic_stress_current_at(
    const int gp, const PointHistory& value)
{
  generalized_maxwell_.histbranchelaststresscurr_->at(gp) = value;
}


void Mat::ViscoElastState::set_branch_stress_current_at(const int gp, const PointHistory& value)
{
  generalized_maxwell_.histbranchstresscurr_->at(gp) = value;
}


bool Mat::ViscoElastState::has_fsls_history() const
{
  return fsls_.histfslsartstresslastall_ != nullptr;
}


int Mat::ViscoElastState::fsls_num_gauss_points() const
{
  if (fsls_.histfslsartstresslastall_ == nullptr) return 0;
  return static_cast<int>(fsls_.histfslsartstresslastall_->size());
}


int Mat::ViscoElastState::fsls_history_size_at(const int gp) const
{
  return static_cast<int>(fsls_.histfslsartstresslastall_->at(gp).size());
}


const Mat::ViscoElastState::StressVector& Mat::ViscoElastState::fsls_history_at(
    const int gp, const int step) const
{
  return fsls_.histfslsartstresslastall_->at(gp).at(step);
}


void Mat::ViscoElastState::set_fsls_current_at(const int gp, const StressVector& value)
{
  fsls_.histfslsartstresscurr_->at(gp) = value;
}


void Mat::ViscoElastState::clear()
{
  iso_rate_.histscgcurr_ = nullptr;
  iso_rate_.histscglast_ = nullptr;
  iso_rate_.histmodrcgcurr_ = nullptr;
  iso_rate_.histmodrcglast_ = nullptr;

  generalized_maxwell_.histbranchstresscurr_ = nullptr;
  generalized_maxwell_.histbranchstresslast_ = nullptr;
  generalized_maxwell_.histbranchelaststresscurr_ = nullptr;
  generalized_maxwell_.histbranchelaststresslast_ = nullptr;

  fsls_.histfslsartstresscurr_ = nullptr;
  fsls_.histfslsartstresslastall_ = nullptr;

  isinitvis_ = false;
}

FOUR_C_NAMESPACE_CLOSE
