// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_mat_elast_visco_generalizedmaxwell.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN

Mat::Elastic::PAR::GeneralizedMaxwell::GeneralizedMaxwell(
    const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      numbranch_(matdata.parameters.get<int>("NUMBRANCH")),
      matids_(matdata.parameters.get<std::vector<int>>("MATIDS")),
      solve_(matdata.parameters.get<std::string>("SOLVE"))
{
  if (numbranch_ <= 0)
    FOUR_C_THROW(
        "Invalid NUMBRANCH={} in VISCO_GeneralizedMaxwell (MAT {}). NUMBRANCH has to be "
        "positive.",
        numbranch_, matdata.id);

  if (static_cast<int>(matids_.size()) != numbranch_)
    FOUR_C_THROW(
        "Invalid VISCO_GeneralizedMaxwell branch declaration in MAT {}: NUMBRANCH={} but "
        "MATIDS has size {}.",
        matdata.id, numbranch_, matids_.size());

  for (int branch_index = 0; branch_index < numbranch_; ++branch_index)
  {
    const int branch_mat = matids_.at(branch_index);
    if (branch_mat <= 0)
      FOUR_C_THROW(
          "Invalid MATIDS[{}]={} in VISCO_GeneralizedMaxwell (MAT {}). Branch material ids "
          "have to be positive.",
          branch_index, branch_mat, matdata.id);
  }

  if (solve_ != "OneStepTheta" && solve_ != "ExponentialTimeDiscretization")
    FOUR_C_THROW(
        "Invalid input for SOLVE='{}' in VISCO_GeneralizedMaxwell (MAT {}). Use "
        "OneStepTheta or ExponentialTimeDiscretization.",
        solve_, matdata.id);
}

Mat::Elastic::GeneralizedMaxwell::GeneralizedMaxwell(Mat::Elastic::PAR::GeneralizedMaxwell* params)
    : params_(params), branchespotsum_(0), branchtau_(0), internalpotsum_(0)
{
  if (params_->numbranch_ <= 0)
    FOUR_C_THROW("Invalid VISCO_GeneralizedMaxwell setup for MAT {}: NUMBRANCH={} is not positive.",
        params_->id(), params_->numbranch_);

  if (static_cast<int>(params_->matids_.size()) != params_->numbranch_)
    FOUR_C_THROW(
        "Invalid VISCO_GeneralizedMaxwell setup for MAT {}: NUMBRANCH={} but MATIDS has size "
        "{}.",
        params_->id(), params_->numbranch_, params_->matids_.size());

  if (Global::Problem::instance()->materials() == nullptr)
    FOUR_C_THROW(
        "Cannot validate VISCO_GeneralizedMaxwell branches for MAT {} because no global material "
        "bundle is available.",
        params_->id());

  if (Global::Problem::instance()->materials()->num() == 0)
    FOUR_C_THROW(
        "Cannot validate VISCO_GeneralizedMaxwell branches for MAT {} because the global "
        "material bundle is empty.",
        params_->id());

  const int probinst = Global::Problem::instance()->materials()->get_read_from_problem();

  // loop over materials of GeneralizedMaxwell (branches)
  std::vector<int>::const_iterator m;
  int branch_index = 0;
  for (m = params_->matids_.begin(); m != params_->matids_.end(); ++m)
  {
    // make sure the summands of the current branch is empty
    internalpotsum_.clear();

    // get parameters of each branch
    const int matid = *m;

    auto* branch_parameter =
        Global::Problem::instance(probinst)->materials()->parameter_by_id(matid);
    if (branch_parameter->type() != Core::Materials::mes_viscobranch)
      FOUR_C_THROW(
          "Invalid branch declaration in VISCO_GeneralizedMaxwell (MAT {}): MATIDS[{}]={} "
          "refers to material type {}. Expected VISCO_GeneralizedMaxwellBranch.",
          params_->id(), branch_index, matid, branch_parameter->type());

    std::shared_ptr<Mat::Elastic::Summand> visco_branch = Mat::Elastic::Summand::factory(matid);
    std::shared_ptr<Mat::Elastic::ViscoBranch> visco_branch_typed =
        std::dynamic_pointer_cast<Mat::Elastic::ViscoBranch>(visco_branch);
    if (visco_branch_typed == nullptr)
      FOUR_C_THROW(
          "Failed to create VISCO_GeneralizedMaxwellBranch summand from MAT {} referenced by "
          "VISCO_GeneralizedMaxwell (MAT {}, MATIDS[{}]).",
          matid, params_->id(), branch_index);

    double tau = -1.0;
    int branchmatid = -1;

    visco_branch_typed->read_material_parameters(tau, branchmatid);

    if (tau <= 0.0)
      FOUR_C_THROW(
          "Invalid TAU={} in VISCO_GeneralizedMaxwellBranch (MAT {}, referenced by "
          "VISCO_GeneralizedMaxwell MAT {}, MATIDS[{}]). TAU has to be positive.",
          tau, matid, params_->id(), branch_index);

    if (branchmatid <= 0)
      FOUR_C_THROW(
          "Invalid MATID={} in VISCO_GeneralizedMaxwellBranch (MAT {}, referenced by "
          "VISCO_GeneralizedMaxwell MAT {}, MATIDS[{}]). MATID has to be positive.",
          branchmatid, matid, params_->id(), branch_index);

    std::shared_ptr<Mat::Elastic::Summand> sum = Mat::Elastic::Summand::factory(branchmatid);
    if (sum == nullptr)
      FOUR_C_THROW(
          "Failed to create branch elasticity summand for MATID {} in "
          "VISCO_GeneralizedMaxwellBranch (MAT {}, referenced by VISCO_GeneralizedMaxwell MAT "
          "{}, MATIDS[{}]).",
          branchmatid, matid, params_->id(), branch_index);

    // write summand in the vector of summands of each branch
    internalpotsum_.push_back(sum);

    // write into vector of summands of the GeneralizedMaxwell material
    branchespotsum_.push_back(internalpotsum_);
    branchtau_.push_back(tau);

    ++branch_index;

  }  // end for-loop over branches
}

void Mat::Elastic::GeneralizedMaxwell::read_material_parameters(
    int& numbranch, const std::vector<int>*& matids, std::string& solve)
{
  numbranch = params_->numbranch_;
  matids = &params_->matids_;
  solve = params_->solve_;
}

// Viscobranch
Mat::Elastic::PAR::ViscoBranch::ViscoBranch(const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      tau_(matdata.parameters.get<double>("TAU")),
      matid_(matdata.parameters.get<int>("MATID"))
{
  if (tau_ <= 0.0)
    FOUR_C_THROW(
        "Invalid TAU={} in VISCO_GeneralizedMaxwellBranch (MAT {}). TAU has to be positive.", tau_,
        matdata.id);

  if (matid_ <= 0)
    FOUR_C_THROW(
        "Invalid MATID={} in VISCO_GeneralizedMaxwellBranch (MAT {}). MATID has to be "
        "positive.",
        matid_, matdata.id);
}

Mat::Elastic::ViscoBranch::ViscoBranch(Mat::Elastic::PAR::ViscoBranch* params) : params_(params) {}

void Mat::Elastic::ViscoBranch::read_material_parameters(double& tau, int& matid)
{
  tau = params_->tau_;
  matid = params_->matid_;
}

FOUR_C_NAMESPACE_CLOSE
