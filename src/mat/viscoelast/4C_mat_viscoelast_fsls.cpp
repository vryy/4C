// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_mat_viscoelast_fsls.hpp"

#include "4C_material_parameter_base.hpp"

#include <cmath>

FOUR_C_NAMESPACE_OPEN

Mat::ViscoElast::PAR::Fsls::Fsls(const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      tau_(matdata.parameters.get<double>("TAU")),
      alpha_(matdata.parameters.get<double>("ALPHA")),
      beta_(matdata.parameters.get<double>("BETA"))
{
  if (tau_ <= 0.0)
    FOUR_C_THROW(
        "Invalid TAU={} in VISCO_FSLS (MAT {}). TAU has to be positive.", tau_, matdata.id);

  if (alpha_ < 0.0 || alpha_ >= 1.0)
    FOUR_C_THROW(
        "Invalid ALPHA={} in VISCO_FSLS (MAT {}). Expected 0 <= ALPHA < 1.", alpha_, matdata.id);
}

Mat::ViscoElast::Fsls::Fsls(Mat::ViscoElast::PAR::Fsls* params) : params_(params) {}

void Mat::ViscoElast::Fsls::read_material_parameters_visco(
    double& tau, double& beta, double& alpha, std::string& solve)
{
  tau = params_->tau_;
  alpha = params_->alpha_;
  beta = params_->beta_;
}


void Mat::ViscoElast::Kernels::evaluate_fsls_kernel(const FslsStressVector& stress,
    const FslsTangentMatrix& cmat, FslsStressVector& q_current_for_history,
    FslsStressVector& q_additive, FslsTangentMatrix& cmatq_additive, const FslsKernelInput& input)
{
  if (input.dt < 0.0)
    FOUR_C_THROW(
        "Invalid time step size dt={} in FSLS kernel evaluation (MAT {}, GP {}, ELE {}). "
        "Expected dt >= 0.",
        input.dt, input.visco_mat_id, input.gp, input.ele_gid);

  if (input.tau <= 0.0)
    FOUR_C_THROW(
        "Invalid TAU={} in FSLS kernel evaluation (MAT {}, GP {}, ELE {}). Expected TAU > 0.",
        input.tau, input.visco_mat_id, input.gp, input.ele_gid);

  if (input.alpha < 0.0 || input.alpha >= 1.0)
    FOUR_C_THROW(
        "Invalid ALPHA={} in FSLS kernel evaluation (MAT {}, GP {}, ELE {}). Expected 0 <= "
        "ALPHA < 1.",
        input.alpha, input.visco_mat_id, input.gp, input.ele_gid);

  if (input.previous_history == nullptr)
    FOUR_C_THROW("Missing FSLS history state in kernel evaluation (MAT {}, GP {}, ELE {}).",
        input.visco_mat_id, input.gp, input.ele_gid);

  const FslsHistory& previous_history = *input.previous_history;
  if (previous_history.empty())
    FOUR_C_THROW("Missing FSLS history state in kernel evaluation (MAT {}, GP {}, ELE {}).",
        input.visco_mat_id, input.gp, input.ele_gid);

  if (input.gp < 0 || input.gp >= static_cast<int>(previous_history.size()))
    FOUR_C_THROW(
        "Invalid Gauss point index GP={} for FSLS history in kernel evaluation (MAT {}, ELE "
        "{}). History container size is {}.",
        input.gp, input.visco_mat_id, input.ele_gid, previous_history.size());

  const auto& fsls_history_at_gp = previous_history.at(input.gp);
  const int hs = fsls_history_at_gp.size();
  if (hs <= 0)
    FOUR_C_THROW(
        "Invalid FSLS history size {} at GP {} in kernel evaluation (MAT {}, ELE {}). "
        "Expected at least one entry.",
        hs, input.gp, input.visco_mat_id, input.ele_gid);

  // calculate artificial history stress Qq with weights b_j
  // Qq = sum[j=1 up to j=n][b_j*Q_(n+1-j)] (short: b*Qj)
  // b_j = (j-1-alpha)/j * b_(j-1), with b_0 = 1
  double bj = 1.;
  double fac = 1.;
  FslsStressVector q_history_sum(Core::LinAlg::Initialization::zero);

  for (int j = 1; j <= hs; j++)
  {
    fac = (j - 1. - input.alpha) / j;
    bj = bj * fac;

    FslsStressVector qj(fsls_history_at_gp.at(hs - j));
    q_history_sum.update(bj, qj, 1.0);
  }

  const double dtalpha = std::pow(input.dt, input.alpha);
  const double taualpha = std::pow(input.tau, input.alpha);
  const double denominator = dtalpha + taualpha;
  if (denominator <= 0.0)
    FOUR_C_THROW(
        "Invalid FSLS update denominator dt^alpha + tau^alpha = {} in kernel evaluation (MAT "
        "{}, GP {}, ELE {}): dt={}, tau={}, alpha={}. Expected a positive denominator.",
        denominator, input.visco_mat_id, input.gp, input.ele_gid, input.dt, input.tau, input.alpha);

  const double lambda_1 = dtalpha / denominator;
  const double lambda_2 = -1. * taualpha / denominator;

  q_current_for_history.update(lambda_1 * input.beta, stress, 0.);
  q_current_for_history.update(lambda_2, q_history_sum, 1.);

  q_additive.update(1.0, q_current_for_history, 0.0);
  q_additive.update(input.beta, stress, -1.);

  cmatq_additive.update(lambda_1 * input.beta, cmat, 0.);
  cmatq_additive.update(input.beta, cmat, -1.);
}

FOUR_C_NAMESPACE_CLOSE
