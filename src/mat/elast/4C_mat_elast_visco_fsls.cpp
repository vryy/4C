// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_mat_elast_visco_fsls.hpp"

#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN

Mat::Elastic::PAR::Fsls::Fsls(const Core::Mat::PAR::Parameter::Data& matdata)
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

Mat::Elastic::Fsls::Fsls(Mat::Elastic::PAR::Fsls* params) : params_(params) {}

void Mat::Elastic::Fsls::read_material_parameters_visco(
    double& tau, double& beta, double& alpha, std::string& solve)
{
  tau = params_->tau_;
  alpha = params_->alpha_;
  beta = params_->beta_;
}
FOUR_C_NAMESPACE_CLOSE
