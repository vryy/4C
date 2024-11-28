// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_mat_elast_visco_coupmyocard.hpp"

#include "4C_material_parameter_base.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN


Mat::Elastic::PAR::CoupMyocard::CoupMyocard(const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata), n_(matdata.parameters.get<double>("N"))
{
}

Mat::Elastic::CoupMyocard::CoupMyocard(Mat::Elastic::PAR::CoupMyocard* params) : params_(params) {}

void Mat::Elastic::CoupMyocard::add_coefficients_visco_principal(
    const Core::LinAlg::Matrix<3, 1>& prinv, Core::LinAlg::Matrix<8, 1>& mu,
    Core::LinAlg::Matrix<33, 1>& xi, Core::LinAlg::Matrix<7, 1>& rateinv,
    Teuchos::ParameterList& params, const int gp, const int eleGID)
{
  // material parameter
  const double eta = params_->n_;

  // get time algorithmic parameters.
  const double dt = params.get<double>("delta time");

  // contribution: \dot{C}
  mu(2) = .5 * eta;

  // contribution: id4sharp_{ijkl} = 1/2 (\delta_{ik}\delta_{jl} + \delta_{il}\delta_{jk})
  xi(2) = eta / dt;
}
FOUR_C_NAMESPACE_CLOSE
