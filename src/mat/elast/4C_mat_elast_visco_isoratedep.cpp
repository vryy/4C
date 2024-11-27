// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_mat_elast_visco_isoratedep.hpp"

#include "4C_material_parameter_base.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN


Mat::Elastic::PAR::IsoRateDep::IsoRateDep(const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata), n_(matdata.parameters.get<double>("N"))
{
}

Mat::Elastic::IsoRateDep::IsoRateDep(Mat::Elastic::PAR::IsoRateDep* params) : params_(params) {}

void Mat::Elastic::IsoRateDep::add_coefficients_visco_modified(
    const Core::LinAlg::Matrix<3, 1>& modinv, Core::LinAlg::Matrix<8, 1>& modmu,
    Core::LinAlg::Matrix<33, 1>& modxi, Core::LinAlg::Matrix<7, 1>& modrateinv,
    Teuchos::ParameterList& params, const int gp, const int eleGID)
{
  const double n = params_->n_;

  // get time algorithmic parameters.
  double dt = params.get<double>("delta time");

  modmu(1) += 2. * n * modrateinv(1);
  modmu(2) += (2. * n * (modinv(0) - 3.)) / dt;

  modxi(1) += (4. * n) / dt;
  modxi(2) += (4. * n * (modinv(0) - 3.)) / (dt * dt);
}
FOUR_C_NAMESPACE_CLOSE
