// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_mat_viscoelast_isoratedep.hpp"

#include "4C_linalg_fixedsizematrix_tensor_products.hpp"
#include "4C_material_parameter_base.hpp"

#include <Teuchos_ParameterList.hpp>

#include <cmath>

FOUR_C_NAMESPACE_OPEN


Mat::ViscoElast::PAR::IsoRateDep::IsoRateDep(const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata), n_(matdata.parameters.get<double>("N"))
{
}

Mat::ViscoElast::IsoRateDep::IsoRateDep(Mat::ViscoElast::PAR::IsoRateDep* params) : params_(params)
{
}

void Mat::ViscoElast::IsoRateDep::add_coefficients_visco_modified(
    const Core::LinAlg::Matrix<3, 1>& modinv, Core::LinAlg::Matrix<8, 1>& modmu,
    Core::LinAlg::Matrix<33, 1>& modxi, Core::LinAlg::Matrix<7, 1>& modrateinv,
    const Teuchos::ParameterList& /*params*/, const double dt, const int gp, const int eleGID)
{
  const double n = params_->n_;

  modmu(1) += 2. * n * modrateinv(1);
  modmu(2) += (2. * n * (modinv(0) - 3.)) / dt;

  modxi(1) += (4. * n) / dt;
  modxi(2) += (4. * n * (modinv(0) - 3.)) / (dt * dt);
}


void Mat::ViscoElast::Kernels::evaluate_mu_xi_kernel(
    const std::vector<std::shared_ptr<Mat::ViscoElast::Summand>>& summands,
    const bool isoprinc_active, const bool isomod_active, Matrix31& prinv, Matrix31& modinv,
    Matrix81& mu, Matrix81& modmu, Matrix331& xi, Matrix331& modxi, Matrix71& rateinv,
    Matrix71& modrateinv, const Teuchos::ParameterList& params, const double dt, const int gp,
    const int ele_gid)
{
  if (isoprinc_active)
  {
    for (const auto& summand : summands)
      summand->add_coefficients_visco_principal(prinv, mu, xi, rateinv, params, dt, gp, ele_gid);
  }

  if (isomod_active)
  {
    for (const auto& summand : summands)
      summand->add_coefficients_visco_modified(
          modinv, modmu, modxi, modrateinv, params, dt, gp, ele_gid);
  }
}


void Mat::ViscoElast::Kernels::evaluate_kin_quant_vis_kernel(const Matrix61& rcg,
    const Matrix61& scg, const Matrix31& prinv, const Matrix61& scg_previous,
    const Matrix61& modrcg_previous, const double dt, Matrix61& modrcg, Matrix61& scgrate,
    Matrix61& modrcgrate, Matrix71& modrateinv, const int visco_mat_id, const int gp)
{
  if (dt <= 0.0)
    FOUR_C_THROW(
        "Invalid time step size dt={} in MAT_ViscoElastHyper (MAT {}, GP {}) for "
        "rate-dependent viscous update. Expected dt > 0.",
        dt, visco_mat_id, gp);

  const double modscale = std::pow(prinv(2), -1. / 3.);
  modrcg.update(modscale, rcg, 0.0);

  scgrate.update(1.0, scg, 0.0);
  scgrate.update(-1.0, scg_previous, 1.0);
  scgrate.scale(1.0 / dt);

  modrcgrate.update(1.0, modrcg, 0.0);
  modrcgrate.update(-1.0, modrcg_previous, 1.0);
  modrcgrate.scale(1.0 / dt);

  modrateinv(1) =
      0.5 * (modrcgrate(0) * modrcgrate(0) + modrcgrate(1) * modrcgrate(1) +
                modrcgrate(2) * modrcgrate(2) + .5 * modrcgrate(3) * modrcgrate(3) +
                .5 * modrcgrate(4) * modrcgrate(4) + .5 * modrcgrate(5) * modrcgrate(5));
}


void Mat::ViscoElast::Kernels::evaluate_iso_visco_principal_kernel(Matrix61& stress, Matrix66& cmat,
    const Matrix81& mu, const Matrix331& xi, const Matrix66& id4sharp, const Matrix61& scgrate)
{
  stress.update(mu(2), scgrate, 1.0);
  cmat.update(xi(2), id4sharp, 1.0);
}


void Mat::ViscoElast::Kernels::evaluate_iso_visco_modified_kernel(Matrix61& stressisomodisovisco,
    Matrix61& stressisomodvolvisco, Matrix66& cmatisomodisovisco, Matrix66& cmatisomodvolvisco,
    const Matrix31& prinv, const Matrix31& modinv, const Matrix81& modmu, const Matrix331& modxi,
    const Matrix61& rcg, const Matrix61& id2, const Matrix61& icg, const Matrix66& id4,
    const Matrix61& modrcgrate)
{
  const double modscale = std::pow(prinv(2), -1. / 3.);

  Matrix61 modstress(Core::LinAlg::Initialization::zero);
  modstress.update(modmu(1), id2);
  modstress.update(modmu(2), modrcgrate, 1.0);

  Matrix66 projection;
  projection.multiply_nt(1. / 3., icg, rcg);
  projection.update(1.0, id4, -1.0);

  stressisomodisovisco.multiply_nn(modscale, projection, modstress, 1.0);

  Matrix66 modcmat(Core::LinAlg::Initialization::zero);
  Matrix66 modcmat2(Core::LinAlg::Initialization::zero);
  modcmat.multiply_nt(modxi(1), id2, modrcgrate);
  modcmat.multiply_nt(modxi(1), modrcgrate, id2, 1.0);
  modcmat.update(modxi(2), id4, 1.0);
  modcmat.scale(std::pow(modinv(2), -4. / 3.));

  modcmat2.multiply_nn(projection, modcmat);
  cmatisomodisovisco.multiply_nt(1.0, modcmat2, projection, 1.0);

  modcmat.clear();
  modcmat.multiply_nt(-1.0 / 3.0, icg, icg);
  Core::LinAlg::FourTensorOperations::add_holzapfel_product(modcmat, icg, 1.0);
  Core::LinAlg::Matrix<1, 1> tracemat;
  tracemat.multiply_tn(2. / 3. * std::pow(modinv(2), -2. / 3.), modstress, rcg);
  cmatisomodisovisco.update(tracemat(0, 0), modcmat, 1.0);
  cmatisomodisovisco.multiply_nt(-2. / 3., icg, stressisomodisovisco, 1.0);
  cmatisomodisovisco.multiply_nt(-2. / 3., stressisomodisovisco, icg, 1.0);

  // no volumetric visco contribution for visco_isoratedep
  stressisomodvolvisco.clear();
  cmatisomodvolvisco.clear();
}

FOUR_C_NAMESPACE_CLOSE
