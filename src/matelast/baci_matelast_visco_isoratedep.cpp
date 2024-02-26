/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of the isochoric contribution of a viscous rate dependent material law,
modified from Pioletti, 1997

\level 2
*/
/*----------------------------------------------------------------------*/

#include "baci_matelast_visco_isoratedep.hpp"

#include "baci_mat_par_material.hpp"

BACI_NAMESPACE_OPEN


MAT::ELASTIC::PAR::IsoRateDep::IsoRateDep(const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : Parameter(matdata), n_(*matdata->Get<double>("N"))
{
}

MAT::ELASTIC::IsoRateDep::IsoRateDep(MAT::ELASTIC::PAR::IsoRateDep* params) : params_(params) {}

void MAT::ELASTIC::IsoRateDep::AddCoefficientsViscoModified(
    const CORE::LINALG::Matrix<3, 1>& modinv, CORE::LINALG::Matrix<8, 1>& modmu,
    CORE::LINALG::Matrix<33, 1>& modxi, CORE::LINALG::Matrix<7, 1>& modrateinv,
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
BACI_NAMESPACE_CLOSE
