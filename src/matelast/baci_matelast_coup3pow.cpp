/*----------------------------------------------------------------------*/
/*! \file
\brief Implementationo of a volumetric general power-type material in terms of the Jacobi
determinant

\level 2
*/
/*----------------------------------------------------------------------*/

#include "baci_matelast_coup3pow.hpp"

#include "baci_mat_par_material.hpp"

BACI_NAMESPACE_OPEN

MAT::ELASTIC::PAR::Coup3Pow::Coup3Pow(const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : Parameter(matdata), c_(*matdata->Get<double>("C")), d_(*matdata->Get<int>("D"))
{
}

MAT::ELASTIC::Coup3Pow::Coup3Pow(MAT::ELASTIC::PAR::Coup3Pow* params) : params_(params) {}

void MAT::ELASTIC::Coup3Pow::AddStrainEnergy(double& psi, const CORE::LINALG::Matrix<3, 1>& prinv,
    const CORE::LINALG::Matrix<3, 1>& modinv, const CORE::LINALG::Matrix<6, 1> glstrain,
    const int eleGID)
{
  // material Constants c and beta
  const double c = params_->c_;
  const int d = params_->d_;

  // add to overall strain energy
  psi += c * pow((pow(prinv(2), 1. / 3.) - 1.), d);
}

void MAT::ELASTIC::Coup3Pow::AddDerivativesPrincipal(CORE::LINALG::Matrix<3, 1>& dPI,
    CORE::LINALG::Matrix<6, 1>& ddPII, const CORE::LINALG::Matrix<3, 1>& prinv, const int gp,
    const int eleGID)
{
  const double c = params_->c_;
  const int d = params_->d_;

  // If d<2
  if (d < 2)
    dserror(
        "The Elast_Coup3Pow - material only works for positive integer exponents, which are larger "
        "than two.");

  dPI(2) += 1. / 3. * c * d * pow(prinv(2), -2. / 3.) * pow((pow(prinv(2), 1. / 3.) - 1.), d - 1.);

  if (d == 2)
    ddPII(2) +=
        -2. / 9. * c * d * pow(prinv(2), -5. / 3.) * pow((pow(prinv(2), 1. / 3.) - 1.), d - 1.) +
        1. / 9. * c * d * (d - 1.) * pow(prinv(2), -4. / 3.);
  else
    ddPII(2) +=
        -2. / 9. * c * d * pow(prinv(2), -5. / 3.) * pow((pow(prinv(2), 1. / 3.) - 1.), d - 1.) +
        1. / 9. * c * d * (d - 1.) * pow(prinv(2), -4. / 3.) *
            pow((pow(prinv(2), 1. / 3.) - 1.), d - 2.);
}
BACI_NAMESPACE_CLOSE
