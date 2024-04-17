/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of the volumetic penalty material as in Roernbauer2008 (student
thesis)

\level 1
*/
/*----------------------------------------------------------------------*/

#include "baci_matelast_volpenalty.hpp"

#include "baci_mat_par_material.hpp"

FOUR_C_NAMESPACE_OPEN


MAT::ELASTIC::PAR::VolPenalty::VolPenalty(const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : Parameter(matdata),
      eps_(*matdata->Get<double>("EPSILON")),
      gam_(*matdata->Get<double>("GAMMA"))
{
  if (eps_ < 0. || gam_ <= 0.)
    dserror("VolPenalty parameters EPSILON and GAMMA have to be greater zero");
}

MAT::ELASTIC::VolPenalty::VolPenalty(MAT::ELASTIC::PAR::VolPenalty* params) : params_(params) {}

void MAT::ELASTIC::VolPenalty::AddStrainEnergy(double& psi, const CORE::LINALG::Matrix<3, 1>& prinv,
    const CORE::LINALG::Matrix<3, 1>& modinv, const CORE::LINALG::Matrix<6, 1>& glstrain,
    const int gp, const int eleGID)
{
  const double eps = params_->eps_;
  const double gam = params_->gam_;

  // strain energy: Psi=\epsilon \left( J^{\gamma} + \frac 1 {J^{\gamma}} -2 \right)
  // add to overall strain energy
  psi += eps * (pow(modinv(2), gam) + pow(modinv(2), -gam) - 2.);
}

void MAT::ELASTIC::VolPenalty::AddDerivativesModified(CORE::LINALG::Matrix<3, 1>& dPmodI,
    CORE::LINALG::Matrix<6, 1>& ddPmodII, const CORE::LINALG::Matrix<3, 1>& modinv, const int gp,
    const int eleGID)
{
  const double eps = params_->eps_;
  const double gam = params_->gam_;

  dPmodI(2) += eps * gam * (pow(modinv(2), gam - 1.) - pow(modinv(2), -gam - 1.));

  ddPmodII(2) +=
      eps * gam * ((gam - 1.) * pow(modinv(2), gam - 2.) + (gam + 1.) * pow(modinv(2), -gam - 2.));
}

void MAT::ELASTIC::VolPenalty::Add3rdVolDeriv(
    const CORE::LINALG::Matrix<3, 1>& modinv, double& d3PsiVolDJ3)
{
  const double eps = params_->eps_;
  const double gam = params_->gam_;
  const double J = modinv(2);
  d3PsiVolDJ3 += eps * (gam * (gam - 1.) * (gam - 2.) * pow(J, gam - 3.) +
                           (-gam) * (-gam - 1.) * (-gam - 2.) * pow(J, -gam - 3.));
}
FOUR_C_NAMESPACE_CLOSE
