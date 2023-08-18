/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of the isochoric contribution of a material to test the isochoric parts of the
Elasthyper-Toolbox.

\level 1
*/
/*----------------------------------------------------------------------*/

#include "baci_matelast_isotestmaterial.H"

#include "baci_mat_par_material.H"

MAT::ELASTIC::PAR::IsoTestMaterial::IsoTestMaterial(const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : Parameter(matdata), c1_(matdata->GetDouble("C1")), c2_(matdata->GetDouble("C2"))
{
}

MAT::ELASTIC::IsoTestMaterial::IsoTestMaterial(MAT::ELASTIC::PAR::IsoTestMaterial* params)
    : params_(params)
{
}

void MAT::ELASTIC::IsoTestMaterial::AddStrainEnergy(double& psi,
    const CORE::LINALG::Matrix<3, 1>& prinv, const CORE::LINALG::Matrix<3, 1>& modinv,
    const CORE::LINALG::Matrix<6, 1>& glstrain, const int gp, const int eleGID)
{
  const double c1 = params_->c1_;
  const double c2 = params_->c2_;
  const double d = c1 + 2. * c2;

  // strain energy: Psi = C (\overline{I}_{\boldsymbol{C}}-3)^2.
  ///   \Psi = C1 (\overline{I}_{\boldsymbol{C}}-3) + 0.5 C1 (\overline{I}_{\boldsymbol{C}}-3)^2
  ///        + C2 (\overline{II}_{\boldsymbol{C}}-3)  + 0.5 C2 (\overline{II}_{\boldsymbol{C}}-3)^2
  ///        + D (\overline{I}_{\boldsymbol{C}}-3) (\overline{II}_{\boldsymbol{C}}-3).

  //  // add to overall strain energy
  psi += c1 * (modinv(0) - 3) + 0.5 * c1 * (modinv(0) - 3) * (modinv(0) - 3) +
         c2 * (modinv(1) - 3) + 0.5 * c2 * (modinv(1) - 3) * (modinv(1) - 3) +
         d * (modinv(0) - 3) * (modinv(1) - 3);
}

void MAT::ELASTIC::IsoTestMaterial::AddDerivativesModified(CORE::LINALG::Matrix<3, 1>& dPmodI,
    CORE::LINALG::Matrix<6, 1>& ddPmodII, const CORE::LINALG::Matrix<3, 1>& modinv, const int gp,
    const int eleGID)
{
  const double c1 = params_->c1_;
  const double c2 = params_->c2_;

  const double d = c1 + 2. * c2;

  dPmodI(0) += c1 + c1 * (modinv(0) - 3.) + d * (modinv(1) - 3.);
  dPmodI(1) += c2 + d * (modinv(0) - 3.) + c2 * (modinv(1) - 3.);

  ddPmodII(0) += c1;
  ddPmodII(1) += c2;
  ddPmodII(5) += d;
}