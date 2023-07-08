/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of the isochoric contribution of a Mooney-Rivlin-type material

\level 1
*/
/*----------------------------------------------------------------------*/

#include "matelast_isomooneyrivlin.H"
#include "mat_par_material.H"


MAT::ELASTIC::PAR::IsoMooneyRivlin::IsoMooneyRivlin(const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : Parameter(matdata), c1_(matdata->GetDouble("C1")), c2_(matdata->GetDouble("C2"))
{
}

MAT::ELASTIC::IsoMooneyRivlin::IsoMooneyRivlin(MAT::ELASTIC::PAR::IsoMooneyRivlin* params)
    : params_(params)
{
}

void MAT::ELASTIC::IsoMooneyRivlin::AddStrainEnergy(double& psi,
    const CORE::LINALG::Matrix<3, 1>& prinv, const CORE::LINALG::Matrix<3, 1>& modinv,
    const CORE::LINALG::Matrix<6, 1>& glstrain, const int gp, const int eleGID)
{
  const double c1 = params_->c1_;
  const double c2 = params_->c2_;

  // strain energy: Psi = C1 (\overline{I}_{\boldsymbol{C}}-3) + C2
  // (\overline{II}_{\boldsymbol{C}}-3). add to overall strain energy
  psi += c1 * (modinv(0) - 3.) + c2 * (modinv(1) - 3.);
}

void MAT::ELASTIC::IsoMooneyRivlin::AddDerivativesModified(CORE::LINALG::Matrix<3, 1>& dPmodI,
    CORE::LINALG::Matrix<6, 1>& ddPmodII, const CORE::LINALG::Matrix<3, 1>& modinv, const int gp,
    const int eleGID)
{
  const double c1 = params_->c1_;
  const double c2 = params_->c2_;

  dPmodI(0) += c1;
  dPmodI(1) += c2;
}