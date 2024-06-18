/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of the isochoric contribution of a material to test the isochoric parts of the
Elasthyper-Toolbox.

\level 1
*/
/*----------------------------------------------------------------------*/

#include "4C_matelast_isotestmaterial.hpp"

#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN

Mat::Elastic::PAR::IsoTestMaterial::IsoTestMaterial(const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      c1_(matdata.parameters.get<double>("C1")),
      c2_(matdata.parameters.get<double>("C2"))
{
}

Mat::Elastic::IsoTestMaterial::IsoTestMaterial(Mat::Elastic::PAR::IsoTestMaterial* params)
    : params_(params)
{
}

void Mat::Elastic::IsoTestMaterial::AddStrainEnergy(double& psi,
    const Core::LinAlg::Matrix<3, 1>& prinv, const Core::LinAlg::Matrix<3, 1>& modinv,
    const Core::LinAlg::Matrix<6, 1>& glstrain, const int gp, const int eleGID)
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

void Mat::Elastic::IsoTestMaterial::add_derivatives_modified(Core::LinAlg::Matrix<3, 1>& dPmodI,
    Core::LinAlg::Matrix<6, 1>& ddPmodII, const Core::LinAlg::Matrix<3, 1>& modinv, const int gp,
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
FOUR_C_NAMESPACE_CLOSE
