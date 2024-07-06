/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of the volumetic penalty material as in Roernbauer2008 (student
thesis)

\level 1
*/
/*----------------------------------------------------------------------*/

#include "4C_matelast_volpenalty.hpp"

#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN


Mat::Elastic::PAR::VolPenalty::VolPenalty(const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      eps_(matdata.parameters.get<double>("EPSILON")),
      gam_(matdata.parameters.get<double>("GAMMA"))
{
  if (eps_ < 0. || gam_ <= 0.)
    FOUR_C_THROW("VolPenalty parameters EPSILON and GAMMA have to be greater zero");
}

Mat::Elastic::VolPenalty::VolPenalty(Mat::Elastic::PAR::VolPenalty* params) : params_(params) {}

void Mat::Elastic::VolPenalty::add_strain_energy(double& psi,
    const Core::LinAlg::Matrix<3, 1>& prinv, const Core::LinAlg::Matrix<3, 1>& modinv,
    const Core::LinAlg::Matrix<6, 1>& glstrain, const int gp, const int eleGID)
{
  const double eps = params_->eps_;
  const double gam = params_->gam_;

  // strain energy: Psi=\epsilon \left( J^{\gamma} + \frac 1 {J^{\gamma}} -2 \right)
  // add to overall strain energy
  psi += eps * (pow(modinv(2), gam) + pow(modinv(2), -gam) - 2.);
}

void Mat::Elastic::VolPenalty::add_derivatives_modified(Core::LinAlg::Matrix<3, 1>& dPmodI,
    Core::LinAlg::Matrix<6, 1>& ddPmodII, const Core::LinAlg::Matrix<3, 1>& modinv, const int gp,
    const int eleGID)
{
  const double eps = params_->eps_;
  const double gam = params_->gam_;

  dPmodI(2) += eps * gam * (pow(modinv(2), gam - 1.) - pow(modinv(2), -gam - 1.));

  ddPmodII(2) +=
      eps * gam * ((gam - 1.) * pow(modinv(2), gam - 2.) + (gam + 1.) * pow(modinv(2), -gam - 2.));
}

void Mat::Elastic::VolPenalty::add3rd_vol_deriv(
    const Core::LinAlg::Matrix<3, 1>& modinv, double& d3PsiVolDJ3)
{
  const double eps = params_->eps_;
  const double gam = params_->gam_;
  const double J = modinv(2);
  d3PsiVolDJ3 += eps * (gam * (gam - 1.) * (gam - 2.) * pow(J, gam - 3.) +
                           (-gam) * (-gam - 1.) * (-gam - 2.) * pow(J, -gam - 3.));
}
FOUR_C_NAMESPACE_CLOSE
