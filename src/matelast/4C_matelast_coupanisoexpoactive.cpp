/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of the coupled contribution of an anisotropic active fiber material

\level 2
*/
/*----------------------------------------------------------------------*/

#include "4C_matelast_coupanisoexpoactive.hpp"

#include "4C_io_linedefinition.hpp"
#include "4C_matelast_aniso_structuraltensor_strategy.hpp"
#include "4C_material_input_base.hpp"

FOUR_C_NAMESPACE_OPEN


MAT::ELASTIC::PAR::CoupAnisoExpoActive::CoupAnisoExpoActive(
    const Teuchos::RCP<CORE::MAT::PAR::Material>& matdata)
    : ParameterAniso(matdata),
      k1_(matdata->Get<double>("K1")),
      k2_(matdata->Get<double>("K2")),
      gamma_(matdata->Get<double>("GAMMA")),
      k1comp_(matdata->Get<double>("K1COMP")),
      k2comp_(matdata->Get<double>("K2COMP")),
      init_(matdata->Get<int>("INIT")),
      adapt_angle_(matdata->Get<bool>("ADAPT_ANGLE")),
      s_(matdata->Get<double>("S")),
      lambdamax_(matdata->Get<double>("LAMBDAMAX")),
      lambda0_(matdata->Get<double>("LAMBDA0")),
      dens_(matdata->Get<double>("DENS"))
{
}

MAT::ELASTIC::CoupAnisoExpoActive::CoupAnisoExpoActive(
    MAT::ELASTIC::PAR::CoupAnisoExpoActive* params)
    : params_(params),
      anisotropy_extension_(params_->init_, params->gamma_, params_->adapt_angle_ != 0,
          params_->structural_tensor_strategy(), {0})
{
  d_p_iact_ = 0.0;
  lambdaact_ = 1.0;
  anisotropy_extension_.register_needed_tensors(
      FiberAnisotropyExtension<1>::FIBER_VECTORS |
      FiberAnisotropyExtension<1>::STRUCTURAL_TENSOR_STRESS |
      FiberAnisotropyExtension<1>::STRUCTURAL_TENSOR);
}

void MAT::ELASTIC::CoupAnisoExpoActive::register_anisotropy_extensions(MAT::Anisotropy& anisotropy)
{
  anisotropy.register_anisotropy_extension(anisotropy_extension_);
}

void MAT::ELASTIC::CoupAnisoExpoActive::PackSummand(CORE::COMM::PackBuffer& data) const
{
  AddtoPack(data, lambdaact_);
  anisotropy_extension_.PackAnisotropy(data);
}

void MAT::ELASTIC::CoupAnisoExpoActive::UnpackSummand(
    const std::vector<char>& data, std::vector<char>::size_type& position)
{
  ExtractfromPack(position, data, lambdaact_);
  anisotropy_extension_.UnpackAnisotropy(data, position);

  d_p_iact_ = EvaluatedPsiActive();
}

void MAT::ELASTIC::CoupAnisoExpoActive::Setup(int numgp, INPUT::LineDefinition* linedef)
{
  // setup first derivative of active fiber potential w.r.t. active fiber stretch (const during the
  // whole simulation)
  lambdaact_ = 1.0;

  d_p_iact_ = EvaluatedPsiActive();
}

void MAT::ELASTIC::CoupAnisoExpoActive::AddStrainEnergy(double& psi,
    const CORE::LINALG::Matrix<3, 1>& prinv, const CORE::LINALG::Matrix<3, 1>& modinv,
    const CORE::LINALG::Matrix<6, 1>& glstrain, const int gp, const int eleGID)
{
  // rigth Cauchy Green in strain like Voigt notation
  CORE::LINALG::Matrix<6, 1> rcg(true);

  for (int i = 0; i < 3; ++i) rcg(i) = 2.0 * glstrain(i) + 1.0;
  rcg(3) = 2.0 * glstrain(3);
  rcg(4) = 2.0 * glstrain(4);
  rcg(5) = 2.0 * glstrain(5);

  double I4 = anisotropy_extension_.get_structural_tensor_stress(gp, 0).Dot(glstrain);

  double k1 = params_->k1_;
  double k2 = params_->k2_;

  if (prinv(0) < 1.0)
  {
    k1 = params_->k1comp_;
    k2 = params_->k2comp_;
  }

  // passive contribution
  psi += (k1 / (2.0 * k2)) * (exp(k2 * (I4 - 1.0) * (I4 - 1.0)) - 1.0);

  // active contribution
  psi += params_->s_ / params_->dens_ *
         (lambdaact_ + (1. / 3.) * (std::pow(params_->lambdamax_ - lambdaact_, 3.0) /
                                       std::pow(params_->lambdamax_ - params_->lambda0_, 2.0)));
}

template <typename T>
void MAT::ELASTIC::CoupAnisoExpoActive::EvaluateFunc(
    T& psi, CORE::LINALG::Matrix<3, 3, T> const& rcg, const int gp, int const eleGID) const
{
  T I4_fad;
  static CORE::LINALG::Matrix<6, 1, T> Av_T(true);
  for (int i = 0; i < 6; ++i)
    Av_T(i) = anisotropy_extension_.get_structural_tensor_stress(gp, 0)(i);
  I4_fad = Av_T(0) * rcg(0, 0) + Av_T(1) * rcg(1, 1) + Av_T(2) * rcg(2, 2) +
           Av_T(3) * (rcg(0, 1) + rcg(1, 0)) + Av_T(4) * (rcg(1, 2) + rcg(2, 1)) +
           Av_T(5) * (rcg(0, 2) + rcg(2, 0));

  T k1 = params_->k1_;
  T k2 = params_->k2_;

  if (I4_fad < 1.0)
  {
    k1 = params_->k1comp_;
    k2 = params_->k2comp_;
  }

  psi += (k1 / (2.0 * k2)) * (exp(k2 * (I4_fad - 1.0) * (I4_fad - 1.0)) - 1.0);
};

template <typename T>
void MAT::ELASTIC::CoupAnisoExpoActive::evaluate_active_stress_cmat_aniso(
    CORE::LINALG::Matrix<3, 3, T> const& CM, CORE::LINALG::Matrix<6, 6, T>& cmat,
    CORE::LINALG::Matrix<6, 1, T>& stress, const int gp, const int eleGID) const
{
  T lambda_sq = 0.0;
  static CORE::LINALG::Matrix<6, 1, T> Av_T(true);
  for (int i = 0; i < 6; ++i)
    Av_T(i) = anisotropy_extension_.get_structural_tensor_stress(gp, 0)(i);
  lambda_sq = Av_T(0) * CM(0, 0) + Av_T(1) * CM(1, 1) + Av_T(2) * CM(2, 2) + Av_T(3) * CM(0, 1) +
              Av_T(4) * CM(1, 2) + Av_T(5) * CM(0, 2) + Av_T(3) * CM(1, 0) + Av_T(4) * CM(2, 1) +
              Av_T(5) * CM(2, 0);

  T dPIact_T = 0.0;
  dPIact_T = d_p_iact_;
  stress.Update(dPIact_T * 1. / lambda_sq, Av_T, 0.0);
  cmat.MultiplyNT(-2.0 * dPIact_T * 1. / (lambda_sq * lambda_sq), Av_T, Av_T, 0.0);
}

void MAT::ELASTIC::CoupAnisoExpoActive::add_active_stress_cmat_aniso(
    CORE::LINALG::Matrix<3, 3> const& CM, CORE::LINALG::Matrix<6, 6>& cmat,
    CORE::LINALG::Matrix<6, 1>& stress, const int gp, const int eleGID) const
{
  double lambda_sq = CM.Dot(anisotropy_extension_.GetStructuralTensor(gp, 0));

  double dPIact_T = 0.0;
  dPIact_T = d_p_iact_;
  stress.Update(
      dPIact_T * 1. / lambda_sq, anisotropy_extension_.get_structural_tensor_stress(gp, 0), 1.0);
  cmat.MultiplyNT(-2.0 * dPIact_T * 1. / (lambda_sq * lambda_sq),
      anisotropy_extension_.get_structural_tensor_stress(gp, 0),
      anisotropy_extension_.get_structural_tensor_stress(gp, 0), 1.0);
}

void MAT::ELASTIC::CoupAnisoExpoActive::evaluate_first_derivatives_aniso(
    CORE::LINALG::Matrix<2, 1>& dPI_aniso, CORE::LINALG::Matrix<3, 3> const& rcg, int gp,
    int eleGID)
{
  double I4 = anisotropy_extension_.GetStructuralTensor(gp, 0).Dot(rcg);

  double k1 = params_->k1_;
  double k2 = params_->k2_;

  if (I4 < 1.0)
  {
    k1 = params_->k1comp_;
    k2 = params_->k2comp_;
  }

  dPI_aniso(0) = k1 * (I4 - 1.0) * exp(k2 * (I4 - 1.0) * (I4 - 1.0));
}

void MAT::ELASTIC::CoupAnisoExpoActive::evaluate_second_derivatives_aniso(
    CORE::LINALG::Matrix<3, 1>& ddPII_aniso, CORE::LINALG::Matrix<3, 3> const& rcg, int gp,
    int eleGID)
{
  double I4 = anisotropy_extension_.GetStructuralTensor(gp, 0).Dot(rcg);

  double k1 = params_->k1_;
  double k2 = params_->k2_;

  if (I4 < 1.0)
  {
    k1 = params_->k1comp_;
    k2 = params_->k2comp_;
  }

  ddPII_aniso(0) =
      (1.0 + 2.0 * k2 * (I4 - 1.0) * (I4 - 1.0)) * k1 * exp(k2 * (I4 - 1.0) * (I4 - 1.0));
}

template <typename T>
void MAT::ELASTIC::CoupAnisoExpoActive::GetDerivativesAniso(
    CORE::LINALG::Matrix<2, 1, T>& dPI_aniso, CORE::LINALG::Matrix<3, 1, T>& ddPII_aniso,
    CORE::LINALG::Matrix<4, 1, T>& dddPIII_aniso, CORE::LINALG::Matrix<3, 3, T> const& rcg,
    const int gp, const int eleGID) const
{
  T I4 = 0.0;
  CORE::LINALG::Matrix<3, 3, T> AM(true);
  for (int i = 0; i < 3; ++i)
    AM(i, i) = anisotropy_extension_.get_structural_tensor_stress(gp, 0)(i);
  AM(0, 1) = AM(1, 0) = anisotropy_extension_.get_structural_tensor_stress(gp, 0)(3);
  AM(2, 1) = AM(1, 2) = anisotropy_extension_.get_structural_tensor_stress(gp, 0)(4);
  AM(0, 2) = AM(2, 0) = anisotropy_extension_.get_structural_tensor_stress(gp, 0)(5);

  I4 = AM.Dot(rcg);

  T k1 = params_->k1_;
  T k2 = params_->k2_;

  if (I4 < 1.0)
  {
    k1 = params_->k1comp_;
    k2 = params_->k2comp_;
  }


  dPI_aniso(0) = k1 * (I4 - 1.0) * exp(k2 * (I4 - 1.0) * (I4 - 1.0));

  ddPII_aniso(0) =
      (1.0 + 2.0 * k2 * (I4 - 1.0) * (I4 - 1.0)) * k1 * exp(k2 * (I4 - 1.0) * (I4 - 1.0));

  dddPIII_aniso(0) = (3.0 + 2.0 * k2 * (I4 - 1.0) * (I4 - 1.0)) * 2.0 * k1 * k2 * (I4 - 1.0) *
                     exp(k2 * (I4 - 1.0) * (I4 - 1.0));
};

void MAT::ELASTIC::CoupAnisoExpoActive::add_stress_aniso_principal(
    const CORE::LINALG::Matrix<6, 1>& rcg, CORE::LINALG::Matrix<6, 6>& cmat,
    CORE::LINALG::Matrix<6, 1>& stress, Teuchos::ParameterList& params, const int gp,
    const int eleGID)
{
  double I4 = anisotropy_extension_.get_structural_tensor_stress(gp, 0).Dot(rcg);

  double k1 = params_->k1_;
  double k2 = params_->k2_;

  if (I4 < 1.0)
  {
    k1 = params_->k1comp_;
    k2 = params_->k2comp_;
  }

  double gamma = 2. * (k1 * (I4 - 1.) * exp(k2 * (I4 - 1.) * (I4 - 1.)));
  stress.Update(gamma, anisotropy_extension_.get_structural_tensor_stress(gp, 0), 1.0);

  double delta =
      2. * (1. + 2. * k2 * (I4 - 1.) * (I4 - 1.)) * 2. * k1 * exp(k2 * (I4 - 1.) * (I4 - 1.));
  cmat.MultiplyNT(delta, anisotropy_extension_.get_structural_tensor_stress(gp, 0),
      anisotropy_extension_.get_structural_tensor_stress(gp, 0), 1.0);
}

void MAT::ELASTIC::CoupAnisoExpoActive::GetFiberVecs(
    std::vector<CORE::LINALG::Matrix<3, 1>>& fibervecs)
{
  // this method does not support Gauss point fibers
  fibervecs.push_back(anisotropy_extension_.GetFiber(BaseAnisotropyExtension::GPDEFAULT, 0));
}

void MAT::ELASTIC::CoupAnisoExpoActive::SetFiberVecs(const double newgamma,
    const CORE::LINALG::Matrix<3, 3>& locsys, const CORE::LINALG::Matrix<3, 3>& defgrd)
{
  anisotropy_extension_.SetFiberVecs(newgamma, locsys, defgrd);
}

double MAT::ELASTIC::CoupAnisoExpoActive::EvaluatedPsiActive() const
{
  return params_->s_ / params_->dens_ *
         (1.0 - std::pow(params_->lambdamax_ - lambdaact_, 2.0) /
                    std::pow(params_->lambdamax_ - params_->lambda0_, 2.0));
}


// explicit instantiation of template functions
template void MAT::ELASTIC::CoupAnisoExpoActive::GetDerivativesAniso<double>(
    CORE::LINALG::Matrix<2, 1, double>&, CORE::LINALG::Matrix<3, 1, double>&,
    CORE::LINALG::Matrix<4, 1, double>&, CORE::LINALG::Matrix<3, 3, double> const&, int,
    const int) const;
template void MAT::ELASTIC::CoupAnisoExpoActive::GetDerivativesAniso<FAD>(
    CORE::LINALG::Matrix<2, 1, FAD>&, CORE::LINALG::Matrix<3, 1, FAD>&,
    CORE::LINALG::Matrix<4, 1, FAD>&, CORE::LINALG::Matrix<3, 3, FAD> const&, int, const int) const;
template void MAT::ELASTIC::CoupAnisoExpoActive::evaluate_active_stress_cmat_aniso<double>(
    CORE::LINALG::Matrix<3, 3, double> const&, CORE::LINALG::Matrix<6, 6, double>&,
    CORE::LINALG::Matrix<6, 1, double>&, int, const int) const;
template void MAT::ELASTIC::CoupAnisoExpoActive::evaluate_active_stress_cmat_aniso<FAD>(
    CORE::LINALG::Matrix<3, 3, FAD> const&, CORE::LINALG::Matrix<6, 6, FAD>&,
    CORE::LINALG::Matrix<6, 1, FAD>&, int, const int) const;
template void MAT::ELASTIC::CoupAnisoExpoActive::EvaluateFunc<double>(
    double&, CORE::LINALG::Matrix<3, 3, double> const&, int, const int) const;
template void MAT::ELASTIC::CoupAnisoExpoActive::EvaluateFunc<FAD>(
    FAD&, CORE::LINALG::Matrix<3, 3, FAD> const&, int, const int) const;

FOUR_C_NAMESPACE_CLOSE
