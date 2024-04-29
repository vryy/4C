/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of the base functionality of an exponential anisotropic summand

\level 1
*/
/*----------------------------------------------------------------------*/

#include "4C_matelast_coupanisoexpobase.hpp"

#include "4C_mat_par_material.hpp"

FOUR_C_NAMESPACE_OPEN


MAT::ELASTIC::PAR::CoupAnisoExpoBase::CoupAnisoExpoBase(
    const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : k1_(matdata->Get<double>("K1")),
      k2_(matdata->Get<double>("K2")),
      gamma_(matdata->Get<double>("GAMMA")),
      k1comp_(matdata->Get<double>("K1COMP")),
      k2comp_(matdata->Get<double>("K2COMP")),
      init_(matdata->Get<int>("INIT"))
{
}

MAT::ELASTIC::PAR::CoupAnisoExpoBase::CoupAnisoExpoBase()
    : k1_(0.0), k2_(0.0), gamma_(0.0), k1comp_(0.0), k2comp_(0.0), init_(0.0)
{
}

MAT::ELASTIC::CoupAnisoExpoBase::CoupAnisoExpoBase(MAT::ELASTIC::PAR::CoupAnisoExpoBase* params)
    : params_(params)
{
}

void MAT::ELASTIC::CoupAnisoExpoBase::AddStrainEnergy(double& psi,
    const CORE::LINALG::Matrix<3, 1>& prinv, const CORE::LINALG::Matrix<3, 1>& modinv,
    const CORE::LINALG::Matrix<6, 1>& glstrain, const int gp, const int eleGID)
{
  // right Cauchy Green
  CORE::LINALG::Matrix<3, 3> C(true);
  for (int i = 0; i < 3; ++i) C(i, i) = 2.0 * glstrain(i) + 1.0;
  C(0, 1) = C(1, 0) = glstrain(3);
  C(1, 2) = C(2, 1) = glstrain(4);
  C(0, 2) = C(2, 0) = glstrain(5);

  EvaluateFunc<double>(psi, C, gp, eleGID);
}

template <typename T>
void MAT::ELASTIC::CoupAnisoExpoBase::EvaluateFunc(
    T& psi, CORE::LINALG::Matrix<3, 3, T> const& C, const int gp, int const eleGID) const
{
  CORE::LINALG::Matrix<3, 3, T> A_T(GetCoupAnisoExpoBaseInterface().GetStructuralTensor(gp).A());
  const double scalarProduct = GetCoupAnisoExpoBaseInterface().GetScalarProduct(gp);

  T I4 = C.Dot(A_T);

  T k1;
  T k2;
  if (I4 < 1.0)
  {
    k1 = params_->k1comp_;
    k2 = params_->k2comp_;
  }
  else
  {
    k1 = params_->k1_;
    k2 = params_->k2_;
  }

  psi += (k1 / (2.0 * k2)) * (std::exp(k2 * (I4 - scalarProduct) * (I4 - scalarProduct)) - 1.0);
}

void MAT::ELASTIC::CoupAnisoExpoBase::EvaluateFirstDerivativesAniso(
    CORE::LINALG::Matrix<2, 1>& dPI_aniso, CORE::LINALG::Matrix<3, 3> const& rcg, int gp,
    int eleGID)
{
  double I4 = GetCoupAnisoExpoBaseInterface().GetStructuralTensor(gp).Dot(rcg);
  const double scalarProduct = GetCoupAnisoExpoBaseInterface().GetScalarProduct(gp);

  double k1 = params_->k1_;
  double k2 = params_->k2_;

  if (I4 < 1.0)
  {
    k1 = params_->k1comp_;
    k2 = params_->k2comp_;
  }

  dPI_aniso(0) =
      k1 * (I4 - scalarProduct) * std::exp(k2 * (I4 - scalarProduct) * (I4 - scalarProduct));
}

void MAT::ELASTIC::CoupAnisoExpoBase::EvaluateSecondDerivativesAniso(
    CORE::LINALG::Matrix<3, 1>& ddPII_aniso, CORE::LINALG::Matrix<3, 3> const& rcg, int gp,
    int eleGID)
{
  double I4 = GetCoupAnisoExpoBaseInterface().GetStructuralTensor(gp).Dot(rcg);
  const double scalarProduct = GetCoupAnisoExpoBaseInterface().GetScalarProduct(gp);

  double k1 = params_->k1_;
  double k2 = params_->k2_;

  if (I4 < 1.0)
  {
    k1 = params_->k1comp_;
    k2 = params_->k2comp_;
  }

  ddPII_aniso(0) = (1.0 + 2.0 * k2 * std::pow((I4 - scalarProduct), 2)) * k1 *
                   std::exp(k2 * std::pow((I4 - scalarProduct), 2));
}

template <typename T>
void MAT::ELASTIC::CoupAnisoExpoBase::GetDerivativesAniso(CORE::LINALG::Matrix<2, 1, T>& dPI_aniso,
    CORE::LINALG::Matrix<3, 1, T>& ddPII_aniso, CORE::LINALG::Matrix<4, 1, T>& dddPIII_aniso,
    CORE::LINALG::Matrix<3, 3, T> const& rcg, const int gp, const int eleGID) const
{
  CORE::LINALG::Matrix<3, 3, T> AM(GetCoupAnisoExpoBaseInterface().GetStructuralTensor(gp).A());
  const double scalarProduct = GetCoupAnisoExpoBaseInterface().GetScalarProduct(gp);

  T I4 = AM.Dot(rcg);

  T k1 = params_->k1_;
  T k2 = params_->k2_;

  if (I4 < 1.0)
  {
    k1 = params_->k1comp_;
    k2 = params_->k2comp_;
  }


  dPI_aniso(0) = k1 * (I4 - scalarProduct) * std::exp(k2 * std::pow((I4 - scalarProduct), 2));

  ddPII_aniso(0) = (1.0 + 2.0 * k2 * std::pow((I4 - scalarProduct), 2)) * k1 *
                   std::exp(k2 * std::pow((I4 - scalarProduct), 2));

  dddPIII_aniso(0) = (3.0 + 2.0 * k2 * std::pow((I4 - scalarProduct), 2)) * 2.0 * k1 * k2 *
                     (I4 - scalarProduct) * std::exp(k2 * std::pow((I4 - scalarProduct), 2));
}

void MAT::ELASTIC::CoupAnisoExpoBase::AddStressAnisoPrincipal(const CORE::LINALG::Matrix<6, 1>& rcg,
    CORE::LINALG::Matrix<6, 6>& cmat, CORE::LINALG::Matrix<6, 1>& stress,
    Teuchos::ParameterList& params, const int gp, const int eleGID)
{
  double I4 = GetCoupAnisoExpoBaseInterface().GetStructuralTensor_stress(gp).Dot(rcg);
  const double scalarProduct = GetCoupAnisoExpoBaseInterface().GetScalarProduct(gp);

  double k1;
  double k2;
  if (I4 < 1.0)
  {
    k1 = params_->k1comp_;
    k2 = params_->k2comp_;
  }
  else
  {
    k1 = params_->k1_;
    k2 = params_->k2_;
  }

  double gamma =
      2. * (k1 * (I4 - scalarProduct) * std::exp(k2 * std::pow((I4 - scalarProduct), 2)));
  stress.Update(gamma, GetCoupAnisoExpoBaseInterface().GetStructuralTensor_stress(gp), 1.0);

  double delta = 2. * (1. + 2. * k2 * std::pow((I4 - scalarProduct), 2)) * 2. * k1 *
                 std::exp(k2 * std::pow((I4 - scalarProduct), 2));
  cmat.MultiplyNT(delta, GetCoupAnisoExpoBaseInterface().GetStructuralTensor_stress(gp),
      GetCoupAnisoExpoBaseInterface().GetStructuralTensor_stress(gp), 1.0);
}

void MAT::ELASTIC::CoupAnisoExpoBase::GetFiberVecs(
    std::vector<CORE::LINALG::Matrix<3, 1>>& fibervecs)
{
  FOUR_C_THROW("Getting the fiber vectors is not implemented in the base version of CoupAnisoExpo");
}

void MAT::ELASTIC::CoupAnisoExpoBase::SetFiberVecs(const double newgamma,
    const CORE::LINALG::Matrix<3, 3>& locsys, const CORE::LINALG::Matrix<3, 3>& defgrd)
{
  FOUR_C_THROW("Setting the fiber vectors is not implemented in the base version of CoupAnisoExpo");
}

void MAT::ELASTIC::CoupAnisoExpoBase::SetFiberVecs(const CORE::LINALG::Matrix<3, 1>& fibervec)
{
  FOUR_C_THROW("Setting the fiber vectors is not implemented in the base version of CoupAnisoExpo");
}


// explicit instantiation of template functions
template void MAT::ELASTIC::CoupAnisoExpoBase::GetDerivativesAniso<double>(
    CORE::LINALG::Matrix<2, 1, double>&, CORE::LINALG::Matrix<3, 1, double>&,
    CORE::LINALG::Matrix<4, 1, double>&, CORE::LINALG::Matrix<3, 3, double> const&, int, int) const;
template void MAT::ELASTIC::CoupAnisoExpoBase::GetDerivativesAniso<FAD>(
    CORE::LINALG::Matrix<2, 1, FAD>&, CORE::LINALG::Matrix<3, 1, FAD>&,
    CORE::LINALG::Matrix<4, 1, FAD>&, CORE::LINALG::Matrix<3, 3, FAD> const&, int, int) const;
template void MAT::ELASTIC::CoupAnisoExpoBase::EvaluateFunc<double>(
    double&, CORE::LINALG::Matrix<3, 3, double> const&, int, int) const;
template void MAT::ELASTIC::CoupAnisoExpoBase::EvaluateFunc<FAD>(
    FAD&, CORE::LINALG::Matrix<3, 3, FAD> const&, int, int) const;

FOUR_C_NAMESPACE_CLOSE
