/*----------------------------------------------------------------------*/
/*! \file
\brief the input line should read
     MAT 1 ELAST_CoupAnisoExpoActive K1 10.0 K2 1.0 GAMMA 35.0 K1COMP 0.0 K2COMP 1.0 INIT 0
ADAPT_ANGLE 0 S 54000 LAMBDAMAX 1.4 LAMBDA0 0.8 DENS 1050

\level 2

\maintainer Amadeus Gebauer
*/

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_coupanisoexpoactive.H"
#include "elast_aniso_structuraltensor_strategy.H"

#include "../drt_mat/matpar_material.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_linedefinition.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::CoupAnisoExpoActive::CoupAnisoExpoActive(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : ParameterAniso(matdata),
      k1_(matdata->GetDouble("K1")),
      k2_(matdata->GetDouble("K2")),
      gamma_(matdata->GetDouble("GAMMA")),
      k1comp_(matdata->GetDouble("K1COMP")),
      k2comp_(matdata->GetDouble("K2COMP")),
      init_(matdata->GetInt("INIT")),
      adapt_angle_(matdata->GetInt("ADAPT_ANGLE")),
      s_(matdata->GetDouble("S")),
      lambdamax_(matdata->GetDouble("LAMBDAMAX")),
      lambda0_(matdata->GetDouble("LAMBDA0")),
      dens_(matdata->GetDouble("DENS"))
{
}


/*----------------------------------------------------------------------*
 |  Constructor                             (public)   fb         07/16 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::CoupAnisoExpoActive::CoupAnisoExpoActive(
    MAT::ELASTIC::PAR::CoupAnisoExpoActive* params)
    : params_(params),
      anisotropyExtension_(params_->init_, params->gamma_, params_->adapt_angle_ != 0,
          params_->StructuralTensorStrategy())
{
  dPIact_ = 0.0;
  lambdaact_ = 1.0;
  anisotropyExtension_.RegisterNeededTensors(FiberAnisotropyExtension::FIBER_VECTORS |
                                             FiberAnisotropyExtension::STRUCTURAL_TENSOR_STRESS |
                                             FiberAnisotropyExtension::STRUCTURAL_TENSOR);
}

void MAT::ELASTIC::CoupAnisoExpoActive::RegisterAnisotropyExtensions(MAT::Anisotropy& anisotropy)
{
  anisotropy.RegisterAnisotropyExtension(anisotropyExtension_);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoExpoActive::PackSummand(DRT::PackBuffer& data) const
{
  AddtoPack(data, dPIact_);
  AddtoPack(data, lambdaact_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoExpoActive::UnpackSummand(
    const std::vector<char>& data, std::vector<char>::size_type& position)
{
  ExtractfromPack(position, data, dPIact_);
  ExtractfromPack(position, data, lambdaact_);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoExpoActive::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  // setup first derivative of active fiber potential w.r.t. active fiber stretch (const during the
  // whole simulation)
  lambdaact_ = 1.0;

  dPIact_ = params_->s_ / params_->dens_ *
            (1.0 - std::pow(params_->lambdamax_ - lambdaact_, 2.0) /
                       std::pow(params_->lambdamax_ - params_->lambda0_, 2.0));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoExpoActive::AddStrainEnergy(double& psi,
    const LINALG::Matrix<3, 1>& prinv, const LINALG::Matrix<3, 1>& modinv,
    const LINALG::Matrix<6, 1>& glstrain, const int eleGID)
{
  // rigth Cauchy Green in strain like Voigt notation
  LINALG::Matrix<6, 1> rcg(true);

  for (int i = 0; i < 3; ++i) rcg(i) = 2.0 * glstrain(i) + 1.0;
  rcg(3) = 2.0 * glstrain(3);
  rcg(4) = 2.0 * glstrain(4);
  rcg(5) = 2.0 * glstrain(5);

  double I4 = anisotropyExtension_.GetStructuralTensor_stress(0, 0).Dot(glstrain);

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


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <typename T>
void MAT::ELASTIC::CoupAnisoExpoActive::EvaluateFunc(
    T& psi, LINALG::Matrix<3, 3, T> const& rcg, int const eleGID) const
{
  T I4_fad;
  static LINALG::Matrix<6, 1, T> Av_T(true);
  for (int i = 0; i < 6; ++i) Av_T(i) = anisotropyExtension_.GetStructuralTensor_stress(0, 0)(i);
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


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <typename T>
void MAT::ELASTIC::CoupAnisoExpoActive::EvaluateActiveStressCmatAniso(
    LINALG::Matrix<3, 3, T> const& CM, LINALG::Matrix<6, 6, T>& cmat,
    LINALG::Matrix<6, 1, T>& stress, const int eleGID) const
{
  T lambda_sq = 0.0;
  static LINALG::Matrix<6, 1, T> Av_T(true);
  for (int i = 0; i < 6; ++i) Av_T(i) = anisotropyExtension_.GetStructuralTensor_stress(0, 0)(i);
  lambda_sq = Av_T(0) * CM(0, 0) + Av_T(1) * CM(1, 1) + Av_T(2) * CM(2, 2) + Av_T(3) * CM(0, 1) +
              Av_T(4) * CM(1, 2) + Av_T(5) * CM(0, 2) + Av_T(3) * CM(1, 0) + Av_T(4) * CM(2, 1) +
              Av_T(5) * CM(2, 0);

  T dPIact_T = 0.0;
  dPIact_T = dPIact_;
  stress.Update(dPIact_T * 1. / lambda_sq, Av_T, 0.0);
  cmat.MultiplyNT(-2.0 * dPIact_T * 1. / (lambda_sq * lambda_sq), Av_T, Av_T, 0.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoExpoActive::AddActiveStressCmatAniso(LINALG::Matrix<3, 3> const& CM,
    LINALG::Matrix<6, 6>& cmat, LINALG::Matrix<6, 1>& stress, const int eleGID) const
{
  double lambda_sq = CM.Dot(anisotropyExtension_.GetStructuralTensor(0, 0));

  double dPIact_T = 0.0;
  dPIact_T = dPIact_;
  stress.Update(
      dPIact_T * 1. / lambda_sq, anisotropyExtension_.GetStructuralTensor_stress(0, 0), 1.0);
  cmat.MultiplyNT(-2.0 * dPIact_T * 1. / (lambda_sq * lambda_sq),
      anisotropyExtension_.GetStructuralTensor_stress(0, 0),
      anisotropyExtension_.GetStructuralTensor_stress(0, 0), 1.0);
}

void MAT::ELASTIC::CoupAnisoExpoActive::EvaluateFirstDerivativesAniso(
    LINALG::Matrix<2, 1>& dPI_aniso, LINALG::Matrix<3, 3> const& rcg, int gp, int eleGID)
{
  double I4 = anisotropyExtension_.GetStructuralTensor(gp, 0).Dot(rcg);

  double k1 = params_->k1_;
  double k2 = params_->k2_;

  if (I4 < 1.0)
  {
    k1 = params_->k1comp_;
    k2 = params_->k2comp_;
  }

  dPI_aniso(0) = k1 * (I4 - 1.0) * exp(k2 * (I4 - 1.0) * (I4 - 1.0));
}

void MAT::ELASTIC::CoupAnisoExpoActive::EvaluateSecondDerivativesAniso(
    LINALG::Matrix<3, 1>& ddPII_aniso, LINALG::Matrix<3, 3> const& rcg, int gp, int eleGID)
{
  double I4 = anisotropyExtension_.GetStructuralTensor(gp, 0).Dot(rcg);

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


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <typename T>
void MAT::ELASTIC::CoupAnisoExpoActive::GetDerivativesAniso(LINALG::Matrix<2, 1, T>& dPI_aniso,
    LINALG::Matrix<3, 1, T>& ddPII_aniso, LINALG::Matrix<4, 1, T>& dddPIII_aniso,
    LINALG::Matrix<3, 3, T> const& rcg, const int eleGID) const
{
  T I4 = 0.0;
  LINALG::Matrix<3, 3, T> AM(true);
  for (int i = 0; i < 3; ++i) AM(i, i) = anisotropyExtension_.GetStructuralTensor_stress(0, 0)(i);
  AM(0, 1) = AM(1, 0) = anisotropyExtension_.GetStructuralTensor_stress(0, 0)(3);
  AM(2, 1) = AM(1, 2) = anisotropyExtension_.GetStructuralTensor_stress(0, 0)(4);
  AM(0, 2) = AM(2, 0) = anisotropyExtension_.GetStructuralTensor_stress(0, 0)(5);

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

  return;
};


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoExpoActive::AddStressAnisoPrincipal(const LINALG::Matrix<6, 1>& rcg,
    LINALG::Matrix<6, 6>& cmat, LINALG::Matrix<6, 1>& stress, Teuchos::ParameterList& params,
    const int eleGID)
{
  double I4 = anisotropyExtension_.GetStructuralTensor_stress(0, 0).Dot(rcg);

  double k1 = params_->k1_;
  double k2 = params_->k2_;

  if (I4 < 1.0)
  {
    k1 = params_->k1comp_;
    k2 = params_->k2comp_;
  }

  double gamma = 2. * (k1 * (I4 - 1.) * exp(k2 * (I4 - 1.) * (I4 - 1.)));
  stress.Update(gamma, anisotropyExtension_.GetStructuralTensor_stress(0, 0), 1.0);

  double delta =
      2. * (1. + 2. * k2 * (I4 - 1.) * (I4 - 1.)) * 2. * k1 * exp(k2 * (I4 - 1.) * (I4 - 1.));
  cmat.MultiplyNT(delta, anisotropyExtension_.GetStructuralTensor_stress(0, 0),
      anisotropyExtension_.GetStructuralTensor_stress(0, 0), 1.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoExpoActive::GetFiberVecs(std::vector<LINALG::Matrix<3, 1>>& fibervecs)
{
  fibervecs.push_back(anisotropyExtension_.GetFiber(BaseAnisotropyExtension::GPDEFAULT, 0));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoExpoActive::SetFiberVecs(
    const double newgamma, const LINALG::Matrix<3, 3>& locsys, const LINALG::Matrix<3, 3>& defgrd)
{
  anisotropyExtension_.SetFiberVecs(newgamma, locsys, defgrd);
}


// explicit instantiation of template functions
template void MAT::ELASTIC::CoupAnisoExpoActive::GetDerivativesAniso<double>(
    LINALG::Matrix<2, 1, double>&, LINALG::Matrix<3, 1, double>&, LINALG::Matrix<4, 1, double>&,
    LINALG::Matrix<3, 3, double> const&, const int) const;
template void MAT::ELASTIC::CoupAnisoExpoActive::GetDerivativesAniso<FAD>(
    LINALG::Matrix<2, 1, FAD>&, LINALG::Matrix<3, 1, FAD>&, LINALG::Matrix<4, 1, FAD>&,
    LINALG::Matrix<3, 3, FAD> const&, const int) const;
template void MAT::ELASTIC::CoupAnisoExpoActive::EvaluateActiveStressCmatAniso<double>(
    LINALG::Matrix<3, 3, double> const&, LINALG::Matrix<6, 6, double>&,
    LINALG::Matrix<6, 1, double>&, const int) const;
template void MAT::ELASTIC::CoupAnisoExpoActive::EvaluateActiveStressCmatAniso<FAD>(
    LINALG::Matrix<3, 3, FAD> const&, LINALG::Matrix<6, 6, FAD>&, LINALG::Matrix<6, 1, FAD>&,
    const int) const;
template void MAT::ELASTIC::CoupAnisoExpoActive::EvaluateFunc<double>(
    double&, LINALG::Matrix<3, 3, double> const&, const int) const;
template void MAT::ELASTIC::CoupAnisoExpoActive::EvaluateFunc<FAD>(
    FAD&, LINALG::Matrix<3, 3, FAD> const&, const int) const;
