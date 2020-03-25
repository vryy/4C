/*----------------------------------------------------------------------*/
/*! \file
\brief the input line should read
  MAT 1 ELAST_CoupAnisoExpo K1 10.0 K2 1.0 GAMMA 35.0  K1COMP 0.0 K2COMP 1.0 [INIT 1] [ADAPT_ANGLE
No]

\level 1

\maintainer Amadeus Gebauer
*/

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_coupanisoexpo.H"
#include "elast_aniso_structuraltensor_strategy.H"

#include "../drt_mat/matpar_material.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_linedefinition.H"


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::CoupAnisoExpo::CoupAnisoExpo(Teuchos::RCP<MAT::PAR::Material> matdata)
    : ParameterAniso(matdata),
      k1_(matdata->GetDouble("K1")),
      k2_(matdata->GetDouble("K2")),
      gamma_(matdata->GetDouble("GAMMA")),
      k1comp_(matdata->GetDouble("K1COMP")),
      k2comp_(matdata->GetDouble("K2COMP")),
      init_(matdata->GetInt("INIT")),
      adapt_angle_(matdata->GetInt("ADAPT_ANGLE"))
{
}


/*----------------------------------------------------------------------*
 |  Constructor                             (public)   st         03/12 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::CoupAnisoExpo::CoupAnisoExpo(MAT::ELASTIC::PAR::CoupAnisoExpo* params)
    : params_(params),
      anisotropyExtension_(params_->init_, params->gamma_, params_->adapt_angle_ != 0,
          params_->StructuralTensorStrategy())
{
  anisotropyExtension_.RegisterNeededTensors(FiberAnisotropyExtension::FIBER_VECTORS |
                                             FiberAnisotropyExtension::STRUCTURAL_TENSOR_STRESS |
                                             FiberAnisotropyExtension::STRUCTURAL_TENSOR);
}

void MAT::ELASTIC::CoupAnisoExpo::RegisterAnisotropyExtensions(MAT::Anisotropy& anisotropy)
{
  anisotropy.RegisterAnisotropyExtension(anisotropyExtension_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoExpo::AddStrainEnergy(double& psi, const LINALG::Matrix<3, 1>& prinv,
    const LINALG::Matrix<3, 1>& modinv, const LINALG::Matrix<6, 1>& glstrain, const int eleGID)
{
  // rigth Cauchy Green in strain-like Voigt notation
  LINALG::Matrix<6, 1> rcg(true);

  for (int i = 0; i < 3; ++i) rcg(i) = 2.0 * glstrain(i) + 1.0;
  rcg(3) = 2.0 * glstrain(3);
  rcg(4) = 2.0 * glstrain(4);
  rcg(5) = 2.0 * glstrain(5);

  double I4 = anisotropyExtension_.GetStructuralTensor_stress(BaseAnisotropyExtension::GPDEFAULT, 0)
                  .Dot(rcg);

  double k1 = params_->k1_;
  double k2 = params_->k2_;

  if (prinv(0) < 1.0)
  {
    k1 = params_->k1comp_;
    k2 = params_->k2comp_;
  }

  psi += (k1 / (2.0 * k2)) * (exp(k2 * (I4 - 1.0) * (I4 - 1.0)) - 1.0);
};


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <typename T>
void MAT::ELASTIC::CoupAnisoExpo::EvaluateFunc(
    T& psi, LINALG::Matrix<3, 3, T> const& rcg, int const eleGID) const
{
  T I4_fad = 0.0;
  static LINALG::Matrix<6, 1, T> Av_T(true);
  for (int i = 0; i < 6; ++i)
    Av_T(i) =
        anisotropyExtension_.GetStructuralTensor_stress(BaseAnisotropyExtension::GPDEFAULT, 0)(i);
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

  return;
};

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <typename T>
void MAT::ELASTIC::CoupAnisoExpo::GetDerivativesAniso(LINALG::Matrix<2, 1, T>& dPI_aniso,
    LINALG::Matrix<3, 1, T>& ddPII_aniso, LINALG::Matrix<4, 1, T>& dddPIII_aniso,
    LINALG::Matrix<3, 3, T> const& rcg, const int eleGID) const
{
  T I4 = 0.0;
  LINALG::Matrix<3, 3, T> AM(true);
  for (int i = 0; i < 3; ++i)
    AM(i, i) =
        anisotropyExtension_.GetStructuralTensor_stress(BaseAnisotropyExtension::GPDEFAULT, 0)(i);
  AM(0, 1) = AM(1, 0) =
      anisotropyExtension_.GetStructuralTensor_stress(BaseAnisotropyExtension::GPDEFAULT, 0)(3);
  AM(2, 1) = AM(1, 2) =
      anisotropyExtension_.GetStructuralTensor_stress(BaseAnisotropyExtension::GPDEFAULT, 0)(4);
  AM(0, 2) = AM(2, 0) =
      anisotropyExtension_.GetStructuralTensor_stress(BaseAnisotropyExtension::GPDEFAULT, 0)(5);

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
void MAT::ELASTIC::CoupAnisoExpo::AddStressAnisoPrincipal(const LINALG::Matrix<6, 1>& rcg,
    LINALG::Matrix<6, 6>& cmat, LINALG::Matrix<6, 1>& stress, Teuchos::ParameterList& params,
    const int eleGID)
{
  double I4 = 0.0;
  I4 = anisotropyExtension_.GetStructuralTensor_stress(BaseAnisotropyExtension::GPDEFAULT, 0)
           .Dot(rcg);

  double k1 = params_->k1_;
  double k2 = params_->k2_;

  if (I4 < 1.0)
  {
    k1 = params_->k1comp_;
    k2 = params_->k2comp_;
  }

  double gamma = 2. * (k1 * (I4 - 1.) * exp(k2 * (I4 - 1.) * (I4 - 1.)));
  stress.Update(gamma,
      anisotropyExtension_.GetStructuralTensor_stress(BaseAnisotropyExtension::GPDEFAULT, 0), 1.0);

  double delta =
      2. * (1. + 2. * k2 * (I4 - 1.) * (I4 - 1.)) * 2. * k1 * exp(k2 * (I4 - 1.) * (I4 - 1.));
  cmat.MultiplyNT(delta,
      anisotropyExtension_.GetStructuralTensor_stress(BaseAnisotropyExtension::GPDEFAULT, 0),
      anisotropyExtension_.GetStructuralTensor_stress(BaseAnisotropyExtension::GPDEFAULT, 0), 1.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoExpo::GetFiberVecs(
    std::vector<LINALG::Matrix<3, 1>>& fibervecs  ///< vector of all fiber vectors
)
{
  fibervecs.push_back(anisotropyExtension_.GetFiber(BaseAnisotropyExtension::GPDEFAULT, 0));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoExpo::SetFiberVecs(
    const double newgamma, const LINALG::Matrix<3, 3>& locsys, const LINALG::Matrix<3, 3>& defgrd)
{
  anisotropyExtension_.SetFiberVecs(newgamma, locsys, defgrd);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoExpo::SetFiberVecs(const LINALG::Matrix<3, 1>& fibervec)
{
  anisotropyExtension_.SetFiberVecs(fibervec);
}


// explicit instantiation of template functions
template void MAT::ELASTIC::CoupAnisoExpo::GetDerivativesAniso<double>(
    LINALG::Matrix<2, 1, double>&, LINALG::Matrix<3, 1, double>&, LINALG::Matrix<4, 1, double>&,
    LINALG::Matrix<3, 3, double> const&, const int) const;
template void MAT::ELASTIC::CoupAnisoExpo::GetDerivativesAniso<FAD>(LINALG::Matrix<2, 1, FAD>&,
    LINALG::Matrix<3, 1, FAD>&, LINALG::Matrix<4, 1, FAD>&, LINALG::Matrix<3, 3, FAD> const&,
    const int) const;
template void MAT::ELASTIC::CoupAnisoExpo::EvaluateFunc<double>(
    double&, LINALG::Matrix<3, 3, double> const&, const int) const;
template void MAT::ELASTIC::CoupAnisoExpo::EvaluateFunc<FAD>(
    FAD&, LINALG::Matrix<3, 3, FAD> const&, const int) const;
