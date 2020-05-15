/*----------------------------------------------------------------------*/
/*! \file
\brief the input line should read MAT 1 ELAST_CoupAnisoExpoTwoCoup A4 18.472 B4 16.026 A6 2.481
B6 11.120 A8 0.216 B8 11.436 GAMMA 0.0 [INIT 1] [FIB_COMP Yes] [ADAPT_ANGLE No]

\maintainer Amadeus Gebauer

\level 2

*/

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_coupanisoexpotwocoup.H"

#include "../drt_mat/matpar_material.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_mat/material_service.H"
#include "../drt_mat/anisotropy_extension.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::CoupAnisoExpoTwoCoup::CoupAnisoExpoTwoCoup(
    const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : ParameterAniso(matdata),
      A4_(matdata->GetDouble("A4")),
      B4_(matdata->GetDouble("B4")),
      A6_(matdata->GetDouble("A6")),
      B6_(matdata->GetDouble("B6")),
      A8_(matdata->GetDouble("A8")),
      B8_(matdata->GetDouble("B8")),
      gamma_(matdata->GetDouble("GAMMA")),
      init_(matdata->GetInt("INIT")),
      fib_comp_(matdata->GetInt("FIB_COMP")),
      adapt_angle_(matdata->GetInt("ADAPT_ANGLE"))
{
}


/*----------------------------------------------------------------------*
 |  Constructor                                                         |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::CoupAnisoExpoTwoCoup::CoupAnisoExpoTwoCoup(
    MAT::ELASTIC::PAR::CoupAnisoExpoTwoCoup* params)
    : params_(params), anisotropyExtension_(params_)
{
  anisotropyExtension_.RegisterNeededTensors(FiberAnisotropyExtension<2>::FIBER_VECTORS |
                                             FiberAnisotropyExtension<2>::STRUCTURAL_TENSOR_STRESS);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoExpoTwoCoup::AddStressAnisoPrincipal(const LINALG::Matrix<6, 1>& rcg,
    LINALG::Matrix<6, 6>& cmat, LINALG::Matrix<6, 1>& stress, Teuchos::ParameterList& params,
    const int eleGID)
{
  int gp = anisotropyExtension_.GetVirtualGaussPoint(params);

  LINALG::Matrix<6, 1> A1 = anisotropyExtension_.GetStructuralTensor_stress(gp, 0);
  LINALG::Matrix<6, 1> A2 = anisotropyExtension_.GetStructuralTensor_stress(gp, 1);
  LINALG::Matrix<6, 1> A1A2 = anisotropyExtension_.GetCoupledStructuralTensor_stress(gp);

  double a1a2 = anisotropyExtension_.GetCoupledScalarProduct(gp);

  double I4 = A1.Dot(rcg);
  double I6 = A2.Dot(rcg);
  double I8;
  I8 = A1A2(0) * rcg(0) + A1A2(1) * rcg(1) + A1A2(2) * rcg(2) + A1A2(3) * rcg(3) +
       A1A2(4) * rcg(4) + A1A2(5) * rcg(5);

  double A4 = params_->A4_;
  double B4 = params_->B4_;
  double A6 = params_->A6_;
  double B6 = params_->B6_;
  double A8 = params_->A8_;
  double B8 = params_->B8_;

  // check if fibers should support compression or not - if not, set the multipliers in front of
  // their strain-energy contribution to zero when the square of their stretches (fiber invariants
  // I4, I6) is smaller than one, respectively - mhv 03/14
  if ((params_->fib_comp_) == 0)
  {
    if (I4 < 1.0) A4 = 0.;
    if (I6 < 1.0) A6 = 0.;
  }

  double gamma = 2.0 * A4 * (I4 - 1.0) * exp(B4 * (I4 - 1.0) * (I4 - 1.0));
  stress.Update(gamma, A1, 1.0);
  gamma = 2.0 * A6 * (I6 - 1.0) * exp(B6 * (I6 - 1.0) * (I6 - 1.0));
  stress.Update(gamma, A2, 1.0);
  gamma = 2.0 * A8 * (I8 - a1a2) * exp(B8 * (I8 - a1a2) * (I8 - a1a2));
  stress.Update(gamma, A1A2, 1.0);

  double delta = 2.0 * (1.0 + 2.0 * B4 * (I4 - 1.0) * (I4 - 1.0)) * 2.0 * A4 *
                 exp(B4 * (I4 - 1.0) * (I4 - 1.0));
  cmat.MultiplyNT(delta, A1, A1, 1.0);
  delta = 2.0 * (1.0 + 2.0 * B6 * (I6 - 1.0) * (I6 - 1.0)) * 2.0 * A6 *
          exp(B6 * (I6 - 1.0) * (I6 - 1.0));
  cmat.MultiplyNT(delta, A2, A2, 1.0);
  delta =
      4.0 * A8 * exp(B8 * (I8 - a1a2) * (I8 - a1a2)) * (1 + 2.0 * B8 * (I8 - a1a2) * (I8 - a1a2));
  cmat.MultiplyNT(delta, A1A2, A1A2, 1.0);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoExpoTwoCoup::GetFiberVecs(
    std::vector<LINALG::Matrix<3, 1>>& fibervecs  ///< vector of all fiber vectors
)
{
  if (params_->init_ == DefaultAnisotropyExtension<2>::INIT_MODE_NODAL_FIBERS)
  {
    // This method expects constant fibers within this element but the init mode is such that
    // fibers are defined on the Gauss points
    // We therefore cannot return sth here.

    // ToDo: This may needs improvements later on if needed!
    return;
  }

  fibervecs.push_back(anisotropyExtension_.GetFiber(BaseAnisotropyExtension::GPDEFAULT, 0));
  fibervecs.push_back(anisotropyExtension_.GetFiber(BaseAnisotropyExtension::GPDEFAULT, 1));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoExpoTwoCoup::SetFiberVecs(
    const double newgamma, const LINALG::Matrix<3, 3>& locsys, const LINALG::Matrix<3, 3>& defgrd)
{
  anisotropyExtension_.SetFiberVecs(newgamma, locsys, defgrd);
}
void MAT::ELASTIC::CoupAnisoExpoTwoCoup::RegisterAnisotropyExtensions(MAT::Anisotropy& anisotropy)
{
  anisotropy.RegisterAnisotropyExtension(anisotropyExtension_);
}

MAT::ELASTIC::CoupAnisoExpoTwoCoupAnisoExtension::CoupAnisoExpoTwoCoupAnisoExtension(
    MAT::ELASTIC::PAR::CoupAnisoExpoTwoCoup* params)
    : DefaultAnisotropyExtension(params->init_, params->gamma_, params->adapt_angle_ != 0,
          params->StructuralTensorStrategy(), {0, 1})
{
}

void MAT::ELASTIC::CoupAnisoExpoTwoCoupAnisoExtension::PackAnisotropy(DRT::PackBuffer& data) const
{
  DefaultAnisotropyExtension::PackAnisotropy(data);

  DRT::ParObject::AddtoPack(data, a1a2_);
  DRT::ParObject::AddtoPack(data, A1A2_);
}

void MAT::ELASTIC::CoupAnisoExpoTwoCoupAnisoExtension::UnpackAnisotropy(
    const std::vector<char>& data, std::vector<char>::size_type& position)
{
  DefaultAnisotropyExtension::UnpackAnisotropy(data, position);

  DRT::ParObject::ExtractfromPack(position, data, a1a2_);
  DRT::ParObject::ExtractfromPack(position, data, A1A2_);
}

void MAT::ELASTIC::CoupAnisoExpoTwoCoupAnisoExtension::OnFibersInitialized()
{
  // Setup structural tensor of the coupling part
  const int fibersperele = GetFibersPerElement();

  A1A2_.resize(fibersperele);
  a1a2_.resize(fibersperele);

  for (int gp = 0; gp < fibersperele; ++gp)
  {
    LINALG::Matrix<3, 1> a1 = GetFiber(gp, 0);
    LINALG::Matrix<3, 1> a2 = GetFiber(gp, 1);
    A1A2_[gp](0) = a1(0) * a2(0);
    A1A2_[gp](1) = a1(1) * a2(1);
    A1A2_[gp](2) = a1(2) * a2(2);
    A1A2_[gp](3) = 0.5 * (a1(0) * a2(1) + a1(1) * a2(0));
    A1A2_[gp](4) = 0.5 * (a1(1) * a2(2) + a1(2) * a2(1));
    A1A2_[gp](5) = 0.5 * (a1(0) * a2(2) + a1(2) * a2(0));

    a1a2_[gp] = 0.0;
    for (int i = 0; i < 3; ++i)
    {
      a1a2_[gp] += a1(i) * a2(i);
    }
  }
}

const LINALG::Matrix<6, 1>&
MAT::ELASTIC::CoupAnisoExpoTwoCoupAnisoExtension::GetCoupledStructuralTensor_stress(int gp) const
{
  return A1A2_[gp];
}

double MAT::ELASTIC::CoupAnisoExpoTwoCoupAnisoExtension::GetCoupledScalarProduct(int gp) const
{
  return a1a2_[gp];
}
