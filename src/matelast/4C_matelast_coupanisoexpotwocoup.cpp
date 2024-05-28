/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of the passive material behaviour of cardiac muscle
according to Holzapfel and Ogden, "Constitutive modelling of passive myocardium", 2009.

\level 2
*/
/*----------------------------------------------------------------------*/

#include "4C_matelast_coupanisoexpotwocoup.hpp"

#include "4C_mat_anisotropy_extension.hpp"
#include "4C_mat_service.hpp"
#include "4C_material_input_base.hpp"

FOUR_C_NAMESPACE_OPEN


MAT::ELASTIC::PAR::CoupAnisoExpoTwoCoup::CoupAnisoExpoTwoCoup(
    const Teuchos::RCP<CORE::MAT::PAR::Material>& matdata)
    : ParameterAniso(matdata),
      A4_(matdata->Get<double>("A4")),
      B4_(matdata->Get<double>("B4")),
      A6_(matdata->Get<double>("A6")),
      B6_(matdata->Get<double>("B6")),
      A8_(matdata->Get<double>("A8")),
      B8_(matdata->Get<double>("B8")),
      gamma_(matdata->Get<double>("GAMMA")),
      init_(matdata->Get<int>("INIT")),
      fib_comp_(matdata->Get<bool>("FIB_COMP")),
      adapt_angle_(matdata->Get<bool>("ADAPT_ANGLE"))
{
}

MAT::ELASTIC::CoupAnisoExpoTwoCoup::CoupAnisoExpoTwoCoup(
    MAT::ELASTIC::PAR::CoupAnisoExpoTwoCoup* params)
    : params_(params), anisotropy_extension_(params_)
{
  anisotropy_extension_.register_needed_tensors(
      FiberAnisotropyExtension<2>::FIBER_VECTORS |
      FiberAnisotropyExtension<2>::STRUCTURAL_TENSOR_STRESS);
}

void MAT::ELASTIC::CoupAnisoExpoTwoCoup::PackSummand(CORE::COMM::PackBuffer& data) const
{
  anisotropy_extension_.PackAnisotropy(data);
}

void MAT::ELASTIC::CoupAnisoExpoTwoCoup::UnpackSummand(
    const std::vector<char>& data, std::vector<char>::size_type& position)
{
  anisotropy_extension_.UnpackAnisotropy(data, position);
}

void MAT::ELASTIC::CoupAnisoExpoTwoCoup::add_stress_aniso_principal(
    const CORE::LINALG::Matrix<6, 1>& rcg, CORE::LINALG::Matrix<6, 6>& cmat,
    CORE::LINALG::Matrix<6, 1>& stress, Teuchos::ParameterList& params, const int gp,
    const int eleGID)
{
  CORE::LINALG::Matrix<6, 1> A1 = anisotropy_extension_.get_structural_tensor_stress(gp, 0);
  CORE::LINALG::Matrix<6, 1> A2 = anisotropy_extension_.get_structural_tensor_stress(gp, 1);
  CORE::LINALG::Matrix<6, 1> A1A2 = anisotropy_extension_.get_coupled_structural_tensor_stress(gp);

  double a1a2 = anisotropy_extension_.get_coupled_scalar_product(gp);

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

void MAT::ELASTIC::CoupAnisoExpoTwoCoup::GetFiberVecs(
    std::vector<CORE::LINALG::Matrix<3, 1>>& fibervecs  ///< vector of all fiber vectors
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

  fibervecs.push_back(anisotropy_extension_.GetFiber(BaseAnisotropyExtension::GPDEFAULT, 0));
  fibervecs.push_back(anisotropy_extension_.GetFiber(BaseAnisotropyExtension::GPDEFAULT, 1));
}

void MAT::ELASTIC::CoupAnisoExpoTwoCoup::SetFiberVecs(const double newgamma,
    const CORE::LINALG::Matrix<3, 3>& locsys, const CORE::LINALG::Matrix<3, 3>& defgrd)
{
  anisotropy_extension_.SetFiberVecs(newgamma, locsys, defgrd);
}

void MAT::ELASTIC::CoupAnisoExpoTwoCoup::register_anisotropy_extensions(MAT::Anisotropy& anisotropy)
{
  anisotropy.register_anisotropy_extension(anisotropy_extension_);
}

MAT::ELASTIC::CoupAnisoExpoTwoCoupAnisoExtension::CoupAnisoExpoTwoCoupAnisoExtension(
    MAT::ELASTIC::PAR::CoupAnisoExpoTwoCoup* params)
    : DefaultAnisotropyExtension(params->init_, params->gamma_, params->adapt_angle_ != 0,
          params->structural_tensor_strategy(), {0, 1})
{
}

void MAT::ELASTIC::CoupAnisoExpoTwoCoupAnisoExtension::PackAnisotropy(
    CORE::COMM::PackBuffer& data) const
{
  DefaultAnisotropyExtension::PackAnisotropy(data);

  CORE::COMM::ParObject::AddtoPack(data, a1a2_);
  CORE::COMM::ParObject::AddtoPack(data, a1_a2_);
}

void MAT::ELASTIC::CoupAnisoExpoTwoCoupAnisoExtension::UnpackAnisotropy(
    const std::vector<char>& data, std::vector<char>::size_type& position)
{
  DefaultAnisotropyExtension::UnpackAnisotropy(data, position);

  CORE::COMM::ParObject::ExtractfromPack(position, data, a1a2_);
  CORE::COMM::ParObject::ExtractfromPack(position, data, a1_a2_);
}

void MAT::ELASTIC::CoupAnisoExpoTwoCoupAnisoExtension::on_fibers_initialized()
{
  // Setup structural tensor of the coupling part
  const int fibersperele = GetFibersPerElement();

  a1_a2_.resize(fibersperele);
  a1a2_.resize(fibersperele);

  for (int gp = 0; gp < fibersperele; ++gp)
  {
    CORE::LINALG::Matrix<3, 1> a1 = GetFiber(gp, 0);
    CORE::LINALG::Matrix<3, 1> a2 = GetFiber(gp, 1);
    a1_a2_[gp](0) = a1(0) * a2(0);
    a1_a2_[gp](1) = a1(1) * a2(1);
    a1_a2_[gp](2) = a1(2) * a2(2);
    a1_a2_[gp](3) = 0.5 * (a1(0) * a2(1) + a1(1) * a2(0));
    a1_a2_[gp](4) = 0.5 * (a1(1) * a2(2) + a1(2) * a2(1));
    a1_a2_[gp](5) = 0.5 * (a1(0) * a2(2) + a1(2) * a2(0));

    a1a2_[gp] = 0.0;
    for (int i = 0; i < 3; ++i)
    {
      a1a2_[gp] += a1(i) * a2(i);
    }
  }
}

const CORE::LINALG::Matrix<6, 1>&
MAT::ELASTIC::CoupAnisoExpoTwoCoupAnisoExtension::get_coupled_structural_tensor_stress(int gp) const
{
  switch (this->get_fiber_location())
  {
    case FiberLocation::ElementFibers:
      return a1_a2_[GPDEFAULT];
    case FiberLocation::GPFibers:
      return a1_a2_[gp];
    default:
      FOUR_C_THROW(
          "You have not specified, whether you want fibers on GP level or on element level.");
  }

  // Can not land here because of the FOUR_C_THROW(). Just here to ensure no compiler warning.
  std::abort();
}

double MAT::ELASTIC::CoupAnisoExpoTwoCoupAnisoExtension::get_coupled_scalar_product(int gp) const
{
  switch (this->get_fiber_location())
  {
    case FiberLocation::ElementFibers:
      return a1a2_[GPDEFAULT];
    case FiberLocation::GPFibers:
      return a1a2_[gp];
    default:
      FOUR_C_THROW(
          "You have not specified, whether you want fibers on GP level or on element level.");
  }

  // Can not land here because of the FOUR_C_THROW(). Just here to ensure no compiler warning.
  std::abort();
}

FOUR_C_NAMESPACE_CLOSE
