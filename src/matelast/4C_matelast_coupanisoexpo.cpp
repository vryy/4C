/*----------------------------------------------------------------------*/
/*! \file
\brief Implentation for an exponential strain energy function for fibers

\level 1
*/
/*----------------------------------------------------------------------*/

#include "4C_matelast_coupanisoexpo.hpp"

#include "4C_matelast_aniso_structuraltensor_strategy.hpp"
#include "4C_material_input_base.hpp"

FOUR_C_NAMESPACE_OPEN


MAT::ELASTIC::CoupAnisoExpoAnisotropyExtension::CoupAnisoExpoAnisotropyExtension(
    const int init_mode, const double gamma, const bool adapt_angle,
    const Teuchos::RCP<ELASTIC::StructuralTensorStrategyBase>& structuralTensorStrategy,
    const int fiber_id)
    : DefaultAnisotropyExtension<1>(
          init_mode, gamma, adapt_angle, structuralTensorStrategy, {fiber_id - 1})
{
}

double MAT::ELASTIC::CoupAnisoExpoAnisotropyExtension::GetScalarProduct(int gp) const
{
  return 1.0;
}

const CORE::LINALG::Matrix<3, 1>& MAT::ELASTIC::CoupAnisoExpoAnisotropyExtension::GetFiber(
    int gp) const
{
  return DefaultAnisotropyExtension<1>::GetFiber(gp, 0);
}

const CORE::LINALG::Matrix<3, 3>&
MAT::ELASTIC::CoupAnisoExpoAnisotropyExtension::GetStructuralTensor(int gp) const
{
  return DefaultAnisotropyExtension<1>::GetStructuralTensor(gp, 0);
}

const CORE::LINALG::Matrix<6, 1>&
MAT::ELASTIC::CoupAnisoExpoAnisotropyExtension::get_structural_tensor_stress(int gp) const
{
  return DefaultAnisotropyExtension<1>::get_structural_tensor_stress(gp, 0);
}

MAT::ELASTIC::PAR::CoupAnisoExpo::CoupAnisoExpo(
    const Teuchos::RCP<CORE::MAT::PAR::Material>& matdata)
    : MAT::PAR::ParameterAniso(matdata),
      MAT::ELASTIC::PAR::CoupAnisoExpoBase(matdata),
      adapt_angle_(matdata->Get<bool>("ADAPT_ANGLE")),
      fiber_id_(matdata->Get<int>("FIBER_ID"))
{
}

MAT::ELASTIC::CoupAnisoExpo::CoupAnisoExpo(MAT::ELASTIC::PAR::CoupAnisoExpo* params)
    : MAT::ELASTIC::CoupAnisoExpoBase(params),
      params_(params),
      anisotropy_extension_(params_->init_, params->gamma_, params_->adapt_angle_ != 0,
          params_->structural_tensor_strategy(), params->fiber_id_)
{
  anisotropy_extension_.register_needed_tensors(
      FiberAnisotropyExtension<1>::FIBER_VECTORS |
      FiberAnisotropyExtension<1>::STRUCTURAL_TENSOR_STRESS |
      FiberAnisotropyExtension<1>::STRUCTURAL_TENSOR);
}

void MAT::ELASTIC::CoupAnisoExpo::register_anisotropy_extensions(MAT::Anisotropy& anisotropy)
{
  anisotropy.register_anisotropy_extension(anisotropy_extension_);
}

void MAT::ELASTIC::CoupAnisoExpo::PackSummand(CORE::COMM::PackBuffer& data) const
{
  anisotropy_extension_.PackAnisotropy(data);
}

void MAT::ELASTIC::CoupAnisoExpo::UnpackSummand(
    const std::vector<char>& data, std::vector<char>::size_type& position)
{
  anisotropy_extension_.UnpackAnisotropy(data, position);
}

void MAT::ELASTIC::CoupAnisoExpo::GetFiberVecs(
    std::vector<CORE::LINALG::Matrix<3, 1>>& fibervecs  ///< vector of all fiber vectors
)
{
  if (anisotropy_extension_.FibersInitialized())
    fibervecs.push_back(anisotropy_extension_.GetFiber(BaseAnisotropyExtension::GPDEFAULT));
}

void MAT::ELASTIC::CoupAnisoExpo::SetFiberVecs(const double newgamma,
    const CORE::LINALG::Matrix<3, 3>& locsys, const CORE::LINALG::Matrix<3, 3>& defgrd)
{
  anisotropy_extension_.SetFiberVecs(newgamma, locsys, defgrd);
}

void MAT::ELASTIC::CoupAnisoExpo::SetFiberVecs(const CORE::LINALG::Matrix<3, 1>& fibervec)
{
  anisotropy_extension_.SetFiberVecs(fibervec);
}

FOUR_C_NAMESPACE_CLOSE
