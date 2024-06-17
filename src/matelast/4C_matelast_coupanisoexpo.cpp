/*----------------------------------------------------------------------*/
/*! \file
\brief Implentation for an exponential strain energy function for fibers

\level 1
*/
/*----------------------------------------------------------------------*/

#include "4C_matelast_coupanisoexpo.hpp"

#include "4C_matelast_aniso_structuraltensor_strategy.hpp"
#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN


Mat::Elastic::CoupAnisoExpoAnisotropyExtension::CoupAnisoExpoAnisotropyExtension(
    const int init_mode, const double gamma, const bool adapt_angle,
    const Teuchos::RCP<Elastic::StructuralTensorStrategyBase>& structuralTensorStrategy,
    const int fiber_id)
    : DefaultAnisotropyExtension<1>(
          init_mode, gamma, adapt_angle, structuralTensorStrategy, {fiber_id - 1})
{
}

double Mat::Elastic::CoupAnisoExpoAnisotropyExtension::GetScalarProduct(int gp) const
{
  return 1.0;
}

const Core::LinAlg::Matrix<3, 1>& Mat::Elastic::CoupAnisoExpoAnisotropyExtension::get_fiber(
    int gp) const
{
  return DefaultAnisotropyExtension<1>::get_fiber(gp, 0);
}

const Core::LinAlg::Matrix<3, 3>&
Mat::Elastic::CoupAnisoExpoAnisotropyExtension::get_structural_tensor(int gp) const
{
  return DefaultAnisotropyExtension<1>::get_structural_tensor(gp, 0);
}

const Core::LinAlg::Matrix<6, 1>&
Mat::Elastic::CoupAnisoExpoAnisotropyExtension::get_structural_tensor_stress(int gp) const
{
  return DefaultAnisotropyExtension<1>::get_structural_tensor_stress(gp, 0);
}

Mat::Elastic::PAR::CoupAnisoExpo::CoupAnisoExpo(const Core::Mat::PAR::Parameter::Data& matdata)
    : Mat::PAR::ParameterAniso(matdata),
      Mat::Elastic::PAR::CoupAnisoExpoBase(matdata),
      adapt_angle_(matdata.parameters.Get<bool>("ADAPT_ANGLE")),
      fiber_id_(matdata.parameters.Get<int>("FIBER_ID"))
{
}

Mat::Elastic::CoupAnisoExpo::CoupAnisoExpo(Mat::Elastic::PAR::CoupAnisoExpo* params)
    : Mat::Elastic::CoupAnisoExpoBase(params),
      params_(params),
      anisotropy_extension_(params_->init_, params->gamma_, params_->adapt_angle_ != 0,
          params_->structural_tensor_strategy(), params->fiber_id_)
{
  anisotropy_extension_.register_needed_tensors(
      FiberAnisotropyExtension<1>::FIBER_VECTORS |
      FiberAnisotropyExtension<1>::STRUCTURAL_TENSOR_STRESS |
      FiberAnisotropyExtension<1>::STRUCTURAL_TENSOR);
}

void Mat::Elastic::CoupAnisoExpo::register_anisotropy_extensions(Mat::Anisotropy& anisotropy)
{
  anisotropy.register_anisotropy_extension(anisotropy_extension_);
}

void Mat::Elastic::CoupAnisoExpo::PackSummand(Core::Communication::PackBuffer& data) const
{
  anisotropy_extension_.pack_anisotropy(data);
}

void Mat::Elastic::CoupAnisoExpo::UnpackSummand(
    const std::vector<char>& data, std::vector<char>::size_type& position)
{
  anisotropy_extension_.unpack_anisotropy(data, position);
}

void Mat::Elastic::CoupAnisoExpo::GetFiberVecs(
    std::vector<Core::LinAlg::Matrix<3, 1>>& fibervecs  ///< vector of all fiber vectors
)
{
  if (anisotropy_extension_.fibers_initialized())
    fibervecs.push_back(anisotropy_extension_.get_fiber(BaseAnisotropyExtension::GPDEFAULT));
}

void Mat::Elastic::CoupAnisoExpo::SetFiberVecs(const double newgamma,
    const Core::LinAlg::Matrix<3, 3>& locsys, const Core::LinAlg::Matrix<3, 3>& defgrd)
{
  anisotropy_extension_.set_fiber_vecs(newgamma, locsys, defgrd);
}

void Mat::Elastic::CoupAnisoExpo::SetFiberVecs(const Core::LinAlg::Matrix<3, 1>& fibervec)
{
  anisotropy_extension_.set_fiber_vecs(fibervec);
}

FOUR_C_NAMESPACE_CLOSE
