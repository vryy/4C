/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of a growth strategy for anisotropic growth

\level 3
*/
/*----------------------------------------------------------------------*/

#include "4C_mixture_growth_strategy_anisotropic.hpp"

#include "4C_linalg_fixedsizematrix_generators.hpp"
#include "4C_linalg_fixedsizematrix_voigt_notation.hpp"
#include "4C_mat_service.hpp"
#include "4C_matelast_aniso_structuraltensor_strategy.hpp"
#include "4C_material_parameter_base.hpp"
#include "4C_mixture_growth_strategy.hpp"

FOUR_C_NAMESPACE_OPEN

MIXTURE::PAR::AnisotropicGrowthStrategy::AnisotropicGrowthStrategy(
    const Teuchos::RCP<CORE::MAT::PAR::Material>& matdata)
    : MIXTURE::PAR::MixtureGrowthStrategy(matdata),
      init_mode_(matdata->Get<int>("INIT")),
      fiber_id_(matdata->Get<int>("FIBER_ID"))
{
}

std::unique_ptr<MIXTURE::MixtureGrowthStrategy>
MIXTURE::PAR::AnisotropicGrowthStrategy::create_growth_strategy()
{
  return std::make_unique<MIXTURE::AnisotropicGrowthStrategy>(this);
}

MIXTURE::AnisotropicGrowthStrategy::AnisotropicGrowthStrategy(
    MIXTURE::PAR::AnisotropicGrowthStrategy* params)
    : params_(params),
      anisotropy_extension_(params_->init_mode_, 0.0, false,
          Teuchos::rcp(new MAT::ELASTIC::StructuralTensorStrategyStandard(nullptr)),
          {params->fiber_id_ - 1})
{
  anisotropy_extension_.register_needed_tensors(
      MAT::FiberAnisotropyExtension<1>::STRUCTURAL_TENSOR);
}

void MIXTURE::AnisotropicGrowthStrategy::pack_mixture_growth_strategy(
    CORE::COMM::PackBuffer& data) const
{
  MixtureGrowthStrategy::pack_mixture_growth_strategy(data);

  anisotropy_extension_.PackAnisotropy(data);
}

void MIXTURE::AnisotropicGrowthStrategy::unpack_mixture_growth_strategy(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  MixtureGrowthStrategy::unpack_mixture_growth_strategy(position, data);

  anisotropy_extension_.UnpackAnisotropy(data, position);
}

void MIXTURE::AnisotropicGrowthStrategy::register_anisotropy_extensions(MAT::Anisotropy& anisotropy)
{
  anisotropy.register_anisotropy_extension(anisotropy_extension_);
}

void MIXTURE::AnisotropicGrowthStrategy::evaluate_inverse_growth_deformation_gradient(
    CORE::LINALG::Matrix<3, 3>& iFgM, const MIXTURE::MixtureRule& mixtureRule,
    double currentReferenceGrowthScalar, int gp) const
{
  const CORE::LINALG::Matrix<3, 3> Id = CORE::LINALG::IdentityMatrix<3>();

  iFgM.Update(1.0 / currentReferenceGrowthScalar - 1.0,
      anisotropy_extension_.GetStructuralTensor(gp, 0), 1.0, Id);
}

void MIXTURE::AnisotropicGrowthStrategy::evaluate_growth_stress_cmat(
    const MIXTURE::MixtureRule& mixtureRule, double currentReferenceGrowthScalar,
    const CORE::LINALG::Matrix<1, 6>& dCurrentReferenceGrowthScalarDC,
    const CORE::LINALG::Matrix<3, 3>& F, const CORE::LINALG::Matrix<6, 1>& E_strain,
    Teuchos::ParameterList& params, CORE::LINALG::Matrix<6, 1>& S_stress,
    CORE::LINALG::Matrix<6, 6>& cmat, const int gp, const int eleGID) const
{
  S_stress.Clear();
  cmat.Clear();
}
FOUR_C_NAMESPACE_CLOSE
