/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of a growth strategy for anisotropic growth

\level 3
*/
/*----------------------------------------------------------------------*/

#include "baci_mixture_growth_strategy_anisotropic.H"
#include "baci_mixture_growth_strategy.H"
#include "baci_mat_service.H"
#include "baci_lib_voigt_notation.H"
#include "baci_mat_par_material.H"
#include "baci_matelast_aniso_structuraltensor_strategy.H"

MIXTURE::PAR::AnisotropicGrowthStrategy::AnisotropicGrowthStrategy(
    const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : MIXTURE::PAR::MixtureGrowthStrategy(matdata),
      init_mode_(matdata->GetInt("INIT")),
      fiber_id_(matdata->GetInt("FIBER_ID"))
{
}

std::unique_ptr<MIXTURE::MixtureGrowthStrategy>
MIXTURE::PAR::AnisotropicGrowthStrategy::CreateGrowthStrategy()
{
  std::unique_ptr<MIXTURE::AnisotropicGrowthStrategy> strategy(
      new MIXTURE::AnisotropicGrowthStrategy(this));
  return std::move(strategy);
}

MIXTURE::AnisotropicGrowthStrategy::AnisotropicGrowthStrategy(
    MIXTURE::PAR::AnisotropicGrowthStrategy* params)
    : params_(params),
      anisotropyExtension_(params_->init_mode_, 0.0, false,
          Teuchos::rcp(new MAT::ELASTIC::StructuralTensorStrategyStandard(nullptr)),
          {params->fiber_id_ - 1})
{
  anisotropyExtension_.RegisterNeededTensors(MAT::FiberAnisotropyExtension<1>::STRUCTURAL_TENSOR);
}

void MIXTURE::AnisotropicGrowthStrategy::PackMixtureGrowthStrategy(DRT::PackBuffer& data) const
{
  MixtureGrowthStrategy::PackMixtureGrowthStrategy(data);

  anisotropyExtension_.PackAnisotropy(data);
}

void MIXTURE::AnisotropicGrowthStrategy::UnpackMixtureGrowthStrategy(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  MixtureGrowthStrategy::UnpackMixtureGrowthStrategy(position, data);

  anisotropyExtension_.UnpackAnisotropy(data, position);
}

void MIXTURE::AnisotropicGrowthStrategy::RegisterAnisotropyExtensions(MAT::Anisotropy& anisotropy)
{
  anisotropy.RegisterAnisotropyExtension(anisotropyExtension_);
}

void MIXTURE::AnisotropicGrowthStrategy::EvaluateInverseGrowthDeformationGradient(
    CORE::LINALG::Matrix<3, 3>& iFgM, const MIXTURE::MixtureRule& mixtureRule,
    double currentReferenceGrowthScalar, int gp) const
{
  CORE::LINALG::Matrix<3, 3> Id(false);
  MAT::IdentityMatrix(Id);

  iFgM.Update(1.0 / currentReferenceGrowthScalar - 1.0,
      anisotropyExtension_.GetStructuralTensor(gp, 0), 1.0, Id);
}

void MIXTURE::AnisotropicGrowthStrategy::EvaluateGrowthStressCmat(
    const MIXTURE::MixtureRule& mixtureRule, double currentReferenceGrowthScalar,
    const CORE::LINALG::Matrix<1, 6>& dCurrentReferenceGrowthScalarDC,
    const CORE::LINALG::Matrix<3, 3>& F, const CORE::LINALG::Matrix<6, 1>& E_strain,
    Teuchos::ParameterList& params, CORE::LINALG::Matrix<6, 1>& S_stress,
    CORE::LINALG::Matrix<6, 6>& cmat, const int gp, const int eleGID) const
{
  S_stress.Clear();
  cmat.Clear();
}