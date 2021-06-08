/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of a growth strategy for anisotropic growth

\level 3
*/
/*----------------------------------------------------------------------*/

#include "mixture_growth_strategy_anisotropic.H"
#include "mixture_growth_strategy.H"
#include "../drt_mat/material_service.H"
#include "../drt_lib/voigt_notation.H"
#include "../drt_mat/matpar_material.H"
#include "../drt_matelast/elast_aniso_structuraltensor_strategy.H"

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
  return strategy;
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
    LINALG::Matrix<3, 3>& iFgM, const MIXTURE::MixtureRule& mixtureRule,
    double currentReferenceGrowthScalar, int gp) const
{
  LINALG::Matrix<3, 3> Id(false);
  MAT::IdentityMatrix(Id);

  iFgM.Update(1.0 + 1.0 / currentReferenceGrowthScalar,
      anisotropyExtension_.GetStructuralTensor(gp, 0), -1.0, Id);
}

void MIXTURE::AnisotropicGrowthStrategy::AddGrowthStressCmat(
    const MIXTURE::MixtureRule& mixtureRule, double currentReferenceGrowthScalar,
    const LINALG::Matrix<3, 3>& F, const LINALG::Matrix<6, 1>& E_strain,
    Teuchos::ParameterList& params, LINALG::Matrix<6, 1>& S_stress, LINALG::Matrix<6, 6>& cmat,
    const int gp, const int eleGID) const
{
}