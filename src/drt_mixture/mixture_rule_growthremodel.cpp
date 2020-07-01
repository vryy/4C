/*----------------------------------------------------------------------*/
/*! \file

\brief Mixture rule for growth and remodeling simulations with homogenized constrained mixtures

\level 3


*/
/*----------------------------------------------------------------------*/
#include "mixture_rule_growthremodel.H"
#include "../drt_mat/material_service.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/elasthyper_service.H"

MIXTURE::PAR::GrowthRemodelMixtureRule::GrowthRemodelMixtureRule(
    const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : MixtureRule(matdata), growth_type_(matdata->GetInt("GROWTH_TYPE"))
{
}

Teuchos::RCP<MIXTURE::MixtureRule> MIXTURE::PAR::GrowthRemodelMixtureRule::CreateRule()
{
  return Teuchos::rcp(new MIXTURE::GrowthRemodelMixtureRule(this));
}

MIXTURE::GrowthRemodelMixtureRule::GrowthRemodelMixtureRule(
    MIXTURE::PAR::GrowthRemodelMixtureRule* params)
    : MixtureRule(params), params_(params)
{
}

void MIXTURE::GrowthRemodelMixtureRule::PackMixtureRule(DRT::PackBuffer& data) const
{
  MixtureRule::PackMixtureRule(data);
}

void MIXTURE::GrowthRemodelMixtureRule::UnpackMixtureRule(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  MixtureRule::UnpackMixtureRule(position, data);
}

void MIXTURE::GrowthRemodelMixtureRule::Setup(Teuchos::ParameterList& params, const int eleGID)
{
  MixtureRule::Setup(params, 0);
}

void MIXTURE::GrowthRemodelMixtureRule::Update(
    LINALG::Matrix<3, 3> const& F, Teuchos::ParameterList& params, const int gp, const int eleGID)
{
  // Update base mixture rule, which also updates the constituents.
  MixtureRule::Update(F, params, gp, eleGID);

  const double dt = params.get<double>("delta time", -1.0);
  if (dt < 0.0)
  {
    dserror("The element does not write the timestep, which is fatal.");
  }

  // Evaluate inverse growth deformation gradient
  static LINALG::Matrix<3, 3> iFg;
  EvaluateInverseGrowthDeformationGradient(iFg, gp);

  for (auto const& constituent : *Constituents())
  {
    constituent->UpdateElasticPart(F, iFg, params, dt, gp, eleGID);
  }
}

void MIXTURE::GrowthRemodelMixtureRule::Evaluate(const LINALG::Matrix<3, 3>& F,
    const LINALG::Matrix<6, 1>& E_strain, Teuchos::ParameterList& params,
    LINALG::Matrix<6, 1>& S_stress, LINALG::Matrix<6, 6>& cmat, const int gp, const int eleGID)
{
  LINALG::Matrix<3, 3> iF_gM;  // growth deformation gradient

  // Evaluate growth kinematics
  EvaluateInverseGrowthDeformationGradient(iF_gM, gp);


  // define temporary matrices
  static LINALG::Matrix<6, 1> cstress;
  static LINALG::Matrix<6, 6> ccmat;

  // Iterate over all constituents and apply their contributions to the stress and linearization
  for (auto const& constituent : *Constituents())
  {
    cstress.Clear();
    ccmat.Clear();
    constituent->EvaluateElasticPart(F, iF_gM, params, cstress, ccmat, gp, eleGID);

    // Add stress contribution to global stress
    S_stress.Update(1.0, cstress, 1.0);
    cmat.Update(1.0, ccmat, 1.0);
  }
}

void MIXTURE::GrowthRemodelMixtureRule::EvaluateInverseGrowthDeformationGradient(
    LINALG::Matrix<3, 3>& iFgM, int gp) const
{
  double current_reference_density = 0.0;
  for (const auto& constituent : *Constituents())
  {
    current_reference_density += constituent->CurrentRefDensity(gp);
  }

  iFgM.Clear();

  switch (params_->growth_type_)
  {
    case PAR::GrowthRemodelMixtureRule::GROWTH_TYPE_ISOTROPIC:
    {
      for (int i = 0; i < 3; ++i)
      {
        iFgM(i, i) = std::pow(current_reference_density / GetMaterialMassDensity(), -1.0 / 3.0);
      }
      break;
    }
    case PAR::GrowthRemodelMixtureRule::GROWTH_TYPE_ANISOTROPIC:
    {
      dserror(
          "The anisotropic growth is not yet implemented. Just implement the growth deformation "
          "gradient here.");
      break;
    }
    default:
    {
      dserror(
          "The growth type %i is unknown. Please use either %i for isotropic growth or %i for "
          "anisotropic growth",
          params_->growth_type_, PAR::GrowthRemodelMixtureRule::GROWTH_TYPE_ISOTROPIC,
          PAR::GrowthRemodelMixtureRule::GROWTH_TYPE_ANISOTROPIC);
      break;
    }
  }
}
