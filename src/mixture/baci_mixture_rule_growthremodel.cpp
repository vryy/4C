/*----------------------------------------------------------------------*/
/*! \file

\brief Mixture rule for growth and remodeling simulations with homogenized constrained mixtures

\level 3


*/
/*----------------------------------------------------------------------*/
#include "baci_mixture_rule_growthremodel.H"
#include <Epetra_ConfigDefs.h>
#include <Epetra_SerialDenseMatrix.h>
#include <cmath>
#include <iosfwd>
#include "baci_linalg_fixedsizematrix.H"
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <algorithm>
#include "baci_utils_exceptions.H"
#include "baci_mat_par_material.H"
#include "baci_mixture_constituent.H"
#include "baci_mixture_growth_strategy.H"

// forward declarations
namespace DRT
{
  class PackBuffer;
}

MIXTURE::PAR::GrowthRemodelMixtureRule::GrowthRemodelMixtureRule(
    const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : MixtureRule(matdata),
      growth_strategy_matid_(matdata->GetInt("GROWTH_STRATEGY")),
      initial_reference_density_(matdata->GetDouble("DENS")),
      mass_fractions_(*matdata->Get<std::vector<double>>("MASSFRAC"))
{
}

std::unique_ptr<MIXTURE::MixtureRule> MIXTURE::PAR::GrowthRemodelMixtureRule::CreateRule()
{
  return std::unique_ptr<MIXTURE::GrowthRemodelMixtureRule>(
      new MIXTURE::GrowthRemodelMixtureRule(this));
}

MIXTURE::GrowthRemodelMixtureRule::GrowthRemodelMixtureRule(
    MIXTURE::PAR::GrowthRemodelMixtureRule* params)
    : MixtureRule(params), params_(params)
{
  if (params->growth_strategy_matid_ <= 0)
  {
    dserror(
        "You have not specified a growth strategy material id. Reference to the material with the "
        "growth strategy.");
  }
  growthStrategy_ = MIXTURE::PAR::MixtureGrowthStrategy::Factory(params->growth_strategy_matid_)
                        ->CreateGrowthStrategy();
}

void MIXTURE::GrowthRemodelMixtureRule::PackMixtureRule(DRT::PackBuffer& data) const
{
  MixtureRule::PackMixtureRule(data);

  growthStrategy_->PackMixtureGrowthStrategy(data);
}

void MIXTURE::GrowthRemodelMixtureRule::UnpackMixtureRule(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  MixtureRule::UnpackMixtureRule(position, data);

  growthStrategy_->UnpackMixtureGrowthStrategy(position, data);
}

void MIXTURE::GrowthRemodelMixtureRule::RegisterAnisotropyExtensions(MAT::Anisotropy& anisotropy)
{
  growthStrategy_->RegisterAnisotropyExtensions(anisotropy);
}

void MIXTURE::GrowthRemodelMixtureRule::Setup(Teuchos::ParameterList& params, const int eleGID)
{
  MixtureRule::Setup(params, 0);
}

void MIXTURE::GrowthRemodelMixtureRule::Update(CORE::LINALG::Matrix<3, 3> const& F,
    Teuchos::ParameterList& params, const int gp, const int eleGID)
{
  // Update base mixture rule, which also updates the constituents.
  MixtureRule::Update(F, params, gp, eleGID);


  // Evaluate inverse growth deformation gradient
  if (growthStrategy_->HasInelasticGrowthDeformationGradient())
  {
    const double dt = params.get<double>("delta time", -1.0);
    if (dt < 0.0)
    {
      dserror("The element does not write the timestep, which is fatal.");
    }

    // Evaluate inverse growth deformation gradient
    static CORE::LINALG::Matrix<3, 3> iFg;
    growthStrategy_->EvaluateInverseGrowthDeformationGradient(
        iFg, *this, ComputeCurrentReferenceGrowthScalar(gp), gp);

    for (const auto& constituent : Constituents())
    {
      constituent->UpdateElasticPart(F, iFg, params, dt, gp, eleGID);
    }
  }
}

void MIXTURE::GrowthRemodelMixtureRule::Evaluate(const CORE::LINALG::Matrix<3, 3>& F,
    const CORE::LINALG::Matrix<6, 1>& E_strain, Teuchos::ParameterList& params,
    CORE::LINALG::Matrix<6, 1>& S_stress, CORE::LINALG::Matrix<6, 6>& cmat, const int gp,
    const int eleGID)
{
  CORE::LINALG::Matrix<3, 3> iF_gM;  // growth deformation gradient

  if (growthStrategy_->HasInelasticGrowthDeformationGradient())
  {
    // Evaluate growth kinematics
    const double currentReferenceGrowthScalar = ComputeCurrentReferenceGrowthScalar(gp);

    growthStrategy_->EvaluateInverseGrowthDeformationGradient(
        iF_gM, *this, currentReferenceGrowthScalar, gp);
  }

  // define temporary matrices
  static CORE::LINALG::Matrix<6, 1> cstress;
  static CORE::LINALG::Matrix<6, 6> ccmat;

  // Iterate over all constituents and apply their contributions to the stress and linearization
  for (std::size_t i = 0; i < Constituents().size(); ++i)
  {
    MixtureConstituent& constituent = *Constituents()[i];
    cstress.Clear();
    ccmat.Clear();
    if (growthStrategy_->HasInelasticGrowthDeformationGradient())
      constituent.EvaluateElasticPart(F, iF_gM, params, cstress, ccmat, gp, eleGID);
    else
      constituent.Evaluate(F, E_strain, params, cstress, ccmat, gp, eleGID);


    // Add stress contribution to global stress
    const double current_ref_constituent_density = params_->initial_reference_density_ *
                                                   params_->mass_fractions_[i] *
                                                   constituent.GetGrowthScalar(gp);

    const CORE::LINALG::Matrix<1, 6> dGrowthScalarDC = constituent.GetDGrowthScalarDC(gp, eleGID);

    S_stress.Update(current_ref_constituent_density, cstress, 1.0);
    cmat.Update(current_ref_constituent_density, ccmat, 1.0);

    cmat.MultiplyNN(2.0 * params_->initial_reference_density_ * params_->mass_fractions_[i],
        cstress, dGrowthScalarDC, 1.0);
  }

  cstress.Clear();
  ccmat.Clear();


  const auto [currentReferenceGrowthScalar, dCurrentReferenceGrowthScalarDC] = std::invoke(
      [&]()
      {
        double growthScalar = 0.0;
        CORE::LINALG::Matrix<1, 6> dGrowthScalarDC(true);

        for (std::size_t i = 0; i < Constituents().size(); ++i)
        {
          MixtureConstituent& constituent = *Constituents()[i];

          growthScalar += params_->mass_fractions_[i] * constituent.GetGrowthScalar(gp);
          dGrowthScalarDC.Update(
              params_->mass_fractions_[i], constituent.GetDGrowthScalarDC(gp, eleGID), 1.0);
        }

        return std::make_tuple(growthScalar, dGrowthScalarDC);
      });

  growthStrategy_->EvaluateGrowthStressCmat(*this, currentReferenceGrowthScalar,
      dCurrentReferenceGrowthScalarDC, F, E_strain, params, cstress, ccmat, gp, eleGID);

  S_stress.Update(1.0, cstress, 1.0);
  cmat.Update(1.0, ccmat, 1.0);
}

double MIXTURE::GrowthRemodelMixtureRule::ComputeCurrentReferenceGrowthScalar(int gp) const
{
  double current_reference_growth_scalar = 0.0;
  for (std::size_t i = 0; i < Constituents().size(); ++i)
  {
    MixtureConstituent& constituent = *Constituents()[i];
    current_reference_growth_scalar +=
        params_->mass_fractions_[i] * constituent.GetGrowthScalar(gp);
  }
  return current_reference_growth_scalar;
}

double MIXTURE::GrowthRemodelMixtureRule::GetConstituentInitialReferenceMassDensity(
    const MixtureConstituent& constituent) const
{
  for (std::size_t i = 0; i < Constituents().size(); ++i)
  {
    MixtureConstituent& cur_constituent = *Constituents()[i];
    if (cur_constituent.Id() == constituent.Id())
    {
      return params_->initial_reference_density_ * params_->mass_fractions_[i];
    }
  }
  dserror("The constituent could not be found!");
  return 0.0;
}

void MIXTURE::GrowthRemodelMixtureRule::RegisterVtkOutputDataNames(
    std::unordered_map<std::string, int>& names_and_size) const
{
  names_and_size[OUTPUT_CURRENT_REFERENCE_DENSITY] = 1;
}

bool MIXTURE::GrowthRemodelMixtureRule::EvaluateVtkOutputData(
    const std::string& name, Epetra_SerialDenseMatrix& data) const
{
  if (name == OUTPUT_CURRENT_REFERENCE_DENSITY)
  {
    for (int gp = 0; gp < NumGP(); ++gp)
    {
      data(gp, 0) = ComputeCurrentReferenceGrowthScalar(gp) * params_->initial_reference_density_;
    }
    return true;
  }
  return false;
}