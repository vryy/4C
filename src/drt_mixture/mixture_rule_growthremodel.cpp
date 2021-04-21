/*----------------------------------------------------------------------*/
/*! \file

\brief Mixture rule for growth and remodeling simulations with homogenized constrained mixtures

\level 3


*/
/*----------------------------------------------------------------------*/
#include "mixture_rule_growthremodel.H"
#include <Epetra_ConfigDefs.h>
#include <Epetra_SerialDenseMatrix.h>
#include <bits/c++config.h>
#include <cmath>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <algorithm>
#include <utility>
#include "../drt_lib/drt_dserror.H"
#include "../drt_mat/matpar_material.H"
#include "mixture_constituent.H"
#include <functional>

// forward declarations
namespace DRT
{
  class PackBuffer;
}

MIXTURE::PAR::GrowthRemodelMixtureRule::GrowthRemodelMixtureRule(
    const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : MixtureRule(matdata),
      growth_type_(matdata->GetInt("GROWTH_TYPE")),
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

  for (const auto& constituent : Constituents())
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
  for (std::size_t i = 0; i < Constituents().size(); ++i)
  {
    MixtureConstituent& constituent = *Constituents()[i];
    cstress.Clear();
    ccmat.Clear();
    constituent.EvaluateElasticPart(F, iF_gM, params, cstress, ccmat, gp, eleGID);

    // Add stress contribution to global stress
    const double current_ref_constituent_density = params_->initial_reference_density_ *
                                                   params_->mass_fractions_[i] *
                                                   constituent.GetGrowthScalar(gp);
    S_stress.Update(current_ref_constituent_density, cstress, 1.0);
    cmat.Update(current_ref_constituent_density, ccmat, 1.0);
  }
}

void MIXTURE::GrowthRemodelMixtureRule::EvaluateInverseGrowthDeformationGradient(
    LINALG::Matrix<3, 3>& iFgM, int gp) const
{
  const double current_reference_density = ComputeCurrentReferenceDensity(gp);

  iFgM.Clear();

  switch (params_->growth_type_)
  {
    case PAR::GrowthRemodelMixtureRule::GROWTH_TYPE_ISOTROPIC:
    {
      for (int i = 0; i < 3; ++i)
      {
        iFgM(i, i) =
            std::pow(current_reference_density / params_->initial_reference_density_, -1.0 / 3.0);
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

double MIXTURE::GrowthRemodelMixtureRule::ComputeCurrentReferenceDensity(int gp) const
{
  double current_reference_density = 0.0;
  for (std::size_t i = 0; i < Constituents().size(); ++i)
  {
    MixtureConstituent& constituent = *Constituents()[i];
    current_reference_density += params_->initial_reference_density_ * params_->mass_fractions_[i] *
                                 constituent.GetGrowthScalar(gp);
  }
  return current_reference_density;
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
      data(gp, 0) = ComputeCurrentReferenceDensity(gp);
    }
    return true;
  }
  return false;
}