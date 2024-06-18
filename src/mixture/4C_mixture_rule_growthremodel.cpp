/*----------------------------------------------------------------------*/
/*! \file

\brief Mixture rule for growth and remodeling simulations with homogenized constrained mixtures

\level 3


*/
/*----------------------------------------------------------------------*/
#include "4C_mixture_rule_growthremodel.hpp"

#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_material_parameter_base.hpp"
#include "4C_mixture_constituent.hpp"
#include "4C_mixture_growth_strategy.hpp"
#include "4C_utils_exceptions.hpp"

#include <Epetra_ConfigDefs.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include <algorithm>
#include <cmath>
#include <iosfwd>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Core::Communication
{
  class PackBuffer;
}

MIXTURE::PAR::GrowthRemodelMixtureRule::GrowthRemodelMixtureRule(
    const Core::Mat::PAR::Parameter::Data& matdata)
    : MixtureRule(matdata),
      growth_strategy_matid_(matdata.parameters.get<int>("GROWTH_STRATEGY")),
      initial_reference_density_(matdata.parameters.get<double>("DENS")),
      mass_fractions_(matdata.parameters.get<std::vector<double>>("MASSFRAC"))
{
}

std::unique_ptr<MIXTURE::MixtureRule> MIXTURE::PAR::GrowthRemodelMixtureRule::create_rule()
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
    FOUR_C_THROW(
        "You have not specified a growth strategy material id. Reference to the material with the "
        "growth strategy.");
  }
  growth_strategy_ = MIXTURE::PAR::MixtureGrowthStrategy::factory(params->growth_strategy_matid_)
                         ->create_growth_strategy();
}

void MIXTURE::GrowthRemodelMixtureRule::pack_mixture_rule(
    Core::Communication::PackBuffer& data) const
{
  MixtureRule::pack_mixture_rule(data);

  growth_strategy_->pack_mixture_growth_strategy(data);
}

void MIXTURE::GrowthRemodelMixtureRule::unpack_mixture_rule(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  MixtureRule::unpack_mixture_rule(position, data);

  growth_strategy_->unpack_mixture_growth_strategy(position, data);
}

void MIXTURE::GrowthRemodelMixtureRule::register_anisotropy_extensions(Mat::Anisotropy& anisotropy)
{
  growth_strategy_->register_anisotropy_extensions(anisotropy);
}

void MIXTURE::GrowthRemodelMixtureRule::setup(Teuchos::ParameterList& params, const int eleGID)
{
  MixtureRule::setup(params, 0);
}

void MIXTURE::GrowthRemodelMixtureRule::update(Core::LinAlg::Matrix<3, 3> const& F,
    Teuchos::ParameterList& params, const int gp, const int eleGID)
{
  // Update base mixture rule, which also updates the constituents.
  MixtureRule::update(F, params, gp, eleGID);


  // Evaluate inverse growth deformation gradient
  if (growth_strategy_->has_inelastic_growth_deformation_gradient())
  {
    const double dt = params.get<double>("delta time", -1.0);
    if (dt < 0.0)
    {
      FOUR_C_THROW("The element does not write the timestep, which is fatal.");
    }

    // Evaluate inverse growth deformation gradient
    static Core::LinAlg::Matrix<3, 3> iFg;
    growth_strategy_->evaluate_inverse_growth_deformation_gradient(
        iFg, *this, compute_current_reference_growth_scalar(gp), gp);

    for (const auto& constituent : constituents())
    {
      constituent->update_elastic_part(F, iFg, params, dt, gp, eleGID);
    }
  }
}

void MIXTURE::GrowthRemodelMixtureRule::evaluate(const Core::LinAlg::Matrix<3, 3>& F,
    const Core::LinAlg::Matrix<6, 1>& E_strain, Teuchos::ParameterList& params,
    Core::LinAlg::Matrix<6, 1>& S_stress, Core::LinAlg::Matrix<6, 6>& cmat, const int gp,
    const int eleGID)
{
  Core::LinAlg::Matrix<3, 3> iF_gM;  // growth deformation gradient

  if (growth_strategy_->has_inelastic_growth_deformation_gradient())
  {
    // Evaluate growth kinematics
    const double currentReferenceGrowthScalar = compute_current_reference_growth_scalar(gp);

    growth_strategy_->evaluate_inverse_growth_deformation_gradient(
        iF_gM, *this, currentReferenceGrowthScalar, gp);
  }

  // define temporary matrices
  Core::LinAlg::Matrix<6, 1> cstress;
  Core::LinAlg::Matrix<6, 6> ccmat;

  // Iterate over all constituents and apply their contributions to the stress and linearization
  for (std::size_t i = 0; i < constituents().size(); ++i)
  {
    MixtureConstituent& constituent = *constituents()[i];
    cstress.clear();
    ccmat.clear();
    if (growth_strategy_->has_inelastic_growth_deformation_gradient())
      constituent.evaluate_elastic_part(F, iF_gM, params, cstress, ccmat, gp, eleGID);
    else
      constituent.evaluate(F, E_strain, params, cstress, ccmat, gp, eleGID);


    // Add stress contribution to global stress
    const double current_ref_constituent_density = params_->initial_reference_density_ *
                                                   params_->mass_fractions_[i] *
                                                   constituent.get_growth_scalar(gp);

    const Core::LinAlg::Matrix<1, 6> dGrowthScalarDC =
        constituent.get_d_growth_scalar_d_cg(gp, eleGID);

    S_stress.Update(current_ref_constituent_density, cstress, 1.0);
    cmat.Update(current_ref_constituent_density, ccmat, 1.0);

    cmat.MultiplyNN(2.0 * params_->initial_reference_density_ * params_->mass_fractions_[i],
        cstress, dGrowthScalarDC, 1.0);
  }

  cstress.clear();
  ccmat.clear();


  const auto [currentReferenceGrowthScalar, dCurrentReferenceGrowthScalarDC] = std::invoke(
      [&]()
      {
        double growthScalar = 0.0;
        Core::LinAlg::Matrix<1, 6> dGrowthScalarDC(true);

        for (std::size_t i = 0; i < constituents().size(); ++i)
        {
          MixtureConstituent& constituent = *constituents()[i];

          growthScalar += params_->mass_fractions_[i] * constituent.get_growth_scalar(gp);
          dGrowthScalarDC.Update(
              params_->mass_fractions_[i], constituent.get_d_growth_scalar_d_cg(gp, eleGID), 1.0);
        }

        return std::make_tuple(growthScalar, dGrowthScalarDC);
      });

  growth_strategy_->evaluate_growth_stress_cmat(*this, currentReferenceGrowthScalar,
      dCurrentReferenceGrowthScalarDC, F, E_strain, params, cstress, ccmat, gp, eleGID);

  S_stress.Update(1.0, cstress, 1.0);
  cmat.Update(1.0, ccmat, 1.0);
}

double MIXTURE::GrowthRemodelMixtureRule::compute_current_reference_growth_scalar(int gp) const
{
  double current_reference_growth_scalar = 0.0;
  for (std::size_t i = 0; i < constituents().size(); ++i)
  {
    MixtureConstituent& constituent = *constituents()[i];
    current_reference_growth_scalar +=
        params_->mass_fractions_[i] * constituent.get_growth_scalar(gp);
  }
  return current_reference_growth_scalar;
}

double MIXTURE::GrowthRemodelMixtureRule::get_constituent_initial_reference_mass_density(
    const MixtureConstituent& constituent) const
{
  for (std::size_t i = 0; i < constituents().size(); ++i)
  {
    MixtureConstituent& cur_constituent = *constituents()[i];
    if (cur_constituent.id() == constituent.id())
    {
      return params_->initial_reference_density_ * params_->mass_fractions_[i];
    }
  }
  FOUR_C_THROW("The constituent could not be found!");
  return 0.0;
}

void MIXTURE::GrowthRemodelMixtureRule::register_output_data_names(
    std::unordered_map<std::string, int>& names_and_size) const
{
  names_and_size[OUTPUT_CURRENT_REFERENCE_DENSITY] = 1;
}

bool MIXTURE::GrowthRemodelMixtureRule::evaluate_output_data(
    const std::string& name, Core::LinAlg::SerialDenseMatrix& data) const
{
  if (name == OUTPUT_CURRENT_REFERENCE_DENSITY)
  {
    for (int gp = 0; gp < num_gp(); ++gp)
    {
      data(gp, 0) =
          compute_current_reference_growth_scalar(gp) * params_->initial_reference_density_;
    }
    return true;
  }
  return false;
}
FOUR_C_NAMESPACE_CLOSE
