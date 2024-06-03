/*----------------------------------------------------------------------*/
/*! \file

\brief Mixture rule for growth and remodeling simulations with homogenized constrained mixtures

\level 3


*/
/*----------------------------------------------------------------------*/
#include "4C_mixture_rule_simple.hpp"

#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_material_parameter_base.hpp"
#include "4C_mixture_constituent.hpp"
#include "4C_utils_exceptions.hpp"

#include <Epetra_ConfigDefs.h>
#include <Teuchos_RCP.hpp>

#include <algorithm>
#include <iosfwd>
#include <numeric>

FOUR_C_NAMESPACE_OPEN

MIXTURE::PAR::SimpleMixtureRule::SimpleMixtureRule(
    const Teuchos::RCP<CORE::MAT::PAR::Material>& matdata)
    : MixtureRule(matdata),
      initial_reference_density_(matdata->Get<double>("DENS")),
      mass_fractions_(matdata->Get<std::vector<double>>("MASSFRAC"))
{
  // check, whether the mass frac sums up to 1
  const double sum = std::accumulate(mass_fractions_.begin(), mass_fractions_.end(), 0.0);

  if (std::abs(1.0 - sum) > 1e-8)
    FOUR_C_THROW("Mass fractions don't sum up to 1, which is unphysical.");
}

std::unique_ptr<MIXTURE::MixtureRule> MIXTURE::PAR::SimpleMixtureRule::CreateRule()
{
  return std::unique_ptr<MIXTURE::SimpleMixtureRule>(new MIXTURE::SimpleMixtureRule(this));
}

MIXTURE::SimpleMixtureRule::SimpleMixtureRule(MIXTURE::PAR::SimpleMixtureRule* params)
    : MixtureRule(params), params_(params)
{
}

void MIXTURE::SimpleMixtureRule::Evaluate(const CORE::LINALG::Matrix<3, 3>& F,
    const CORE::LINALG::Matrix<6, 1>& E_strain, Teuchos::ParameterList& params,
    CORE::LINALG::Matrix<6, 1>& S_stress, CORE::LINALG::Matrix<6, 6>& cmat, const int gp,
    const int eleGID)
{
  // define temporary matrices
  CORE::LINALG::Matrix<6, 1> cstress;
  CORE::LINALG::Matrix<6, 6> ccmat;

  // This is the simplest mixture rule
  // Just iterate over all constituents and add all stress/cmat contributions
  for (std::size_t i = 0; i < constituents().size(); ++i)
  {
    MixtureConstituent& constituent = *constituents()[i];
    cstress.Clear();
    ccmat.Clear();
    constituent.Evaluate(F, E_strain, params, cstress, ccmat, gp, eleGID);

    // Add stress contribution to global stress
    // In this basic mixture rule, the mass fractions do not change
    double constituent_density = params_->initial_reference_density_ * params_->mass_fractions_[i];
    S_stress.Update(constituent_density, cstress, 1.0);
    cmat.Update(constituent_density, ccmat, 1.0);
  }
}
FOUR_C_NAMESPACE_CLOSE
