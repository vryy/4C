// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_mixture_rule_simple.hpp"

#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_material_parameter_base.hpp"
#include "4C_mixture_constituent.hpp"
#include "4C_utils_exceptions.hpp"

#include <Epetra_ConfigDefs.h>

#include <algorithm>
#include <iosfwd>
#include <memory>
#include <numeric>

FOUR_C_NAMESPACE_OPEN

Mixture::PAR::SimpleMixtureRule::SimpleMixtureRule(const Core::Mat::PAR::Parameter::Data& matdata)
    : MixtureRule(matdata),
      initial_reference_density_(matdata.parameters.get<double>("DENS")),
      mass_fractions_(matdata.parameters.get<std::vector<double>>("MASSFRAC"))
{
  // check, whether the mass frac sums up to 1
  const double sum = std::accumulate(mass_fractions_.begin(), mass_fractions_.end(), 0.0);

  if (std::abs(1.0 - sum) > 1e-8)
    FOUR_C_THROW("Mass fractions don't sum up to 1, which is unphysical.");
}

std::unique_ptr<Mixture::MixtureRule> Mixture::PAR::SimpleMixtureRule::create_rule()
{
  return std::unique_ptr<Mixture::SimpleMixtureRule>(new Mixture::SimpleMixtureRule(this));
}

Mixture::SimpleMixtureRule::SimpleMixtureRule(Mixture::PAR::SimpleMixtureRule* params)
    : MixtureRule(params), params_(params)
{
}

void Mixture::SimpleMixtureRule::evaluate(const Core::LinAlg::Matrix<3, 3>& F,
    const Core::LinAlg::Matrix<6, 1>& E_strain, Teuchos::ParameterList& params,
    Core::LinAlg::Matrix<6, 1>& S_stress, Core::LinAlg::Matrix<6, 6>& cmat, const int gp,
    const int eleGID)
{
  // define temporary matrices
  Core::LinAlg::Matrix<6, 1> cstress;
  Core::LinAlg::Matrix<6, 6> ccmat;

  // This is the simplest mixture rule
  // Just iterate over all constituents and add all stress/cmat contributions
  for (std::size_t i = 0; i < constituents().size(); ++i)
  {
    MixtureConstituent& constituent = *constituents()[i];
    cstress.clear();
    ccmat.clear();
    constituent.evaluate(F, E_strain, params, cstress, ccmat, gp, eleGID);

    // Add stress contribution to global stress
    // In this basic mixture rule, the mass fractions do not change
    double constituent_density = params_->initial_reference_density_ * params_->mass_fractions_[i];
    S_stress.update(constituent_density, cstress, 1.0);
    cmat.update(constituent_density, ccmat, 1.0);
  }
}
FOUR_C_NAMESPACE_CLOSE
