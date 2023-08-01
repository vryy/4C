/*----------------------------------------------------------------------*/
/*! \file

\brief Mixture rule for growth and remodeling simulations with homogenized constrained mixtures

\level 3


*/
/*----------------------------------------------------------------------*/
#include "baci_mixture_rule_simple.H"
#include <Epetra_ConfigDefs.h>
#include <Teuchos_RCP.hpp>
#include <algorithm>
#include "baci_utils_exceptions.H"
#include "baci_mat_par_material.H"
#include "baci_mixture_constituent.H"
#include <iosfwd>
#include "baci_linalg_fixedsizematrix.H"

MIXTURE::PAR::SimpleMixtureRule::SimpleMixtureRule(const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : MixtureRule(matdata),
      initial_reference_density_(matdata->GetDouble("DENS")),
      mass_fractions_(*matdata->Get<std::vector<double>>("MASSFRAC"))
{
  // check, whether the mass frac sums up to 1
  double sum = 0.0;
  for (double massfrac : mass_fractions_) sum += massfrac;

  if (std::abs(1.0 - sum) > 1e-8) dserror("Mass fractions don't sum up to 1, which is unphysical.");
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
  static CORE::LINALG::Matrix<6, 1> cstress;
  static CORE::LINALG::Matrix<6, 6> ccmat;

  // This is the simplest mixture rule
  // Just iterate over all constituents and add all stress/cmat contributions
  for (std::size_t i = 0; i < Constituents().size(); ++i)
  {
    MixtureConstituent& constituent = *Constituents()[i];
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